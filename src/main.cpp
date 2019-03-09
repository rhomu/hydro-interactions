//
// Multiparticle simulation of stochastic particles with hydrodynamic
// interactions
//
// Romain Mueller (c) 2018

#include "header.hpp"
#include "random.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// ===================================================================
// parameters

// temperature
double kBT = 1e-2;
// number of particles
int npart = 10;
// number of types
int ntypes = 1;
// total number of time steps
int nsteps = 100000;
// number of steps between writes
int ninfo = 100;
// time step
double tau = 1e-4;
// verbosity level
int verbose = 1;
// width of the output
int width = 61;
// intput file
string input;
// output directory
string output;

// strength of potential
double eps = 1e-6;
// radius of colloids
double rad = 1;
// implement Ramin's correction?
bool correction = false;

constexpr int dim = 3;
using vect = Eigen::VectorXd;
using matr = Eigen::MatrixXd;
using vec3 = Eigen::Vector3d;
using ldlt = Eigen::LDLT<matr, Eigen::Lower>;
using chol = Eigen::LLT<matr, Eigen::Lower>;

//constexpr double pi = 4*std::atan(1);
//constexpr double sqrt_6pi = sqrt(6*pi);

// =============================================================================

// Declare and parse all program options
void parse_options(int ac, char **av)
{
  // options allowed only in the command line
  opt::options_description generic("Generic options");
  generic.add_options()
    ("help,h", "produce help message")
    ("input,i", opt::value<string>(&input), "input file directory")
    ("output,o", opt::value<string>(&output), "output directory")
    ("verbose", opt::value<int>(&verbose)->implicit_value(2), "verbosity level (0, 1, 2, default=1)")
    ;

  // options allowed only in the config file
  opt::options_description config("Configuration options");
  config.add_options()
    ("npart", opt::value<int>(&npart), "number of particles")
    ("nsteps", opt::value<int>(&nsteps), "total number of time steps")
    ("ninfo", opt::value<int>(&ninfo), "number of time steps between two analyses")
    ("kBT", opt::value<double>(&kBT), "temperature")
    ("eps", opt::value<double>(&eps), "stength of harmonic potential")
    ("tau", opt::value<double>(&tau), "time step")
    ("rad", opt::value<double>(&rad), "radii of the colloids")
    ("correction", opt::value<bool>(&correction), "enable correction term?")
    ;

  // command line options
  opt::options_description cmdline_options;
  cmdline_options.add(generic);
  opt::options_description config_file_options;
  config_file_options.add(config);

  // first unnamed argument is the input directory
  opt::positional_options_description p;
  p.add("input", 1);
  p.add("output", 1);

  // the variables map
  opt::variables_map vm;

  // parse first the cmd line
  opt::store(
      opt::command_line_parser(ac, av)
      .options(cmdline_options)
      .positional(p)
      .run(), vm);
  opt::notify(vm);

  // print help msg
  if(vm.count("help"))
  {
    cout << cmdline_options << endl;
    exit(0);
  }

  // parse input file (values are not erased, such that cmd line args are 'stronger')
  if(!input.empty())
  {
    std::fstream file(input.c_str(), std::fstream::in);
    if(!file.good())
      throw inline_str("error while opening runcard file ", input);
    opt::store(opt::parse_config_file(file, config_file_options), vm);
    opt::notify(vm);
  }
  else throw inline_str("please specify an input file");

  // check output
  if(output.empty())
    throw inline_str("please specify an output directory");

  // print the simulation parameters
  if(verbose)
  {
    cout << "Run parameters" << endl;
    cout << string(width, '=') << endl;
    print_vm(vm, width);
  }
}

// write all particles positions to file
void write_frame(int t, const vect& p)
{
  ofstream file(inline_str(output, "/frame", t, ".dat"),
                ios::out | ios::binary);

  for(int a=0; a<npart; ++a)
    for(int i=0; i<dim; ++i)
      write_binary(file, p(a+i));

  file.close();
}

// write histogram
void write_histogram(const vector<double>& hist)
{
  ofstream file(inline_str(output, "/histogram.dat"),
                ios::out | ios::binary);

  for(auto& h : hist) write_binary(file, h);

  file.close();
}

// =============================================================================
// simulation

void simulate()
{
  // ---------------------------------------------------------------------------
  // init

  // total size
  const int size = dim*npart;
  // self explanatory
  const double sqrt_2kBTtau = sqrt(2*kBT*tau);
  // particles
  vect p(size);
  // forces
  vect f(size);
  // noise
  vect n(size);
  // mobility matrix
  matr M = matr::Zero(size, size);
  // derivatives of the mobility matrix
  vector<matr> dM(size, matr::Zero(size, size)),
               dE(size, matr::Zero(size, size));
  // diffusion matrix from Cholesky dec
  chol E;
  matr Einv;
  // Ramin's term
  vect R = vect::Zero(size);

  // init
  for(int A=0; A<size; ++A) p(A) = random_normal();

  // ---------------------------------------------------------------------------
  // the algo

  for(int time=0; time<=nsteps; ++time)
  {
    // store step
    if(time%ninfo == 0)
    {
      // print status
      cout << "t = " << time <<  " / " << nsteps << endl;

      // store and clear
      write_frame(time, p);
    }

    // compute mobility matrix and noise
    for(int b=0; b<size; b+=dim)
    for(int j=0; j<dim ; ++j)
    {
      for(int a=b+dim; a<size; a+=dim)
      {
        const vec3 r  { p(a + 0) - p(b + 0),
                        p(a + 1) - p(b + 1),
                        p(a + 2) - p(b + 2) };
        const auto r2 = r.dot(r);
        const auto rr = sqrt(r2);

        for(int i=0; i<dim; ++i)
        {
          M(a+i, b+j) = 0.75*r(i)*r(j)/r2/rr;

          if(not correction) continue;

          for(int k=0; k<dim; ++k)
          {
            dM[a+k](a+i, b+j) = -3/r2/r2/rr*0.75*r(i)*r(j)*r(k);
            dM[b+k](a+i, b+j) = +3/r2/r2/rr*0.75*r(i)*r(j)*r(k);
          }
          dM[a+i](a+i, b+j) += 0.75*r(j)/r2/rr;
          dM[b+i](a+i, b+j) -= 0.75*r(j)/r2/rr;
          dM[a+j](a+i, b+j) += 0.75*r(i)/r2/rr;
          dM[b+j](a+i, b+j) -= 0.75*r(i)/r2/rr;
        }

        // correct for delta(i, j)
        M(a+j, b+j) += 0.75/rr;

        if(not correction) continue;

        for(int k=0; k<dim; ++k)
        {
          dM[a+k](a+j, b+j) -= 0.75*r(k)/r2/rr;
          dM[b+k](a+j, b+j) += 0.75*r(k)/r2/rr;
        }
      }

      // A == B
      M(b+j, b+j) = 1./rad;

      // noise
      n(b+j) = sqrt_2kBTtau*random_normal();
    }

    // compute cholesky decomposition
    E.compute(M);

    // compute Ramin's term
    if(correction)
    {
      Einv = E.matrixL();
      Einv = Einv.inverse();

      // dE in terms of dM
      for(int A=0; A<size; ++A)
      {
        dE[A] = Einv*dM[A].selfadjointView<Eigen::Lower>()*Einv.transpose();

        // phi operator
        dE[A] = dE[A].selfadjointView<Eigen::Lower>();
        for(int B=0; B<size; ++B)
          dE[A](B, B) *= .5;
        dE[A] = E.matrixL()*dE[A];
      }

      // compute div E = d_i E_{ij}
      for(int A=0; A<size; ++A)
      {
        R(A) = 0.;

        for(int B=0; B<size; ++B)
            R(A) += dE[B](B, A);

        R(A) *= kBT;
      }
    }

    // compute potential force
    f = -eps*p;

    // update positions
    p += tau*M.selfadjointView<Eigen::Lower>()*f + E.matrixL()*(n + tau*R);
  }

  // write results
  //write_histogram(hist);
}

// =============================================================================

int main(int argc, char *argv[])
{
  cout << "HYDRO : Multiparticle stochastic dynamics with hydrodynamical\n"
       << "        interactions                Romain Mueller (c) 2018-9\n"
       << string(width, '=') << endl;

  try
  {
    // ========================================
    // Initialization

    // parse options (may throw)
    parse_options(argc, argv);

    if(verbose) cout << endl << "Initialization"
                     << endl << string(width, '=')
                     << endl;

    init_random();

    // ========================================
    // Running

    if(verbose) cout << endl << "Run" << endl << string(width, '=') << endl;

    // do the job
    simulate();
  }
  // error messages
  catch(const string& s) {
    cerr << argv[0] << ": " << s << endl;
    return 1;
  }
  // all the rest (mainly from boost)
  catch(const exception& e) {
    cerr << argv[0] << ": " << e.what() << endl;
    return 1;
  }

  return 0;
}
