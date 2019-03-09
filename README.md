# Hydro interactions

Simulating hydrodynamic interactions the right way.

## Building

We use cmake. Typically you would type:
```
mkdir build
cd build
cmake ..
make
```

We rely on the `boost::program_options` which must be installed prior to
building the program. We also use modern C++ features, such that you will
require a modern compiler (tested with g++-4.9).

## Running

Run examples: in main directory type
```
cd example
../build/hydro . -o .
```
Program interaction is fairly limited but you can change the simulation
parameters in `example/parameters`.
