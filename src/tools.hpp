// tools.hpp
// some utility functions

#ifndef TOOLS_HPP_
#define TOOLS_HPP_

#include <iostream>
#include <cmath>
#include <boost/program_options.hpp>

/** Print variables from variables_map
  *
  * from: https://gist.github.com/gesquive/8673796
  */
void print_vm(const boost::program_options::variables_map& vm, unsigned padding);

namespace detail
{
  /** Convert to strig and catenate arguments */
  template<class Head>
  void inline_str_add_args(std::ostream& stream, Head&& head)
  {
    stream << std::forward<Head>(head);
  }
  /** Convert to strig and catenate arguments */
  template<class Head, class... Tail>
  void inline_str_add_args(std::ostream& stream, Head&& head, Tail&&... tail)
  {
    stream << std::forward<Head>(head);
    inline_str_add_args(stream, std::forward<Tail>(tail)...);
  }
} // namespace detail

/** Convert any number of arguments to string and catenate
 *
 * It does pretty much what is advertised. Look at the code if you want to learn
 * some pretty neat modern C++.
 * */
template<class... Args>
std::string inline_str(Args&&... args)
{
  std::stringstream s;
  detail::inline_str_add_args(s, std::forward<Args>(args)...);
  return s.str();
}

// write single value to binary stream
template<class T>
inline std::ostream& write_binary(std::ostream& stream, const T& value)
{
  return stream.write((char*)&value, sizeof(T));
}

#endif//TOOLS_HPP_
