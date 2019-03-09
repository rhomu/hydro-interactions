#ifndef RANDOM_HPP_
#define RANDOM_HPP_

#include <random>
#include "ziggurat_inline.hpp"

// -----------------------------------------------------------------------------
// init

void init_random()
{
  // init
  std::random_device rd;
  zigset(rd(), rd(), rd(), rd());
  //zigset(23, 2, 1, 377);
}

// -----------------------------------------------------------------------------
// return random real, uniform distribution

inline float random_real()
{
  return r4_uni_value();
}

inline float random_real(float min, float max)
{
  return min+random_real()*(max-min);
}

// -----------------------------------------------------------------------------
// return random real, normally distributed

inline float random_normal()
{
  return r4_nor_value();
}

inline float random_normal(float mu, float sigma)
{
  return mu + sigma*random_normal();
}

// -----------------------------------------------------------------------------
// return random uint

inline uint32_t random_uint32()
{
  return kiss_value();
}

inline uint32_t random_uint32(uint32_t lower, uint32_t upper)
{
  // this is ok when the range is small compared to UINT_MAX!!!
  return ((random_uint32() % (upper-lower+1)) + lower);
}

#endif//RANDOM_HPP_
