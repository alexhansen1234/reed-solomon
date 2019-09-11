#include <iostream>
#include <iomanip>
#include <array>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include <cassert>
#include <random>
#include <chrono>
#include "galois_field.hpp"
#include "reed_solomon.hpp"
#include "rs_util.hpp"
#include "polynomial.hpp"

// constexpr size_t nbits(size_t value)
// {
//   return value == 0 ? 0 : 1 + nbits(value>>1);
// }
//
// constexpr size_t msb(size_t value)
// {
//   return value == 0 ? 0 : 1 << (nbits(value) - 1);
// }

//#define BBC

// constexpr int n = 15;
// constexpr int k = 11;
// constexpr int generator = 19;

constexpr int n = 255;
constexpr int k = 223;
constexpr int generator = 285;

constexpr int two_t = n-k;

typedef GaloisField<n,generator> gf;

int main(int argc, char ** argv)
{
  // gf::show_addition_table();
  // gf::show_multiplication_table();
  // gf::show_division_table();
  gf::show_exponentiation_table();
  gf::show_logarithm_table();

  // std::cout << "\n";

  std::default_random_engine rng;
  std::uniform_int_distribution<int> distribution(0,n);
  auto seedval = std::chrono::system_clock::now().time_since_epoch().count();

  rng.seed( seedval );

  std::cout << "Seed value: " << seedval << "\n";

  Polynomial<gf> message;
  Polynomial<gf> original;
  Polynomial<gf> encoder({1});
  Polynomial<gf> f, e, lambda, omega, dividend;

  gf primitive = 2;

  for(int i=0; i < k; ++i)
  {
    #ifdef BBC
    message.push_back( i+1 );
    #else
    message.push_back( distribution(rng) );
    #endif
  }
  for(int i=0; i < n-k; ++i)
    message.push_back(0);

  std::cout << "\nOriginal Message:\n" << message << "[" << message.coefficients.size() << "]" << "\n";

  for(int i=0; i < two_t; ++i)
  {
    encoder *= Polynomial<gf>({1, primitive.pow(i)});
  }

  std::cout << "\nEncoder Polynomial:\n" << encoder << "[" << encoder.coefficients.size() << "]" << "\n";

  auto val = message / encoder;

  std::cout << "\nError Correction:\n" << val.second << "[" << val.second.coefficients.size() << "]" << "\n";

  message += val.second;

  std::cout << "\nEncoded Message:\n" << message << "[" << message.coefficients.size() << "]" << "\n";

  original = message;

#ifdef BBC
  message[5] = 11;
  message[12] = 1;
#else
  for(int i=0; i < (n-k)/2; ++i)
  {
    int j = distribution(rng);
    if( message[ rs_util::mod(j,n) ] != j )
      message[ rs_util::mod(j,n) ] = j;
    else
      continue;
  }
#endif

  std::cout << "\nCorrupted Message:\n" << message << "[" << message.coefficients.size() << "]" << "\n";

  Polynomial<gf> syndrome;

  for(int i=two_t-1; i >= 0; --i)
  {
    Polynomial<gf> term({1, primitive.pow(i)});
    auto res = message / term;
    syndrome.push_back( res.second );
  }

  std::cout << "\nSyndromes:\n" << syndrome << "[" << syndrome.coefficients.size() << "]" << "\n";

  if( syndrome.is_zero() )
  {
    f.push_back(0);
    goto end;
  }

  dividend.push_back(1);

  for(int i=0; i < two_t; ++i)
  {
    dividend.push_back(0);
  }

  std::tie(lambda, omega) = ReedSolomonCodec<n,k,generator>::gcd(dividend, syndrome);

  std::cout << "\nLambda:\n" << lambda << "[" << lambda.coefficients.size() << "]" << "\n";
  std::cout << "\nOmega:\n" << omega << "[" << omega.coefficients.size() << "]" << "\n";

  e = ReedSolomonCodec<n,k,generator>::get_error_locator(lambda);

  std::cout << "\nError Locator:\n" << e << "[" << e.coefficients.size() << "]" << "\n";

  f = ReedSolomonCodec<n,k,generator>::get_error_polynomial(e, lambda, omega);

  std::cout << "\nError Magnitude:\n" << f << "[" << f.coefficients.size() << "]" << "\n";

  std::cout << f + message << "\n";

  std::cout << original << "\n";

end:

  auto g = f + message - original;

  std::cout << "\nTotal Error:\n";

  std::cout << g.remove_leading_zeros() << "\n";

  return !g.is_zero();
}
