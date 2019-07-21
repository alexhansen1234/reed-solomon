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
#include "rs-util.hpp"
#include "polynomial.hpp"

constexpr size_t nbits(size_t value)
{
  return value == 0 ? 0 : 1 + nbits(value>>1);
}

constexpr size_t msb(size_t value)
{
  return value == 0 ? 0 : 1 << (nbits(value) - 1);
}

constexpr int mod(int a, int b)
{
  return ((a%b)+b)%b;
}

constexpr size_t reduce_galois(size_t value, size_t polynomial)
{
  return value & msb(polynomial) ? value ^ polynomial : value;
}

constexpr size_t galois_element(size_t index, size_t element, size_t end_index, size_t polynomial)
{
  return index == end_index ? element : galois_element(index+1, reduce_galois(element<<1, polynomial), end_index, polynomial);
}

template<size_t n, size_t generator_polynomial>
constexpr int gf_exponential(int index)
{
  return galois_element(0, 1, index, generator_polynomial);
}

template<size_t n, std::array exponential, size_t search>
constexpr int gf_logarithm(int element)
{
  return element == 0 ? -1 : (element == exponential[search] ? search : gf_logarithm<n, exponential, mod(search+1, n+1)>(element));
}

template<size_t n, size_t generator_polynomial, class Function, size_t... indices>
constexpr auto gf_elements_helper(Function f, std::index_sequence<indices...>)
-> std::array<typename std::result_of<Function(int)>::type, n>
{
  return {{ f(indices)... }};
}

template<size_t n, size_t generator_polynomial, class Function>
constexpr auto gf_elements(Function f)
-> std::array<typename std::result_of<Function(int)>::type, n>
{
  return gf_elements_helper<n,generator_polynomial>(f, std::make_index_sequence<n>{});
}

template<typename T, size_t n>
std::vector<T> dft( std::vector<T> in )
{
  std::vector<T> out;

  out.resize(n, T(-1));
  in.resize(n, T(-1));

  for(size_t i=0; i < n; ++i)
  {
    for(size_t j=0; j < n; ++j)
    {
      out[i] += T(i*j) * in[j];
    }
  }

  return out;
}

template<typename T, size_t n>
std::vector<T> idft( std::vector<T> in )
{
  std::vector<T> out;

  out.resize(n, T(-1));
  in.resize(n, T(-1));

  for(size_t i=0; i < n; ++i)
  {
    for(size_t j=0; j < n; ++j)
    {
      out[i] += T(n - mod(i*j,n)) * in[j];
    }
  }

  return out;
}

template <size_t n, size_t generator>
class GaloisField
{
public:

  static_assert( rs_util::is_power_of_two(n+1), "n must be one less than a power of two." );
  static_assert( generator < 2*(n+1), "generator must have order at most m, where 2^(m-1)-1=n.");

  GaloisField()
  {
    this->value = 0;
  }

  GaloisField(int element)
  {
    assert(static_cast<size_t>(element) < n+1);
    assert(element >= 0);

    this->value = element;
  }

  GaloisField( const GaloisField<n,generator> & rhs )
  {
    this->value = rhs.value;
  }

  GaloisField<n,generator> operator= (size_t & rhs)
  {
    GaloisField<n,generator> ret(rhs);

    return ret;
  }

  GaloisField<n,generator> operator+ (const GaloisField<n,generator> & rhs) const
  {
    if( rhs.value == 0 )
      return *this;
    else
      return this->value ^ rhs.value;
  }

  GaloisField<n,generator> operator- (const GaloisField<n,generator> & rhs) const
  {
    if( rhs.value == 0 )
      return *this;
    else
      return this->value ^ rhs.value;
  }

  GaloisField<n,generator> operator* (const GaloisField<n,generator> & rhs) const
  {
    if( this->value == 0 || rhs.value == 0)
      return 0;
    else
    {
      auto lhs_index = GaloisField<n,generator>::logarithm[this->value];
      auto rhs_index = GaloisField<n,generator>::logarithm[rhs.value];
      return GaloisField<n,generator>::exponential[mod(lhs_index + rhs_index,n)];
    }
  }

  GaloisField<n,generator> operator/ (const GaloisField<n,generator> & rhs) const
  {
    assert( rhs.value != 0 );

    if( this->value == 0 )
      return 0;
    else
    {
      auto lhs_index = GaloisField<n,generator>::logarithm[this->value];
      auto rhs_index = GaloisField<n,generator>::logarithm[rhs.value];
      auto result_index = mod(lhs_index-rhs_index,n);
      return GaloisField<n,generator>( GaloisField<n,generator>::exponential[result_index] );
    }
  }

  GaloisField<n,generator> & operator += (const GaloisField<n,generator> & rhs)
  {
    *this = *this + rhs;
    return *this;
  }

  GaloisField<n,generator> & operator -= (const GaloisField<n,generator> & rhs)
  {
    *this = *this - rhs;
    return *this;
  }

  GaloisField<n,generator> & operator *= (const GaloisField<n,generator> & rhs)
  {
    *this = *this * rhs;
    return *this;
  }

  GaloisField<n,generator> & operator /= (const GaloisField<n,generator> & rhs)
  {
    *this = *this / rhs;
    return *this;
  }

  bool operator== ( const GaloisField<n,generator> & rhs ) const
  {
    return this->value == rhs.value;
  }

  bool operator!= ( const GaloisField<n,generator> & rhs ) const
  {
    return this->value != rhs.value;
  }

  GaloisField<n,generator> pow(const int m) const
  {
    auto this_index = GaloisField<n,generator>::logarithm[this->value];
    auto new_index = mod(this_index * m, n);
    auto new_value = GaloisField<n,generator>::exponential[new_index];
    return GaloisField<n,generator>(new_value);
  }

  template<size_t m, size_t gen>
  friend std::ostream& operator<<(std::ostream &out, const GaloisField<m,gen> & c);

  static void show_addition_table()
  {
    typedef GaloisField<n,generator> gf;

    std::cout << "GF(" << n+1 << ") addition table with generator " << generator << ":\n\n";

    for(size_t i=0; i < n+1; ++i)
    {
      gf a = i;
      for(size_t j=0; j < n+1; ++j)
      {
        gf b = j;
        std::cout << std::setw( static_cast<int>(log10(n)+1) ) << a + b << " ";
      }
      std::cout << "\n\n";
    }
  }

  static void show_multiplication_table()
  {
    typedef GaloisField<n,generator> gf;

    std::cout << "GF(" << n+1 << ") product table with generator " << generator << ":\n\n";

    for(size_t i=0; i < n+1; ++i)
    {
      gf a = i;

      for(size_t j=0; j < n+1; ++j)
      {
        gf b = j;
        std::cout << std::setw( static_cast<int>(log10(n)+1) ) << a * b << " ";
      }
      std::cout << "\n\n";
    }
  }

  static void show_exponentiation_table()
  {
    std::cout << "Exponential Lookup for GF(" << n+1 << "), primitive polynomial " << generator << ":\n";
    for(size_t i=0; i < n; ++i)
      std::cout << "a^" << i << " = " << GaloisField<n,generator>::exponential[i] << "\n";
  }

  static void show_logarithm_table()
  {
    std::cout << "Logarithm Lookup for GF(" << n+1 << "), primitive polynomial " << generator << ":\n";
    for(size_t i=1; i < n+1; ++i)
      std::cout << "log(" << i << ") = " << GaloisField<n,generator>::logarithm[i] << "\n";
  }

  constexpr static std::array<int, n> exponential = gf_elements<n, generator>( gf_exponential<n, generator> );
  constexpr static std::array<int, n+1> logarithm = gf_elements<n+1, generator>( gf_logarithm<n+1, exponential, 0> );

private:

  int value;

protected:

};

template<size_t n, size_t generator>
std::ostream & operator<<(std::ostream & out, const GaloisField<n,generator> & c)
{
     out << c.value;
     return out;
}

template<size_t N, size_t K, size_t generator>
class ReedSolomonCodec
{
public:
  typedef GaloisField<N,generator> base_field;

  ReedSolomonCodec()
  {
      static_assert( rs_util::is_power_of_two(N+1) , "Template parameter N must be of the form 2^m-1");
      static_assert( (N-K) % 2 == 0, "N-K must be an even number");
      static_assert(  N-K > 0, "N-K must be positive");
  }

  static std::vector<uint8_t> encode( std::string input )
  {
    std::vector<uint8_t> ret;
    return ret;
  }

  static std::string decode( std::vector<uint8_t> input )
  {
    std::string ret;
    return ret;
  }

  static std::pair<Polynomial<base_field>,Polynomial<base_field>>
    gcd( const Polynomial<base_field> & a, const Polynomial<base_field> & b )
  {
    Polynomial<base_field> dividend = a;
    Polynomial<base_field> divisor = b;
    Polynomial<base_field> u0({0}), u1({0}), v0({1}), v1({0});

    while(1)
    {
      Polynomial<base_field> quotient;
      Polynomial<base_field> remainder;
      std::tie(quotient,remainder) = dividend / divisor;
      dividend = divisor;
      divisor = remainder;
      u1 = v0;
      v1 += quotient*v0;

      if( divisor.order() < (N-K)/2 )
        break;
      else
      {
        std::swap(u1,u0);
        std::swap(v1,v0);
      }
    }

    return std::pair( v1, divisor );
  }

  static Polynomial<base_field> diff( const Polynomial<base_field> & a )
  {
    auto b = a;
    for(int i = b.coefficients.size()-1,j=0; i >= 0; --i,++j)
    {
      if( i % 2 == 0 )
      {
        b[j] = 0;
      }
    }
    auto c = b / Polynomial<base_field>({1,0});
    return c.first;
  }

  static Polynomial<base_field> get_error_locator( Polynomial<base_field> & lambda )
  {
    Polynomial<base_field> ret;
    base_field primitive = 2;

    for(int i=1; i <= static_cast<int>(N); ++i)
      ret.push_back( lambda.eval( primitive.pow( i ) ) );

    return ret;
  }

  static Polynomial<base_field> get_error_polynomial(
    Polynomial<base_field> & error_locations,
    Polynomial<base_field> & lambda,
    Polynomial<base_field> & omega)
  {
    Polynomial<base_field> dlambda = diff(lambda);
    Polynomial<base_field> errors;
    base_field primitive = 2;
    int order = N-1;

    std::cout << "\nd(lambda):\n" << dlambda << "[" << dlambda.coefficients.size() << "]" << "\n";

    for(int i=0; i < static_cast<int>(error_locations.order()+1); ++i)
    {
      if( error_locations[i] == 0 )
      {
        base_field q = primitive.pow( mod(order-i, N) );
        base_field p = primitive.pow( mod(i-order, N) );
        errors.push_back( q * omega.eval(p) / dlambda.eval(p) );
      }
      else
      {
        errors.push_back(0);
      }
    }

    return errors;
  }
};


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
  gf::show_addition_table();
  gf::show_multiplication_table();
  gf::show_exponentiation_table();
  gf::show_logarithm_table();

  std::cout << "\n";

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
    if( message[mod(j,n)] != j )
      message[mod(j,n)] = j;
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
