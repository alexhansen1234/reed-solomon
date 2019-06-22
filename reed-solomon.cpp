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

// template<typename T>
// int gf_test_zero_vec( std::vector<T> & a )
// {
//   for(auto & b : a)
//     if( b != T(-1) )
//       return 0;
//   return 1;
// }
//
// template<typename T>
// auto conv( std::vector<T> a, std::vector<T> b )
// -> std::vector<T>
// {
//   std::vector<T> ret;
//
//   size_t len = a.size() + b.size() - 1;
//
//   a.resize(len, T(-1));
//   b.resize(len, T(-1));
//   ret.resize(len, T(-1));
//
//   for(size_t n=0; n < len; ++n)
//   {
//     for(size_t m=0; m < len; ++m)
//     {
//       ret[n] += a[m] * b[mod(n - m, len)];
//     }
//   }
//
//   return ret;
// }
//
// template<typename T>
// auto poly_div( std::vector<T> dividend, std::vector<T> divisor )
// -> std::pair<std::vector<T>, std::vector<T>>
// {
//   size_t first_nonzero_coeff;
//
//   // Find the first non-zero coefficient, searching from left to right
//   for(first_nonzero_coeff=0; first_nonzero_coeff < divisor.size(); ++first_nonzero_coeff)
//   {
//     if( divisor[first_nonzero_coeff] != T(-1) )
//       break;
//   }
//
//   // If the dividend is lower order than the divisor, return the zero polynomial
//   // as quotient and dividend as remainder
//   if( dividend.size() < divisor.size() - first_nonzero_coeff )
//   {
//     std::vector<T> zerovec = {T(-1)};
//     return std::pair<std::vector<T>, std::vector<T>> ( zerovec, dividend );
//   }
//
//   // Store the quotient here
//   std::vector<T> quotient;
//
//
//   for(size_t i=0; i < dividend.size() - divisor.size() + first_nonzero_coeff + 1; ++i)
//   {
//     // Find the leading coefficient for the next round of division, which is
//     // the remaining term in the dividend and the first term in the divisor
//     T lead = dividend[i] / divisor[first_nonzero_coeff];
//
//     quotient.push_back(lead);
//
//     // Find the next row in the division algorithm
//     for(size_t j=first_nonzero_coeff, l=0; j < divisor.size(); ++j, ++l)
//     {
//       dividend[i+l] += lead * divisor[j];
//     }
//   }
//
//   size_t start;
//
//   for(start=0; start < dividend.size()-1; ++start)
//   {
//     if( dividend[start] != T(-1) )
//       break;
//   }
//
//   std::vector<T> remainder(dividend.begin() + start, dividend.end());
//
//   return std::pair<std::vector<T>, std::vector<T>>(quotient,remainder);
// }
//
// template<typename T>
// std::vector<T> gf_poly_sub(std::vector<T> & a, std::vector<T> & b)
// {
//   // In Galois Fields of characteristic 2, subtraction and addition
//   // are equivalent
//   std::vector<T> ret;
//
//   if( a.size() > b.size() )
//   {
//     for( ; a.size() != b.size(); )
//       b.insert( b.begin() , T(-1) );
//   }
//   else if ( b.size() > a.size() )
//   {
//     for( ; b.size() != a.size(); )
//       a.insert( a.begin() , T(-1) );
//   }
//
//   ret.resize(a.size(), T(-1));
//
//   for(size_t i=0; i < a.size(); ++i)
//   {
//     ret[i] = a[i] - b[i];
//   }
//
//   return ret;
// }
//
// template<typename T>
// constexpr auto gf_poly_add = gf_poly_sub<T>;
//
// template<typename T>
// size_t gf_poly_order(std::vector<T> & a)
// {
//   size_t order = a.size()-1;
//   for(size_t i=0; i < a.size(); ++i)
//   {
//     if(a[i] == T(-1))
//       order--;
//     else
//       break;
//   }
//   return order;
// }
//
// template<typename T>
// auto gf_poly_gcd(std::vector<T> & a, std::vector<T> & b)
// {
//   std::vector<T> t0 = {T(-1)}, t1 = {T(0)};
//   std::vector<T> r = b, newr = a;
//   std::vector<T> temp;
//   std::vector<T> quot;
//   std::vector<T> rem;
//   std::vector<T> cov;
//   size_t terminate = gf_poly_order(a) - gf_poly_order(b) + 1;
//
//   do
//   {
//     auto tup = poly_div<T>(newr, r);
//     quot = tup.first;
//     rem  = tup.second;
//
//     cov = conv(t1, quot);
//
//     if( gf_poly_order( rem ) < terminate )
//     {
//       break;
//     }
//
//     temp = t1;
//     t1 = gf_poly_sub(t0, cov);
//     t0 = temp;
//
//     newr = r; r = rem;
//
//   } while( 1 );
//
//   return std::pair<std::vector<T>,std::vector<T>>(gf_poly_sub(t0,cov), rem);
// }
//
// template<typename T>
// std::vector<T> gf_poly_diff( std::vector<T> & a )
// {
//   std::vector<T> ret;
//
//   int order = static_cast<int>(a.size() - 1);
//
//   for(size_t i=0; i < a.size(); ++i)
//   {
//     if( (order-i) % 2 == 0 )
//       ret.push_back(T(-1));
//     else
//       ret.push_back(a[i]);
//   }
//   std::rotate(ret.begin(), ret.end()-1, ret.end());
//
//   return ret;
// }
//
// template<typename T>
// T gf_eval_poly( T input, std::vector<T> & poly )
// {
//   int order = static_cast<int>(poly.size() - 1);
//   T val = T(-1);
//   for(int i=0; i <= order; ++i)
//   {
//     val += poly[i] * T( input.log() * (order - i) );
//   }
//   return val;
// }
//
// template<typename T, size_t n>
// std::vector<T> error_locator( std::vector<T> & lambda,
//   std::vector<T> & omega, std::vector<T> error_message )
// {
//   int message_order = static_cast<int>(error_message.size() - 1);
//
//   std::vector<T> error_location;
//
//   std::vector<T> error_polynomial;
//
//   std::vector<T> lambda_diff = gf_poly_diff( lambda );
//
//   printvec("lambda:", lambda);
//   printvec("d(lamdba):", lambda_diff);
//
//   for(int i=message_order ; i >= 0 ; --i)
//   {
//     auto temp = gf_eval_poly( T( mod( -i, n ) ) , lambda );
//     error_location.push_back( temp );
//   }
//
//   printvec("error locator:", error_location);
//
//   for(int i=0; i <= message_order; ++i)
//   {
//     if( error_location[i] == T(-1) )
//     {
//       std::cout << message_order << "\n";
//       std::cout << "error order: " << -message_order + i << "\n";
//       T x = T( -message_order + i );
//       T y = gf_eval_poly( x, omega );
//       y /= gf_eval_poly( x, lambda_diff );
//       y *= T( message_order - i );
//       error_polynomial.push_back(y);
//     }
//     else
//     {
//       error_polynomial.push_back(T(-1));
//     }
//   }
//   return error_polynomial;
// }

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
    Polynomial<base_field> errors({0});
    base_field primitive = 2;
    int order = N-1;

    std::cout << "\nd(lambda):\n" << dlambda << "\n";

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


// constexpr int n = 255;
// constexpr int k = 223;
// constexpr int generator = 285;

constexpr int n = 255;
constexpr int k = 223;
constexpr int generator = 285;

constexpr int two_t = n-k;

typedef GaloisField<n,generator> gf;

int main(int argc, char ** argv)
{
  // gf::show_addition_table();
  // gf::show_multiplication_table();

  std::default_random_engine rng;
  std::uniform_int_distribution<int> distribution(0,n);
  auto seedval = std::chrono::system_clock::now().time_since_epoch().count();
  rng.seed( seedval );

  std::cout << "Seed value: " << seedval << "\n";

  Polynomial<gf> message;
  Polynomial<gf> original;
  Polynomial<gf> encoder({1});

  gf primitive = 2;

  for(int i=0; i < k; ++i)
    message.push_back( distribution(rng) );
  for(int i=0; i < n-k; ++i)
    message.push_back(0);


  std::cout << "\nOriginal Message:\n" << message << "\n";

  for(int i=0; i < two_t; ++i)
  {
    encoder *= Polynomial<gf>({1, primitive.pow(i)});
  }

  std::cout << "\nEncoder Polynomial:\n" << encoder << "\n";

  auto val = message / encoder;

  std::cout << "\nError Correction:\n" << val.second << "\n";

  message += val.second;

  std::cout << "\nEncoded Message:\n" << message << "\n";

  original = message;

  for(int i=0; i < (n-k)/2; ++i)
  {
    int j = distribution(rng);
    if( message[mod(j,n)] != j )
      message[mod(j,n)] = j;
    else
      continue;
  }

  std::cout << "\nCorrupted Message:" << message << "\n";

  Polynomial<gf> syndrome;

  for(int i=two_t-1; i >= 0; --i)
  {
    Polynomial<gf> term({1, primitive.pow(i)});
    auto res = message / term;
    syndrome.push_back( res.second );
  }

  std::cout << "\nSyndromes:\n" << syndrome << "\n";

  Polynomial<gf> dividend({1});

  for(int i=0; i < two_t; ++i)
  {
    dividend.push_back(0);
  }

  // std::cout << dividend << "\n";

  Polynomial<gf> lambda, omega;

  std::tie(lambda, omega) = ReedSolomonCodec<n,k,generator>::gcd(dividend, syndrome);

  std::cout << "\nLambda:\n" << lambda << "\n";
  std::cout << "\nOmega:\n" << omega << "\n";

  auto e = ReedSolomonCodec<n,k,generator>::get_error_locator(lambda);

  std::cout << "\nError Locator:\n" << e << "\n";

  auto f = ReedSolomonCodec<n,k,generator>::get_error_polynomial(e, lambda, omega);

  std::cout << "\nError Magnitude:\n" << f << "\n";

  // std::cout << f + message << "\n";
  //
  // std::cout << original << "\n";

  auto g = f + message - original;

  std::cout << "\nTotal Error:\n";
  std::cout << g.remove_leading_zeros() << "\n";

  return g.is_zero();
}
