#include <iostream>
#include <iomanip>
#include <array>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include "rs-util.hpp"

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

  enum gf_init_type
  {
    INDEX,
    ELEMENT
  };

  static_assert( rs_util::is_power_of_two(n+1), "n must be one less than a power of two." );
  static_assert( generator < 2*(n+1), "generator must have order at most m, where 2^(m-1)-1=n.");

  GaloisField(int index)
  {
    this->index = index;

    if( index == -1 )
      this->value = 0;
    else
      this->value = exponential[mod(index,n)];
  }

  GaloisField(int in, gf_init_type t)
  {
    if( t == INDEX )
    {
      this->index = in;
      if( index == -1 )
        this->value = 0;
      else
        this->value = exponential[mod(in,n)];
    }
    else if( t == ELEMENT )
    {
      this->value = in;
      this->index = logarithm[mod(in,n+1)];
    }
  }

  size_t log()
  {
    return this->index;
  }

  GaloisField<n,generator> operator+ (const GaloisField<n,generator> & rhs)
  {
    return GaloisField<n,generator>( logarithm[ this->value ^ rhs.value ] );
  }

  GaloisField<n,generator> operator- (const GaloisField<n,generator> & rhs)
  {
    return GaloisField<n,generator>( logarithm[ this->value ^ rhs.value ] );
  }

  GaloisField<n,generator> operator* (const GaloisField<n,generator> & rhs)
  {
    if( this->index == -1 || rhs.index == -1)
      return GaloisField<n,generator>( -1, INDEX );
    else
      return GaloisField<n,generator>( mod(this->index + rhs.index,n) );
  }

  GaloisField<n,generator> operator/ (const GaloisField<n,generator> & rhs)
  {
    if( this->index == -1 || rhs.index == -1)
      return GaloisField<n,generator>( -1, INDEX );
    else
      return GaloisField<n,generator>( mod(this->index - rhs.index,n) );
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

  bool operator== (const GaloisField<n, generator> & rhs)
  {
    return this->value == rhs.value;
  }

  bool operator!= (const GaloisField<n, generator> & rhs)
  {
    return this->value != rhs.value;
  }

  template<size_t m, size_t gen>
  friend std::ostream& operator<<(std::ostream &out, const GaloisField<m,gen>& c);

  static void show_addition_table()
  {
    typedef GaloisField<n,generator> gf;

    std::cout << "GF(" << n+1 << ") addition table with generator " << generator << ":\n\n";

    for(size_t i=0; i < n+1; ++i)
    {
      gf a(i, gf::ELEMENT);
      for(size_t j=0; j < n+1; ++j)
      {
        gf b(j, gf::ELEMENT);
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
      gf a(i, gf::ELEMENT);
      for(size_t j=0; j < n+1; ++j)
      {
        gf b(j, gf::ELEMENT);
        std::cout << std::setw( static_cast<int>(log10(n)+1) ) << a * b << " ";
      }
      std::cout << "\n\n";
    }
  }

  constexpr static std::array<int, n> exponential = gf_elements<n, generator>( gf_exponential<n, generator> );
  constexpr static std::array<int, n+1> logarithm = gf_elements<n+1, generator>( gf_logarithm<n+1, exponential, 0> );

private:

  int value;
  int index;

protected:
};

template<size_t n, size_t generator_polynomial>
std::ostream & operator<<(std::ostream &out, const GaloisField<n,generator_polynomial>& c)
{
     out << c.value;
     return out;
}

typedef GaloisField<15,19> gf;

void printvec(std::string str, std::vector<gf> in)
{
  std::cout << str << "\n";
  for(auto & a : in)
    std::cout << a << " ";
  std::cout << "\n";
}

template<typename T>
int test_zero_vec( std::vector<T> & a )
{
  for(auto & b : a)
    if( b != T(-1) )
      return 0;
  return 1;
}

template<typename T>
auto conv( std::vector<T> a, std::vector<T> b )
-> std::vector<T>
{
  std::vector<T> ret;

  size_t len = a.size() + b.size() - 1;

  a.resize(len, T(-1));
  b.resize(len, T(-1));
  ret.resize(len, T(-1));

  for(size_t n=0; n < len; ++n)
  {
    for(size_t m=0; m < len; ++m)
    {
      ret[n] += a[m] * b[mod(n - m, len)];
    }
  }

  return ret;
}

template<typename T>
auto poly_div( std::vector<T> dividend, std::vector<T> divisor )
-> std::pair<std::vector<T>, std::vector<T>>
{
  size_t first_nonzero_coeff;

  // Find the first non-zero coefficient, searching from left to right
  for(first_nonzero_coeff=0; first_nonzero_coeff < divisor.size(); ++first_nonzero_coeff)
  {
    if( divisor[first_nonzero_coeff] != T(-1) )
      break;
  }

  // If the dividend is lower order than the divisor, return the zero polynomial
  // as quotient and dividend as remainder
  if( dividend.size() < divisor.size() - first_nonzero_coeff )
  {
    std::vector<T> zerovec = {T(-1)};
    return std::pair<std::vector<T>, std::vector<T>> ( zerovec, dividend );
  }

  // Store the quotient here
  std::vector<T> quotient;


  for(size_t i=0; i < dividend.size() - divisor.size() + first_nonzero_coeff + 1; ++i)
  {

    // Find the leading coefficient for the next round of division, which is
    // the remaining term in the dividend and the first term in the divisor
    T lead = dividend[i] / divisor[first_nonzero_coeff];
    quotient.push_back(lead);

    // Find the next row in the division algorithm
    for(size_t j=first_nonzero_coeff, l=0; j < divisor.size(); ++j, ++l)
    {
      dividend[i+l] += lead * divisor[j];
    }

  }
  return std::pair<std::vector<T>, std::vector<T>>(quotient, dividend);
}

template<typename T>
std::vector<T> gf_poly_sub(std::vector<T> & a, std::vector<T> & b)
{
  // In Galois Fields of characteristic 2, subtraction and addition
  // are equivalent
  std::vector<T> ret;

  if( a.size() > b.size() )
  {
    for( ; a.size() != b.size(); )
      b.insert( b.begin() , T(-1) );
  }
  else if ( b.size() > a.size() )
  {
    for( ; b.size() != a.size(); )
      a.insert( a.begin() , T(-1) );
  }

  ret.resize(a.size(), T(-1));

  for(size_t i=0; i < a.size(); ++i)
  {
    ret[i] = a[i] - b[i];
  }

  return ret;
}

template<typename T>
constexpr auto gf_poly_add = gf_poly_sub<T>;

template<typename T>
size_t gf_poly_order(std::vector<T> & a)
{
  size_t order = a.size()-1;
  for(size_t i=0; i < a.size(); ++i)
  {
    if(a[i] == T(-1))
      order--;
    else
      break;
  }
  return order;
}

template<typename T>
auto gf_gcd(std::vector<T> & a, std::vector<T> & b)
{
  std::vector<T> t0 = {T(-1)}, t1 = {T(0)};
  std::vector<T> r = b, newr = a;
  std::vector<T> temp;
  std::vector<T> quot;
  std::vector<T> rem;
  std::vector<T> cov;
  size_t terminate = gf_poly_order(a) - gf_poly_order(b) + 1;

  do
  {
    auto tup = poly_div<T>(newr, r);
    quot = tup.first;
    rem  = tup.second;

    cov = conv(t1, quot);

    if( gf_poly_order( rem ) < terminate )
    {
      break;
    }

    temp = t1;
    t1 = gf_poly_sub(t0, cov);
    t0 = temp;

    newr = r; r = rem;

  } while( 1 );

  return std::pair<std::vector<T>,std::vector<T>>(gf_poly_sub(t0,cov), rem);
}

template<typename T>
std::vector<T> gf_diff( std::vector<T> & a )
{
  std::vector<T> ret;

  int order = static_cast<int>(a.size() - 1);

  for(size_t i=0; i < a.size(); ++i)
  {
    if( (order-i) % 2 == 0 )
      ret.push_back(T(-1));
    else
      ret.push_back(a[i]);
    std::cout << "a[" << i <<"] = " << a[i] << "\n";
  }
  std::rotate(ret.begin(), ret.end()-1, ret.end());
  return ret;
}

template<typename T>
T gf_eval_poly( T input, std::vector<T> & poly )
{
  int order = static_cast<int>(poly.size() - 1);
  T val = T(-1);
  for(int i=0; i <= order; ++i)
  {
    val += poly[i] * gf( input.log() * (order - i) );
  }
  return val;
}

template<typename T, size_t n>
std::vector<T> error_locator( std::vector<T> & lambda, std::vector<T> & omega, std::vector<T> error_message )
{
  int message_order = static_cast<int>(error_message.size() - 1);

  std::vector<T> error_location;

  std::vector<T> error_polynomial;

  std::vector<T> lambda_diff = gf_diff( lambda );

  printvec("lambda:", lambda);
  printvec("lamdba':", lambda_diff);

  for(int i=message_order ; i >= 0 ; --i)
  {
    auto temp = gf_eval_poly( gf( mod( -i, n ) ) , lambda );
    error_location.push_back( temp );
  }

  printvec("error locator:", error_location);

  for(int i=0; i <= message_order; ++i)
  {
    if( error_location[i] == T(-1) )
    {
      std::cout << message_order << "\n";
      std::cout << "error order: " << -message_order + i << "\n";
      T X = gf( -message_order + i );
      T y = gf_eval_poly( X, omega );
      y /= gf_eval_poly( X, lambda_diff );
      y *= gf( message_order - i );
      error_polynomial.push_back(y);
    }
    else
    {
      error_polynomial.push_back(T(-1));
    }
  }
  return error_polynomial;
}

template<size_t N, size_t K, size_t generator>
class ReedSolomonCodec
{
  ReedSolomonCodec()
  {
      static_assert( rs-util::is_power_of_two(N+1) , "Template parameter N must be of the form 2^m-1");
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

private:

  GaloisField<N,generator> base_field;

}

int main(int argc, char ** argv)
{
  gf::show_addition_table();
  gf::show_multiplication_table();

  // Create terms of message encoder polynomial

  std::vector<gf> v0 = {gf(0), gf(0)};
  std::vector<gf> v1 = {gf(0), gf(1)};
  std::vector<gf> v2 = {gf(0), gf(2)};
  std::vector<gf> v3 = {gf(0), gf(3)};
  std::vector<gf> vc;

  // Convolve monomials to produce generator polynomial
  // RS encoding works by finding the remainder after division
  // of this polynomial by some message polynomial padded on
  // on the right with zeros, then adding the remainder to
  // the message to make the remainder zero.

  vc = conv<gf>(v0,v1);
  vc = conv<gf>(vc,v2);
  vc = conv<gf>(vc,v3);

  // Stores the message to be encoded

  std::vector<gf> message;

  std::vector<gf> original_message;

  // message = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 0, 0, 0}

  for(int i=1; i < 12; ++i)
    message.push_back(gf(i, gf::ELEMENT));
  for(int i=0; i < 4; ++i)
    message.push_back(gf(-1));

  auto val = poly_div( message, vc );

  printvec("Quotient", val.first);
  printvec("Remainder", val.second);

  for(size_t i=0; i < message.size(); ++i)
  {
    message[i] = message[i] + val.second[i];
  }

  printvec("Encoded message:", message);

  original_message = message;

  // Introduce errors, for a RS(15,11) code, n-k=2t,
  // n being the block size and k being the message size,
  // can corrected t errors.
  //
  // 15-11 = 2(2) -> two errors can be detected and corrected

  message[5] = gf(11, gf::ELEMENT);
  message[12] = gf(1, gf::ELEMENT);

  printvec("Transmitted message:", message);

  // Generate the syndrome polynomial by dividing the
  // received message by each individual factor
  // in the generator.

  auto s0 = poly_div(message, v0).second;
  auto s1 = poly_div(message, v1).second;
  auto s2 = poly_div(message, v2).second;
  auto s3 = poly_div(message, v3).second;

  // Take dividend, X^(2t) and divide by the syndrome

  std::vector<gf> syndrome = {s3.back(), s2.back(), s1.back(), s0.back()};
  std::vector<gf> dividend = {gf(0), gf(-1), gf(-1), gf(-1), gf(-1)};

  printvec("syndrome:", syndrome);

  // Find the lambda and omega polynomials
  // Lambda is the error locator polynomial,
  // and omega is used in the Forney algorithm
  // along with the derivative of lambda

  auto q = gf_gcd<gf>(dividend, syndrome);

  printvec("q:", q.first);
  printvec("p:", q.second);

  auto error = error_locator<gf, 15>(q.first, q.second, message);

  printvec("error:", error);

  std::vector<gf> recovered;

  for(size_t i=0; i < message.size(); ++i)
  {
    recovered.push_back(message[i] + error[i]);
    if(recovered[i] != original_message[i])
      std::cout << "error at position " << i << "\n";
  }

  printvec("recovered message:", recovered);

  return 0;
}
