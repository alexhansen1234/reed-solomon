#ifndef REED_SOLOMON_HPP
#define REED_SOLOMON_HPP

#include "polynomial.hpp"
#include "galois_field.hpp"


namespace rs_util
{
  constexpr int mod(int a, int b)
  {
    return ((a%b)+b)%b;
  }

  constexpr int is_power_of_two(int a)
  {
    return (a > 0) && (a & (a-1)) == 0;
  }
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
        base_field q = primitive.pow( rs_util::mod(order-i, N) );
        base_field p = primitive.pow( rs_util::mod(i-order, N) );
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

#endif
