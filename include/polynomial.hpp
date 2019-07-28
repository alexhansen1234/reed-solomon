#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <iostream>
#include <vector>
#include <initializer_list>
#include <type_traits>
#include <cassert>

template<typename base_field>
class Polynomial
{
public:
  Polynomial() {}

  Polynomial( base_field a )
  {
    this->coefficients.push_back(a);
  }

  Polynomial( std::initializer_list<base_field> list )
  {
    for( auto & elem : list )
      coefficients.push_back(elem);
  }

  size_t order(void) const
  {
    if( this->coefficients.size() == 1 && this->coefficients[0] == 0 )
      return 0;
    else
      return this->coefficients.size()-1;
  }

  int mod(int a, int b)
  {
    return ((a%b)+b)%b;
  }


  Polynomial<base_field> operator + (const Polynomial<base_field> & rhs )
  {
    Polynomial<base_field> ret;

    size_t result_order = std::max<size_t>( this->order(), rhs.order() );

    if( result_order == 0 && (this->coefficients[0] == 0 || rhs.coefficients[0] == 0))
    {
      ret.coefficients.resize(1,0);
      return ret;
    }
    else
    {
        ret.coefficients.resize(result_order+1, 0);
    }

    auto ret_riter = ret.coefficients.rbegin();
    auto lhs_riter = this->coefficients.rbegin();
    auto rhs_riter = rhs.coefficients.rbegin();
    auto lhs_rend = this->coefficients.rend();
    auto rhs_rend = rhs.coefficients.rend();

    for( ;lhs_riter != lhs_rend; ++lhs_riter, ++ret_riter)
    {
      *ret_riter += *lhs_riter;
    }

    ret_riter = ret.coefficients.rbegin();

    for( ;rhs_riter != rhs_rend; ++rhs_riter, ++ret_riter)
    {
      *ret_riter += *rhs_riter;
    }

    return ret;
  }

  Polynomial<base_field> operator - (const Polynomial<base_field> & rhs )
  {
    Polynomial<base_field> ret;

    size_t result_order = std::max<size_t>( this->order(), rhs.order() );

    if( result_order == 0 && (this->coefficients[0] == 0 || rhs.coefficients[0] == 0))
    {
      ret.coefficients.resize(1,0);
      return ret;
    }
    else
    {
      ret.coefficients.resize(result_order+1, 0);
    }

    auto ret_riter = ret.coefficients.rbegin();
    auto lhs_riter = this->coefficients.rbegin();
    auto rhs_riter = rhs.coefficients.rbegin();
    auto lhs_rend = this->coefficients.rend();
    auto rhs_rend = rhs.coefficients.rend();

    for( ;lhs_riter != lhs_rend; ++lhs_riter, ++ret_riter)
    {
      *ret_riter += *lhs_riter;
    }

    ret_riter = ret.coefficients.rbegin();

    for( ;rhs_riter != rhs_rend; ++rhs_riter, ++ret_riter)
    {
      *ret_riter -= *rhs_riter;
    }

    return ret;
  }

  Polynomial<base_field> operator * ( Polynomial<base_field> & rhs )
  {
    Polynomial<base_field> ret;

    size_t result_order = this->order() + rhs.order();

    if( ( this->order() == 0 || rhs.order() == 0) && ( this->coefficients.at(0) == 0 || rhs.coefficients.at(0) == 0 ))
    {
      ret.coefficients.resize(1,0);
      return ret;
    }
    else
    {
      ret.coefficients.resize(result_order+1, 0);
    }

    size_t outer_bounds = ret.order()+1;

    for( size_t n=0; n < outer_bounds; ++n )
    {
      for( size_t m=0; m < n+1; ++m)
      {
        base_field temp;

        try
        {
          temp = this->coefficients.at(m) * rhs.coefficients.at(n-m);
        }
        catch( const std::out_of_range & oor)
        {
          temp = 0;
        }

        ret.coefficients.at(n) += temp;

      }
    }

    return ret;
  }

  std::pair<Polynomial<base_field>,Polynomial<base_field>>
  operator / ( const Polynomial<base_field> & rhs ) const
  {
    Polynomial<base_field> quotient;
    Polynomial<base_field> remainder(*this);

    size_t divisor_start = 0;

    for( ; divisor_start < rhs.coefficients.size(); ++divisor_start )
      if( rhs.coefficients[divisor_start] != 0 )
        break;

    for(size_t i=0; i <= remainder.coefficients.size() - (rhs.coefficients.size() - divisor_start); ++i)
    {
      base_field lead = remainder.coefficients[i] / rhs.coefficients[divisor_start];
      quotient.coefficients.push_back( lead );

      for(size_t j=i,k=divisor_start; k < rhs.coefficients.size(); ++j, ++k)
      {
        remainder.coefficients[j] -= lead * rhs.coefficients[k];
      }
    }

    quotient.remove_leading_zeros();
    remainder.remove_leading_zeros();

    return std::pair(quotient, remainder);
  }

  Polynomial<base_field> & operator += ( Polynomial<base_field> & rhs )
  {
    *this = *this + rhs;
    return *this;
  }

  Polynomial<base_field> & operator += ( Polynomial<base_field> && rhs )
  {
    *this = *this + rhs;
    return *this;
  }


  Polynomial<base_field> & operator *= ( Polynomial<base_field> & rhs )
  {
    *this = *this * rhs;
    return *this;
  }

  Polynomial<base_field> & operator *= ( Polynomial<base_field> && rhs )
  {
    *this = *this * rhs;
    return *this;
  }

  base_field & operator[] ( size_t index )
  {
    return this->coefficients.at(index);
  }

  Polynomial<base_field> & push_back(base_field a)
  {
    this->coefficients.push_back(a);
    return *this;
  }

  Polynomial<base_field> & push_back(base_field & a)
  {
    this->coefficients.push_back(a);
    return *this;
  }

  Polynomial<base_field> & push_back( Polynomial<base_field> & a)
  {
    for(size_t i=0; i < a.coefficients.size(); ++i)
    {
      this->coefficients.push_back(a.coefficients.at(i));
    }

    return *this;
  }

  base_field eval( const base_field & a )
  {
    base_field ret = 0;
    int order = this->order();
    for(int i=0; i <= order; ++i)
      ret += this->coefficients[i] * a.pow(order-i-1);

    return ret;
  }

  template<typename b>
  friend std::ostream & operator<< (std::ostream & out, const Polynomial<b> & c);

  std::vector<base_field> coefficients;

  auto & remove_leading_zeros(void)
  {
    size_t start=0;

    for( ; start < this->coefficients.size()-1; ++start )
      if( this->coefficients[start] != 0 )
        break;

    auto begin = this->coefficients.begin();
    auto end   = this->coefficients.end();

    this->coefficients = std::vector<base_field>(begin+start,end);
    return *this;
  }

  int is_zero(void) const
  {
    for(auto & a : this->coefficients)
      if(a != 0)
        return 0;
    return 1;
  }

protected:
private:

};

template<typename base_field>
std::ostream & operator<<(std::ostream & out, const Polynomial<base_field> & c)
{
    out << "[";
    for(size_t i=0; i < c.coefficients.size()-1; ++i)
    {
      out << c.coefficients[i] << ",";
    }
    out << c.coefficients.back() << "]";

    return out;
}

#endif
