#ifndef GALOIS_FIELD_HPP
#define GALOIS_FIELD_HPP

namespace gf_util
{
  constexpr int is_power_of_two(int a)
  {
    return (a > 0) && (a & (a-1)) == 0;
  }

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
}

constexpr size_t reduce_galois(size_t value, size_t polynomial)
{
  return value & gf_util::msb(polynomial) ? value ^ polynomial : value;
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
  return element == 0 ? -1 : (element == exponential[search] ? search : gf_logarithm<n, exponential, gf_util::mod(search+1, n+1)>(element));
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
      out[i] += T(n - gf_util::mod(i*j,n)) * in[j];
    }
  }

  return out;
}

template <size_t n, size_t generator>
class GaloisField
{
public:

  static_assert( gf_util::is_power_of_two(n+1), "n must be one less than a power of two." );
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
      return GaloisField<n,generator>::exponential[ gf_util::mod(lhs_index + rhs_index,n) ];
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
      auto result_index = gf_util::mod(lhs_index-rhs_index,n);
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
    auto new_index = gf_util::mod(this_index * m, n);
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

  static void show_division_table()
  {
    typedef GaloisField<n,generator> gf;

    std::cout << "GF(" << n+1 << ") quotient table with generator " << generator << ":\n\n";

    for(size_t i=0; i < n+1; ++i)
    {
      gf a = i;

      for(size_t j=1; j < n+1; ++j)
      {
        gf b = j;
        std::cout << std::setw( static_cast<int>(log10(n)+1) ) << a / b << " ";
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

#endif
