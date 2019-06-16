#include <iostream>
#include <array>
#include <utility>

int test_func(int a, int b)
{
  return a + b;
}

template <size_t n>
constexpr int msb(void)
{
  if(n == 0 || n == 1)
    return 0;
  else
  {
    int c=1;
    do
    {
      ++c;
    } while(n>>1 != 1);
    return 1<<c;
  }
}

template <class Function, size_t n, size_t irreducible_polynomial, size_t... indices>
constexpr std::array<int, n-1> generate_field_helper(Function f)
{
  return {{ f(indices, irreducible_polynomial)... }};
}

template <size_t n, size_t irreducible_polynomial>
constexpr std::array<int, n-1> generate_field(void)
{
  return generate_field_helper<test_func,n,irreducible_polynomial,std::make_index_sequence<n>{}>();
}

int main(int argc, char ** argv)
{
  auto arr = generate_field<16, 19>();
  return 0;
}
