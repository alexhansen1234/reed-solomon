namespace rs_util
{
  constexpr int is_power_of_two(int a)
  {
    return (a > 0) && (a & (a-1)) == 0;
  }

  template <size_t c>
  constexpr int get_msb()
  {
    return 1 + get_msb<(c>>1)>();
  }

  template <>
  constexpr int get_msb<1>()
  {
    return 0;
  }

  template <>
  constexpr int get_msb<0>()
  {
    return 0;
  }

  template <size_t c>
  constexpr int msb()
  {
    return get_msb<c>();
  }

  constexpr int mod(int a, int b)
  {
    return (((a % b) + b) % b);
  }

  void test_msb()
  {
    std::cout << msb<19>() << "\n";
  }

}
