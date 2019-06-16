#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

typedef struct
{
  size_t size;
  int * index_to_code;
  int * code_to_index;
} rs_table;

int is_power_of_two(int);
int msb(int);
int mod(int,int);

void init_galois_field( rs_table * tab, int n, int generator_polynomial )
{
  assert( is_power_of_two( n ) );

  assert( n * 2 > generator_polynomial );

  tab->size = n-1;

  tab->index_to_code = (int *)malloc( sizeof(int) * tab->size );

  tab->code_to_index = (int *)malloc( sizeof(int) * tab->size );

  int p = generator_polynomial;

  int lhs = msb( p );

  int rhs = p - lhs;

  int val = 1;

  for(int i=0; i < n-1; ++i)
  {
    tab->index_to_code[i] = val;
    tab->code_to_index[val] = i;

    val <<= 1;

    if( val & lhs )
    {
      val = val ^ lhs ^ rhs;
    }
  }
}

int gf_mult( rs_table * tab, int a, int b )
{
  if( a == 0 || b == 0 )
    return 0;

  int ap=a,bp=b;

  a = mod(a-1, tab->size) + 1;
  b = mod(b-1, tab->size) + 1;

  return tab->index_to_code[mod( tab->code_to_index[a] + tab->code_to_index[b], tab->size )];
}

int gf_div( rs_table * tab, int a, int b )
{
  if( a == 0 || b == 0 )
    return 0;

  a = mod(a-1, tab->size) + 1;
  b = mod(b-1, tab->size) + 1;

  return tab->index_to_code[mod( tab->code_to_index[a] - tab->code_to_index[b], tab->size )];
}

void delete_galois_field( rs_table * tab )
{
  free( tab->index_to_code );
  free( tab->code_to_index );
}

int msb( int n )
{
  int c=0;
  if( n == 0 )
    return 0;
  if( n == 1 )
    return 1;
  while( (n >>= 1) != 1 )
    c++;
  return pow(2, c+1 );
}

int is_power_of_two( int a )
{
  return (a > 0 && (a & (a - 1)) == 0);
}

int mod(int a, int b)
{
  while( a < 0 )
    a += b;
  return a % b;
}

int main(int argc, char ** argv)
{
  rs_table tab;

  init_galois_field( &tab, 256, 285 );

  for(int i=0; i < tab.size; i++)
  {
    printf("%3d\n", tab.index_to_code[i]);
  }

  delete_galois_field( &tab );
  return 0;
}
