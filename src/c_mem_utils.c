#include <stddef.h>

void *cMemPointsToDbl ( double *p, int elemSize, int shift )
{
  return p + elemSize * shift;
}

size_t cGetMemAddrs ( void *p )
{
  return (size_t) p;
}

