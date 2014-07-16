/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "CGLA.h"

namespace
{
  unsigned int seed = 1;
  unsigned int current_rand = 1;
}

namespace CGLA
{
  void gel_srand(unsigned int s)
  {
    seed = current_rand = s;
  }

  unsigned int gel_rand(unsigned int k)
  {
    unsigned int b = 3125;
    unsigned int c = 49;
    unsigned int result = seed;

    for (;k > 0;k>>=1)
    {
      if (k & 1) result = result * b + c;
      c += b * c;
      b *= b;
    }
    return result;
  }

  unsigned int gel_rand()
  {
    const unsigned int b = 3125;
    const unsigned int c = 49;

    current_rand = current_rand*b + c;
    return current_rand;
  }
}