/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/CGLA/Random.h>

#include <cassert>
#include <random>

namespace CGLA::Random
{
using namespace detail;

static_assert(std::uniform_random_bit_generator<GelPrngBase<>>);
static_assert(!std::uniform_random_bit_generator<GelPrng>);

thread_local Pcg32Random pcg32_local = {0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL};

void gel_srand(unsigned int s)
{
    pcg32_local.pcg32_seed(s, s);
}

unsigned int gel_rand(unsigned int k)
{
    auto rng = Pcg32Random(0x853c49e6748fea9bULL, k);
    return rng.pcg32_rand();
}

unsigned int gel_rand()
{
    return pcg32_local.pcg32_rand();
}
}
