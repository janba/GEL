#include "luminaire.h"

using namespace CGLA;

luminaire::luminaire(int samples) : samples_(samples)
{
    trs_.identity();
}

luminaire::~luminaire(void) {}


const Mat4x4f& luminaire::transform(void) const
{
    return trs_;
}

void luminaire::set_transform(const Mat4x4f& t)
{
    trs_ = t;
}

int luminaire::samples(void) const
{
    return samples_;
}

void luminaire::set_samples(int n)
{
    samples_ = n;
}

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
