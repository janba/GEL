#include "material.h"

material::material(float ior, float ext) : ior_(ior), extinction_(ext)
{}

material::~material(void)
{}

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
