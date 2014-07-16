/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "AABox.h"
#include "BoundingLNode.h"

namespace Geometry
{
  template class BoundingLNode<AABox>;
  template class BoundingLNode<OBox>;
}
