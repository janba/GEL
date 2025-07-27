/** @file ExceptionStandard.h
 * @brief Exceptions are not much used in CGLA, but these classes define what we throw.
 */

/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * @brief Exceptions are not much used in CGLA, but these classes define what we throw.
 * ----------------------------------------------------------------------- */

#ifndef CGLA_EXCEPTIONSTANDARD_H
#define CGLA_EXCEPTIONSTANDARD_H

#include <string_view>
#include <string>
#include <ostream>

namespace CGLA
{
/// Base class of all CGLA exceptions
class CGLAMotherException {
    std::string str;

public:
    constexpr explicit CGLAMotherException(const std::string_view s) : str(s) {}

    void print(std::ostream& os) const
    {
        os << str << "\n";
    }
};

#define CGLA_DERIVEEXCEPTION(nameoe)                                               \
class nameoe: public CGLAMotherException {                                         \
public:                                                                            \
    constexpr explicit nameoe(const std::string_view s): CGLAMotherException(s) {} \
};

}

#endif
