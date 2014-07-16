/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/*!
\file LinAlgIO.h
\brief I/O functions for the LinAlg types
*/
#ifndef LINALGIO_H__HAA_AGUST_2001
#define LINALGIO_H__HAA_AGUST_2001

#include "Matrix.h"
#include <iostream>

namespace LinAlg
{

/*!
\name I/O from MatLab
Reads or writes the type to or from a m-file with a somwhat 
strict format. MatLab can 'load' the matrix by running the m-file
in question, and write to the format via the 'LinAlg.m' function.

Since MatLab uses double as the working percision, this is the only 
value of T where this interface to MatLab works. For other types, 
this functionality is good for loading, saving and debuging.

\todo write the LinAlg.m fuction.
\todo make a template overloading, of these functions.
\author Henrik Aanæs
\version Aug 2001
*/
//@{
void ToMatlab(const CMatrix& M,const std::string& VarName,const std::string& FileName="c:\\test.m",const bool append=true,const std::string& Comment = "");
void FromMatlab(CMatrix& M,const std::string& VarName,const std::string& FileName="c:\\test.m");
//@}






}
#endif // !defined(LINALGIO_H__HAA_AGUST_2001)
