/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file string_utils.h
 * @brief Split a string into pieces.
 */

#ifndef __UTIL_STRING_UTILS_H
#define __UTIL_STRING_UTILS_H

#include <string>
#include <list>

namespace Util
{
  std::string trim(const std::string& s, const std::string& wspaces);
  std::string trim(const std::string& s);
  void split(const std::string& s, std::list<std::string>& result, const std::string& delim);
  void split(const std::string& s, std::list<std::string>& result);
  void trim_split(const std::string& s, std::list<std::string>& result, const std::string& delim);
  void trim_split(const std::string& s, std::list<std::string>& result);
  void get_first(std::string& s, std::string& first);
  void get_last(std::string& s, std::string& last);
}

#endif // STRING_UTILS_H