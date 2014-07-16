/**
 * @file ArgExtracter.h
 * @brief Simple functionality for getting command line parameters.
 */
/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#ifndef __UTIL_ARG_EXTRACTER__
#define __UTIL_ARG_EXTRACTER__
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <list>
#include <string>
#include <cassert>

namespace Util
{
    template<class T>
    inline T string_convert(const std::string& x)
    {
        assert(0);
        T t;
        return t;
    }
    template<> 
    inline int string_convert(const std::string& x){ 
        return std::atoi(x.c_str());}
    template<> 
    inline float string_convert(const std::string& x){ 
        return static_cast<float>(std::atof(x.c_str()));}
    template<> 
    inline double string_convert(const std::string& x){ 
        return std::atof(x.c_str());}
    template<> 
    inline std::string string_convert(const std::string& x){ 
        return x;}


    struct UpCase {void operator()(char& x) {x=toupper(x);}};

    inline void up_case_string(std::string& s)
    {
        std::for_each(s.begin(), s.end(), UpCase());
    }



    /// Extract arguments from command line parameters.
    class ArgExtracter
    {	
        std::list<std::string> avec;
        typedef std::list<std::string>::iterator LSI;

        bool extract(const std::string& argname, LSI& iter)
        {
            for(iter = avec.begin();iter != avec.end(); ++iter)
            {
                if((*iter)==argname)
                {
                    iter = avec.erase(iter);
                    return true;
                }
            }
            return false;
        }

    public:

        ArgExtracter(int argc, char **argv)
        {
            for(int i=0;i<argc; ++i)
                avec.push_back(std::string(argv[i]));
        }

        bool extract(const std::string& argname)
        {
            LSI iter;
            return extract(argname, iter);
        }

        template<class T>
        bool extract(const std::string& argname, T& val)
        {
            LSI iter;
            if(extract(argname, iter))
            {
                val = string_convert<T>(iter->c_str());
                avec.erase(iter);
                return true;
            }
            return false;
        }

        size_t no_remaining_args() const
        {
            return avec.size();
        }

        const std::string& get_last_arg() const
        {
            static const std::string emptystring("");
            if(no_remaining_args() > 0)
                return avec.back();
            return emptystring;
        }

        void get_all_args(std::vector<std::string>& args)
        {
            args = std::vector<std::string>(avec.begin(), avec.end());
        }
    };

}


#endif
