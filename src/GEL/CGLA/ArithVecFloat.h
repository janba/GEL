/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file ArithVecFloat.h
 * @brief Abstract 2D floating point vector class
 */

#ifndef __CGLA__ARITHVECFLOAT_H__
#define __CGLA__ARITHVECFLOAT_H__

#include "ArithVec.h"

namespace CGLA {
    
    template<class T, class V, unsigned int N>
    class ArithVecFloat: public ArithVec<T,V,N>
    {
    public:
        
        ArithVecFloat() 
        {
#ifndef NDEBUG
            std::fill_n(this->data.begin(), N, CGLA_INIT_VALUE);
#endif
        }
        
        ArithVecFloat(T a): 
        ArithVec<T,V,N>(a) {}
        
        ArithVecFloat(T a, T b): 
        ArithVec<T,V,N>(a,b) {}
        
        ArithVecFloat(T a, T b, T c): 
        ArithVec<T,V,N>(a,b,c) {}
        
        ArithVecFloat(T a, T b, T c, T d): 
        ArithVec<T,V,N>(a,b,c,d) {}
        
        /// Compute Euclidean length.
        T length() const 
        {
            return sqrt(sqr_length(*this));
        }
        
        /// Normalize vector.
        void normalize() 
        {
            (*this) /= this->length();
        }
        
        /// Conditionally normalize vector. The condition being that the vector has non-zero length
        void cond_normalize() 
        {
            T sql = sqr_length(*this);
            if(sql > 0)
                (*this) /= sqrt(sql);
        }
        
    };
    
    /// Returns length of vector
    template<class T, class V, unsigned int N>
    inline T length(const ArithVecFloat<T,V,N>& v) 
    {
        return v.length();
    }
	
	
    /// Returns normalized vector
    template<class T, class V, unsigned int N>
    inline V normalize(const ArithVecFloat<T,V,N>& v) 
    {
        return v/v.length();
    }
    
    /// Returns normalized vector if the vector has non-zero length - otherwise the 0 vector.
    template<class T, class V, unsigned int N>
    inline V cond_normalize(const ArithVecFloat<T,V,N>& v) 
    {
        T sql = sqr_length(v);
        if(sql > 0)
            return v/sqrt(sql);
        return v*1.0;
    }
    
}

#endif

