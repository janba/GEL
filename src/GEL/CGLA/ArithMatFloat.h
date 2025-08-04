/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/** @file ArithMatFloat.h
 * Abstract matrix class
 */

#ifndef CGLA_ARITHMATFLOAT_H
#define CGLA_ARITHMATFLOAT_H

#include <numeric>
#include <GEL/CGLA/Vec.h>

namespace CGLA
{
    
    /** \brief Basic class template for matrices.
     
     In this template a matrix is defined as an array of vectors. This may
     not in all cases be the most efficient but it has the advantage that
     it is possible to use the double subscripting notation:
     
     T x = m[i][j]
     
     This template should be used through inheritance just like the
     vector template */
    template <class VVT, class HVT, class MT, unsigned int ROWS>
    class ArithMatFloat
    {
    public:
        
        /// Horizontal vector type
        using HVectorType = HVT;
        
        /// Vertical vector type
        using VVectorType = VVT;
        
        /// The type of a matrix element
        using ScalarType = typename HVT::ScalarType;
		
    protected:
        
        /// The actual contents of the matrix.
        HVT data[ROWS];
        
    protected:
        
        /// Construct 0 matrix
        constexpr ArithMatFloat()
        {
            std::fill_n(data, ROWS, HVT(CGLA_INIT_VALUE));
        }
        
        /// Construct a matrix where all entries are the same.
        explicit constexpr ArithMatFloat(ScalarType x)
        {
            std::fill_n(data, ROWS, HVT(x));
        }
        
        /// Construct a matrix where all rows are the same.
        explicit constexpr ArithMatFloat(HVT _a)
        {
            std::fill_n(data, ROWS, _a);
        }
        
        /// Construct a matrix with two rows.
        constexpr ArithMatFloat(HVT _a, HVT _b)
        {
            static_assert(ROWS==2);
            data[0] = _a;
            data[1] = _b;
        }
        
        /// Construct a matrix with three rows.
        constexpr ArithMatFloat(HVT _a, HVT _b, HVT _c)
        {
            static_assert(ROWS==3);
            data[0] = _a;
            data[1] = _b;
            data[2] = _c;
        }
        
        /// Construct a matrix with four rows.
        constexpr ArithMatFloat(HVT _a, HVT _b, HVT _c, HVT _d)
        {
            static_assert(ROWS==4);
            data[0] = _a;
            data[1] = _b;
            data[2] = _c;
            data[3] = _d;
        }
		
    public:

        /// @name C API
        /// @{

        /// Get vertical dimension of matrix
        static constexpr unsigned int get_v_dim() {return VVT::get_dim();}
        
        /// Get horizontal dimension of matrix
        static constexpr unsigned int get_h_dim() {return HVT::get_dim();}
        
        
        /** Get const pointer to data array.
         This function may be useful when interfacing with some other API
         such as OpenGL (TM). */
        constexpr const ScalarType* get() const
        {
            return data[0].get();
        }
        
        /** Get pointer to data array.
         This function may be useful when interfacing with some other API
         such as OpenGL (TM). */
        constexpr ScalarType* get()
        {
            return data[0].get();
        }

        /// @}
        /// @name Index operators
        /// @{
        
        /// Const index operator. Returns i'th row of matrix.
        constexpr const HVT& operator [] ( unsigned int i ) const
        {
            assert(i<ROWS);
            return data[i];
        }
        
        /// Non-const index operator. Returns i'th row of matrix.
        constexpr HVT& operator [] ( unsigned int i )
        {
            assert(i<ROWS);
            return data[i];
        }

        /// @}
        /// @name Partial equality
        /// @{
        
        /// Equality operator.
        constexpr bool operator==(const MT& v) const
        {
            return std::inner_product(data, &data[ROWS], &v[0], true,
                                      std::logical_and<>(), std::equal_to<HVT>());
        }
        
        /// Inequality operator.
        constexpr bool operator!=(const MT& v) const
        {
            return !(*this==v);
        }
        
        /// @}
        /// @name Arithmetic operations
        /// @{
        
        /// Multiply scalar onto matrix. All entries are multiplied by scalar.
        constexpr MT operator *(ScalarType k) const
        {
            MT v_new;
            for(std::size_t i = 0; i < ROWS; ++i)
                v_new.data[i] = data[i] * k;
            return v_new;
        }
        
        /// Divide all entries in matrix by scalar.
        constexpr MT operator / (ScalarType k) const
        {
            MT v_new;
            for(std::size_t i = 0; i < ROWS; ++i)
                v_new.data[i] = data[i] / k;
            return v_new;
        }
        
        /// Assignment multiplication of matrix by scalar.
        constexpr const MT& operator *=(ScalarType k)
        {
            for(auto& x: data)
                x *= k;
            return static_cast<const MT&>(*this);
        }
        
        /// Assignment division of matrix by scalar.
        constexpr const MT& operator /=(ScalarType k)
        {
            for(auto& x: data)
                x /= k;
            return static_cast<const MT&>(*this);
        }
        
        //----------------------------------------------------------------------
        
        /// Add two matrices.
        constexpr MT operator + (const MT& m1) const
        {
            MT v_new;
            std::transform(data, &data[ROWS], &m1[0], &v_new[0], std::plus<HVT>());
            return v_new;
        }
        
        /// Subtract two matrices.
        constexpr MT operator - (const MT& m1) const
        {
            MT v_new;
            std::transform(data, &data[ROWS], &m1[0], &v_new[0], std::minus<HVT>());
            return v_new;
        }
        
        /// Assigment addition of matrices.
        constexpr const MT& operator +=(const MT& v)
        {
            std::transform(data, &data[ROWS], &v[0], data, std::plus<HVT>());
            return static_cast<const MT&>(*this);
        }
        
        /// Assigment subtraction of matrices.
        constexpr const MT& operator -=(const MT& v)
        {
            std::transform(data, &data[ROWS], &v[0], data, std::minus<HVT>());
            return static_cast<const MT&>(*this);
        }
        
        //----------------------------------------------------------------------
        
        /// Negate matrix.
        constexpr MT operator - () const
        {
            MT v_new;
            std::transform(data, &data[ROWS], &v_new[0], std::negate<HVT>());
            return v_new;
        }

        /// @}
    };

    /// @}
    /// @name Arithmetic with scalars
    /// @related ArithMatFloat
    /// @{
    
    /// Multiply scalar onto matrix
    template <class VVT, class HVT, class MT, unsigned int ROWS>
    constexpr MT operator * (double k, const ArithMatFloat<VVT,HVT,MT,ROWS>& v)
    {
        return v * k;
    }
    
    /// Multiply scalar onto matrix
    template <class VVT, class HVT, class MT, unsigned int ROWS>
    constexpr MT operator * (float k, const ArithMatFloat<VVT,HVT,MT,ROWS>& v)
    {
        return v * k;
    }
    
    /// Multiply scalar onto matrix
    template <class VVT, class HVT, class MT, unsigned int ROWS>
    constexpr MT operator * (int k, const ArithMatFloat<VVT,HVT,MT,ROWS>& v)
    {
        return v * k;
    }
    
    /// Multiply vector onto matrix
    template <class VVT, class HVT, class MT, unsigned int ROWS>
    constexpr VVT operator*(const ArithMatFloat<VVT,HVT,MT,ROWS>& m,const HVT& v)
    {
        VVT v2;
        for(unsigned int i=0;i<ROWS;i++) v2[i] = dot(m[i], v);
            return v2;
    }

    /// @}
    /// @name Matrix operations
    /// @related ArithMatFloat
    /// @{
    
#ifndef WIN32
    /// @brief Multiply two arbitrary matrices.
    /// @details In principle, this function could return a matrix, but in general
    /// the new matrix will be of a type that is different from either of
    /// the two matrices that are multiplied together. We do not want to
    /// return an ArithMatFloat - so it seems best to let the return value be
    /// a reference arg.
    ///
    /// This template can only be instantiated if the dimensions of the
    /// matrices match -- i.e. if the multiplication can actually be
    /// carried out. This is more type safe than the win32 version below.
    /// @related ArithMatFloat
    template <class VVT, class HVT,
    class HV1T, class VV2T,
    class MT1, class MT2, class MT,
    unsigned int ROWS1, unsigned int ROWS2>
    constexpr void mul(const ArithMatFloat<VVT,HV1T,MT1,ROWS1>& m1,
                    const ArithMatFloat<VV2T,HVT,MT2,ROWS2>& m2,
                    ArithMatFloat<VVT,HVT,MT,ROWS1>& m)
    {
        constexpr auto cols = ArithMatFloat<VVT,HVT,MT,ROWS1>::get_h_dim();
        for(unsigned int i=0;i<ROWS1;i++)
            for(unsigned int j=0;j<cols;j++)
            {
                m[i][j] = 0;
                for(unsigned int k=0;k<ROWS2;k++)
                    m[i][j] += m1[i][k] * m2[k][j];
            }
    }
    
    
    /// Transpose. See the discussion on mul if you are curious as to why
    /// I don't simply return the transpose.
    /// @sa CGLA::mul
    /// @related ArithMatFloat
    template <class VVT, class HVT, class M1T, class M2T, unsigned int ROWS, unsigned int COLS>
    constexpr void transpose(const ArithMatFloat<VVT,HVT,M1T,ROWS>& m,
                          ArithMatFloat<HVT,VVT,M2T,COLS>& m_new)
    {
        for(unsigned int i=0;i<M2T::get_v_dim();++i)
            for(unsigned int j=0;j<M2T::get_h_dim();++j)
                m_new[i][j] = m[j][i];
    }
    
#else
    
    //----------------- win32 -------------------------------
    // Visual studio is not good at deducing the args. to these template functions.
    // This means that you can call the two functions below with
    // matrices of wrong dimension.
    
    template <class M1, class M2, class M>
    constexpr void mul(const M1& m1, const M2& m2, M& m)
    {
        constexpr auto cols = M::get_h_dim();
        constexpr auto rows1 = M1::get_v_dim();
        constexpr auto rows2 = M2::get_v_dim();
        
        for(unsigned int i=0;i<rows1;++i)
            for(unsigned int j=0;j<cols;++j)
            {
                m[i][j] = 0;
                for(unsigned int k=0;k<rows2;++k)
                    m[i][j] += m1[i][k] * m2[k][j];
            }
    }
    
    
    /** Transpose. See the discussion on mul if you are curious as to why
     I don't simply return the transpose. */
    template <class M1, class M2>
    constexpr void transpose(const M1& m1, M2& m2)
    {
        for(unsigned int i=0;i<M2::get_v_dim();++i)
            for(unsigned int j=0;j<M2::get_h_dim();++j)
                m2[i][j] = m1[j][i];
    }
    
#endif
    
    /** Compute the outer product of a and b: a * transpose(b). This is
     a matrix with a::rows and b::columns. */
    template <class VVT, class HVT, class MT, unsigned int ROWS>
    constexpr void outer_product(const VVT& a, const HVT& b,
                       ArithMatFloat<VVT,HVT,MT,ROWS>& m)
    {
        constexpr auto R = VVT::get_dim();
        constexpr auto C = HVT::get_dim();
        for(unsigned int i=0;i<R;++i)
            for(unsigned int j=0;j<C;++j)
            {
                m[i][j] = a[i] * b[j];
            }
    }
    
    
    /** Compute the outer product of a and b using an arbitrary
     binary operation: op(a, transpose(b)). This is
     a matrix with a::rows and b::columns. */
    template <class VVT, class HVT, class MT, int ROWS, class BinOp>
    constexpr void outer_product(const VVT& a, const HVT& b,
                       ArithMatFloat<VVT,HVT,MT,ROWS>& m, BinOp op)
    {
        int R = VVT::get_dim();
        int C = HVT::get_dim();
        for(int i=0;i<R;++i)
            for(int j=0;j<C;++j)
            {
                m[i][j] = op(a[i], b[j]);
            }
    }
    
    /** Copy a matrix to another matrix, cell by cell.
     This conversion that takes a const matrix as first argument
     (source) and a non-const matrix as second argument
     (destination). The contents of the first matrix is simply copied
     to the second matrix. 
     
     However, if the first matrix is	larger than the second,
     the cells outside the range of the destination are simply not
     copied. If the destination is larger, the cells outside the 
     range of the source matrix are not touched.
     
     An obvious use of this function is to copy a 3x3 rotation matrix
     into a 4x4 transformation matrix.
     */
    template <class M1, class M2>
    constexpr void copy_matrix(const M1& inmat, M2& outmat)
    {
        const unsigned int R = std::min(inmat.get_v_dim(), outmat.get_v_dim());
        const unsigned int C = std::min(inmat.get_h_dim(), outmat.get_h_dim());
        for(unsigned int i=0;i<R;++i)
            for(unsigned int j=0;j<C;++j)
                outmat[i][j] = inmat[i][j];
    }

    /// @}
    /// @name stdio operations
    /// @related ArithMatFloat
    /// @{

    /// Put to operator
    template <class VVT, class HVT, class MT, unsigned int ROWS>
    std::ostream&
    operator<<(std::ostream&os, const ArithMatFloat<VVT,HVT,MT,ROWS>& m)
    {
        os << "[\n";
        for(unsigned int i=0;i<ROWS;i++) os << "  " << m[i] << "\n";
        os << "]\n";
        return os;
        }
        
    /// Get from operator
    template <class VVT, class HVT, class MT, unsigned int ROWS>
    std::istream& operator>>(std::istream&is,
                                    const ArithMatFloat<VVT,HVT,MT,ROWS>& m)
    {
        for(unsigned int i=0;i<ROWS;i++) is>>m[i];
        return is;
    }

    /// @}
}
#endif
