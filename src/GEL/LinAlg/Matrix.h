/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#if !defined(MATRIX_H_HAA_AGUST_2001)
#define MATRIX_H_HAA_AGUST_2001

#include <cassert>
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "../CGLA/ArithMatFloat.h"

namespace LinAlg
{

	/*!
		\file Matrix.h
		\brief The Matrix type
	*/


	/*!
		\brief The Matrix type.

		This matrix teplate is one of the basic linear algebra types. In 
		principle it an overloaded Array of type T. The typedef Matrix with 
		T set to double

		typedef CMatrixType<double> CMatrix;

		is the delclaration intended to be used as the matrix type outside 
		the this and related files. The reson for making a template is 
		twofold, first it allows for an esay change of precision in this
		linear algebre package. Secondly it allows for other MAtrix types 
		in special cases.

		The data structure consists of the array of type T, T* Data, 
		with the Ilife vector, T** Ilife, refering to the begining 
		of each of the columns. Hereby a coherent data area is avalable and 
		the concept of 2-dimensionality is acheived by the Ilife vector. The 
		pointer to the data area can be obtained by operator [] (int i) 
		with i=0, e.g. 

		Matrix A;

		T* pData=A[0];  

		will  make pData point to the data.

		It should be noted that much of the functionality associated with this
		Matrix type and the associated Vector type is located in the 
		LapackFunc.h file.

		\see LapackFunc.h
		\see CVectorType  
		\see Additional matrix operators (Located in this file Matrix.h).
		\see CMatrix
		\see CVector
		\see CVec2Type
		\see CVec3Type
		\see CVec2
		\see CVec3

		\author Henrik Aanæs
		\version Aug 2001
	*/
	template <class T>
	class CMatrixType 
	{
		friend class CVectorType<T>;

	private:
		///The number of rows in the matrix
		int nRows;
		///The number of columns in the matrix
		int nCols;
		/// The number of elements in the matrix, i.e. nElems=nRows*nCols.
		int nElems;
		///The pointer to the data containing the elements of the matrix.
		T* Data;
		///The Ilife vector containing pointers to the start of the columns. That is pointers to elements in Data.
		T** Ilife;
	
		///Sets all the elements in the matrix to Scalar.
		void SetScalar(const T& Scalar) 
		{
			T* Itt=Data;
			T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					*Itt=Scalar;
					++Itt;
				}
		} 
	
		///Clears the allocated memory.
		void CleanUp(void)
		{
			if(Data!=NULL) 
				delete [] Data; 
			if(Ilife!=NULL) 
				delete [] Ilife;
		};
	
		///Allocates memory. Note that nRows,nCols and nElems should be set first.
		void Allocate(void)
		{ 
			Data= new T[nElems]; 
			Ilife=new T*[nRows];	
			for(int cRow=0;cRow<nRows;cRow++) 
				Ilife[cRow]=&Data[cRow*nCols];
		};
	
	public:
		///	\name The Constructors and Destructors of the CMatrixType class.
		//@{
		CMatrixType<T>() : nRows(0),nCols(0),nElems(0),Data(NULL),Ilife(NULL) {};

		CMatrixType<T>(const int Row,const int Col) : nRows(Row),nCols(Col),nElems(Row*Col) {Allocate();}

		CMatrixType<T>(const int Row,const int Col,const T& Scalar) : nRows(Row),nCols(Col),nElems(Row*Col) {Allocate(); SetScalar(Scalar);};

		CMatrixType<T>(const CMatrixType<T>& Rhs) : nRows(Rhs.nRows),nCols(Rhs.nCols),nElems(Rhs.nRows*Rhs.nCols) 
		{ 
			Allocate();
			memcpy(Data,Rhs.Data,nElems*sizeof(T));
		}

		template <class VVT, class HVT, class MT, unsigned int ROWS>
		CMatrixType<T>(const CGLA::ArithMatFloat<VVT,HVT,MT,ROWS>& Rhs): 
			nRows(Rhs.get_v_dim()),nCols(Rhs.get_h_dim()),
			nElems(nRows*nCols) 
		{ 
			Allocate();
			for(int i=0;i<nElems;++i)
				Data[i] = Rhs.get()[i];
		}

		CMatrixType<T>(const CVectorType<T>& Rhs) : nRows(Rhs.Length()),nCols(1),nElems(Rhs.Length()) 
		{ 
			Allocate();
			memcpy(Data,Rhs.Data,nElems*sizeof(T));
		}
		virtual ~CMatrixType<T>(){CleanUp();};
		//@}

		/*!
			\name	Asignment opertors.
			The assignment operators of the CMatrixType class.
			Note that the matrix is resized if it does not have the correct size.
		*/
		//@{
		CMatrixType<T>& operator=(const T& Rhs){ SetScalar(Rhs);return *this;}

    CMatrixType<T>& operator=(const CMatrixType<T>& Rhs)
		{ 
			Resize(Rhs.nRows,Rhs.nCols); 
			memcpy(Data,Rhs.Data,nElems*sizeof(T));
			return *this;
		}

		template <class VVT, class HVT, class MT, unsigned int ROWS>
		CMatrixType<T>& operator=(const CGLA::ArithMatFloat<VVT,HVT,MT,ROWS>& Rhs)
		{ 
			Resize(Rhs.get_v_dim(),Rhs.get_h_dim()); 
			for(int i=0;i<nElems;++i)
				Data[i] = Rhs.get()[i];
			return *this;
		}

		CMatrixType<T>& operator=(const CVectorType<T>&Rhs)
		{
			Resize(Rhs.Lenght(),1);
			memcpy(Data,Rhs.Dat,nElems*sizeof(T));
			return *this;
		}
		//@}


		/*!
			\name Dimension functions.
			Functions dealing with the dimension of the matrix.
		*/
		//@{
		/// Returns the number of rows in the matrix
		const int Rows(void) const {return nRows;};
		///Retruns the number of columns in the matrix. 
		const int Cols(void) const {return nCols;};
		/// Resizes the matrix IF it does not already have the desired dimensions.
		void Resize(const int Row,const int Col)
		{
			if(Row!=nRows || Col!=nCols)
				{
					CleanUp();
					nRows=Row;
					nCols=Col;
					nElems=Row*Col;
					Data= new T[nElems];
					Ilife=new T*[nRows];
					for(int cRow=0;cRow<nRows;cRow++)
						{
							Ilife[cRow]=&Data[cRow*nCols];
						}
				}
		}
		//@}
	
		/*!
			\name Acces functions and operators.
			The [] operators are the ones intaended for usual use. The get 
			and set functions are added to allow for genral Matrix function 
			templates.
		*/
		//@{
		T* operator[](const int Row) 
		{
			assert(Row>=0);
			assert(Row<nRows);
		
			return Ilife[Row];
		}
	
		const T* operator[](const int Row) const 
		{
			assert(Row>=0);
			assert(Row<nRows);
		
			return Ilife[Row];
		}

		const T& get(const int Row,const int Col) const
		{
			assert(Row>=0 && Row<nRows);
			assert(Col>=0 && Col<nCols);
			return Ilife[Row][Col];
		}

		void set(const int Row,const int Col,const T val) 
		{
			assert(Row>=0 && Row<nRows);
			assert(Col>=0 && Col<nCols);
			Ilife[Row][Col]=val;
		}
		//@}
	
		/*! 
			\name Arithmetic operators.
			The aritmic operators of the CMatrixType class.
		*/
		//@{
		CMatrixType<T> operator+(const T& Rhs) const {CMatrixType<T>Ret(*this); return Ret+=Rhs;};
		CMatrixType<T> operator-(const T& Rhs) const {CMatrixType<T>Ret(*this); return Ret-=Rhs;};
		CMatrixType<T> operator*(const T& Rhs) const {CMatrixType<T>Ret(*this); return Ret*=Rhs;};
		CMatrixType<T> operator/(const T& Rhs) const {CMatrixType<T>Ret(*this); return Ret/=Rhs;};
		CMatrixType<T>& operator+=(const T& Rhs)
		{

			T* Itt=Data;
			T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					(*Itt)+=Rhs;
					++Itt;
				}
		
			return *this; 
		}
	
		CMatrixType<T>& operator-=(const T& Rhs)
		{
			T* Itt=Data;
			T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					(*Itt)-=Rhs;
					++Itt;
				}
		
			return *this; 
		}
	
		CMatrixType<T>& operator*=(const T& Rhs)
		{
			T* Itt=Data;
			T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					(*Itt)*=Rhs;
					++Itt;
				}
		
			return *this; 
		}
	
		CMatrixType<T>& operator/=(const T& Rhs)
		{
			T* Itt=Data;
			T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					(*Itt)/=Rhs;
					++Itt;
				}
		
			return *this; 
		}

		CMatrixType<T> operator+(const CMatrixType<T>& Rhs) const {CMatrixType<T>Ret(*this); return Ret+=Rhs;};
		CMatrixType<T> operator-(const CMatrixType<T>& Rhs) const {CMatrixType<T>Ret(*this); return Ret-=Rhs;};

		CMatrixType<T>& operator+=(const CMatrixType<T>& Rhs)
		{
			assert(nRows==Rhs.nRows);
			assert(nCols==Rhs.nCols);

			T* IttRhs=Rhs.Data;
			T* Itt=Data;
			T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					(*Itt)+=(*IttRhs);
					++Itt;
					++IttRhs;
				}

			return *this; 
		}
	
		CMatrixType<T>& operator-=(const CMatrixType<T>& Rhs)
		{
			assert(nRows==Rhs.nRows);
			assert(nCols==Rhs.nCols);
			T* IttRhs=Rhs.Data;
			T* Itt=Data;
			T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					(*Itt)-=(*IttRhs);
					++Itt;
					++IttRhs;
				}
			return *this; 
		}

		CMatrixType<T> operator* (const CMatrixType<T>& Rhs) const
		{
			assert(nCols==Rhs.nRows);
			CMatrixType<T> Ret(nRows,Rhs.nCols);

			for(int cRow=0;cRow<nRows;cRow++)
				{
					for(int cCol=0;cCol<Rhs.nCols;cCol++)
						{
							T* Ret_ij=&(Ret[cRow][cCol]);
							*Ret_ij=0;
							for(int cI=0;cI<nCols;cI++)
								{
									(*Ret_ij)+=(Ilife[cRow])[cI]*Rhs[cI][cCol];
								}
						}
				}
			return Ret;
		}

		CVectorType<T> operator* (const CVectorType<T>& Rhs) const 
		{
			assert(nCols==Rhs.nElems);
			CVectorType<T> Ret(nRows,0);

			T* ThisItt=Data;
			T* Stop=&(Data[nElems]);
			T* RhsItt;
			T* RetItt=Ret.Data;
			while(ThisItt!=Stop)
				{
					RhsItt=Rhs.Data;
					for(int cCol=nCols;cCol>0;--cCol)
						{
							(*RetItt)+=(*RhsItt)*(*ThisItt);
							++ThisItt;
							++RhsItt;
						}
					++RetItt;
				}
			return Ret;
		}
		//@}
	
		/*!
			\name Matrix Transpose.
			Functions to transpose the matrix
		*/
		//@{
		///Transposes the matrix itself.
		CMatrixType<T>& Transpose(void)
		{
			T* NewData =new T[nElems];
			T*RowP;
			int cRow;
			for(cRow=0;cRow<nRows;cRow++)
				{
					RowP=Ilife[cRow];
					for(int cCol=0;cCol<nCols;cCol++)
						{
							NewData[cCol*nRows+cRow]=RowP[cCol];
						}
				}
		
			CleanUp();
		
			Data=NewData;
			int temp=nRows;
			nRows=nCols;
			nCols=temp;
			Ilife=new T*[nRows];	
			for(cRow=0;cRow<nRows;cRow++) 
				{
					Ilife[cRow]=&Data[cRow*nCols];
				}
		
			return *this;
		}

		///Returns the Transpose of the matrix is returned in AT.
		void Transposed(CMatrixType<T>& AT) const
		{
			AT.Resize(nCols,nRows);
			for(int cRow=nRows-1;cRow>=0;--cRow)
				{
					for(int cCol=nCols-1;cCol>=0;--cCol)
						{
							AT[cCol][cRow]=Ilife[cRow][cCol];
						}
				}
		}

		///Returns the Transpose of the matrix is returned.
		CMatrixType<T> Transposed(void) const
		{
			CMatrixType<T> Ret(nCols,nRows);
			Transposed(Ret);
			return Ret;
		}
		//@}


		/*!
			\name Elementwise functions.
			Functions the perform some elementary operation on each 
			of the elements of the matrix.
		*/
		//@{
		///Pairwise multipling the elemtns of the two matrices. Equvivalent to MatLAb's .*
		void ElemMult(const CMatrixType<T>& Rhs)
		{
			T*Itt=Data;
			T*Stop=&(Data[nElems]);
			T*RhsItt=Rhs.Data;
			while(Itt!=Stop)
				{
					(*Itt)*=(*RhsItt);
					++Itt;
					++RhsItt;
				}
		}

		///Pairwise dividnig the the elemtns of the two matrices. Equvivalent to MatLAb's ./
		void ElemDiv(const CMatrixType<T>& Rhs)
		{
			T*Itt=Data;
			T*Stop=&(Data[nElems]);
			T*RhsItt=Rhs.Data;
			while(Itt!=Stop)
				{
					(*Itt)/=(*RhsItt);
					++Itt;
					++RhsItt;
				}
		}

		///The elements of the matrix squared. Equvivalent to MatLAb's .^2
		void ElemSqr(void)
		{
			T*Itt=Data;
			T*Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					(*Itt)*=(*Itt);
					++Itt;
				}
		}

		///The square root of the elements. Equvivalent to MatLAb's .^0.5
		void ElemSqrt(void)
		{
			T*Itt=Data;
			T*Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					(*Itt)=sqrt(*Itt);
					++Itt;
				}
		}
		//@}

		/*!
			\name Norms.
			Computes various norms of the matrix.
		*/
		//@{
		///The Max norm returns the maximum absolute value
		const T NormMax(void) const
		{
			T Ret=0;
			T Abs;
			const T*Itt=Data;
			const T*Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					Abs=(*Itt)>0?(*Itt):-(*Itt);
					if(Ret<Abs)
						Ret=Abs;

					++Itt;
				}
			return Ret;
		}

		/// The Frobenius norm returns the square root of the element sum of squares.
		const T NormFrobenius(void) const
		{
			T Ret=0;
			const T*Itt=Data;
			const T*Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					Ret+=(*Itt)*(*Itt);
					++Itt;
				}
			return sqrt(Ret);
		}
		//@}


		/*!
			\name Sub-part acces functions
			These functions are used to acces sub parts of the matrix.
		*/
		//@{
		///Get the cRow'th row and set it equal to Row. Equvivalent to Row=this(cRow,:) in Matlab.
		void GetRow(CVectorType<T>& Row,const int cRow) const
		{
			Row.Resize(nCols);
			memcpy(Row.Data,Ilife[cRow],nCols*sizeof(T));
		}

		///Sets the cRow'th row of the matrix equal to Row. Equvivalent to this(cRow,:)=Row in Matlab.
		void SetRow(CVectorType<T>& Row,const int cRow)
		{
			assert(Row.Length()==nCols);
			memcpy(Ilife[cRow],Row.Data,nCols*sizeof(T));
		}

		///Copies row 'from' of this matrix to row 'to' of matrix A. Equvivalent to A(to,:)=this(from,:) in Matlab.
		void RowTransfere(CMatrixType<T>& A, const int from,const int to) const
		{
			assert(A.Cols()==nCols);
			memcpy(A.Ilife[to],Ilife[from],nCols*sizeof(T));
		}

		void GetCol(CVectorType<T>&Col, const int nCol) const
		{
			Col.Resize(nRows);
			const T* Itt=&Data[nCol];
			T* VecItt=&Col[0];
			for(int cRow=0;cRow<nRows;cRow++)
				{
					*VecItt=*Itt;
					++VecItt;
					Itt+=nCols;
				}
		}

		void SetCol(CVectorType<T>&Col, const int nCol)
		{
			assert(nRows==Col.Length());
			assert(nCol<nCols);

			T* Itt=&Data[nCol];
			T* VecItt=&Col[0];
			for(int cRow=0;cRow<nRows;cRow++)
				{
					*Itt=*VecItt;
					++VecItt;
					Itt+=nCols;
				}
		}

		void GetSubMatrix(CMatrixType<T>& A, const int StartRow, const int StartCol, const int RowDim, const int ColDim) const
		{
			assert(StartRow+RowDim<=nRows);
			assert(StartCol+ColDim<=nCols);
			A.Resize(RowDim,ColDim);
			const T* Itt=&(Ilife[StartRow][StartCol]);

			for(int cRow=0;cRow<A.Rows();cRow++)
				{
					memcpy(A[cRow],Itt,A.Cols()*sizeof(T));
					Itt+=nCols;
				}
		}

		void SetSubMatrix(CMatrixType<T>& A, const int StartRow, const int StartCol)
		{
			assert(A.Rows()+StartRow<=nRows);
			assert(A.Cols()+StartCol<=nCols);
			T* Itt=&(Ilife[StartRow][StartCol]);

			for(int cRow=0;cRow<A.Rows();cRow++)
				{
					memcpy(Itt,A[cRow],A.Cols()*sizeof(T));
					Itt+=nCols;
				}
		}

		//@}



	};

	/*!
		\name Additional matrix operators
		These are operators heavily associated with CMAtrixType, but not included
		in the class definition it self.
	*/
	//@{
	template<class T>
	inline CMatrixType<T> operator+(const T& Lhs,const CMatrixType<T>&Rhs ) 
	{
		return Rhs+Lhs;
	}

	template<class T>
	inline CMatrixType<T> operator*(const T& Lhs,const CMatrixType<T>&Rhs ) 
	{
		return Rhs*Lhs;
	}

	template <class T>
	std::ostream& operator<<(std::ostream &s, const CMatrixType<T> &A)
	{
    int nRows=A.Rows();
    int nCols=A.Cols();
	
    for (int cRow=0; cRow<nRows; cRow++)
			{
        for (int cCol=0; cCol<nCols; cCol++)
					{
            s << A[cRow][cCol] << " ";
					}
        s << "\n";
			}
    return s;
	}

	template <class T>
	std::istream& operator>>(std::istream &s, CMatrixType<T> &A)
	{
    int nRows;
		int nCols;

    s >> nRows >> nCols;

    if ( nRows!=A.Rows() || nCols!=A.Cols() )
			A.Resize(nRows,nCols);

		T* RowP;


    for (int cRow=0; cRow<nRows;cRow++)
			{
				RowP=A[cRow];
        for (int cCol=0;cCol<nCols;cCol++)
					{
            s >>  RowP[cCol];
					}
			}

    return s;
	}
	//@}

	/// The Matrix annotation intended for use.
	typedef CMatrixType<double> CMatrix;
}
#endif // !defined(MATRIX_H_HAA_AGUST_2001)

