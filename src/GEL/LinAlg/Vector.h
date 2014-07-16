/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/*!
	\file Vector.h
	\brief The Vector type
*/
#if !defined(VECTOR_H__HAA_AGUST_2001)
#define VECTOR_H__HAA_AGUST_2001

#include <cassert>
#include <iostream>
#include <cmath>
#include "../CGLA/Vec2f.h"
#include "../CGLA/Vec3f.h"

namespace LinAlg
{
	template <class T> class CMatrixType;  

	/*!
		\brief The Vector type.
	
		This vector teplate is one of the basic linear algebra types. In 
		principle it an overloaded Array of type T. The typedef Vector with 
		T set to double

		typedef CVectorType<double> CVector;
	
		is the delclaration intended to be used as the matrix type outside 
		the this and related files. The reson for making a template is 
		twofold, first it allows for an esay change of precision in this
		linear algebre package. Secondly it allows for other Vector types 
		in special cases.
	
		The pointer to the data area can be obtained by operator [] (int i) 
		with i=0, e.g. 
	
		CVector A;
	
		T* pData=A[0];  
	
		will  make pData point to the data.
	
		\see LapackFunc.h
		\see CMatrixType  
		\see Additional vector operators (located in this file Vector.h)
		\see CMatrix
		\see CVector
		\author Henrik Aanæs
		\version Aug 2001
	*/

	template <class T>
	class CVectorType  
	{

		friend class CMatrixType<T>;

	private:

		///The number of elements in the vector
		int nElems;

		///The pointer to the vector data
		T*Data;

		///Sets all the elements in the vector to a scalar
		void SetScalar(const T& Scalar) 
		{
			T* Itt=Data;
			T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					*Itt=Scalar;
					++Itt;
				}
		}; 

		///Frees the used memory
		void CleanUp(void)  {if(Data!=NULL) delete [] Data;};
		///Allocate the needed memory. Note that nElems should be initialized correctly first.

		void Allocate(void) { Data= new T[nElems];};


	public:

		/// \name Constructurs and Destructors
		//@{
		///Creats a vector a size 0.
		CVectorType<T>():nElems(0),Data(NULL) {};

		///Creates a vector of the specified Length
		explicit CVectorType<T>(const int Length): nElems(Length) {Allocate();}

		///Creates a vector of the specified Length, and sets all the values to a scalar.
		CVectorType<T>(const int Length,const T& Scalar): nElems(Length) {Allocate();SetScalar(Scalar);}

		///Copy constructor
		CVectorType<T>(const CVectorType<T>& Rhs): nElems(Rhs.nElems)
		{
			Allocate();
			memcpy(Data,Rhs.Data,nElems*sizeof(T));
		}

		///Copy constructor from Vec2Type
		template <class TT, class V, unsigned int N> 
		CVectorType<T>(const CGLA::ArithVec<TT,V,N>& Rhs):nElems(N)
		{
			Allocate(); 
			for(int i=0;i<N;++i)
				Data[i]=Rhs[i];
		}

		///
		virtual ~CVectorType<T>(){ CleanUp();};
		//@}


		/*!\name Asignment operators. 
			Note that the vector is resized if it does not have the correct size.
		*/
		//@{
		CVectorType<T>& operator=(const T& Rhs) {SetScalar(Rhs);return *this;}

		template <class TT, class V, unsigned int N> 
		CVectorType<T>& operator=(const CGLA::ArithVec<TT,V,N>& Rhs)
		{
			Resize(N);
			for(int i=0;i<N;++i)
				Data[i]=Rhs[i];
		}


		CVectorType<T>& operator=(const CVectorType<T>& Rhs)
		{
			Resize(Rhs.nElems);
			memcpy(Data,Rhs.Data,nElems*sizeof(T));
			return *this;
		}
		//@}

		/// \name Dimension functions
		//@{
		///Returns the length of the vector.
		const int Length(void) const {return nElems;}
		///Realocates memory IF the vector does not already hve the correct size.
		void Resize(const int Length)
		{
			if(Length!=nElems)
				{
					CleanUp();
					nElems=Length;
					Data= new T[nElems];
				}
		}
		//@}

		/*!
			\name Acces functions and operators.
			The [] operators are the ones intaended for usual use. The get 
			and set functions are added to allow for genral vector function 
			templates.
		*/
		//@{
		const T& get(const int i)const 
		{
			assert(i>=0);
			assert(i<nElems);

			return Data[i];
		}

		void set(const int i,const T& val)
		{
			assert(i>=0);
			assert(i<nElems);

			Data[i]=val;
		}

		T& operator[](const int i) 
		{
			assert(i>=0);
			assert(i<nElems);
		
			return Data[i];
		}

		const T& operator[](const int i) const 
		{
			assert(i>=0);
			assert(i<nElems);
		
			return Data[i];
		}
		//@}

		/// \name aritmethic operators
		//@{
		CVectorType<T> operator+(const CVectorType<T>& Rhs) const {CVectorType<T>Ret(*this); return Ret+=Rhs;};
		CVectorType<T> operator-(const CVectorType<T>& Rhs) const {CVectorType<T>Ret(*this); return Ret-=Rhs;};
		CVectorType<T> operator+(const T& Rhs) const {CVectorType<T>Ret(*this); return Ret+=Rhs;};
		CVectorType<T> operator-(const T& Rhs) const {CVectorType<T>Ret(*this); return Ret-=Rhs;};
		CVectorType<T> operator*(const T& Rhs) const {CVectorType<T>Ret(*this); return Ret*=Rhs;};
		CVectorType<T> operator/(const T& Rhs) const {CVectorType<T>Ret(*this); return Ret/=Rhs;};

		CVectorType<T>& operator+=(const CVectorType<T>& Rhs)
		{
			assert(nElems==Rhs.nElems);

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

		CVectorType<T>& operator-=(const CVectorType<T>& Rhs)
		{
			assert(nElems==Rhs.nElems);

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

		CVectorType<T>& operator+=(const T& Rhs)
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

		CVectorType<T>& operator-=(const T& Rhs)
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

		CVectorType<T>& operator*=(const T& Rhs)
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

		CVectorType<T>& operator/=(const T& Rhs)
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

		T operator*(const CVectorType<T>& Rhs) const
		{
			assert(nElems==Rhs.nElems);
			T Ret=0;

			T* IttRhs=Rhs.Data;
			T* Itt=Data;
			T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					Ret+=(*Itt)*(*IttRhs);
					++Itt;
					++IttRhs;
				}
			return Ret; 
		}

		CVectorType<T> operator*(const CMatrixType<T>& Rhs) const
		{
			assert(nElems==Rhs.Rows());
			CVectorType<T> Ret(Rhs.Cols(),0);

			T* ThisItt=Data;
			const T* RhsItt=Rhs.Data;
			const T* Stop=&(Rhs.Data[Rhs.nElems]);
			T* RetItt;

			while(RhsItt!=Stop)
				{
					RetItt=Ret.Data;
					for(int cCol=Rhs.Cols()-1;cCol>=0;--cCol)
						{
							(*RetItt)+=(*ThisItt)*(*RhsItt);
							++RetItt;
							++RhsItt;
						}
					++ThisItt;
				}
			return Ret;
		}
		//@}

		/// \name Norms
		//@{
		const T Norm(void) const {return sqrt((*this)*(*this));}
		const T NormMax(void) const
		{
			T Ret=0;
			T Abs;
			const T* Itt=Data;
			const T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					Abs=(*Itt)>0?(*Itt):-(*Itt);
					if(Ret<Abs)
						Ret=Abs;

					++Itt;
				}
			return Ret;
		}
		//@}

		/// Sum of all the elements in the vector
		const T Sum(void) const
		{
			T Ret=0;
			const T* Itt=Data;
			const T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					Ret+=(*Itt);
					++Itt;
				}
		}


		/*!
			\name Elementwise functions.
			Functions the perform some elementary operation on each 
			of the elements of the vector.
		*/
		//@{
		///The elements of the vector squared. Equvivalent to MatLAb's .^2
		void ElemSqr(void)
		{
			T* Itt=Data;
			T* Stop=&(Data[nElems]);
			while(Itt!=Stop)
				{
					(*Itt)*=(*Itt);
					++Itt;
				}
		}


	};


	/*!
		\name Additional vector operators
		These are operators heavily associated with CVectorType, but not included
		in the class definition it self.
	*/
	//@{
	template<class T>
	inline CVectorType<T> operator+(const T& Lhs,const CVectorType<T>&Rhs ) 
	{
		return Rhs+Lhs;
	}

	template<class T>
	inline CVectorType<T> operator*(const T& Lhs,const CVectorType<T>&Rhs ) 
	{
		return Rhs*Lhs;
	}

	template <class T>
	std::ostream& operator<<(std::ostream &s, const CVectorType<T> &a)
	{
		int nElems=a.Length();
	
		for (int cElem=0; cElem<nElems; cElem++)
			{
				s << a[cElem] << " ";
			}
		s << "\n";
	
		return s;
	}


	template <class T>
	std::istream& operator>>(std::istream &s, CVectorType<T> &a)
	{
		int nElems;

		s >> nElems;

		if ( nElems!=a.Length())
			a.Resize(nElems);


		for (int cElem=0; cElem<a.Length();cElem++)
			{
				s >>  a[cElem];
			}

		return s;
	}
	//@}

	/// The Vector annotation intended for use.
	typedef CVectorType<double> CVector;
}
#endif // !defined(VECTOR_H__HAA_AGUST_2001)
