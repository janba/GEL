/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "LinAlgIO.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

namespace
{

		bool Comp(istreambuf_iterator<char>& Itt,const string str)
		{
		
				size_t cChar=0;
				bool Ret=true;
				while(cChar<str.size() && Ret)
				{
						Ret=(*Itt==str[cChar]);
						++Itt;
						cChar++;
				}
		
				return Ret;
		}
	
		bool FindVar(istreambuf_iterator<char>& Itt,
								 istreambuf_iterator<char> End,
								 const string VarName)
		{
				bool bFound=false;
		
				while(!bFound)
				{
						Itt=find(Itt,End,VarName[0]);
						if(Itt==End)
								return false;
						bFound=Comp(Itt,VarName);
				}
		
				return true;
		}
	

}

namespace LinAlg
{
		void ToMatlab(const CMatrix& A,
									const std::string& VarName,
									const std::string& FileName,
									const bool append,
									const std::string& Comment)
		{
				std::ofstream Mfile(FileName.c_str(), 
														append ? std::ios::app : std::ios::trunc);
																				 
	
	
				if(!Comment.empty())
						Mfile << "% " << Comment.c_str() << "\n";
				Mfile << VarName.c_str() << "=[";
	

				const double* pRow;

				for(int cRow=0;cRow<A.Rows();cRow++)
				{
						pRow=A[cRow];
						for(int cCol=0;cCol<A.Cols();cCol++)
						{
								Mfile << "\t" << pRow[cCol];
						}
						Mfile << ";\n";
				}
	
				Mfile << "];\n\n";
		}


		void FromMatlab(CMatrix& A,
										const std::string& VarName,
										const std::string& FileName)
		{
				ifstream file(FileName.c_str());
				istreambuf_iterator<char> Itt(file);
				istreambuf_iterator<char> End;

				string VarId(VarName);
				VarId+="=[";

				FindVar(Itt,End,VarId);

				vector<double> data;
				double elem;
				int nCols=0;

				while((*Itt)!=';')	
				{
						Itt=(file >> elem );
						data.push_back(elem);
						nCols++;
				}

				Itt++;
				Itt++;
				int nElems=nCols;

				while((*Itt)!=']')
				{
						while((*Itt)!=';')	
						{
								Itt=(file >> elem );
								data.push_back(elem);
								nElems++;
						}
						Itt++;
						Itt++;
				}

				A.Resize(nElems/nCols,nCols);

				memcpy(A[0],&(data[0]),sizeof(double)*data.size());


/*	for(int cI=0;cI<20;cI++)
		{
		cout << *Itt;
		++Itt;
		}
*/
/*	for(int cI=0;cI<data.size();cI++)
		cout << data[cI] << "\t";

		cout << endl << nCols << endl << nElems << endl;
*/
	
		}

}
