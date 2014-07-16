/**
	 Small test program for KDTree class. The idea is to insert a lot of points
	 and find them again. Next a range search is used to test the find point
	 inside sphere function. 

	 Send bug reports or improvements to jab@imm.dtu.dk
*/

#include <cstdlib>
#include <GEL/Geometry/KDTree.h>
#include <GEL/CGLA/Vec3f.h>

using namespace GEO;
using namespace std;
using namespace CGLA;

void make_ran_point(Vec3f& p0)
{
	p0=Vec3f(10.0f*gel_rand()/GEL_RAND_MAX,
					 10.0f*gel_rand()/GEL_RAND_MAX,
					 10.0f*gel_rand()/GEL_RAND_MAX);
}

int main()
{
	cout << "\n\nTest 1: Insert and find " << endl;
	int sizes[15] = {1,2,3,4,7,21,134,256,589,678,751,892,1023,1024,1025};
	for(int i=0;i<15;++i)
		{
			int N=sizes[i];
			cout << "Tree size " << N << endl;
			gel_srand(0);	
			cout << "Insert and find points " << endl;
			KDTree<Vec3f,int> tree;
			for(int i=0;i<N;++i)
				{
					Vec3f p0;
					make_ran_point(p0);
					tree.insert(p0, i);
				}
			tree.build();
			gel_srand(0);
			for(int i=0;i<N;++i)
				{
					Vec3f p0;
					make_ran_point(p0);
					Vec3f pos;
					float d = 0.0000001f;
					int x;
					if(!tree.closest_point(p0, d, pos, x))
						{
							cout << " test failed " << endl;
							exit(1);
						}
					assert(x==i);
				}
		}

	cout << "Test passed " << endl;

	cout << "\n\nTest 2: Points in range " << endl;
	KDTree<Vec3f,int> tree;
	int K = 312;
	for(int i=0;i<K;++i)
		{
			Vec3f p0;
			make_ran_point(p0);
			tree.insert(p0, i);
		}
	
	cout << "Building tree of size = " << K <<endl;
	tree.build();

	Vec3f p0;
	make_ran_point(p0);
	std::vector<Vec3f> keys;
	std::vector<int> vals;

	float range = 1.95f;
	int N = tree.in_sphere(p0, range , keys, vals);

	cout << "Points in range (" << range << ") of " << p0 << " : " << N << endl;
	for(int i=0;i<N;++i)
		{
			cout << keys[i] << (keys[i]-p0).length() << " - " << vals[i] << endl;
		}
}
