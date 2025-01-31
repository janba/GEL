#include "GEL/HMesh/RsR.h"

int main() {
	read_config("C:/Users/ruicu/Desktop/Gel_integration/test/example_config.txt");
	HMesh::Manifold output;
	reconstruct_single(output);
	return 0;
}