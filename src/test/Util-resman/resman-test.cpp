#include <iostream>
#include <GEL/Util/ResourceManager.h>

using namespace std;
using namespace Util;

struct Blob
{
	double x;

	int step(int iter)
	{
		x += .5 / ( (iter-.75)*(iter-.25) );
		return 0;
	}

	~Blob() {cout << "Blob destroyed" << endl;}
};


double x=0;
int step(int iter)
{
	x += .5 / ( (iter-.75)*(iter-.25) );				
	return 0;
}

int (*step_fun)(int) = step;


ResourcePtr<Blob> g_ptr;
int main()
{
 	Blob b;
 	b.x = 312312;

	// Test 1: Register resource and verify that there is 1 resource
	register_static_resource<Blob>("myblob", &b);
	cout << "no resources (1) " << get_no_resources<Blob>() << endl;

	// Test 1.1: In debug mode, uncommenting the line below will make the 
	// program fail since a static resource cannot be removed by calling 
	// delete. In release mode, calling this function on a static resource
	// has no effect.
	// get_resource_ptr<Blob>("myblob").remove_when_unused();


	// Test 2:
	// Create a new resources and mark it as "remove when unused".
	// it goes away immediately since nothing points to it.
	cout << "Should print \"Blob destroyed\" below " << endl;
	register_dynamic_resource<Blob>("myblob2", new Blob).remove_when_unused();
	cout << "no resources (still 1) " << get_no_resources<Blob>() << endl;

	// Test 3:
	// Create a new resource and get a ResourcePointer. Mark as remove
	// when unused and verify that it goes away when ResourcePtr 
	// overwritten by a null resource pointer.
	ResourcePtr<Blob> b_ptr_x = register_dynamic_resource<Blob>("myblob3", new Blob);
	b_ptr_x.remove_when_unused();
	{
		ResourcePtr<Blob> gni;
		cout << "Should print \"Blob destroyed\" below " << endl;
		b_ptr_x = gni;
	}
	cout << "no resources (still 1) " << get_no_resources<Blob>() << endl;

	// Test 4:
	// Print value to ascertain the pointer is valid
	ResourcePtr<Blob> b_ptr = get_resource_ptr<Blob>("myblob");
	cout << "Value should be " << b.x << " : " << (b_ptr->x) << endl;

	// Test 5:
	// Print usage count
	cout << "Usage count (should be 1) " << b_ptr.usage() << endl;

	// Test 6:
	// Test that assigning to a new ResourcePtr increases usage count
	ResourcePtr<Blob> b_ptr2;
	b_ptr2 = get_resource_ptr<Blob>("myblob");
	cout << "Usage count (should be 2) : " << b_ptr.usage() << endl;


	// Test 7:
	// Create a local ResourcePtr
	{
		ResourcePtr<Blob> b_ptr3 = b_ptr;
		cout << "Usage count (should be 3) : " << b_ptr.usage() << endl;
	}
	
	// Test 8:
	// Test copy constructor
 	ResourcePtr<Blob> b_ptr4(b_ptr);
	cout << "Usage count (should be 3) : " << b_ptr.usage() << endl;

	// Test 9:
	// Get and throw away a resource ptr to test that usage is not increased.
	get_resource_ptr<Blob>("myblob");
	cout << "Usage count (should be 3) : " << b_ptr.usage() << endl;

	// Test 10:
	// Get raw pointer
	Blob* blob_raw_ptr = get_resource_ptr<Blob>("myblob").get_raw_ptr();
	cout << "Access with raw pointer " << blob_raw_ptr->x << endl;
	cout << "Usage count (should be 3) : " << b_ptr.usage() << endl;

	// --------------
	
	Blob b2;
	b2.x = 666;
 	register_static_resource<Blob>("my other blob", &b2);

	// Test 10:
	// Assign another resource to b_ptr4, check that usage count decreases
 	b_ptr4 = get_resource_ptr<Blob>("my other blob");
	cout << "Usage count (should be 2) : " << b_ptr.usage() << endl;
	cout << "Usage count (should be 1) : " << b_ptr4.usage() << endl;

	// Test 11:
	// Check that self assignment does not do anything
	b_ptr4 = b_ptr4;
	cout << "Usage count (should be 1) : " << b_ptr4.usage() << endl;

	// Test 12:
	// Assign the null resource ptr
 	ResourcePtr<Blob> gni;
	cout << "Usage count for 0 resource (should be -1) : " << gni.usage() 
			 << endl;
	
	// Test 13:
	// Test that overwriting a resource ptr with a new resource ptr to the
	// same resource does not increase usage count.
 	b_ptr4 = get_resource_ptr<Blob>("my other blob");
	cout << "Usage count (should be 1) : " << b_ptr4.usage() << endl;

	// Test 14:
	// Efficiency test. This test just computes an approximation of pi
	// using a function in a resource.
	b_ptr4->x = 0;
	for(int i=1;i<5000000;++i)
		{
			b_ptr4->step(i);
		}
	cout << " PI " << b_ptr4->x << endl;

// 	x = 0;
// 	for(int i=1;i<500000000;++i)
// 		{
// 			step(i);
// 		}
// 	cout << " PI " << x << endl;

	// Test 15: 
	// Global resource ptr test. If the relinquish line is commented out,
	// the program will give a warning upon termination and fail in debug
	// mode.
	g_ptr = b_ptr4;
	g_ptr.relinquish_resource();

	// At the end of the program, the two local blobs are destroyed.
	cout << "Two Blobs are now destroyed " << endl;
}

