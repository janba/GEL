/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file ResourceManager.h
 * @brief Resource manager that is not used much in GEL but potentially useful for larger projects.
 */

#ifndef __UTIL_RESOURCE_MANAGER_H
#define __UTIL_RESOURCE_MANAGER_H

#include <cassert>
#include <string>
#include <list>
#include <typeinfo>

#define CLEAN_SHUTDOWN 1

namespace Util
{

	typedef unsigned char FlagByte;
	const FlagByte REMOVE_WHEN_UNUSED = 0x01;
	const FlagByte STATIC_RESOURCE    = 0x02;

	/** \brief This class template represents a resource record. 

	    There is a pointer to 
			the actual resource (of type RES) and a usage counter which should never
			be less than 0. */
	template<class RES>
	class ResourceRecord
	{
		/// Name of the resource (const)
    const std::string name;

		/// A pointer to the resource (the pointer not the pointee is const)
		RES * const res;
		
		/// Usage counter.
		int usage;

		/// Flags
		FlagByte flags;

	public:
	
		/// Construct a null record.
		ResourceRecord(): res(0), usage(0), flags(0) {}

		/// Construct a resource record with a name and a pointer.
		ResourceRecord(const std::string& _name, RES* _res, bool static_res=false): 
			name(_name), 
			res(_res), 
			usage(0), 
			flags(static_res ? STATIC_RESOURCE : 0) 
		{
			assert(res != 0);
		}

		~ResourceRecord() 
		{
#if CLEAN_SHUTDOWN
			assert(usage==0);		
#endif
		}
	
		void erase_resource()
		{
			assert(usage==0);
			if(!(flags&STATIC_RESOURCE)) 
				{
					delete res;
				}
		}

		/// Increment the usage counter
		void increment_usage() 
		{
			assert(usage>=0);
			++usage;
		}

		/// Decrement the usage counter. assert that counter is >0
		void decrement_usage() 
		{
			assert(usage>0);
			--usage;
		}

		/// Return the usage count (mostly for debugging)
		int get_usage() const { return usage; }

		/// Get the name of a resource
		const std::string& get_name() const {return name;}

		/// Get a pointer to the resource (const pointer)
		RES * const get_ptr() 
		{
			return res;
		}

		/// Get a resource pointer (const pointer and pointee)
		const RES * const get_ptr() const 
		{
			return res;
		}

		void remove_when_unused() 
		{
/* 			assert(!(flags&STATIC_RESOURCE)); */
/* 			if(!(flags&STATIC_RESOURCE)) */
				flags = flags|REMOVE_WHEN_UNUSED;
		}

		FlagByte get_flags() const { return flags;}
	
	};


	template<class RES> class ResourceManager;


	/** \brief Template for a pointer to a reference counted resource.

      The ResourcePtr class template is a template for a reference
			counted resource pointer. It is a smart pointer that actually points
			to a ResourceRecord and not to the actual resource. Since the record
			contains a reference count, we can increase the reference count when
			the pointer is copied and decrease it in the destructor. Only the
			ResourceManager can create ResourcePtr's directly from a raw
			ResourceRecord. */
	template<class RES>
	class ResourcePtr
	{
		typedef ResourceRecord<RES> RR;
		typedef std::list<RR> RRList;
		typedef typename RRList::iterator RRListIter;

		static RRListIter get_null_rrlist_iter()
		{
			static RRList l;
			return l.end();
		}

#define NULL_RRLIST_ITER ResourcePtr<RES>::get_null_rrlist_iter() 		

		friend class ResourceManager<RES>;

	private:		
		RRListIter rr;
		RES* res;
	
		/** Explicit constructor with value. The constructor can only be called
				by a friend, i.e. the resource manager. */
		explicit ResourcePtr(const RRListIter& _rr): rr(_rr), res(0)
		{
			if(rr != NULL_RRLIST_ITER)
				{
					rr->increment_usage();
					res = rr->get_ptr();
				}
		}

		void decrement_usage();
	
	public:
	
		ResourcePtr(): rr(NULL_RRLIST_ITER), res(0) {}
	
		ResourcePtr(const ResourcePtr& r2): rr(r2.rr), res(0)
		{
			if(rr != NULL_RRLIST_ITER)
				{
					rr->increment_usage();
					res = rr->get_ptr();
				}
		}

		const ResourcePtr& operator=(const ResourcePtr& r2) 
		{
			// Guard against self-assignment
			if (r2.rr != this->rr)
				{
					if(rr != NULL_RRLIST_ITER)	decrement_usage();

					rr  = r2.rr;
					res = 0;

					if(rr != NULL_RRLIST_ITER)
						{
							res = rr->get_ptr();
							rr->increment_usage();
						}
				}
			return *this;
		}

		~ResourcePtr() {decrement_usage();}

		RES& operator*() const 
		{
			assert(rr != NULL_RRLIST_ITER);
			assert(res != 0);
			return *res;
		}

		RES* const operator->() const 
		{
			assert(rr != NULL_RRLIST_ITER);
			assert(res !=0);
			return res;
		}

		RES* const get_raw_ptr() const 
		{
			assert(rr != NULL_RRLIST_ITER);
			assert(res != 0);
			return res;
		}

		int usage() const
		{
			if(rr != NULL_RRLIST_ITER)
				return rr->get_usage();
			return -1;
		}

		bool is_valid() const {return rr != NULL_RRLIST_ITER;}

		void relinquish_resource()
		{
			ResourcePtr p;
			*this = p;
		}

		/** Calling this function sets the REMOVE_WHEN_UNUSED flag. If this 
				flag is set, the resource record is removed from the manager */
		void remove_when_unused() {	rr->remove_when_unused(); }
	};


	/** \brief Resource manager class.

      The ResourceManager is a singleton, and it uses Scott Meyers construct
			on first use idiom, i.e. the static function get_instance() will
			construct it when first called. There may be a problem with this
			idion if resources depend on each other. However, that is only
			when the program shuts down.  See modern C++ design by
			Alexandrescu for more information.

			I'll try to make everything in this class private and accessed
			by a simple interface of friend functions.

	*/

	template<class RES>
	class ResourceManager
	{
		typedef ResourceRecord<RES> RR;
		typedef std::list<RR> RRList;
		typedef typename RRList::iterator RRListIter;
		RRList resources;
	
		ResourceManager(): resources(0) {}
		ResourceManager(const ResourceManager&);
		ResourceManager& operator=(const ResourceManager&);

	public:


		~ResourceManager()
		{
#if CLEAN_SHUTDOWN
			RRListIter i = resources.begin(); 
			while(i != resources.end())
				{
					if(i->get_usage()==0)
						{
							RRListIter tmp = i;
							++i;
							erase_resource(tmp);
						}
					else
						{
							std::cout << "Warning, ResourceManager:\n\n"
												<< (typeid(this).name())
												<< "\n\nis shutting down, and resource \n\n" 
												<< i->get_name() << "\n\nhas usage: " << i->get_usage()
												<< std::endl;
							std::cout <<
 								"In other words, this resource is not unused at this\n"
								"point. That is unfortunate because then it cannot\n"
								"be deleted (since it might be used by some other\n"
								"part of the program during shutdown). Please ensure\n"
								"that all resources are unused at program\n"
								"termination. If you use a global ResourcePtr, this is\n"
								"done by calling relinquish_resource just before\n"
								"exiting.\n"
												<< std::endl;
							++i;
						}
				}
#endif
		}

		int get_no_resources() const
		{
			return resources.size();
		}

		static ResourceManager& get_instance() 
		{
			static ResourceManager instance;
			return instance;
		}

		void erase_resource(RRListIter iter)
		{
			iter->erase_resource();
			resources.erase(iter);
		}

		/** Find a resource and return a ResourcePtr to it. Returns the 
				0 ResourcePtr if no string found. */
		ResourcePtr<RES> get_resource_ptr(const std::string& str)
		{
			for(RRListIter i = resources.begin(); i != resources.end(); ++i)
				if((*i).get_name() == str) 
					return ResourcePtr<RES>(i);
			return ResourcePtr<RES>(NULL_RRLIST_ITER);
		}

		/** Add a resource with name passed as first argument and a pointer
				passed as 2nd argument. A ResourcePtr to the just inserted
				element is returned. */
		ResourcePtr<RES> register_resource(const std::string& str, 
																			 RES* res, bool static_resource)
		{
			for(RRListIter i = resources.begin(); i != resources.end(); ++i)
				if(i->get_name() == str)
					return (ResourcePtr<RES>(NULL_RRLIST_ITER));

			resources.push_front(RR(str, res, static_resource));
			ResourcePtr<RES> ptr = ResourcePtr<RES>(resources.begin());
			return ptr;
		}

/*		friend class ResourcePtr<RES>;
		friend int get_no_resources<RES>();
		friend ResourcePtr<RES> get_resource_ptr<RES>(const std::string& str); 
		friend ResourcePtr<RES> register_static_resource<RES>(const std::string&,
																													RES*); 
		friend ResourcePtr<RES> register_dynamic_resource<RES>(const std::string&,
																													 RES*);*/ 
	};

	template<class RES>
	inline void ResourcePtr<RES>::decrement_usage()
	{
		assert( (rr == NULL_RRLIST_ITER) || (res != 0));
		if(rr != NULL_RRLIST_ITER)
			{
				rr->decrement_usage();
				if(rr->get_usage() == 0 && (rr->get_flags()&REMOVE_WHEN_UNUSED))
					{
						ResourceManager<RES>& man = ResourceManager<RES>::get_instance();
						man.erase_resource(rr);
					}
			}
	}

	template<class RES>
	inline int get_no_resources()
	{
		ResourceManager<RES>& man = ResourceManager<RES>::get_instance();
		return man.get_no_resources();
	}

	template<class RES>
	inline ResourcePtr<RES> get_resource_ptr(const std::string& str)
	{
		ResourceManager<RES>& man = ResourceManager<RES>::get_instance();
		return man.get_resource_ptr(str);
	}

	template<class RES>
	inline ResourcePtr<RES> register_dynamic_resource(const std::string& str,RES* res)
	{
		ResourceManager<RES>& man = ResourceManager<RES>::get_instance();
		ResourcePtr<RES> ptr = man.register_resource(str, res, false);
		return ptr;
	}

	template<class RES>
	inline ResourcePtr<RES> register_static_resource(const std::string& str,RES* res)
	{
		ResourceManager<RES>& man = ResourceManager<RES>::get_instance();
		ResourcePtr<RES> ptr = man.register_resource(str, res, true);
		return ptr;
	}

}


#endif
