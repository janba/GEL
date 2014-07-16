/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file HashTable.h
 * @brief Hash table class.
 */

#ifndef __UTIL_HASHTABLE_H
#define __UTIL_HASHTABLE_H

#include "HashKey.h"

#include <functional>
#include <iostream>

namespace Util
{

	namespace
	{
		const int HASH_EXPAND_LOAD = 80;
	};

	/** \brief Hashtable class template. 

	This class template has taken for too long to get right. That is because
	it is flexible. It can be used either as a simple (key, value) database
	or as a reference counted database.

	There are two template arguments: KEY_T is the hash key type. This type
	must provide operators == and != and a hash function which takes an 
	integer denoting the table size and returns another integer - the table
	argument.
		
	VAL_T is the value type. There are few restrictions on this type but it 
	must have default and copy constructors and operator=. The HashTable 
	template is intended for use with simple types.
	*/
	template<class KEY_T, class VAL_T>
	class HashTable
	{
		/** struct containing hashtable entries.
				Apart from the value, this struct also contains a reference
				count and the hash key */
		struct EntryType
		{
			// Hash key
			KEY_T key;

			/// Table value corresponding to key
			VAL_T val;

			/// Reference counter
			signed char ref;
		
		public:
			EntryType(): ref(-1) {}
		};
	
		/// The table itself (array of entries)
		EntryType* table;   

	private:

		/// Number of elements in the table
		int items;          
  
		/// Size of table
		int size;       

		/// Resize the table (new size is argument)
		void resize(int new_size);

		/** Private find entry function. Returns true iff entry was found 
				and reference non-zero. the second argument is set to the table
				index corresponding to key. */
		bool find_entry_priv(const KEY_T&, int&) const;
  
	public:

		/// Construct an empty hash table
		HashTable(): items(0), size(1) 
		{
			table = new EntryType[size];
		}

		/// Destruct it
		~HashTable() { delete [] table; }

		/** Create entry in hash table.
				The first argument is the key, and the second is the value
				to be inserted. The third argument is the amount that the 
				reference count should be increased.
			
				If third argument is zero: The key,value pair is inserted
				if key is not found. In either case the reference count is 
				set to 1. function returns true if key was not found (or 
				reference count=0 which counts as key begin absent)

				If third argument is non-zero: Behaviour is as above except 
				that the reference count is incremented or set to value of 
				third argument depending on whether the key was already present. */
		bool create_entry(const KEY_T&, const VAL_T&,int=0);

		bool create_entry(const KEY_T& kv, VAL_T*& val, int incr=0);



		/** Removes an entry from database.
				This function returns false if the key was not found. Otherwise
				the reference count is set to zero and the function returns true
				having first reduced the size of the table, if the size was 
				below threshold after removal. */
		bool remove_entry(const KEY_T&);

		/** Find entry in database.
				The second argument is set to the value that corresponds to the key
				passed as first argument. The function returns true if the key
				was found (and reference non-zero) */
		bool find_entry(const KEY_T&, VAL_T&) const;

		/** Decrease reference count. THe key is first argument, the value
				is returned in second argument, and the third argument contains 
				the amount of decrease, 1 is default. The function returns true
				if the key was found and reference non-zero. False otherwise.

				if the reference count is zero after being decreased, the function
				may decrease the size of the table. 
		*/
		bool decr_reference(const KEY_T&, VAL_T&, int=1);


		/** Change Entry. Does not affect reference count. Returns true
				if the entry exists - false otherwise. */
		bool change_entry(const KEY_T&, const VAL_T&);

		/// Print various statistics for the table 
		void print() const;

		/** Get the first entry in HashTable. All arguments are passed by 
				reference. The first argument is the index. Upon return the 
				first argument will contain the first valid index in the table.
				This index should be passed to the get_next function. (see below).
				Subsequent arguments are the key, the value and the reference 
				count. These are set to the values that correspond to the index. 
				The function returns false if there are no valid indices - i.e. 
				the table is empty 
		*/
		bool get_first(int& idx, KEY_T&,VAL_T& val,int&) const;

		/** Get the next entry in the HashTable. All arguments are passed
				by reference.
				This first argument is a valid index to an element in the table 
				(from a previous call to either this function or get_first)
				get_next finds the next valid entry and returns key, value and 
				reference count in the next three arguments. get_next returns 
				false if there is no next valid index. */
		bool get_next(int& idx, KEY_T&,VAL_T& val,int&) const;

		/// Returns the sum of the reference counts that are > 0
		int integrate() const; 


		/// Returns no of items in table.
		int no_items() const {return items;}

		/// Returns true if number of items==0
		bool empty() const {return items==0;}

		/// Returns total size of table (in bytes)
		int total_size() const {return sizeof(*this) + sizeof(EntryType)*size;}

		/// Returns the full size of an entry (key + value + ref.count).
		static int entry_size() {return sizeof(EntryType);}
	
	};


	template<class KEY_T, class VAL_T>
	bool HashTable<KEY_T, VAL_T>::find_entry_priv(const KEY_T& kv, int& k) const
	{
		k = kv.hash(size);
		int fk = k;
		while(table[k].ref != -1)
			{
				if (table[k].key == kv) 
					return (table[k].ref > 0) ? true : false;
				k = (k+1)&(size-1);
				if (k==fk) break;
			}
		return false;
	}

	template<class KEY_T, class VAL_T>
	bool HashTable<KEY_T, VAL_T>::find_entry(const KEY_T& kv, VAL_T& val) const
	{
		int k;
		if(find_entry_priv(kv,k))
			{
				val = table[k].val;
				return true;
			}
		return false;
	}

	template<class KEY_T, class VAL_T>
	void HashTable<KEY_T, VAL_T>::resize(int new_size)
	{ 
		int j=0;
		EntryType* tmp_table = new EntryType[new_size];
		for(int i=0; i<size;i++)
			if (table[i].ref >0)
				{
					j++;
					int k = table[i].key.hash(new_size);
				
					while(tmp_table[k].ref != -1) 
						k=(k+1)&(new_size-1); 
				
					tmp_table[k].key = table[i].key;
					tmp_table[k].val = table[i].val;
					tmp_table[k].ref = table[i].ref;
				}
      
		delete table;
		table = tmp_table;
		size = new_size;
	}

	template<class KEY_T, class VAL_T>
	bool HashTable<KEY_T, VAL_T>::create_entry(const KEY_T& kv, 
																						 const VAL_T& val,
																						 int incr)
	{
		if (items*100 >= HASH_EXPAND_LOAD * size) 
			resize(size*2);

		int k;
		if(!find_entry_priv(kv,k))
			{
				++items;
				table[k].ref = (incr==0? 1: incr);
				table[k].key = kv;
				table[k].val = val;
				return true;
			}
		table[k].ref += incr;
		table[k].val = val;
		return false;
	}

	template<class KEY_T, class VAL_T>
	bool HashTable<KEY_T, VAL_T>::create_entry(const KEY_T& kv, 
																						 VAL_T*& val,
																						 int incr)
	{
		if (items*100 >= HASH_EXPAND_LOAD * size) 
			resize(size*2);

		int k;
		if(!find_entry_priv(kv,k))
			{
				++items;
				table[k].ref = (incr==0? 1: incr);
				table[k].key = kv;
				val = &table[k].val;
				return true;
			}
		table[k].ref += incr;
		val=&table[k].val;
		return false;
	}


	template<class KEY_T, class VAL_T>
	bool HashTable<KEY_T, VAL_T>::remove_entry(const KEY_T& kv)
	{
		int k;
		if(find_entry_priv(kv,k))
			{
				table[k].ref = 0;
				--items;
				if(items*100 < (HASH_EXPAND_LOAD * size) / 2)
					resize(size/2);
				return true;
			}
		return false;
	}



	template<class KEY_T, class VAL_T>
	bool HashTable<KEY_T, VAL_T>::decr_reference(const KEY_T& kv,  
																							 VAL_T& val, 
																							 int decr)
	{
		int k;
		if(find_entry_priv(kv,k))
			{
				table[k].ref -= decr;
				val = table[k].val;
				if(table[k].ref == 0)
					{
						--items;
						if(items*100 < (HASH_EXPAND_LOAD * size) / 2)
							resize(size/2);
					}
				return true;
			}
		return false;
	}

	template<class KEY_T, class VAL_T>
	bool HashTable<KEY_T, VAL_T>::change_entry(const KEY_T& kv,  
																						 const VAL_T& val)
	{
		int k;
		if(find_entry_priv(kv,k))
			{
				table[k].val = val;
				return true;
			}
		return false;
	}


	template<class KEY_T, class VAL_T>
	void HashTable<KEY_T, VAL_T>::print() const
	{
		std::cout << "Entry size " << sizeof(EntryType) << " bytes\n"
							<< "Max number of items " << size	<< "\n" 
							<< "Number of items " << items << std::endl;
	}

	template<class KEY_T, class VAL_T>
	int HashTable<KEY_T, VAL_T>::integrate() const 
	{
		int integral =0;
		for(int i=0;i<size; i++)
			if(table[i].ref>0) integral += table[i].ref;
		return integral;
	}

	template<class KEY_T, class VAL_T>
	bool HashTable<KEY_T, VAL_T>::get_first(int& idx, 
																					KEY_T& key,
																					VAL_T& val,
																					int& r) const 
	{
		if(size==0) return false;
	
		idx=0;
		while(table[idx].ref <= 0)
			{
				++idx;
				if(idx==size) return false;
			}
		key = table[idx].key;
		val = table[idx].val;
		r   = table[idx].ref;
		return true;
	}

	template<class KEY_T, class VAL_T>
	bool HashTable<KEY_T, VAL_T>::get_next(int& idx, 
																				 KEY_T& key,
																				 VAL_T& val,
																				 int& r) const
	{
		idx++;
		if(idx==size) return false;
		while(table[idx].ref <= 0)
			{
				++idx;
				if(idx==size) return false;
			}
		key = table[idx].key;
		val = table[idx].val;
		r   = table[idx].ref;
		return true;
	}

}


#endif
