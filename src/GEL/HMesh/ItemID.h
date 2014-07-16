/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file ItemID.h
 * @brief Base class for the integer IDs used to refer to mesh entities.
 */


#ifndef __HMESH_ITEMID_H__
#define __HMESH_ITEMID_H__

#include <iostream>

namespace HMesh
{
    /** The ItemID class is simply a wrapper around an index. This class associates a type
     with the index. */
    template<typename T>
    class ItemID
    {
    public:
        typedef size_t IndexType;

        ItemID(): index(INVALID_INDEX){}

        bool operator ==(const ItemID& other) const { return index == other.index; }
        bool operator !=(const ItemID& other) const { return index != other.index; }
        bool operator <(const ItemID& other) const { return index < other.index; }
        
        IndexType get_index() const { return index;}
		
    private:
        IndexType index;
        static const IndexType INVALID_INDEX =  -1;

        explicit ItemID(IndexType _index): index(_index){}

        friend class ConnectivityKernel;
            
        template<typename ITEM> friend class ItemVector;
        template<typename ITEM, typename ITEMID> friend class AttributeVector;
        template<typename X>
        friend std::ostream& operator<<(std::ostream& os, const ItemID<X>&);

    };

	template<typename T>
	inline std::ostream& operator<<(std::ostream& os, const ItemID<T>& iid)
	{
		return (os << iid.index);
	}
	
}

#endif