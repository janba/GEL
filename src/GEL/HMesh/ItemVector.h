/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */
/**
 * @file ItemVector.h
 * @brief Contains the vector type used for the mesh entities.
 */

#ifndef __HMESH_ITEMVECTOR_H__
#define __HMESH_ITEMVECTOR_H__

#include <cassert>
#include <vector>
#include "ItemID.h"

namespace HMesh
{
    /** The ItemVector is a vector of mesh entities.
        An ItemVector is a layer on top of the regular vector class which keeps
        track of the unused elements in a vector. This allows for garbage collection.
        ItemVector is used for storing vectors of Face, HalfEdge, and Vertex type.
     */
    template<typename ITEM>
    class ItemVector
    {
    public:
        typedef ItemID<ITEM> IDType;
        
        
        /// default constructor
        ItemVector(size_t _size = 0, ITEM i = ITEM()); 
        
        /// Get a reference to item i from kernel
        ITEM& get(IDType i);
        /// Get a const reference to item i from kernel
        const ITEM& get(IDType i) const;
        /// Get a reference to item i from kernel
        ITEM& operator[](IDType i);
        /// Get a const reference to item i from kernel
        const ITEM& operator[](IDType i) const;

        /// Add an entity to the kernel
        IDType add(const ITEM& i);

        /// remove an entity from kernel - entity is NOT erased!
        void remove(IDType i);

        /// erase unused entities from the kernel
        void cleanup();

        /// active size of vector
        size_t size() const;

        /// total size of vector
        size_t allocated_size() const;

        /// Clear the kernel
        void clear();

        /// Check if entity i is used
        bool in_use(IDType i) const;

        /// get starting index (default: skip to first active index)
        IDType index_begin(bool skip = true) const;

        /// get one past the end index of vector
        IDType index_end() const;

        /// get the next index (default: skip to first active index)
        IDType index_next(IDType index, bool skip = true) const;

        /// get the previous index (default: skip to first active index)
        IDType index_prev(IDType index, bool skip = true) const;

    private:

        size_t size_active;
        std::vector<ITEM> items;

        /// Memory consideration - objects flagged as unused should be remembered for future use (unless purged)
        std::vector<bool> active_items;
    };

    template<typename ITEM>
    inline ItemVector<ITEM>::ItemVector(size_t _size, ITEM i) 
        :   size_active(_size), 
            items(_size, i), 
            active_items(_size, true){}

    template<typename ITEM>
    inline ITEM& ItemVector<ITEM>::get(IDType id)
    {
        assert(id.index < items.size());
        return items[id.index];
    }

    template<typename ITEM>
    inline const ITEM& ItemVector<ITEM>::get(IDType id) const
    {
        assert(id.index < items.size());
        return items[id.index];
    }

    template<typename ITEM>
    inline ITEM& ItemVector<ITEM>::operator [](IDType id)
    {
        assert(id.index < items.size());
        return items[id.index];
    } 

    template<typename ITEM>
    inline const ITEM& ItemVector<ITEM>::operator [](IDType id) const
    {
        assert(id.index < items.size());
        return items[id.index];
    }

    template<typename ITEM>
    inline typename ItemVector<ITEM>::IDType ItemVector<ITEM>::add(const ITEM& item)
    {
        items.push_back(item);
        active_items.push_back(true);
        ++size_active;
        return IDType(items.size() - 1);
    }

    template<typename ITEM>
    inline void ItemVector<ITEM>::remove(typename ItemVector<ITEM>::IDType id)
    {
        if(active_items[id.index]){
            --size_active;
            active_items[id.index] = false;
        }
    }

    template<typename ITEM>
    inline void ItemVector<ITEM>::cleanup()
    {
        std::vector<ITEM> new_items;
        for(size_t i = 0; i < items.size(); ++i){
            if(active_items[i]) 
                new_items.push_back(items[i]);          
        }
        std::swap(items, new_items);
        active_items = std::vector<bool>(items.size(), true);
        size_active = items.size();
    }

    template<typename ITEM>
    inline size_t ItemVector<ITEM>::size() const
    { return size_active; }

    template<typename ITEM>
    inline size_t ItemVector<ITEM>::allocated_size() const
    { return items.size(); }

    template<typename ITEM>
    inline void ItemVector<ITEM>::clear()
    {
        items.clear();
        active_items.clear();
        size_active = 0;
    }

    template<typename ITEM>
    inline bool ItemVector<ITEM>::in_use(typename ItemVector<ITEM>::IDType id) const
    {
        if(id.index < active_items.size())
			return active_items[id.index];
		return false;
    }

    template<typename ITEM>
    inline typename ItemVector<ITEM>::IDType ItemVector<ITEM>::index_begin(bool skip) const
    {
        size_t i = 0;

        if(!skip)
            return IDType(i);
        
        while(i < active_items.size() && !active_items[i])
            ++i;


        return IDType(i);
    }

    template<typename ITEM>
    inline typename ItemVector<ITEM>::IDType ItemVector<ITEM>::index_end() const
    { return IDType(items.size()); }

    template<typename ITEM>
    inline typename ItemVector<ITEM>::IDType ItemVector<ITEM>::index_next(typename ItemVector<ITEM>::IDType id, bool skip) const
    {
        if(id.index < items.size())
            ++id.index;

        if(!skip)
            return IDType(id);

        while(id.index < items.size() && !active_items[id.index])
           ++id.index;

        return id;
    }

    template<typename ITEM>
    inline typename ItemVector<ITEM>::IDType ItemVector<ITEM>::index_prev(typename ItemVector<ITEM>::IDType id, bool skip) const
    {
        if(id.index > 0)
            --id.index;

        if(!skip)
            return id;

        while(!active_items[id.index] && id.index > 0)
            --id.index;

        return id;
    }

}

#endif 
