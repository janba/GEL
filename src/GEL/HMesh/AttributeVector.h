/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
  @file AttributeVector.h
 Contains an abstract class template for attribute vectors for HMesh entities. 
 Also contains templates for attribute vectors for concrete entities.
 */

#ifndef __HMESH_ATTRIBUTEVECTOR_H__
#define __HMESH_ATTRIBUTEVECTOR_H__

#include <cassert>
#include <vector>
#include <map>

namespace HMesh 
{
    /** Abstract class for HMesh entity attribute vectors. This class is used for storing all attributes
     associated with any type of mesh entity - be it Vertex, HalfEdge, or Face. Also the position attribute
     is stored in an AttributeVector class.      
     */
    template<typename ITEM, typename ITEMID>
    class AttributeVector
    {
    public:
        /// Construct from optional size and item (size should be identical to number of entities in associated container
        AttributeVector(size_t _size = 0, ITEM item = ITEM());

        /// Add an item
       // ITEMID add(const ITEM& item);

        /// const reference to item given by ID
        const ITEM& get(ITEMID id) const;

        /// reference to item given by ID
        ITEM& get(ITEMID id);

        /// const reference to item given by ID
        const ITEM& operator [](ITEMID id) const;

        /// reference to item given by ID
        ITEM& operator [](ITEMID id);

        /// resize the vector (may be necessary if associated container size grows)
        void resize(size_t _size, ITEM item = ITEM());

        /// number of attribute items in vector
        size_t size() const;

        /// clear the vector
        void clear();

        /// clenup unused items from the vector, given by remap from associated container
        void cleanup(const std::map<ITEMID, ITEMID>& map);

    private:
        std::vector<ITEM> items;
    };

    template<typename ITEM>
    class VertexAttributeVector : public AttributeVector<ITEM, VertexID>
    {
    public:
        VertexAttributeVector(size_t _size = 0, ITEM item = ITEM());
    };

    template<typename ITEM>
    class FaceAttributeVector : public AttributeVector<ITEM, FaceID>
    {
    public:
        FaceAttributeVector(size_t _size = 0, ITEM item = ITEM());
    };

    template<typename ITEM>
    class HalfEdgeAttributeVector : public AttributeVector<ITEM, HalfEdgeID>
    {
    public:
        HalfEdgeAttributeVector(size_t _size = 0, ITEM item = ITEM());
    };

    template<typename ITEM>
    inline VertexAttributeVector<ITEM>::VertexAttributeVector(size_t _size, ITEM item) : AttributeVector<ITEM, VertexID>(_size, item){}

    template<typename ITEM>
    inline FaceAttributeVector<ITEM>::FaceAttributeVector(size_t _size, ITEM item) : AttributeVector<ITEM, FaceID>(_size, item){}

    template<typename ITEM>
    inline HalfEdgeAttributeVector<ITEM>::HalfEdgeAttributeVector(size_t _size, ITEM item) : AttributeVector<ITEM, HalfEdgeID>(_size, item){}


    template<typename ITEM, typename ITEMID>
    inline AttributeVector<ITEM, ITEMID>::AttributeVector(size_t _size, ITEM item) : items(_size, item){}

    template<typename ITEM, typename ITEMID>
    inline void AttributeVector<ITEM, ITEMID>::clear()
    { items.clear(); }

    template<typename ITEM, typename ITEMID>
    inline void AttributeVector<ITEM, ITEMID>::cleanup(const std::map<ITEMID, ITEMID>& remap)
    {
        std::vector<ITEM> new_items(remap.size());
        for(typename std::map<ITEMID, ITEMID>::const_iterator it = remap.begin(); it != remap.end(); ++it){
            assert(it->second.index < remap.size());
            new_items[it->second.index] = items[it->first.index];
        }
        std::swap(items, new_items);
    }


    template<typename ITEM, typename ITEMID>
    inline void AttributeVector<ITEM, ITEMID>::resize(size_t _size, ITEM item)
    { items.resize(_size, item); }

    template<typename ITEM, typename ITEMID>
    inline size_t AttributeVector<ITEM, ITEMID>::size() const
    { return items.size(); }

    // just returning should be ok; manifold and attribs should always be in sync.
    // const context means manifold and attribs should be const, hence in sync.
    template<typename ITEM, typename ITEMID>
    inline const ITEM& AttributeVector<ITEM, ITEMID>::get(ITEMID id) const
    { 
        assert(id.index < items.size());
        return items[id.index]; 
    }

    template<typename ITEM, typename ITEMID>
    inline ITEM& AttributeVector<ITEM, ITEMID>::get(ITEMID id)
    {
        if(id.index >= items.size())
            items.resize(id.index + 1);
        return items[id.index]; 
    }

    template<typename ITEM, typename ITEMID>
    inline const ITEM& AttributeVector<ITEM, ITEMID>::operator [](ITEMID id) const
    { 
        assert(id.index < items.size());
        return items[id.index]; 
    }

    template<typename ITEM, typename ITEMID>
    inline ITEM& AttributeVector<ITEM, ITEMID>::operator [](ITEMID id)
    {
        if(id.index >= items.size())
            items.resize(id.index + 1);
        return items[id.index]; 
    }
}

#endif