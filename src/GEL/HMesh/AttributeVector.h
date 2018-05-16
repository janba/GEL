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
        /// Construct from  size and optional default value (size should be identical to number of entities.
        AttributeVector(size_t _size, ITEM _default_value) :
        items(_size, _default_value),
        default_value(_default_value) {}

        /// Construct from optional default value.
        AttributeVector(ITEM _default_value):
        default_value(_default_value) {}

        /// const reference to item given by ID
        const ITEM& get(ITEMID id) const     {
            assert(id.index < items.size());
            return items[id.index];
        }

        /// reference to item given by ID
        ITEM& get(ITEMID id)     {
            if(id.index >= items.size())
                items.resize(id.index + 1, default_value);
            return items[id.index];
        }


        /// const reference to item given by ID
        const ITEM& operator [](ITEMID id) const {
            return get(id);
        }


        /// reference to item given by ID
        ITEM& operator [](ITEMID id) {
            return get(id);
        }

        /// resize the vector (may be necessary if associated container size grows)
        void resize(size_t _size, ITEM _default_value = ITEM()) {
            default_value = _default_value;
            items.resize(_size, default_value);
        }


        /// number of attribute items in vector
        size_t size() const {
            return items.size();
        }


        /// clear the vector
        void clear() {
            items.clear();
        }

        /// cleanup unused items from the vector, given by remap from associated container
        void cleanup(const std::map<ITEMID, ITEMID>& remap) {
            std::vector<ITEM> new_items(remap.size());
            for(const auto& mapping : remap){
                assert(mapping.second.index < remap.size());
                if(mapping.first.index < items.size())
                    new_items[mapping.second.index] = items[mapping.first.index];
            }
            std::swap(items, new_items);
        }

    private:
        std::vector<ITEM> items;
        ITEM default_value;
    };

    template<typename ITEM>
    class VertexAttributeVector : public AttributeVector<ITEM, VertexID>
    {
    public:
        VertexAttributeVector(size_t _size, ITEM _default_value):
        AttributeVector<ITEM, VertexID>(_size, _default_value){}

        VertexAttributeVector(ITEM _default_value = ITEM()):
        AttributeVector<ITEM, VertexID>(_default_value){}
    };

    template<typename ITEM>
    class FaceAttributeVector : public AttributeVector<ITEM, FaceID>
    {
    public:
        FaceAttributeVector(size_t _size, ITEM _default_value):
        AttributeVector<ITEM, FaceID>(_size, _default_value){}

        FaceAttributeVector(ITEM _default_value = ITEM()):
        AttributeVector<ITEM, FaceID>(_default_value){}


    };

    template<typename ITEM>
    class HalfEdgeAttributeVector : public AttributeVector<ITEM, HalfEdgeID>
    {
    public:
        HalfEdgeAttributeVector(size_t _size, ITEM _default_value):
        AttributeVector<ITEM, HalfEdgeID>(_size, _default_value){}

        HalfEdgeAttributeVector(ITEM _default_value = ITEM()):
        AttributeVector<ITEM, HalfEdgeID>(_default_value){}

    };
}

#endif
