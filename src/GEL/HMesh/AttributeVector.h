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
#include <GEL/Util/Serialization.h>

namespace HMesh 
{
    /** Abstract class for HMesh entity attribute vectors. This class is used for storing all attributes
     associated with any type of mesh entity - be it Vertex, HalfEdge, or Face. Also the position attribute
     is stored in an AttributeVector class.      
     */
    template<typename ITEM, typename ITEMID>
    class AttributeVector
    {
        std::vector<ITEM> items;
        ITEM default_value;
        
    public:
        /// Construct from  size and optional default value (size should be identical to number of entities.
        AttributeVector(size_t _size, ITEM _default_value) :
        items(_size, _default_value),
        default_value(_default_value) {}

        /// Construct from optional default value.
        AttributeVector(ITEM _default_value):
        default_value(_default_value) {}
        
        /// Return raw pointer to data
        ITEM* data() {return items.data();}
        
        /// Return raw const pointer to data
        const ITEM* data() const {return items.data();}

        /// Returns  value corresponding to index.
        ITEM operator [](ITEMID idx) const {
            return items.at(idx.index);
        }

        /// reference to item given by ID. Note that this function may resize the vector.
        ITEM& operator [](ITEMID idx) {
            if (idx.index >= items.size())
                resize(idx.index+1);
            return items.at(idx.index);
        }

        /** resize the vector (may be necessary if associated container size grows). This
         version of the function takes a new default value as argument. This will override and
         overwrite the previous default value. Rarely useful since the default value can be set
         in the constructor. */
        void resize(size_t _size, ITEM _default_value) {
            default_value = _default_value;
            items.resize(_size, default_value);
        }

        /// resize the vector (may be necessary if associated container size grows)
        void resize(size_t _size) {
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
        
        void serialize(Util::Serialization& ser) const {
            ser.write(items);
            ser.write(default_value);
        }
        
        void deserialize(Util::Serialization& ser) {
            ser.read(items);
            ser.read(default_value);
        }

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
