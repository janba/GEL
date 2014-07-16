/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Iterators.h
 * @brief Contains class for iterating over mesh entities in a HMesh.
 */
#ifndef __HMESH_ITERATORS_H__
#define __HMESH_ITERATORS_H__

#include "ItemVector.h"

namespace HMesh
{
    /** Traverse the entities of an HMesh in the order they are stored in the
      data structure. */
    template<typename ITEM>
    class IDIterator
    {
    public:
        typedef ItemID<ITEM> ID;
        
        // typedefs to accommodiate stl compliance
        typedef std::bidirectional_iterator_tag iterator_category;
        typedef ptrdiff_t difference_type;
        typedef ID value_type;
        typedef value_type reference;
        typedef value_type* pointer;

        /// constructor (default: skipping enabled)
        IDIterator(const ItemVector<ITEM>& _item_vector, ID _id, bool _skip = true);

        /// prefix increment 
        IDIterator& operator ++();		
        /// postfix increment
        IDIterator& operator ++(int);
        /// prefix decrement
        IDIterator& operator --();
        /// postfix decrement
        IDIterator& operator --(int);

        /// equal to
        bool operator ==(const IDIterator& other) const;
        /// not equal to
        bool operator !=(const IDIterator& other) const;

        /// indirection
        reference operator *();
        /// member by pointer
        pointer operator ->();
        /// cast
        //operator VertexID() const;

    private:
        const ItemVector<ITEM>* item_vector;
        ID id;
        bool skip;
    };
    
    template<typename ITEM>
    class IDIteratorPair
    {
    public:
        IDIteratorPair(IDIterator<ITEM> _f, IDIterator<ITEM> _l): f(_f), l(_l) {}
        IDIterator<ITEM> begin() {return f;}
        IDIterator<ITEM> end() {return l;}
    private:
        IDIterator<ITEM> f,l;
    };
    


    /*-----------------------------------------
     * IDIterator template implementation
     *-----------------------------------------*/

    template<typename ITEM>
    inline IDIterator<ITEM>::IDIterator(const ItemVector<ITEM>& _item_vector, ID _id, bool _skip) 
        : item_vector(&_item_vector), id(_id), skip(_skip){}

    template<typename ITEM>
    inline IDIterator<ITEM>& IDIterator<ITEM>::operator ++(int)
    { return ++(*this); }

    template<typename ITEM>
    inline IDIterator<ITEM>& IDIterator<ITEM>::operator --(int)
    { return --(*this); }

    template<typename ITEM>
    inline bool IDIterator<ITEM>::operator ==(const IDIterator<ITEM>& other) const
    { return item_vector == other.item_vector && id == other.id; }

    template<typename ITEM>
    inline bool IDIterator<ITEM>::operator !=(const IDIterator<ITEM>& other) const
    { return item_vector != other.item_vector || id != other.id; }

    template<typename ITEM>
    inline ItemID<ITEM> IDIterator<ITEM>::operator *()
    { return id; }

    template<typename ITEM>
    inline ItemID<ITEM>* IDIterator<ITEM>::operator ->()
    { return &id; }

    template<typename ITEM>
    inline IDIterator<ITEM>& IDIterator<ITEM>::operator ++()
    {
        id = item_vector->index_next(id, skip);
        return *this;
    }

    template<typename ITEM>
    inline IDIterator<ITEM>& IDIterator<ITEM>::operator --()
    {
        id = item_vector->index_prev(id, skip);
        return *this;
    }
}

#endif