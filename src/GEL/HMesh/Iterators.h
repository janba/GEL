/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Iterators.h
 * @brief Contains class for iterating over mesh entities in a HMesh.
 */
#ifndef HMESH_ITERATORS_H
#define HMESH_ITERATORS_H

#include <GEL/HMesh/ItemVector.h>
#include <iterator>
#include <ranges>

namespace HMesh {
    /** Traverse the entities of an HMesh in the order they are stored in the
      data structure. */
    template<typename ITEM> class IDIterator {
    public:
        typedef ItemVector<ITEM> vector_type;

        // typedefs to accommodate stl compliance
        typedef ptrdiff_t difference_type;
        typedef std::bidirectional_iterator_tag iterator_category;
        typedef ItemID<ITEM> value_type;
        typedef value_type reference;
        typedef value_type* pointer;

        // default constructor required by certain range concepts
        IDIterator() = default;
        /// constructor (default: skipping enabled)
        IDIterator(const vector_type& _item_vector, value_type _id, bool _skip = true);

        /// prefix increment
        IDIterator& operator++();
        /// postfix increment
        IDIterator operator++(int);
        /// prefix decrement
        IDIterator& operator--();
        /// postfix decrement
        IDIterator operator--(int);

        /// equal to
        bool operator==(const IDIterator& other) const;
        /// not equal to
        bool operator!=(const IDIterator& other) const;

        /// indirection
        value_type operator*() const;

    private:
        const vector_type* item_vector;
        value_type id;
        bool skip = true;
    };

    template<typename IterType> class IteratorPair {
        IterType f, l;

    public:
        IteratorPair(IterType _f, IterType _l): f(_f), l(_l) {}
        [[nodiscard]] IterType begin() const { return f; }
        [[nodiscard]] IterType end() const { return l; }
    };

    /*-----------------------------------------
     * IDIterator template implementation
     *-----------------------------------------*/

    template<typename ITEM>
    IDIterator<ITEM>::IDIterator(const vector_type& _item_vector, value_type _id, bool _skip)
        : item_vector(&_item_vector), id(_id), skip(_skip)
    {}

    template<typename ITEM>
    IDIterator<ITEM> IDIterator<ITEM>::operator++(int)
    {
        auto tmp(*this);
        ++(*this);
        return tmp;
    }

    template<typename ITEM>
    IDIterator<ITEM> IDIterator<ITEM>::operator--(int)
    {
        auto tmp(*this);
        --(*this);
        return tmp;
    }

    template<typename ITEM>
    bool IDIterator<ITEM>::operator==(const IDIterator<ITEM>& other) const
    {
        return item_vector == other.item_vector && id == other.id;
    }

    template<typename ITEM>
    bool IDIterator<ITEM>::operator!=(const IDIterator<ITEM>& other) const
    {
        return item_vector != other.item_vector || id != other.id;
    }

    template<typename ITEM>
    typename IDIterator<ITEM>::value_type IDIterator<ITEM>::operator*() const
    {
        return id;
    }

    template<typename ITEM>
    IDIterator<ITEM>& IDIterator<ITEM>::operator++()
    {
        id = item_vector->index_next(id, skip);
        return *this;
    }

    template<typename ITEM>
    IDIterator<ITEM>& IDIterator<ITEM>::operator--()
    {
        id = item_vector->index_prev(id, skip);
        return *this;
    }

    static_assert(std::bidirectional_iterator<IDIterator<int>>);
    static_assert(std::ranges::viewable_range<const IteratorPair<IDIterator<int>>&>);
} // namespace HMesh

template<typename T> constexpr bool ::std::ranges::enable_borrowed_range<HMesh::IteratorPair<T>> = true;

#endif
