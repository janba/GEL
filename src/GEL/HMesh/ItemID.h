/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file ItemID.h
 * @brief Base class for the integer IDs used to refer to mesh entities.
 */

#ifndef HMESH_ITEMID_H
#define HMESH_ITEMID_H

#include <iostream>

namespace HMesh {
    /** The ItemID class is simply a wrapper around an index. This class associates a type
     with the index. */
    template<typename T>
    class ItemID {
    public:
        using IndexType = size_t;
        using EntityType = T;

        constexpr ItemID(): index(INVALID_INDEX) {}

        bool operator==(const ItemID& other) const { return index == other.index; }
        bool operator!=(const ItemID& other) const { return index != other.index; }
        bool operator<(const ItemID& other) const { return index < other.index; }

        IndexType get_index() const { return index; }
        template<typename ITEM, typename ITEMID> friend class AttributeVector;
        //    private:
        IndexType index;
        static constexpr IndexType INVALID_INDEX = -1;

        explicit constexpr ItemID(IndexType _index): index(_index) {}

        friend class ConnectivityKernel;

        template<typename ITEM> friend class ItemVector;

        template<typename X> friend std::ostream& operator<<(std::ostream& os, const ItemID<X>&);
    };
    static_assert(sizeof(ItemID<size_t>) == sizeof(size_t));

    template<typename T> std::ostream& operator<<(std::ostream& os, const ItemID<T>& iid)
    {
        return (os << iid.index);
    }

} // namespace HMesh

// Injecting into std::hash is explicitly allowed
template<typename T> struct std::hash<HMesh::ItemID<T>> {
    std::size_t operator()(const HMesh::ItemID<T>& item_id) const noexcept
    {
        constexpr std::hash<size_t> hasher;
        return hasher(item_id.index);
    }
};

#endif
