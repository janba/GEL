//
// Created by Cem Akarsubasi on 8/7/25.
//

/// @file AssociativeContainers.h
/// @brief Type definitions for sets and maps

#ifndef GEL_ASSOCIATIVECONTAINERS_H
#define GEL_ASSOCIATIVECONTAINERS_H

#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>

/// FIXME: Include these when adding phmap
//#include <parallel_hashmap/btree.h>
//#include <parallel_hashmap/phmap.h>

/// FIXME: This is a facade that is meant to be removed when we include phmap as a dependency
/// @private
namespace phmap
{
template <typename Key>
using Hash = std::hash<Key>;

template <typename T>
using EqualTo = std::equal_to<T>;

template <typename T>
using Less = std::less<T>;
}

#include <GEL/Util/Assert.h>

namespace Util
{
/// @name Associative Containers
/// This contains associative containers that are meant to be used from the rest of GEL. It is currently a facade
/// as we want to defer adding this as a dependency and just move over some algorithms to use these
/// @{

template <typename Key, typename Value, typename Hash = phmap::Hash<Key>, typename Eq = phmap::EqualTo<Key>>
using HashMap = std::unordered_map<Key, Value, Hash, Eq>;
//phmap::flat_hash_map<Key, Value, Hash, Eq>;

template <typename Key, typename Value, typename Compare = phmap::Less<Key>>
    requires std::is_nothrow_copy_constructible_v<Compare>
using BTreeMap = std::map<Key, Value, Compare>;
//phmap::btree_map<Key, Value, Compare>;

template <typename Key, typename Value, typename Compare = phmap::Less<Key>>
    requires std::is_nothrow_copy_constructible_v<Compare>
using BTreeMultiMap = std::multimap<Key, Value, Compare>;
//phmap::btree_multimap<Key, Value, Compare>;

template <typename Key, typename Hash = phmap::Hash<Key>, typename Eq = phmap::EqualTo<Key>>
using HashSet = std::unordered_set<Key, Hash, Eq>;
//phmap::flat_hash_set<Key, Hash, Eq>;

template <typename Key, typename Compare = phmap::Less<Key>>
    requires std::is_nothrow_copy_constructible_v<Compare>
using BTreeSet = std::set<Key, Compare>;
//phmap::btree_set<Key, Compare>;

template <typename Key, typename Compare = phmap::Less<Key>>
    requires std::is_nothrow_copy_constructible_v<Compare>
using BTreeMultiSet = std::multiset<Key, Compare>;
//phmap::btree_multiset<Key, Compare>;

//using phmap::erase_if;

/// @}

// TODO: Might consider creating a header for common concepts used in GEL
namespace Concepts
{
    /// Checks if a range of type Range yields elements of type T
    template <typename Range, typename T>
    concept RangeYields = std::ranges::range<Range> && std::is_same_v<
        std::decay_t<std::ranges::range_value_t<Range>>, T>;
}
}

#endif //GEL_ASSOCIATIVECONTAINERS_H
