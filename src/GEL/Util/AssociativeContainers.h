//
// Created by Cem Akarsubasi on 8/7/25.
//

/// @file AssociativeContainers.h
/// @brief Type definitions for sets and maps

#ifndef GEL_ASSOCIATIVECONTAINERS_H
#define GEL_ASSOCIATIVECONTAINERS_H

#include <parallel_hashmap/btree.h>
#include <parallel_hashmap/phmap.h>

#include <GEL/Util/Assert.h>

namespace Util
{
template <typename Key, typename Value, typename Hash = phmap::Hash<Key>, typename Eq = phmap::EqualTo<Key>>
using HashMap = phmap::flat_hash_map<Key, Value, Hash, Eq>;

template <typename Key, typename Value, typename Compare = phmap::Less<Key>>
    requires std::is_nothrow_copy_constructible_v<Compare>
using BTreeMap = phmap::btree_map<Key, Value, Compare>;

template <typename Key, typename Value, typename Compare = phmap::Less<Key>>
    requires std::is_nothrow_copy_constructible_v<Compare>
using BTreeMultiMap = phmap::btree_multimap<Key, Value, Compare>;

template <typename Key, typename Hash = phmap::Hash<Key>, typename Eq = phmap::EqualTo<Key>>
using HashSet = phmap::flat_hash_set<Key, Hash, Eq>;

template <typename Key, typename Compare = phmap::Less<Key>>
    requires std::is_nothrow_copy_constructible_v<Compare>
using BTreeSet = phmap::btree_set<Key, Compare>;

template <typename Key, typename Compare = phmap::Less<Key>>
    requires std::is_nothrow_copy_constructible_v<Compare>
using BTreeMultiSet = phmap::btree_multiset<Key, Compare>;

using phmap::erase_if;

// TODO: Might consider creating a header for common concepts used in GEL
namespace Concepts
{
    /// Checks if a range of type Range yields elements of type T
    template<typename Range, typename T>
    concept RangeYields = std::ranges::range<Range> && std::is_same_v<std::decay_t<std::ranges::range_value_t<Range>>, T>;
}


/// A priority queue sorted by value. Sorted in ascending order by default.
///
/// @details This class generalizes a priority queue by allowing key ordering relations to be stored indirectly
/// This means that keys can be removed arbitrarily and not just from the front. Existing keys can be updated by
/// reinserting them to the queue.
///
/// Note that the memory footprint of this is substantially larger than an ordinary vector-based min/max heap. In
/// addition, the internal queue is always kept sorted, meaning there is a potentially higher up-front cost to
/// constructing this. So make sure you really need the mutability or really benefit from it when using this over
/// a generational priority queue.
template <typename Key,
          typename Value,
          typename Hash = phmap::Hash<Key>,
          typename Eq = phmap::EqualTo<Key>,
          typename Compare = phmap::Less<Value>>
class MutablePriorityQueue {
    Util::HashMap<Key, Value, Hash, Eq> m_distance_map;

    struct CompareFunctor {
        const MutablePriorityQueue& self;

        bool operator()(const Key& key1, const Key& key2) const
        {
            Compare cmp;
            return cmp(self.m_distance_map.at(key1), self.m_distance_map.at(key2));
        }
    };

    Util::BTreeSet<Key, CompareFunctor> m_queue;

public:
    using value_type = std::pair<Key, Value>;
    using size_type = decltype(m_queue)::size_type;

    MutablePriorityQueue() : m_queue(CompareFunctor{*this}) {};

    /// Remove the first value in the queue
    value_type pop_front()
    {
        GEL_ASSERT(!m_queue.empty());
        const auto key = *m_queue.begin();
        m_queue.erase(key);
        const auto value = m_distance_map.at(key);
        m_distance_map.erase(key);
        return std::make_pair(key, value);
    }

    [[nodiscard]]
    size_type size() const
    {
        return m_queue.size();
    }

    [[nodiscard]]
    bool empty() const
    {
        return m_queue.empty();
    }

    /// Inserts or updates the given keys
    /// @details For n keys inside the queue and m keys in the given range, runs in O(m log m+n) time
    template <std::ranges::input_range KVs> requires Concepts::RangeYields<KVs, value_type>
    void insert_range(KVs&& kvs)
    {
        // remove keys
        for (const auto& [key, _] : kvs) {
            if (m_distance_map.contains(key))
                m_queue.erase(key);
        }
        // write new distances
        for (const auto& [key, value] : kvs) {
            m_distance_map[key] = value;
        }
        // reinsert every key
        for (const auto& [key, _] : kvs) {
            m_queue.insert(key);
        }
    }

    /// Removes the given keys
    /// @details For n keys inside the queue and m keys in the given range, runs in O(m log n) time
    template <std::ranges::input_range KeyRange> requires Concepts::RangeYields<KeyRange, Key>
    void remove_range(KeyRange&& keys)
    {
        for (const auto& key : keys) {
            if (m_distance_map.contains(key))
                m_queue.erase(key);
        }
        for (const auto& key : keys) {
            m_distance_map.erase(key);
        }
    }
};
}

#endif //GEL_ASSOCIATIVECONTAINERS_H
