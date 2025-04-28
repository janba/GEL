//
//  AttribVec.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 01/03/2016.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#ifndef UTIL__AttribVec_h
#define UTIL__AttribVec_h

#include <assert.h>
#include <map>
#include "Serialization.h"
namespace Util {


template<typename IndexT>
inline std::enable_if_t<std::is_integral_v<IndexT>, size_t> get_index(IndexT i) {return i;}

template<typename IndexT>
inline std::enable_if_t<!std::is_integral_v<IndexT>, size_t> get_index(IndexT i) {return i.index;}


template<typename IndexT, typename ValT>
class AttribVec
{
    std::vector<ValT> items;
    ValT default_value;
    
public:
    /// Construct from  size and optional default value (size should be identical to number of entities.
    AttribVec(size_t _size, ValT _default_value) :
    items(_size, _default_value),
    default_value(_default_value) {}
    
    /// Construct from optional default value.
    AttribVec(ValT _default_value):
    default_value(_default_value) {}
    
    /// Default constructor
    AttribVec() {}
    
    /// Return raw pointer to data
    ValT* data() {return items.data();}
    
    /// Return raw const pointer to data
    const ValT* data() const {return items.data();}
    
    /** Returns  value corresponding to index. Note that an object and not a
     const reference is returned. This is because the vector is resizable and the
     reference might be invalidated at surprising times. If you really need a reference
     call the at function. */
    ValT operator [](IndexT idx) const {
        return items.at(get_index(idx));
    }
    
    /** This at function returns a const reference. Use at own risk. */
    const ValT& at(IndexT idx) const {
        return items.at(get_index(idx));
    }
    
    /// reference to item given by ID. Note that this function may resize the vector.
    ValT& operator [](IndexT idx) {
        if (get_index(idx) >= items.size())
            resize(get_index(idx)+1);
        return items.at(get_index(idx));
    }
    
    /** resize the vector (may be necessary if associated container size grows). This
     version of the function takes a new default value as argument. This will override and
     overwrite the previous default value. Rarely useful since the default value can be set
     in the constructor. */
    void resize(size_t _size, ValT _default_value) {
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
    void cleanup(const std::map<IndexT, IndexT>& remap) {
        std::vector<ValT> new_items(remap.size());
        for(const auto& mapping : remap){
            assert(get_index(mapping.second) < remap.size());
            if(get_index(mapping.first) < items.size())
                new_items[get_index(mapping.second)] = items[get_index(mapping.first)];
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

}


#endif /* AttribVec_h */
