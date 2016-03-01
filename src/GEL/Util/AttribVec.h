//
//  AttribVec.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 01/03/2016.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#ifndef UTIL__AttribVec_h
#define UTIL__AttribVec_h

namespace Util {
    
    template<typename IndexT, typename ValT>
    class AttribVec
    {
        std::vector<ValT> items;
        
    public:
        /// Construct from optional size and item (size should be identical to number of entities in associated container
        AttribVec(size_t _size, ValT item = ValT()): items(_size, item) {}
        
        AttribVec() {}
        
        /// const reference to item given by ID
        const ValT& operator [](IndexT id) const {
            assert(id<items.size());
            return items[id];
        }
        
        /// reference to item given by ID
        ValT& operator [](IndexT id) {
            if(id>= items.size())
                items.resize(id+1);
            return items[id];
        }
        
        /// number of attribute items in vector
        size_t size() const {
            return items.size();
        }
        
        /// clear the vector
        void clear() {
            items.clear();
        }
    };
}


#endif /* AttribVec_h */
