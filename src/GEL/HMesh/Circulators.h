//
//  Circulators.h
//  GEL
//
//  Created by Andreas Bærentzen on 24/10/2020.
//  Copyright © 2020 J. Andreas Bærentzen. All rights reserved.
//
#include <iterator>

#include "ConnectivityKernel.h"
#include "Iterators.h"


#ifndef HMesh_Circulators_h
#define HMesh_Circulators_h

namespace HMesh {

    template<typename  ITEM>
    class VertexCirculator {
        const ConnectivityKernel* ck;
        HalfEdgeID he;
        bool begun;
    public:
        using value_type = ItemID<ITEM>;

        VertexCirculator(const ConnectivityKernel& _ck, HalfEdgeID _he, bool _begun=false):
            ck(&_ck), he(_he), begun(_begun) {}
        
        value_type operator*() const;

        VertexCirculator& operator++() {
            begun = true;
            he = ck->opp(ck->prev(he));
            return *this;
        }
        
        VertexCirculator& operator++(int) {
            auto tmp(*this);
            ++(*this);
            return tmp;
        }
        
        bool operator !=(const VertexCirculator& c2) {
            if (he != c2.he)
                return true;
            if(begun != c2.begun)
                return true;
            return false;
        }
        
        bool operator ==(const VertexCirculator& c2) {
            if (he == c2.he && begun == c2.begun)
                return true;
            return false;
        }
    };

    template<>
    inline VertexCirculator<Vertex>::value_type VertexCirculator<Vertex>::operator*() const {
        return ck->vert(he);
    }

    template<>
    inline VertexCirculator<HalfEdge>::value_type VertexCirculator<HalfEdge>::operator*() const {
        return he;
    }

    template<>
    inline VertexCirculator<Face>::value_type VertexCirculator<Face>::operator*() const {
        return ck->face(he);
    }



    template<typename  ITEM>
    class FaceCirculator {
        const ConnectivityKernel* ck;
        HalfEdgeID he;
        bool begun;
    public:
        using value_type = ItemID<ITEM>;

        FaceCirculator(const ConnectivityKernel& _ck, HalfEdgeID _he, bool _begun=false):
            ck(&_ck), he(_he), begun(_begun) {}
        
        value_type operator*() const {
            return ck->vert(he);
        }
        
        FaceCirculator& operator++() {
            begun = true;
            he = ck->next(he);
            return *this;
        }
        
        FaceCirculator& operator++(int) {
            auto tmp(*this);
            ++(*this);
            return tmp;
        }
        
        bool operator !=(const FaceCirculator& c2) {
            if (he != c2.he)
                return true;
            if(begun != c2.begun)
                return true;
            return false;
        }
        
        bool operator ==(const FaceCirculator& c2) {
            if (he == c2.he && begun == c2.begun)
                return true;
            return false;
        }
    };

    template<>
    inline FaceCirculator<Vertex>::value_type FaceCirculator<Vertex>::operator*() const {
        return ck->vert(he);
    }

    template<>
    inline FaceCirculator<HalfEdge>::value_type FaceCirculator<HalfEdge>::operator*() const {
        return he;
    }

    template<>
    inline FaceCirculator<Face>::value_type FaceCirculator<Face>::operator*() const {
        return ck->face(he);
    }


}


#endif /* Circulator_h */
