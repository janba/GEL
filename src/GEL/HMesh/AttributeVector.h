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

#include <GEL/HMesh/ConnectivityKernel.h>
#include <GEL/Util/AttribVec.h>

namespace HMesh 
{
    /** HMesh entity attribute vectors. This class is used for storing all attributes
     associated with any type of mesh entity - be it Vertex, HalfEdge, or Face. Also the position attribute
     is stored in an AttributeVector class.      
     */

    template<typename ITEM>
    class VertexAttributeVector : public Util::AttribVec<VertexID, ITEM>
    {
    public:
        VertexAttributeVector(size_t _size, ITEM _default_value):
        Util::AttribVec<VertexID,ITEM>(_size, _default_value){}

        VertexAttributeVector(ITEM _default_value = ITEM()):
        Util::AttribVec<VertexID,ITEM>(_default_value){}
    };

    template<typename ITEM>
    class FaceAttributeVector : public Util::AttribVec<FaceID, ITEM>
    {
    public:
        FaceAttributeVector(size_t _size, ITEM _default_value):
        Util::AttribVec<FaceID,ITEM>(_size, _default_value){}

        FaceAttributeVector(ITEM _default_value = ITEM()):
        Util::AttribVec<FaceID,ITEM>(_default_value){}


    };

    template<typename ITEM>
    class HalfEdgeAttributeVector : public Util::AttribVec<HalfEdgeID, ITEM>
    {
    public:
        HalfEdgeAttributeVector(size_t _size, ITEM _default_value):
        Util::AttribVec<HalfEdgeID,ITEM>(_size, _default_value){}

        HalfEdgeAttributeVector(ITEM _default_value = ITEM()):
        Util::AttribVec<HalfEdgeID,ITEM>(_default_value){}

    };
}

#endif
