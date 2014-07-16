/*
 *  harmonics.h
 *  GEL
 *
 *  Created by J. Andreas Bærentzen on 01/09/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __MESHEDIT_HARMONICS_H__
#define __MESHEDIT_HARMONICS_H__

#include "../CGLA/Vec3d.h"
#include "../HMesh/Manifold.h"
#include "../HMesh/AttributeVector.h"
#include "../LinAlg/Matrix.h"
#include "../LinAlg/Vector.h"

namespace  HMesh {
    
    class Harmonics
    {
        HMesh::Manifold& mani;
        HMesh::VertexAttributeVector<int> vtouched;
        
        int maximum_eigenvalue;
        
        bool is_initialized;
        
        std::vector<CGLA::Vec3d> proj;
        
        LinAlg::CMatrix Q;
        LinAlg::CVector qnorm;
        LinAlg::CVector V;
        LinAlg::CVector S;
        
        std::vector<float> max_eig_values;
        
        void make_laplace_operator();
        void make_laplace_operator_sparse();
        
    public:
        
        std::vector<CGLA::Vec3d> analyze_signal(const HMesh::VertexAttributeVector<CGLA::Vec3d>& sig);
        
        HMesh::VertexAttributeVector<CGLA::Vec3d> reconstruct_signal(const std::vector<CGLA::Vec3d>& sig_proj, int);
        
        /// Initial analysis of harmonics
        Harmonics(HMesh::Manifold& mani);
        
        /// Add a frequency to mesh reconstruction
        void add_frequency(int f, float scale = 1.0f);
        
        /// Reset the shape to use 0 eigenvalues
        void reset_shape();
        
        /// Do a partial reconstruct with an interval of eigenvalues
        void partial_reconstruct(int E0, int E1, float scale=1.0f);
        
        double compute_adf(HMesh::VertexAttributeVector<double>& adf, double t, double fiedler_boost=0);

        double compute_esum(HMesh::VertexAttributeVector<double>& adf, int e0, int e1);
        
    };
    
}
#endif

