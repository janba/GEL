/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file mesh_optimization.h
 * @brief Functions for energy minimization based mesh optimization.
 */

#ifndef __HMESH_MESH_OPTIMIZATION_H
#define __HMESH_MESH_OPTIMIZATION_H

#include <algorithm>
#include "Manifold.h"
#include "../CGLA/Vec3d.h"
#include "../Geometry/Implicit.h"

namespace HMesh
{
    // forward declarations
    //class Manifold;
    //class HalfEdgeID;

    /// This class represents the energy of an edge. It is used in optimization schemes where edges are swapped (aka flipped). */
    class EnergyFun
    {
    public:
        virtual double delta_energy(const Manifold& m, HalfEdgeID h) const = 0;
        virtual double energy(const Manifold& m, HalfEdgeID h) const {return 0;}
    };
	
	class MinAngleEnergy: public EnergyFun
	{
		double min_angle(const CGLA::Vec3d& v0, const CGLA::Vec3d& v1, const CGLA::Vec3d& v2) const;
		double thresh;
		
	public:
		
		MinAngleEnergy(double _thresh): thresh(_thresh) {}
		
		double delta_energy(const HMesh::Manifold& m, HMesh::HalfEdgeID h) const;
	};
	
	class DihedralEnergy: public EnergyFun
	{
		const double gamma;
		const bool use_alpha;
		
		double cos_ang(const CGLA::Vec3d& n1, const CGLA::Vec3d& n2) const
		{ 
			return (std::max)(-1.0, (std::min)(1.0, CGLA::dot(n1, n2)));
		}
		
		double edge_alpha_energy(CGLA::Vec3d v1, CGLA::Vec3d v2, double ca) const
		{
			return pow(CGLA::length(v1-v2)*(acos(ca)), 1.0f/gamma); 
		}
		
		double edge_c_energy(CGLA::Vec3d v1, CGLA::Vec3d v2, double ca) const
		{
			return pow(length(v1-v2)*(1-ca), 1.0f/gamma); 
		}

		void compute_angles(const HMesh::Manifold & m, HMesh::HalfEdgeID h) const;
		


		mutable double ab_12;
		mutable double ab_a1;
		mutable double ab_b1;
		mutable double ab_2c;
		mutable double ab_2d;
		
		mutable double aa_12;
		mutable double aa_b1;
		mutable double aa_c1;
		mutable double aa_2a;
		mutable double aa_2d;
		
	public:
		
		DihedralEnergy(double _gamma = 4.0, bool _use_alpha=false): 
		gamma(_gamma), use_alpha(_use_alpha) {}
		
		double energy(const HMesh::Manifold& m, HMesh::HalfEdgeID h) const;
		
		double delta_energy(const HMesh::Manifold& m, HMesh::HalfEdgeID h) const;
	
		double min_angle(const HMesh::Manifold& m, HMesh::HalfEdgeID h) const
		{
			compute_angles(m, h);
			return (std::min)((std::min)((std::min)((std::min)(aa_12, aa_b1), aa_c1), aa_2a), aa_2d);
		}		
	};
	
	class CurvatureEnergy: public EnergyFun
	{
        mutable std::vector<CGLA::Vec3d> va_ring_bef;
        mutable std::vector<CGLA::Vec3d> va_ring_aft;
        mutable std::vector<CGLA::Vec3d> vb_ring_bef;
        mutable std::vector<CGLA::Vec3d> vb_ring_aft;
        mutable std::vector<CGLA::Vec3d> vc_ring_bef;
        mutable std::vector<CGLA::Vec3d> vc_ring_aft;
        mutable std::vector<CGLA::Vec3d> vd_ring_bef;
        mutable std::vector<CGLA::Vec3d> vd_ring_aft;

		double abs_mean_curv(const CGLA::Vec3d& v, const std::vector<CGLA::Vec3d>& ring) const;
	public:
		double delta_energy(const HMesh::Manifold& m, HMesh::HalfEdgeID h) const;
	};		
    
    class ValencyEnergy: public EnergyFun
	{
	public:
		double delta_energy(const Manifold& m, HalfEdgeID h) const;
	};


    /// Optimize in a greedy fashion.
    void priority_queue_optimization(Manifold& m, const EnergyFun& efun);

    /// Optimize with simulated annealing. Avoids getting trapped in local minima
    void simulated_annealing_optimization(Manifold& m, const EnergyFun& efun, int max_iter=10000);

    /// Minimize the angle between adjacent triangles. Almost the same as mean curvature minimization 
    void minimize_dihedral_angle(Manifold& m, int max_iter=10000, bool anneal=false, bool alpha=false, double gamma=4.0);

    /// Minimizes mean curvature. This is really the same as dihedral angle optimization except that we weight by edge length 
    void minimize_curvature(Manifold& m, bool anneal=false);

    /// Minimizes gaussian curvature. Probably less useful than mean curvature.
    void minimize_gauss_curvature(Manifold& m, bool anneal=false);

    /// Maximizes the minimum angle of triangles. Makes the mesh more Delaunay.
    void maximize_min_angle(Manifold& m, float thresh, bool anneal=false);

    /// Tries to achieve valence 6 internally and 4 along edges.
    void optimize_valency(Manifold& m, bool anneal=false);

    /// Make radom flips. Useful for generating synthetic test cases.
    void randomize_mesh(Manifold& m, int max_iter);

    /// Perform many operations in order to equalize edge lengths.
    void edge_equalize(HMesh::Manifold& m, const Geometry::Implicit& imp, int max_iter);

}


#endif
