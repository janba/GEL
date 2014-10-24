/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "curvature.h"

#include <iostream>
#include "../CGLA/eigensolution.h"
#include "../CGLA/Vec2d.h"
#include "../CGLA/Vec3d.h"
#include "../CGLA/Mat3x3d.h"
#include "../CGLA/Mat2x2d.h"
#include "../CGLA/Mat2x3d.h"

#include "Manifold.h"
#include "AttributeVector.h"
#include "x3d_save.h"
#include "x3d_load.h"
#include "obj_load.h"
#include "mesh_optimization.h"

#include "../LinAlg/Matrix.h"
#include "../LinAlg/Vector.h"
#include "../LinAlg/LapackFunc.h"

using namespace std;

using namespace LinAlg;
using namespace CGLA;
using namespace HMesh;

namespace HMesh
{
    namespace 
    {
        //double scal = 0.001;
        //double vector_scal = 0.001;

        template<class T> 
        void smooth_something_on_mesh(const Manifold& m, VertexAttributeVector<T>& vec, int smooth_steps)
        {
            for(int iter=0;iter<smooth_steps;++iter){
                VertexAttributeVector<T> new_vec(m.allocated_vertices());
                for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
                    new_vec[*v] = vec[*v];
                    for(Walker w = m.walker(*v); !w.full_circle(); w = w.circulate_vertex_cw()){
                        new_vec[*v] += vec[w.vertex()];
                    }
                    new_vec[*v] /= (valency(m, *v) + 1.0);
                }
                swap(vec,new_vec);
            }		
        }
    }

    double mixed_area(const Manifold& m, VertexID v)
    {
        double area_mixed = 0;
        //For each triangle T from the 1-ring neighborhood of x
        for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw()){
            double f_area = area(m, w.face());

            Vec3d v0(m.pos(v));
            Vec3d v1(m.pos(w.vertex()));
            Vec3d v2(m.pos(w.next().vertex()));

            double a0 = acos(dot(v1-v0, v2-v0)/(length(v1-v0)*length(v2-v0)));
            double a1 = acos(dot(v2-v1, v0-v1)/(length(v2-v1)*length(v0-v1)));
            double a2 = acos(dot(v0-v2, v1-v2)/(length(v0-v2)*length(v1-v2)));

            if(a0>(M_PI/2.0) && a1>(M_PI/2.0) && a2>(M_PI/2.0)) // f is non-obtuse
            {
                // Add Voronoi formula (see Section 3.3)
                area_mixed += (1.0/8) * 
                    ((1.0/tan(a1)) * sqr_length(v2-v0) + 
                    (1.0/tan(a2)) * sqr_length(v1-v0));
            }
            else // Voronoi inappropriate
            {
                // Add either area(f)/4 or area(f)/2
                area_mixed += f_area/3;
            }
        }
        return area_mixed;
    }

    double barycentric_area(const Manifold& m, VertexID v)
    {
        double barea = 0;
        //For each triangle T from the 1-ring neighborhood of x
        for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw()){
            barea += area(m, w.face())/3.0;
        }
        return barea;
    }

    void unnormalized_mean_curvature_normal(const Manifold& m, VertexID v, Vec3d& curv_normal, double& w_sum)
    {
        if(boundary(m, v))
            return;

        Vec3d vertex(m.pos(v));
        curv_normal = Vec3d(0);
        w_sum = 0;
        for(Walker walker = m.walker(v); !walker.full_circle(); walker = walker.circulate_vertex_ccw()){
            Vec3d nbr(m.pos(walker.vertex()));
            Vec3d left(m.pos(walker.next().vertex()));
            Vec3d right(m.pos(walker.opp().next().vertex()));

            double d_left = dot(cond_normalize(nbr-left),cond_normalize(vertex-left));
            double d_right = dot(cond_normalize(nbr-right),cond_normalize(vertex-right));
            double a_left  = acos(min(1.0, max(-1.0, d_left)));
            double a_right = acos(min(1.0, max(-1.0, d_right)));

//            double w = 1.0/(1e-300+tan(a_left));
//            w += 1.0/(1e-300+tan(a_right));
            double w = sin(a_left + a_right) / (1e-300 + sin(a_left)*sin(a_right));
            
//            double wl = dot(vertex-left, nbr-left)/length(cross(vertex-left, nbr-left));
//            double wr = dot(vertex-right, nbr-right)/length(cross(vertex-right, nbr-right));
//            double w = wl + wr;
            curv_normal += w * (nbr-vertex);
            w_sum += w;
        }

    }

    Vec3d mean_curvature_normal(const Manifold& m, VertexID v)
    {
        Vec3d curv_normal;
        double w_sum;
        unnormalized_mean_curvature_normal(m, v, curv_normal, w_sum);

        return curv_normal / (4*mixed_area(m, v));
    }

    double sum_curvatures(const Manifold& m, VertexAttributeVector<double>& curvature)
    {
        double sum = 0;
        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
            if(boundary(m, *v))
                continue;	
            sum += curvature[*v] * mixed_area(m, *v);
        }
        return sum;
    }


    double gaussian_curvature_angle_defect(const Manifold& m, VertexID v)
    {
        if(boundary(m, v))
            return 0;

        Vec3d vertex(m.pos(v));
        vector<Vec3d> edges;
        for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw()){
            Vec3d e(normalize(m.pos(w.vertex()) - vertex));
            edges.push_back(e);
        }
        size_t N=edges.size();
        double angle_sum = 0;
        for(size_t i = 0; i < N; ++i)
        {
            double dot_prod = 
                std::max(-1.0, std::min(1.0, dot(edges[i],edges[(i+1)%N])));
            angle_sum += acos(dot_prod);
        }
        return (2*M_PI - angle_sum)/mixed_area(m, v);

    }

    Mat3x3d curvature_tensor(const Manifold& m, HalfEdgeID h)
    {
        if(boundary(m, h))
            return Mat3x3d(0);

        Walker w = m.walker(h);
        Vec3d edge(m.pos(w.vertex()) - m.pos(w.opp().vertex()));
        double edge_len = length(edge);
        edge /= edge_len;

        Vec3d h_norm(normal(m, w.face()));
        Vec3d h_opp_norm(normal(m, w.opp().face()));

        Vec3d nc = cross(h_norm, h_opp_norm);

        double sign = (dot(nc, edge) >= 0) ? 1 : -1;
        double beta = asin(nc.length());

        Mat3x3d mat;
        outer_product(edge, edge, mat);
        return sign * edge_len * beta * mat;
    }

    Mat3x3d curvature_tensor_from_edges(const Manifold& m, VertexID v)
    {
        Mat3x3d curv_tensor(0);

        if(boundary(m, v))
            return curv_tensor;

        for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw())
            curv_tensor += 0.5*curvature_tensor(m, w.halfedge());

        curv_tensor /= mixed_area(m, v);

        return curv_tensor;
    }


    void curvature_tensor_paraboloid(const Manifold& m, VertexID v, Mat2x2d& curv_tensor, Mat3x3d& frame)
    {
        if(boundary(m, v))
            return;
        // First estimate the normal and compute a transformation matrix
        // which takes us into tangent plane coordinates.
        Vec3d Norm = Vec3d(normal(m, v));
        Vec3d X,Y;
        orthogonal(Norm,X,Y);
        frame = Mat3x3d(X,Y,Norm);
        Vec3d centre(m.pos(v));

        vector<Vec3d> points;
        for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw())
            points.push_back(Vec3d(m.pos(w.vertex())));

        int N = int(points.size());

        CVector b(N);
        // Compute the matrix of parameter values
        CMatrix PMat(N, 3);
        for(int i = 0; i < N; ++i){
            Vec3d p = frame * (points[i]-centre);
            b[i] = p[2];

            PMat.set(i,0,0.5*sqr(p[0]));
            PMat.set(i,1,p[0]*p[1]);
            PMat.set(i,2,0.5*sqr(p[1]));
        }

        // Compute the coefficients of the polynomial surface
        CVector x(3);
        x = LinearLSSolve(PMat,b);
        if(isnan(x[0])) cout << __LINE__ << " " << PMat << b << endl ;

        // Finally compute the shape tensor from the coefficients
        // using the first and second fundamental forms.
        curv_tensor = - Mat2x2d(x[0],x[1],x[1],x[2]);

    }

    void curvature_tensors_from_edges(const Manifold& m, VertexAttributeVector<Mat3x3d>& curvature_tensors)
    {
        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v)
            curvature_tensors[*v] = curvature_tensor_from_edges(m, *v);
    }

    void smooth_curvature_tensors(const Manifold& m, VertexAttributeVector<Mat3x3d>& curvature_tensors)
    {
        assert(curvature_tensors.size() == m.allocated_vertices());
        VertexAttributeVector<Mat3x3d> tmp_curvature_tensors(m.allocated_vertices());
        double tmp_area;

        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
            if(boundary(m, *v))
                continue;
            double a = mixed_area(m, *v);
            tmp_curvature_tensors[*v] = curvature_tensors[*v] * a;
            tmp_area = a;
            for(Walker w = m.walker(*v); !w.full_circle(); w = w.circulate_vertex_cw()){
                if(!boundary(m, w.vertex())){
                    double a = mixed_area(m, w.vertex());
                    tmp_curvature_tensors[*v] += curvature_tensors[w.vertex()]*a;
                    tmp_area += a;
                }
                tmp_curvature_tensors[*v] /= tmp_area;
            }
        }
        curvature_tensors = move(tmp_curvature_tensors);
    }

    void gaussian_curvature_angle_defects(const Manifold& m, VertexAttributeVector<double>& curvature, int smooth_steps)
    {
        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v)
            curvature[*v] = gaussian_curvature_angle_defect(m, *v);

        smooth_something_on_mesh(m, curvature, smooth_steps);
    }

    void mean_curvatures(const Manifold& m, VertexAttributeVector<double>& curvature, int smooth_steps)
    {
        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v)
			if(!boundary(m,*v))
			{
				Vec3d N = -mean_curvature_normal(m, *v);
				curvature[*v] = length(N) * sign(dot(N,Vec3d(normal(m, *v))));
			}	
        smooth_something_on_mesh(m, curvature, smooth_steps);	
    }


    void curvature_paraboloids( const Manifold& m, 
                                VertexAttributeVector<Vec3d>& min_curv_direction, 
                                VertexAttributeVector<Vec3d>& max_curv_direction,
                                VertexAttributeVector<Vec2d>& curvature)
    {

        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
            Mat2x2d tensor;
            Mat3x3d frame;
            curvature_tensor_paraboloid(m, *v, tensor, frame);

            Mat2x2d Q,L;
            int s = power_eigensolution(tensor, Q, L);

            if(s < 2)	
                cout << tensor << Q << L << endl;

            int max_idx = 0;
            int min_idx = 1;

            if(L[max_idx][max_idx]<L[min_idx][min_idx]) swap(max_idx, min_idx);

            Mat3x3d frame_t = transpose(frame);

            max_curv_direction[*v] = cond_normalize(frame_t * Vec3d(Q[max_idx][0], Q[max_idx][1], 0));

            min_curv_direction[*v] = cond_normalize(frame_t * Vec3d(Q[min_idx][0], Q[min_idx][1], 0));

            curvature[*v][0] = L[min_idx][min_idx];
            curvature[*v][1] = L[max_idx][max_idx];
        }
    }


    void curvature_from_tensors(const Manifold& m,
                                const VertexAttributeVector<Mat3x3d>& curvature_tensors,
                                VertexAttributeVector<Vec3d>& min_curv_direction,
                                VertexAttributeVector<Vec3d>& max_curv_direction,
                                VertexAttributeVector<Vec2d>& curvature)
    {
        assert(curvature_tensors.size() == m.allocated_vertices());

        double max_val = -1e30;
        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
            Mat3x3d C,Q,L;
            C = curvature_tensors[*v];
            int s = power_eigensolution(C, Q, L);
            Vec3d dmin, dmax;
            if(s == 0)
            {
                Vec3d n(normal(m, *v));
                orthogonal(n, dmin, dmax);
                curvature[*v] = Vec2d(0);
                cout << " rank 0 " << endl;
            }
            else if(s == 1)
            {
                Vec3d n(normal(m, *v));
                dmin = normalize(Q[0]);
                dmax = cross(n, dmin);
                curvature[*v] = Vec2d(0);
                cout << " rank 1 " << endl;
            }
            else
            {
                Vec2d l(fabs(L[0][0]), fabs(L[1][1]));

                int max_idx = 0;
                int min_idx = 1;

                if(l[max_idx] < l[min_idx]) swap(max_idx, min_idx);
                
                // Yes - the biggest eigenvalue corresponds to the min direction
                // and vice versa.
                dmin = normalize(Q[max_idx]);
                dmax = normalize(Q[min_idx]);

                curvature[*v][0] = L[min_idx][min_idx];
                curvature[*v][1] = L[max_idx][max_idx];

            }
            min_curv_direction[*v] = dmin;
            max_curv_direction[*v] = dmax;
            max_val = max(fabs(curvature[*v][1]), max_val);

        }
        //scal = 1.0/max_val;
    }
}

