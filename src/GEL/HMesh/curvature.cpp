/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/HMesh/curvature.h>

#include <iostream>
#include <GEL/CGLA/CGLA.h>

#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/AttributeVector.h>
#include <GEL/HMesh/x3d_save.h>
#include <GEL/HMesh/x3d_load.h>
#include <GEL/HMesh/obj_load.h>
#include <GEL/HMesh/mesh_optimization.h>

using namespace std;
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
                for(auto v: m.vertices()){
                    new_vec[v] = vec[v];
                    for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw()){
                        new_vec[v] += vec[w.vertex()];
                    }
                    new_vec[v] /= (valency(m, v) + 1.0);
                }
                swap(vec,new_vec);
            }		
        }
    

    }

    void smooth_vectors_on_mesh(const Manifold& m, VertexAttributeVector<Vec3d>& vec, int smooth_steps)
    {
        for(int iter=0; iter<smooth_steps; ++iter){
            VertexAttributeVector<Vec3d> new_vec = vec;
            for(auto v: m.vertices()){
                for(auto vn: m.incident_vertices(v)) {
                    double sgn = dot(vec[vn], vec[v]) < 0 ? -1: 1;
                    new_vec[v] += sgn*vec[vn];
                }
                new_vec[v] /= (valency(m, v) + 1.0);
            }
            swap(vec,new_vec);
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

            if(a0<(M_PI/2.0) && a1<(M_PI/2.0) && a2<(M_PI/2.0)) // f is non-obtuse
                area_mixed += (1.0/8) * 
                    ((1.0/tan(a1)) * sqr_length(v2-v0) + 
                    (1.0/tan(a2)) * sqr_length(v1-v0));
            else if (a0>=(M_PI/2.0)) // a0 is the obtuse angle
                area_mixed += f_area/2;
            else
                area_mixed += f_area/4;
        }
        return area_mixed;
    }

    double barycentric_area(const Manifold& m, VertexID v)
    {
        double barea = 0;
        //For each triangle T from the 1-ring neighborhood of x
        for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw())
        if(w.face() != InvalidFaceID) {
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
        return curv_normal / (4.0*mixed_area(m, v));
    }

    double sum_curvatures(const Manifold& m, VertexAttributeVector<double>& curvature)
    {
        double sum = 0;
        for(auto v: m.vertices()){
            if(boundary(m, v))
                continue;	
            sum += curvature[v] * mixed_area(m, v);
        }
        return sum;
    }

    double angle_defect(const Manifold& m, VertexID v)
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
        return (2*M_PI - angle_sum);
        
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
        Vec3d centre(m.pos(v));

        vector<Vec3d> points;
        for (auto vn: m.incident_vertices(v))
            points.push_back(m.pos(vn));

        int N = int(points.size());
        vector<Vec3d> A(N);
        vector<double> b(N);
        for(int n = 0; n < N; ++n){
            Vec3d p = (points[n]-centre);
            double x = dot(p,X);
            double y = dot(p,Y);
            A[n] = Vec3d(0.5*x*x, x*y, 0.5*y*y);
            b[n] = dot(p,Norm);
        }
        try {
            Vec3d x = ls_solve(A,b);
            // Finally compute the shape tensor from the coefficients
            // using the first and second fundamental forms.
            curv_tensor = - Mat2x2d(x[0],x[1],x[1],x[2]);
        }
        catch (const Mat3x3fSingular& e) {
            cout << "Caught exception" <<endl;
        }
        curv_tensor = Mat2x2d(0.0);
    }

    PrincipalCurvatures principal_curvatures( const Manifold& m, VertexID v) {

        Mat2x2d tensor;
        Mat3x3d frame;
        curvature_tensor_paraboloid(m, v, tensor, frame);
        
        Mat2x2d Q, L;
        int s = power_eigensolution(tensor, Q, L);
        
        int max_idx = 0;
        int min_idx = 1;
        
        if(abs(L[max_idx][max_idx])<abs(L[min_idx][min_idx])) 
            swap(max_idx, min_idx);
        
        Mat3x3d frame_t = transpose(frame);
        
        PrincipalCurvatures pc = {
            cond_normalize(frame_t * Vec3d(Q[min_idx][0], Q[min_idx][1], 0)),
            cond_normalize(frame_t * Vec3d(Q[max_idx][0], Q[max_idx][1], 0)),
            L[min_idx][min_idx],
            L[max_idx][max_idx]
        };

        return pc;
    }

    void curvature_tensors_from_edges(const Manifold& m, VertexAttributeVector<Mat3x3d>& curvature_tensors)
    {
        for(auto v: m.vertices())
            curvature_tensors[v] = curvature_tensor_from_edges(m, v);
    }

    void smooth_curvature_tensors(const Manifold& m, VertexAttributeVector<Mat3x3d>& curvature_tensors)
    {
        assert(curvature_tensors.size() == m.allocated_vertices());
        VertexAttributeVector<Mat3x3d> tmp_curvature_tensors;
        double tmp_area;

        for(auto v: m.vertices()){
            if(boundary(m, v))
                continue;
            double a = mixed_area(m, v);
            tmp_curvature_tensors[v] = curvature_tensors[v] * a;
            tmp_area = a;
            for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw()){
                if(!boundary(m, w.vertex())){
                    double a = mixed_area(m, w.vertex());
                    tmp_curvature_tensors[v] += curvature_tensors[w.vertex()]*a;
                    tmp_area += a;
                }
                tmp_curvature_tensors[v] /= tmp_area;
            }
        }
        curvature_tensors = std::move(tmp_curvature_tensors);
    }

    void gaussian_curvature_angle_defects(const Manifold& m, VertexAttributeVector<double>& curvature, int smooth_steps)
    {
        for(auto v: m.vertices())
            curvature[v] = gaussian_curvature_angle_defect(m, v);

        smooth_something_on_mesh(m, curvature, smooth_steps);
    }

    void mean_curvatures(const Manifold& m, VertexAttributeVector<double>& curvature, int smooth_steps)
    {
        for(auto v: m.vertices())
			if(!boundary(m,v))
			{
				Vec3d N = -mean_curvature_normal(m, v);
				curvature[v] = length(N) * sign(dot(N,Vec3d(normal(m, v))));
			}	
        smooth_something_on_mesh(m, curvature, smooth_steps);	
    }


    void curvature_paraboloids( const Manifold& m, 
                                VertexAttributeVector<Vec3d>& min_curv_direction, 
                                VertexAttributeVector<Vec3d>& max_curv_direction,
                                VertexAttributeVector<Vec2d>& curvature)
    {

        for(auto v: m.vertices()){
           
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
        for(auto v: m.vertices()){
            Mat3x3d C,Q,L;
            C = curvature_tensors[v];
            int s = power_eigensolution(C, Q, L);
            Vec3d dmin, dmax;
            if(s == 0)
            {
                Vec3d n(normal(m, v));
                orthogonal(n, dmin, dmax);
                curvature[v] = Vec2d(0);
                cout << " rank 0 " << endl;
            }
            else if(s == 1)
            {
                Vec3d n(normal(m, v));
                dmin = normalize(Q[0]);
                dmax = cross(n, dmin);
                curvature[v] = Vec2d(0);
                cout << " rank 1 " << endl;
            }
            else
            {
                Vec2d l(fabs(L[0][0]), fabs(L[1][1]));

                int max_idx = 0;
                int min_idx = 1;

                if(abs(l[max_idx]) < abs(l[min_idx])) swap(max_idx, min_idx);
                
                // Yes - the biggest eigenvalue corresponds to the min direction
                // and vice versa.
                dmin = normalize(Q[max_idx]);
                dmax = normalize(Q[min_idx]);

                curvature[v][0] = L[min_idx][min_idx];
                curvature[v][1] = L[max_idx][max_idx];

            }
            min_curv_direction[v] = dmin;
            max_curv_direction[v] = dmax;
            max_val = max(fabs(curvature[v][1]), max_val);

        }
        //scal = 1.0/max_val;
    }
}

