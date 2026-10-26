/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/HMesh/curvature.h>
#include <algorithm>
#include <random>
#include <iostream>
#include <cmath>
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
            if(w.face() != InvalidFaceID) {
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


    Vec3d mean_curvature_normal(const Manifold& m, VertexID v)
    {
        Vec3d curv_normal(0);
        if(!boundary(m, v)) {
            Vec3d vertex(m.pos(v));
            for(Walker walker = m.walker(v); !walker.full_circle(); walker = walker.circulate_vertex_ccw()){
                Vec3d nbr(m.pos(walker.vertex()));
                Vec3d left(m.pos(walker.next().vertex()));
                Vec3d right(m.pos(walker.opp().next().vertex()));
                
                double d_left = dot(cond_normalize(nbr-left),cond_normalize(vertex-left));
                double d_right = dot(cond_normalize(nbr-right),cond_normalize(vertex-right));
                double a_left  = acos(min(1.0, max(-1.0, d_left)));
                double a_right = acos(min(1.0, max(-1.0, d_right)));
                double w = sin(a_left + a_right) / (1e-300 + sin(a_left)*sin(a_right));
                curv_normal += w * (nbr-vertex);
            }
        }
        return -curv_normal / (4.0*mixed_area(m, v));
    }

    double mean_curvature(const Manifold& m, VertexID v)
    {
        Vec3d N = mean_curvature_normal(m, v);
        return length(N) * sign(dot(N,Vec3d(normal(m, v))));
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


    double gaussian_curvature(const Manifold& m, VertexID v)
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

    void gaussian_curvature(const Manifold& m, VertexAttributeVector<double>& curvature, int smooth_steps)
    {
        for(auto v: m.vertices())
            curvature[v] = gaussian_curvature(m, v);
        smooth_something_on_mesh(m, curvature, smooth_steps);
    }

    void mean_curvature(const Manifold& m, VertexAttributeVector<double>& curvature, int smooth_steps)
    {
        for(auto v: m.vertices())
            curvature[v] = mean_curvature(m, v);
        smooth_something_on_mesh(m, curvature, smooth_steps);	
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
    
    
    void curvature_tensor_paraboloid(const Manifold& m, VertexID v, Mat2x2d& curv_tensor, Mat3x3d& frame,
                                      bool use_opposite_vertices)
    {
        // Tangent frame of the mesh normal. Rows map world coordinates into this frame.
        Vec3d norm = Vec3d(normal(m, v));
        Vec3d X(norm[2]+0.721, norm[1]+0.163, norm[0]+0.542);
        Vec3d Y=normalize(cross(norm, X));
        X = normalize(cross(Y, norm));
        Vec3d centre(m.pos(v));
        frame = Mat3x3d(X,Y,norm);

        // The 1-ring, plus the vertex across each link edge when it exists.
        // Boundary edges, unused ids, and repeated vertices contribute nothing.
        vector<Vec3d> points;
        vector<VertexID> used;
        auto consider = [&](VertexID id) {
            if (id == v || !m.in_use(id))
                return;
            if (std::find(used.begin(), used.end(), id) != used.end())
                return;
            used.push_back(id);
            points.push_back(m.pos(id));
        };
        for (auto vn: m.incident_vertices(v))
            consider(vn);
        // Width of the Gaussian: mean squared distance of the direct 1-ring.
        // Opposite vertices are measured against this same width.
        double ring_sq = 0;
        const int ring_n = int(points.size());
        for (int n = 0; n < ring_n; ++n) {
            const Vec3d d = points[n] - centre;
            ring_sq += dot(d, d);
        }
        const double mean_sq = ring_n > 0 ? ring_sq / double(ring_n) : 0.0;
        if (use_opposite_vertices && !used.empty()) {
            for (Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw()) {
                if (!m.in_use(w.halfedge()) || w.face() == InvalidFaceID)
                    continue;
                const Walker link = w.next();
                if (!m.in_use(link.halfedge()))
                    continue;
                // opp() keeps the parent's circulation counter, so start a new
                // walk on the opposite face. A boundary edge has no face there.
                const HalfEdgeID across_id = link.opp().halfedge();
                if (!m.in_use(across_id))
                    continue;
                const Walker across = m.walker(across_id);
                if (across.face() == InvalidFaceID)
                    continue;
                const VertexID end_a = link.vertex();
                const VertexID end_b = across.vertex();
                for (Walker face = across; !face.full_circle(); face = face.circulate_face_ccw()) {
                    const VertexID opposite = face.vertex();
                    if (opposite == end_a || opposite == end_b)
                        continue;
                    consider(opposite);
                }
            }
        }
        const int N = int(points.size());

        // z = su x + sv y + (1/2) cu x^2 + cuv x y + (1/2) cv y^2.
        // The linear terms travel with the quadratic ones. A separate plane fit
        // absorbs curvature wherever the 1-ring is not centrosymmetric.
        vector<VecN<5>> rows(N);
        vector<double> rhs(N);
        for (int n = 0; n < N; ++n) {
            const Vec3d p = frame * (points[n] - centre);
            const double x = p[0];
            const double y = p[1];
            // Weighted least squares: row scale sqrt(w) makes the solver
            // minimize sum w r^2 with w = exp(-dist^2 / mean 1-ring dist^2).
            const double dist2 = dot(p, p);
            const double weight = mean_sq > 0.0 ? std::exp(-dist2 / mean_sq) : 1.0;
            const double scale = std::sqrt(weight);
            rows[n] = VecN<5>(0);
            rows[n][0] = scale * x;
            rows[n][1] = scale * y;
            rows[n][2] = scale * 0.5 * x * x;
            rows[n][3] = scale * x * y;
            rows[n][4] = scale * 0.5 * y * y;
            rhs[n] = scale * p[2];
        }
        const VecN<5> coef = ls_solve(rows, rhs);
        const double su = coef[0];
        const double sv = coef[1];
        const double cu = coef[2];
        const double cuv = coef[3];
        const double cv = coef[4];

        // First and second fundamental forms of the graph at the vertex.
        // The shape operator below uses the outward-sphere sign of the old
        // pure-quadratic fit (a convex cap with the normal pointing out is positive).
        const double E = 1.0 + su * su;
        const double F = su * sv;
        const double G = 1.0 + sv * sv;
        const double metric = std::sqrt(E * G - F * F);
        const double e = -cu / metric;
        const double f = -cuv / metric;
        const double g = -cv / metric;

        // Cholesky factor of the first fundamental form, I = R^T R, and the
        // symmetric shape operator in the orthonormal tangent coordinates w = R v.
        const double R00 = std::sqrt(E);
        const double R01 = F / R00;
        const double R11 = metric / R00;
        const Mat2x2d Rinv(1.0 / R00, -R01 / (R00 * R11),
                           0.0, 1.0 / R11);
        const Mat2x2d second(e, f, f, g);
        curv_tensor = transpose(Rinv) * second * Rinv;

        const Vec3d r_x(1.0, 0.0, su);
        const Vec3d r_y(0.0, 1.0, sv);
        const Vec3d n_local = cond_normalize(Vec3d(-su, -sv, 1.0));
        const Vec2d v0(Rinv[0][0], Rinv[1][0]);
        const Vec2d v1(Rinv[0][1], Rinv[1][1]);
        const Vec3d t0_local = v0[0] * r_x + v0[1] * r_y;
        const Vec3d t1_local = v1[0] * r_x + v1[1] * r_y;

        const Mat3x3d to_world = transpose(frame);
        const Vec3d n = cond_normalize(to_world * n_local);
        const Vec3d t0 = cond_normalize(to_world * t0_local);
        const Vec3d t1 = cond_normalize(cross(n, t0));
        frame = Mat3x3d(t0, t1, n);
    }

    PrincipalCurvatures principal_curvatures(const Manifold& m, VertexID v, bool use_opposite_vertices) {

        Mat2x2d tensor;
        Mat3x3d frame;
        curvature_tensor_paraboloid(m, v, tensor, frame, use_opposite_vertices);
        
        Mat2x2d Q, L;
        eigensolution_2x2(tensor, Q, L);
        
        int max_idx = 0;
        int min_idx = 1;
        
        if((L[max_idx][max_idx])<(L[min_idx][min_idx])) 
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

