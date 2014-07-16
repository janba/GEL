/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "../GL/glew.h"

#include "../CGLA/Mat4x4f.h"
#include "../CGLA/Vec3d.h"
#include "../HMesh/Manifold.h"
#include "../Geometry/TriMesh.h"

#include "draw.h"
#include "SinglePassWireframeRenderer.h"
#include "IDBufferWireFrameRenderer.h"
#include "SOIL.h"

namespace GLGraphics
{
    using namespace CGLA;
    using namespace Geometry;
    using namespace HMesh;
    using namespace std;
    
    namespace
    {
        void set_material(const Geometry::Material& material)
        {
            if(material.has_texture && material.tex_id >=0)
            {
                glEnable(GL_TEXTURE_2D);
                glBindTexture(GL_TEXTURE_2D, material.tex_id);
                glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
            }
            else
                glDisable(GL_TEXTURE_2D);
            
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material.ambient);
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material.diffuse);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material.specular);
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material.shininess);
        }
    }
    
    void draw_ball() {
        static TriMesh mesh;
        if(mesh.geometry.no_vertices()==0)
        {
            GLfloat sphere_verts[] = {0.91358,0,0,-0.91358,0,0,0,0.91358,0,0,-0.91358,0,0,0,0.91358,0,0,-0.91358,0.54321,-0.54321,0.54321,0.54321,0.54321,-0.54321,-0.54321,0.54321,0.54321,-0.54321,-0.54321,-0.54321,0.54321,0.54321,0.54321,-0.54321,0.54321,-0.54321,-0.54321,-0.54321,0.54321,0.54321,-0.54321,-0.54321,0.444444,0,0.814815,0,-0.814815,0.444444,0.814815,-0.444444,0,0.444444,0,-0.814815,0,0.814815,-0.444444,0.814815,0.444444,0,-0.444444,0,0.814815,0,0.814815,0.444444,-0.814815,0.444444,0,-0.444444,0,-0.814815,0,-0.814815,-0.444444,-0.814815,-0.444444,0,0,0.444444,0.814815,0.814815,0,0.444444,0.444444,0.814815,0,0,0.444444,-0.814815,-0.814815,0,-0.444444,-0.444444,0.814815,0,0,-0.444444,0.814815,-0.814815,0,0.444444,-0.444444,-0.814815,0,0,-0.444444,-0.814815,0.814815,0,-0.444444,0.444444,-0.814815,0};
            GLuint sphere_indices[] = {28,11,15,27,5,15,33,7,15,33,13,16,35,4,16,38,7,16,38,14,17,37,1,17,28,7,17,37,14,18,36,6,18,30,8,18,30,12,19,32,3,19,29,8,19,29,11,20,28,1,20,37,8,20,34,13,21,33,5,21,27,9,21,27,11,22,29,3,22,32,9,22,32,12,23,31,2,23,34,9,23,31,12,24,30,6,24,36,10,24,36,14,25,38,4,25,35,10,25,35,13,26,34,2,26,31,10,26,22,9,27,21,5,27,15,11,27,15,7,28,17,1,28,20,11,28,20,8,29,19,3,29,22,11,29,19,8,30,18,6,30,24,12,30,24,10,31,26,2,31,23,12,31,23,9,32,22,3,32,19,12,32,16,7,33,15,5,33,21,13,33,21,9,34,23,2,34,26,13,34,26,10,35,25,4,35,16,13,35,25,10,36,24,6,36,18,14,36,18,8,37,20,1,37,17,14,37,17,7,38,16,4,38,25,14,38};
            
            for(int i=0;i<38;++i) {
                int idx0 = 3*i;
                mesh.geometry.add_vertex(Vec3f(sphere_verts[idx0],sphere_verts[idx0+1],sphere_verts[idx0+2]));
            }
            for(int i=0;i<72;++i) {
                int idx0 = 3*i;
                mesh.geometry.add_face(Vec3i(sphere_indices[idx0]-1,sphere_indices[idx0+1]-1,sphere_indices[idx0+2]-1));
            }
            mesh.compute_normals();
        }
        draw(mesh);
    }

    
    void draw(const Geometry::IndexedFaceSet& geometry)
    {
        glBegin(GL_TRIANGLES);
        for(int i=0;i<geometry.no_faces();i++)
        {
            Vec3i g_face = geometry.face(i);
            Vec3f vert0 = geometry.vertex(g_face[0]);
            Vec3f vert1 = geometry.vertex(g_face[1]);
            Vec3f vert2 = geometry.vertex(g_face[2]);
            Vec3f norm = normalize(cross(vert1-vert0, vert2-vert0));
            glNormal3fv(norm.get());
            glVertex3fv(vert0.get());
            glVertex3fv(vert1.get());
            glVertex3fv(vert2.get());
        }
        glEnd();
    }
    
    void draw(const Geometry::TriMesh& tm, bool per_vertex_norms)
    {
        int old_mat_idx = -1;
        glBegin(GL_TRIANGLES);
        for(int i=0;i<tm.geometry.no_faces();i++)
        {
            int new_mat_idx = i<static_cast<int>(tm.mat_idx.size()) ? tm.mat_idx[i] : -1;
            if(new_mat_idx != old_mat_idx)
            {
                glEnd();
                set_material(tm.materials[tm.mat_idx[i]]);
                glBegin(GL_TRIANGLES);
                old_mat_idx = new_mat_idx;
            }
            Vec3i n_face = tm.normals.face(i);
            Vec3i g_face = tm.geometry.face(i);
            Vec3i t_face = tm.texcoords.face(i);
            
            if(!per_vertex_norms)
            {
                Vec3f vert0 = tm.geometry.vertex(g_face[0]);
                Vec3f vert1 = tm.geometry.vertex(g_face[1]);
                Vec3f vert2 = tm.geometry.vertex(g_face[2]);
                Vec3f norm = normalize(cross(vert1-vert0, vert2-vert0));
                glNormal3fv(norm.get());
            }
            for(int j=0;j<3;j++)
            {
                if(per_vertex_norms && n_face != Geometry::NULL_FACE)
                {
                    Vec3f norm = tm.normals.vertex(n_face[j]);
                    glNormal3fv(norm.get());
                }
                if(t_face != Geometry::NULL_FACE)
                {
                    Vec3f texc = tm.texcoords.vertex(t_face[j]);
                    glTexCoord2fv(texc.get());
                }
                Vec3f vert = tm.geometry.vertex(g_face[j]);
                glVertex3fv(vert.get());
            }
        }
        glEnd();
        glDisable(GL_TEXTURE_2D);
    }
    
    void load_textures(Geometry::TriMesh& tm)
    {
        for(unsigned int i=0;i<tm.materials.size(); ++i)
        {
            Geometry::Material& mat = tm.materials[i];
            if(mat.tex_name != "")
            {
                string name = mat.tex_path + mat.tex_name;
                mat.tex_id = SOIL_load_OGL_texture(name.data(), 0, 0, SOIL_FLAG_TEXTURE_REPEATS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_POWER_OF_TWO);
            }
        }
    }

	/// Draw an object of type T which contains only triangles as wireframe. In practice T = Manifold or TriMesh.
	template<typename T>
    void draw_triangles_in_wireframe(T& m, bool per_vertex_norms, const CGLA::Vec3f& line_color)
	{
		static SinglePassWireframeRenderer swr;
		swr.enable(line_color);
		draw(m, per_vertex_norms);
		swr.disable();
	}
    
	template
    void draw_triangles_in_wireframe(HMesh::Manifold& m, bool per_vertex_norms, const CGLA::Vec3f& line_color);
    
	template
    void draw_triangles_in_wireframe(Geometry::TriMesh& m, bool per_vertex_norms, const CGLA::Vec3f& line_color);
    
    template<class T>
    void draw_wireframe_oldfashioned(const T& m, bool per_vertex_norms, const Vec3f& line_color)
    {
        // Store state that we change
        glPushAttrib(GL_POLYGON_BIT);
        GLboolean lights_on;
        glGetBooleanv(GL_LIGHTING, &lights_on);
        Vec4f current_color;
        glGetFloatv(GL_CURRENT_COLOR, &current_color[0]);
        
        // Draw filled
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        draw(m, per_vertex_norms);
        
        // Draw lines
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDisable(GL_LIGHTING);
        glEnable(GL_POLYGON_OFFSET_LINE);
        glPolygonOffset(0,-5);
        glColor3fv(line_color.get());
        draw(m, per_vertex_norms);
        
        // Put back old state
        glColor3fv(current_color.get());
        if(lights_on) glEnable(GL_LIGHTING);
        glPopAttrib();
    }
    template
    void draw_wireframe_oldfashioned(const HMesh::Manifold& m, bool per_vertex_norms, const Vec3f& line_color);
    
    template
    void draw_wireframe_oldfashioned(const Geometry::TriMesh& m, bool per_vertex_norms, const Vec3f& line_color);
    
    
    void draw(const Manifold& m, bool per_vertex_norms)
    {
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
            if(!per_vertex_norms)
                glNormal3dv(normal(m, *f).get());
            if(no_edges(m, *f)== 3)
                glBegin(GL_TRIANGLES);
            else
                glBegin(GL_POLYGON);
            
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                Vec3d n = normal(m, w.vertex());
                if(per_vertex_norms)
                    glNormal3dv(n.get());
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }
    }
    
    
    bool depth_pick(int x, int y, float& depth)
    {
        float new_depth;
        
        // Get the minimum and maximum depth values.
        float minmax_depth[2];
        glGetFloatv(GL_DEPTH_RANGE, minmax_depth);
        glPushAttrib(GL_PIXEL_MODE_BIT);
        glReadBuffer(GL_FRONT);
        // Read a single pixel at the position of the mouse cursor.
        glReadPixels(x, y, 1,1, GL_DEPTH_COMPONENT,
                     GL_FLOAT, (void*) &new_depth);
        glPopAttrib();
        
        // If the depth corresponds to the far plane, we clicked on the
        // background.
        if(new_depth == minmax_depth[1])
            return false;
        
        depth = new_depth;
        return true;
    }
    
    Vec3d screen2world(int x, int y, float depth)
    {
        // Enquire about the viewport dimensions
        GLint viewport[4];
        glGetIntegerv(GL_VIEWPORT, viewport);
        // Copy modelview matrix.
        double mvmat[16];
        glGetDoublev(GL_MODELVIEW_MATRIX, mvmat);
        
        // Copy the projection matrix. We assume it is unchanged.
        double prjmat[16];
        glGetDoublev(GL_PROJECTION_MATRIX, prjmat);
        
        // Now unproject the point from screen to world coordinates.
        double ox, oy, oz;
        gluUnProject(x, y,depth,
                     mvmat,prjmat,viewport,
                     &ox, &oy, &oz);
        
        return Vec3d(ox,oy,oz);
    }
    
    Vec3d world2screen(const Vec3d& p)
    {
        // Enquire about the viewport dimensions
        GLint viewport[4];
        glGetIntegerv(GL_VIEWPORT, viewport);
        // Copy modelview matrix.
        double mvmat[16];
        glGetDoublev(GL_MODELVIEW_MATRIX, mvmat);
        
        // Copy the projection matrix. We assume it is unchanged.
        double prjmat[16];
        glGetDoublev(GL_PROJECTION_MATRIX, prjmat);
        
        // Now unproject the point from screen to world coordinates.
        double ox, oy, oz;
        gluProject(p[0], p[1], p[2],
                   mvmat,prjmat,viewport,
                     &ox, &oy, &oz);
        
        return Vec3d(ox,oy,oz);
    }


    
    
    void draw(const Geometry::AABox& box)
    {
        glBegin(GL_QUADS);
        Vec3f norm_neg[] = {Vec3f(0,0,-1), Vec3f(-1,0,0), Vec3f(0,-1,0)};
        Vec3f norm_pos[] = {Vec3f(0,0, 1), Vec3f( 1,0,0), Vec3f(0, 1,0)};
        for(int j=0;j<3;++j)
        {
            glNormal3fv(norm_neg[j].get());
            Vec3f p = box.get_pmin();
            glVertex3f(p[0], p[1], p[2]);
            p[(j+1)%3] = box.get_pmax()[(j+1)%3];
            glVertex3f(p[0], p[1], p[2]);
            p[j] = box.get_pmax()[j];
            glVertex3f(p[0], p[1], p[2]);
            p[(j+1)%3] = box.get_pmin()[(j+1)%3];
            glVertex3f(p[0], p[1], p[2]);
        }
        glEnd();
        glBegin(GL_QUADS);
        for(int j=0;j<3;++j)
        {
            glNormal3fv(norm_pos[j].get());
            Vec3f p = box.get_pmax();
            glVertex3f(p[0], p[1], p[2]);
            p[j] = box.get_pmin()[j];
            glVertex3f(p[0], p[1], p[2]);
            p[(j+1)%3] = box.get_pmin()[(j+1)%3];
            glVertex3f(p[0], p[1], p[2]);
            p[j] = box.get_pmax()[j];
            glVertex3f(p[0], p[1], p[2]);
        }
        glEnd();
    }
    
    void draw(const Geometry::OBox& box)
    {
        Mat4x4f m = identity_Mat4x4f();
        copy_matrix(box.get_rotation(), m);
        glPushMatrix();
        glMultMatrixf(m.get());
        draw(box.get_aabox());
        glPopMatrix();
    }
    
    /** Draw the tree. The first argument is the level counter, the second
     argument is the level at which to stop drawing. */
    template <class BoxType>
    void draw(const Geometry::BoundingINode<BoxType>& node, int level, int max_level)
    {
        if(level == max_level)
        {
            draw(node);
            return;
        }
        node->left->draw(level + 1, max_level);
        node->right->draw(level + 1, max_level);
    }
    
    template <class BoxType>
    void draw(const Geometry::BoundingLNode<BoxType>& node, int level, int max_level)
    {
#if USE_LEAF_BOXES
        draw(node); 
#endif
    }
    
    template <class BoxType>
    void draw(const Geometry::BoundingTree<BoxType>& tree, int max_level)
    {
        draw(*tree.root, 0, max_level);
    }
}
