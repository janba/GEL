/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "TriMesh.h"
#include "../CGLA/Vec3f.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace CGLA;

namespace Geometry {
    
	namespace
	{
		string get_path(const string& _filename)
		{
			// Make sure we only have slashes and not backslashes
			string filename = _filename;
			replace(filename.begin(),filename.end(),'\\','/');
			
			// Find the last occurrence of slash.
			// Everything before is path.
			size_t n = filename.rfind("/");
			if(n > filename.size())
				return "./";
			string pathname(filename,0,n);
			pathname.append("/");
			return pathname;
		}
	}
	class TriMeshObjLoader
	{
		TriMesh *mesh;
		std::string pathname;
        
		int get_vert(int i) {
			assert(i!=0);
			if (i<0) {
				return mesh->geometry.no_vertices()+i;
			} else
				return i-1;
		}
        
		int get_normal(int i) {
			if (i<0) {
				return mesh->normals.no_vertices()+i;
			} else
				return i-1;
		}
        
		int get_texcoord(int i) {
			if (i<0) {
				return mesh->texcoords.no_vertices()+i;
			} else
				return i-1;
		}

    void read_material_library(const string& filename, vector<Material>& materials);

	public:
        
		TriMeshObjLoader(TriMesh *_mesh): mesh(_mesh) {}
		
		void load(const std::string& filename);
        
    void load_material_library(const string& filename, vector<Material>& materials)
    {
      pathname = get_path(filename);
      read_material_library(filename, materials);
    }        
	};
	
	void TriMeshObjLoader::read_material_library(const string& filename, vector<Material>& materials)
	{
    string fn = pathname + filename;
    ifstream material_file(fn.data());  
    if(!material_file)
    {
      cerr << "Could not open " << filename << endl;
      return;
    }
    materials.resize(1);

    string buf;
    unsigned int nummaterials=0;
    while(material_file >> buf)
    {
      switch(buf[0])
      {
      case 'n':
        ++nummaterials;
        materials.push_back(Material());
        material_file >> materials[nummaterials].name;
        break;
      case 'N':
        switch(buf[1])
        {
        case 's':
          material_file >> materials[nummaterials].shininess;
          materials[nummaterials].shininess *= 128.0f/1000.0f;
          break;
        case 'i':
          material_file >> materials[nummaterials].ior;
          break;
        }
        break;
      case 'K':
        switch(buf[1])
        {
        case 'd':
          material_file >> materials[nummaterials].diffuse[0]
                        >> materials[nummaterials].diffuse[1]
                        >> materials[nummaterials].diffuse[2];
          break;
        case 's':
          material_file >> materials[nummaterials].specular[0]
                        >> materials[nummaterials].specular[1]
                        >> materials[nummaterials].specular[2];
          break;
        case 'a':
          material_file >> materials[nummaterials].ambient[0]
                        >> materials[nummaterials].ambient[1]
                        >> materials[nummaterials].ambient[2];
          break;
        }
        break;
      case 'T':
        switch(buf[1])
        {
        case 'f':
          material_file >> materials[nummaterials].transmission[0]
                        >> materials[nummaterials].transmission[1]
                        >> materials[nummaterials].transmission[2];
          break;
        }
        break;
      case 'i':
        material_file >> materials[nummaterials].illum;
        break;
      case 'm': // Map ... all maps are treated equally.
        material_file >> materials[nummaterials].tex_name;
        materials[nummaterials].has_texture = true;
        materials[nummaterials].tex_path = pathname;
        break;
      case '#':
      default:
        material_file.ignore(1024, '\n');
        break;
      }
    }
  }
    
  void TriMeshObjLoader::load(const std::string& filename)
  {
    pathname = get_path(filename);
    ifstream obj_file(filename.data());
    if(!obj_file)
    {
			cerr << "File " << filename << " does not exist" << endl;
      exit(0);
    }

    mesh->materials.resize(1);
    string buf;
    int current_material=0;
    while(obj_file >> buf)
    {
      switch(buf[0])
      {
      case 'v': // v, vn, vt
        {
          if(buf == "v")
          {
            Vec3f v_geo;
            obj_file >> v_geo;
            mesh->geometry.add_vertex(v_geo);
          }
          else if(buf == "vn")
          {
            Vec3f v_norm;
            obj_file >> v_norm;
            mesh->normals.add_vertex(v_norm);
          }
          else if(buf == "vt")
          {
            Vec3f v_tex;
            obj_file >> v_tex[0] >> v_tex[1];
            v_tex[2]=1;
            mesh->texcoords.add_vertex(v_tex);
          }
        }
        break;
      case 'f':
        {
          char line[1024];
                        
          vector<int> v_indices;
          vector<int> t_indices;
          vector<int> n_indices;
                        
          obj_file.getline(line, 1022);
          char* pch = strtok(line, " \t");
          while(pch != 0)
          {
            int v,t,n;
            if(sscanf(pch, "%d/%d/%d", &v, &t, &n)==3)
            {
              v_indices.push_back(get_vert(v));
              t_indices.push_back(get_texcoord(t));
              n_indices.push_back(get_normal(n));
            }
            else if(sscanf(pch, "%d//%d", &v, &n)==2)
            {
              v_indices.push_back(get_vert(v));
              n_indices.push_back(get_normal(n));
            }
            else if(sscanf(pch, "%d/%d", &v, &t)==2)
            {
              v_indices.push_back(get_vert(v));
              t_indices.push_back(get_texcoord(t));
            }
            else if(sscanf(pch, "%d", &v)==1)
              v_indices.push_back(get_vert(v));
            pch = strtok(0, " \t");
          }
          if(v_indices.size()>=3)
            for(size_t i=0;i<=(v_indices.size()-3);++i)
            {
              int idx = mesh->geometry.add_face(Vec3i(v_indices[0], v_indices[i+1], v_indices[i+2]));
              if(t_indices.size())
              mesh->texcoords.add_face(Vec3i(t_indices[0], t_indices[i+1], t_indices[i+2]), idx);
              if(n_indices.size())
              mesh->normals.add_face(Vec3i(n_indices[0], n_indices[i+1], n_indices[i+2]), idx);
              mesh->mat_idx.push_back(current_material);
            }
        }
        break;
      case 'm':
        obj_file >> buf;
        read_material_library(buf, mesh->materials);
        break;
      case 'u':
        obj_file >> buf;
        current_material = mesh->find_material(buf);
        break;
      case '#':
      default:
        obj_file.ignore(1024, '\n');
        break;                        
      }
    }
  }
    
	void obj_load(const string& filename, TriMesh &mesh)
	{
		TriMeshObjLoader loader(&mesh);
		loader.load(filename);
	}
    
  /// Load materials from an MTL file
  void mtl_load(const std::string& filename, std::vector<Material>& materials)
  {
    TriMeshObjLoader loader(0);
    loader.load_material_library(filename, materials);
  }    
}
