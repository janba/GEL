/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/cleanup.h>
#include <GEL/HMesh/load.h>
#include <GEL/HMesh/refine_edges.h>
#include <GEL/HMesh/stl_load.h>
#include <fstream>
#include <sstream>

using namespace std;
using namespace CGLA;

namespace HMesh {
    using std::string;

    void right_trim_stl(std::string& str) {
        if (!str.empty())
            while (isspace(str.back())) str.pop_back();
    }

    istream& get_multi_line_stl(istream& is, string& buf) {
        getline(is, buf);
        right_trim_stl(buf);
        if (!buf.empty())
            while (buf.back() == '\\') {
                buf.pop_back();
                string continuation;
                getline(is, continuation);
                right_trim_stl(continuation);
                buf += continuation;
            }
        return is;
    }

    bool stl_load_ascii(
        const std::string& filename, Manifold& m,
        VertexAttributeVector<int>& orig_vertex_indices) {
        ifstream stl_file(filename.data());

        if (stl_file) {
            string        buf;
            vector<Vec3d> vertices;
            vector<int>   faces;
            vector<int>   indices;
            Vec3d         normal;
            while (get_multi_line_stl(stl_file, buf)) {
                istringstream iss(std::move(buf));
                string        code;
                if (iss >> code) {
                    if (code == "facet") {
                        if (iss >> code && code == "normal") { iss >> normal; }
                    } else if (code == "vertex") {
                        Vec3d v;
                        iss >> v;
                        vertices.push_back(v);
                    } else if (code == "endfacet") {
                        // Determine the order of indices to match the normal
                        // direction
                        Vec3d n(normal);
                        Vec3d e1 = vertices[vertices.size() - 2] -
                                   vertices[vertices.size() - 3];
                        Vec3d e2 = vertices[vertices.size() - 1] -
                                   vertices[vertices.size() - 3];
                        Vec3d n2 = cross(e1, e2);

                        if (dot(n, n2) > 0) {
                            indices.push_back(vertices.size() - 3);
                            indices.push_back(vertices.size() - 2);
                            indices.push_back(vertices.size() - 1);
                        } else {
                            indices.push_back(vertices.size() - 3);
                            indices.push_back(vertices.size() - 1);
                            indices.push_back(vertices.size() - 2);
                        }
                        faces.push_back(3);
                    }
                }
            }
            // cout << "Loaded " << filename << " : " << vertices.size()
            //      << " vertices and " << faces.size() << " faces" << endl;
            m.clear();

            orig_vertex_indices = build(
                m, vertices.size(), reinterpret_cast<double*>(&vertices[0]),
                faces.size(), &faces[0], &indices[0]);

            return true;
        }
        return false;
    }

    bool stl_load_binary(
        const std::string& filename, Manifold& m,
        VertexAttributeVector<int>& orig_vertex_indices) {
        ifstream stl_file(filename.data(), ios::binary);
        if (stl_file) {

            char header[80];
            stl_file.read(header, 80);

            uint32_t num_triangles;
            stl_file.read(reinterpret_cast<char*>(&num_triangles), 4);

            vector<Vec3d> vertices;
            vector<int>   faces;
            vector<int>   indices;
            for (uint32_t i = 0; i < num_triangles; ++i) {
                Vec3f normal;
                stl_file.read(reinterpret_cast<char*>(&normal), 12);

                Vec3f v[3];
                for (int j = 0; j < 3; ++j) {
                    stl_file.read(reinterpret_cast<char*>(&v[j]), 12);
                    vertices.push_back(Vec3d(v[j]));
                }
                uint16_t attr;
                stl_file.read(reinterpret_cast<char*>(&attr), 2);

                // Determine the order of indices to match the normal direction
                Vec3d n(normal);
                Vec3d e1 = vertices[vertices.size() - 2] -
                           vertices[vertices.size() - 3];
                Vec3d e2 = vertices[vertices.size() - 1] -
                           vertices[vertices.size() - 3];
                Vec3d n2 = cross(e1, e2);
                if (dot(n, n2) > 0) {
                    indices.push_back(vertices.size() - 3);
                    indices.push_back(vertices.size() - 2);
                    indices.push_back(vertices.size() - 1);
                } else {
                    indices.push_back(vertices.size() - 3);
                    indices.push_back(vertices.size() - 1);
                    indices.push_back(vertices.size() - 2);
                }
                faces.push_back(3);
            }

            // cout << "Loaded " << filename << " : " << vertices.size()
            //      << " vertices and " << faces.size() << " faces" << endl;
            m.clear();

            orig_vertex_indices = build(
                m, vertices.size(), reinterpret_cast<double*>(&vertices[0]),
                faces.size(), &faces[0], &indices[0]);

            return true;
        }
        return false;
    }

    bool stl_load(
        const std::string& filename, Manifold& m,
        VertexAttributeVector<int>& orig_vertex_indices) {

        // Open the file
        ifstream stl_file(filename.data(), ios::binary);
        if (!stl_file) {
            cerr << "Could not open file " << filename << endl;
            return false;
        }

        // Read the first 5 bytes to determine if the file is ASCII or binary
        char header[5];
        stl_file.read(header, 5);
        stl_file.close();

        // Check if the file is ASCII or binary (ASCII files start with "solid")
        bool is_ascii = (header[0] == 's' || header[0] == 'S') &&
                        (header[1] == 'o' || header[1] == 'O') &&
                        (header[2] == 'l' || header[2] == 'L') &&
                        (header[3] == 'i' || header[3] == 'I') &&
                        (header[4] == 'd' || header[4] == 'D');

        // Call the appropriate loader
        bool status = false;
        if (is_ascii) {
            status = stl_load_ascii(filename, m, orig_vertex_indices);
        } else {
            status = stl_load_binary(filename, m, orig_vertex_indices);
        }

        if (!status) {
            cerr << "Could not load file " << filename << endl;
            return false;
        }

        // The current state of the mesh is not valid, so we need to clean it up
        // since it is just a triangle soup.
        double e_len = HMesh::average_edge_length(m);

        int missing_edges = stitch_mesh(m, 1e-6 * e_len);

        if (missing_edges > 0) {
            cerr << "Could not stitch " << missing_edges << " edges" << endl;
            return false;
        }

        m.cleanup();
        return true;
    }

    bool stl_load(const string& filename, Manifold& m) {
        VertexAttributeVector<int> orig_vertex_indices;
        return stl_load(filename, m, orig_vertex_indices);
    }

} // namespace HMesh
