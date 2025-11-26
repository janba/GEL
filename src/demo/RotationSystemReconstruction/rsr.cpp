#include <sstream>
#include <GEL/HMesh/RsR.h>
#include <GEL/HMesh/obj_save.h>

using Point = CGLA::Vec3d;
using Vector = CGLA::Vec3d;

struct PointCloud {
    std::vector<Point> vertices;
    std::vector<Vector> normals;
};

void read_obj(std::string file_path, PointCloud& pc) {
    std::ifstream file(file_path);
    std::string line;
    if (!file.is_open()) {
        throw std::runtime_error("Error: Unable to open file " + file_path);
    }
    while (std::getline(file, line))
    {
        std::vector<std::string> info;
        int pos = 0;
        while ((pos = line.find(" ")) != std::string::npos) {
            info.push_back(line.substr(0, pos));
            line.erase(0, pos + 1);
        }
        info.push_back(line);
        if (info.size() == 0) {
            continue;
        }
        if (info.at(0) == "v") {
            Point vertex(std::stof(info.at(1)),
                std::stof(info.at(2)), std::stof(info.at(3)));
            pc.vertices.push_back(vertex);
        }
        if (info.at(0) == "vn") {
            Vector normal(std::stof(info.at(1)),
                std::stof(info.at(2)), std::stof(info.at(3)));
            pc.normals.push_back(normal);
        }
    }
    file.close();
    return;
}

void read_ply(std::string file_path, PointCloud& pc) {
    std::ifstream file(file_path, std::ios::binary);

    if (!file.is_open()) {
        throw std::runtime_error("Error: Unable to open file " + file_path);
    }

    std::string line;
    bool headerEnded = false;
    size_t vertexCount = 0;
    bool hasNormals = false;
    bool isBinary = false;

    // Read the header
    while (std::getline(file, line)) {
        if (line.find("format ascii") != std::string::npos) {
            isBinary = false;
        }
        else if (line.find("format binary_little_endian") != std::string::npos) {
            isBinary = true;
        }
        else if (line.find("format binary_big_endian") != std::string::npos) {
            throw std::runtime_error("Error: Big-endian format not supported.");
        }
        else if (line.find("element vertex") != std::string::npos) {
            std::istringstream iss(line);
            std::string element, vertex;
            iss >> element >> vertex >> vertexCount;
        }
        else if (line.find("property float nx") != std::string::npos) {
            hasNormals = true;
        }
        else if (line == "end_header") {
            headerEnded = true;
            break;
        }
    }

    if (!headerEnded) {
        throw std::runtime_error("Error: PLY file does not have a valid header.");
    }

    // Read the body
    if (isBinary) {
        // Binary PLY reading
        for (size_t i = 0; i < vertexCount; ++i) {
            float x, y, z, nx = 0, ny = 0, nz = 0;

            // Read vertex coordinates
            file.read(reinterpret_cast<char*>(&x), sizeof(float));
            file.read(reinterpret_cast<char*>(&y), sizeof(float));
            file.read(reinterpret_cast<char*>(&z), sizeof(float));
            pc.vertices.push_back({ x, y, z });

            // Read normals if present
            if (hasNormals) {
                file.read(reinterpret_cast<char*>(&nx), sizeof(float));
                file.read(reinterpret_cast<char*>(&ny), sizeof(float));
                file.read(reinterpret_cast<char*>(&nz), sizeof(float));
                pc.normals.push_back({ nx, ny, nz });
            }
        }
    }
    else {
        // ASCII PLY reading
        size_t count = 0;
        while (std::getline(file, line) && count < vertexCount) {
            std::istringstream iss(line);
            float x, y, z, nx = 0, ny = 0, nz = 0;

            iss >> x >> y >> z;
            pc.vertices.push_back({ x, y, z });

            if (hasNormals) {
                iss >> nx >> ny >> nz;
                pc.normals.push_back({ nx, ny, nz });
            }

            ++count;
        }
    }

    file.close();
    return;
}

int main() {
    // {    
    //     PointCloud input;
    //     // Test on genus-0 shape
    //     read_obj("../../../data/PointClouds/owl-little.obj", input);
    //     HMesh::Manifold output;
    //     reconstruct_single(output, input.vertices,
    //     input.normals, false);
    //     HMesh::obj_save("owl-little-out.obj", output);
    // }

    // {
    //     PointCloud input;
    //     // Test on high-genus shape
    //     read_obj("../../../data/PointClouds/capital_A.obj", input);
    //     HMesh::Manifold output;
    //     reconstruct_single(output, input.vertices,
    //         input.normals, 
    //         true,   // Use Euclidean distance
    //         -1,     // Genus auto-detect
    //         30,     // Neighborhood size
    //         20,     // Radius for local operations
    //         60,     // Angle threshold in degrees
    //         40);    // Sample count

    //     HMesh::obj_save("capital_A-out.obj", output);
    // }

    // For new normal computation
    {
       PointCloud input;
       // Test on genus-0 shape
       //read_obj("../../../data/PointClouds/asn.obj", input);
       read_obj("../../../data/as.obj", input);
       std::cout << HMesh::DUMMY << std::endl;
       HMesh::Manifold output;
       HMesh::reconstruct_single(output, input.vertices, input.normals, true, -1, 10, 20, 60, 40);
       HMesh::obj_save("as-neighbor15-out.obj", output);
    }
	return 0;
}