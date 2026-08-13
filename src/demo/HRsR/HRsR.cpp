#include <sstream>
#include <GEL/HMesh/RsR.h>
#include <GEL/HMesh/HierarchicalReconstruction.h>
#include <GEL/HMesh/RSRExperimental.h>
#include <GEL/HMesh/obj_save.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Util/ThreadPool.h>

#include <filesystem>
#include <chrono>
#include <string_view>
#include <string>
#include <optional>
#include <numeric>
#include <cstring>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <fstream>
#include <iostream>

using Point = CGLA::Vec3d;
using Vector = CGLA::Vec3d;
namespace fs = std::filesystem;


struct PointCloud{
    std::vector<Point> vertices;
    std::vector<Vector> normals;
};

constexpr auto IS_EUCLIDEAN = false;
constexpr auto K_PARAM = 70;
constexpr auto GENUS = -1;
constexpr auto R_PARAM = 20;
constexpr auto THETA = 60;
constexpr auto N_PARAM = 20;


struct Property {
    std::string name;
    std::string type;
    size_t offset;   // for binary
    size_t size;
    bool isList = false;
};

size_t type_size(const std::string& type)
{
    if (type == "char" || type == "uchar" || type == "int8" || type == "uint8") return 1;
    if (type == "short" || type == "ushort" || type == "int16" || type == "uint16") return 2;
    if (type == "int" || type == "uint" || type == "int32" || type == "uint32" || type == "float") return 4;
    if (type == "double" || type == "float64") return 8;
    throw std::runtime_error("Unknown type: " + type);
}

auto test_options()
{
    HMesh::RSR::RSROpts opts;
    opts.dist = IS_EUCLIDEAN ? HMesh::RSR::Distance::Euclidean : HMesh::RSR::Distance::Tangent;
    opts.num_neighbors = K_PARAM;
    opts.genus = GENUS;
    opts.max_neighbor_dist = R_PARAM;
    opts.max_normal_ang = THETA;
    opts.max_handle_dist = N_PARAM;
    return opts;
}

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

void read_ply(const std::string& path, PointCloud& pc)
{
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open file");

    std::string line;
    bool isBinary = false;
    size_t vertexCount = 0;

    std::vector<Property> properties;
    bool inVertex = false;

    // -------- Parse header --------
    while (std::getline(file, line))
    {
        if (line.find("format ascii") != std::string::npos)
            isBinary = false;
        else if (line.find("format binary_little_endian") != std::string::npos)
            isBinary = true;
        else if (line.find("binary_big_endian") != std::string::npos)
            throw std::runtime_error("Big endian not supported");

        else if (line.rfind("element", 0) == 0)
        {
            std::istringstream iss(line);
            std::string tmp, name;
            iss >> tmp >> name;

            if (name == "vertex") {
                iss >> vertexCount;
                inVertex = true;
            }
            else {
                inVertex = false;
            }
        }
        else if (inVertex && line.rfind("property", 0) == 0)
        {
            std::istringstream iss(line);
            std::string tmp, type, name;
            iss >> tmp >> type;

            if (type == "list") {
                // Skip list property completely
                Property p;
                p.isList = true;
                properties.push_back(p);
                continue;
            }

            iss >> name;

            Property p;
            p.name = name;
            p.type = type;
            p.size = type_size(type);
            p.isList = false;
            properties.push_back(p);
        }
        else if (line.find("end_header") != std::string::npos)
        {
            break;
        }
    }

    if (vertexCount == 0)
        throw std::runtime_error("No vertices found");

    // -------- Precompute binary offsets --------
    size_t stride = 0;
    for (auto& p : properties)
    {
        if (!p.isList) {
            p.offset = stride;
            stride += p.size;
        }
    }

    // Map property indices
    std::unordered_map<std::string, int> propIndex;
    for (int i = 0; i < (int)properties.size(); ++i)
        if (!properties[i].isList)
            propIndex[properties[i].name] = i;

    auto get = [&](const std::string& name) -> int {
        return propIndex.count(name) ? propIndex[name] : -1;
        };

    int ix = get("x"), iy = get("y"), iz = get("z");
    int inx = get("nx"), iny = get("ny"), inz = get("nz");

    // -------- Read data --------
    std::vector<char> buffer(stride);

    for (size_t i = 0; i < vertexCount; ++i)
    {
        double x = 0, y = 0, z = 0;
        double nx = 0, ny = 0, nz = 0;

        if (isBinary)
        {
            file.read(buffer.data(), stride);

            auto read_from_buffer = [&](int idx) -> double {
                if (idx < 0) return 0.0;
                const Property& p = properties[idx];
                const char* ptr = buffer.data() + p.offset;

                if (p.type == "float") {
                    float v; std::memcpy(&v, ptr, 4); return v;
                }
                if (p.type == "double") {
                    double v; std::memcpy(&v, ptr, 8); return v;
                }
                if (p.type == "uchar" || p.type == "uint8") {
                    uint8_t v; std::memcpy(&v, ptr, 1); return v;
                }
                if (p.type == "int" || p.type == "int32") {
                    int32_t v; std::memcpy(&v, ptr, 4); return v;
                }

                return 0.0;
                };

            x = read_from_buffer(ix);
            y = read_from_buffer(iy);
            z = read_from_buffer(iz);

            nx = read_from_buffer(inx);
            ny = read_from_buffer(iny);
            nz = read_from_buffer(inz);
        }
        else
        {
            std::getline(file, line);
            std::istringstream iss(line);

            for (size_t j = 0; j < properties.size(); ++j)
            {
                if (properties[j].isList) {
                    int count;
                    iss >> count;
                    for (int k = 0; k < count; ++k) {
                        std::string tmp;
                        iss >> tmp;
                    }
                    continue;
                }

                std::string token;
                iss >> token;
                double v = std::stod(token);

                if ((int)j == ix) x = v;
                else if ((int)j == iy) y = v;
                else if ((int)j == iz) z = v;
                else if ((int)j == inx) nx = v;
                else if ((int)j == iny) ny = v;
                else if ((int)j == inz) nz = v;
            }
        }

        // Keep even if normals missing
        if (std::isfinite(x) && std::isfinite(y) && std::isfinite(z))
        {
            pc.vertices.emplace_back(x, y, z);

            if (std::isfinite(nx) && std::isfinite(ny) && std::isfinite(nz))
                pc.normals.emplace_back(nx, ny, nz);
        }
    }

    std::cout << "Loaded " << pc.vertices.size() << " vertices\n";
}

void read_pc_ply(fs::path file_path,
    PointCloud& pc_out)
{
    PointCloud pc;
    read_ply(file_path.string(), pc); // <-- use the robust reader we wrote

    const auto& vertex_pos = pc.vertices;
    const auto& pc_normals = pc.normals;

    // -------- index array --------
    std::vector<int> idx_array(vertex_pos.size());
    std::iota(idx_array.begin(), idx_array.end(), 0);

    // -------- optional downsampling --------
    bool downSample = true;
    if (downSample)
    {
        const int num_sample = 20000000;

        if (vertex_pos.size() > num_sample)
        {
            std::mt19937 g(12345);
            std::shuffle(idx_array.begin(), idx_array.end(), g);
            idx_array.resize(num_sample);
        }
    }

    // -------- main copy --------
    for (int idx : idx_array)
    {
        const auto& p = vertex_pos[idx];

        if (std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]))
        {
            pc_out.vertices.emplace_back(
                float(p[0]), float(p[1]), float(p[2])
            );

            // normals (only if available and aligned)
            if (pc_normals.size() == vertex_pos.size())
            {
                const auto& n = pc_normals[idx];

                if (std::isfinite(n[0]) &&
                    std::isfinite(n[1]) &&
                    std::isfinite(n[2]))
                {
                    pc_out.normals.emplace_back(
                        float(n[0]), float(n[1]), float(n[2])
                    );
                }
            }
        }
    }
}

auto reconstruct_serial(const std::string_view file_name, const HMesh::RSR::RSROpts& opts) {
    PointCloud input;
    fs::path path(file_name);

    // Get extension and normalize
    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == ".obj")
    {
        read_obj(path.string().c_str(), input);
    }
    else if (ext == ".ply")
    {
        read_pc_ply(path.string().c_str(), input);
    }
    else
    {
        throw std::runtime_error("Unsupported file format: " + ext);
    }

    HMesh::Manifold output;
    reconstruct_single(output,
        input.vertices,
        input.normals,
        opts.dist == HMesh::RSR::Distance::Euclidean, opts.genus, opts.num_neighbors,
        opts.max_neighbor_dist, opts.max_normal_ang, opts.max_handle_dist);
    std::cout << output.positions.size() << "\n";

    return output;
}

auto reconstruct_parallel(const std::string_view file_name, const HMesh::RSR::RSROpts& opts) -> std::optional<HMesh::Manifold> {
    std::cout << "======================\n"
        << "Begin new function\n";
    PointCloud input;
    fs::path path(file_name);

    // Get extension and normalize
    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == ".obj")
    {
        read_obj(path.string().c_str(), input);
    }
    else if (ext == ".ply")
    {
        read_pc_ply(path.string().c_str(), input);
    }
    else
    {
        throw std::runtime_error("Unsupported file format: " + ext);
    }

    std::cout << "obj vertices: " << input.vertices.size() << "\n";
    std::cout << "obj normals: " << input.normals.size() << "\n";

    HMesh::Manifold output;
    point_cloud_to_mesh(input.vertices, {}, opts, output);
    std::cout << output.positions.size() << "\n";

    return output;
}

auto reconstruct_collapse_reexpand(const std::string_view file_name, const HMesh::RSR::CollapseOpts& collapse_opts,
    const HMesh::RSR::RSROpts& rsr_opts, const HMesh::RSR::ReexpandOpts& reexpand) -> std::optional<HMesh::Manifold>
{
    std::cout << "======================\n"
        << "Begin new function\n";
    PointCloud input;
    fs::path path(file_name);

    // Get extension and normalize
    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == ".obj")
    {
        read_obj(path.string().c_str(), input);
    }
    else if (ext == ".ply")
    {
        read_pc_ply(path.string().c_str(), input);
    }
    else
    {
        throw std::runtime_error("Unsupported file format: " + ext);
    }

    std::cout << "obj vertices: " << input.vertices.size() << "\n";
    std::cout << "obj normals: " << input.normals.size() << "\n";

    std::cout << collapse_opts.max_collapses << std::endl;

    HMesh::Manifold output;
    std::cout << reexpand.refinement_iterations << std::endl;
    point_cloud_collapse_reexpand(input.vertices, input.normals, collapse_opts, rsr_opts,
        reexpand, output);
    std::cout << output.positions.size() << "\n";

    return output;
}

void iterate_reconstruct_serial(fs::path input_dir, fs::path output_dir) {
    const auto opts = test_options();
    {
        fs::create_directories(output_dir);

        // Open log file
        std::ofstream log_file(output_dir / "timing.txt");
        if (!log_file.is_open())
        {
            std::cerr << "Failed to open log file.\n";
            return;
        }

        for (const auto& entry : fs::directory_iterator(input_dir))
        {
            if (!entry.is_regular_file())
                continue;

            fs::path input_path = entry.path();

            if (input_path.extension() != ".obj" && input_path.extension() != ".ply")
                continue;

            std::cout << "Processing: " << input_path << std::endl;

            // ---- Start timing ----
            auto start = std::chrono::high_resolution_clock::now();

            std::optional<HMesh::Manifold> manifold =
                reconstruct_serial(input_path.string(), opts);

            // ---- End timing ----
            auto end = std::chrono::high_resolution_clock::now();
            double elapsed_sec =
                std::chrono::duration<double>(end - start).count();

            if (manifold.has_value())
            {
                fs::path output_path = output_dir / input_path.stem();
                HMesh::obj_save(output_path.string() + ".obj", *manifold);

                std::cout << "Saved: " << output_path << std::endl;

                log_file << input_path.filename().string()
                    << " SUCCESS "
                    << elapsed_sec << " sec\n";
            }
            else
            {
                std::cout << "Failed: " << input_path << std::endl;

                log_file << input_path.filename().string()
                    << " FAILED "
                    << elapsed_sec << " sec\n";
            }
        }

        log_file.close();
    }
}

void iterate_reconstruct_parallel(fs::path input_dir, fs::path output_dir) {
    const auto opts = test_options();
    {
        fs::create_directories(output_dir);

        // Open log file
        std::ofstream log_file(output_dir / "timing.txt");
        if (!log_file.is_open())
        {
            std::cerr << "Failed to open log file.\n";
            return;
        }

        for (const auto& entry : fs::directory_iterator(input_dir))
        {
            if (!entry.is_regular_file())
                continue;

            fs::path input_path = entry.path();

            if (input_path.extension() != ".obj" && input_path.extension() != ".ply")
                continue;

            std::cout << "Processing: " << input_path << std::endl;

            // ---- Start timing ----
            auto start = std::chrono::high_resolution_clock::now();

            std::optional<HMesh::Manifold> manifold = reconstruct_parallel(input_path.string(), opts);

            // ---- End timing ----
            auto end = std::chrono::high_resolution_clock::now();
            double elapsed_sec =
                std::chrono::duration<double>(end - start).count();

            if (manifold.has_value())
            {
                fs::path output_path = output_dir / input_path.stem();
                HMesh::obj_save(output_path.string() + ".obj", *manifold);

                std::cout << "Saved: " << output_path << std::endl;

                log_file << input_path.filename().string()
                    << " SUCCESS "
                    << elapsed_sec << " sec\n";
            }
            else
            {
                std::cout << "Failed: " << input_path << std::endl;

                log_file << input_path.filename().string()
                    << " FAILED "
                    << elapsed_sec << " sec\n";
            }
        }

        log_file.close();
    }
}

void iterate_reconstruct_collapse(fs::path input_dir, fs::path output_dir, int iteration) {
    const auto opts = test_options();
    {
        fs::create_directories(output_dir);

        // Open log file
        std::ofstream log_file(output_dir / "timing.txt");
        if (!log_file.is_open())
        {
            std::cerr << "Failed to open log file.\n";
            return;
        }

        for (const auto& entry : fs::directory_iterator(input_dir))
        {
            if (!entry.is_regular_file())
                continue;

            fs::path input_path = entry.path();

            if (input_path.extension() != ".obj" && input_path.extension() != ".ply")
                continue;

            std::cout << "Processing: " << input_path << std::endl;

            // ---- Start timing ----
            auto start = std::chrono::high_resolution_clock::now();

            HMesh::RSR::CollapseOpts collapseopts = HMesh::RSR::CollapseOpts();
            collapseopts.max_iterations = iteration;
            HMesh::RSR::ReexpandOpts expand_opts;
            expand_opts.refinement_iterations = 0;
            std::optional<HMesh::Manifold> manifold = reconstruct_collapse_reexpand(input_path.string(), collapseopts,
                opts, expand_opts);

            // ---- End timing ----
            auto end = std::chrono::high_resolution_clock::now();
            double elapsed_sec =
                std::chrono::duration<double>(end - start).count();

            if (manifold.has_value())
            {
                fs::path output_path = output_dir / input_path.stem();
                HMesh::obj_save(output_path.string() + ".obj", *manifold);

                std::cout << "Saved: " << output_path << std::endl;

                log_file << input_path.filename().string()
                    << " SUCCESS "
                    << elapsed_sec << " sec\n";
            }
            else
            {
                std::cout << "Failed: " << input_path << std::endl;

                log_file << input_path.filename().string()
                    << " FAILED "
                    << elapsed_sec << " sec\n";
            }
        }

        log_file.close();
    }
}



int main() {
    fs::path input_dir = "/home/ruicu/HRsR/data/WNNC";
    fs::path output_dir = "/home/ruicu/HRsR/result/HRsR/iter_3/WNNC";
    iterate_reconstruct_collapse(input_dir, output_dir, 3);

	return 0;
}