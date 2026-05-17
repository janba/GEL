#include <sstream>
#include <GEL/HMesh/RsR.h>
#include <GEL/HMesh/HierarchicalReconstruction.h>
#include <GEL/HMesh/RSRExperimental.h>
#include <GEL/HMesh/obj_save.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Util/ThreadPool.h>

#include <filesystem>
#include <chrono>
#include <ranges>
#include <array>
#include <string_view>
#include <sys/resource.h>
#include <sys/wait.h>

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

//static constexpr auto FILE_AS = "../../../../data/arma_comp_new.obj";
//static constexpr auto FILE_AS = "../../../../data/bunny_w_normal.obj";
//static constexpr auto FILE_AS = "../../../../data/point_on_square.obj";
//static constexpr auto FILE_AS = "../../../../data/one_ring.obj";
//static constexpr auto FILE_AS = "../../../../data/point_on_square.obj";
//static constexpr auto FILE_AS = "../../../../data/bun_complete.obj";
static constexpr auto FILE_AS = "C:/Users/ruicu/Desktop/SGP26/data/TanksTemple/Barn.ply";


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

double read_binary_value(std::ifstream& file, const std::string& type)
{
    if (type == "float") {
        float v; file.read(reinterpret_cast<char*>(&v), 4); return v;
    }
    if (type == "double") {
        double v; file.read(reinterpret_cast<char*>(&v), 8); return v;
    }
    if (type == "uchar" || type == "uint8") {
        uint8_t v; file.read(reinterpret_cast<char*>(&v), 1); return v;
    }
    if (type == "int" || type == "int32") {
        int32_t v; file.read(reinterpret_cast<char*>(&v), 4); return v;
    }
    if (type == "uint" || type == "uint32") {
        uint32_t v; file.read(reinterpret_cast<char*>(&v), 4); return v;
    }
    if (type == "short" || type == "int16") {
        int16_t v; file.read(reinterpret_cast<char*>(&v), 2); return v;
    }
    if (type == "ushort" || type == "uint16") {
        uint16_t v; file.read(reinterpret_cast<char*>(&v), 2); return v;
    }

    throw std::runtime_error("Unsupported binary type: " + type);
}

auto test_options()
{
    HMesh::RSR::RSROpts opts;
    opts.dist = IS_EUCLIDEAN ? HMesh::RSR::Distance::Euclidean : HMesh::RSR::Distance::Tangent;
    opts.k = K_PARAM;
    opts.genus = GENUS;
    opts.r = R_PARAM;
    opts.theta = THETA;
    opts.n = N_PARAM;
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

//void read_ply(std::string file_path, PointCloud& pc) {
//    std::ifstream file(file_path, std::ios::binary);
//
//    if (!file.is_open()) {
//        throw std::runtime_error("Error: Unable to open file " + file_path);
//    }
//
//    std::string line;
//    bool headerEnded = false;
//    size_t vertexCount = 0;
//    bool hasNormals = false;
//    bool isBinary = false;
//
//    // Read the header
//    while (std::getline(file, line)) {
//        if (line.find("format ascii") != std::string::npos) {
//            isBinary = false;
//        }
//        else if (line.find("format binary_little_endian") != std::string::npos) {
//            isBinary = true;
//        }
//        else if (line.find("format binary_big_endian") != std::string::npos) {
//            throw std::runtime_error("Error: Big-endian format not supported.");
//        }
//        else if (line.find("element vertex") != std::string::npos) {
//            std::istringstream iss(line);
//            std::string element, vertex;
//            iss >> element >> vertex >> vertexCount;
//        }
//        else if (line.find("property float nx") != std::string::npos) {
//            hasNormals = true;
//        }
//        else if (line == "end_header") {
//            headerEnded = true;
//            break;
//        }
//    }
//
//    if (!headerEnded) {
//        throw std::runtime_error("Error: PLY file does not have a valid header.");
//    }
//
//    // Read the body
//    if (isBinary) {
//        // Binary PLY reading
//        for (size_t i = 0; i < vertexCount; ++i) {
//            float x, y, z, nx = 0, ny = 0, nz = 0;
//
//            // Read vertex coordinates
//            file.read(reinterpret_cast<char*>(&x), sizeof(float));
//            file.read(reinterpret_cast<char*>(&y), sizeof(float));
//            file.read(reinterpret_cast<char*>(&z), sizeof(float));
//            pc.vertices.push_back({ x, y, z });
//
//            // Read normals if present
//            if (hasNormals) {
//                file.read(reinterpret_cast<char*>(&nx), sizeof(float));
//                file.read(reinterpret_cast<char*>(&ny), sizeof(float));
//                file.read(reinterpret_cast<char*>(&nz), sizeof(float));
//                pc.normals.push_back({ nx, ny, nz });
//            }
//        }
//    }
//    else {
//        // ASCII PLY reading
//        size_t count = 0;
//        while (std::getline(file, line) && count < vertexCount) {
//            std::istringstream iss(line);
//            float x, y, z, nx = 0, ny = 0, nz = 0;
//
//            iss >> x >> y >> z;
//            pc.vertices.push_back({ x, y, z });
//
//            if (hasNormals) {
//                iss >> nx >> ny >> nz;
//                pc.normals.push_back({ nx, ny, nz });
//            }
//
//            ++count;
//        }
//    }
//
//    file.close();
//    return;
//}

inline size_t ply_type_size(const std::string& t)
{
    if (t == "float" || t == "float32") return 4;
    if (t == "double" || t == "float64") return 8;
    if (t == "uchar" || t == "uint8") return 1;
    if (t == "char" || t == "int8") return 1;
    if (t == "ushort" || t == "uint16") return 2;
    if (t == "short" || t == "int16") return 2;
    if (t == "uint" || t == "uint32") return 4;
    if (t == "int" || t == "int32") return 4;
    throw std::runtime_error("Unsupported PLY type: " + t);
}

template<typename T>
T read_binary(std::ifstream& file)
{
    T v;
    file.read(reinterpret_cast<char*>(&v), sizeof(T));
    return v;
}

double read_property_binary(std::ifstream& file, const std::string& type)
{
    if (type == "float" || type == "float32") return read_binary<float>(file);
    if (type == "double" || type == "float64") return read_binary<double>(file);
    if (type == "uchar" || type == "uint8") return read_binary<uint8_t>(file);
    if (type == "char" || type == "int8") return read_binary<int8_t>(file);
    if (type == "ushort" || type == "uint16") return read_binary<uint16_t>(file);
    if (type == "short" || type == "int16") return read_binary<int16_t>(file);
    if (type == "uint" || type == "uint32") return read_binary<uint32_t>(file);
    if (type == "int" || type == "int32") return read_binary<int32_t>(file);

    throw std::runtime_error("Unsupported binary type: " + type);
}

double read_property_ascii(std::istringstream& iss, const std::string&)
{
    double v;
    iss >> v;
    return v;
}

bool is_valid(double v)
{
    return std::isfinite(v) && std::abs(v) < 1e10;
}

//void read_ply(const std::string& file_path, PointCloud& pc)
//{
//    std::ifstream file(file_path, std::ios::binary);
//    if (!file) throw std::runtime_error("Cannot open file");
//
//    std::string line;
//    bool isBinary = false;
//    size_t vertexCount = 0;
//
//    std::vector<Property> properties;
//    bool inVertexElement = false;
//
//    // -------- Parse header --------
//    while (std::getline(file, line))
//    {
//        if (line.find("format ascii") != std::string::npos)
//            isBinary = false;
//        else if (line.find("format binary_little_endian") != std::string::npos)
//            isBinary = true;
//        else if (line.find("binary_big_endian") != std::string::npos)
//            throw std::runtime_error("Big endian not supported");
//
//        else if (line.find("element vertex") != std::string::npos)
//        {
//            std::istringstream iss(line);
//            std::string tmp;
//            iss >> tmp >> tmp >> vertexCount;
//            inVertexElement = true;
//        }
//        else if (line.find("element") != std::string::npos)
//        {
//            inVertexElement = false;
//        }
//        else if (inVertexElement && line.find("property") != std::string::npos)
//        {
//            std::istringstream iss(line);
//            std::string tmp, type, name;
//            iss >> tmp >> type >> name;
//            properties.push_back({ type, name });
//        }
//        else if (line == "end_header")
//        {
//            break;
//        }
//    }
//
//    if (vertexCount == 0)
//        throw std::runtime_error("No vertices found");
//
//    // -------- Read body --------
//    for (size_t i = 0; i < vertexCount; ++i)
//    {
//        double x = 0, y = 0, z = 0;
//        double nx = 0, ny = 0, nz = 0;
//
//        if (isBinary)
//        {
//            for (const auto& prop : properties)
//            {
//                double v = read_property_binary(file, prop.type);
//
//                if (prop.name == "x") x = v;
//                else if (prop.name == "y") y = v;
//                else if (prop.name == "z") z = v;
//                else if (prop.name == "nx") nx = v;
//                else if (prop.name == "ny") ny = v;
//                else if (prop.name == "nz") nz = v;
//            }
//        }
//        else
//        {
//            std::getline(file, line);
//            std::istringstream iss(line);
//
//            for (const auto& prop : properties)
//            {
//                double v = read_property_ascii(iss, prop.type);
//
//                if (prop.name == "x") x = v;
//                else if (prop.name == "y") y = v;
//                else if (prop.name == "z") z = v;
//                else if (prop.name == "nx") nx = v;
//                else if (prop.name == "ny") ny = v;
//                else if (prop.name == "nz") nz = v;
//            }
//        }
//
//        // -------- Sanity filter (critical) --------
//        if (!is_valid(x) || !is_valid(y) || !is_valid(z))
//            continue;
//
//        pc.vertices.emplace_back(x, y, z);
//
//        if (is_valid(nx) && is_valid(ny) && is_valid(nz))
//            pc.normals.emplace_back(nx, ny, nz);
//    }
//}

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
        //const int num_sample = 1000000;

        if (vertex_pos.size() > num_sample)
        {
            std::random_device rd;
            //std::mt19937 g(rd());
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
        opts.dist == HMesh::RSR::Distance::Euclidean, opts.genus, opts.k,
        opts.r, opts.theta, opts.n);
    // k: 70 is too large
    // r: needs isEuclidean false
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
    // k: 70 is too large
    // r: needs isEuclidean false
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
    // k: 70 is too large
    // r: needs isEuclidean false
    std::cout << output.positions.size() << "\n";

    return output;
}

void iterate_reconstruct_serial_linux(fs::path input_dir, fs::path output_dir) {
    const auto opts = test_options();

    fs::create_directories(output_dir);

    std::ofstream log_file(output_dir / "timing.txt");
    if (!log_file.is_open())
    {
        std::cerr << "Failed to open log file.\n";
        return;
    }
    log_file.close();

    for (const auto& entry : fs::directory_iterator(input_dir))
    {
        if (!entry.is_regular_file())
            continue;

        fs::path input_path = entry.path();

        if (input_path.extension() != ".obj" && input_path.extension() != ".ply")
            continue;
        log_file.open(output_dir / "timing.txt", std::ios::app);
        if (!log_file.is_open())
        {
            std::cerr << "Failed to open log file.\n";
            return;
        }
        std::cout << "Processing: " << input_path << std::endl;

        // ---- Start timing (parent) ----
        auto start = std::chrono::high_resolution_clock::now();

        pid_t pid = fork();

        if (pid == 0)
        {
            // ===== CHILD PROCESS =====
            std::optional<HMesh::Manifold> manifold =
                reconstruct_serial(input_path.string(), opts);

            if (manifold.has_value())
            {
                fs::path output_path = output_dir / input_path.stem();
                HMesh::obj_save(output_path.string() + ".obj", *manifold);
                exit(0); // success
            }
            else
            {
                exit(1); // failure
            }
        }
        else if (pid > 0)

        {
            // ===== PARENT PROCESS =====
            int status;
            struct rusage usage;

            wait4(pid, &status, 0, &usage);

            auto end = std::chrono::high_resolution_clock::now();
            double elapsed_sec =
                std::chrono::duration<double>(end - start).count();

            // ru_maxrss is in KB on Linux
            double peak_mb = usage.ru_maxrss / 1024.0;

            bool success = WIFEXITED(status) && (WEXITSTATUS(status) == 0);

            if (success)
            {
                std::cout << "Saved: " << input_path.stem() << std::endl;

                log_file << input_path.filename().string()
                    << " SUCCESS "
                    << elapsed_sec << " sec "
                    << peak_mb << " MB\n";
            }
            else
            {
                std::string fail_reason;

                if (WIFEXITED(status))
                {
                    int code = WEXITSTATUS(status);
                    fail_reason = "EXIT_CODE=" + std::to_string(code);
                }
                else if (WIFSIGNALED(status))
                {
                    int sig = WTERMSIG(status);
                    fail_reason = "SIGNALED=" + std::to_string(sig);

                    // Optional: human-readable signal
                    fail_reason += " (" + std::string(strsignal(sig)) + ")";
                }
                else
                {
                    fail_reason = "UNKNOWN_FAILURE";
                }

                std::cout << "Failed: " << input_path << " | " << fail_reason << std::endl;

                log_file << input_path.filename().string()
                    << " FAILED "
                    << fail_reason << " "
                    << elapsed_sec << " sec "
                    << peak_mb << " MB\n";
            }
        }
        else
        {
            std::cerr << "Fork failed!\n";
        }
        log_file.close();
    }

    //log_file.close();
}

void iterate_reconstruct_parallel_linux(fs::path input_dir, fs::path output_dir) {
    const auto opts = test_options();

    fs::create_directories(output_dir);

    std::ofstream log_file(output_dir / "timing.txt");
    if (!log_file.is_open())
    {
        std::cerr << "Failed to open log file.\n";
        return;
    }
    log_file.close();

    for (const auto& entry : fs::directory_iterator(input_dir))
    {
        if (!entry.is_regular_file())
            continue;

        fs::path input_path = entry.path();

        if (input_path.extension() != ".obj" && input_path.extension() != ".ply")
            continue;
        log_file.open(output_dir / "timing.txt", std::ios::app);
        if (!log_file.is_open())
        {
            std::cerr << "Failed to open log file.\n";
            return;
        }
        std::cout << "Processing: " << input_path << std::endl;

        // ---- Start timing (parent) ----
        auto start = std::chrono::high_resolution_clock::now();

        pid_t pid = fork();

        if (pid == 0)
        {
            // ===== CHILD PROCESS =====
            std::optional<HMesh::Manifold> manifold = reconstruct_parallel(input_path.string(), opts);

            if (manifold.has_value())
            {
                fs::path output_path = output_dir / input_path.stem();
                HMesh::obj_save(output_path.string() + ".obj", *manifold);
                exit(0); // success
            }
            else
            {
                exit(1); // failure
            }
        }
        else if (pid > 0)
        {
            // ===== PARENT PROCESS =====
            int status;
            struct rusage usage;

            wait4(pid, &status, 0, &usage);

            auto end = std::chrono::high_resolution_clock::now();
            double elapsed_sec =
                std::chrono::duration<double>(end - start).count();

            // ru_maxrss is in KB on Linux
            double peak_mb = usage.ru_maxrss / 1024.0;

            bool success = WIFEXITED(status) && (WEXITSTATUS(status) == 0);

            if (success)
            {
                std::cout << "Saved: " << input_path.stem() << std::endl;

                log_file << input_path.filename().string()
                    << " SUCCESS "
                    << elapsed_sec << " sec "
                    << peak_mb << " MB\n";
            }
            else
            {
                std::string fail_reason;

                if (WIFEXITED(status))
                {
                    int code = WEXITSTATUS(status);
                    fail_reason = "EXIT_CODE=" + std::to_string(code);
                }
                else if (WIFSIGNALED(status))
                {
                    int sig = WTERMSIG(status);
                    fail_reason = "SIGNALED=" + std::to_string(sig);

                    // Optional: human-readable signal
                    fail_reason += " (" + std::string(strsignal(sig)) + ")";
                }
                else
                {
                    fail_reason = "UNKNOWN_FAILURE";
                }

                std::cout << "Failed: " << input_path << " | " << fail_reason << std::endl;

                log_file << input_path.filename().string()
                    << " FAILED "
                    << fail_reason << " "
                    << elapsed_sec << " sec "
                    << peak_mb << " MB\n";
            }
        }
        else
        {
            std::cerr << "Fork failed!\n";
        }
        log_file.close();
    }

    //log_file.close();
}

void iterate_reconstruct_collapse_linux(fs::path input_dir, fs::path output_dir, int iteration) {
    const auto opts = test_options();

    fs::create_directories(output_dir);

    std::ofstream log_file(output_dir / "timing.txt");
    if (!log_file.is_open())
    {
        std::cerr << "Failed to open log file.\n";
        return;
    }
    log_file.close();

    for (const auto& entry : fs::directory_iterator(input_dir))
    {
        if (!entry.is_regular_file())
            continue;

        fs::path input_path = entry.path();

        if (input_path.extension() != ".obj" && input_path.extension() != ".ply")
            continue;
        log_file.open(output_dir / "timing.txt", std::ios::app);
        if (!log_file.is_open())
        {
            std::cerr << "Failed to open log file.\n";
            return;
        }
        std::cout << "Processing: " << input_path << std::endl;

        // ---- Start timing (parent) ----
        auto start = std::chrono::high_resolution_clock::now();

        pid_t pid = fork();

        if (pid == 0)
        {
            // ===== CHILD PROCESS =====
            HMesh::RSR::CollapseOpts collapseopts = HMesh::RSR::CollapseOpts();
            collapseopts.max_iterations = iteration;
            //collapseopts.reduction_per_iteration = 0.7;
            //collapseopts.initial_neighbors = 70;
            collapseopts.distance = HMesh::RSR::Distance::Tangent;
            HMesh::RSR::ReexpandOpts expand_opts;
            //expand_opts.refinement_iterations = 0;
            std::optional<HMesh::Manifold> manifold = reconstruct_collapse_reexpand(input_path.string(), collapseopts,
                opts, expand_opts);

            if (manifold.has_value())
            {
                fs::path output_path = output_dir / input_path.stem();
                HMesh::obj_save(output_path.string() + ".obj", *manifold);
                exit(0); // success
            }
            else
            {
                exit(1); // failure
            }
        }
        else if (pid > 0)

        {
            // ===== PARENT PROCESS =====
            int status;
            struct rusage usage;

            wait4(pid, &status, 0, &usage);

            auto end = std::chrono::high_resolution_clock::now();
            double elapsed_sec =
                std::chrono::duration<double>(end - start).count();

            // ru_maxrss is in KB on Linux
            double peak_mb = usage.ru_maxrss / 1024.0;

            bool success = WIFEXITED(status) && (WEXITSTATUS(status) == 0);

            if (success)
            {
                std::cout << "Saved: " << input_path.stem() << std::endl;

                log_file << input_path.filename().string()
                    << " SUCCESS "
                    << elapsed_sec << " sec "
                    << peak_mb << " MB\n";
            }
            else
            {
                std::string fail_reason;

                if (WIFEXITED(status))
                {
                    int code = WEXITSTATUS(status);
                    fail_reason = "EXIT_CODE=" + std::to_string(code);
                }
                else if (WIFSIGNALED(status))
                {
                    int sig = WTERMSIG(status);
                    fail_reason = "SIGNALED=" + std::to_string(sig);

                    // Optional: human-readable signal
                    fail_reason += " (" + std::string(strsignal(sig)) + ")";
                }
                else
                {
                    fail_reason = "UNKNOWN_FAILURE";
                }

                std::cout << "Failed: " << input_path << " | " << fail_reason << std::endl;

                log_file << input_path.filename().string()
                    << " FAILED "
                    << fail_reason << " "
                    << elapsed_sec << " sec "
                    << peak_mb << " MB\n";
            }
        }
        else
        {
            std::cerr << "Fork failed!\n";
        }
        log_file.close();
    }

    //log_file.close();
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
            //std::cout << "enter loop" << std::endl;
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

            //std::optional<HMesh::Manifold> manifold = reconstruct_parallel(input_path.string(), opts);


            //HMesh::RSR::CollapseOpts collapseopts = HMesh::RSR::CollapseOpts();
            //collapseopts.max_iterations = 1;
            //std::optional<HMesh::Manifold> manifold = reconstruct_collapse_reexpand(input_path.string(), collapseopts,
                //opts, HMesh::RSR::ReexpandOpts());

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
            //std::cout << "enter loop" << std::endl;
            if (!entry.is_regular_file())
                continue;

            fs::path input_path = entry.path();

            if (input_path.extension() != ".obj" && input_path.extension() != ".ply")
                continue;

            std::cout << "Processing: " << input_path << std::endl;

            // ---- Start timing ----
            auto start = std::chrono::high_resolution_clock::now();

            std::optional<HMesh::Manifold> manifold = reconstruct_parallel(input_path.string(), opts);

            //HMesh::RSR::CollapseOpts collapseopts = HMesh::RSR::CollapseOpts();
            //collapseopts.max_iterations = 1;
            //std::optional<HMesh::Manifold> manifold = reconstruct_collapse_reexpand(input_path.string(), collapseopts,
                //opts, HMesh::RSR::ReexpandOpts());

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

void iterate_reconstruct_collapse(fs::path input_dir, fs::path output_dir) {
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
            //std::cout << "enter loop" << std::endl;
            if (!entry.is_regular_file())
                continue;

            fs::path input_path = entry.path();

            if (input_path.extension() != ".obj" && input_path.extension() != ".ply")
                continue;

            std::cout << "Processing: " << input_path << std::endl;

            // ---- Start timing ----
            auto start = std::chrono::high_resolution_clock::now();

            HMesh::RSR::CollapseOpts collapseopts = HMesh::RSR::CollapseOpts();
            collapseopts.max_iterations = 1;
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
    /*fs::path input_dir = "C:\\Users\\ruicu\\Desktop\\SGP26\\data\\TanksTemple";
    fs::path output_dir = "C:\\Users\\ruicu\\Desktop\\SGP26\\Results\\RsR\\Tanks";*/
    if(true){
        Util::detail::ImmediatePool pool;
        std::cout << pool.size() << std::endl;

        //fs::path input_dir = "/home/ruicu/HRsR/data/TanksTemple";
        //fs::path input_dir = "/home/ruicu/HRsR/data/DTU";
        fs::path input_dir = "/home/ruicu/HRsR/data/WNNC";
        //fs::path input_dir = "/home/ruicu/HRsR/data/Stanford";
        //fs::path input_dir = "/home/ruicu/HRsR/data/topology";
        fs::path output_dir = "/home/ruicu/HRsR/result/RsR/Stanford";

        //iterate_reconstruct_serial_linux(input_dir, output_dir);

        output_dir = fs::path("/home/ruicu/HRsR/result/PRsR/Stanford");

        //iterate_reconstruct_parallel_linux(input_dir, output_dir);

        output_dir = fs::path("/home/ruicu/HRsR/result/HRsR/iter_1/topology");

        //iterate_reconstruct_collapse_linux(input_dir, output_dir, 1);

        output_dir = fs::path("/home/ruicu/HRsR/result/HRsR/iter_2/Stanford");

        //iterate_reconstruct_collapse_linux(input_dir, output_dir, 2);

        output_dir = fs::path("/home/ruicu/HRsR/result/HRsR/iter_3/WNNC");
        iterate_reconstruct_collapse_linux(input_dir, output_dir, 3);

    }
    
    /*Serial*/
   /* std::optional<HMesh::Manifold> manifold = reconstruct_serial(FILE_AS, opts);
    if (manifold.has_value())
        HMesh::obj_save("bun_complete-out_s.obj", *manifold);*/

    /*Parallel*/
    /*const auto opts = test_options();
    auto manifold = reconstruct_parallel(FILE_AS, opts);
    if (manifold.has_value())
        HMesh::obj_save("Barn.obj", *manifold);*/

    /*collapse*/
    //std::optional<HMesh::Manifold> manifold = reconstruct_parallel(FILE_AS, opts);
    /*const auto opts = test_options();
    HMesh::RSR::CollapseOpts collapseopts = HMesh::RSR::CollapseOpts();
    collapseopts.max_iterations = 1;
    collapseopts.distance = HMesh::RSR::Distance::Tangent;
    std::optional<HMesh::Manifold> manifold = reconstruct_collapse_reexpand(FILE_AS, collapseopts,
        opts, HMesh::RSR::ReexpandOpts());

    if(manifold.has_value())
        HMesh::obj_save("Barn.obj", *manifold);*/
        //HMesh::obj_save("bun_complete-out.obj", *manifold);

   /* {
        HMesh::Manifold one_ring;
        HMesh::obj_load(FILE_AS, one_ring);

        HMesh::VertexID center_id = HMesh::VertexID(0);
        std::cout << one_ring.pos(center_id) << std::endl;
        HMesh::HalfEdgeID h1, h2;
        for (const auto he : one_ring.incident_halfedges(center_id)) {
            const auto& walker = one_ring.walker(he);
            if (walker.vertex() == center_id) {
                h1 = he;
                h2 = walker.next().halfedge();
                std::cout << "h1: " << walker.opp().vertex() << " and " << walker.vertex() << std::endl;
                std::cout << "h2: " << walker.next().opp().vertex() << " and " << walker.next().vertex() << std::endl;
                break;
            }
            if (walker.opp().vertex() == center_id) {
                h1 = walker.opp().halfedge();
                h2 = walker.opp().next().halfedge();
                std::cout << "h1: " << walker.opp().vertex() << " and " << walker.vertex() << std::endl;
                std::cout << "h2: " << walker.next().opp().vertex() << " and " << walker.next().vertex() << std::endl;
                break;
            }
        }
        auto V_new = one_ring.slit_one_ring(h1, h2);


        one_ring.positions[V_new] = CGLA::Vec3d{ 0.5, 0.5, 0.5 };

        HMesh::obj_save("one_ring_before_closing_hole.obj", one_ring);

        const auto& walker_vn = one_ring.walker(V_new);
        HMesh::HalfEdgeID h_id = walker_vn.halfedge();
        const auto f = one_ring.close_hole(one_ring.kernel.out(V_new));
        std::cout << one_ring.walker(one_ring.kernel.out(V_new)).vertex() << " and "
            << one_ring.walker(one_ring.kernel.out(V_new)).opp().vertex() << std::endl;
        one_ring.split_face_by_edge(f, center_id, V_new);

        std::cout << "Sanity check" << std::endl;
        {
            const auto& walker = one_ring.walker(center_id);
            auto he = walker.halfedge();
            HMesh::HalfEdgeID h = he;
            h = walker.circulate_vertex_ccw().halfedge();
            while (h != he) {
                const auto& h_walker = one_ring.walker(h);
                std::cout << h_walker.opp().vertex() << " and " << h_walker.vertex() << std::endl;
                h = h_walker.circulate_vertex_ccw().halfedge();
            }
        }

        HMesh::obj_save("one_ring.obj", one_ring);
    }
    */

	return 0;
}