//
// Created by Cem Akarsubasi on 5/22/25.
//

#include <string> // MSVC: getline, stod
#include <charconv>
#include <fstream>

#include <GEL/Util/RawObj.h>

#include <GEL/Util/Assert.h>

namespace Util
{
using namespace Combinators;

std::ostream& operator<<(std::ostream& os, const RawObj& obj)
{
    for (auto const& v : obj.vertices) {
        os << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
    }
    for (auto const& vn : obj.normals) {
        os << "vn" << vn[0] << " " << vn[1] << " " << vn[2] << "\n";
    }
    for (auto const& vt : obj.texture_coordinates) {
        os << "vt" << vt[0] << " " << vt[1] << "\n";
    }
    return os;
}

std::istream& operator>>(std::istream& is, RawObj& obj)
{
    obj.vertices.clear();
    std::string line;
    while (auto& more = std::getline(is, line)) {
        std::string_view line_view = line;
        if (auto pos = parse_prefix_then_float_triplet("v", line_view)) {
            obj.vertices.emplace_back(*pos);
        }
        if (auto pos = parse_prefix_then_float_triplet("vn", line_view)) {
            obj.normals.emplace_back(*pos);
        }
        if (auto pos = parse_prefix_then_float_triplet("vt", line_view)) {
            obj.normals.emplace_back(*pos);
        }
    }
    return is;
}

namespace Combinators
{
    bool parse_one(const char c, std::string_view& s)
    {
        if (s.empty()) {
            return false;
        } else if (s.front() == c) {
            s = s.substr(1);
            return true;
        }
        return false;
    }

    bool parse_string(const std::string_view& fragment, std::string_view& s)
    {
        if (s.starts_with(fragment)) {
            s = s.substr(fragment.length());
            return true;
        } else {
            return false;
        }
    }

    bool ignore_spaces(std::string_view& s)
    {
        const auto pos = s.find_first_not_of(' ');
        if (pos != std::string::npos) {
            s = s.substr(pos);
        } else {
            s = s.substr(s.length());
        }
        return true;
    }

    /// In a far away not-so-dark future, we will eventually have a floating point parser in the standard library
    /// With sane behavior (no ignoring whitespace or accepting leading + signs),
    /// Proper error reporting (no errno or exceptions)
    /// Supported by all compilers.
    ///
    /// Today is not that day, and as it turns out, while from_chars support can be tested,
    /// some compilers can claim to support it while not supporting it for floating point values.

#if !defined(__clang__) && !defined(_MSC_VER)
    std::optional<double> parse_float(std::string_view& s)
    {
        double out;
        const auto begin = s.begin();
        const auto end_pos = s.find(' ');
        const auto end = (end_pos == std::string::npos) ? s.end() : s.begin() + end_pos;
        const auto [p, ec] = std::from_chars(begin, end, out);
        if (ec != std::errc()) {
            return std::nullopt;
        } else {
            const auto offset = p - begin;
            s = s.substr(offset);
            return out;
        }
    }
    /// FIXME: Check if we live in that future where this code can be deleted
#else
    std::optional<double> parse_float(std::string_view& s) {
        const auto begin = s.begin();
        const auto end_pos = s.find(' ');
        if (end_pos == 0) {
            return std::nullopt;
        }
        const auto end = (end_pos == std::string::npos) ? s.end() : s.begin() + end_pos;
        // we need to make a copy because stod and strtod need a null terminated string
        // and a string view is not null terminated
        const auto copy = std::string(begin, end);
        try {
            size_t end_point = 0;
            double out = std::stod(copy, &end_point);
            if (end_point == 0) {
                return std::nullopt;
            }
            s = s.substr(end_point);
            return out;
        } catch (std::exception& ex) {
            return std::nullopt;
        }
    }
#endif

    std::optional<std::uint64_t> parse_uint(std::string_view& s)
    {
        throw std::runtime_error("unimplemented");
        return std::nullopt;
    }

    std::optional<double> parse_float_ws(std::string_view& s)
    {
        if (const auto d = parse_float(s)) {
            ignore_spaces(s);
            return d;
        } else {
            return std::nullopt;
        }
    }

    std::optional<CGLA::Vec2d> parse_float_doublet_ws(std::string_view& s)
    {
        auto s_temp = s;
        const auto out1 = parse_float_ws(s_temp);
        if (!out1) return std::nullopt;
        const auto out2 = parse_float_ws(s_temp);
        if (!out2) return std::nullopt;
        s = s_temp;
        return CGLA::Vec2d{*out1, *out2};
    }

    std::optional<CGLA::Vec3d> parse_float_triplet_ws(std::string_view& s)
    {
        auto s_temp = s;
        const auto out1 = parse_float_ws(s_temp);
        if (!out1) return std::nullopt;
        const auto out2 = parse_float_ws(s_temp);
        if (!out2) return std::nullopt;
        const auto out3 = parse_float_ws(s_temp);
        if (!out3) return std::nullopt;
        s = s_temp;
        return CGLA::Vec3d{*out1, *out2, *out3};
    }

    std::optional<CGLA::Vec2d> parse_prefix_then_float_doublet(const std::string_view prefix, std::string_view& s)
    {
        auto s_temp = s;
        if (!parse_string(prefix, s_temp)) return std::nullopt;
        ignore_spaces(s_temp);
        auto pos = parse_float_doublet_ws(s_temp);
        if (pos) {
            s = s_temp;
            return pos;
        } else {
            return std::nullopt;
        }
    }

    std::optional<CGLA::Vec3d> parse_prefix_then_float_triplet(const std::string_view prefix, std::string_view& s)
    {
        auto s_temp = s;
        if (!parse_string(prefix, s_temp)) return std::nullopt;
        ignore_spaces(s_temp);
        auto pos = parse_float_triplet_ws(s_temp);
        if (pos) {
            s = s_temp;
            return pos;
        } else {
            return std::nullopt;
        }
    }

    std::optional<RawObj::FaceElement> parse_face_element(std::string_view& s)
    {
        auto s_temp = s;
        if (!parse_string("f", s_temp)) return std::nullopt;
        ignore_spaces(s_temp);
        // TODO

        return std::nullopt;
    }
}

std::optional<RawObj> read_raw_obj(const std::filesystem::path& file_path)
{
    std::ifstream file(file_path);
    if (!file.is_open()) return std::nullopt;
    RawObj obj;
    file >> obj;
    if (file.bad()) return std::nullopt;
    return obj;
}


void write_raw_obj(const std::filesystem::path& file_path, const RawObj& obj)
{
    // TODO: error reporting
    std::ofstream file(file_path);
    file << obj << std::endl;
    file.close();
}
}
