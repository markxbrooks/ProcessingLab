#include "mosflm/io/Json.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <system_error>

#include <fcntl.h>
#include <unistd.h>

namespace mosflm::io {
namespace {

void escape(std::string& out, const std::string& s) {
    out += '"';
    for (const char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof buf, "\\u%04x", static_cast<unsigned char>(c));
                    out += buf;
                } else {
                    out += c;  // UTF-8 passes through
                }
        }
    }
    out += '"';
}

const Json::Array kEmptyArray;
const Json::Object kEmptyObject;
const std::string kEmptyString;

}  // namespace

Json& Json::operator[](const std::string& key) {
    if (is_null()) {
        v_ = Object{};
    }
    auto& obj = std::get<Object>(v_);
    for (auto& [k, v] : obj) {
        if (k == key) {
            return v;
        }
    }
    obj.emplace_back(key, Json());
    return obj.back().second;
}

const Json* Json::find(const std::string& key) const {
    if (!is_object()) {
        return nullptr;
    }
    for (const auto& [k, v] : std::get<Object>(v_)) {
        if (k == key) {
            return &v;
        }
    }
    return nullptr;
}

void Json::push_back(Json v) {
    if (is_null()) {
        v_ = Array{};
    }
    std::get<Array>(v_).push_back(std::move(v));
}

const Json::Array& Json::items() const {
    return is_array() ? std::get<Array>(v_) : kEmptyArray;
}

Json::Array& Json::elements() {
    if (is_null()) {
        v_ = Array{};
    }
    return std::get<Array>(v_);
}

const Json::Object& Json::members() const {
    return is_object() ? std::get<Object>(v_) : kEmptyObject;
}

double Json::number() const {
    if (const auto* d = std::get_if<double>(&v_)) {
        return *d;
    }
    if (const auto* i = std::get_if<std::int64_t>(&v_)) {
        return static_cast<double>(*i);
    }
    return 0.0;
}

const std::string& Json::string() const {
    const auto* s = std::get_if<std::string>(&v_);
    return s != nullptr ? *s : kEmptyString;
}

bool Json::boolean() const {
    const auto* b = std::get_if<bool>(&v_);
    return b != nullptr && *b;
}

std::string Json::dump(int indent) const {
    std::string out;
    dump_to(out, indent, 0);
    out += '\n';
    return out;
}

void Json::dump_to(std::string& out, int indent, int depth) const {
    const auto newline = [&](int d) {
        out += '\n';
        out.append(static_cast<std::size_t>(indent * d), ' ');
    };
    if (std::holds_alternative<std::monostate>(v_)) {
        out += "null";
    } else if (const auto* b = std::get_if<bool>(&v_)) {
        out += *b ? "true" : "false";
    } else if (const auto* i = std::get_if<std::int64_t>(&v_)) {
        out += std::to_string(*i);
    } else if (const auto* d = std::get_if<double>(&v_)) {
        if (!std::isfinite(*d)) {
            out += "null";
        } else {
            char buf[32];
            std::snprintf(buf, sizeof buf, "%.10g", *d);
            out += buf;
        }
    } else if (const auto* s = std::get_if<std::string>(&v_)) {
        escape(out, *s);
    } else if (const auto* a = std::get_if<Array>(&v_)) {
        if (a->empty()) {
            out += "[]";
            return;
        }
        // Arrays of plain numbers stay on one line (matrices, cells).
        bool flat = true;
        for (const auto& e : *a) {
            flat = flat && e.is_number();
        }
        out += '[';
        for (std::size_t k = 0; k < a->size(); ++k) {
            if (k > 0) {
                out += flat ? ", " : ",";
            }
            if (!flat) {
                newline(depth + 1);
            }
            (*a)[k].dump_to(out, indent, depth + 1);
        }
        if (!flat) {
            newline(depth);
        }
        out += ']';
    } else if (const auto* o = std::get_if<Object>(&v_)) {
        if (o->empty()) {
            out += "{}";
            return;
        }
        out += '{';
        for (std::size_t k = 0; k < o->size(); ++k) {
            if (k > 0) {
                out += ',';
            }
            newline(depth + 1);
            escape(out, (*o)[k].first);
            out += ": ";
            (*o)[k].second.dump_to(out, indent, depth + 1);
        }
        newline(depth);
        out += '}';
    }
}

void write_file_atomic(const std::filesystem::path& path, const std::string& text) {
    const std::filesystem::path tmp = path.string() + ".tmp";
    const int fd = ::open(tmp.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) {
        throw std::runtime_error("cannot write " + tmp.string());
    }
    const char* p = text.data();
    std::size_t n = text.size();
    while (n > 0) {
        const ssize_t w = ::write(fd, p, n);
        if (w <= 0) {
            ::close(fd);
            throw std::runtime_error("write failed: " + tmp.string());
        }
        p += w;
        n -= static_cast<std::size_t>(w);
    }
    ::fsync(fd);
    ::close(fd);
    std::error_code ec;
    std::filesystem::rename(tmp, path, ec);  // atomic on POSIX
    if (ec) {
        throw std::runtime_error("cannot rename " + tmp.string() + ": " + ec.message());
    }
}

}  // namespace mosflm::io
