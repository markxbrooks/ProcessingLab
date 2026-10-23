#include "mosflm/io/Json.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string_view>
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

/// Recursive-descent reader for RFC 8259 JSON.
class Parser {
public:
    explicit Parser(std::string_view text) : s_(text) {}

    Json document() {
        Json v = value();
        skip_ws();
        if (p_ != s_.size()) {
            fail("trailing characters");
        }
        return v;
    }

private:
    [[noreturn]] void fail(const std::string& what) const {
        throw std::runtime_error("JSON: " + what + " at offset " + std::to_string(p_));
    }

    void skip_ws() {
        while (p_ < s_.size() && (s_[p_] == ' ' || s_[p_] == '\t' || s_[p_] == '\n' || s_[p_] == '\r')) {
            ++p_;
        }
    }

    bool consume(std::string_view word) {
        if (s_.substr(p_, word.size()) == word) {
            p_ += word.size();
            return true;
        }
        return false;
    }

    void expect(char c) {
        skip_ws();
        if (p_ >= s_.size() || s_[p_] != c) {
            fail(std::string("expected '") + c + "'");
        }
        ++p_;
    }

    Json value() {
        skip_ws();
        if (p_ >= s_.size()) {
            fail("unexpected end");
        }
        const char c = s_[p_];
        if (c == '{') return object();
        if (c == '[') return array();
        if (c == '"') return Json(string());
        if (consume("true")) return Json(true);
        if (consume("false")) return Json(false);
        if (consume("null")) return Json();
        return number();
    }

    Json object() {
        ++p_;
        Json o = Json::object();
        skip_ws();
        if (p_ < s_.size() && s_[p_] == '}') {
            ++p_;
            return o;
        }
        for (;;) {
            skip_ws();
            if (p_ >= s_.size() || s_[p_] != '"') {
                fail("expected a member name");
            }
            std::string key = string();
            expect(':');
            o[key] = value();
            skip_ws();
            if (p_ < s_.size() && s_[p_] == ',') {
                ++p_;
                continue;
            }
            expect('}');
            return o;
        }
    }

    Json array() {
        ++p_;
        Json a = Json::array();
        skip_ws();
        if (p_ < s_.size() && s_[p_] == ']') {
            ++p_;
            return a;
        }
        for (;;) {
            a.push_back(value());
            skip_ws();
            if (p_ < s_.size() && s_[p_] == ',') {
                ++p_;
                continue;
            }
            expect(']');
            return a;
        }
    }

    unsigned hex4() {
        if (p_ + 4 > s_.size()) {
            fail("short \\u escape");
        }
        unsigned v = 0;
        for (int k = 0; k < 4; ++k) {
            const char c = s_[p_++];
            v <<= 4;
            if (c >= '0' && c <= '9') v |= static_cast<unsigned>(c - '0');
            else if (c >= 'a' && c <= 'f') v |= static_cast<unsigned>(c - 'a' + 10);
            else if (c >= 'A' && c <= 'F') v |= static_cast<unsigned>(c - 'A' + 10);
            else fail("bad \\u escape");
        }
        return v;
    }

    static void utf8(std::string& out, unsigned cp) {
        if (cp < 0x80) {
            out += static_cast<char>(cp);
        } else if (cp < 0x800) {
            out += static_cast<char>(0xC0 | (cp >> 6));
            out += static_cast<char>(0x80 | (cp & 0x3F));
        } else if (cp < 0x10000) {
            out += static_cast<char>(0xE0 | (cp >> 12));
            out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
            out += static_cast<char>(0x80 | (cp & 0x3F));
        } else {
            out += static_cast<char>(0xF0 | (cp >> 18));
            out += static_cast<char>(0x80 | ((cp >> 12) & 0x3F));
            out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
            out += static_cast<char>(0x80 | (cp & 0x3F));
        }
    }

    std::string string() {
        ++p_;  // opening quote
        std::string out;
        while (p_ < s_.size() && s_[p_] != '"') {
            char c = s_[p_++];
            if (c != '\\') {
                out += c;  // UTF-8 passes through
                continue;
            }
            if (p_ >= s_.size()) {
                fail("unterminated escape");
            }
            c = s_[p_++];
            switch (c) {
                case '"': out += '"'; break;
                case '\\': out += '\\'; break;
                case '/': out += '/'; break;
                case 'b': out += '\b'; break;
                case 'f': out += '\f'; break;
                case 'n': out += '\n'; break;
                case 'r': out += '\r'; break;
                case 't': out += '\t'; break;
                case 'u': {
                    unsigned cp = hex4();
                    if (cp >= 0xD800 && cp < 0xDC00 && consume("\\u")) {
                        cp = 0x10000 + ((cp - 0xD800) << 10) + (hex4() - 0xDC00);
                    }
                    utf8(out, cp);
                    break;
                }
                default: fail("bad escape");
            }
        }
        if (p_ >= s_.size()) {
            fail("unterminated string");
        }
        ++p_;  // closing quote
        return out;
    }

    Json number() {
        const std::size_t start = p_;
        bool real = false;
        if (p_ < s_.size() && s_[p_] == '-') ++p_;
        while (p_ < s_.size()) {
            const char c = s_[p_];
            if (c >= '0' && c <= '9') {
                ++p_;
            } else if (c == '.' || c == 'e' || c == 'E' || c == '+' || c == '-') {
                real = true;
                ++p_;
            } else {
                break;
            }
        }
        const std::string tok(s_.substr(start, p_ - start));
        if (tok.empty() || tok == "-") {
            fail("unexpected character");
        }
        std::size_t used = 0;
        try {
            if (!real) {
                const long long v = std::stoll(tok, &used);
                if (used == tok.size()) {
                    return Json(v);
                }
            }
            const double d = std::stod(tok, &used);
            if (used == tok.size()) {
                return Json(d);
            }
        } catch (const std::exception&) {
        }
        fail("bad number '" + tok + "'");
    }

    std::string_view s_;
    std::size_t p_{0};
};

}  // namespace

Json Json::parse(std::string_view text) {
    return Parser(text).document();
}

Json read_json_file(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("cannot read " + path.string());
    }
    std::ostringstream s;
    s << in.rdbuf();
    try {
        return Json::parse(s.str());
    } catch (const std::runtime_error& ex) {
        throw std::runtime_error(path.string() + ": " + ex.what());
    }
}

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
