#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace mosflm::io {

/// Minimal JSON value for ProcessingLab's run records (``input.json``,
/// ``process.json``): null, bool, number, string, array and object, with
/// objects keeping insertion order so records read top to bottom.
class Json {
public:
    using Array = std::vector<Json>;
    using Object = std::vector<std::pair<std::string, Json>>;

    Json() = default;
    Json(std::nullptr_t) {}
    Json(bool v) : v_(v) {}
    Json(double v) : v_(v) {}
    Json(int v) : v_(static_cast<std::int64_t>(v)) {}
    Json(long v) : v_(static_cast<std::int64_t>(v)) {}
    Json(long long v) : v_(static_cast<std::int64_t>(v)) {}
    Json(unsigned long v) : v_(static_cast<std::int64_t>(v)) {}
    Json(unsigned long long v) : v_(static_cast<std::int64_t>(v)) {}
    Json(const char* v) : v_(std::string(v)) {}
    Json(std::string v) : v_(std::move(v)) {}
    Json(Array v) : v_(std::move(v)) {}
    Json(Object v) : v_(std::move(v)) {}

    [[nodiscard]] static Json object() { return Json(Object{}); }
    [[nodiscard]] static Json array() { return Json(Array{}); }

    [[nodiscard]] bool is_null() const { return std::holds_alternative<std::monostate>(v_); }
    [[nodiscard]] bool is_object() const { return std::holds_alternative<Object>(v_); }
    [[nodiscard]] bool is_array() const { return std::holds_alternative<Array>(v_); }
    [[nodiscard]] bool is_string() const { return std::holds_alternative<std::string>(v_); }
    [[nodiscard]] bool is_number() const {
        return std::holds_alternative<double>(v_) || std::holds_alternative<std::int64_t>(v_);
    }

    /// Object member, inserted (as null) if absent. Turns a null into an object.
    Json& operator[](const std::string& key);
    /// Object member or nullptr.
    [[nodiscard]] const Json* find(const std::string& key) const;
    /// Append to an array (a null becomes an array).
    void push_back(Json v);

    [[nodiscard]] const Array& items() const;
    /// Mutable array elements (a null becomes an empty array).
    [[nodiscard]] Array& elements();
    [[nodiscard]] const Object& members() const;
    [[nodiscard]] double number() const;  ///< 0 if not a number
    [[nodiscard]] const std::string& string() const;  ///< "" if not a string
    [[nodiscard]] bool boolean() const;  ///< false if not a bool

    /// Parse JSON text (objects keep their member order).
    ///
    /// :raises std::runtime_error: With the offset, if ``text`` is not JSON.
    [[nodiscard]] static Json parse(std::string_view text);

    /// Serialise; non-finite numbers are written as ``null``.
    [[nodiscard]] std::string dump(int indent = 2) const;

private:
    void dump_to(std::string& out, int indent, int depth) const;
    std::variant<std::monostate, bool, std::int64_t, double, std::string, Array, Object> v_;
};

/// Read and parse a JSON file.
///
/// :raises std::runtime_error: If it cannot be read or is not JSON.
[[nodiscard]] Json read_json_file(const std::filesystem::path& path);

/// Write ``text`` to ``path`` atomically: a temporary file in the same
/// directory is written, flushed to disk and renamed over ``path``, so a
/// reader never sees a half-written file.
///
/// :raises std::runtime_error: If the file cannot be written.
void write_file_atomic(const std::filesystem::path& path, const std::string& text);

}  // namespace mosflm::io
