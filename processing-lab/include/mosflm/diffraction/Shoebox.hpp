#pragma once

#include "mosflm/diffraction/SummationIntegrator.hpp"
#include "mosflm/geometry/DetectorGeometry.hpp"
#include "mosflm/image/DiffractionImage.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <list>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mosflm::diffraction {

/// LRU of loaded diffraction frames so PROCESS 3D does not reopen HDF5 per
/// spot. Returned images stay alive via ``shared_ptr`` even if later evicted.
class FrameCache {
public:
    explicit FrameCache(std::size_t capacity = 16);

    /// Load or return a cached frame. ``frame`` is 1-based (ImageIO).
    [[nodiscard]] std::shared_ptr<const image::DiffractionImage> get(
        const std::filesystem::path& path,
        int frame
    );

    /// Insert an already-loaded image (e.g. PROCESS's first-frame peek).
    void put(
        const std::filesystem::path& path,
        int frame,
        image::DiffractionImage image
    );

    [[nodiscard]] std::size_t size() const;
    [[nodiscard]] std::size_t capacity() const { return capacity_; }

private:
    using Key = std::pair<std::string, int>;
    struct KeyHash {
        std::size_t operator()(const Key& k) const noexcept;
    };
    using ImagePtr = std::shared_ptr<image::DiffractionImage>;

    void touch_unlocked(const Key& key);
    void evict_if_needed_unlocked();

    std::size_t capacity_;
    mutable std::mutex mutex_;
    std::list<Key> order_;
    std::unordered_map<Key, std::pair<ImagePtr, std::list<Key>::iterator>, KeyHash>
        items_;
};

/// Extracted (x, y, φ) volume around a predicted reflection.
///
/// Pixels cover the measurement box including the background rim:
/// ``side = 2·outer_half+1`` with ``outer_half = peak_half + bg_rim``.
/// Indexing: ``iphi`` in ``[0, nphi)``, ``dr``/``dc`` relative to the peak
/// centre in ``[-outer_half, outer_half]``. The predicted φ bin is always
/// ``nphi/2`` (odd ``nphi``).
struct Shoebox {
    int h{};
    int k{};
    int l{};
    int col0{};
    int row0{};
    int peak_half{};
    int outer_half{};
    int nphi{};
    int peak_frame{};
    int first_frame{};
    int width{};
    int height{};
    double x_mm{};
    double y_mm{};
    double phi_calc_deg{};
    /// ``nphi * side * side`` voxels.
    std::vector<std::int32_t> pixels;
    /// 1 if the voxel is on-detector and the φ-slice loaded.
    std::vector<std::uint8_t> valid;
    bool ok{false};

    [[nodiscard]] int side() const { return 2 * outer_half + 1; }

    [[nodiscard]] std::size_t index(int iphi, int dr, int dc) const;

    [[nodiscard]] std::int32_t at(int iphi, int dr, int dc) const;

    [[nodiscard]] bool is_valid(int iphi, int dr, int dc) const;
};

/// Inputs for :func:`extract_shoebox`.
struct ShoeboxRequest {
    int h{};
    int k{};
    int l{};
    double x_mm{};
    double y_mm{};
    double phi_calc_deg{};
    int peak_frame{};
    int nphi{5};
    int peak_half_px{5};
    int bg_rim_px{2};
    /// Local COM refine of the peak centre on the φ-centre slice (pixels).
    /// 0 disables.
    double centroid_max_shift_px{8.0};
};

/// Map predicted millimetres to integer peak-centre pixels (same convention
/// as MeasurementBox / StandardProfile).
void shoebox_peak_pixel(
    const DetectorGeometry& detector,
    const IntegrationSettings& settings,
    double x_mm,
    double y_mm,
    int& col,
    int& row
);

/// Force shoebox depth to a positive odd integer (NADD).
[[nodiscard]] int odd_nadd(int nphi);

/// Extract a shoebox centred on ``req.peak_frame`` spanning ``nphi`` frames.
/// Missing edge frames (outside ``[first_allowed, last_allowed]`` or unreadable)
/// are marked invalid. ``last_allowed <= 0`` means no upper clamp.
[[nodiscard]] Shoebox extract_shoebox(
    FrameCache& cache,
    const std::function<std::filesystem::path(int frame)>& path_of,
    const ShoeboxRequest& req,
    const DetectorGeometry& detector,
    const IntegrationSettings& settings,
    int first_allowed_frame = 1,
    int last_allowed_frame = 0
);

}  // namespace mosflm::diffraction
