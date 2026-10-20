#include "mosflm/diffraction/Shoebox.hpp"

#include "mosflm/geometry/LabFrame.hpp"
#include "mosflm/image/ImageIO.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace mosflm::diffraction {
namespace {

void to_pixel(
    const DetectorGeometry& detector,
    const IntegrationSettings& settings,
    double x_mm,
    double y_mm,
    int& col,
    int& row
) {
    if (settings.mosflm_detector_axes) {
        double xdet = 0.0;
        double ydet = 0.0;
        mosflm::cpp_mm_to_mosflm_xdet_yscale(
            x_mm,
            y_mm,
            settings.image_beam_x_mm,
            settings.image_beam_y_mm,
            detector.pixel_size_x_mm,
            detector.pixel_size_y_mm,
            detector.yscale,
            xdet,
            ydet
        );
        col = static_cast<int>(std::lround(xdet));
        row = static_cast<int>(std::lround(ydet));
        return;
    }
    SummationIntegrator::mm_to_image_pixel(detector, x_mm, y_mm, col, row);
}

}  // namespace

std::size_t FrameCache::KeyHash::operator()(const Key& k) const noexcept {
    const std::size_t h1 = std::hash<std::string>{}(k.first);
    const std::size_t h2 = std::hash<int>{}(k.second);
    return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
}

FrameCache::FrameCache(std::size_t capacity)
    : capacity_(capacity < 1 ? 1 : capacity) {}

void FrameCache::touch_unlocked(const Key& key) {
    auto it = items_.find(key);
    if (it == items_.end()) {
        return;
    }
    order_.erase(it->second.second);
    order_.push_front(key);
    it->second.second = order_.begin();
}

void FrameCache::evict_if_needed_unlocked() {
    while (items_.size() >= capacity_ && !order_.empty()) {
        const Key last = order_.back();
        order_.pop_back();
        items_.erase(last);
    }
}

std::shared_ptr<const image::DiffractionImage> FrameCache::get(
    const std::filesystem::path& path,
    int frame
) {
    const Key key{path.string(), frame};
    {
        std::lock_guard<std::mutex> lock(mutex_);
        auto it = items_.find(key);
        if (it != items_.end()) {
            touch_unlocked(key);
            return it->second.first;
        }
    }
    auto loaded = std::make_shared<image::DiffractionImage>(
        image::read_diffraction_image(path, frame)
    );
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = items_.find(key);
    if (it != items_.end()) {
        touch_unlocked(key);
        return it->second.first;
    }
    evict_if_needed_unlocked();
    order_.push_front(key);
    auto [ins, _] = items_.emplace(
        key, std::make_pair(std::move(loaded), order_.begin())
    );
    return ins->second.first;
}

void FrameCache::put(
    const std::filesystem::path& path,
    int frame,
    image::DiffractionImage image
) {
    const Key key{path.string(), frame};
    auto loaded = std::make_shared<image::DiffractionImage>(std::move(image));
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = items_.find(key);
    if (it != items_.end()) {
        it->second.first = std::move(loaded);
        touch_unlocked(key);
        return;
    }
    evict_if_needed_unlocked();
    order_.push_front(key);
    items_.emplace(key, std::make_pair(std::move(loaded), order_.begin()));
}

std::size_t FrameCache::size() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return items_.size();
}

int odd_nadd(int nphi) {
    int n = nphi < 1 ? 1 : nphi;
    if ((n % 2) == 0) {
        ++n;
    }
    return n;
}

void shoebox_peak_pixel(
    const DetectorGeometry& detector,
    const IntegrationSettings& settings,
    double x_mm,
    double y_mm,
    int& col,
    int& row
) {
    to_pixel(detector, settings, x_mm, y_mm, col, row);
}

std::size_t Shoebox::index(int iphi, int dr, int dc) const {
    const int s = side();
    return static_cast<std::size_t>(iphi) * static_cast<std::size_t>(s) *
               static_cast<std::size_t>(s) +
           static_cast<std::size_t>(dr + outer_half) * static_cast<std::size_t>(s) +
           static_cast<std::size_t>(dc + outer_half);
}

std::int32_t Shoebox::at(int iphi, int dr, int dc) const {
    if (iphi < 0 || iphi >= nphi || dr < -outer_half || dr > outer_half ||
        dc < -outer_half || dc > outer_half) {
        return 0;
    }
    const std::size_t idx = index(iphi, dr, dc);
    if (idx >= pixels.size()) {
        return 0;
    }
    return pixels[idx];
}

bool Shoebox::is_valid(int iphi, int dr, int dc) const {
    if (iphi < 0 || iphi >= nphi || dr < -outer_half || dr > outer_half ||
        dc < -outer_half || dc > outer_half) {
        return false;
    }
    const std::size_t idx = index(iphi, dr, dc);
    if (idx >= valid.size()) {
        return false;
    }
    return valid[idx] != 0;
}

Shoebox extract_shoebox(
    FrameCache& cache,
    const std::function<std::filesystem::path(int frame)>& path_of,
    const ShoeboxRequest& req,
    const DetectorGeometry& detector,
    const IntegrationSettings& settings,
    int first_allowed_frame,
    int last_allowed_frame
) {
    Shoebox box;
    box.h = req.h;
    box.k = req.k;
    box.l = req.l;
    box.x_mm = req.x_mm;
    box.y_mm = req.y_mm;
    box.phi_calc_deg = req.phi_calc_deg;
    box.peak_frame = req.peak_frame;
    box.nphi = odd_nadd(req.nphi);
    box.peak_half = std::max(1, req.peak_half_px);
    box.outer_half = box.peak_half + std::max(0, req.bg_rim_px);
    const int half_phi = box.nphi / 2;
    box.first_frame = req.peak_frame - half_phi;

    shoebox_peak_pixel(
        detector, settings, req.x_mm, req.y_mm, box.col0, box.row0
    );

    const auto fill = [&]() -> bool {
        const int side = box.side();
        const std::size_t nvox =
            static_cast<std::size_t>(box.nphi) * static_cast<std::size_t>(side) *
            static_cast<std::size_t>(side);
        box.pixels.assign(nvox, 0);
        box.valid.assign(nvox, 0);
        bool centre_ok = false;
        for (int iphi = 0; iphi < box.nphi; ++iphi) {
            const int frame = box.first_frame + iphi;
            if (frame < first_allowed_frame) {
                continue;
            }
            if (last_allowed_frame > 0 && frame > last_allowed_frame) {
                continue;
            }
            std::shared_ptr<const image::DiffractionImage> img;
            try {
                img = cache.get(path_of(frame), frame);
            } catch (const std::exception&) {
                continue;
            }
            if (!img) {
                continue;
            }
            box.width = img->metadata.width;
            box.height = img->metadata.height;
            const int w = box.width;
            const int h = box.height;
            if (w <= 0 || h <= 0) {
                continue;
            }
            if (iphi == half_phi) {
                centre_ok = true;
            }
            const auto& pix = img->pixels;
            const std::size_t nimg = pix.size();
            for (int dr = -box.outer_half; dr <= box.outer_half; ++dr) {
                for (int dc = -box.outer_half; dc <= box.outer_half; ++dc) {
                    const int col = box.col0 + dc;
                    const int row = box.row0 + dr;
                    if (col < 0 || row < 0 || col >= w || row >= h) {
                        continue;
                    }
                    const std::size_t src =
                        static_cast<std::size_t>(row) *
                            static_cast<std::size_t>(w) +
                        static_cast<std::size_t>(col);
                    if (src >= nimg) {
                        continue;
                    }
                    const std::size_t dst = box.index(iphi, dr, dc);
                    box.pixels[dst] = pix[src];
                    box.valid[dst] = 1;
                }
            }
        }
        return centre_ok && box.width > 0 && box.height > 0;
    };

    box.ok = fill();
    if (box.ok && req.centroid_max_shift_px > 0.0) {
        const int cphi = box.nphi / 2;
        const int search = std::min(
            box.outer_half,
            std::max(
                box.peak_half,
                static_cast<int>(std::ceil(req.centroid_max_shift_px))
            )
        );
        double rim = 0.0;
        int nrim = 0;
        for (int dr = -box.outer_half; dr <= box.outer_half; ++dr) {
            for (int dc = -box.outer_half; dc <= box.outer_half; ++dc) {
                if (std::abs(dc) <= box.peak_half &&
                    std::abs(dr) <= box.peak_half) {
                    continue;
                }
                if (!box.is_valid(cphi, dr, dc)) {
                    continue;
                }
                rim += static_cast<double>(box.at(cphi, dr, dc));
                ++nrim;
            }
        }
        const double bg = (nrim > 0) ? (rim / static_cast<double>(nrim)) : 0.0;
        double peak_v = 0.0;
        int peak_dc = 0;
        int peak_dr = 0;
        double wsum = 0.0;
        double wc = 0.0;
        double wr = 0.0;
        for (int dr = -search; dr <= search; ++dr) {
            for (int dc = -search; dc <= search; ++dc) {
                if (!box.is_valid(cphi, dr, dc)) {
                    continue;
                }
                const double v =
                    std::max(0.0, static_cast<double>(box.at(cphi, dr, dc)) - bg);
                if (v > peak_v) {
                    peak_v = v;
                    peak_dc = dc;
                    peak_dr = dr;
                }
                wsum += v;
                wc += v * static_cast<double>(dc);
                wr += v * static_cast<double>(dr);
            }
        }
        const double cap = req.centroid_max_shift_px;
        int dcol = 0;
        int drow = 0;
        if (peak_v > 3.0 * std::max(1.0, std::sqrt(std::max(bg, 1.0))) &&
            std::hypot(static_cast<double>(peak_dc), static_cast<double>(peak_dr)) <=
                cap) {
            dcol = peak_dc;
            drow = peak_dr;
        } else if (wsum > 10.0) {
            const double fc = wc / wsum;
            const double fr = wr / wsum;
            if (std::hypot(fc, fr) <= cap) {
                dcol = static_cast<int>(std::lround(fc));
                drow = static_cast<int>(std::lround(fr));
            }
        }
        if (dcol != 0 || drow != 0) {
            box.col0 += dcol;
            box.row0 += drow;
            box.ok = fill();
        }
    }
    return box;
}

}  // namespace mosflm::diffraction
