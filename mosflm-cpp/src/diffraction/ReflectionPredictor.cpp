#include "mosflm/diffraction/ReflectionPredictor.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace mosflm::diffraction {

std::vector<PredictedReflection> ReflectionPredictor::predict(
    const PredictRequest& request
) const {
    (void)request;
    // Vertical-slice placeholder: real PREDICT port fills this list.
    return {};
}

void ReflectionPredictor::write_dump(
    const std::filesystem::path& path,
    const std::vector<PredictedReflection>& spots
) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error(
            "ReflectionPredictor::write_dump: cannot open " + path.string()
        );
    }
    out << "PREDICT_DUMP n=" << spots.size() << '\n';
    for (const auto& s : spots) {
        out << s.h << ' ' << s.k << ' ' << s.l << ' ' << s.x_mm << ' ' << s.y_mm
            << ' ' << s.phi_begin_deg << ' ' << s.phi_end_deg << ' '
            << (s.full ? 1 : 0) << '\n';
    }
}

std::vector<PredictedReflection> ReflectionPredictor::read_dump(
    const std::filesystem::path& path
) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error(
            "ReflectionPredictor::read_dump: cannot open " + path.string()
        );
    }
    std::string header;
    if (!std::getline(in, header)) {
        throw std::runtime_error("ReflectionPredictor::read_dump: empty file");
    }
    if (header.rfind("PREDICT_DUMP", 0) != 0) {
        throw std::runtime_error(
            "ReflectionPredictor::read_dump: missing PREDICT_DUMP header"
        );
    }

    std::vector<PredictedReflection> spots;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        std::istringstream iss(line);
        PredictedReflection s;
        int full_flag = 1;
        if (!(iss >> s.h >> s.k >> s.l >> s.x_mm >> s.y_mm >> s.phi_begin_deg >>
              s.phi_end_deg >> full_flag)) {
            throw std::runtime_error(
                "ReflectionPredictor::read_dump: bad line: " + line
            );
        }
        s.full = (full_flag != 0);
        spots.push_back(s);
    }
    return spots;
}

}  // namespace mosflm::diffraction
