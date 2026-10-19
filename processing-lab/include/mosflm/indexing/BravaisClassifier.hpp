#pragma once

#include <Eigen/Dense>

#include <array>
#include <string>
#include <vector>

namespace mosflm::indexing {

/// One DPS lattice-character recommendation (characters 1–44).
struct BravaisCandidate {
    int character{};          ///< Lattice character index (1–44)
    std::string system_id;    ///< Two-char code, e.g. "cP", "oC"
    float criterion{};        ///< DPS goof / FOM (lower is better)
    Eigen::Matrix3d orientation{Eigen::Matrix3d::Identity()};  ///< Reciprocal OM
};

/// Thin C++ façade over the Cambridge 5-argument DPS ``bravais()`` API.
///
/// Links the existing C implementation (`index/bravais.c` via `dps_index.a`);
/// does not reimplement lattice-character logic. CCTBX may validate transformed
/// cell metrics separately (see ``docs/CCTBX_EVAL.md``).
class BravaisClassifier {
public:
    /// Classify from a reciprocal orientation matrix (DPS convention:
    /// rows are reciprocal basis vectors in lab frame, Å⁻¹).
    ///
    /// \param reciprocal_om  3×3 matrix; row ``i`` is reciprocal axis ``i``.
    /// \param output_type    Passed through to ``bravais`` (1 = fill arrays).
    [[nodiscard]] std::vector<BravaisCandidate> classify(
        const Eigen::Matrix3d& reciprocal_om,
        int output_type = 1
    ) const;

    /// Convenience: build a diagonal reciprocal OM for an orthogonal cell
    /// (``diag(1/a, 1/b, 1/c)``) and classify.
    [[nodiscard]] std::vector<BravaisCandidate> classify_orthogonal_cell(
        double a_A,
        double b_A,
        double c_A
    ) const;
};

}  // namespace mosflm::indexing
