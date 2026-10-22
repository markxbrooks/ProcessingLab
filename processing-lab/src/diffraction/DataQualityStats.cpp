#include "mosflm/diffraction/DataQualityStats.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <numeric>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

namespace mosflm::diffraction {
namespace {

/// FLAG bit for a negative profile-fitted intensity (ProcessPipeline). A
/// negative I is a valid measurement of a weak reflection, not a misfit:
/// rejecting it empties the high-resolution shells and biases them upward.
constexpr int kFlagNegProfile = 4;

using Miller = std::tuple<int, int, int>;

struct Obs {
    double intensity{0.0};
    double sigma{0.0};
    double d_A{0.0};
};

int shell_index(double d_A, double dmax_A, double dmin_A, int n_shells) {
    if (!(dmax_A > dmin_A) || n_shells < 1) {
        return 0;
    }
    // Equal width in 1/d (XDS-like); shell 0 = lowest resolution.
    const double s_low = 1.0 / dmax_A;
    const double s_high = 1.0 / dmin_A;
    const double s = 1.0 / d_A;
    if (s <= s_low) {
        return 0;
    }
    if (s >= s_high) {
        return n_shells - 1;
    }
    const double t = (s - s_low) / (s_high - s_low);
    int idx = static_cast<int>(std::floor(t * static_cast<double>(n_shells)));
    if (idx < 0) {
        idx = 0;
    }
    if (idx >= n_shells) {
        idx = n_shells - 1;
    }
    return idx;
}

double shell_d_limit(int shell, double dmax_A, double dmin_A, int n_shells) {
    const double s_low = 1.0 / dmax_A;
    const double s_high = 1.0 / dmin_A;
    const double s =
        s_low + (s_high - s_low) *
                    (static_cast<double>(shell + 1) /
                     static_cast<double>(n_shells));
    return 1.0 / s;
}

double pearson_cc(
    const std::vector<double>& a,
    const std::vector<double>& b
) {
    const std::size_t n = a.size();
    if (n < 3 || b.size() != n) {
        return 0.0;
    }
    double mean_a = 0.0;
    double mean_b = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        mean_a += a[i];
        mean_b += b[i];
    }
    mean_a /= static_cast<double>(n);
    mean_b /= static_cast<double>(n);
    double num = 0.0;
    double da = 0.0;
    double db = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double xa = a[i] - mean_a;
        const double xb = b[i] - mean_b;
        num += xa * xb;
        da += xa * xa;
        db += xb * xb;
    }
    if (!(da > 0.0) || !(db > 0.0)) {
        return 0.0;
    }
    return num / std::sqrt(da * db);
}

/// Rough 0.1% two-sided significance for Pearson r (Karplus-style mark).
bool cc_significant(double cc, int n_pairs) {
    if (n_pairs < 3) {
        return false;
    }
    // z ≈ 3.29 for two-sided p=0.001; r_crit ≈ tanh(z/sqrt(n-3)).
    const double df = static_cast<double>(std::max(n_pairs - 3, 1));
    const double r_crit = std::tanh(3.291 / std::sqrt(df));
    return std::abs(cc) >= r_crit;
}

void accumulate_group_stats(
    const std::vector<Obs>& group,
    double& sum_abs_dev,
    double& sum_rmeas_num,
    double& sum_i_abs,
    double& sum_sigma,
    int& compared,
    double& sum_merged_isigma,
    int& n_merged_isigma,
    std::vector<double>& half_a,
    std::vector<double>& half_b,
    std::mt19937& rng
) {
    const int n = static_cast<int>(group.size());
    if (n < 1) {
        return;
    }
    double sum_i = 0.0;
    double sum_w = 0.0;
    double sum_wi = 0.0;
    for (const auto& o : group) {
        sum_i += o.intensity;
        const double w =
            (o.sigma > 0.0) ? (1.0 / (o.sigma * o.sigma)) : 0.0;
        sum_w += w;
        sum_wi += w * o.intensity;
    }
    const double mean =
        (sum_w > 0.0) ? (sum_wi / sum_w) : (sum_i / static_cast<double>(n));
    const double sigma_mean =
        (sum_w > 0.0) ? (1.0 / std::sqrt(sum_w)) : group[0].sigma;

    if (sigma_mean > 0.0 && std::isfinite(mean) && std::isfinite(sigma_mean)) {
        sum_merged_isigma += mean / sigma_mean;
        ++n_merged_isigma;
    }

    if (n < 2) {
        return;
    }

    double abs_dev = 0.0;
    double abs_i = 0.0;
    double sig = 0.0;
    for (const auto& o : group) {
        abs_dev += std::abs(o.intensity - mean);
        abs_i += std::abs(o.intensity);
        sig += o.sigma;
    }
    sum_abs_dev += abs_dev;
    sum_i_abs += abs_i;
    sum_sigma += sig;
    compared += n;
    // R-meas numerator: sqrt(n/(n-1)) * Σ|I_i − Ī|
    sum_rmeas_num +=
        std::sqrt(static_cast<double>(n) / static_cast<double>(n - 1)) *
        abs_dev;

    // Random half-datasets for CC(1/2).
    std::vector<int> order(static_cast<std::size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);
    const int n1 = n / 2;
    if (n1 < 1 || n - n1 < 1) {
        return;
    }
    double s1 = 0.0;
    double s2 = 0.0;
    for (int i = 0; i < n; ++i) {
        if (i < n1) {
            s1 += group[static_cast<std::size_t>(order[static_cast<std::size_t>(i)])]
                      .intensity;
        } else {
            s2 += group[static_cast<std::size_t>(order[static_cast<std::size_t>(i)])]
                      .intensity;
        }
    }
    half_a.push_back(s1 / static_cast<double>(n1));
    half_b.push_back(s2 / static_cast<double>(n - n1));
}

void finalize_shell(
    ResolutionShellStats& s,
    double sum_abs_dev,
    double sum_rmeas_num,
    double sum_i_abs,
    double sum_sigma,
    int compared,
    double sum_merged_isigma,
    int n_merged_isigma,
    const std::vector<double>& half_a,
    const std::vector<double>& half_b
) {
    if (s.possible > 0) {
        s.completeness_pct =
            100.0 * static_cast<double>(s.unique) /
            static_cast<double>(s.possible);
    }
    s.compared = compared;
    if (n_merged_isigma > 0) {
        s.mean_isigma =
            sum_merged_isigma / static_cast<double>(n_merged_isigma);
    } else {
        s.mean_isigma = -99.0;
    }
    if (compared > 0 && sum_i_abs > 0.0) {
        s.r_obs_pct = 100.0 * sum_abs_dev / sum_i_abs;
        // Expected R from Σσ / Σ|I| (sample-statistics estimate).
        s.r_exp_pct = 100.0 * sum_sigma / sum_i_abs;
        s.rmeas_pct = 100.0 * sum_rmeas_num / sum_i_abs;
    } else {
        s.r_obs_pct = -99.9;
        s.r_exp_pct = -99.9;
        s.rmeas_pct = -99.9;
    }
    if (half_a.size() >= 3) {
        const double cc = pearson_cc(half_a, half_b);
        s.cc_half_pct = 100.0 * cc;
        s.cc_significant =
            cc_significant(cc, static_cast<int>(half_a.size()));
    } else {
        s.cc_half_pct = 0.0;
        s.cc_significant = false;
    }
    // Anomalous not computed in PROCESS (FRIEDEL'S_LAW assumed).
    s.anomal_corr_pct = 0.0;
    s.sigano = 0.0;
    s.nano = 0;
}

void count_possible(
    const gemmi::UnitCell& cell,
    const gemmi::SpaceGroup* sg,
    double dmax_A,
    double dmin_A,
    int n_shells,
    std::vector<int>& possible_per_shell,
    int& possible_total
) {
    possible_per_shell.assign(static_cast<std::size_t>(n_shells), 0);
    possible_total = 0;
    if (sg == nullptr || !(dmax_A > dmin_A) || !(dmin_A > 0.0)) {
        return;
    }
    const gemmi::ReciprocalAsu asu(sg);
    const gemmi::GroupOps gops = sg->operations();
    const int hmax =
        static_cast<int>(std::ceil(cell.a / dmin_A)) + 2;
    const int kmax =
        static_cast<int>(std::ceil(cell.b / dmin_A)) + 2;
    const int lmax =
        static_cast<int>(std::ceil(cell.c / dmin_A)) + 2;
    for (int h = -hmax; h <= hmax; ++h) {
        for (int k = -kmax; k <= kmax; ++k) {
            for (int l = -lmax; l <= lmax; ++l) {
                if (h == 0 && k == 0 && l == 0) {
                    continue;
                }
                const gemmi::Op::Miller m{h, k, l};
                if (!asu.is_in(m)) {
                    continue;
                }
                if (gops.is_systematically_absent(m)) {
                    continue;
                }
                const double d = cell.calculate_d(m);
                if (!(d >= dmin_A) || !(d <= dmax_A)) {
                    continue;
                }
                const int sh = shell_index(d, dmax_A, dmin_A, n_shells);
                ++possible_per_shell[static_cast<std::size_t>(sh)];
                ++possible_total;
            }
        }
    }
}

/// ``# `` is two columns, matching the indent on each data row, so the
/// labels share the field widths of :func:`print_shell_row`.
void print_shell_header(std::ostream& os) {
    os << "# " << std::setw(9) << "limit_A"
       << ' ' << std::setw(8) << "observed"
       << ' ' << std::setw(7) << "unique"
       << ' ' << std::setw(8) << "possible"
       << ' ' << std::setw(7) << "compl%"
       << ' ' << std::setw(7) << "R_obs%"
       << ' ' << std::setw(7) << "R_exp%"
       << ' ' << std::setw(8) << "compared"
       << ' ' << std::setw(7) << "I/sigma"
       << ' ' << std::setw(7) << "Rmeas%"
       << ' ' << std::setw(7) << "CC(1/2)"
       << ' ' << std::setw(3) << "sig"
       << ' ' << std::setw(6) << "Anomal"
       << ' ' << std::setw(7) << "SigAno"
       << ' ' << std::setw(6) << "Nano"
       << '\n';
}

void print_shell_row(std::ostream& os, const ResolutionShellStats& s, bool total) {
    const double r_obs = s.compared > 0 ? s.r_obs_pct : -99.9;
    const double r_exp = s.compared > 0 ? s.r_exp_pct : -99.9;
    const double rmeas = s.compared > 0 ? s.rmeas_pct : -99.9;
    os << std::fixed << "  ";
    if (total) {
        os << std::setw(9) << "total";
    } else {
        os << std::setw(9) << std::setprecision(2) << s.d_limit_A;
    }
    os << ' ' << std::setw(8) << s.observed
       << ' ' << std::setw(7) << s.unique
       << ' ' << std::setw(8) << s.possible
       << ' ' << std::setw(7) << std::setprecision(1) << s.completeness_pct
       << ' ' << std::setw(7) << std::setprecision(1) << r_obs
       << ' ' << std::setw(7) << std::setprecision(1) << r_exp
       << ' ' << std::setw(8) << s.compared
       << ' ' << std::setw(7) << std::setprecision(2) << s.mean_isigma
       << ' ' << std::setw(7) << std::setprecision(1) << rmeas
       << ' ' << std::setw(7) << std::setprecision(1) << s.cc_half_pct
       << ' ' << std::setw(3) << (s.cc_significant ? 1 : 0)
       << ' ' << std::setw(6) << std::setprecision(0) << s.anomal_corr_pct
       << ' ' << std::setw(7) << std::setprecision(3) << s.sigano
       << ' ' << std::setw(6) << s.nano
       << '\n';
}

}  // namespace

DataQualityReport compute_data_quality(
    const crystal::CrystalMetadata& crystal,
    const std::vector<crystal::MosflmReflection>& reflections,
    double dmax_A,
    double dmin_A,
    int n_shells,
    double min_isigma
) {
    DataQualityReport report;
    if (n_shells < 1) {
        n_shells = 1;
    }
    if (!(dmax_A > dmin_A)) {
        dmax_A = 200.0;
        dmin_A = 1.45;
    }

    report.shells.resize(static_cast<std::size_t>(n_shells));
    for (int i = 0; i < n_shells; ++i) {
        report.shells[static_cast<std::size_t>(i)].d_limit_A =
            shell_d_limit(i, dmax_A, dmin_A, n_shells);
    }

    std::vector<int> possible_shells;
    int possible_total = 0;
    count_possible(
        crystal.cell,
        crystal.spacegroup,
        dmax_A,
        dmin_A,
        n_shells,
        possible_shells,
        possible_total
    );
    for (int i = 0; i < n_shells; ++i) {
        report.shells[static_cast<std::size_t>(i)].possible =
            possible_shells[static_cast<std::size_t>(i)];
    }
    report.total.possible = possible_total;
    report.total.d_limit_A = dmin_A;

    // Collect accepted observations by ASU HKL, per shell.
    std::vector<std::map<Miller, std::vector<Obs>>> by_shell(
        static_cast<std::size_t>(n_shells)
    );
    std::map<Miller, std::vector<Obs>> by_all;
    double sum_sysabs_isigi = 0.0;
    int n_sysabs_isigi = 0;

    for (const auto& row : reflections) {
        if ((row.flag & ~kFlagNegProfile) != 0) {
            ++report.n_rejected_misfits;
            continue;
        }
        if (!(row.sigma > 0.0f) || !std::isfinite(row.intensity) ||
            !std::isfinite(row.sigma)) {
            continue;
        }
        int h = row.h;
        int k = row.k;
        int l = row.l;
        int misym = 0;
        crystal::MtzAdapter::reduce_to_asu(crystal.spacegroup, h, k, l, misym);
        const gemmi::Op::Miller m{h, k, l};
        if (crystal.spacegroup != nullptr) {
            const gemmi::GroupOps gops = crystal.spacegroup->operations();
            if (gops.is_systematically_absent(m)) {
                ++report.n_systematic_absent;
                sum_sysabs_isigi +=
                    static_cast<double>(row.intensity / row.sigma);
                ++n_sysabs_isigi;
                continue;
            }
        }
        const double isigi =
            static_cast<double>(row.intensity / row.sigma);
        if (isigi < min_isigma) {
            continue;
        }
        const double d = crystal.cell.calculate_d(m);
        if (!(d >= dmin_A) || !(d <= dmax_A)) {
            continue;
        }
        Obs o;
        o.intensity = static_cast<double>(row.intensity);
        o.sigma = static_cast<double>(row.sigma);
        o.d_A = d;
        const Miller key{h, k, l};
        const int sh = shell_index(d, dmax_A, dmin_A, n_shells);
        by_shell[static_cast<std::size_t>(sh)][key].push_back(o);
        by_all[key].push_back(o);
        ++report.n_in_subset;
        ++report.n_accepted_obs;
        ++report.shells[static_cast<std::size_t>(sh)].observed;
        ++report.total.observed;
    }
    if (n_sysabs_isigi > 0) {
        report.mean_isigma_systematic_absent =
            sum_sysabs_isigi / static_cast<double>(n_sysabs_isigi);
    }
    report.n_unique_accepted = static_cast<int>(by_all.size());
    report.total.unique = report.n_unique_accepted;

    std::mt19937 rng(20261018u);

    auto process_map =
        [&](std::map<Miller, std::vector<Obs>>& groups,
            ResolutionShellStats& stats) {
            stats.unique = static_cast<int>(groups.size());
            double sum_abs_dev = 0.0;
            double sum_rmeas_num = 0.0;
            double sum_i_abs = 0.0;
            double sum_sigma = 0.0;
            int compared = 0;
            double sum_merged_isigma = 0.0;
            int n_merged_isigma = 0;
            std::vector<double> half_a;
            std::vector<double> half_b;
            for (auto& kv : groups) {
                accumulate_group_stats(
                    kv.second,
                    sum_abs_dev,
                    sum_rmeas_num,
                    sum_i_abs,
                    sum_sigma,
                    compared,
                    sum_merged_isigma,
                    n_merged_isigma,
                    half_a,
                    half_b,
                    rng
                );
            }
            finalize_shell(
                stats,
                sum_abs_dev,
                sum_rmeas_num,
                sum_i_abs,
                sum_sigma,
                compared,
                sum_merged_isigma,
                n_merged_isigma,
                half_a,
                half_b
            );
        };

    for (int i = 0; i < n_shells; ++i) {
        process_map(
            by_shell[static_cast<std::size_t>(i)],
            report.shells[static_cast<std::size_t>(i)]
        );
    }
    process_map(by_all, report.total);
    report.total.possible = possible_total;
    if (report.total.possible > 0) {
        report.total.completeness_pct =
            100.0 * static_cast<double>(report.total.unique) /
            static_cast<double>(report.total.possible);
    }

    return report;
}

void write_data_quality_table(
    std::ostream& os,
    const DataQualityReport& report
) {
    os << "\n          COMPLETENESS AND QUALITY OF DATA SET\n"
       << "          ------------------------------------\n\n"
       << " NOTE: computed from UNSCALED integrated intensities (no batch\n"
       << "       scale, absorption or decay correction, no outlier\n"
       << "       rejection). R-factors are therefore higher than after\n"
       << "       scaling (AIMLESS / XSCALE / CORRECT); use them to compare\n"
       << "       integration runs, not as final merging statistics.\n\n"
       << " R-FACTOR\n"
       << " observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))\n"
       << " expected = expected R-FACTOR derived from Sigma(I)\n\n"
       << " COMPARED = number of reflections used for calculating R-FACTOR\n"
       << " I/SIGMA  = mean of intensity/Sigma(I) of unique reflections\n"
       << "            (after merging symmetry-related observations)\n"
       << " Sigma(I) = standard deviation of reflection intensity I\n"
       << "            estimated from sample statistics\n\n"
       << " R-meas   = redundancy independent R-factor (intensities)\n"
       << "            Diederichs & Karplus (1997), Nature Struct. Biol. 4,"
          " 269-275.\n\n"
       << " CC(1/2)  = percentage of correlation between intensities from\n"
       << "            random half-datasets. Correlation significant at\n"
       << "            the 0.1% level has sig=1 in the table.\n"
       << "            Karplus & Diederichs (2012), Science 336, 1030-33\n"
       << " Anomal   = percentage of correlation between random half-sets\n"
       << "  Corr      of anomalous intensity differences. Correlation\n"
       << "            significant at the 0.1% level is marked.\n"
       << " SigAno   = mean anomalous difference in units of its estimated\n"
       << "            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)\n"
       << "            are structure factor estimates obtained from the\n"
       << "            merged intensity observations in each parity class.\n"
       << "  Nano    = Number of unique reflections used to calculate\n"
       << "            Anomal_Corr & SigAno. At least two observations\n"
       << "            for each (+ and -) parity are required.\n\n\n"
       << " SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION"
          " OF RESOLUTION (UNSCALED)\n"
       << "RESOLUTION_BINS nshells=" << report.shells.size()
       << " columns=d_limit_A,observed,unique,possible,completeness_pct,"
          "r_obs_pct,r_exp_pct,compared,mean_isigma,rmeas_pct,cc_half_pct,"
          "cc_sig,anomal_corr_pct,sigano,nano\n"
       << "# d_limit_A is the high-resolution edge (Å). The last row is"
          " the whole range (total).\n"
       << "# R_obs, R_exp, and Rmeas are -99.9 when compared is 0."
          " Percent columns are already in percent.\n";
    print_shell_header(os);
    for (const auto& shell : report.shells) {
        print_shell_row(os, shell, /*total=*/false);
    }
    print_shell_row(os, report.total, /*total=*/true);
    os << "RESOLUTION_BINS_END\n";

    os << "\n\n"
       << " NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES "
       << std::setw(7) << report.n_in_subset << "\n"
       << " NUMBER OF REJECTED MISFITS                         "
       << std::setw(7) << report.n_rejected_misfits << "\n"
       << " NUMBER OF SYSTEMATIC ABSENT REFLECTIONS            "
       << std::setw(7) << report.n_systematic_absent << "\n";
    if (report.n_systematic_absent > 0 &&
        report.mean_isigma_systematic_absent > -90.0) {
        os << " MEAN I/SIGMA OF SYSTEMATIC ABSENCES                "
           << std::setw(7) << std::setprecision(2) << std::fixed
           << report.mean_isigma_systematic_absent << "\n";
        // CORRECT-lite: elevated sys-abs I/σ suggests wrong SG / screw.
        if (report.mean_isigma_systematic_absent > 3.0 &&
            report.n_systematic_absent >= 10) {
            os << " *** WARNING **   SYSTEMATIC ABSENCES LOOK STRONG"
                  " — check SPACE GROUP / indexing\n";
        }
    }
    os << " NUMBER OF ACCEPTED OBSERVATIONS                    "
       << std::setw(7) << report.n_accepted_obs << "\n"
       << " NUMBER OF UNIQUE ACCEPTED REFLECTIONS              "
       << std::setw(7) << report.n_unique_accepted << "\n\n"
       << std::flush;
}

}  // namespace mosflm::diffraction
