#include "mosflm/indexing/LocalIndexAssign.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

namespace mosflm::indexing {
namespace {

struct Edge {
    int i{0};
    int j{0};
    double len{0.0};
    int dh{0};
    int dk{0};
    int dl{0};
};

double branch_length(
    const Eigen::Matrix3d& Ainv,
    const Eigen::Vector3d& ri,
    const Eigen::Vector3d& rj,
    double epsilon,
    int delta_max,
    int& dh,
    int& dk,
    int& dl
) {
    const Eigen::Vector3d d = Ainv * (rj - ri);
    dh = static_cast<int>(std::lround(d[0]));
    dk = static_cast<int>(std::lround(d[1]));
    dl = static_cast<int>(std::lround(d[2]));
    const double e0 = std::fabs(d[0] - static_cast<double>(dh));
    const double e1 = std::fabs(d[1] - static_cast<double>(dk));
    const double e2 = std::fabs(d[2] - static_cast<double>(dl));
    const int a0 = std::abs(dh);
    const int a1 = std::abs(dk);
    const int a2 = std::abs(dl);
    if (a0 > delta_max || a1 > delta_max || a2 > delta_max) {
        return 1.0 + e0 + e1 + e2;
    }
    if (e0 <= epsilon && e1 <= epsilon && e2 <= epsilon) {
        return e0 + e1 + e2;
    }
    return 1.0 + e0 + e1 + e2;
}

}  // namespace

LocalIndexResult assign_local_indices(
    const Eigen::Matrix3d& amat,
    const std::vector<float>& x,
    const std::vector<float>& y,
    const std::vector<float>& z,
    const LocalIndexSettings& settings
) {
    LocalIndexResult result;
    const std::size_t n_in = x.size();
    if (n_in == 0 || y.size() != n_in || z.size() != n_in) {
        return result;
    }
    Eigen::FullPivLU<Eigen::Matrix3d> lu(amat);
    if (!lu.isInvertible()) {
        return result;
    }
    const Eigen::Matrix3d Ainv = lu.inverse();

    const std::size_t n = std::min(n_in, settings.max_spots);
    result.spots.assign(n_in, LocalIndexSpot{});

    std::vector<Eigen::Vector3d> r(n);
    for (std::size_t i = 0; i < n; ++i) {
        r[i] = Eigen::Vector3d(
            static_cast<double>(x[i]),
            static_cast<double>(y[i]),
            static_cast<double>(z[i])
        );
    }

    if (n < 2) {
        if (n == 1) {
            const Eigen::Vector3d hkl = Ainv * r[0];
            result.spots[0].h = static_cast<int>(std::lround(hkl[0]));
            result.spots[0].k = static_cast<int>(std::lround(hkl[1]));
            result.spots[0].l = static_cast<int>(std::lround(hkl[2]));
            result.spots[0].indexed = true;
            result.spots[0].subtree = 1;
            result.n_indexed = 1;
            result.n_largest_subtree = 1;
        }
        return result;
    }

    std::vector<Edge> edges;
    edges.reserve(n * (n - 1) / 2);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            Edge e;
            e.i = static_cast<int>(i);
            e.j = static_cast<int>(j);
            e.len = branch_length(
                Ainv, r[i], r[j], settings.epsilon, settings.delta_max, e.dh,
                e.dk, e.dl
            );
            edges.push_back(e);
        }
    }
    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        return a.len < b.len;
    });

    // Kruskal MST.
    std::vector<int> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int a) {
        while (parent[static_cast<std::size_t>(a)] != a) {
            parent[static_cast<std::size_t>(a)] =
                parent[static_cast<std::size_t>(
                    parent[static_cast<std::size_t>(a)]
                )];
            a = parent[static_cast<std::size_t>(a)];
        }
        return a;
    };
    auto unite = [&](int a, int b) {
        a = find(a);
        b = find(b);
        if (a != b) {
            parent[static_cast<std::size_t>(a)] = b;
        }
    };

    std::vector<Edge> mst;
    mst.reserve(n - 1);
    for (const Edge& e : edges) {
        if (find(e.i) == find(e.j)) {
            continue;
        }
        unite(e.i, e.j);
        mst.push_back(e);
        if (mst.size() + 1 >= n) {
            break;
        }
    }

    // Adjacency for tree walk.
    std::vector<std::vector<std::pair<int, Edge>>> adj(n);
    for (const Edge& e : mst) {
        adj[static_cast<std::size_t>(e.i)].push_back({e.j, e});
        Edge rev = e;
        rev.dh = -e.dh;
        rev.dk = -e.dk;
        rev.dl = -e.dl;
        std::swap(rev.i, rev.j);
        adj[static_cast<std::size_t>(e.j)].push_back({e.i, rev});
    }

    // Root = node with most MST neighbours.
    int root = 0;
    std::size_t best_deg = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (adj[i].size() > best_deg) {
            best_deg = adj[i].size();
            root = static_cast<int>(i);
        }
    }

    const Eigen::Vector3d h0 = Ainv * r[static_cast<std::size_t>(root)];
    std::vector<int> hh(n, 0);
    std::vector<int> kk(n, 0);
    std::vector<int> ll(n, 0);
    std::vector<int> subtree(n, 0);
    std::vector<char> seen(n, 0);
    hh[static_cast<std::size_t>(root)] = static_cast<int>(std::lround(h0[0]));
    kk[static_cast<std::size_t>(root)] = static_cast<int>(std::lround(h0[1]));
    ll[static_cast<std::size_t>(root)] = static_cast<int>(std::lround(h0[2]));
    subtree[static_cast<std::size_t>(root)] = 1;
    seen[static_cast<std::size_t>(root)] = 1;

    std::vector<int> stack;
    stack.push_back(root);
    int next_sub = 1;
    while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        for (const auto& [v, e] : adj[static_cast<std::size_t>(u)]) {
            if (seen[static_cast<std::size_t>(v)]) {
                continue;
            }
            seen[static_cast<std::size_t>(v)] = 1;
            hh[static_cast<std::size_t>(v)] =
                hh[static_cast<std::size_t>(u)] + e.dh;
            kk[static_cast<std::size_t>(v)] =
                kk[static_cast<std::size_t>(u)] + e.dk;
            ll[static_cast<std::size_t>(v)] =
                ll[static_cast<std::size_t>(u)] + e.dl;
            if (e.len >= settings.ell_min) {
                ++next_sub;
                subtree[static_cast<std::size_t>(v)] = next_sub;
            } else {
                subtree[static_cast<std::size_t>(v)] =
                    subtree[static_cast<std::size_t>(u)];
            }
            stack.push_back(v);
        }
    }

    // Count subtree sizes; keep only the largest (+ ties with size ≥ 50% of max).
    std::vector<std::size_t> sub_count(
        static_cast<std::size_t>(next_sub + 1), 0
    );
    for (std::size_t i = 0; i < n; ++i) {
        if (seen[i] && subtree[i] > 0) {
            ++sub_count[static_cast<std::size_t>(subtree[i])];
        }
    }
    std::size_t max_sub = 0;
    for (std::size_t c : sub_count) {
        max_sub = std::max(max_sub, c);
    }
    const std::size_t keep_cut =
        std::max<std::size_t>(3, max_sub / 2);

    for (std::size_t i = 0; i < n; ++i) {
        if (!seen[i]) {
            // Global round fallback for isolates.
            const Eigen::Vector3d hkl = Ainv * r[i];
            result.spots[i].h = static_cast<int>(std::lround(hkl[0]));
            result.spots[i].k = static_cast<int>(std::lround(hkl[1]));
            result.spots[i].l = static_cast<int>(std::lround(hkl[2]));
            result.spots[i].indexed = true;
            result.spots[i].subtree = 0;
            ++result.n_indexed;
            continue;
        }
        const std::size_t sc =
            sub_count[static_cast<std::size_t>(subtree[i])];
        if (sc < keep_cut && sc < max_sub) {
            continue;  // alien subtree
        }
        result.spots[i].h = hh[i];
        result.spots[i].k = kk[i];
        result.spots[i].l = ll[i];
        result.spots[i].indexed = true;
        result.spots[i].subtree = subtree[i];
        ++result.n_indexed;
    }
    result.n_largest_subtree = max_sub;

    // Remaining spots beyond MST cap: global round.
    for (std::size_t i = n; i < n_in; ++i) {
        const Eigen::Vector3d ri(
            static_cast<double>(x[i]),
            static_cast<double>(y[i]),
            static_cast<double>(z[i])
        );
        const Eigen::Vector3d hkl = Ainv * ri;
        result.spots[i].h = static_cast<int>(std::lround(hkl[0]));
        result.spots[i].k = static_cast<int>(std::lround(hkl[1]));
        result.spots[i].l = static_cast<int>(std::lround(hkl[2]));
        result.spots[i].indexed = true;
        ++result.n_indexed;
    }
    return result;
}

}  // namespace mosflm::indexing
