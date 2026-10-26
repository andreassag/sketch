#pragma once

#include "sketch.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <limits>

namespace sketch {

struct DistanceResult {
    std::string query_name;
    std::string ref_name;
    size_t shared_hashes = 0;
    size_t query_hashes = 0;
    size_t ref_hashes = 0;
    double jaccard = 0.0;
    double containment = 0.0;
    double mash_dist = 1.0;
    double ani = 0.0;
};

inline size_t count_shared(const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) noexcept {
    size_t i = 0, j = 0;
    size_t shared = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) {
            shared++;
            i++;
            j++;
        } else if (a[i] < b[j]) {
            i++;
        } else {
            j++;
        }
    }
    return shared;
}

inline DistanceResult compute_distance(const GenomeSketch& query, const GenomeSketch& ref) {
    DistanceResult res;
    res.query_name = query.name;
    res.ref_name = ref.name;
    res.query_hashes = query.hashes.size();
    res.ref_hashes = ref.hashes.size();

    if (query.hashes.empty() || ref.hashes.empty()) {
        return res;
    }

    res.shared_hashes = count_shared(query.hashes, ref.hashes);
    size_t union_size = res.query_hashes + res.ref_hashes - res.shared_hashes;

    if (union_size > 0) {
        res.jaccard = static_cast<double>(res.shared_hashes) / static_cast<double>(union_size);
    }
    if (res.query_hashes > 0) {
        res.containment = static_cast<double>(res.shared_hashes) / static_cast<double>(res.query_hashes);
    }

    uint32_t k = std::min(query.params.k, ref.params.k);
    if (res.jaccard <= 0.0) {
        res.mash_dist = 1.0;
        res.ani = 0.0;
    } else if (res.jaccard >= 1.0) {
        res.mash_dist = 0.0;
        res.ani = 1.0;
    } else {
        // Mash distance: D = -1/k * ln(2J / (1 + J))
        double term = (2.0 * res.jaccard) / (1.0 + res.jaccard);
        res.mash_dist = -std::log(term) / static_cast<double>(k);
        if (res.mash_dist > 1.0) res.mash_dist = 1.0;
        res.ani = std::max(0.0, 1.0 - res.mash_dist);
    }

    return res;
}

// Neighbor-Joining Tree Builder
inline std::string build_neighbor_joining_newick(const std::vector<std::string>& taxa,
                                                const std::vector<std::vector<double>>& raw_dist) {
    size_t n = taxa.size();
    if (n == 0) return "();";
    if (n == 1) return "(" + taxa[0] + ":0.0);";
    if (n == 2) {
        double d = raw_dist[0][1] / 2.0;
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(4);
        ss << "(" << taxa[0] << ":" << d << "," << taxa[1] << ":" << d << ");";
        return ss.str();
    }

    // Active nodes
    std::vector<std::string> nodes = taxa;
    std::vector<std::vector<double>> D = raw_dist;

    while (nodes.size() > 2) {
        size_t m = nodes.size();
        std::vector<double> r(m, 0.0);
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < m; ++j) {
                r[i] += D[i][j];
            }
        }

        // Q-matrix: Q(i,j) = (m - 2)*D(i,j) - r_i - r_j
        double min_q = std::numeric_limits<double>::infinity();
        size_t best_i = 0, best_j = 1;

        for (size_t i = 0; i < m; ++i) {
            for (size_t j = i + 1; j < m; ++j) {
                double q = static_cast<double>(m - 2) * D[i][j] - r[i] - r[j];
                if (q < min_q) {
                    min_q = q;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        // Branch lengths from best_i and best_j to new internal node
        double dij = D[best_i][best_j];
        double limb_i = (dij + (r[best_i] - r[best_j]) / static_cast<double>(m - 2)) / 2.0;
        double limb_j = dij - limb_i;
        if (limb_i < 0.0) limb_i = 0.0;
        if (limb_j < 0.0) limb_j = 0.0;

        // Form new node string
        std::ostringstream new_label;
        new_label << std::fixed << std::setprecision(4);
        new_label << "(" << nodes[best_i] << ":" << limb_i << "," << nodes[best_j] << ":" << limb_j << ")";

        // New distance row to other nodes: D(u, k) = (D(best_i, k) + D(best_j, k) - dij) / 2
        std::vector<double> new_row;
        new_row.reserve(m - 1);
        std::vector<std::string> new_nodes;
        new_nodes.reserve(m - 1);

        for (size_t k = 0; k < m; ++k) {
            if (k != best_i && k != best_j) {
                double duk = (D[best_i][k] + D[best_j][k] - dij) / 2.0;
                if (duk < 0.0) duk = 0.0;
                new_row.push_back(duk);
                new_nodes.push_back(nodes[k]);
            }
        }
        new_nodes.push_back(new_label.str());

        // Reconstruct D matrix
        size_t new_m = m - 1;
        std::vector<std::vector<double>> new_D(new_m, std::vector<double>(new_m, 0.0));
        for (size_t i = 0; i < new_m - 1; ++i) {
            size_t old_i = (i < best_i) ? i : ((i < best_j - 1) ? (i + 1) : (i + 2));
            for (size_t j = 0; j < new_m - 1; ++j) {
                size_t old_j = (j < best_i) ? j : ((j < best_j - 1) ? (j + 1) : (j + 2));
                new_D[i][j] = D[old_i][old_j];
            }
            new_D[i][new_m - 1] = new_row[i];
            new_D[new_m - 1][i] = new_row[i];
        }

        nodes = std::move(new_nodes);
        D = std::move(new_D);
    }

    // Connect final 2 nodes
    double final_d = D[0][1] / 2.0;
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(4);
    ss << "(" << nodes[0] << ":" << final_d << "," << nodes[1] << ":" << final_d << ");";
    return ss.str();
}

} // namespace sketch
