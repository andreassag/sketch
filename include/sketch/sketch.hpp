#pragma once

#include "nthash.hpp"
#include "monotonic_queue.hpp"
#include "radix_sort.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <cmath>

namespace sketch {

struct SketchParams {
    uint32_t k = 21;
    uint32_t s = 11;
    uint32_t t = 5;       // (k - s) / 2
    uint32_t scale = 1;   // FracMinHash scale factor (1 = keep all syncmers, >1 = fraction)
    bool use_syncmer = true;
};

class GenomeSketch {
public:
    std::string name;
    SketchParams params;
    uint64_t total_bases = 0;
    std::vector<uint64_t> hashes;

    GenomeSketch() = default;

    static GenomeSketch compute(const std::string& name, const std::string& seq, const SketchParams& params) {
        GenomeSketch sk;
        sk.name = name;
        sk.params = params;
        sk.total_bases = seq.size();

        if (seq.size() < params.k) return sk;

        uint64_t max_hash = (params.scale <= 1) ? UINT64_MAX : (UINT64_MAX / params.scale);

        if (!params.use_syncmer) {
            // Pure FracMinHash on canonical k-mers
            NtHasher hasher(params.k);
            for (char c : seq) {
                auto h = hasher.roll(c);
                if (h.has_value() && *h <= max_hash) {
                    sk.hashes.push_back(*h);
                }
            }
        } else {
            // Open Syncmers: (k, s, t)
            // A k-mer of length k contains w = (k - s + 1) s-mers.
            // The s-mer at offset t in the k-mer window must be the minimum s-mer.
            uint32_t w = params.k - params.s + 1;
            NtHasher s_hasher(params.s);
            NtHasher k_hasher(params.k);

            // Ring buffer of s-mer hashes within current k-mer
            std::vector<uint64_t> s_hashes;
            s_hashes.reserve(seq.size());

            // First pass: extract all canonical s-mer hashes along sequence with positions
            std::vector<std::pair<uint32_t, uint64_t>> valid_smers; // (pos, hash)
            for (size_t i = 0; i < seq.size(); ++i) {
                auto h = s_hasher.roll(seq[i]);
                if (h.has_value()) {
                    valid_smers.push_back({static_cast<uint32_t>(i + 1 - params.s), *h});
                }
            }

            // Also compute k-mer hashes
            // A k-mer starting at pos has s-mers starting at pos, pos+1, ..., pos + w - 1.
            // If all w s-mers are present and the minimum is at pos + t, it's an open syncmer.
            // Monotonic queue over valid_smers
            MonotonicMinQueue<uint64_t, uint32_t> mq;
            size_t s_idx = 0;

            for (size_t i = 0; i < seq.size(); ++i) {
                auto kh = k_hasher.roll(seq[i]);
                if (kh.has_value()) {
                    uint32_t kmer_start = static_cast<uint32_t>(i + 1 - params.k);
                    // Add s-mers up to kmer_start + w - 1
                    while (s_idx < valid_smers.size() && valid_smers[s_idx].first <= kmer_start + w - 1) {
                        mq.push(valid_smers[s_idx].second, valid_smers[s_idx].first);
                        s_idx++;
                    }
                    // Pop s-mers before kmer_start
                    mq.pop_before(kmer_start);

                    if (!mq.empty()) {
                        auto min_entry = mq.min();
                        if (min_entry.index == kmer_start + params.t) {
                            if (*kh <= max_hash) {
                                sk.hashes.push_back(*kh);
                            }
                        }
                    }
                }
            }
        }

        sort_and_dedup(sk.hashes);
        return sk;
    }

    void save(const std::string& filepath) const {
        std::ofstream out(filepath, std::ios::binary);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot write to sketch file: " + filepath);
        }

        const char magic[4] = {'S', 'K', 'C', 'H'};
        uint16_t version = 1;
        out.write(magic, 4);
        out.write(reinterpret_cast<const char*>(&version), sizeof(version));

        uint32_t name_len = static_cast<uint32_t>(name.size());
        out.write(reinterpret_cast<const char*>(&name_len), sizeof(name_len));
        out.write(name.data(), name_len);

        out.write(reinterpret_cast<const char*>(&params.k), sizeof(params.k));
        out.write(reinterpret_cast<const char*>(&params.s), sizeof(params.s));
        out.write(reinterpret_cast<const char*>(&params.t), sizeof(params.t));
        out.write(reinterpret_cast<const char*>(&params.scale), sizeof(params.scale));
        uint8_t use_sync = params.use_syncmer ? 1 : 0;
        out.write(reinterpret_cast<const char*>(&use_sync), sizeof(use_sync));
        out.write(reinterpret_cast<const char*>(&total_bases), sizeof(total_bases));

        uint64_t num_hashes = hashes.size();
        out.write(reinterpret_cast<const char*>(&num_hashes), sizeof(num_hashes));
        if (num_hashes > 0) {
            out.write(reinterpret_cast<const char*>(hashes.data()), num_hashes * sizeof(uint64_t));
        }
    }

    static GenomeSketch load(const std::string& filepath) {
        std::ifstream in(filepath, std::ios::binary);
        if (!in.is_open()) {
            throw std::runtime_error("Cannot read sketch file: " + filepath);
        }

        char magic[4];
        in.read(magic, 4);
        if (std::memcmp(magic, "SKCH", 4) != 0) {
            throw std::runtime_error("Invalid sketch magic bytes in: " + filepath);
        }

        uint16_t version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));
        if (version != 1) {
            throw std::runtime_error("Unsupported sketch version: " + std::to_string(version));
        }

        GenomeSketch sk;
        uint32_t name_len = 0;
        in.read(reinterpret_cast<char*>(&name_len), sizeof(name_len));
        sk.name.resize(name_len);
        in.read(&sk.name[0], name_len);

        in.read(reinterpret_cast<char*>(&sk.params.k), sizeof(sk.params.k));
        in.read(reinterpret_cast<char*>(&sk.params.s), sizeof(sk.params.s));
        in.read(reinterpret_cast<char*>(&sk.params.t), sizeof(sk.params.t));
        in.read(reinterpret_cast<char*>(&sk.params.scale), sizeof(sk.params.scale));
        uint8_t use_sync = 0;
        in.read(reinterpret_cast<char*>(&use_sync), sizeof(use_sync));
        sk.params.use_syncmer = (use_sync != 0);
        in.read(reinterpret_cast<char*>(&sk.total_bases), sizeof(sk.total_bases));

        uint64_t num_hashes = 0;
        in.read(reinterpret_cast<char*>(&num_hashes), sizeof(num_hashes));
        sk.hashes.resize(num_hashes);
        if (num_hashes > 0) {
            in.read(reinterpret_cast<char*>(sk.hashes.data()), num_hashes * sizeof(uint64_t));
        }
        return sk;
    }
};

} // namespace sketch
