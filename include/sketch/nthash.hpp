#pragma once

#include "dna.hpp"
#include <cstdint>
#include <bit>
#include <array>
#include <string_view>
#include <optional>
#include <vector>

namespace sketch {

// 64-bit random seeds for A, C, G, T
// Generated using 64-bit splitmix64 constants
inline constexpr std::array<uint64_t, 4> NTHASH_SEEDS = {
    0x3c8bf092713723bfULL, // A
    0x91849a6e2345ef01ULL, // C
    0xa4b1928424a1b029ULL, // G
    0x573b2a89364a2731ULL  // T
};

// Seed for reverse complement base:
// comp_seed(0) = T seed, comp_seed(1) = G seed, comp_seed(2) = C seed, comp_seed(3) = A seed
inline constexpr uint64_t get_seed(uint8_t code) noexcept {
    return (code < 4) ? NTHASH_SEEDS[code] : 0ULL;
}

inline constexpr uint64_t get_comp_seed(uint8_t code) noexcept {
    return (code < 4) ? NTHASH_SEEDS[3 - code] : 0ULL;
}

class NtHasher {
public:
    explicit NtHasher(uint32_t k) : k_(k) {}

    void reset() noexcept {
        h_f_ = 0;
        h_r_ = 0;
        window_.clear();
    }

    uint32_t k() const noexcept { return k_; }

    // Computes hash from scratch for a string of length k
    static std::optional<uint64_t> hash_kmer(std::string_view kmer) {
        if (kmer.size() == 0) return std::nullopt;
        uint32_t k = static_cast<uint32_t>(kmer.size());
        uint64_t hf = 0;
        uint64_t hr = 0;
        for (uint32_t i = 0; i < k; ++i) {
            uint8_t c = base_to_2bit(kmer[i]);
            if (c >= 4) return std::nullopt;
            uint64_t s_in = get_seed(c);
            uint64_t s_comp = get_comp_seed(c);
            hf = std::rotl(hf, 1) ^ s_in;
            hr ^= std::rotl(s_comp, i);
        }
        return std::min(hf, hr);
    }

    // Rolls the hasher by one character.
    // Returns canonical hash if a full valid k-mer window is present.
    std::optional<uint64_t> roll(char base) {
        uint8_t code = base_to_2bit(base);
        if (code >= 4) {
            // Invalid base (e.g. 'N'): reset window
            reset();
            return std::nullopt;
        }

        if (window_.size() < k_) {
            window_.push_back(code);
            uint64_t s_in = get_seed(code);
            uint64_t s_comp = get_comp_seed(code);
            h_f_ = std::rotl(h_f_, 1) ^ s_in;
            hr_init_ ^= std::rotl(s_comp, static_cast<int>(window_.size() - 1));

            if (window_.size() == k_) {
                h_r_ = hr_init_;
                return std::min(h_f_, h_r_);
            }
            return std::nullopt;
        }

        // Window is already of size k: roll out window_[0], roll in code
        uint8_t out_code = window_[window_head_];
        window_[window_head_] = code;
        window_head_ = (window_head_ + 1) % k_;

        uint64_t s_in = get_seed(code);
        uint64_t s_out = get_seed(out_code);
        uint64_t sc_in = get_comp_seed(code);
        uint64_t sc_out = get_comp_seed(out_code);

        // Forward roll
        h_f_ = std::rotl(h_f_, 1) ^ std::rotl(s_out, static_cast<int>(k_)) ^ s_in;

        // Reverse complement roll
        // h_r' = rotr(h_r, 1) ^ rotr(sc_out, 1) ^ rotl(sc_in, k - 1)
        h_r_ = std::rotr(h_r_, 1) ^ std::rotr(sc_out, 1) ^ std::rotl(sc_in, static_cast<int>(k_ - 1));

        return std::min(h_f_, h_r_);
    }

private:
    uint32_t k_;
    uint64_t h_f_ = 0;
    uint64_t h_r_ = 0;
    uint64_t hr_init_ = 0;
    std::vector<uint8_t> window_;
    size_t window_head_ = 0;
};

} // namespace sketch
