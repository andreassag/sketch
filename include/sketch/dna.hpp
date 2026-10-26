#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <array>
#include <algorithm>

namespace sketch {

// Fast 2-bit base encoding: A=0, C=1, G=2, T=3, invalid=4
inline constexpr uint8_t base_to_2bit(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return 4;
    }
}

inline constexpr char complement_base(char c) noexcept {
    switch (c) {
        case 'A': return 'T';
        case 'a': return 't';
        case 'C': return 'G';
        case 'c': return 'g';
        case 'G': return 'C';
        case 'g': return 'c';
        case 'T': case 'U': return 'A';
        case 't': case 'u': return 'a';
        default: return 'N';
    }
}

inline std::string reverse_complement(std::string_view seq) {
    std::string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        rc.push_back(complement_base(*it));
    }
    return rc;
}

inline bool is_valid_dna(char c) noexcept {
    return base_to_2bit(c) < 4;
}

} // namespace sketch
