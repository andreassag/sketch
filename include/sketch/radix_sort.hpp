#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include <algorithm>

namespace sketch {

inline void radix_sort_u64(std::vector<uint64_t>& vec) {
    if (vec.size() <= 1) return;

    std::vector<uint64_t> temp(vec.size());
    std::array<size_t, 256> count;

    for (int shift = 0; shift < 64; shift += 8) {
        count.fill(0);

        for (uint64_t v : vec) {
            uint8_t byte = static_cast<uint8_t>((v >> shift) & 0xFF);
            count[byte]++;
        }

        std::array<size_t, 256> offset;
        offset[0] = 0;
        for (size_t i = 1; i < 256; ++i) {
            offset[i] = offset[i - 1] + count[i - 1];
        }

        for (uint64_t v : vec) {
            uint8_t byte = static_cast<uint8_t>((v >> shift) & 0xFF);
            temp[offset[byte]++] = v;
        }

        vec.swap(temp);
    }
}

inline void sort_and_dedup(std::vector<uint64_t>& vec) {
    radix_sort_u64(vec);
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

} // namespace sketch
