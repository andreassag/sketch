#pragma once

#include <cstdint>
#include <deque>
#include <utility>

namespace sketch {

template <typename T, typename IndexType = uint32_t>
class MonotonicMinQueue {
public:
    struct Entry {
        T value;
        IndexType index;
    };

    void push(T val, IndexType idx) {
        while (!deque_.empty() && deque_.back().value >= val) {
            deque_.pop_back();
        }
        deque_.push_back({val, idx});
    }

    void pop_before(IndexType min_valid_idx) {
        while (!deque_.empty() && deque_.front().index < min_valid_idx) {
            deque_.pop_front();
        }
    }

    const Entry& min() const {
        return deque_.front();
    }

    bool empty() const noexcept {
        return deque_.empty();
    }

    void clear() noexcept {
        deque_.clear();
    }

private:
    std::deque<Entry> deque_;
};

} // namespace sketch
