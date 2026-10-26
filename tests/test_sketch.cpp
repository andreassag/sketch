#include "sketch/dna.hpp"
#include "sketch/nthash.hpp"
#include "sketch/monotonic_queue.hpp"
#include "sketch/radix_sort.hpp"
#include "sketch/sketch.hpp"
#include "sketch/distance.hpp"
#include <iostream>
#include <cassert>
#include <random>
#include <sstream>

void test_nthash_rolling_and_invariance() {
    std::string seq = "ACGTACGTTTGCAACCGTTAGCTAGCTAGCTGATCGATCGATCGATC";
    uint32_t k = 15;

    sketch::NtHasher hasher(k);
    std::vector<uint64_t> rolling_hashes;

    for (char c : seq) {
        auto h = hasher.roll(c);
        if (h.has_value()) {
            rolling_hashes.push_back(*h);
        }
    }

    assert(rolling_hashes.size() == seq.size() - k + 1);

    // Verify rolling matches scratch hash
    for (size_t i = 0; i < rolling_hashes.size(); ++i) {
        std::string_view kmer(seq.data() + i, k);
        auto direct_h = sketch::NtHasher::hash_kmer(kmer);
        assert(direct_h.has_value());
        assert(*direct_h == rolling_hashes[i]);
    }

    // Verify strand invariance: reverse complement k-mer has same canonical hash
    for (size_t i = 0; i < rolling_hashes.size(); ++i) {
        std::string kmer = seq.substr(i, k);
        std::string rc = sketch::reverse_complement(kmer);
        auto direct_f = sketch::NtHasher::hash_kmer(kmer);
        auto direct_rc = sketch::NtHasher::hash_kmer(rc);
        assert(direct_f == direct_rc);
    }

    std::cout << "[PASS] test_nthash_rolling_and_invariance\n";
}

void test_monotonic_queue() {
    sketch::MonotonicMinQueue<int, uint32_t> mq;
    std::vector<int> vals = {5, 3, 8, 2, 9, 1, 4, 7};

    // Push first 3 (5, 3, 8) -> min is 3 at idx 1
    mq.push(5, 0);
    mq.push(3, 1);
    mq.push(8, 2);
    assert(mq.min().value == 3 && mq.min().index == 1);

    // Push 2 at idx 3 -> min is 2 at idx 3
    mq.push(2, 3);
    assert(mq.min().value == 2 && mq.min().index == 3);

    // Pop before idx 2 (window [2, 3]) -> 2 is still min
    mq.pop_before(2);
    assert(mq.min().value == 2 && mq.min().index == 3);

    // Pop before idx 4 -> pop 2
    mq.push(9, 4);
    mq.pop_before(4);
    assert(mq.min().value == 9 && mq.min().index == 4);

    std::cout << "[PASS] test_monotonic_queue\n";
}

void test_radix_sort() {
    std::mt19937_64 rng(42);
    std::vector<uint64_t> vec(1000);
    for (auto& v : vec) v = rng();

    auto std_vec = vec;
    std::sort(std_vec.begin(), std_vec.end());

    sketch::radix_sort_u64(vec);
    assert(vec == std_vec);

    std::cout << "[PASS] test_radix_sort\n";
}

void test_sketch_ani() {
    // Generate synthetic genome of 5000 bp
    std::mt19937 rng(1234);
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::string genome1;
    for (int i = 0; i < 5000; ++i) {
        genome1 += bases[rng() % 4];
    }

    // Identical sketch should give ANI = 1.0 (100%)
    sketch::SketchParams params;
    params.k = 21;
    params.s = 11;
    params.t = 5;
    params.use_syncmer = true;

    auto sk1 = sketch::GenomeSketch::compute("g1", genome1, params);
    auto sk1_copy = sketch::GenomeSketch::compute("g1_copy", genome1, params);

    auto res_id = sketch::compute_distance(sk1, sk1_copy);
    assert(res_id.shared_hashes == sk1.hashes.size());
    assert(res_id.jaccard == 1.0);
    assert(res_id.ani == 1.0);

    // Introduce 2% random mutations into genome2
    std::string genome2 = genome1;
    for (size_t i = 0; i < genome2.size(); i += 50) {
        genome2[i] = bases[(sketch::base_to_2bit(genome2[i]) + 1) % 4];
    }

    auto sk2 = sketch::GenomeSketch::compute("g2", genome2, params);
    auto res_mut = sketch::compute_distance(sk1, sk2);

    assert(res_mut.shared_hashes > 0);
    assert(res_mut.jaccard < 1.0 && res_mut.jaccard > 0.3);
    assert(res_mut.ani >= 0.95 && res_mut.ani <= 0.99);

    std::cout << "[PASS] test_sketch_ani: mutated ANI = " << (res_mut.ani * 100.0) << "%\n";
}

void test_serialization() {
    std::string genome = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC";
    sketch::SketchParams params;
    params.k = 15;
    params.s = 9;
    params.t = 3;

    auto sk = sketch::GenomeSketch::compute("test_sample", genome, params);
    std::string tmp_file = "/tmp/test_sketch_io.skh";
    sk.save(tmp_file);

    auto loaded = sketch::GenomeSketch::load(tmp_file);
    assert(loaded.name == sk.name);
    assert(loaded.params.k == sk.params.k);
    assert(loaded.params.s == sk.params.s);
    assert(loaded.hashes == sk.hashes);

    std::cout << "[PASS] test_serialization\n";
}

void test_nj_tree() {
    std::vector<std::string> taxa = {"A", "B", "C"};
    std::vector<std::vector<double>> D = {
        {0.0, 0.2, 0.4},
        {0.2, 0.0, 0.4},
        {0.4, 0.4, 0.0}
    };
    std::string newick = sketch::build_neighbor_joining_newick(taxa, D);
    assert(!newick.empty() && newick.back() == ';');
    assert(newick.find("A:") != std::string::npos);
    assert(newick.find("B:") != std::string::npos);
    assert(newick.find("C:") != std::string::npos);

    std::cout << "[PASS] test_nj_tree: " << newick << "\n";
}

int main() {
    std::cout << "Running sketch unit tests...\n";
    test_nthash_rolling_and_invariance();
    test_monotonic_queue();
    test_radix_sort();
    test_sketch_ani();
    test_serialization();
    test_nj_tree();
    std::cout << "All sketch unit tests passed successfully!\n";
    return 0;
}
