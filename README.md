# sketch: Rolling ntHash, Open Syncmers & FracMinHash ANI Engine

[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Zero Dependencies](https://img.shields.io/badge/dependencies-none-brightgreen.svg)]()

`sketch` is a high-performance, zero-dependency C++20 sequence sketching and Average Nucleotide Identity (ANI) estimation engine. It implements single-cycle 64-bit canonical **ntHash** rolling hashing, bounded-context **Open Syncmers** via an $O(1)$ amortized monotonic queue, 8-bit LSD radix sorting, and Mash/Poisson ANI estimation with Saitou-Nei Neighbor-Joining phylogenetic tree reconstruction.

---

## Key Features

1. **24-Byte Single-Cycle Canonical ntHash**:
   - $O(1)$ rolling hash per base using `std::rotl` and `std::rotr` with precomputed 64-bit random seeds.
   - Guaranteed **100% exact strand invariance**: canonical hash $\min(h_f, h_r)$ is identical for forward sequence and reverse complement without re-hashing.
2. **Open Syncmers ($k, s, t$) via Monotonic Queue**:
   - Fixed-size window minimum $s$-mer detection in $O(1)$ amortized time.
   - Preserves conserved core $k$-mers across mutated genomes better than standard minimizers.
3. **Cache-Friendly 8-Bit LSD Radix Sort**:
   - Sorts 64-bit hashes in 8 byte-level passes without heap thrashing, followed by linear deduplication.
4. **Distance & Phylogeny Engine**:
   - Two-pointer merge scan for Jaccard index $J$, Containment $C$, and Mash distance:
     $$\text{ANI} = \max\left(0.0, 1.0 + \frac{1}{k} \ln\left(\frac{2J}{1 + J}\right)\right)$$
   - Full all-vs-all distance matrix export (TSV) and Saitou-Nei Neighbor-Joining phylogenetic tree reconstruction in Newick format.
5. **Compact Binary Serialization (`.skh`)**:
   - Portable binary sketch file format with magic bytes (`SKCH\x01\x00`).

---

## Architectural Overview

```
FASTA Sequence ---> [ NtHasher (k, s) ] ---> [ MonotonicMinQueue ]
                            |                        |
                   Rolling 64-bit Hash        Open Syncmer Check (min at t)
                            |                        |
                            v                        v
                   [ FracMinHash Filter ] <--- [ Valid Syncmers ]
                            |
                   [ 8-Bit Radix Sort + Dedup ]
                            |
                            v
                   Binary Sketch File (.skh)
                            |
           +----------------+----------------+
           |                                 |
    Two-Pointer Merge                 Distance Matrix
           |                                 |
    Jaccard / Containment / ANI       Neighbor-Joining Tree (Newick)
```

---

## Build & Test

### Prerequisites
- Modern C++ compiler supporting C++20 (`g++` >= 11 or `clang++` >= 13)
- `make` (or `cmake` >= 3.20)

### Building via Makefile
```bash
# Build CLI binary
make all

# Run unit tests
make test
```

### Building via CMake
```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
ctest --output-on-failure
```

---

## CLI Usage

### 1. Compute Genome Sketch
```bash
./bin/sketch sketch -i genome.fa -o genome.skh -k 21 -s 11 --scale 1
```

### 2. Pairwise Distance & ANI
```bash
./bin/sketch dist strain_A.skh strain_B.skh
```
Output:
```text
Query:       Strain_A (1840 hashes)
Reference:   Strain_B (1832 hashes)
Shared:      1433
Jaccard:     0.640018
Containment: 0.778804
Mash Dist:   0.011801
ANI (%):     98.819909 %
```

### 3. All-vs-All ANI Matrix
```bash
./bin/sketch matrix strain_A.skh strain_B.skh strain_C.skh
```

### 4. Build Phylogenetic Tree (Newick)
```bash
./bin/sketch tree strain_A.skh strain_B.skh strain_C.skh > tree.nwk
```

---

## Verification & Benchmarks

Run automated verification:
```bash
./scripts/test_example.sh
./scripts/benchmark_compare.sh
```

### Performance on 100 kbp Genome:
- **Sketching Throughput**: ~10 ms total elapsed time (~10 Mbp/s single-threaded).
- **Pairwise Distance Query**: < 50 microseconds per comparison.
- **Memory Footprint**: Peak RSS < 7 MB.

---

## License

MIT License. Copyright (c) 2026 Andreas Sagen.
