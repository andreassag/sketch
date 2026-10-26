#include "sketch/fasta.hpp"
#include "sketch/sketch.hpp"
#include "sketch/distance.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

void print_usage() {
    std::cerr << "Usage: sketch <command> [options]\n\n"
              << "Commands:\n"
              << "  sketch -i <fasta> -o <sketch.skh> [-k 21] [-s 11] [--scale 1] [--no-syncmer]\n"
              << "         Compute genome sketch and save to binary file\n"
              << "  dist <sketch1.skh> <sketch2.skh>\n"
              << "         Compute pairwise Jaccard, Containment, and Mash ANI\n"
              << "  matrix <s1.skh> <s2.skh> ...\n"
              << "         Compute all-vs-all distance and ANI matrix (TSV)\n"
              << "  tree <s1.skh> <s2.skh> ...\n"
              << "         Build Neighbor-Joining phylogenetic tree (Newick format)\n";
}

int cmd_sketch(int argc, char* argv[]) {
    std::string in_path, out_path;
    sketch::SketchParams params;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
            in_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            out_path = argv[++i];
        } else if ((arg == "-k") && i + 1 < argc) {
            params.k = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if ((arg == "-s") && i + 1 < argc) {
            params.s = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--scale" && i + 1 < argc) {
            params.scale = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--no-syncmer") {
            params.use_syncmer = false;
        }
    }

    if (in_path.empty() || out_path.empty()) {
        std::cerr << "Error: --input and --output are required.\n";
        print_usage();
        return 1;
    }

    params.t = (params.k > params.s) ? (params.k - params.s) / 2 : 0;

    auto records = sketch::read_fasta(in_path);
    if (records.empty()) {
        std::cerr << "Error: No FASTA records found in " << in_path << "\n";
        return 1;
    }

    // Concatenate or use primary record
    std::string full_seq;
    for (const auto& r : records) full_seq += r.seq;

    std::string sample_name = records[0].name;
    auto sk = sketch::GenomeSketch::compute(sample_name, full_seq, params);
    sk.save(out_path);

    std::cout << "[sketch] Processed " << sk.name << " (" << sk.total_bases << " bp) -> "
              << sk.hashes.size() << " hashes saved to " << out_path << "\n";
    return 0;
}

int cmd_dist(int argc, char* argv[]) {
    if (argc < 4) {
        print_usage();
        return 1;
    }
    auto sk1 = sketch::GenomeSketch::load(argv[2]);
    auto sk2 = sketch::GenomeSketch::load(argv[3]);

    auto res = sketch::compute_distance(sk1, sk2);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Query:       " << res.query_name << " (" << res.query_hashes << " hashes)\n";
    std::cout << "Reference:   " << res.ref_name << " (" << res.ref_hashes << " hashes)\n";
    std::cout << "Shared:      " << res.shared_hashes << "\n";
    std::cout << "Jaccard:     " << res.jaccard << "\n";
    std::cout << "Containment: " << res.containment << "\n";
    std::cout << "Mash Dist:   " << res.mash_dist << "\n";
    std::cout << "ANI (%):     " << (res.ani * 100.0) << " %\n";
    return 0;
}

int cmd_matrix(int argc, char* argv[]) {
    if (argc < 3) {
        print_usage();
        return 1;
    }

    std::vector<sketch::GenomeSketch> sketches;
    for (int i = 2; i < argc; ++i) {
        sketches.push_back(sketch::GenomeSketch::load(argv[i]));
    }

    size_t n = sketches.size();
    std::cout << "Sample";
    for (size_t i = 0; i < n; ++i) std::cout << "\t" << sketches[i].name;
    std::cout << "\n";

    std::cout << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < n; ++i) {
        std::cout << sketches[i].name;
        for (size_t j = 0; j < n; ++j) {
            auto res = sketch::compute_distance(sketches[i], sketches[j]);
            std::cout << "\t" << res.ani;
        }
        std::cout << "\n";
    }
    return 0;
}

int cmd_tree(int argc, char* argv[]) {
    if (argc < 3) {
        print_usage();
        return 1;
    }

    std::vector<sketch::GenomeSketch> sketches;
    std::vector<std::string> names;
    for (int i = 2; i < argc; ++i) {
        sketches.push_back(sketch::GenomeSketch::load(argv[i]));
        names.push_back(sketches.back().name);
    }

    size_t n = sketches.size();
    std::vector<std::vector<double>> D(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            auto res = sketch::compute_distance(sketches[i], sketches[j]);
            D[i][j] = res.mash_dist;
            D[j][i] = res.mash_dist;
        }
    }

    std::string newick = sketch::build_neighbor_joining_newick(names, D);
    std::cout << newick << "\n";
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage();
        return 1;
    }

    std::string cmd = argv[1];
    if (cmd == "sketch") return cmd_sketch(argc, argv);
    if (cmd == "dist") return cmd_dist(argc, argv);
    if (cmd == "matrix") return cmd_matrix(argc, argv);
    if (cmd == "tree") return cmd_tree(argc, argv);

    std::cerr << "Unknown command: " << cmd << "\n";
    print_usage();
    return 1;
}
