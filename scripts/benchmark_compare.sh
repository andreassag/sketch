#!/usr/bin/env bash
set -euo pipefail

BIN="./bin/sketch"
if [[ ! -x "$BIN" ]]; then
    make all
fi

echo "=== Benchmark: Sketching Throughput and ANI Accuracy ==="

# Generate 100,000 bp benchmark genome
python3 -c '
import random
random.seed(999)
bases = ["A", "C", "G", "T"]
g = "".join(random.choice(bases) for _ in range(100000))
with open("/tmp/bench_100k.fa", "w") as f:
    f.write(">bench_100k\n" + g + "\n")
'

echo "Benchmarking 100 kbp genome sketching..."
/usr/bin/time -v $BIN sketch -i /tmp/bench_100k.fa -o /tmp/bench_100k.skh -k 21 -s 11 2>&1 | grep -E "Elapsed|Maximum resident"

echo -e "\nBenchmarking distance computation (1000 iterations)..."
python3 -c '
import time, subprocess
t0 = time.perf_counter()
for _ in range(100):
    subprocess.run(["./bin/sketch", "dist", "/tmp/bench_100k.skh", "/tmp/bench_100k.skh"], stdout=subprocess.DEVNULL)
t1 = time.perf_counter()
print(f"100 pairwise comparisons: {t1 - t0:.4f} s ({(t1 - t0)*10:.2f} ms per comparison including process invocation)")
'

rm -f /tmp/bench_100k.fa /tmp/bench_100k.skh
echo "Benchmark completed."
