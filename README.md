# BioLSHasher : Genomic Hash Function Testing Framework

<!-- ![Logo](https://github.com/bpanda-dev/biolshasher/blob/main/BioLSHasherLogo.png?raw=true) -->

![Logo](BioLSHasherLogo.png?raw=true)

## Overview

BioLSHasher is a C++ testing framework for evaluating locality-sensitive hash (LSH) functions in the context of biological sequence analysis. Built on top of the [SMHasher3](https://gitlab.com/fwojcik/smhasher3) testing harness, it provides a complete workflow to add and register custom hash implementations, run mutation-aware benchmarks (collision-curve analysis and threshold based near neighbour search performance evaluation), and produce CSV results along with interactive plots.

Use BioLSHasher to measure how well LSH families preserve biological similarity, compare retrieval quality across (b,r) LSH parameterizations, and select practical hashing strategies for genomic applications.

<!-- 
BioLSHasher is a specialized Locality Sensitive Hash function testing framework designed with biological context in mind. It is built upon the hash function testing framework [SMHasher3](https://gitlab.com/fwojcik/smhasher3). BioLSHasher currently supports two tests Collision Curve Test and Threshold-based Similarity Search Test, both of which support and-or amplification. -->

> **Full documentation:** [biolshasher-guide.readthedocs.io](https://biolshasher-guide.readthedocs.io/en/latest/)


## Getting Started

To get started with BioLSHasher, follow these steps:

1. **Clone the repository**:

    ```bash
    git clone https://github.com/bpanda-dev/BioLSHasher.git
    cd BioLSHasher
    ```

2. **Set up conda environment**: (BioLSHasher is bundled with multiple python scripts for plotting and for aiding users to connect their hash function to the BioLSHasher framework.)
   1. Prerequisites: [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/download)
   2. Create the conda environment from the provided `environment.yaml`:

      ```bash
       conda env create -f environment.yaml
      ```

      This creates a conda environment name `biolshasher`.
   3. Activate the environment:

      ```bash
       conda activate biolshasher
       # to deactivate biolshasher, use "conda deactivate"
      ```

3. **Build BioLSHasher**:

   ```bash
    mkdir build
    cd build
    cmake ..
    make -j$(nproc)
   ```

4. **Run an example test**:

   ```bash
    # A sample test run of the two analyses: collision and similarity search on the included
    # example hash function (exampleHash),which preserves Hamming distance.
    ./BioLSHasher --test=LSHCollision,LSHApproxNearestNeighbour exampleHash --ncpu=4
   ```

    The run should complete in under a minute. If everything works correctly, two output files are written to the `results/` directory under `BioLSHasher`:
    - `collisionResults_Hamming.csv`
    - `ApproxNearestNeighbourResults_Hamming.csv`

    > **Note:**
    > 1. The filename reflects the similarity metric the hash function operates on, rather than the hash function name.
    > 2. These files use a custom CSV format and require post-processing before they can be used as standard CSVs.

5. **Process the output files** to generate standard CSVs, plots, and an interactive visualisation:

    ```bash
    python ../analysis/generate_report.py \
        --coll=../results/collisionResults_Hamming.csv \
        --ann=../results/approxNearestNeighbourResults_Hamming.csv
    ```

    This writes all plots and processed CSVs to the `results/` directory under `Hamming` name. An HTML file is also generated, which aggregates all plots into a single interactive visualisation. If you use firefox, you can visualize it using ```firefox ../results/Hamming/Hamming_dashboard.html```
    
   
---

## License

This project builds upon the SMHasher3 framework. Please refer to individual source files for specific licensing information.

## References

- [SMHasher3](https://gitlab.com/fwojcik/smhasher3): Original hash testing framework biolshasher is based on.
- Hash functions in genomic sequence analysis - Ke Chen, Xiang Li, Qian Shi, Mingfu Shao, Paul Medvedev
- Shrivastava, Anshumali, and Ping Li. "In defense of minhash over simhash." Artificial intelligence and statistics. PMLR, 2014.
