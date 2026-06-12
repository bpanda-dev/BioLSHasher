# BioLSHasher: Genomic Testing Framework for Locality Sensitive Hash Functions

![Logo](BioLSHasherLogo.png?raw=true)

## Overview

BioLSHasher is a C++ benchmark suite for evaluating how well a practical hash function fits the definition of a locality-senstive hash function in the context of genomics. 
A locality senstive hash (LSH) function is, informally, a hash function that maps similar elements to the same bucket with high probability and maps dissimilar elements to different buckets with low probability. 
Many hash functions attempt to mimick this property even if they are not provably LSH functions.
This benchmark suite gives developers and users a way to evaluate the extent to which their hash function exhibits LSH properties in practice. 

BioLSHasher by default has inbuilt support for five different similarity metrics: Hamming, Jaccard, angular, cosine and edit; but users can also define and implement custom metrics. Rather than testing LSH properties on random strings, BioLSHasher tests it on evolutionarily related strings, as is common in genomics. 
Concretely, it generates mutated strings from a random source string by introducing substitutions, deletions, and insertions, and 
then measures whether the hash function exhibits LSH properties on such string pairs.

To measure LSH properties, BioLSHasher performs two types of analysis. 
The first is a collision curve analysis, which shows the observed collision rates as a function of similarity.
The second is an analysis of false positive and false negative rates when used as part of nearest neighbour similarity search. 
Both of these analyses are combined with the effect of AND-OR amplification.
These analyses generate a dashboard of results, an example of which can be viewed [here](link needed). 

This README provides a brief guide to installing and running a simple test. The full documentation is provided at [biolshasher-guide.readthedocs.io](https://biolshasher-guide.readthedocs.io/en/latest/).

<!-- 
BioLSHasher is a specialized Locality Sensitive Hash function testing framework designed with biological context in mind. It is built upon the hash function testing framework [SMHasher3](https://gitlab.com/fwojcik/smhasher3). BioLSHasher currently supports two tests Collision Curve Test and Threshold-based Similarity Search Test, both of which support and-or amplification. -->

## Prerequisites
* [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/download)
  
## Getting Started

1. Clone the repository.

    ```bash
    git clone https://github.com/bpanda-dev/BioLSHasher.git
    cd BioLSHasher
    ```

2. Create conda environment from the provided `environment.yaml`.

      ```bash
       conda env create -f environment.yaml
      ```

      This creates a conda environment name `biolshasher`.

3. Activate the environment.

      ```bash
       conda activate biolshasher
      ```
    To later deactivate the biolshasher envirenoment, use "conda deactivate".
   
5. Build BioLSHasher.

   ```bash
    mkdir build
    cd build
    cmake ..
    make -j$(nproc)
   ```

6. Run an example test.

    The following analyzes a simple hash function (called exampleHash), which is an LSH for Hamming similarity.
    The run should complete in under a minute.
      
   ```bash
      ./BioLSHasher --test=LSHCollision,LSHApproxNearestNeighbour exampleHash --ncpu=4
   ```

     It should produce to output files, written to the `results/` directory under `BioLSHasher`:
    - `collisionResults_Hamming.outdata`
    - `ApproxNearestNeighbourResults_Hamming.outdata`


8. Generate standard CSVs, plots, and an interactive visualisation dahsboard.

    ```bash
        python ../analysis/generate_report.py \
        --coll=../results/collisionResults_Hamming.outdata \
        --ann=../results/approxNearestNeighbourResults_Hamming.outdata
    ```

   This writes all plots and processed CSVs to the `results/Hamming` directory.
   An HTML file is also generated, which aggregates all plots into a single interactive visualisation.

9. Open the dashboard to see the results.

   For example, if you use firefox as your browser, run
   
   ```bash
       firefox ../results/Hamming/Hamming_dashboard.html
   ```

  You should see a dashboard similar to the one [here](tODO)
  
## Acknowledgments

- BioLSHasher is built on top of the [SMHasher3](https://gitlab.com/fwojcik/smhasher3) suite. 
- Logo design: ChatGPT
- The project is inspired by the [paper](https://doi.org/10.1101/gr.281453.125) of Chen, Li, Shi, Shao, and Medvedev on "Hash functions in nucleotide sequence analysis".

## Citation

Paper describing BioLSHasher is forthcoming. For now, please cite
- Bikram Kumar Panda and Paul Medvedev, Benchmarking Locality-Sensitive Hash Functions for Genome Sequence Analysis using BioLSHasher, in progress.

## License

This project builds upon the SMHasher3 framework. Please refer to individual source files for specific licensing information.
