# BioLSHasher : Genomic Hash Function Testing Framework

## Overview

BioLSHasher is a specialized Locality Sensitive Hash function testing framework designed with biological context in mind. It is built upon the hash function testing framework [SMHasher3](https://gitlab.com/fwojcik/smhasher3). BioLSHasher currently supports two tests Collision Curve Test and Threshold-based Similarity Search Test, both of which support and-or amplification.

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

<!-- ## A brief description of tests included in BioLSHasher

### Test A:  Collision Curve

The Collision Curve test measures how often a hash function produces the same output for pairs of sequences at varying levels of similarity. BioLSHasher generates many sequence pairs, mutates one sequence in each pair at a controlled rate, and then checks whether the hash value of both sequences in the pair collide. By plotting the collision rate against the true similarity, we get a curve that shows how well the hash function distinguishes similar sequences from dissimilar ones. A good LSH should collide frequently for highly similar pairs and rarely for dissimilar pairs. The test also supports AND-OR amplification, which reshapes this curve by combining `b` hashes per band (AND-ing) and `r` independent bands (OR-ing). This lets us trace how different (b, r) configurations shift the effective similarity threshold and steepen the transition between the high-collision and low-collision regimes.

<img src="documentation/CollisionCurve.png" alt="Collision Curve" width="600"/>


### Test B:  c-ANN Test

The c-Approximate Nearest Neighbour test evaluates how well an LSH function performs when used as a hashing scheme for similarity search. BioLSHasher builds a database of reference sequences, creates mutated query sequences with known nearest neighbours, and then uses the hash function to retrieve candidates from the database. For each `(b, r)` configuration, it measures how many true neighbours were found (Recall) and how many of the returned candidates were incorrect (False Positive Rate). This tells you how the hash function will behave in a real similarity-search pipeline and helps you pick the right `(b, r)` parameters for your desired balance of accuracy and efficiency.
> **Why is this test needed?**:  While collision curve analysis characterizes the sensitivity of a hash family to pairwise similarity variations, it does not capture factors that matter in practice like the ability to correctly retrieve true neighbors from a database while minimizing false positive candidates. Minimizing false positives is important because real similarity-search pipelines include a re-scoring step using computationally expensive algorithms (e.g., full sequence alignment), and an inflated candidate set directly increases that cost.

<img src="documentation/c-ANN%20concept.png" alt="c-ANN concept" width="600"/>
todo:caption

<img src="documentation/ANN_flowchart.png" alt="ANN flowchart" width="600"/>
todo:caption -->

## Usage Guide for Adding a Novel Hash function for testing

This section walks through the complete workflow of: **adding a new hash function**, **running tests**, and **generating plots** from the results.

### Part 1 : Adding a New Hash Function

There are two approaches to add a new hash function to BioLSHasher: using the **interactive template generator script**, or **manually** using the [EXAMPLE_TEMPLATE.cpp](hashes/EXAMPLE_TEMPLATE.cpp) file.

#### Option A: Using `createHashTemplate.py` (Recommended)

BioLSHasher ships with an interactive Python script that scaffolds a new hash `.cpp` file and registers it in the build system automatically. Run it from the repository root:

```bash
python3 createHashTemplate.py
```

The script walks you through **10 guided steps**:


| Step | Prompt              | What it sets                                                                  |
| ---- | ------------------- | ----------------------------------------------------------------------------- |
| 1    | Hash Name           | C++ function name, `REGISTER_HASH` identifier, and output filename            |
| 2    | Author Name         | Copyright header                                                              |
| 3    | License             | License text in file header (MIT default, 8 options)                          |
| 4    | Family Name         | `REGISTER_FAMILY(...)` grouping (defaults to hash name)                       |
| 5    | Description         | Human-readable description in `REGISTER_HASH`                                 |
| 6    | Output Bits Size    | 32, 64, 128, 256, 512, 1024 or custom (multiple allowed)                      |
| 7    | LSH Candidacy       | Confirms the hash is an LSH candidate; exits if not (BioLSHasher is LSH-only)   |
| 8    | Similarity Name     | Built-in (`Hamming`, `Jaccard`, `Cosine`, `Angular`, `Edit`) or custom name   |
| 9    | Similarity Function | Auto-set for built-in metrics; prompts for a C++ function name if custom      |
| 10   | Hash Parameters     | (Optional) Add configurable parameters (e.g., k-mer size, window size)         |

Every input is validated for naming rules and C++ keyword checks. The script will re-prompt if bad input is provided. If the hash is not an LSH candidate,the script stops with a message that BioLSHasher only supports LSH-related tests.

<!-- **Parameters (Step 10):**
- Optional step to add hash function parameters.
- Each parameter should have unique name.
- Parameters are automatically generated as `#define` macros at the top of the file (e.g., `#define k 10`)
- The `$.parameterValues` in `REGISTER_HASH` references the macro names, allowing easy customization by editing the `#define` values at the top of the file.
-->

**It generates:**
1. A compilable C++ template file at `hashes/<hashname>.cpp` containing:
    - Copyright header with your chosen license
    - A hash function stub for each selected output bit size
    - The similarity function implementation (included automatically for built-in metrics like Hamming or Edit; a stub for custom metrics)
    - A default equality checker for hash outputs, defining the condition under which two outputs produced by the hash function are considered identical.
    - Macros defined for parameters of the hash function.
2. An updated `hashes/Hashsrc.cmake` with the new hash registered.

> **Full documentation:** See [`createHashTemplate.md`](documentation/createHashTemplate.md) for the complete reference including validation rules, example sessions, and troubleshooting.

#### Option B: Manual Creation

[TODO]

#### After Creating the Template

Regardless of the methods used above:

**1. Implement your hash logic** : open `hashes/<hashname>.cpp` and replace the placeholder body:

```cpp

static void MyHash( const void * in, const size_t len, const seed_t seed, void * out ) {
    // Your hash implementation goes here
    // 'in' = pointer to input data (genomic sequence, unencoded)
    // 'len' = input length in bytes
    // 'seed' = seed value
    uint32_t hash = 0;
    // ... your logic ...
    PUT_U32(hash, (uint8_t *)out, 0);
}
```

**2. Set the correct hash flags** in `REGISTER_HASH(...)`:

| Flag                                     | When to use                                           | FLAG TYPE |
| ---------------------------------------- | ----------------------------------------------------- | --------- |
| `FLAG_HASH_LOCALITY_SENSITIVE`           | Your hash is an LSH function                          | HASH      |
| `FLAG_IMPL_SLOW`                         | Hash is computationally expensive                     | IMPL      |
| `FLAG_IMPL_VERY_SLOW`                    | Very slow (reduces LSH test parameters automatically) | IMPL      |
<!-- | `FLAG_IMPL_SMALL_SEQUENCE_LENGTH`        | Uses small sequences (40 bases instead of 512)        | IMPL      | -->

- **IMPL flags** are for controlling **test execution behaviour** : they tell the testing module how to adjust parameters (e.g., fewer iterations, shorter sequences) based on the computational cost of your hash implementation. They do not affect hash semantics.
- **HASH flags** are for declaring the **mathematical properties** of your hash function. They tell the testing module whether the hash function is an LSH candidate or not.

#### Example: Adding Flags in `REGISTER_HASH(...)`

If your hash is a locality-sensitive MinHash that preserves Jaccard similarity using overlapping 5-mers, and it's computationally expensive:

```cpp
REGISTER_HASH(MyMinHash_64,
   $.desc            = "My custom MinHash (64-bit, Jaccard, overlapping k-mers)",
   $.hash_flags      = FLAG_HASH_LOCALITY_SENSITIVE,
   $.impl_flags      = FLAG_IMPL_VERY_SLOW,
   // ... other registration fields ...
);
```
Note that we use ORing of the flags to combine them.

> **Tip:** If your hash is extremely slow (minutes per test), use `FLAG_IMPL_VERY_SLOW` which _significantly_ reduces the number of sequence pairs to work on.

**3. Build and verify:**

```bash
cd build
cmake ..
make 

# Confirm your hash is registered. The following command should return the name of your Hash.
./BioLSHasher --list | grep MyHash
```

---

### Part 2 : Running Tests

#### CLI Reference

```bash
./BioLSHasher [options] [<hashname>]
```

**Key arguments:**

| Arguments               | Description                                                   | Status in BioLSHasher     |
|-------------------------| ------------------------------------------------------------- |-------------------------|
| `--list`                | List all registered hashes with descriptions                  | Active                  |
| `--listnames`           | List just hash names (useful for scripting)                   | Active                  |
| `--tests`               | Print all valid test suite names                              | Active                  |
| `--test=<name>[,...]`   | Enable **only** these tests (disables all others first)       | Active                  |
| `--notest=<name>[,...]` | Disable specific tests                                        | In Progress             |
| `--ncpu=N`              | Number of threads for parallel tests                          | Active |
| `--verbose`             | Verbose output with more stats and diagrams                   | Not Active              |
| `--version`             | Print version string                                          | Active                  |

#### Test 1 : LSH Collision Test (`--test=LSHCollision`)
The Collision Curve test measures how often a hash function produces the same output for pairs of sequences at varying levels of similarity. 
<!--- The `(AND, OR)` grid is configured in [`lib/LSHGlobals.cpp`](lib/LSHGlobals.cpp) via `g_ANN_start_B`, `g_ANN_MAX_B`, `g_ANN_start_R`, `g_ANN_MAX_R` variables. --->

> **Full documentation:** See [`Collisiontest.md`](documentation/CollisionTest.md) for the complete reference including AND-OR basics, internal pipeline, pseudocode, all configurable parameters, and caveats.

```bash
# Run LSH collision test (multi-threaded recommended)
./BioLSHasher --test=LSHCollision exampleHash --ncpu=4
```

**Output:** `results/collisionResults_<similarityname>.csv`  (`similarityname` is the similarity metric the `<hashname>` is tested against.)

<!--
**Test parameters** (automatically adjusted based on hash flags):

| Parameter                 | Normal Hash | `FLAG_IMPL_VERY_SLOW`                            |
| ------------------------- | ----------- | ------------------------------------------------ |
| Sequence pairs per bin    | 50,000      | 5,000                                            |
| Hash repetitions per pair | 2,000       | 2,000                                            |
| Sequence length           | 512 bases   | 512 (or 40 if `FLAG_IMPL_SMALL_SEQUENCE_LENGTH`) |
-->
---

#### Test 2 : Approximate Nearest Neighbour Test (`--test=LSHApproxNearestNeighbour`)

This test evaluates the hash function as an indexing scheme for Similarity search applications, which is an important use case in many genomic LSH pipelines.

<!--
It:

1. Generates a reference database of random sequences.
2. Samples query sequences from the reference, then mutates them to target 90–100% similarity.
3. Computes brute-force ground truth (the true nearest neighbour for each query).
4. For each `(b, r)` configuration, builds an LSH index with `r` tables (each using a `b`-hash band signature), inserts all reference sequences, queries with each mutated sequence, and compares results against ground truth.
5. Reports **Recall**, **Precision**, **FPR**, and **F1** averaged over multiple independent runs.

This helps you select the optimal `(b, r)` parameters for your application by directly measuring retrieval quality.
-->

> **Full documentation:** See [ApproximateNearestNeighbour.md](documentation/ApproximateNearestNeighbour.md) for the complete reference including the 5-phase pipeline internals, all configurable parameters, evaluation metrics, and caveats.

```bash
./BioLSHasher --test=LSHApproxNearestNeighbour exampleHash --ncpu=4
```

**Output:** `results/ApproxNearestNeighbourResults_<similarityname>.csv` (`<similarityname>` is the similarity metric the `<hashname>` is tested against.)

<!--
**Output format:**

```
:1: LSH Approx Nearest Neighbour Summary
:2: Hashname, SequenceLength, TokenLength, Distance Metric, Mutation Model, Mutation Expression
:3: <values>
:4.1:/:4.2:/:4.3: Hash-specific parameters
:5: b, r, Avg_Recall, Avg_Precision, Avg_FPR, Avg_F1_Score   (column headers)
:6: <one data row per (b,r) pair>
```
--->

> **Note**: Each run **appends** to the CSV if it already exists, so you can accumulate results across multiple token lengths or configurations.

---

#### Configuring the Mutation Model

The mutation model is set at compile time in [`lib/LSHGlobals.cpp`](lib/LSHGlobals.cpp):

| Variable                     | Options                                                                                                                                                                            | Default         |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------- |
| `g_mutation_model`           | `MUTATION_MODEL_SIMPLE_SNP_ONLY` (0), `MUTATION_MODEL_GEOMETRIC_MUTATOR` (1)                                                                                                       | `1` (Geometric) |
| `g_mutation_expression_type` | `MUTATION_EXPRESSION_BALANCED` (0), `MUTATION_EXPRESSION_SUB_ONLY` (1), `MUTATION_EXPRESSION_DEL_LITE` (2), `MUTATION_EXPRESSION_INS_LITE` (3), `MUTATION_EXPRESSION_SUB_LITE` (4) | `0` (Balanced)  |

> **Note:** If a hash has hamming similarity metric, please ensure that you only use substitution only mutation model, i.e. `MUTATION_MODEL_SIMPLE_SNP_ONLY`. 

See [`MutationModels.md`](documentation/MutationModels.md) for full documentation on what each mutation expression does.

After changing these values, rebuild:

```bash
cd build && make -j$(nproc)
```

---

### Part 3 : Generating Reports and Visualisations

To process your test results, generate plots, and create an interactive HTML dashboard, use the `generate_report.py` script:

```bash
# Ensure that you have activated the biolshasher conda environment and are in the root directory of BioLSHasher. At least one of --coll or --ann args must be provided to generate the plots and visualisations.
python analysis/generate_report.py \
    --coll=results/collisionResults_<similarityname>.csv \
    --ann=results/ApproxNearestNeighbourResults_<similarityname>.csv
```

<!-- All plots are saved in the current working directory.

| Filename Pattern              | Description                                 |
| ----------------------------- | ------------------------------------------- |
| `*_binaveraged.png`           | Bin-averaged collision curves (overlay)     |
| `*_snpcolor_multiplot.png`    | SNP-rate colored scatter subplots           |
| `*_boxplot_multiplot.png`     | Box plot subplots per similarity bin        |
| `*_monocolor_multiplot.png`   | Single-color scatter subplots               |
| `*_verificationCurves.png`    | Mutation parameter distributions (sanity check) |

#### 3b. ANN Result Plots (`analysis/plot_ANN.py`)

Visualises the Recall vs FPR trade-off across `(b, r)` configurations.

```bash
  # ensure that you have activate biolshasher conda environment and are in the root directory of biolshasher.  
  python analysis/plot_ANN.py results/ApproxNearestNeighbourResults_<hashname>.csv
```

Produces **6 plots** (3 linear + 3 log-scale):

| Filename Pattern                       | Description                                              |
| -------------------------------------- | -------------------------------------------------------- |
| `*_fpr_vs_recall.png` / `*_log.png`   | FPR vs Recall per `b` value, annotated with `r`          |
| `*_best_fpr_per_recall_bin.png` / `*_log.png` | Pareto-optimal (min FPR) per recall bin           |
| `*_all_points.png` / `*_log.png`      | All `(b,r)` configs as an annotated scatter              | -->

---

### Plotting Error

If you encounter an error similar to the following while generating reports and plots:

```bash
--- Processing Collision Curve data: ../results/collisionResults_Hamming.csv ---
Reading: ../results/collisionResults_Hamming.csv
Traceback (most recent call last):
  File "~/Hashing/BioLSHasher/analysis/plot_collisioncurve.py", line 1208, in <module>
    main()
  File "~/Hashing/BioLSHasher/analysis/plot_collisioncurve.py", line 1079, in main
    sections, df = read_collision_data_complete(args.csvfile)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "~/Hashing/BioLSHasher/analysis/plot_collisioncurve.py", line 102, in read_collision_data_complete
    collision_rates = [float(x.strip()) for x in line_content.split(',')]
                       ^^^^^^^^^^^^^^^^
ValueError: could not convert string to float: ''
Error occurred while plotting collision curve.
```

this usually indicates that the output CSV file has been corrupted or partially overwritten during either the current run or during previous experiment(s). This is applicable for both: Collision Test and Similarity Search Test.

**Solution:** Delete the existing results CSV file and rerun the experiment so that a clean output file is generated.

For example:

```bash
rm ../results/collisionResults_Hamming.csv
```

Then rerun the experiment and regenerate the plots.

---

## License

This project builds upon the SMHasher3 framework. Please refer to individual source files for specific licensing information.

## References

- [SMHasher3](https://gitlab.com/fwojcik/smhasher3): Original hash testing framework biolshasher is based on.
- Hash functions in genomic sequence analysis - Ke Chen, Xiang Li, Qian Shi, Mingfu Shao, Paul Medvedev
- Shrivastava, Anshumali, and Ping Li. "In defense of minhash over simhash." Artificial intelligence and statistics. PMLR, 2014.
