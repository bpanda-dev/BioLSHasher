#include "specifics.h"
#include "HashInfo.h"
#include "TestGlobals.h"
#include "Instantiate.h"
#include "coolstuff.h"

#include "BioDataGeneration.h"
#include "LSHGlobals.h"
#include "ApproxNearestNeighbour.h"

#include <random>
#include <fstream>
#include <iostream>
#include <filesystem>

#include <numeric>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <set>



// #if defined(HAVE_THREADS)
// #include <atomic>
// #endif
#include <atomic>
#include <chrono>
#include <thread>
#include "assertMsg.h"


// ============================================================================
// ANN LSH Index Data Structure
// ============================================================================
class ANN_LSH_Index {
private:
    uint32_t numTables_r;
    uint32_t hashesPerTable_b;
    HashFn hashfn;

    // Seeds: Outer vector is tables (r), inner vector is hash functions (b)
    std::vector<std::vector<seed_t>> table_seeds;

    // Hash Tables: Outer vector corresponds to the r tables
    std::vector<std::unordered_map<uint64_t, std::vector<uint32_t>>> hash_tables;

    // This is a robust hash combiner to replace bitwise AND, which can
    // heavily bias the signature distribution toward 0.
    static inline uint64_t combine_signatures(uint64_t seed, uint64_t hash_val) {
        // Based on golden ratio mix (similar to boost::hash_combine but 64-bit) https://stackoverflow.com/a/2595226
        return seed ^ (hash_val + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
    }

    [[nodiscard]] uint64_t compute_band_signature(const std::string& inp_string, uint32_t table_idx) const {
        // uint64_t band_signature = 0; // Initial seed for the combiner
        // starting with a non zero seed. 
        uint64_t band_signature = static_cast<uint64_t>(table_idx + 1) * 0x9e3779b97f4a7c15ULL; // Just to start with a non zero value.
        const auto* strPtr = reinterpret_cast<const uint8_t *>(inp_string.c_str());
        const auto strLen = inp_string.length();
        for (uint32_t k = 0; k < hashesPerTable_b; k++) {
            uint64_t current_hash = 0;//TODO: make the output datatype flexible based on the hash function's output type.
            hashfn(strPtr, strLen, table_seeds[table_idx][k], &current_hash);
            band_signature = combine_signatures(band_signature, current_hash);
        }
        return band_signature;
    }


public:
    ANN_LSH_Index(uint32_t r, uint32_t b, HashFn hash_fn, SeedGenerator& seedGen)
        : numTables_r(r), hashesPerTable_b(b), hashfn(hash_fn) {
        table_seeds.resize(r, std::vector<seed_t>(b));  // Allocate seeds for r tables and b hash functions per table
        hash_tables.resize(r);  // Allocate r hash tables (initially empty)

        // Initialize b independent hash seeds for each of the L tables
        for (uint32_t l = 0; l < r; l++) {
            for (uint32_t k = 0; k < b; k++) {
                table_seeds[l][k] = seedGen.nextSeed();
            }
        }
    }

    void insert(const std::string& inp_string, uint32_t position) {
        for (uint32_t l = 0; l < numTables_r; l++) {
            uint64_t band_signature = compute_band_signature(inp_string, l);
            hash_tables[l][band_signature].push_back(position);
        }
    }

    // [[nodiscard]] std::vector<uint32_t> query(const std::string& query_seq) const {
    //     std::set<uint32_t> unique_candidates;

    //     for (uint32_t l = 0; l < numTables_r; l++) {
    //         uint64_t band_signature = compute_band_signature(query_seq, l);
    //         auto it = hash_tables[l].find(band_signature);
    //         if (it != hash_tables[l].end()) {
    //             for (uint32_t pos : it->second) {
    //                 unique_candidates.insert(pos);
    //             }
    //         }
    //     }
    //     return {unique_candidates.begin(), unique_candidates.end()};
    // }

    [[nodiscard]] std::vector<uint32_t> query(const std::string& query_seq) const {
        std::vector<uint32_t> candidates;

        for (uint32_t l = 0; l < numTables_r; l++) {
            uint64_t band_signature = compute_band_signature(query_seq, l);
            auto it = hash_tables[l].find(band_signature);
            if (it != hash_tables[l].end()) {
                // Just dump all candidates in, no dedup here
                candidates.insert(candidates.end(), it->second.begin(), it->second.end());
            }
        }
        return candidates;  // may contain duplicates — caller handles it
    }

    [[nodiscard]] size_t index_size_bytes() const {
        size_t total = 0;
        for (uint32_t l = 0; l < numTables_r; l++) {
            // Each bucket: one uint64_t key + vector of uint32_t positions
            for (const auto& bucket : hash_tables[l]) {
                total += sizeof(uint64_t);                              // key  //TODO: make it flexible based on the hash function's output type.
                total += bucket.second.size() * sizeof(uint32_t);      // stored positions
                total += sizeof(std::vector<uint32_t>);                 // vector overhead
            }
            // unordered_map overhead (approximate)
            total += hash_tables[l].bucket_count() * sizeof(void*);
        }
        // Seeds
        total += numTables_r * hashesPerTable_b * sizeof(seed_t);
        return total;
    }
    
    [[maybe_unused]] void print_stats() const {
        printf("\n--- LSH Index Statistics ---\n");
        printf("Tables (L): %u, Hashes per Table (K): %u\n", numTables_r, hashesPerTable_b);
        for (uint32_t l = 0; l < numTables_r; l++) {
            size_t total_items = 0;
            for(const auto& bucket : hash_tables[l]) {
                total_items += bucket.second.size();
            }
            printf("Table %u: %zu distinct buckets, %zu total items\n", l, hash_tables[l].size(), total_items);
        }
        printf("----------------------------\n\n");
    }
    [[maybe_unused]] void print_tables(size_t limit_per_table = 1000) const {
        printf("\n--- LSH Index Contents ---\n");
        for (uint32_t l = 0; l < numTables_r; ++l) {
            printf("--- Table %u ---\n", l);
            if (hash_tables[l].empty()) {
                printf("  (empty)\n");
            } else {
                size_t count = 0;
                for (const auto& pair : hash_tables[l]) {
                    if (count >= limit_per_table) {
                        printf("... and %zu more buckets\n", hash_tables[l].size() - count);
                        break;
                    }
                    printf("  Hash: 0x%016llx -> Positions: [", static_cast<unsigned long long>(pair.first));
                    for (size_t i = 0; i < pair.second.size(); ++i) {
                        printf("%u", pair.second[i]);
                        if (i < pair.second.size() - 1) {
                            printf(", ");
                        }
                    }
                    printf("]\n");
                    count++;
                }
            }
            printf("\n");
        }
        printf("------------------------\n\n");
    }
};

static sim_bins_struct LSHInnerAgg(const HashInfo * hinfo,SequenceRecordsWithMetadataStruct &sequenceRecordsForAgg, SeedGenerator &seedGen, std::ofstream &out_file);

void writeSequencesToFile(const SequenceRecordsWithMetadataStruct& sequenceRecords, const std::string& outputFilename, int flag);

static void run_progress_bar(std::atomic<uint32_t>& counter, uint32_t total);

//----------------------------------------------------------------------------//
template <typename hashtype>
bool LSHApproxNearestNeighbourTest(const HashInfo *hinfo, flags_t flags) {
    printf("[[[ Approximate Nearest Neighbour Test ]]]\n\n");

    // BIOLSHASHER_ASSERT(hinfo->hashfn == nullptr , "Error: hash function pointer is null!");
    HashFn hash = hinfo->hashFn();
    if (!hash) {
        printf("Error: hash function pointer is null!\n");
        exit(EXIT_FAILURE);
    }

	BIOLSHASHER_ASSERT(hinfo->isLocalSensitive(),"Flag FLAG_HASH_LOCALITY_SENSITIVE not found. Please ensure that your LSH hash function is defined with FLAG_HASH_LOCALITY_SENSITIVE tag. \n\t\t It may be possible that the provided hash is not an LSH function.");
    BIOLSHASHER_ASSERT(hinfo->similarityfn != nullptr, "LSH function should be defined with a similarity metric.");
    BIOLSHASHER_ASSERT((std::strlen(hinfo->similarity_name) > 0), "The similarity metric should be named.");


    bool result = true; //TODO 

    // Create results directory if it doesn't exist
	std::string results_dir_str = "../results";
	std::filesystem::path results_dir = results_dir_str;
	if (!std::filesystem::exists(results_dir)) {
		std::filesystem::create_directories(results_dir);
	}

	std::string filename = results_dir_str + "/approxNearestNeighbourResults_" + hinfo->similarity_name + ".csv";

	std::ios_base::openmode mode = std::ios::trunc;  // Default: replace
	if (std::filesystem::exists(filename)) {
		mode = std::ios::app; // File exists: append
	}
	std::ofstream out_file(filename, mode);
	if (!out_file.is_open()) {
		throw std::runtime_error("Error: Could not open output file");
	}

    uint32_t sequenceLength = g_SequenceLength;

    SetIsTestActive(true);
    SeedGenerator seedGen(g_GoldenRatio ^ std::chrono::system_clock::now().time_since_epoch().count());

    // -----------------------------------------------------------------------------------------------------
    // 1. Generate reference sequences.
    SequenceRecordsWithMetadataStruct ReferenceSequenceRecords;
    ReferenceSequenceRecords.KeyCount = g_Nseq_in_Database;
    ReferenceSequenceRecords.OriginalSequenceLength = sequenceLength;
    ReferenceSequenceRecords.isBasesDrawnFromUniformDist = true;
    ReferenceSequenceRecords.DataGenSeed = seedGen.nextSeed();
    ReferenceSequenceRecords.DataMutateSeed = 0; // not applicable

    // side-effect: populates ReferenceSequenceRecords
    [[maybe_unused]] SequenceDataGenerator dataGenReference(&ReferenceSequenceRecords);

    // print the reference sequence to the output file  //TODO: Add a random number as suffix. and proper output file positing.
    // writeSequencesToFile(ReferenceSequenceRecords, "../results/reference_sequences.txt", 0);

    // -----------------------------------------------------------------------------------------------------
    // 2. Generate query sequences by randomly sampling from the reference sequences.

    const uint32_t numQueries = g_numQueriesForApproxNNTest;
    BIOLSHASHER_ASSERT(numQueries < g_Nseq_in_Database, "Number of queries exceeds number of available reference sequences"); // Ensure we don't sample more queries than available sequences

    std::mt19937 rngQuery(seedGen.nextSeed());
    std::uniform_int_distribution<uint32_t> distribution(0, g_Nseq_in_Database - 1);    // [a,b], both boundaries are inclusive, thus subtracted 1.

    SequenceRecordsWithMetadataStruct QuerySequenceRecords;
    QuerySequenceRecords.KeyCount = numQueries;
    QuerySequenceRecords.OriginalSequenceLength = sequenceLength; // Sampled Queries which have not yet mutated are same length as reference.
    QuerySequenceRecords.isBasesDrawnFromUniformDist = g_isBasesDrawnFromUniformDistribution;
    QuerySequenceRecords.DataGenSeed = 0; // Not applicable, as we are sampling from reference sequences rather than generating new random sequences.
    QuerySequenceRecords.DataMutateSeed = seedGen.nextSeed();

    // Manual population of query sequences by sampling from unique k-mers
    QuerySequenceRecords.Records.resize(numQueries);
    for (uint32_t i = 0; i < numQueries; i++) {
        uint32_t idx = distribution(rngQuery);
        QuerySequenceRecords.Records[i].SeqASCIIOrg = ReferenceSequenceRecords.Records[idx].SeqASCIIOrg;
        QuerySequenceRecords.Records[i].OriginalLength = sequenceLength;
    }

    // Print the query sequences to the output file
    // writeSequencesToFile(QuerySequenceRecords, "../results/query_sequences.txt", 0);


    // -----------------------------------------------------------------------------------------------------
    // 3. Perform aggegration to get the bins from simx% to simy% similarity (note: simy>simx).

    SequenceRecordsWithMetadataStruct AggSequenceRecord;
    AggSequenceRecord.KeyCount = g_NAggCasesApproxNNTest;
    AggSequenceRecord.OriginalSequenceLength = sequenceLength; // Sampled Queries which have not yet mutated are same length as reference.
    AggSequenceRecord.isBasesDrawnFromUniformDist = g_isBasesDrawnFromUniformDistribution;
    AggSequenceRecord.DataGenSeed = 0; // Not applicable, as we are sampling from reference sequences rather than generating new random sequences.
    AggSequenceRecord.DataMutateSeed = seedGen.nextSeed();

    // side-effect: populates AggSequenceRecord
    [[maybe_unused]] SequenceDataGenerator dataGenAgg(&AggSequenceRecord);
    sim_bins_struct sim_bins = LSHInnerAgg(hinfo,AggSequenceRecord, seedGen, out_file);

    // remove the agg sequence records to free up memory, as we no longer need them.
    // We have already extracted the binned error parameters and their statistics
    // that we need for the next steps.
    AggSequenceRecord.Records.clear();
    AggSequenceRecord.Records.shrink_to_fit();

    // print bin means and stddevs using
    // if(REPORT(VERBOSE, flags)) {
    //     for (size_t bin_idx = 0; bin_idx < sim_bins.bin_error_parameters_mean.size(); bin_idx++) {
    //       printf("Bin %zu: Count %d, Mean = %0.2f, Stddev = %0.2f\n", bin_idx, sim_bins.bin_fill_count[bin_idx], sim_bins.bin_error_parameters_mean[bin_idx], sim_bins.bin_error_parameters_stddev[bin_idx]);
    //     }
    // }

    // Reset mutation-related fields in each record
    for (uint32_t i = 0; i < QuerySequenceRecords.KeyCount; i++) {
        QuerySequenceRecords.Records[i].SeqASCIIMut.clear();
        QuerySequenceRecords.Records[i].similarity = 0.0;
        QuerySequenceRecords.Records[i].foundationalParameter = 0.0;
        QuerySequenceRecords.Records[i].snpRate = 0.0;
        QuerySequenceRecords.Records[i].delRate = 0.0;
        QuerySequenceRecords.Records[i].stayRate = 0.0;
        QuerySequenceRecords.Records[i].insmean = 0.0;
        QuerySequenceRecords.Records[i].insRate = 0.0;
    }

    // Generate new error parameters targeting similarity range [low, high]
    double target_sim_low = g_simThresholdForApproxNNTest; // Only two digits are allowed after decimal.
    constexpr double target_sim_high = 1.0; // Only two digits are allowed after decimal.

    // Convert similarity values to bin indices using the same formula as the binning step
    // bin_index = (sim - binstart) / binsize, clamped to [0, bincount-1]
    // This makes the intent explicit and stays consistent with how bins are defined.
    auto sim_to_bin_idx = [&](double sim) -> int32_t {
        auto idx = static_cast<int32_t>((sim - AggSequenceRecord.binstart) / AggSequenceRecord.binsize);
        return std::clamp(idx, 0, static_cast<int32_t>(AggSequenceRecord.bincount) - 1);
    };

    const int32_t target_bin_low  = sim_to_bin_idx(target_sim_low);   // = 95
    const int32_t target_bin_high = sim_to_bin_idx(target_sim_high);  // = 99, not 100

    Randbin rng_bin_sampler(seedGen.nextSeed());
    Randbin rng_bin_params_sampler(seedGen.nextSeed());

    uint32_t skipped_sequences = 0;
    for (uint32_t idx = 0; idx < QuerySequenceRecords.KeyCount; idx++) {
        uint32_t bin_fill_count = 0;
        int sampled_binid = -1;
        uint32_t attempts = 0;
        constexpr uint32_t max_attempts = 1000;

        while (bin_fill_count == 0 && attempts < max_attempts) {
            sampled_binid = rng_bin_sampler.rand_custom_range(target_bin_low, target_bin_high);
            // printf("Attempt %u: Sampled bin ID = %d\n", attempts + 1,sampled_binid);
            bin_fill_count = sim_bins.bin_fill_count[sampled_binid];  // get the number of entries in this bin. We get this from the agg computation step.
            attempts++;
        }

        if (bin_fill_count == 0) {
            printf("Warning: Could not find non-empty bin after %u attempts\n", max_attempts);
            // sequenceRecordsforTest.Records[idx].snpRate = 1.0; // Assign a default value TODO: Need to handle this.
            skipped_sequences++;
            continue; // Skip this sequence or handle error
        }

        uint32_t rand_param_idx = rng_bin_params_sampler.rand_range(bin_fill_count);
        double sampled_error_param = sim_bins.bin_error_parameters[sampled_binid][rand_param_idx];
        QuerySequenceRecords.Records[idx].foundationalParameter = sampled_error_param;
    }

    if (skipped_sequences > 0) {
        printf("Skipped %u sequences due to empty bins.\n", skipped_sequences);
    }

    if (g_mutation_model == MUTATION_MODEL_SIMPLE_SNP_ONLY) {
        // Side-effect: populates SeqASCIIMut and similarity fields in sequenceRecordsForTest
        [[maybe_unused]] SequenceDataMutatorSubstitutionOnly dataMutTest(&QuerySequenceRecords, hinfo->similarityfn);
        // printf("Completed mutation using simple SNP only model.\n");
    } else if (g_mutation_model == MUTATION_MODEL_GEOMETRIC_MUTATOR) {
        // Side-effect: populates SeqASCIIMut and similarity fields in sequenceRecordsForTest
        [[maybe_unused]] SequenceDataMutatorGeometric dataMutTest(&QuerySequenceRecords, hinfo->similarityfn);
        // printf("Completed mutation using geometric mutator model.\n");
    }

    printf("Re-sampled mutated %u sequences targeting similarity in [%.2f, %.2f].\n", QuerySequenceRecords.KeyCount, target_sim_low, target_sim_high);

    // writeSequencesToFile(QuerySequenceRecords, "../results/query_sequences.txt", 2);

    // print the query, the mutated sequence and the similarity values to terminal
    // printf("\n\n");
    // for (uint32_t i = 0; i < QuerySequenceRecord.KeyCount; i++) {
    //   printf("Foundational Parameter for record %u: %f\n", i, QuerySequenceRecord.Records[i].foundationalParameter);
    //   printf("Similarity for record %u: %f\n", i, QuerySequenceRecord.Records[i].similarity);
    //   printf("Query Sequence: %s\n", QuerySequenceRecord.Records[i].SeqASCIIOrg.c_str());
    //   printf("Mutated Sequence: %s\n", QuerySequenceRecord.Records[i].SeqASCIIMut.c_str());
    // }

    // Debug section: Print similarity statistics across all query records
    double sum_sim = 0.0;
    double min_sim = 1.0;
    double max_sim = 0.0;
    for (uint32_t i = 0; i < QuerySequenceRecords.KeyCount; i++) {
        double s = QuerySequenceRecords.Records[i].similarity;
        sum_sim += s;
        if (s < min_sim){
            min_sim = s;
        }
        if (s > max_sim){
            max_sim = s;
        }
    }
    double mean_sim = sum_sim / QuerySequenceRecords.KeyCount;

    double sq_diff_sum = 0.0;   //variance
    for (uint32_t i = 0; i < QuerySequenceRecords.KeyCount; i++) {
        double diff = QuerySequenceRecords.Records[i].similarity - mean_sim;
        sq_diff_sum += diff * diff;
    }
    double stddev_sim = sqrt(sq_diff_sum / QuerySequenceRecords.KeyCount);

    printf("Similarity stats (%u records): Mean = %f, Min = %f, Max = %f, Stddev = %f\n", QuerySequenceRecords.KeyCount, mean_sim, min_sim, max_sim, stddev_sim);


    // -----------------------------------------------------------------------------------------------------
    // 5: Compute Ground Truth for c-ANN: For each mutated query, compute similarity to all reference sequences
    // and store entires above the minimum c threshold floor. At evaluation time, for each c value,
    // the true neighbour set = all entries with sim >= c * target_sim_low.
    const double min_c = *std::min_element(g_cValuesApproxNNTest.begin(), g_cValuesApproxNNTest.end());
    const double sim_floor = min_c * target_sim_low; // Only storing entries above this.

    // Note: In our case we are only dealing with c=1. That all. so we are using a threshold based similarity search.


    // struct GroundTruthEntry {
    //     uint32_t position; // Position in the reference database
    //     double similarity;
    // };

    // std::vector<std::vector<GroundTruthEntry>> groundTruthAll(QuerySequenceRecords.KeyCount);

    std::vector<std::vector<uint32_t>> groundTruthAll(QuerySequenceRecords.KeyCount);

    printf("Computing ground truth (c-ANN) for %u queries across %u reference sequences (sim floor = %.4f)...\n", QuerySequenceRecords.KeyCount, g_Nseq_in_Database, sim_floor);

    std::atomic<uint32_t> gt_completed(0);

    std::vector<uint32_t> valid_Q = std::vector<uint32_t>(QuerySequenceRecords.KeyCount, 0);    //Important when performing search. Skip the ones with the value 0 after the GROUND TRUTH.

    // Returns [start, end) range for thread i over a total of n items
    auto make_range = [](uint32_t i, uint32_t n_threads, uint32_t n_items) -> std::pair<uint32_t, uint32_t> {
        uint32_t block = n_items / n_threads;
        uint32_t start = i * block;
        uint32_t end   = (i == n_threads - 1) ? n_items : start + block;
        return {start, end};
    };

    auto gt_worker = [&](uint32_t start, uint32_t end) {
        for (uint32_t q_idx = start; q_idx < end; q_idx++) {
            const std::string &querySeq = QuerySequenceRecords.Records[q_idx].SeqASCIIMut;
            const uint32_t mut_len = QuerySequenceRecords.Records[q_idx].MutatedLength;
            groundTruthAll[q_idx].reserve(std::max(96u, g_Nseq_in_Database / 1000)); // tune based on expected hit rate

            for (uint32_t ref_pos = 0; ref_pos < g_Nseq_in_Database; ++ref_pos) {
                const std::string &refSeq = ReferenceSequenceRecords.Records[ref_pos].SeqASCIIOrg;
                double sim = hinfo->similarityfn(querySeq, refSeq, mut_len, sequenceLength);
                if (sim >= sim_floor) {
                    groundTruthAll[q_idx].push_back(ref_pos); // position only
                    valid_Q[q_idx] |= 1;
                }
            }

            // No sort needed — eval phase uses unordered_set lookup, order irrelevant
            // std::sort(groundTruthAll[q_idx].begin(), groundTruthAll[q_idx].end(),
            //           [](const GroundTruthEntry &a, const GroundTruthEntry &b) {
            //               return a.similarity > b.similarity;
            //           });

            gt_completed.fetch_add(1, std::memory_order_relaxed);
        }
    };

    if (g_NCPU == 1) {
        std::thread gt_progress(run_progress_bar, std::ref(gt_completed), numQueries);
        gt_worker(0, numQueries);
        gt_progress.join();
    }
    else {
        #if defined(HAVE_THREADS)
            std::vector<std::thread> gt_threads(g_NCPU);
            
            std::thread gt_progress(run_progress_bar, std::ref(gt_completed), numQueries);
            
            for (uint32_t i = 0; i < g_NCPU; i++) {
                auto [s, e] = make_range(i, g_NCPU, numQueries);
                gt_threads[i] = std::thread(gt_worker, s, e);
            }

            for (auto &t : gt_threads) t.join();
            
            gt_progress.join();
        #endif
    }

    uint32_t total_valid = std::accumulate(valid_Q.begin(), valid_Q.end(), 0);
    printf("Ground truth computation completed for %u queries.\nValid queries with at least one neighbor above sim floor: %u\n", QuerySequenceRecords.KeyCount, total_valid);
    
    if(total_valid == 0){
        printf("Error: No valid queries with neighbors above similarity floor. Cannot proceed with ANN evaluation.\n");
        return false; // TODO:  handle as appropriate
    }
  
    // Write ground truth to output file
    // out_file << ":15:Ground Truth Top-" << TOP_K << " Nearest Sequences\n";
    // for(uint32_t q_idx = 0; q_idx < QuerySequenceRecord.KeyCount; q_idx++){
    //   out_file << "Q" << q_idx << ":";
    //   for(size_t k = 0; k < groundTruthNearest[q_idx].size(); k++){
    //     const auto& entry = groundTruthNearest[q_idx][k];
    //     out_file << entry.position << "," << entry.similarity;
    //     if(k < groundTruthNearest[q_idx].size() - 1) out_file << ";";
    //   }
    //   out_file << "\n";
    // }

    // -----------------------------------------------------------------------------------------------------
    // 6. LSH Index Construction and Querying

    // Define LSH parameters and number of runs
    const uint32_t NUM_RUNS = g_avgRunsForApproxNN;
    std::vector<std::pair<uint32_t, uint32_t>> br_pairs;
    // Programmatically fill br_pairs from b=1 to 10 and r=1 to 10
    BIOLSHASHER_ASSERT(g_ANN_start_B > 0 && g_ANN_start_R > 0 && g_ANN_MAX_B >= g_ANN_start_B && g_ANN_MAX_R >= g_ANN_start_R, "Invalid LSH parameters.");

    for (uint32_t b_val = g_ANN_start_B; b_val <= g_ANN_MAX_B; ++b_val) {
        for (uint32_t r_val = g_ANN_start_R; r_val <= g_ANN_MAX_R; ++r_val) {
          br_pairs.push_back({b_val, r_val});
      }
    }
    printf("Generated %zu (b,r) pairs for testing (b from %u to %u, r from %u to %u).\n", br_pairs.size(), g_ANN_start_B, g_ANN_MAX_B, g_ANN_start_R, g_ANN_MAX_R);


    // File header
    out_file << "\n:1:LSH Approx Nearest Neighbour Summary\n";
    out_file << ":2:" << "Hashname," << "SequenceLength," << "Distance Metric," << "Mutation Model," << "Mutation Expression" << std::endl;
    out_file << ":3:" << hinfo->name << "," << sequenceLength << "," << hinfo->similarity_name << "," << g_mutation_model << "," << g_mutation_expression_type << std::endl;

    // Print hash function ultra-specific parameters to the output file
    hinfo->printParameters(out_file);

    // Header — update this line
    out_file << ":5:b,r,c,Avg_Recall,Avg_Precision,Avg_FPR,Avg_Query_Time_ms,Avg_Index_Size_MB\n";

    const size_t num_c = g_cValuesApproxNNTest.size();

    for (const auto& br_pair : br_pairs) {
        uint32_t b = br_pair.first;  // Hashes per table
        uint32_t r = br_pair.second; // Number of tables

        printf("\n==================================================================\n");
        printf("Running experiment for b=%u, r=%u over %u runs\n", b, r, NUM_RUNS);
        printf("==================================================================\n");

        // Per-c accumulators: for each c value, accumulate metrics across runs
        std::vector<std::vector<double>> all_runs_avg_recall(num_c, std::vector<double>(NUM_RUNS));
        std::vector<std::vector<double>> all_runs_avg_precision(num_c, std::vector<double>(NUM_RUNS));
        std::vector<std::vector<double>> all_runs_avg_fpr(num_c, std::vector<double>(NUM_RUNS));

        std::vector<double> all_runs_query_time_ms(NUM_RUNS, 0.0);   // avg ms per query
        std::vector<double> all_runs_index_size_bytes(NUM_RUNS, 0.0); // total index size in bytes

        for (uint32_t run = 0; run < NUM_RUNS; ++run) {
            // printf("\n--- Run %u/%u for b=%u, r=%u ---\n", run + 1, NUM_RUNS, b, r);
            // printf("\n--- Run %u/%u ---\n", run + 1, NUM_RUNS);

            // Initialise the LSH index class object
            ANN_LSH_Index lsh_index(r, b, hash, seedGen);

            // Populate Index from the reference database
            auto build_start = std::chrono::high_resolution_clock::now();
            for (uint32_t i = 0; i < g_Nseq_in_Database; ++i) {
                lsh_index.insert(ReferenceSequenceRecords.Records[i].SeqASCIIOrg, i);
            }
            auto build_end = std::chrono::high_resolution_clock::now();
            double build_time_ms = std::chrono::duration<double, std::milli>(build_end - build_start).count();
            printf("  Index build time: %.2f ms\n", build_time_ms);

            double index_bytes = static_cast<double>(lsh_index.index_size_bytes());
            all_runs_index_size_bytes[run] = index_bytes;
            printf("  Index size: %.2f MB\n", index_bytes / (1024.0 * 1024.0));
            
            // ============================================================================
            // Query and Evaluation Phase for this run
            // ============================================================================
            // Per-c, per-query metric vectors

            // [c_idx][q_idx]
            std::vector<std::vector<double>> recall_arr   (num_c, std::vector<double>(numQueries, 0.0));
            std::vector<std::vector<double>> precision_arr(num_c, std::vector<double>(numQueries, 0.0));
            std::vector<std::vector<double>> fpr_arr      (num_c, std::vector<double>(numQueries, 0.0));
            std::vector<double> query_time_arr(numQueries, 0.0); // ms per query

            std::atomic<uint32_t> eval_completed(0);

            // lsh_index.query() is const and stateless — safe to call concurrently
            // auto eval_worker = [&](uint32_t start, uint32_t end) {
            //     for (uint32_t q_idx = start; q_idx < end; q_idx++) {
            //         uint32_t is_valid = valid_Q[q_idx];
            //         if (is_valid) {
            //             const std::string &querySeq  = QuerySequenceRecords.Records[q_idx].SeqASCIIMut;
            //             //TODO: MAKE IT EFFICIENT.
            //             std::vector<uint32_t> lsh_results = lsh_index.query(querySeq);
            //             std::set<uint32_t> lsh_results_set(lsh_results.begin(), lsh_results.end());

            //             const auto &ground_truth_full = groundTruthAll[q_idx];

            //             // Evaluate for each c value
            //             for (size_t c_idx = 0; c_idx < num_c; c_idx++) {    // we have kept it in even though no c parameter is effectively used.
            //                 double c              = g_cValuesApproxNNTest[c_idx];
            //                 double sim_threshold  = c * target_sim_low;

            //                 // Build ground truth set: all entries with sim >= c * target_sim_low
            //                 std::set<uint32_t> ground_truth_positions;
            //                 for (const auto &entry : ground_truth_full) {
            //                     if (entry.similarity >= sim_threshold) {
            //                         ground_truth_positions.insert(entry.position);
            //                     }
            //                 }

            //                 uint32_t true_positives  = 0;
            //                 uint32_t false_positives = 0;

            //                 for (uint32_t pos : lsh_results_set) {
            //                     if (ground_truth_positions.count(pos)) {
            //                         true_positives++;
            //                     } else {
            //                         false_positives++;
            //                     }
            //                 }

            //                 int false_negatives = static_cast<int>(ground_truth_positions.size()) - true_positives;
            //                 BIOLSHASHER_ASSERT(false_negatives >= 0, "Internal issue: False negatives cannot be negative");

            //                 double recall    = (true_positives + false_negatives > 0) ? static_cast<double>(true_positives) / (true_positives + false_negatives) : 0.0;
            //                 double precision = (true_positives + false_positives > 0) ? static_cast<double>(true_positives) / (true_positives + false_positives) : 0.0;

            //                 uint32_t total_negatives = g_Nseq_in_Database - static_cast<uint32_t>(ground_truth_positions.size());
            //                 double fpr = (total_negatives > 0) ? static_cast<double>(false_positives) / total_negatives : 0.0;

            //                 // Disjoint write — thread owns this q_idx column
            //                 recall_arr[c_idx][q_idx]    = recall;
            //                 precision_arr[c_idx][q_idx] = precision;
            //                 fpr_arr[c_idx][q_idx]       = fpr;
            //             }
            //         }
            //         eval_completed.fetch_add(1, std::memory_order_relaxed);
            //     }
            // };

            // This new method will NOT WORK with multiple c values. only c=1.0.
            auto eval_worker = [&](uint32_t start, uint32_t end) {
                for (uint32_t q_idx = start; q_idx < end; q_idx++) {
                    if (valid_Q[q_idx]) {
                        const std::string &querySeq = QuerySequenceRecords.Records[q_idx].SeqASCIIMut;

                        // Time just the query — not the eval logic
                        auto q_start = std::chrono::high_resolution_clock::now();
                        std::vector<uint32_t> lsh_results = lsh_index.query(querySeq);
                        auto q_end = std::chrono::high_resolution_clock::now();
                        query_time_arr[q_idx] = std::chrono::duration<double, std::milli>(q_end - q_start).count();

                        std::unordered_set<uint32_t> lsh_results_set(lsh_results.begin(), lsh_results.end());

                        // Built once per query — no per-c rebuild needed since c=1.0
                        // means sim_threshold == sim_floor, which is exactly what GT was filtered by
                        const std::unordered_set<uint32_t> query_ground_truth_set(
                            groundTruthAll[q_idx].begin(),
                            groundTruthAll[q_idx].end()
                        );
                        const uint32_t gt_size = static_cast<uint32_t>(query_ground_truth_set.size());

                        for (size_t c_idx = 0; c_idx < num_c; c_idx++) {
                            // c loop kept for structural consistency, but ground truth
                            // set is identical for all c when c=1.0 and sim_floor = target_sim_low

                            uint32_t true_positives  = 0;
                            uint32_t false_positives = 0;

                            for (uint32_t pos : lsh_results_set) {
                                if (query_ground_truth_set.count(pos)) {
                                    true_positives++;
                                } else {
                                    false_positives++;
                                }
                            }

                            int false_negatives = static_cast<int>(gt_size) - true_positives;
                            BIOLSHASHER_ASSERT(false_negatives >= 0, "Internal issue: False negatives cannot be negative");

                            double recall    = (gt_size > 0)
                                            ? static_cast<double>(true_positives) / gt_size
                                            : 0.0;
                            double precision = (true_positives + false_positives > 0)
                                            ? static_cast<double>(true_positives) / (true_positives + false_positives)
                                            : 0.0;

                            uint32_t total_negatives = g_Nseq_in_Database - gt_size;
                            double fpr = (total_negatives > 0)
                                    ? static_cast<double>(false_positives) / total_negatives
                                    : 0.0;

                            recall_arr[c_idx][q_idx]    = recall;
                            precision_arr[c_idx][q_idx] = precision;
                            fpr_arr[c_idx][q_idx]       = fpr;
                        }
                    }
                    eval_completed.fetch_add(1, std::memory_order_relaxed);
                }
            };

            if (g_NCPU == 1) {
                std::thread eval_progress(run_progress_bar, std::ref(eval_completed), numQueries);
                eval_worker(0, numQueries);
                eval_progress.join();
            }
            else {
                #if defined(HAVE_THREADS)
                    std::vector<std::thread> eval_threads(g_NCPU);
                    std::thread eval_progress(run_progress_bar, std::ref(eval_completed), numQueries);
                    for (uint32_t i = 0; i < g_NCPU; i++) {
                        auto [s, e] = make_range(i, g_NCPU, numQueries);
                        eval_threads[i] = std::thread(eval_worker, s, e);
                    }
                    for (auto &t : eval_threads) t.join();
                    eval_progress.join();
                #endif
            }

            //-------------------------------------------------------------------------------------
            // Calculate averaged metrics for this single run, per c
            for (size_t i = 0; i < g_cValuesApproxNNTest.size(); ++i) {
                double c = g_cValuesApproxNNTest[i];

                double run_avg_recall = 0.0;
                double run_avg_precision = 0.0;
                double run_avg_fpr = 0.0;
                double run_avg_query_time_ms = 0.0;

                printf("Metrics for Run %u/%u:", run + 1, NUM_RUNS);
                
                for (double r_val : recall_arr[i])
                    run_avg_recall += r_val;
                run_avg_recall /= total_valid; //QuerySequenceRecords.KeyCount;

                for (double p_val : precision_arr[i])
                    run_avg_precision += p_val;
                run_avg_precision /= total_valid; //QuerySequenceRecords.KeyCount;

                for (double fpr_val : fpr_arr[i])
                    run_avg_fpr += fpr_val;
                run_avg_fpr /= total_valid; //QuerySequenceRecords.KeyCount;

                for (uint32_t q_idx = 0; q_idx < numQueries; q_idx++) {
                    run_avg_query_time_ms += query_time_arr[q_idx];
                }
                run_avg_query_time_ms /= total_valid;

                printf(" Recall=%.4f, Precision=%.4f, FPR=%.4f \n", run_avg_recall, run_avg_precision, run_avg_fpr);

                // Store metrics for this run
                all_runs_avg_recall[i][run] = run_avg_recall;
                all_runs_avg_precision[i][run] = run_avg_precision;
                all_runs_avg_fpr[i][run] = run_avg_fpr;
                all_runs_query_time_ms[run] = run_avg_query_time_ms;
            }
        }

        // Average the metrics across all runs for the current (b, r) pair, per c
        for (size_t c_idx = 0; c_idx < g_cValuesApproxNNTest.size(); ++c_idx) {
            double c = g_cValuesApproxNNTest[c_idx];
            double final_avg_recall = 0.0;
            double final_avg_precision = 0.0;
            double final_avg_fpr = 0.0;
            double final_avg_query_time_ms = 0.0;
            double final_avg_index_size_bytes = 0.0;
            

            for (double val : all_runs_avg_recall[c_idx])
                final_avg_recall += val;
            final_avg_recall /= NUM_RUNS;

            for (double val : all_runs_avg_precision[c_idx])
                final_avg_precision += val;
            final_avg_precision /= NUM_RUNS;

            for (double val : all_runs_avg_fpr[c_idx])
                final_avg_fpr += val;
            final_avg_fpr /= NUM_RUNS;

            for (double val : all_runs_query_time_ms)
                final_avg_query_time_ms += val;
            final_avg_query_time_ms /= NUM_RUNS;

            for (double val : all_runs_index_size_bytes)
                final_avg_index_size_bytes += val;
            final_avg_index_size_bytes /= NUM_RUNS;


            printf("\n--- Average Metrics for b=%u, r=%u (over %u runs) ---\n", b, r, NUM_RUNS);
            printf("Average Recall:    %.4f\n", final_avg_recall);
            printf("Average Precision: %.4f\n", final_avg_precision);
            printf("Average FPR:       %.4f\n", final_avg_fpr);
            printf("Avg Query Time:   %.4f ms\n", final_avg_query_time_ms);
            printf("Avg Index Size:   %.2f MB\n", final_avg_index_size_bytes / (1024.0 * 1024.0));
            printf("-----------------------------------------------------\n");

          // Write summary for this (b, r, c) triple to file
          out_file << ":6:" << b << "," << r << "," << c << ","
            << final_avg_recall << "," << final_avg_precision << ","
            << final_avg_fpr << ","
            << final_avg_query_time_ms << ","
            << (final_avg_index_size_bytes / (1024.0 * 1024.0)) << "\n";
        }
    }

    // Cleanup: Reset LSH global variables after test completion
    SetIsTestActive(false);

    // printf("%s\n", result ? "" : g_failstr);
    out_file.close();
    return result;
}


INSTANTIATE(LSHApproxNearestNeighbourTest, HASHTYPELIST);

static sim_bins_struct LSHInnerAgg(const HashInfo * hinfo,SequenceRecordsWithMetadataStruct &sequenceRecordsForAgg, SeedGenerator &seedGen, std::ofstream &out_file) {
    // printf("Inside get_params_aggregated_in_bins\n");

    sim_bins_struct sim_bins_agg; // Bin to store the error parameters

    std::vector<double> rand_error_param(sequenceRecordsForAgg.KeyCount, 0.0);
    std::vector<double> similarity_values(sequenceRecordsForAgg.KeyCount, 0.0);

    // printf("\n-----------------Start AGG----------------\n");
        if (g_mutation_model == MUTATION_MODEL_SIMPLE_SNP_ONLY) {
            // Side-effect: populates rand_error_param and similarity fields in sequenceRecordsForAgg
            // printf("Using simple SNP only model for data mutation in aggregation phase.\n");
            [[maybe_unused]] SequenceDataMutatorSubstitutionOnly dataMutAgg(&sequenceRecordsForAgg, &rand_error_param, hinfo->similarityfn);
    } else if (g_mutation_model == MUTATION_MODEL_GEOMETRIC_MUTATOR) {
            // Side-effect: populates rand_error_param and similarity fields in sequenceRecordsForAgg
            printf("Using geometric mutator model for data mutation in aggregation phase.\n");
            [[maybe_unused]] SequenceDataMutatorGeometric dataMutAgg(&sequenceRecordsForAgg, &rand_error_param, hinfo->similarityfn);
    }
    // printf("\n-----------------End AGG----------------");

    // Print Similarity values
    // out_file << ":13:";
    // for (size_t i = 0; i < sequenceRecordsForAgg.KeyCount; i++) {
    //   if (i == sequenceRecordsForAgg.KeyCount - 1)
    //     out_file << rand_error_param[i] << "\n";
    //   else
    //     out_file << rand_error_param[i] << ",";
    // }

    // Extract similarity values
    for (uint32_t idx = 0; idx < sequenceRecordsForAgg.KeyCount; idx++) {
        similarity_values[idx] = sequenceRecordsForAgg.Records[idx].similarity;
    }

    // Perform binning on similarity values
    for (size_t i = 0; i < sequenceRecordsForAgg.KeyCount; i++) {
        const double sim_value = similarity_values[i];
        auto bin_index = static_cast<uint32_t>((sim_value - sequenceRecordsForAgg.binstart) / sequenceRecordsForAgg.binsize);
        if (bin_index >= sequenceRecordsForAgg.bincount) {
          bin_index = sequenceRecordsForAgg.bincount - 1; // Clamp to last bin
        }

        const uint32_t bin_fill_count = sim_bins_agg.bin_fill_count[bin_index];

        if (bin_fill_count == g_bincount_full) {
          continue; // To avoid excessive memory usage, we limit to x entries per bin.
        }

        sim_bins_agg.bin_error_parameters[bin_index][bin_fill_count] = rand_error_param[i];
        sim_bins_agg.bin_fill_count[bin_index]++;
    }

    // Compute mean and stddev for each bin
    for (size_t bin_idx = 0; bin_idx < sequenceRecordsForAgg.bincount; bin_idx++) {
        const auto &params = sim_bins_agg.bin_error_parameters[bin_idx];
        if (!params.empty()) {
          // Compute mean
          double sum = 0.0;
          for (uint32_t param_idx = 0; param_idx < sim_bins_agg.bin_fill_count[bin_idx]; param_idx++) {
            sum += params[param_idx];
          }
        //   printf("Bin %zu: Sum = %f, Count = %u\n", bin_idx, sum, sim_bins_agg.bin_fill_count[bin_idx]);
          double mean = sum / sim_bins_agg.bin_fill_count[bin_idx];
          sim_bins_agg.bin_error_parameters_mean[bin_idx] = mean;

          // Compute stddev
          double sq_sum = 0.0;
          for (uint32_t param_idx = 0; param_idx < sim_bins_agg.bin_fill_count[bin_idx]; param_idx++) {
            sq_sum += (params[param_idx] - mean) * (params[param_idx] - mean);
          }
          double stddev = std::sqrt(sq_sum / sim_bins_agg.bin_fill_count[bin_idx]);
          sim_bins_agg.bin_error_parameters_stddev[bin_idx] = stddev;
        }
        else {
          sim_bins_agg.bin_error_parameters_mean[bin_idx] = 0.0;
          sim_bins_agg.bin_error_parameters_stddev[bin_idx] = 0.0;
        }
    }

    // // Print each bin's fill count
    // for (size_t bin_idx = 0; bin_idx < sequenceRecordsForAgg.bincount;
    // bin_idx++) { 	printf("Bin %zu: Fill Count = %u", bin_idx,
    // sim_bins_agg.bin_fill_count[bin_idx]);
    // 	// Print mean and stddev
    // 	printf(", Mean = %0.2f, Stddev = %0.2f\n",
    // sim_bins_agg.bin_error_parameters_mean[bin_idx],
    // sim_bins_agg.bin_error_parameters_stddev[bin_idx]);
    // 	// // Print error parameters in this bin
    // 	// printf(" | Error Parameters: ");
    // 	// for (size_t j = 0; j < sim_bins_agg.bin_fill_count[bin_idx]; j++) {
    // 	// 	printf("%f ", sim_bins_agg.bin_error_parameters[bin_idx][j]);
    // 	// }
    // 	// printf("\n");
    // }

    //--------------------------------------------	//
    // Fill the partially filled bins using the 	//
    // mean and stddev values computed above. 		//
    //--------------------------------------------	//
    // Now, I want to go through each bin, and if its fill count is less than
    // g_bincount_full, I want to fill the remaining entries with random samples
    // from a normal distribution defined by the mean and stddev of that bin.

    std::mt19937 gen(seedGen.nextSeed()); // Seeded RNG for reproducibility

    for (size_t bin_idx = 0; bin_idx < sequenceRecordsForAgg.bincount; bin_idx++) {
        uint32_t current_fill = sim_bins_agg.bin_fill_count[bin_idx];

        if (current_fill < g_bincount_full && current_fill > 0) {
            // Get the mean and stddev for this bin
            double mean = sim_bins_agg.bin_error_parameters_mean[bin_idx];
            double stddev = sim_bins_agg.bin_error_parameters_stddev[bin_idx];

            // If stddev is 0 or very small, use a small default to avoid degenerate
            // distribution
            if (stddev < 1e-9) {
                stddev = 0.01;
            }

            std::normal_distribution<double> normal_dist(mean, stddev);

            // Fill remaining entries
            for (uint32_t fill_idx = current_fill; fill_idx < g_bincount_full; fill_idx++) {
                double sampled_value = normal_dist(gen);

                // Clamp to valid range [0.0, 1.0] since these are error rates
                sampled_value = std::max(0.0, std::min(1.0, sampled_value));

                sim_bins_agg.bin_error_parameters[bin_idx][fill_idx] = sampled_value;
                sim_bins_agg.bin_fill_count[bin_idx]++;
            }
        }
        else if (current_fill == 0) {
            // Todo: Add the part where in angular similarity, we geniuinely have
            // empty bins from 0 to 49

            // // Bin is completely empty - use bin center as mean with small stddev
            // double bin_center = sequenceRecordsForAgg.binstart + (bin_idx + 0.5) *
            // sequenceRecordsForAgg.binsize; double default_stddev =
            // sequenceRecordsForAgg.binsize / 4.0;  // Quarter of bin width

            // std::normal_distribution<double> normal_dist(bin_center,
            // default_stddev);

            // for (uint32_t fill_idx = 0; fill_idx < g_bincount_full; fill_idx++) {
            // 	double sampled_value = normal_dist(gen);
            // 	sampled_value = std::max(0.0, std::min(1.0, sampled_value));

            // 	sim_bins_agg.bin_error_parameters[bin_idx][fill_idx] =
            // sampled_value; 	sim_bins_agg.bin_fill_count[bin_idx]++;
            // }

            // // Update mean and stddev for this bin
            // sim_bins_agg.bin_error_parameters_mean[bin_idx] = bin_center;
            // sim_bins_agg.bin_error_parameters_stddev[bin_idx] = default_stddev;

            sim_bins_agg.bin_fill_count[bin_idx] = 0;
            sim_bins_agg.bin_error_parameters_mean[bin_idx] = 0.0;
            sim_bins_agg.bin_error_parameters_stddev[bin_idx] = 0.0;
        }
    }
    // Clean up
    // sequenceRecordsForAgg.Records.clear();
    // sequenceRecordsForAgg.Records.shrink_to_fit();  // Actually releases the
    // memory

    return sim_bins_agg;
}

void writeSequencesToFile(const SequenceRecordsWithMetadataStruct& sequenceRecords, const std::string& outputFilename, int flag) {
    // Open the file in write mode
    std::ofstream outFile(outputFilename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << outputFilename << " for writing.\n";
        return;
    }

    if(flag == 0){
      // Iterate over the records and write each original sequence to the file
      for (const auto& record : sequenceRecords.Records) {
          outFile << record.SeqASCIIOrg << "\n";
      }
    }
    else if(flag == 1){
      // Iterate over the records and write each mutated sequence to the file
      for (const auto& record : sequenceRecords.Records) {
          outFile << record.SeqASCIIMut << "\n";
      }
    }
    else if(flag == 2){
      // Iterate over the records and write each original and mutated sequence to the file
      for (const auto& record : sequenceRecords.Records) {
          outFile << "Original: " << record.SeqASCIIOrg << "\n";
          outFile << "Mutated: " << record.SeqASCIIMut << "\n";
          outFile << "Similarity: " << record.similarity << "\n";
          outFile << "-------------------------\n";
      }
    }
    

    // Close the file
    outFile.close();
    printf("Sequences written to %s\n", outputFilename.c_str());
}

static void run_progress_bar(std::atomic<uint32_t>& counter, uint32_t total) {
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> qdist(0, SCIENTIST_QUOTES.size() - 1);
    std::string current_msg = ""; //"\"" + SCIENTIST_QUOTES[qdist(rng)].quote + "\" — " + SCIENTIST_QUOTES[qdist(rng)].author;
    // int ticks = 0;

    const auto& q = SCIENTIST_QUOTES[qdist(rng)];
    current_msg = "\"" + q.quote + "\" — " + q.author;

    while (true) {
        int done = counter.load();
        float pct = static_cast<float>(done) / total * 100.0f;
        int filled = static_cast<int>(pct / 2);

        // if (ticks % 1000 == 0) {
        // 	const auto& q = SCIENTIST_QUOTES[qdist(rng)];
        // 	current_msg = "\"" + q.quote + "\" — " + q.author;
        // }

        // Clear line, print bar, then quote on next line
        printf("\033[2K\r[%-50s] %.1f%% (%d/%d)\n\033[2K%s\033[1A",
            std::string(filled, '#').c_str(), pct, done, total,
            current_msg.c_str());
        fflush(stdout);
        
        if (done >= total){ 
            // Overwrite the two lines (bar + quote) with the clean final state
            printf("\033[2K\r[%-50s] 100.0%% (%d/%d)\n\n\033[1A\033[2K",
                std::string(50, '#').c_str(), total, total);
            fflush(stdout);
            break;
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        // ticks++;
    }
}
