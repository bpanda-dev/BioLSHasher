
#include "HashInfo.h"
#include "TestGlobals.h"
#include "Instantiate.h"
#include "coolstuff.h"

#include "BioDataGeneration.h"
#include "LSHGlobals.h"
#include "LSHCollision.h"
#include "fstream"
#include <iostream>
#include <filesystem>
#include "assertMsg.h"

#if defined(HAVE_THREADS)
  #include <atomic>
#endif

/*-------------------------------------------------------------------------------*/
/*									Collision Test		 						 */
/*-------------------------------------------------------------------------------*/


static void writeSequencesToFile(const SequenceRecordsWithMetadataStruct& sequenceRecords,const std::vector<double>& averageCollisions,const std::string& outputFilename, int flag) {
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
		outFile << "-------------------------\n";
		int i = 0;
		for (const auto& record : sequenceRecords.Records) {
			if( record.similarity >= 0.95 ){
				outFile << "O: " << record.SeqASCIIOrg << "\n";
				outFile << "M: " << record.SeqASCIIMut << "\n";
				outFile << "S: " << record.similarity << "\n";
				outFile << "P: " << averageCollisions[i] << "\n";
				outFile << "CountD:" << record.count_del << "\n";
				outFile << "CountS:" << record.count_sub << "\n";
				outFile << "CountI:" << record.count_ins << "\n";
				outFile << "CountIL:" << record.count_ins_length << "\n";
				outFile << "-------------------------\n";
			}
			i++;
		}
    }
    

    // Close the file
    outFile.close();
    printf("Sequences written to %s\n", outputFilename.c_str());
}



struct common_params_struct{
	uint32_t seqLen;
	seed_t DataGenSeed;
	seed_t DataMutateSeed;
	bool areBasesDrawnFromUniformDist;
};

static void run_progress_bar(std::atomic<uint32_t>& counter, uint32_t total);

template <typename hashtype>
static void LSHCollisionTestInnerInnerWorker(const HashInfo * hinfo, uint32_t N_hash, HashFn hash, seed_t HashSeed, common_params_struct &common_params, SequenceRecordsWithMetadataStruct &sequenceRecordsforTest, double *AverageCollision,int start, int end, uint32_t and_param, uint32_t or_param,  std::atomic<uint32_t>* completed){

	std::vector<hashtype> hash_val_org = std::vector<hashtype>(and_param*or_param);
	std::vector<hashtype> hash_val_mut = std::vector<hashtype>(and_param*or_param);

	std::vector<seed_t> andor_seeds = std::vector<seed_t>(and_param*or_param, 0);

	// This is for reducing atomic contention.
	uint32_t local_count = 0;  // Local counter to reduce atomic operations
    constexpr uint32_t batch_size = 10;  // Update shared counter every 10 sequences

	for(int rec_idx = start; rec_idx < end; rec_idx++){
		SequenceRecordUnit &record = sequenceRecordsforTest.Records[rec_idx];
		uint32_t collision_count = 0;
		
		// --- OPTIMIZATION START ---
        // We prepare the data ONCE per sequence, not N_hash times.
        // This prevents millions of memory allocations and lock contentions.
		const uint8_t* ptrOrg = reinterpret_cast<const uint8_t *>(record.SeqASCIIOrg.c_str());
		const size_t lenOrg = record.OriginalLength;
		const uint8_t* ptrMut = reinterpret_cast<const uint8_t *>(record.SeqASCIIMut.c_str());
		const size_t lenMut = record.MutatedLength;

		for(uint32_t hash_idx = 0; hash_idx < N_hash; hash_idx++){
			std::vector<uint32_t> and_flags = std::vector<uint32_t>(or_param, 0);
			for(uint32_t or_i = 0; or_i < or_param; or_i++){
				and_flags[or_i] = 1;
				for(uint32_t and_i = 0; and_i < and_param; and_i++){
					andor_seeds[or_i * and_param + and_i] = HashSeed + hash_idx + (or_i * and_param + and_i)*1337;

					hash(ptrOrg, lenOrg, andor_seeds[or_i * and_param + and_i], &hash_val_org[or_i * and_param + and_i]);
					hash(ptrMut, lenMut, andor_seeds[or_i * and_param + and_i], &hash_val_mut[or_i * and_param + and_i]);
				}
				// ANDing
				for(uint32_t i = 0; i < and_param; i++){
					if(hinfo->check_equality_fn(&hash_val_org[or_i * and_param + i], &hash_val_mut[or_i * and_param + i]) == false){
						and_flags[or_i] = 0;
						break;
					}
					// if(hash_val_org[or_i * and_param + i] != hash_val_mut[or_i * and_param + i]){
					// 	and_flags[or_i] = 0;
					// 	break;
					// }
				}
			}

			//AND-ORING
			uint32_t flag = 0;
			for(uint32_t i = 0; i < or_param; i++){
				if(and_flags[i] == 1){
					flag = 1;
					break;
				}
			}
			collision_count += flag;
		}

		// Compute average and stddev of collisions for this sequence pair.
		double avg_collision = static_cast<double>(collision_count) / static_cast<double>(N_hash);
		AverageCollision[rec_idx] = avg_collision;
		
		// Batch counter updates to reduce atomic contention
        local_count++;
        if (local_count >= batch_size) {
            completed->fetch_add(local_count, std::memory_order_relaxed);
            local_count = 0;
        }
		// // single increment per sequence, not too much atomic contention
        // completed->fetch_add(1, std::memory_order_relaxed);
	}
	// Flush remaining count at the end
    if (local_count > 0) {
        completed->fetch_add(local_count, std::memory_order_relaxed);
    }
}


template <typename hashtype>
static bool LSHCollisionTestInnerInnerParallel(const HashInfo * hinfo, uint32_t N_seq, uint32_t N_hash, HashFn hash, seed_t HashSeed, common_params_struct &common_params, sim_bins_struct &sim_bins, SeedGenerator &seedGen, std::ofstream &out_file){
	
	// printf("Inside LSHCollisionTestInnerInner\n");

	SequenceRecordsWithMetadataStruct sequenceRecordsforTest;
	sequenceRecordsforTest.OriginalSequenceLength = common_params.seqLen;
	sequenceRecordsforTest.areBasesDrawnFromUniformDist = common_params.areBasesDrawnFromUniformDist;
	sequenceRecordsforTest.DataGenSeed = seedGen.nextSeed();	// A seed that is different from the seed that is used in AGG step.
	sequenceRecordsforTest.DataMutateSeed = seedGen.nextSeed();
	sequenceRecordsforTest.KeyCount = N_seq;
	sequenceRecordsforTest.HashCount = N_hash;
	
	[[maybe_unused]] SequenceDataGenerator dataGenTest(&sequenceRecordsforTest);	// Generates the test data

	// For each of the sequence pair, draw a random mutation parameter from the bin statistics. Then store it in the mutation record. and mutate it.
	seed_t bin_sampling_seed = seedGen.nextSeed();
	seed_t bin_params_sampling_seed = seedGen.nextSeed();
	
	Randbin rng_bin_sampler(bin_sampling_seed);
	Randbin rng_bin_params_sampler(bin_params_sampling_seed);

	uint32_t n_skipped_sequences = 0;
	for(uint32_t idx = 0; idx < N_seq; idx++){
		uint32_t bin_fill_count = 0;
		int sampled_binid = -1;
		uint32_t attempts = 0;
		constexpr uint32_t max_attempts = 1000;

		while(bin_fill_count == 0 && attempts < max_attempts){
			sampled_binid = rng_bin_sampler.rand_range(100);
			bin_fill_count = sim_bins.bin_fill_count[sampled_binid];
			attempts++;
		}
		
		if(bin_fill_count == 0){
			printf("Warning: Could not find non-empty bin after %u attempts\n", max_attempts);
			// sequenceRecordsforTest.Records[idx].snpRate = 1.0; // Assign a default value TODO: Figure out what will work here.
			n_skipped_sequences++;
			continue;  // Skip this sequence or handle error
		}
		
		uint32_t rand_param_idx = rng_bin_params_sampler.rand_range(bin_fill_count);
		double sampled_error_param = sim_bins.bin_error_parameters[sampled_binid][rand_param_idx];
		sequenceRecordsforTest.Records[idx].foundationalParameter = sampled_error_param;
	}

	if(n_skipped_sequences > 0){
		printf("NOTE:Skipped %u sequences due to empty bins.\n", n_skipped_sequences);
	}

	if(g_mutation_model == MUTATION_MODEL_SIMPLE_SNP_ONLY){
		[[maybe_unused]] SequenceDataMutatorSubstitutionOnly dataMutTest(&sequenceRecordsforTest, hinfo->similarityfn);
		printf("Completed mutation using simple SNP only model.\n");
	}
	else if(g_mutation_model == MUTATION_MODEL_GEOMETRIC_MUTATOR){
		[[maybe_unused]] SequenceDataMutatorGeometric dataMutTest(&sequenceRecordsforTest, hinfo->similarityfn);
		printf("Completed mutation using geometric mutator model.\n");
	}

	// Print Similarity values
	out_file << ":5:";
	for (size_t i = 0; i < N_seq; i++) {
		if (i == N_seq - 1)
			out_file << sequenceRecordsforTest.Records[i].similarity << "\n";
		else
			out_file << sequenceRecordsforTest.Records[i].similarity << ",";
	}
	// Print param values
	out_file << ":6:";
	for (size_t i = 0; i < N_seq; i++) {
		if (i == N_seq - 1)
			out_file << sequenceRecordsforTest.Records[i].snpRate << "\n";
		else
			out_file << sequenceRecordsforTest.Records[i].snpRate << ",";
	}
	
	if(g_mutation_model == MUTATION_MODEL_GEOMETRIC_MUTATOR){	
		out_file << ":7:";
		for (size_t i = 0; i < N_seq; i++) {
			if (i == N_seq - 1)
				out_file << sequenceRecordsforTest.Records[i].delRate << "\n";
			else
				out_file << sequenceRecordsforTest.Records[i].delRate << ",";
		}
		out_file << ":8:";
		for (size_t i = 0; i < N_seq; i++) {
			if (i == N_seq - 1)
				out_file << sequenceRecordsforTest.Records[i].insmean << "\n";
			else
				out_file << sequenceRecordsforTest.Records[i].insmean << ",";
		}
		out_file << ":9:";
		for (size_t i = 0; i < N_seq; i++) {
			if (i == N_seq - 1)
				out_file << sequenceRecordsforTest.Records[i].stayRate << "\n";
			else
				out_file << sequenceRecordsforTest.Records[i].stayRate << ",";
		}
		out_file << ":14:";
		for (size_t i = 0; i < N_seq; i++) {
			if (i == N_seq - 1)
				out_file << sequenceRecordsforTest.Records[i].insRate << "\n";
			else
				out_file << sequenceRecordsforTest.Records[i].insRate << ",";
		}
	}

	alignas(64) std::vector<double> AverageCollision(N_seq, 0.0);
	
	std::vector<std::pair<uint32_t, uint32_t>> AND_OR_params;

	// Programmatically fill ANDOR_params reusing the ANN B/R range globals
	BIOLSHASHER_ASSERT(g_start_B > 0 && g_start_R > 0 && g_MAX_B >= g_start_B && g_MAX_R >= g_start_R, 
			"Invalid AND/OR parameter ranges. Please ensure g_start_B, g_start_R, g_MAX_B, and g_MAX_R are set correctly.");
	

	for (uint32_t and_val = g_start_B; and_val <= g_MAX_B; ++and_val) {
		for (uint32_t or_val = g_start_R; or_val <= g_MAX_R; ++or_val) {
			AND_OR_params.push_back({and_val, or_val});
		}
	}
	printf("Generated %zu (AND,OR) pairs for testing (AND from %u to %u, OR from %u to %u).\n", AND_OR_params.size(), g_start_B, g_MAX_B, g_start_R, g_MAX_R);

	// std::vector<std::pair<uint32_t, uint32_t>> ANDOR_params = {{2,3}};//{{1,1},{2,1},{3,2},{1,2},{2,3},{2,2}}; //{{5,5}, {11,11}, {21,21}};
	
	for(auto andor_param : AND_OR_params){
		uint32_t and_param = andor_param.first;
		uint32_t or_param = andor_param.second;
		
		if (and_param == 0 || or_param == 0) {
			printf("Invalid AND/OR parameters: AND=%u, OR=%u\n", and_param, or_param);
			return false;
		}

		std::atomic<uint32_t> completed(0);

		printf("Running collision computation for AND=%u, OR=%u...\n", and_param, or_param);
		
		if (g_NCPU == 1) {	// Single threaded fallback
			printf("Starting collision computation with %u thread...\n", g_NCPU);

			auto start_seq = std::chrono::high_resolution_clock::now();

			std::thread progress_thread(run_progress_bar, std::ref(completed),N_seq);
			LSHCollisionTestInnerInnerWorker<hashtype>(hinfo, N_hash, hash, HashSeed, common_params, sequenceRecordsforTest, AverageCollision.data(), 0, N_seq, and_param, or_param, &completed);

			auto end_seq = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> elapsed_seq = end_seq - start_seq;
			printf("Sequential collision computation took %.3f seconds.\n", elapsed_seq.count());
		}
		else{
			#if defined(HAVE_THREADS)
				printf("Starting collision computation with %u threads...\n", g_NCPU);
				std::vector<std::thread> threads(g_NCPU);
				int block_size = N_seq / g_NCPU;
				// int block_size = std::max(1, N_seq / (g_NCPU * 8));  // 8x finer granularity
				printf("Block size per thread: %d\n", block_size);

				auto start_par = std::chrono::high_resolution_clock::now();

				std::thread progress_thread(run_progress_bar, std::ref(completed), N_seq);
			
				for(uint32_t i = 0; i < g_NCPU; i++){
					int start = i * block_size;
					int end = (i == g_NCPU - 1) ? N_seq : (i + 1) * block_size;

					threads[i] = std::thread(LSHCollisionTestInnerInnerWorker<hashtype>, hinfo, N_hash, hash, HashSeed, std::ref(common_params), std::ref(sequenceRecordsforTest), AverageCollision.data(), start, end, and_param, or_param, &completed);
				}

				for (auto& t : threads) {
					t.join();
				}
				progress_thread.join();
				
				auto end_par = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed_par = end_par - start_par;
				printf("Parallel collision computation took %.3f seconds.\n", elapsed_par.count());
			#endif
		}

		out_file << ":10:" << "AND," << "OR" << std::endl;
		out_file << ":11:" << and_param << "," << or_param << std::endl;
		// Print Average Collision values
		out_file << ":12:"; // Starting  line with collision curves
		for (size_t i = 0; i < N_seq; i++) {
			if (i == N_seq - 1)
				out_file << AverageCollision[i] << "\n";
			else
				out_file << AverageCollision[i] << ",";
		}
	}

	// Write the sequences to a new output file.

	// writeSequencesToFile(sequenceRecordsforTest, AverageCollision, "../results/lsh_collision_sequences.txt", 2);

	return true;
}

static sim_bins_struct LSHCollisionTestInnerAgg(const HashInfo * hinfo, uint32_t N_agg, common_params_struct &common_params, SeedGenerator &seedGen, std::ofstream &out_file){
	
	printf("Inside Mutation configuration Aggregation Step\n");

	SequenceRecordsWithMetadataStruct sequenceRecordsForAgg;

	sequenceRecordsForAgg.OriginalSequenceLength = common_params.seqLen;
	sequenceRecordsForAgg.areBasesDrawnFromUniformDist = common_params.areBasesDrawnFromUniformDist;
	sequenceRecordsForAgg.DataGenSeed = common_params.DataGenSeed;
	sequenceRecordsForAgg.DataMutateSeed = common_params.DataMutateSeed;
	sequenceRecordsForAgg.KeyCount = N_agg;

	[[maybe_unused]] SequenceDataGenerator dataGenAgg(&sequenceRecordsForAgg);

	sim_bins_struct sim_bins_agg;	//Bin to store the error parameters

    std::vector<double> rand_error_param(N_agg, 0.0);
	std::vector<double> similarity_values(N_agg, 0.0);
	
	if(g_mutation_model == MUTATION_MODEL_SIMPLE_SNP_ONLY){
		[[maybe_unused]] SequenceDataMutatorSubstitutionOnly dataMutAgg(&sequenceRecordsForAgg, &rand_error_param, hinfo->similarityfn);
	}
	else if(g_mutation_model == MUTATION_MODEL_GEOMETRIC_MUTATOR){
		printf("Using geometric mutator model for data mutation in aggregation phase.\n");
		[[maybe_unused]] SequenceDataMutatorGeometric dataMutAgg(&sequenceRecordsForAgg, &rand_error_param, hinfo->similarityfn);
	}
	
	// Print Similarity values
	out_file << ":13:";
	for (size_t i = 0; i < N_agg; i++) {
		if (i == N_agg - 1)
			out_file << rand_error_param[i] << "\n";
		else
			out_file << rand_error_param[i] << ",";
	}
	
	// Extract similarity values
	for(uint32_t idx = 0; idx < N_agg; idx++){
		similarity_values[idx] = sequenceRecordsForAgg.Records[idx].similarity;
	}

	// Perform binning on similarity values
	for (size_t i = 0; i < N_agg; i++) {
		float sim_value = similarity_values[i];
		uint32_t bin_index = static_cast<uint32_t>((sim_value - sequenceRecordsForAgg.binstart) / sequenceRecordsForAgg.binsize);
		if (bin_index >= sequenceRecordsForAgg.bincount) {
			bin_index = sequenceRecordsForAgg.bincount - 1; // Clamp to last bin
		}
		
		uint32_t bin_fill_count = sim_bins_agg.bin_fill_count[bin_index];

		if(bin_fill_count == g_bincount_full){
			continue;	// To avoid excessive memory usage, we limit to x entries per bin.
		}
		
		sim_bins_agg.bin_error_parameters[bin_index][bin_fill_count] = rand_error_param[i];
		sim_bins_agg.bin_fill_count[bin_index]++;
	}

	//Compute mean and stddev for each bin
	for (size_t bin_idx = 0; bin_idx < sequenceRecordsForAgg.bincount; bin_idx++) {
		const auto & params = sim_bins_agg.bin_error_parameters[bin_idx];
		if (!params.empty()) {
			// Compute mean
			double sum = 0.0;
			for(uint32_t param_idx = 0; param_idx < sim_bins_agg.bin_fill_count[bin_idx]; param_idx++){
				sum += params[param_idx];
			}
			// printf("Bin %zu: Sum = %f, Count = %u\n", bin_idx, sum, sim_bins_agg.bin_fill_count[bin_idx]);
			double mean = sum / sim_bins_agg.bin_fill_count[bin_idx];
			sim_bins_agg.bin_error_parameters_mean[bin_idx] = mean;

			// Compute stddev
			double sq_sum = 0.0;
			for (uint32_t param_idx = 0; param_idx < sim_bins_agg.bin_fill_count[bin_idx]; param_idx++) {
				sq_sum += (params[param_idx] - mean) * (params[param_idx] - mean);
			}
			double stddev = std::sqrt(sq_sum / sim_bins_agg.bin_fill_count[bin_idx]);
			sim_bins_agg.bin_error_parameters_stddev[bin_idx] = stddev;
			
		} else {
			sim_bins_agg.bin_error_parameters_mean[bin_idx] = 0.0;
			sim_bins_agg.bin_error_parameters_stddev[bin_idx] = 0.0;
		}
	}

	
	// // Print each bin's fill count
	// for (size_t bin_idx = 0; bin_idx < sequenceRecordsForAgg.bincount; bin_idx++) {
	// 	printf("Bin %zu: Fill Count = %u", bin_idx, sim_bins_agg.bin_fill_count[bin_idx]);
	// 	// Print mean and stddev
	// 	printf(", Mean = %0.2f, Stddev = %0.2f\n", sim_bins_agg.bin_error_parameters_mean[bin_idx], sim_bins_agg.bin_error_parameters_stddev[bin_idx]);
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
	// Now, I want to go through each bin, and if its fill count is less than g_bincount_full,
	// I want to fill the remaining entries with random samples from a normal distribution
	// defined by the mean and stddev of that bin.
	
	std::mt19937 gen(seedGen.nextSeed());  // Seeded RNG for reproducibility
	
	for (size_t bin_idx = 0; bin_idx < sequenceRecordsForAgg.bincount; bin_idx++) {
		uint32_t current_fill = sim_bins_agg.bin_fill_count[bin_idx];
		
		if (current_fill < g_bincount_full && current_fill > 0) {
			// Get the mean and stddev for this bin
			double mean = sim_bins_agg.bin_error_parameters_mean[bin_idx];
			double stddev = sim_bins_agg.bin_error_parameters_stddev[bin_idx];
			
			// If stddev is 0 or very small, use a small default to avoid degenerate distribution
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

			//Todo: Add the part where in angular similarity, we geniuinely have empty bins from 0 to 49


			// // Bin is completely empty - use bin center as mean with small stddev
			// double bin_center = sequenceRecordsForAgg.binstart + (bin_idx + 0.5) * sequenceRecordsForAgg.binsize;
			// double default_stddev = sequenceRecordsForAgg.binsize / 4.0;  // Quarter of bin width
			
			// std::normal_distribution<double> normal_dist(bin_center, default_stddev);
			
			// for (uint32_t fill_idx = 0; fill_idx < g_bincount_full; fill_idx++) {
			// 	double sampled_value = normal_dist(gen);
			// 	sampled_value = std::max(0.0, std::min(1.0, sampled_value));
				
			// 	sim_bins_agg.bin_error_parameters[bin_idx][fill_idx] = sampled_value;
			// 	sim_bins_agg.bin_fill_count[bin_idx]++;
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
	sequenceRecordsForAgg.Records.clear();
	sequenceRecordsForAgg.Records.shrink_to_fit();  // Actually releases the memory
	
	return sim_bins_agg;
}

template <typename hashtype>
static bool LSHCollisionTestInner( const HashInfo * hinfo, const uint32_t seqLen, flags_t flags, std::ofstream &out_file, SeedGenerator &seedGen) {

	bool result = true;	//TODO: Update this based on test results.

	HashFn hash = hinfo->hashFn();

	printf("Running LSH Collision Test for hash: %s\n", hinfo->name);
	printf("Sequence Length: %u\n", seqLen);
	printf("Distance Metric: %s\n", hinfo->similarity_name);
	printf("Mutation Model: %s\n", getMutationModelName(g_mutation_model));
	//TODO: Print mutation expression only when valid.
	if(g_mutation_model == MUTATION_MODEL_GEOMETRIC_MUTATOR){
		printf("Mutation Expression Type: %s\n", getMutationExpressionTypeName(g_mutation_expression_type));
	}
	printf("----------------------------------------------\n");

	// File header
	out_file << "\n:1:LSH Collision Test Results\n";
	out_file << ":2:" << "Hashname," << "SequenceLength,"<< "Distance Metric," << "Mutation Model,"<< "Mutation Expression" << std::endl;
	out_file << ":3:" << hinfo->name << "," << seqLen << "," << hinfo->similarity_name << "," << g_mutation_model << "," << g_mutation_expression_type << std::endl;
	
	hinfo->printParameters(out_file); //~TODO:PARAMETERS IN THE test FILE. //DONE!

	seed_t DataGenSeed = seedGen.nextSeed();
	seed_t DataMutateSeed = seedGen.nextSeed();
	seed_t HashSeed = hinfo->Seed(seedGen.nextSeed());

	common_params_struct common_params;
	common_params.seqLen = seqLen;
	common_params.DataGenSeed = DataGenSeed;
	common_params.DataMutateSeed = DataMutateSeed;
	common_params.areBasesDrawnFromUniformDist = g_areBasesDrawnFromUniformDistribution;

	uint32_t N_agg = 0;

	uint32_t N_seq = 0;		// Number of sequences to generate for testing
	uint32_t N_hash = 0;	// Number of hashes to compute per sequence

	//--------------------------------------------//
	if(hinfo->isVerySlow()){
		printf("Hash %s is marked as very slow. Limiting test parameters for practicality.\n", hinfo->name);
		N_agg = g_verySlowNAggCases;	// Number of sequences to generate for testing

		N_seq = g_verySlowNSeq;		// Number of sequences to generate for testing
		N_hash = g_verySlowNHashes;	// Number of hashes to compute per sequence
	}
	else if (hinfo->isSlow()) {
		printf("Hash %s is marked as slow. Limiting test parameters for practicality.\n", hinfo->name);
		N_agg = g_SlowNAggCases;	// Number of sequences to generate for testing

		N_seq = g_SlowNSeq;		// Number of sequences to generate for testing
		N_hash = g_SlowNHashes;	// Number of hashes to compute per sequence
	}
	else{
		N_agg = g_NAggCases;	// Number of sequences to generate for testing

		N_seq = g_NSeq;		// Number of sequences to generate for testing
		N_hash = g_NHashes;	// Number of hashes to compute per sequence
	}

	sim_bins_struct sim_bins = LSHCollisionTestInnerAgg(hinfo, N_agg, common_params, seedGen, out_file);

	//print bin means and stddevs using
	// if (!REPORT(VERBOSE, flags)) {
	// 	for (size_t bin_idx = 0; bin_idx < sim_bins.bin_error_parameters_mean.size(); bin_idx++) {
	// 		printf("Bin %zu: Count %d, Mean = %0.2f, Stddev = %0.2f\n", bin_idx, sim_bins.bin_fill_count[bin_idx], sim_bins.bin_error_parameters_mean[bin_idx], sim_bins.bin_error_parameters_stddev[bin_idx]);
	// 	}
	// }

	LSHCollisionTestInnerInnerParallel<hashtype>(hinfo, N_seq, N_hash, hash, HashSeed, common_params, sim_bins, seedGen, out_file);

    return result;	//TODO: For now, the result is always true. We need to add logic to find where the test fails.
}

//----------------------------------------------------------------------------//
template <typename hashtype>
bool LSHCollisionTest( const HashInfo * hinfo, flags_t flags) {
	printf("[[[ LSH Collision Tests ]]]\n\n");

	BIOLSHASHER_ASSERT(hinfo->isLocalSensitive(),"Flag FLAG_HASH_LOCALITY_SENSITIVE not found. Please ensure that your LSH hash function is defined with FLAG_HASH_LOCALITY_SENSITIVE tag. \n\t\t It may be possible that the provided hash is not an LSH function.");
	BIOLSHASHER_ASSERT(hinfo->similarityfn != nullptr, "LSH function should be defined with a similarity metric.");
	BIOLSHASHER_ASSERT((std::strlen(hinfo->similarity_name) > 0) , "The similarity metric should be named.");

	bool result = true;

	// Create results directory if it doesn't exist
	std::string results_dir_str = "../results";
	std::filesystem::path results_dir = results_dir_str;
	if (!std::filesystem::exists(results_dir)) {
		std::filesystem::create_directories(results_dir);
	}

	// std::string filename = results_dir_str + "/collisionResults_" + std::string(hinfo->name) + ".csv";
	std::string filename = results_dir_str + "/collisionResults_" + hinfo->similarity_name + ".csv";

	std::ios_base::openmode mode = std::ios::trunc;  // Default: replace
	if (std::filesystem::exists(filename)) {
		mode = std::ios::app; // File exists: append
	}
	std::ofstream out_file(filename, mode);
	if (!out_file.is_open()) {
		throw std::runtime_error("Error: Could not open output file");
	}

	uint32_t sequenceLength = g_SequenceLength;

	SeedGenerator seedGen(g_GoldenRatio ^ std::chrono::system_clock::now().time_since_epoch().count());

	SetIsTestActive(true);

	// printf("\nTesting Hash: %s with keybits: %u (sequence length: %u)\n", hinfo->name, (sequenceLength * 8), sequenceLength);
	result &= LSHCollisionTestInner<hashtype>(hinfo, sequenceLength, flags, out_file, seedGen);

    SetIsTestActive(false); // Cleanup: Reset LSH global variables after test completion

    printf("%s\n", result ? "" : g_failstr);
	out_file.close();

	return result;
}

INSTANTIATE(LSHCollisionTest, HASHTYPELIST);


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

        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        // ticks++;
    }
}
