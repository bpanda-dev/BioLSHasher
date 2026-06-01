
#include "BioDataGeneration.h"
#include <random>
#include <set>

#include "HashInfo.h"
#include "assertMsg.h"


bool simulateSNP(SequenceRecordUnit &record, const uint32_t pos, std::mt19937& rng_nuc) {
	
	char currentBase = record.SeqASCIIMut[pos];
	uint8_t currentBaseid = (currentBase >> 1) & 0x03; // Map ASCII to 2-bit (A=0,C=1,T=2,G=3)

	// printf("SNP at pos %u: current base %c (id %u)\n", pos, currentBase, currentBaseid);

	// Randomly select a new base
	// Generate a different base (ensure it's not the same as original)

	std::uniform_int_distribution<uint32_t> dist(0, 2); // TODO: Check if (0,3) works or (0,2)? Study.
	uint32_t newBaseid_int = dist(rng_nuc);
    auto newBaseid = static_cast<uint8_t>(newBaseid_int & 0x03);
	if (currentBaseid <= newBaseid)
		newBaseid += 1;

	// printf("newBaseid_int = %u\n", newBaseid);
	BIOLSHASHER_ASSERT(currentBaseid != newBaseid,"New base should be different from current base");
	BIOLSHASHER_ASSERT(newBaseid <= 3 , "Base ID should be between 0 and 3");

	// Apply SNP
	record.SeqASCIIMut[pos] = bases[newBaseid];
	
	return true;
}

// double ComputeHammingSimilarity(const std::string& seq1, const std::string& seq2, const uint32_t length) {
//     // Hamming distance only defined for equal-length sequences
//     uint32_t similar = 0;
// 	// Count mismatches in overlapping region
// 	for (uint32_t i = 0; i < length; i++) {
// 		if (seq1[i] == seq2[i]) {
// 			similar++;
// 		}
// 	}
// 	return (static_cast<float>(similar) / static_cast<float>(length));
// }

// std::set<std::string> get_kmers(const std::string& sequence, int k) {
//     std::set<std::string> kmers;
//
//     // Ensure the sequence is long enough for at least one k-mer
//     if (sequence.length() < static_cast<size_t>(k)) {
//         return kmers;
//     }
//
//     // Slide the window across the sequence
//     for (size_t i = 0; i <= sequence.length() - k; ++i) {
//         kmers.insert(sequence.substr(i, k));
//     }
//
//     return kmers;
// }

// double ComputeJaccardSimilarity(const std::string& seq1, const std::string& seq2, int k) {
//     // 1. Tokenize sequences into sets of k-mers
//     std::set<std::string> set1 = get_kmers(seq1, k);
//     std::set<std::string> set2 = get_kmers(seq2, k);
//
// 	// printf("Computing Jaccard Similarity for k=%d\n", k);
//
// 	// printf("Set1 size: %zu, Set2 size: %zu\n", set1.size(), set2.size());
// 	// //print the sets for debugging
// 	// // --- Printing Section ---
//     // std::cout << "--- K-mer Analysis (k=" << k << ") ---" << std::endl;
//     // std::cout << "Set1 (" << set1.size() << " unique): ";
//     // for (const auto& kmer : set1) std::cout << kmer << " ";
//
//     // std::cout << "\nSet2 (" << set2.size() << " unique): ";
//     // for (const auto& kmer : set2) std::cout << kmer << " ";
//     // std::cout << "\n" << std::endl;
//
//
//     if (set1.empty() || set2.empty()) return 0.0;
//
//     // 2. Find the intersection
//     std::vector<std::string> intersect;
//     std::set_intersection(set1.begin(), set1.end(),
//                           set2.begin(), set2.end(),
//                           std::back_inserter(intersect));
//
//     // 3. Use the Inclusion-Exclusion principle for the Union size
//     // |A ∪ B| = |A| + |B| - |A ∩ B|
//     size_t intersection_size = intersect.size();
//     size_t union_size = set1.size() + set2.size() - intersection_size;
//
//     return static_cast<double>(intersection_size) / union_size;
// }

// double ComputeCosineSimilarity(const std::string& seq1, const std::string& seq2, int k) {
// 	// 1. Tokenize sequences into sets of k-mers
// 	std::set<std::string> set1 = get_kmers(seq1, k);
// 	std::set<std::string> set2 = get_kmers(seq2, k);
//
// 	if (set1.empty() || set2.empty()) return 0.0;
//
// 	// 2. Find the intersection
// 	std::vector<std::string> intersect;
// 	std::set_intersection(set1.begin(), set1.end(),
// 						  set2.begin(), set2.end(),
// 						  std::back_inserter(intersect));
//
// 	size_t intersection_size = intersect.size();
//
// 	// Cosine similarity = |A ∩ B| / sqrt(|A| * |B|)
// 	return static_cast<double>(intersection_size) /
// 		   std::sqrt(static_cast<double>(set1.size()) * static_cast<double>(set2.size()));
// }

// double ComputeAngularSimilarity(const std::string& seq1, const std::string& seq2, int k) {
// 	double cosine_sim = ComputeCosineSimilarity(seq1, seq2, k);
// 	// Angular similarity = 1 - (angle / π) where angle = arccos(cosine_similarity)
// 	return 1 - (std::acos(cosine_sim) / M_PI);
// }

// UnionBitVectorsStruct CreateUnionBitVectors(const std::string& seq1, const std::string& seq2, int k) {
// 	// Tokenize sequences into sets of k-mers
// 	std::set<std::string> kmer_set_a = get_kmers(seq1, k);
// 	std::set<std::string> kmer_set_b = get_kmers(seq2, k);
//
// 	if (kmer_set_a.empty() || kmer_set_b.empty()) return {std::vector<char>(), std::vector<char>(), std::vector<std::string>()};	//return empty vectors if either set is empty
//
//     // 1. Create the Universe (Union of both sets)
//     // std::set_union requires sorted ranges. std::set is already sorted.
//     std::vector<std::string> union_kmers;
//     std::set_union(kmer_set_a.begin(), kmer_set_a.end(),
//                    kmer_set_b.begin(), kmer_set_b.end(),
//                    std::back_inserter(union_kmers));
//
//     // 2. Pre-allocate vectors
//     size_t n = union_kmers.size();
// 	std::vector<char> vec_a(n, '0');
// 	std::vector<char> vec_b(n, '0');
//
//     // 3. Fill the vectors
//     for (size_t i = 0; i < n; ++i) {
//         const std::string& kmer = union_kmers[i];
//
//         // Use set::count to check for existence (O(log N))
//         if (kmer_set_a.count(kmer)) {
//             vec_a[i] = '1';
//         }
//         if (kmer_set_b.count(kmer)) {
//             vec_b[i] = '1';
//         }
//     }
//
//     return {std::move(vec_a), std::move(vec_b), std::vector<std::string>()};
// }

// int min1(int x,int y,int z){
// 	int min_xy = (x<y)? x:y;
//     if(min_xy>z)
//         return z;
//     else
//         return min_xy;
// }

// double ComputeEditSimilarity(const std::string& seq1, const std::string& seq2) {
//
// 	//write your code here
// 	const int n = seq1.size();
// 	const int m = seq2.size();
//
// 	std::vector<std::vector<int>> D(n + 1, std::vector<int>(m + 1, 0));
//
// 	// Initialize base cases
//     for (int i = 0; i <= n; i++) {
//         D[i][0] = static_cast<int>(i);
//     }
//     for (int j = 0; j <= m; j++) {
//         D[0][j] = static_cast<int>(j);
//     }
//
// 	// Fill DP table
//     for (int i = 1; i <= n; i++) {
//         for (int j = 1; j <= m; j++) {
//             int insert_cost = D[i][j - 1] + 1;
//             int delete_cost = D[i - 1][j] + 1;
//             int match_cost  = D[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? 0 : 1);
//
//             D[i][j] = std::min({insert_cost, delete_cost, match_cost});
//         }
//     }
//
// 	int edit_dist = D[n][m];
// 	double max_len = static_cast<double>(std::max(n, m));
//
//     return 1.0 - (static_cast<double>(edit_dist) / max_len);
// }


/*#-----------------------------------------------------#*/
/*#					Data Generation						#*/
/*#-----------------------------------------------------#*/

SequenceDataGenerator::SequenceDataGenerator(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata){
	BIOLSHASHER_ASSERT(sequenceRecordsWithMetadata != nullptr, "Input pointer for sequence records with metadata cannot be null");
	BIOLSHASHER_ASSERT(((sequenceRecordsWithMetadata->A_percentage + sequenceRecordsWithMetadata->C_percentage  + sequenceRecordsWithMetadata->G_percentage + sequenceRecordsWithMetadata->T_percentage) == 1), "Sum of nucleotide probs. should be 1.");
	BIOLSHASHER_ASSERT((sequenceRecordsWithMetadata->OriginalSequenceLength > 0), "Sequence length must be positive");
	BIOLSHASHER_ASSERT(sequenceRecordsWithMetadata->KeyCount>0, "Key count must be positive");

	// Initialise seed and reserve memory.
    std::mt19937 rng(sequenceRecordsWithMetadata->DataGenSeed);

	if(sequenceRecordsWithMetadata->areBasesDrawnFromUniformDist == false){
		// If not drawn from uniform distribution, we can set custom base percentages here.
		sequenceRecordsWithMetadata->A_percentage = 1; //0.50;
		sequenceRecordsWithMetadata->C_percentage = 0.00;
		sequenceRecordsWithMetadata->G_percentage = 0.00;
		sequenceRecordsWithMetadata->T_percentage = 0.00;
	}
	else{	// Need to mention this path of branch too, since we have an assertion for sum of percentages.
		sequenceRecordsWithMetadata->A_percentage = 0.25;
		sequenceRecordsWithMetadata->C_percentage = 0.25;
		sequenceRecordsWithMetadata->G_percentage = 0.25;
		sequenceRecordsWithMetadata->T_percentage = 0.25;
	}

    sequenceRecordsWithMetadata->Records.resize(sequenceRecordsWithMetadata->KeyCount);	// Allocate memory for sequence records

    uint32_t total_bases = sequenceRecordsWithMetadata->OriginalSequenceLength * sequenceRecordsWithMetadata->KeyCount;
    std::vector<uint8_t> rand_data(total_bases, 0);

	size_t i = 0;
	while (i + 4 <= total_bases) {
		uint32_t val = rng();          // 32-bit output
		memcpy(&rand_data[i], &val, 4);
		i += 4;
	}
	if (i < total_bases) {
		uint32_t val = rng();
		memcpy(&rand_data[i], &val, total_bases - i);
	}
	
	if(sequenceRecordsWithMetadata->areBasesDrawnFromUniformDist == 1){
		for (uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {
			auto& record = sequenceRecordsWithMetadata->Records[rec_idx];

			// setting the size of the sequences
			record.OriginalLength = sequenceRecordsWithMetadata->OriginalSequenceLength;

			record.SeqASCIIOrg.resize(sequenceRecordsWithMetadata->OriginalSequenceLength, 'A');

			for (uint32_t pos = 0; pos < sequenceRecordsWithMetadata->OriginalSequenceLength; pos++) {
				
				uint8_t base = rand_data[(rec_idx * sequenceRecordsWithMetadata->OriginalSequenceLength) + pos] % 4;
				
				record.SeqASCIIOrg[pos] = bases[base];
			}

			//this was for testing
			// for(uint32_t pos = 0; pos < sequenceRecordsWithMetadata->OriginalSequenceLength; pos++) {
			// 	record.SeqASCIIOrg = "TCGAATGCTCCTACACACTGTTTCGAACATGATCCGCCCGGAGGC";
			// }

			/* Debugging Statements! Please never delete!*/
			// Print first 10 sequences for verification
			// if (rec_idx < 10) {
			// 	std::cout << "Record " << rec_idx << " Original Sequence: " << record.SeqASCIIOrg << std::endl;
			// }
		}
	}
	else{
		// Precompute cumulative probabilities
		std::vector<double> cumulative_probs = {
			sequenceRecordsWithMetadata->A_percentage,
			sequenceRecordsWithMetadata->A_percentage + sequenceRecordsWithMetadata->C_percentage,
			sequenceRecordsWithMetadata->A_percentage + sequenceRecordsWithMetadata->C_percentage + sequenceRecordsWithMetadata->G_percentage,
			1.0 // T percentage
		};

		for (uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {
			auto& record = sequenceRecordsWithMetadata->Records[rec_idx];

			// setting the size of the sequences
			record.OriginalLength = sequenceRecordsWithMetadata->OriginalSequenceLength;

			record.SeqASCIIOrg.resize(sequenceRecordsWithMetadata->OriginalSequenceLength, 'A');

			for (uint32_t pos = 0; pos < sequenceRecordsWithMetadata->OriginalSequenceLength; pos++) {
				
				double rand_val = static_cast<double>(rand_data[(rec_idx * sequenceRecordsWithMetadata->OriginalSequenceLength) + pos]) / 256.0;	//CONVERTED INTO 0 TO 1

				uint8_t base = 0;
				if (rand_val < cumulative_probs[0]) {
					base = 0; // A
				} else if (rand_val < cumulative_probs[1]) {
					base = 1; // C
				} else if (rand_val < cumulative_probs[2]) {
					base = 3; // G
				} else {
					base = 2; // T
				}
				
				record.SeqASCIIOrg[pos] = bases[base];
				
			}

			/* Debugging Statements! Please never delete!*/
			// Print first 10 sequences for verification
			// if (rec_idx < 10) {
			// 	std::cout << "Record " << rec_idx << " Original Sequence: " << record.SeqASCIIOrg << std::endl;
			// }
		}
	}
    

	// std::cout << "DataGeneration: SequenceLength = " << sequenceRecordsWithMetadata->OriginalSequenceLength << ", KeyCount = " << sequenceRecordsWithMetadata->KeyCount << std::endl;
}

void SequenceDataGenerator::decodeSequencesToASCII(const std::vector<uint8_t>& seqTwoBit, std::string& seqASCII, const uint32_t sequenceLength) {
    BIOLSHASHER_ASSERT(seqASCII.size() == sequenceLength, "ASCII sequence size does not match expected length");
    for (uint32_t pos = 0; pos < sequenceLength; pos++) {
        uint32_t byte_idx = (pos * two_bit_per_base) / 8;
        uint32_t bit_offset = (pos * two_bit_per_base) % 8;
        uint8_t base = (seqTwoBit[byte_idx] >> bit_offset) & 0x03;
        seqASCII[pos] = bases[base];
    }
}

/*#-----------------------------------------------------#*/
/*#					Data Mutation						#*/
/*#-----------------------------------------------------#*/

// For Agg
SequenceDataMutatorSubstitutionOnly::SequenceDataMutatorSubstitutionOnly(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata, std::vector<double> *rand_error_param, SimilarityFn similarity_fn){
	BIOLSHASHER_ASSERT(sequenceRecordsWithMetadata != nullptr, "Input pointer for sequence records with metadata cannot be null");
	BIOLSHASHER_ASSERT(!sequenceRecordsWithMetadata->Records.empty(), "Sequence records with metadata cannot be empty");
	BIOLSHASHER_ASSERT(similarity_fn != nullptr, "Similarity metric is required for mutation. Please Add Appropriate Similarity metric in your LSH file.");

	// Initialise seed .
	std::mt19937 rng_snp_rate(sequenceRecordsWithMetadata->DataMutateSeed);
    std::mt19937 rng_mut(sequenceRecordsWithMetadata->DataMutateSeed + 3);
	std::mt19937 rng_nuc(sequenceRecordsWithMetadata->DataMutateSeed + 7);

	std::uniform_int_distribution<uint32_t> dist_snp_rate(0, 1000);
	std::uniform_int_distribution<uint32_t> dist_mut(0, 1000);


	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {

		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
		
		uint32_t snp_value = dist_snp_rate(rng_snp_rate);;

		double snp_rate = static_cast<double>(snp_value) / static_cast<double>(1000);	// snp or indel
		rand_error_param->at(rec_idx) = snp_rate;

		sequenceRecordsWithMetadata->Records[rec_idx].foundationalParameter = snp_rate;
		sequenceRecordsWithMetadata->Records[rec_idx].snpRate = snp_rate;
		sequenceRecordsWithMetadata->Records[rec_idx].delRate = 0.0;	// No indels in this mutator
		sequenceRecordsWithMetadata->Records[rec_idx].stayRate = 1 - snp_rate;	// No indels in this mutator
		sequenceRecordsWithMetadata->Records[rec_idx].insmean = 0.0;	// No indels in this mutator
		sequenceRecordsWithMetadata->Records[rec_idx].insRate = 0.0;	// No indels in this mutator
		
		

		// Initialising the mutated sequences with original sequences
		record.MutatedLength = record.OriginalLength;
		record.SeqASCIIMut.resize(record.MutatedLength, 'A');
		record.SeqASCIIMut = record.SeqASCIIOrg;	// Copy original ASCII sequence

		// NOTE: Use a while loop or track position carefully since length changes
		for (int pos = static_cast<int>(record.OriginalLength) - 1; pos >= 0; --pos) {

			uint32_t rand_val = dist_mut(rng_mut);

			double sample = static_cast<double>(rand_val) / 1000.0;	// snp or indel
			if (sample < record.snpRate) {
                simulateSNP(record, pos, rng_nuc);
            }
        }
	}

	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {
		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
		record.similarity = similarity_fn(record.SeqASCIIOrg, record.SeqASCIIMut, record.OriginalLength, record.MutatedLength);
	}

	// Debugging output
    // for(uint32_t rec_idx = 0; rec_idx < std::min(sequenceRecordsWithMetadata->KeyCount, 10u); rec_idx++) {
    //     auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
    //     if (rec_idx < 10) {
    //         std::cout << "\n=== Record " << rec_idx << " ===" << std::endl;
    //         std::cout << "Original Length: " << record.OriginalLength << std::endl;
    //         std::cout << "Mutated Length:  " << record.MutatedLength << std::endl;
    //         std::cout << "Original: " << record.SeqASCIIOrg << std::endl;
    //         std::cout << "Mutated:  " << record.SeqASCIIMut << std::endl;
	// 		std::cout << "Foundational Rate: " << record.foundationalParameter << std::endl;
	// 		std::cout << "Similarity: " << record.similarity << std::endl;
    //     }
    // }	
}

// For test
SequenceDataMutatorSubstitutionOnly::SequenceDataMutatorSubstitutionOnly(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata, SimilarityFn similarity_fn){
	BIOLSHASHER_ASSERT(sequenceRecordsWithMetadata != nullptr, "Input pointer for sequence records with metadata cannot be null");
	BIOLSHASHER_ASSERT(!sequenceRecordsWithMetadata->Records.empty(), "Sequence records with metadata cannot be empty");
	BIOLSHASHER_ASSERT(similarity_fn != nullptr, "Similarity metric is required for mutation. Please Add Appropriate Similarity metric in your LSH file.");

	// Initialise seed and reserve memory.
    std::mt19937 rng_mut(sequenceRecordsWithMetadata->DataMutateSeed);
	std::mt19937 rng_nuc(sequenceRecordsWithMetadata->DataMutateSeed + 7);

	std::uniform_int_distribution<uint32_t> dist_mut(0, 1000);

	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {

		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
		
		double snp_rate = record.foundationalParameter;
		record.snpRate = snp_rate;
		record.delRate = 0.0;	// No indels in this mutator
		record.stayRate = 1 - snp_rate;	// No indels in this mutator
		record.insmean = 0.0;	// No indels in this mutator
		record.insRate = 0.0;	// No indels in this mutator

		// Initialising the mutated sequences with original sequences
		record.MutatedLength = record.OriginalLength;
		record.SeqASCIIMut.resize(record.MutatedLength, 'A');
		record.SeqASCIIMut = record.SeqASCIIOrg;	// Copy original ASCII sequence

		// NOTE: Use a while loop or track position carefully since length changes
		for (int pos = static_cast<int>(record.OriginalLength) - 1; pos >= 0; --pos) {
			uint32_t rand_val = dist_mut(rng_mut);
			double sample = static_cast<double>(rand_val) / 1000.0;	// snp or indel
			if (sample < record.snpRate) {
                simulateSNP(record, pos, rng_nuc);
            }
        }
	}

	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {
		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
		record.similarity = similarity_fn(record.SeqASCIIOrg, record.SeqASCIIMut, record.OriginalLength, record.MutatedLength);
	}
	// std::cout << "DataMutation: Applied SNP mutations only." << std::endl;

	// Debugging output
    // for(uint32_t rec_idx = 0; rec_idx < std::min(sequenceRecordsWithMetadata->KeyCount, 10u); rec_idx++) {
    //     auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
    //     if (rec_idx < 10) {
    //         std::cout << "\n=== Record " << rec_idx << " ===" << std::endl;
    //         std::cout << "Original Length: " << record.OriginalLength << std::endl;
    //         std::cout << "Mutated Length:  " << record.MutatedLength << std::endl;
    //         std::cout << "Original: " << record.SeqASCIIOrg << std::endl;
    //         std::cout << "Mutated:  " << record.SeqASCIIMut << std::endl;
	// 		std::cout << "Foundational Parameter: " << record.foundationalParameter << std::endl;
	// 		std::cout << "Similarity: " << record.similarity << std::endl;
    //     }
    // }	
}


// For Agg
SequenceDataMutatorGeometric::SequenceDataMutatorGeometric(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata, std::vector<double> *rand_error_param, SimilarityFn similarity_fn){
	BIOLSHASHER_ASSERT(sequenceRecordsWithMetadata != nullptr, "Input pointer for sequence records with metadata cannot be null");
	BIOLSHASHER_ASSERT(rand_error_param != nullptr, "Input pointer for random error parameters cannot be null");
	BIOLSHASHER_ASSERT(!sequenceRecordsWithMetadata->Records.empty(), "Sequence records with metadata cannot be empty");
	BIOLSHASHER_ASSERT(similarity_fn != nullptr, "Similarity metric is required for mutation. Please Add Appropriate Similarity metric in your LSH file.");

	// Initialise seed and reserve memory.
    std::mt19937 rng_gmean(sequenceRecordsWithMetadata->DataMutateSeed);	// For drawing the geometric mean.
	std::mt19937 rng_sub_del((sequenceRecordsWithMetadata->DataMutateSeed + 11));	// For drawing the substitution and deletion marks.
	std::mt19937 gen(sequenceRecordsWithMetadata->DataMutateSeed  + 21);	// For drawing the length of the insertion.

	std::uniform_int_distribution<uint32_t> dist_gmean(0, 100); //TODO: Make the upper boundary visible outside. Sample a mean value between 0 and 0.1 for geometric mutator. We can change this range as needed.
	std::uniform_int_distribution<uint32_t> dist_point_mut(0,1000);
	std::uniform_int_distribution<uint8_t> dist_base(0, 3);

	// Mutations to be sampled here.
	uint32_t mutation_expression_type = g_mutation_expression_type; 	// Change the expression type here as needed.

	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {
		double g_mean = static_cast<double>(dist_gmean(rng_gmean))/1000.0;
		double P_sub, P_del;
		mutation_expression(g_mean, mutation_expression_type, &P_sub, &P_del);

		// printf("Record %u: %f Initial P_sub = %f, P_del = %f\n", rec_idx,g_mean, P_sub, P_del);

		while(is_valid_mutation_parameters(P_sub, P_del) == false){
			// Redraw
			g_mean = static_cast<double>(dist_gmean(rng_gmean))/1000.0;
			mutation_expression(g_mean, mutation_expression_type, &P_sub, &P_del);
		}
		rand_error_param->at(rec_idx) = g_mean;
		sequenceRecordsWithMetadata->Records[rec_idx].foundationalParameter = g_mean;
		sequenceRecordsWithMetadata->Records[rec_idx].snpRate = P_sub;
		sequenceRecordsWithMetadata->Records[rec_idx].delRate = P_del;
		sequenceRecordsWithMetadata->Records[rec_idx].stayRate = 1 - P_sub - P_del;
		sequenceRecordsWithMetadata->Records[rec_idx].insmean = g_mean;
		sequenceRecordsWithMetadata->Records[rec_idx].insRate = (1-(1/(1+g_mean)))*(1/(1+g_mean));  // Probability of insertion being of length 1 drawn from a geometric distribution with mean g_mean.

		uint32_t count_sub = 0;
		uint32_t count_del = 0;
		uint32_t count_ins = 0;
		uint32_t count_stay = 0;
		uint32_t count_ins_length = 0;


		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
		
		record.MutatedLength = record.OriginalLength;	// Initialising the mutated sequences with original sequences
		record.SeqASCIIMut.resize(record.MutatedLength, 'A');
		// record.SeqASCIIMut = record.SeqASCIIOrg;	// Copy original ASCII sequence
		
		// create an array of length of original sequence to mark positions for mutation(1=P_sub,2=P_del,0=P_stay)
		std::vector<uint32_t> mutation_marks(record.OriginalLength, 0);
		std::vector<uint32_t> ins_lengths(record.OriginalLength, 0); // To store insertion lengths at each position

		// Precompute cumulative probabilities
		std::vector<double> cumulative_probs = {
			record.snpRate,	// P_sub
			record.snpRate + record.delRate,	// P_del
			1.0 // P_stay
		};


		// Draw and store if its P_sub, P_del or P_stay for each position
		for (uint32_t pos = 0; pos < record.OriginalLength; pos++) {
			double rand_val = static_cast<double>(dist_point_mut(rng_sub_del)) / 1000.0;	// Convert to 0 to 1
			if (rand_val < cumulative_probs[0]) {
				mutation_marks[pos] = 1; // P_sub
				count_sub++;
			} else if (rand_val < cumulative_probs[1]) {
				mutation_marks[pos] = 2; // P_del
				count_del++;
			} else {
				mutation_marks[pos] = 0; // P_stay
				count_stay++;
			}
		}

		if(g_mutation_expression_type == MUTATION_EXPRESSION_SUB_ONLY){
			for (uint32_t pos = 0; pos < record.OriginalLength; pos++) {
				ins_lengths[pos] = 0;
			}

			count_ins_length = 0;
			count_ins = 0;
		}
		else{
			// Draw insertion lengths for each position based on geometric distribution
			double p = 1.0 / (record.insmean + 1.0);
			std::geometric_distribution<uint32_t> geom_dist(p);
			for (uint32_t pos = 0; pos < record.OriginalLength; pos++) {
				uint32_t ins_length = geom_dist(gen);
				ins_lengths[pos] = ins_length;

				count_ins_length += ins_length;
				count_ins += (ins_length > 0) ? 1 : 0;
				// printf("Record %u, Position %u: Insertion Length = %u\n", rec_idx, pos, ins_length);
			}
		}
		

		// Compute the memory required.
		record.MutatedLength = record.OriginalLength;  // Start with original
		for (uint32_t pos = 0; pos < record.OriginalLength; pos++) {
			if (mutation_marks[pos] == 2) { // Deletion
				record.MutatedLength -= 1;
			} 
			// Insertions
			record.MutatedLength += ins_lengths[pos];
		}
		record.SeqASCIIMut.resize(record.MutatedLength, 'A');	// Resize mutated sequence based on computed length
		
		// if(record.MutatedLength > record.OriginalLength){
		// 	std::cout << "Record " << rec_idx << ": Original Length = " << record.OriginalLength << ", Mutated Length = " << record.MutatedLength << std::endl;
		// }

		// Store the counts in the record for reference
		record.count_sub = count_sub;
		record.count_del = count_del;
		record.count_ins = count_ins;
		record.count_stay = count_stay;
		record.count_ins_length = count_ins_length;


		// Apply mutations based on precomputed marks
		uint32_t mpos = 0;	// Position in mutated sequence
		for(uint32_t pos = 0; pos < record.OriginalLength; pos++) {
			
			if(mutation_marks[pos]==1){
				//SNP
				record.SeqASCIIMut[mpos] = record.SeqASCIIOrg[pos];
				simulateSNP(record, mpos, rng_sub_del);
				mpos++;
			}
			else if(mutation_marks[pos]==2){
				//Deletion
				//Skip this base in mutated sequence
				//mpos is not incremented
			}
			else{	//mutation_marks[pos]==0
				record.SeqASCIIMut[mpos] = record.SeqASCIIOrg[pos];
				mpos++;
			}

			for(uint32_t ins_idx=0; ins_idx<ins_lengths[pos]; ins_idx++){
				//Randomly select a base for insertion
				uint8_t base = dist_base(rng_sub_del);
				record.SeqASCIIMut[mpos] = bases[base];
				mpos++;
			}
		}
		BIOLSHASHER_ASSERT(mpos == record.MutatedLength, "Mismatch between computed and actual mutated length");

	}

	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {
		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
		record.similarity = similarity_fn(record.SeqASCIIOrg, record.SeqASCIIMut, record.OriginalLength, record.MutatedLength);
	}
	std::cout << "DataMutation: Applied Geometric mutations only." << std::endl;

	// Debugging output
    // for(uint32_t rec_idx = 0; rec_idx < std::min(sequenceRecordsWithMetadata->KeyCount, 10u); rec_idx++) {
    //     auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
    //     if (rec_idx < 10) {
    //         std::cout << "\n=== Record " << rec_idx << " ===" << std::endl;
    //         std::cout << "Original Length: " << record.OriginalLength << std::endl;
    //         std::cout << "Mutated Length:  " << record.MutatedLength << std::endl;
    //         std::cout << "Original: " << record.SeqASCIIOrg << std::endl;
    //         std::cout << "Mutated:  " << record.SeqASCIIMut << std::endl;
	// 		std::cout << "Foundational Parameter: " << record.foundationalParameter << std::endl;
	// 		std::cout << "Similarity: " << record.similarity << std::endl;
    //     }
    // }	
}


// For test
SequenceDataMutatorGeometric::SequenceDataMutatorGeometric(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata, SimilarityFn similarity_fn){
	BIOLSHASHER_ASSERT(sequenceRecordsWithMetadata != nullptr, "Input pointer for sequence records with metadata cannot be null");
	BIOLSHASHER_ASSERT(!sequenceRecordsWithMetadata->Records.empty(), "Sequence records with metadata cannot be empty");
	BIOLSHASHER_ASSERT(similarity_fn != nullptr, "Similarity metric is required for mutation. Please Add Appropriate Similarity metric in your LSH file.");

	// Initialise seed and reserve memory.
    std::mt19937 rng_sub_del(sequenceRecordsWithMetadata->DataMutateSeed + 11);
	std::mt19937 gen(sequenceRecordsWithMetadata->DataMutateSeed + 21);	// For geometric distribution

	uint32_t mutation_expression_type = g_mutation_expression_type; 	// Change the expression type here as needed.

	std::uniform_int_distribution<uint32_t> dist_point_mut(0,1000);
	std::uniform_int_distribution<uint8_t> dist_base(0, 3);

	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {

		uint32_t count_sub = 0;
		uint32_t count_del = 0;
		uint32_t count_ins = 0;
		uint32_t count_stay = 0;
		uint32_t count_ins_length = 0;

		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];

		double g_mean = record.foundationalParameter;
		double P_sub, P_del;
		mutation_expression(g_mean, mutation_expression_type, &P_sub, &P_del);

		sequenceRecordsWithMetadata->Records[rec_idx].snpRate = P_sub;
		sequenceRecordsWithMetadata->Records[rec_idx].delRate = P_del;
		sequenceRecordsWithMetadata->Records[rec_idx].stayRate = 1 - P_sub - P_del;
		sequenceRecordsWithMetadata->Records[rec_idx].insmean = g_mean;
		sequenceRecordsWithMetadata->Records[rec_idx].insRate = (1-(1/(1+g_mean)))*(1/(1+g_mean));  // Probability of insertion being of length 1 drawn from a geometric distribution with mean g_mean.
		

		// Initialising the mutated sequences with original sequences
		record.MutatedLength = record.OriginalLength;
		record.SeqASCIIMut.resize(record.MutatedLength, 'A');
		// record.SeqASCIIMut = record.SeqASCIIOrg;	// Copy original ASCII sequence

		// create an array of length of original sequence to mark positions for mutation(1=P_sub,2=P_del,0=P_stay)
		std::vector<uint32_t> mutation_marks(record.OriginalLength, 0);
		std::vector<uint32_t> ins_lengths(record.OriginalLength, 0); // To store insertion lengths at each position

		// Precompute cumulative probabilities
		std::vector<double> cumulative_probs = {
			record.snpRate,	// P_sub
			record.snpRate + record.delRate,	// P_del
			1.0 // P_stay
		};

		// Draw and store if its P_sub, P_del or P_stay for each position
		for (uint32_t pos = 0; pos < record.OriginalLength; pos++) {
			double rand_val = static_cast<double>(dist_point_mut(rng_sub_del)) / 1000.0;	// Convert to 0 to 1
			if (rand_val < cumulative_probs[0]) {
				mutation_marks[pos] = 1; // P_sub
				count_sub++;
			} else if (rand_val < cumulative_probs[1]) {
				mutation_marks[pos] = 2; // P_del
				count_del++;
			} else {
				mutation_marks[pos] = 0; // P_stay
				count_stay++;
			}
		}

		if(g_mutation_expression_type == MUTATION_EXPRESSION_SUB_ONLY){
			for (uint32_t pos = 0; pos < record.OriginalLength; pos++) {
				ins_lengths[pos] = 0;
			}
			count_ins_length = 0;
			count_ins = 0;
		}
		else{
			// Draw insertion lengths for each position based on geometric distribution
			double p = 1.0 / (record.insmean + 1.0);
			std::geometric_distribution<uint32_t> geom_dist(p);
			for (uint32_t pos = 0; pos < record.OriginalLength; pos++) {
				uint32_t ins_length = geom_dist(gen);
				ins_lengths[pos] = ins_length;

				count_ins_length += ins_lengths[pos];
				count_ins += (ins_lengths[pos] > 0) ? 1 : 0;
				// printf("Record %u, Position %u: Insertion Length = %u\n", rec_idx, pos, ins_length);
			}
		}
		
		// Compute the memory required.
		record.MutatedLength = record.OriginalLength;  // Start with original
		for (uint32_t pos = 0; pos < record.OriginalLength; pos++) {
			if (mutation_marks[pos] == 2) { // Deletion
				record.MutatedLength -= 1;
			} 
			// Insertions
			record.MutatedLength += ins_lengths[pos];
		}

		record.SeqASCIIMut.resize(record.MutatedLength, 'A');	// Resize mutated sequence based on computed length
		
		// Store the counts in the record for reference
		record.count_sub = count_sub;
		record.count_del = count_del;
		record.count_ins = count_ins;
		record.count_stay = count_stay;
		record.count_ins_length = count_ins_length;

		// Apply mutations based on precomputed marks
		uint32_t mpos = 0;	// Position in mutated sequence
		for(uint32_t pos = 0; pos < record.OriginalLength; pos++) {
			
			if(mutation_marks[pos]==1){
				//SNP
				record.SeqASCIIMut[mpos] = record.SeqASCIIOrg[pos];
				simulateSNP(record, mpos, rng_sub_del);
				mpos++;
			}
			else if(mutation_marks[pos]==2){
				//Deletion
				//Skip this base in mutated sequence
				//mpos is not incremented
			}
			else{	//mutation_marks[pos]==0
				record.SeqASCIIMut[mpos] = record.SeqASCIIOrg[pos];
				mpos++;
			}

			for(uint32_t ins_idx=0; ins_idx<ins_lengths[pos]; ins_idx++){
				//Randomly select a base for insertion
				uint8_t base = dist_base(rng_sub_del);
				record.SeqASCIIMut[mpos] = bases[base];
				mpos++;
			}
		}
		BIOLSHASHER_ASSERT(mpos == record.MutatedLength, "Mismatch between computed and actual mutated length");

	}
	
	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {
		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
		
		// //TODO: Temporary. Please remove. Just added for debugging.!!!! Purpose is to check implementation in tests when no mutations is applied.
		// sequenceRecordsWithMetadata->Records[rec_idx].SeqASCIIMut = sequenceRecordsWithMetadata->Records[rec_idx].SeqASCIIOrg;
		// sequenceRecordsWithMetadata->Records[rec_idx].MutatedLength = sequenceRecordsWithMetadata->Records[rec_idx].OriginalLength;

		
		record.similarity = similarity_fn(record.SeqASCIIOrg, record.SeqASCIIMut, record.OriginalLength, record.MutatedLength);
	}
	std::cout << "DataMutation: Applied Geometric mutations only." << std::endl;

	// Debugging output
    // for(uint32_t rec_idx = 0; rec_idx < std::min(sequenceRecordsWithMetadata->KeyCount, 10u); rec_idx++) {
    //     auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
    //     if (rec_idx < 10) {
    //         std::cout << "\n=== Record " << rec_idx << " ===" << std::endl;
    //         std::cout << "Original Length: " << record.OriginalLength << std::endl;
    //         std::cout << "Mutated Length:  " << record.MutatedLength << std::endl;
    //         std::cout << "Original: " << record.SeqASCIIOrg << std::endl;
    //         std::cout << "Mutated:  " << record.SeqASCIIMut << std::endl;
	// 		std::cout << "Foundational Parameter: " << record.foundationalParameter << std::endl;
	// 		std::cout << "Similarity: " << record.similarity << std::endl;
    //     }
    // }	
}




// For Agg
// SequenceDataMutatorIndependentRateSID::SequenceDataMutatorIndependentRateSID(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata, std::vector<double> *rand_error_param, SimilarityFn similarity_fn){
// 	BIOLSHASHER_ASSERT(sequenceRecordsWithMetadata != nullptr, "Input pointer for sequence records with metadata cannot be null");
// 	BIOLSHASHER_ASSERT(!sequenceRecordsWithMetadata->Records.empty(), "Sequence records with metadata cannot be empty");
// 	BIOLSHASHER_ASSERT(similarity_fn != nullptr, "Similarity metric is required for mutation. Please Add Appropriate Similarity metric in your LSH file.");

// 	// Initialise seed .
// 	std::mt19937 rng_snp_rate(sequenceRecordsWithMetadata->DataMutateSeed);
//     std::mt19937 rng_mut(sequenceRecordsWithMetadata->DataMutateSeed + 3);
// 	std::mt19937 rng_nuc(sequenceRecordsWithMetadata->DataMutateSeed + 7);

// 	std::uniform_int_distribution<uint32_t> dist_snp_rate(0, 1000);
// 	std::uniform_int_distribution<uint32_t> dist_mut(0, 1000);


// 	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {

// 		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
		
// 		uint32_t snp_value = dist_snp_rate(rng_snp_rate);;

// 		double snp_rate = static_cast<double>(snp_value) / static_cast<double>(1000);	// snp or indel
// 		rand_error_param->at(rec_idx) = snp_rate;

// 		sequenceRecordsWithMetadata->Records[rec_idx].foundationalParameter = snp_rate;
// 		sequenceRecordsWithMetadata->Records[rec_idx].snpRate = snp_rate;
// 		sequenceRecordsWithMetadata->Records[rec_idx].delRate = 0.0;	// No indels in this mutator
// 		sequenceRecordsWithMetadata->Records[rec_idx].stayRate = 1 - snp_rate;	// No indels in this mutator
// 		sequenceRecordsWithMetadata->Records[rec_idx].insmean = 0.0;	// No indels in this mutator
// 		sequenceRecordsWithMetadata->Records[rec_idx].insRate = 0.0;	// No indels in this mutator
		
		

// 		// Initialising the mutated sequences with original sequences
// 		record.MutatedLength = record.OriginalLength;
// 		record.SeqASCIIMut.resize(record.MutatedLength, 'A');
// 		record.SeqASCIIMut = record.SeqASCIIOrg;	// Copy original ASCII sequence

// 		// NOTE: Use a while loop or track position carefully since length changes
// 		for (int pos = static_cast<int>(record.OriginalLength) - 1; pos >= 0; --pos) {

// 			uint32_t rand_val = dist_mut(rng_mut);

// 			double sample = static_cast<double>(rand_val) / 1000.0;	// snp or indel
// 			if (sample < record.snpRate) {
//                 simulateSNP(record, pos, rng_nuc);
//             }
//         }
// 	}

// 	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {
// 		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
// 		record.similarity = similarity_fn(record.SeqASCIIOrg, record.SeqASCIIMut, record.OriginalLength, record.MutatedLength);
// 	}

// 	// Debugging output
//     for(uint32_t rec_idx = 0; rec_idx < std::min(sequenceRecordsWithMetadata->KeyCount, 10u); rec_idx++) {
//         auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
//         if (rec_idx < 10) {
//             std::cout << "\n=== Record " << rec_idx << " ===" << std::endl;
//             std::cout << "Original Length: " << record.OriginalLength << std::endl;
//             std::cout << "Mutated Length:  " << record.MutatedLength << std::endl;
//             std::cout << "Original: " << record.SeqASCIIOrg << std::endl;
//             std::cout << "Mutated:  " << record.SeqASCIIMut << std::endl;
// 			std::cout << "Foundational Rate: " << record.foundationalParameter << std::endl;
// 			std::cout << "Similarity: " << record.similarity << std::endl;
//         }
//     }	
// }

// For test
// SequenceDataMutatorIndependentRateSID::SequenceDataMutatorSubstitutionOnly(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata, SimilarityFn similarity_fn){
// 	BIOLSHASHER_ASSERT(sequenceRecordsWithMetadata != nullptr, "Input pointer for sequence records with metadata cannot be null");
// 	BIOLSHASHER_ASSERT(!sequenceRecordsWithMetadata->Records.empty(), "Sequence records with metadata cannot be empty");
// 	BIOLSHASHER_ASSERT(similarity_fn != nullptr, "Similarity metric is required for mutation. Please Add Appropriate Similarity metric in your LSH file.");

// 	// Initialise seed and reserve memory.
//     std::mt19937 rng_mut(sequenceRecordsWithMetadata->DataMutateSeed);
// 	std::mt19937 rng_nuc(sequenceRecordsWithMetadata->DataMutateSeed + 7);

// 	std::uniform_int_distribution<uint32_t> dist_mut(0, 1000);

// 	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {

// 		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
		
// 		double snp_rate = record.foundationalParameter;
// 		record.snpRate = snp_rate;
// 		record.delRate = 0.0;	// No indels in this mutator
// 		record.stayRate = 1 - snp_rate;	// No indels in this mutator
// 		record.insmean = 0.0;	// No indels in this mutator
// 		record.insRate = 0.0;	// No indels in this mutator

// 		// Initialising the mutated sequences with original sequences
// 		record.MutatedLength = record.OriginalLength;
// 		record.SeqASCIIMut.resize(record.MutatedLength, 'A');
// 		record.SeqASCIIMut = record.SeqASCIIOrg;	// Copy original ASCII sequence

// 		// NOTE: Use a while loop or track position carefully since length changes
// 		for (int pos = static_cast<int>(record.OriginalLength) - 1; pos >= 0; --pos) {
// 			uint32_t rand_val = dist_mut(rng_mut);
// 			double sample = static_cast<double>(rand_val) / 1000.0;	// snp or indel
// 			if (sample < record.snpRate) {
//                 simulateSNP(record, pos, rng_nuc);
//             }
//         }
// 	}

// 	for(uint32_t rec_idx = 0; rec_idx < sequenceRecordsWithMetadata->KeyCount; rec_idx++) {
// 		auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
// 		record.similarity = similarity_fn(record.SeqASCIIOrg, record.SeqASCIIMut, record.OriginalLength, record.MutatedLength);
// 	}
// 	std::cout << "DataMutation: Applied SNP mutations only." << std::endl;

// 	// Debugging output
//     for(uint32_t rec_idx = 0; rec_idx < std::min(sequenceRecordsWithMetadata->KeyCount, 10u); rec_idx++) {
//         auto& record = sequenceRecordsWithMetadata->Records[rec_idx];
//         if (rec_idx < 10) {
//             std::cout << "\n=== Record " << rec_idx << " ===" << std::endl;
//             std::cout << "Original Length: " << record.OriginalLength << std::endl;
//             std::cout << "Mutated Length:  " << record.MutatedLength << std::endl;
//             std::cout << "Original: " << record.SeqASCIIOrg << std::endl;
//             std::cout << "Mutated:  " << record.SeqASCIIMut << std::endl;
// 			std::cout << "Foundational Parameter: " << record.foundationalParameter << std::endl;
// 			std::cout << "Similarity: " << record.similarity << std::endl;
//         }
//     }	
// }

