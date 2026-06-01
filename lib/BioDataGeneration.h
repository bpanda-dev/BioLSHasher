
#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <cassert>
#include <unordered_set>
#include <fstream>
#include <iomanip>

#include "HashInfo.h"
#include "LSHGlobals.h"

//Not used
#define two_bit_per_base 2 // 2 bits per base for A,C,G,T

static constexpr char bases[4] = {'A', 'C', 'T', 'G'};
static constexpr size_t NUM_BASES = std::size(bases);

typedef uint64_t  seed_t; // Seed type for random number generators

struct UnionBitVectorsStruct {
    std::vector<char> vec_a;
    std::vector<char> vec_b;
    std::vector<std::string> universe;
};

// Store k-mers with their positions for nearest neighbour lookups
struct KmerEntry {
	std::string kmer;
	uint32_t position;  // Start position in the reference sequence
};

struct NearestKmerEntry {
	std::string kmer;
	uint32_t position;     // Position in the reference sequence
	double similarity;     // Similarity to the query
};

// double ComputeHammingSimilarity(const std::string& seq1, const std::string& seq2, const uint32_t length);
// double ComputeJaccardSimilarity(const std::string& seq1, const std::string& seq2, int k);
// double ComputeCosineSimilarity(const std::string& seq1, const std::string& seq2, int k);
// double ComputeAngularSimilarity(const std::string& seq1, const std::string& seq2, int k);
// double ComputeEditSimilarity(const std::string& seq1, const std::string& seq2);


UnionBitVectorsStruct CreateUnionBitVectors(const std::string& seq1, const std::string& seq2, int k);

struct SequenceRecordUnit{
	uint32_t OriginalLength;
	uint32_t MutatedLength;

	std::string SeqASCIIOrg;
	std::string SeqASCIIMut;

	std::vector<uint8_t> SeqTwoBitOrg;	//Not used
	std::vector<uint8_t> SeqTwoBitMut;	//Not used

	double similarity;

	uint32_t LocalAlignmentScore;
	uint32_t GlobalAlignmentScore;

	std::string LocalCIGARString;
	std::string GlobalCIGARString;

	// Mutation Parameters
    double foundationalParameter;
	double snpRate;
	double delRate;
	double stayRate;
	double insmean;
	double insRate; //For 1 length.

	uint32_t count_ins;
	uint32_t count_del;
	uint32_t count_sub;
	uint32_t count_stay;
	uint32_t count_ins_length;	// Total length of insertions (sum of lengths of all insertions)
};

struct sim_bins_struct{
	std::vector<uint32_t> bin_fill_count = std::vector<uint32_t>(g_bincount_full, 0); // Count of entries in each bin
	std::vector<std::vector<double>> bin_error_parameters = std::vector<std::vector<double>>(100, std::vector<double>(g_bincount_full));
	std::vector<double> bin_error_parameters_mean = std::vector<double>(100, 0);
	std::vector<double> bin_error_parameters_stddev = std::vector<double>(100, 0);
};

struct SequenceRecordsWithMetadataStruct{
	uint32_t KeyCount;                    // Number of sequence pairs
    uint32_t OriginalSequenceLength;      // Original length (before mutations)

	seed_t DataGenSeed;                          // For reproducibility, we have base seed
	seed_t DataMutateSeed;                        // For reproducibility, we have mutation seed
	seed_t HashSeed;		// Seed for hash family generation
	
	uint32_t HashCount;		// Number of hash functions to be used for LSH

	// Content of sequences
	bool areBasesDrawnFromUniformDist;
	double A_percentage;	// Percentage of A bases in the original sequence
	double C_percentage;	// Percentage of C bases in the original sequence
	double G_percentage;	// Percentage of G bases in the original sequence
	double T_percentage;	// Percentage of T bases in the original sequence

	float binstart = 0.0f;	// Bin start
	float binend = 1.0f;		// Bin end
	float binsize = 0.01f;	// Bin size for distance metrics
	uint32_t bincount = 100;	// Number of bins	// This is similarity bins

	// Array of sequence records (AoS)
    std::vector<SequenceRecordUnit> Records;
	

	SequenceRecordsWithMetadataStruct() : KeyCount(0), OriginalSequenceLength(0), DataGenSeed(0), DataMutateSeed(0),
		A_percentage(0.25), C_percentage(0.25), G_percentage(0.25), T_percentage(0.25) {}
};


class SeedGenerator {
private:
    std::mt19937_64 rng;  // 64-bit Mersenne Twister
    std::uniform_int_distribution<seed_t> dist;

public:
    explicit SeedGenerator(seed_t baseSeed) : rng(baseSeed), dist(0, UINT64_MAX) {}

    seed_t nextSeed() {
        return dist(rng);
    }
};

class Randbin {
private:
    std::mt19937 rng;  // Mersenne Twister 19937 generator
    std::uniform_int_distribution<uint32_t> dist;

public:
    // explicit Randbin(seed_t seed) : rng(seed), dist(0, UINT32_MAX) {}

	// To make the 64 bit seed visible to 32 bit mt19937, we can XOR-fold the upper and lower halves of the seed.
	explicit Randbin(seed_t seed) : rng(static_cast<uint32_t>(seed ^ (seed >> 32))),  // XOR-fold upper/lower halves
      	dist(0, UINT32_MAX) {}

    // Generate a random number in the range [0, range)
    uint32_t rand_range(uint32_t range) {
        return dist(rng) % range;
    }
	
	// Generate a uniform random int in [low, high]
    int rand_custom_range(int low, int high) {
        std::uniform_int_distribution<int> int_dist(low, high);
        return int_dist(rng);
    }
};


class SequenceDataGenerator {
	public:
		// Constructor to initialise and generate random sequences.
		SequenceDataGenerator(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata);
		
		
		static void decodeSequencesToASCII(const std::vector<uint8_t>& seqTwoBit, std::string& seqASCII, uint32_t sequenceLength);
};

class SequenceDataMutatorSubstitutionOnly{
	private:
		// uint8_t GetBaseAtPosition(const std::vector<uint8_t>& seq, uint32_t pos);	// Helper: Get base at position
		// void SetBaseAtPosition(std::vector<uint8_t>& seq, uint32_t pos, uint8_t base);	// Helper: Set base at position
		// bool simulateSNP(SequenceRecordUnit &record, const uint32_t pos, std::mt19937 &randGen);
	public:
		SequenceDataMutatorSubstitutionOnly(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata, std::vector<double> *rand_error_param, SimilarityFn similarity_fn);
		SequenceDataMutatorSubstitutionOnly(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata, SimilarityFn similarity_fn);
};

class SequenceDataMutatorGeometric{
	private:
		// uint8_t GetBaseAtPosition(const std::vector<uint8_t>& seq, uint32_t pos);	// Helper: Get base at position
		// void SetBaseAtPosition(std::vector<uint8_t>& seq, uint32_t pos, uint8_t base);	// Helper: Set base at position
		// bool simulateSNP(SequenceRecordUnit &record, const uint32_t pos, std::mt19937 &randGen);
	public:
		SequenceDataMutatorGeometric(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata, std::vector<double> *rand_error_param, SimilarityFn similarity_fn);
		SequenceDataMutatorGeometric(SequenceRecordsWithMetadataStruct *sequenceRecordsWithMetadata, SimilarityFn similarity_fn);
};
