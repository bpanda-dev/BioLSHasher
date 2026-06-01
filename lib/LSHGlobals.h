#pragma once

#include "HashInfo.h"
#include <vector>
#include <cstdint>


#define MUTATION_MODEL_SIMPLE_SNP_ONLY 0
#define MUTATION_MODEL_GEOMETRIC_MUTATOR 1



#define MUTATION_EXPRESSION_BALANCED 0
#define MUTATION_EXPRESSION_SUB_ONLY 1
#define MUTATION_EXPRESSION_DEL_LITE 2
#define MUTATION_EXPRESSION_INS_LITE 3
#define MUTATION_EXPRESSION_SUB_LITE 4


extern bool g_areBasesDrawnFromUniformDistribution;

extern uint32_t g_SequenceLength;


extern uint32_t g_verySlowNAggCases;
extern uint32_t g_verySlowNSeq;
extern uint32_t g_verySlowNHashes;

extern uint32_t g_SlowNAggCases;
extern uint32_t g_SlowNSeq;
extern uint32_t g_SlowNHashes;

extern uint32_t g_NAggCases;
extern uint32_t g_NSeq;
extern uint32_t g_NHashes;



extern uint32_t g_start_B;
extern uint32_t g_start_R;
extern uint32_t g_MAX_B;
extern uint32_t g_MAX_R;

extern double g_simThresholdForCollTest;

extern uint32_t g_Nseq_in_Database;
extern uint32_t g_numQueriesForApproxNNTest;

extern uint32_t g_avgRunsForApproxNN;
extern uint32_t g_ANN_start_B;
extern uint32_t g_ANN_start_R;
extern uint32_t g_ANN_MAX_B;
extern uint32_t g_ANN_MAX_R;

extern double g_simThresholdForApproxNNTest;

extern std::vector<double> g_cValuesApproxNNTest;


// #define MUTATION_EXPRESSION_DEL_EQUAL_SUB 0
// #define MUTATION_EXPRESSION_DEL_HALF_SUB 1
// #define MUTATION_EXPRESSION_DEL_DOUBLE_SUB 2
// #define MUTATION_EXPRESSION_DEL_ZERO 3

// Global variables for runtime communication
// extern uint32_t g_TokenLength;
// extern uint32_t g_window_length;

// extern uint32_t g_NumSignatures;
// extern bool g_IsTestActive;
extern const uint32_t g_GoldenRatio;	

extern const uint32_t g_bincount_full;


extern  uint32_t g_mutation_model; // 0: Simple SNP only, 1: Geometric Mutator
extern  uint32_t g_mutation_expression_type;



// extern const uint32_t g_subseqHash1_subseq_len; // Default subsequence length for SubseqHash1 : this is k
// extern const uint32_t g_subseqHash1_d;           // Default 'p' value for SubseqHash1	: this is d

// Setter functions (called by LSHCollision test)
// void SetTokenLength(uint32_t length);
// void SetWindowLength(uint32_t length);

// void SetNumSignatures(uint32_t num);
void SetIsTestActive(bool active);

// Getter functions (called by hash functions)
// uint32_t GetTokenLength();
// uint32_t GetWindowLength();
// uint32_t GetNumSignatures();
bool IsTestActive();




void mutation_expression(double g_mean, uint32_t expression_type, double *P_sub_out, double *P_del_out);

bool is_valid_mutation_parameters(double P_sub, double P_del);

// uint32_t setDistanceClassForHashInfo(const uint64_t hash_flags);


const char* getMutationModelName(uint32_t model);

const char* getMutationExpressionTypeName(uint32_t type);