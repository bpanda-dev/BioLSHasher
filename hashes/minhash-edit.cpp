/*
 * MINHASH
 * Copyright (C) 2025 IISc
 *
 * A brief intro about minhash:
 * MinHash is a locality-sensitive hashing (LSH) technique used to estimate the similarity between
 * two sets. It is particularly useful for large datasets where exact similarity computation can be
 * computationally expensive. 
 * The core idea of MinHash is to create a compact representation (signature)
 * of each set, such that the similarity between the sets can be approximated by comparing their signatures.
 * Here we use the base version of minhash where we take the minimum hash value of the k-mers (tokens) in the sequence as the signature.
 *
 * 
 */

#include "assertMsg.h"
#include "specifics.h"
#include "Hashlib.h"
#include <random>
#include <vector>
#include <cstdint>

#include <vector>
#include <unordered_map>

#include "LSHGlobals.h"


#define TOKEN_LENGTH 32





static void ExtractTokensMinhash(const uint8_t *src, size_t startBase, size_t baseLen, uint8_t *dst){
      size_t startByte = startBase; // Since each base is represented by 1 byte in this context
      for(size_t i = 0; i < baseLen; i++){
            dst[i] = src[startByte + i];
      }
}

//-----------------------------------------------------------------------------
#define _mersenne_prime uint64_t((uint64_t(1) << 61) - 1) // A large prime number: 2^61 - 1 : Mersenne prime

#define _max_hash_32 uint64_t(uint64_t(uint64_t(1) << 32) - 1)  // 2^32 - 1 : Maximum hash value for unsigned 32-bit hash functions
#define _hash_range_32 uint64_t(1 << 32)  // 2^32     : Total Number of buckets. 


static double EditSimilarity(const std::string& seq1, const std::string& seq2, const uint32_t in1_len, const uint32_t in2_len) {

	//write your code here
	const int n = in1_len;
	const int m = in2_len;

	std::vector<std::vector<int>> D(n + 1, std::vector<int>(m + 1, 0));

	// Initialize base cases
	for (int i = 0; i <= n; i++) {
		D[i][0] = static_cast<int>(i);
	}
	for (int j = 0; j <= m; j++) {
		D[0][j] = static_cast<int>(j);
	}

	// Fill DP table
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= m; j++) {
			int insert_cost = D[i][j - 1] + 1;
			int delete_cost = D[i - 1][j] + 1;
			int match_cost  = D[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? 0 : 1);

			D[i][j] = std::min({insert_cost, delete_cost, match_cost});
		}
	}

	int edit_dist = D[n][m];
	double max_len = static_cast<double>(std::max(n, m));

	return 1.0 - (static_cast<double>(edit_dist) / max_len);
}


// NOTE: Only valid for 32 and 64 bit hash output.  
template<typename hashT>
static bool check_equality(void* inp1, void* inp2){
      if(inp1 == nullptr || inp2 == nullptr){
            return false;
      }
      return (*static_cast<hashT*>(inp1) == *static_cast<hashT*>(inp2));
}


static void minHash32(const void* in , const size_t len, const seed_t seed, void* out) {
      
      //We will use the inbuilt FNV1a hash function
      const HashInfo * Fnv1aInfo = findHash("FNV-1a-32");  //Accessing the 32-bit FNV1a hash function that is available in biolshasher.
      if (!Fnv1aInfo) {
            // Fallback to another hash or throw error. For now we will throw an error.
            fprintf(stderr, "FNV1a hash function not found!\n");
            out = nullptr;
            return;
      }
      seed_t SeedFNV = g_GoldenRatio;      //A constant seed for FNV1a hash function which will be used by all the hashes in this hash family.

      const uint8_t* data = static_cast<const uint8_t*>(in);      // input is a neucleotide sequence in char format (ASCII) and length is in bytes.

      if(in==nullptr || len==0){
            PUT_U32(0, (uint8_t *)out, 0);
            return;
      }

      uint32_t BasesInToken = TOKEN_LENGTH;
      uint32_t BytesInToken = BasesInToken; 
      
      BIOLSHASHER_ASSERT(len >= BasesInToken ,"Input length is less than token length"); // Since each base is represented by 1 byte in this context, the number of bytes in the token is equal to the number of bases in the token.

      std::vector<uint8_t> token(BytesInToken, 0); // Token storage.

      //Generate random number,
      // Generate two random numbers 'a' and 'b' for the permutation hash function.
      // We will use a simple linear transformation as the permutation function: h(x) = (a*x + b) mod p, 
      // where p is a large prime number. note: a should be non-zero and both a and b should be less than p.
      // We will use the Mersenne prime 2^61 - 1 the large prime number. 

      std::mt19937 rng(seed);
      
      uint32_t a = rng();
      uint32_t b = rng();

      // I didnt use assert here, because I dont want to stop the function from running. 
      if(a==b){
            printf("Please Note! the parameter values are equal!");
      }
      if(a==0){

            printf("Please Note! the parameter 'a' is zero!");
      }

      uint32_t MinHashValue = UINT32_MAX; // Initialize min hash value to maximum possible.
      for(int i = 0; i < (int)(len - BasesInToken + 1); i++){
            uint32_t initHash = 0;

            std::fill(token.begin(), token.end(), 0);
            uint32_t TokenStartPos = i;

            // Extract the k-mer (token) from the input data
            ExtractTokensMinhash(data, TokenStartPos, BasesInToken, token.data());  //focus on this => check if this actually correct.

            Fnv1aInfo->hashfn(token.data(), BasesInToken, SeedFNV, &initHash);

            uint32_t TokenHash = 0;
            TokenHash = (((((uint64_t)a * initHash) + b) % _mersenne_prime) & _max_hash_32); // <-- Applying the permutation hash function

            if(TokenHash < MinHashValue){
                  MinHashValue = TokenHash;
            }
      }

      PUT_U32(MinHashValue, (uint8_t*)out, 0);
}
//------------------------------------------------------------

REGISTER_FAMILY(MinHashEdit,
   $.src_url    = "https://github.com/bpanda-dev/BioLSHasher",
   $.src_status = HashFamilyInfo::SRC_ACTIVE
);

REGISTER_HASH(MinHashEdit_32,
   $.desc            = "32-bit MinHash with Edit version",
   $.hash_flags      = FLAG_HASH_LOCALITY_SENSITIVE,
   $.impl_flags      = FLAG_IMPL_VERY_SLOW,
   $.bits            = 32,
   $.hashfn          = minHash32,
   $.parameterNames  = {"k"},
   $.parameterDescriptions = {"Kmer Length"},
   $.parameterValues = {TOKEN_LENGTH},
   $.similarity_name = "Edit",
   $.similarityfn    = EditSimilarity,
   $.check_equality_fn = check_equality<uint32_t>
);

//------------------------------------------------------------
