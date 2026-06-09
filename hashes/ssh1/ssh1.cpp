// This file implements ssh1 and integrates it into the biolshasher3 framework.


// This code is from subseqhash2 src. but is corrected. Contains TB logic too,

#include "specifics.h"
#include "Hashlib.h"
#include "LSHGlobals.h"
#include <iostream>
#include <vector>

#include <random>
#include <map>
#include <cstring>
#include "assertMsg.h"
#include <algorithm>


#define ALPHABETSIZE 4
const char ALPHABET[ALPHABETSIZE] = {'A', 'C', 'G', 'T'};

#define SUBSEQUENCE_LENGTH 19
#define D_PARAM 19

#define MAXK 100
#define MAXD 100


typedef uint64_t kmer;
const double INF = 1e15;

// TODO: This struct gives more flexibility to the output of the hash function,
// but we need to make sure it works well with the smhasher framework. 
// We might need to adjust the way we return the hash value and other information.
// //this needs to go to the global lsh header file, but for now we can keep it here.
struct genomics_seed
{
    int64_t hashval;
    kmer str;		// Remember, the kmer size is maxed out at 64, so it can fit in a uint64_t. If we want to support longer kmers, we might need to change this to a different data structure.
    size_t start, end;
    uint64_t index = 0;	
};


struct DPCell{
    double f_max, f_min;	// values
    int g_max, g_min;	// directions for traceback
};


static bool check_equality_64(void* inp1, void* inp2){
	uint64_t* data1 = static_cast<uint64_t*>(inp1);
	uint64_t* data2 = static_cast<uint64_t*>(inp2);

	if(data1 == nullptr || data2 == nullptr){
		return false;
	}
	return (*data1 == *data2);
}


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


// static std::set<std::string> get_kmers(const std::string& sequence, int k) {
//     std::set<std::string> kmers;
    
//     // Ensure the sequence is long enough for at least one k-mer
//     if (sequence.length() < static_cast<size_t>(k)) {
//         return kmers;
//     }

//     // Slide the window across the sequence
//     for (size_t i = 0; i <= sequence.length() - k; ++i) {
//         kmers.insert(sequence.substr(i, k));
//     }
    
//     return kmers;
// }

// static double JaccardSimilarity(const std::string& seq1, const std::string& seq2, const uint32_t in1_len, const uint32_t in2_len) {

//       uint32_t k = (uint32_t)SUBSEQUENCE_LENGTH; // k-mer length

//       // 1. Tokenize sequences into sets of k-mers
//       std::set<std::string> set1 = get_kmers(seq1, k);
//       std::set<std::string> set2 = get_kmers(seq2, k);

//       if (set1.empty() || set2.empty()) return 0.0;

//       // 2. Find the intersection
//       std::vector<std::string> intersect;
//       std::set_intersection(set1.begin(), set1.end(),
//                               set2.begin(), set2.end(),
//                               std::back_inserter(intersect));

//       // 3. Use the Inclusion-Exclusion principle for the Union size
//       // |A ∪ B| = |A| + |B| - |A ∩ B|
//       size_t intersection_size = intersect.size();
//       size_t union_size = set1.size() + set2.size() - intersection_size;

//       return static_cast<double>(intersection_size) / union_size;
// }


static inline int alphabetIndex(char c)
{
    return 3 & ((c>>2) ^ (c>>1));
}

static inline int dpIndex(int d1, int d2, int d3)
{
    // Assuming dimensions from the original context
    int k = SUBSEQUENCE_LENGTH;
    int d = D_PARAM;
    return d1 * (k+1) * d + d2 * d + d3;
}

char* decode(const kmer enc, const int k, char* str)	// For decoding the strings :debugging
{
    if(str == NULL)
		str = (char*)malloc(sizeof *str *k);

    kmer enc_copy = enc;

    for(int i=k-1; i>=0; i-=1)
    {
		str[i] = ALPHABET[enc_copy & 3];
		enc_copy >>= 2;
    }

    return str;
}

static void SubseqHash64(const void* in, const size_t len, const seed_t rand_seed, void* out) {

	const char* s;
	s = static_cast<const char *>(in);

	int n = len;
    int k = SUBSEQUENCE_LENGTH;
	int d = D_PARAM;
    int dim1 = (n+1) * (k+1) * d;
    // int dim2 = d;

	// printf("dim1: %d\n", dim1);
	// printf("dim2: %d\n", (WINDOW_SIZE+1) * (k+1) * d);
	
	BIOLSHASHER_ASSERT(n >= k, "Sequence length must be >= subsequence length");
	BIOLSHASHER_ASSERT(d < MAXD, "D_PARAM must be less than MAXD");
	BIOLSHASHER_ASSERT(k < MAXK, "SUBSEQUENCE_LENGTH must be less than MAXK");


	//-------------------------------------------------------------------------------------//
	double best_val = 0;


	genomics_seed gen_seed;

	gen_seed.start = 0;
	gen_seed.end = 0;
	gen_seed.index = 0;
	gen_seed.hashval = 0;
	gen_seed.str = 0;

	//-------------------------------------------------------------------------------------//
	double A[MAXK][ALPHABETSIZE][MAXD];
	
	int B1[MAXK][ALPHABETSIZE][MAXD];
	int B2[MAXK][ALPHABETSIZE][MAXD];

	int C[MAXK][ALPHABETSIZE];

	std::vector<int> pos;
	std::vector<int> possign;

    for(int i = 0; i < 4; i++)
		possign.push_back(i);

	for(int i = 0; i < d; i++)
		pos.push_back(i);
        
	if(d < 4)
		for(int i = d; i < 4; i++)
			pos.push_back(i%d);

	std::default_random_engine generator(rand_seed);
	std::uniform_real_distribution<double> distribution((int64_t)1<<30, (int64_t)1<<31);

    for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < 4; j++)
			for(int q = 0; q < d; q++)
				A[i][j][q] = distribution(generator);

		seed_t internal_seed_1 = rand_seed + 7*i + 7; // simple way to change seed for different shuffles
		std::shuffle(pos.begin(), pos.end(), std::default_random_engine(internal_seed_1));

		for(int j = 0; j < 4; j++)
			C[i][j] = pos[j%d];

		for(int q = 0; q < d; q++)
		{	
			seed_t internal_seed_2 = rand_seed + 13*q + 13; // simple way to change seed for different shuffles
			std::shuffle(possign.begin(), possign.end(), std::default_random_engine(internal_seed_2));

			for(int j = 0; j < 4; j++)
			{
				B1[i][j][q] = (possign[j] % 2) ? 1: -1;
				B2[i][j][q] = (possign[j] / 2) ? 1: -1;
			}
		}
	}

    //-------------------------------------------------------------------------------------//

	DPCell* dp;
    int* h;

    dp = (DPCell*) malloc(sizeof *dp * (dim1));
	h = (int*) malloc(sizeof *h * dim1);
    
    // If we are processing the input as a single window, n should be its length.
    n = len;
    dim1 = (n+1) * (k+1) * d;

    int dp_index;
    int del = n - k; // 'del' is the max number of characters we can skip

    memset(h, 0, sizeof(*h) * dim1);

    for(int i = 0; i <= n; i++)
	{
	    dp_index = dpIndex(i, 0, 0); //[i][0][0]

	    dp[dp_index].f_max = 0;
	    dp[dp_index].f_min = 0;
	    
	    dp_index += 1;
	    for(int j = 1; j < d; ++j, ++dp_index)
	    {
			dp[dp_index].f_max = -INF;
			dp[dp_index].f_min = INF;
	    }
	}

    // The main loop over windows is removed. We process the sequence as a single window.
    // Let's assume st = 0.
    int st = 0;

    bool skip = 0;
	for(int i = st; i < st + n; i++)
		if(s[i] == 'N')
		{
			skip = 1;
			break;
		}	

    if(!skip)
    {
		for(int i = 0; i <= n; i++)
		{
			dp_index = dpIndex(i, 0, 0);
		    h[dp_index] = st + 1;
		}
	
		int d1 = C[0][alphabetIndex(s[st])];
		dp_index = dpIndex(1, 1, d1);

		dp[dp_index].f_min = dp[dp_index].f_max = B2[0][alphabetIndex(s[st])][d1] * A[0][alphabetIndex(s[st])][d1];
		dp[dp_index].g_min = dp[dp_index].g_max = 0;
		h[dp_index] = st + 1;

		//printf("%d %d %d %.2lf %.2lf %d %d\n", st, d1, dp_index, dp[dp_index].f_max, dp[dp_index].f_min, dp[dp_index].g_max, dp[dp_index].g_min);
		double v1, v2;

		for(int i = 2; i <= n; i++)
		{
		    int minj = std::max(1, i - del);
		    int maxj = std::min(i, k);

		    for(int j = minj; j <= maxj; j++)
		    {
				int now = alphabetIndex(s[st + i - 1]); //Current char
				int v = C[j-1][now];

				for(int q = 0; q < d; q++)
				{
				    int z = (q + v) % d; // q: previous d, z: current d
				    
				    int idx = dpIndex(i, j, z);
				    
				    int laidx1 = dpIndex(i-1, j, z);// = index_cal(i-1,j,z);
				    int laidx2 = dpIndex(i-1, j-1, q);// = index_cal(i-1,j-1,q);

				    if(h[idx] != st + 1)
				    {
						dp[idx].f_min = INF;
						dp[idx].f_max = -INF;
				    }
				    
				    if(h[laidx1] == st + 1)
				    {
						if(dp[laidx1].f_min < dp[idx].f_min)
						{
						    dp[idx].f_min = dp[laidx1].f_min;
						    dp[idx].g_min = dp[laidx1].g_min;
						}
						if(dp[laidx1].f_max > dp[idx].f_max)
						{
						    dp[idx].f_max = dp[laidx1].f_max;
						    dp[idx].g_max = dp[laidx1].g_max;
						}
						h[idx] = st + 1;
				    }

				    if(h[laidx2] == st + 1)
				    {
						if(B1[j-1][now][z] == -1)
						{
						    v1 = -dp[laidx2].f_min;
						    v2 = -dp[laidx2].f_max;

						    if(B2[j-1][now][z] == -1)
						    {
								v1 -= A[j-1][now][z];
								v2 -= A[j-1][now][z];
						    }
						    else
						    {
								v1 += A[j-1][now][z];
								v2 += A[j-1][now][z];
						    }

						    if(v1 > dp[idx].f_max)
						    {
								dp[idx].f_max = v1;
								dp[idx].g_max = i-1;
						    }
						    if(v2 < dp[idx].f_min)
						    {
								dp[idx].f_min = v2;
								dp[idx].g_min = i-1;
						    }
						}

						else
						{
						    v1 = dp[laidx2].f_min;
						    v2 = dp[laidx2].f_max;

						    if(B2[j-1][now][z] == -1)
						    {
								v1 -= A[j-1][now][z];
								v2 -= A[j-1][now][z];
						    }
						    else
						    {
								v1 += A[j-1][now][z];
								v2 += A[j-1][now][z];
						    }

						    if(v1 < dp[idx].f_min)
						    {
								dp[idx].f_min = v1;
								dp[idx].g_min = i-1;
						    }
						    if(v2 > dp[idx].f_max)
						    {
								dp[idx].f_max = v2;
								dp[idx].g_max = i-1;
						    }
						}
						
						h[idx] = st + 1;
				    }			

				    //printf("%d %d %d %d %d %.2lf %.2lf %d %d\n", i,j,z,idx,h[idx],dp[idx].f_max, dp[idx].f_min, dp[idx].g_max, dp[idx].g_min);
				}
			}
		}
        
        genomics_seed tmp;
		int mod;
		dp_index = dpIndex(n, k, 0);

		int i = 0;
		for(i = 0; i < d; i++){
			if(h[dp_index + i] == st + 1)
			{
				if(fabs(dp[dp_index + i].f_min) > dp[dp_index + i].f_max)
					tmp.hashval = uint64_t(fabs(dp[dp_index + i].f_min) * 32768);
				else
					tmp.hashval = uint64_t(dp[dp_index + i].f_max * 32768);

				mod = i;
                break;
			}
		}
		best_val = tmp.hashval;


		
		// dp_index = dpIndex(n, k, 0);
		// int i = 0;
		// for(i = 0; i < d; i++){
		// 	if(h[dp_index + i] == st + 1)
		// 	{
		// 		best_val = std::max(best_val, fabs(dp[dp_index + i].f_min));
		// 		best_val = std::max(best_val, dp[dp_index + i].f_max);
		// 		mod = i;
		// 		break;
		// 	}
		// }
			
		// If no valid hash was found, we can't proceed. Added by bpanda-dev
		if (i == d) {
			printf("No valid hash found for the given input sequence.\n");
			return; // TODO: Need to handle it with much better way. 
		}

		kmer hashval = 0;
		std::vector<size_t> index;

		int x = n;
		int y = k;
		size_t nextpos;
		int nextval;

		bool z = 1;
		if(tmp.hashval < 0)
			z = 0;

		while(y > 0)
		{
			dp_index = dpIndex(x, y, mod);
			if(z == 1)
			{
				nextpos = st + dp[dp_index].g_max;
				nextval = alphabetIndex(s[nextpos]);

				index.push_back(nextpos);
				hashval = (hashval<<2) | nextval;

				if(B1[y-1][nextval][mod] == -1)
					z = 0;

				x = dp[dp_index].g_max;
				y--;
				mod = (mod + d - C[y][nextval]) % d;
			}
			else
			{
				nextpos = st + dp[dp_index].g_min;
				nextval = alphabetIndex(s[nextpos]);

				index.push_back(nextpos);
				hashval = (hashval<<2) | nextval;

				if(B1[y-1][nextval][mod] == -1)
					z = 1;

				x = dp[dp_index].g_min;
				y--;
				mod = (mod + d - C[y][nextval]) % d;
			}
		}

		tmp.start = index.back();
		tmp.end = index[0];

	    for(size_t a: index){
			tmp.index |= ((uint64_t)1)<<(a - tmp.start);
		}
	    	
		// The logic for comparing with a vector of seeds is removed.
        
		tmp.str = hashval;
        
		//disabled.
		gen_seed = tmp;	// We assign the result to the output parameter 'out'. Note that we need to make sure the caller of this function knows how to interpret the output as a 'genomics_seed' struct.
		
		// To return seed, we can return the gen_seed datastructure.
	}
    
	free(dp);
    free(h);

	//-------------------------------------------------------------------------------------//
	// Copy double bits directly to output
	uint64_t hash;
	memcpy(&hash, &best_val, sizeof(double));

	PUT_U64(hash, (uint8_t*)out, 0);

}

REGISTER_FAMILY(SubseqHash,
   $.src_url    = "https://github.com/Shao-Group/subseqhash",
   $.src_status = HashFamilyInfo::SRC_STABLEISH
);

REGISTER_HASH(SubseqHash_64,
   $.desc            = "Subsequence hash with edit distance tolerance (64-bit)",
   $.hash_flags      = FLAG_HASH_LOCALITY_SENSITIVE,
   $.impl_flags      = FLAG_IMPL_VERY_SLOW,
   $.bits            = 64,
   $.hashfn			 = SubseqHash64,
   $.parameterNames  = {"k", "d"},
   $.parameterDescriptions  = {"Subsequence Length 'k'", "Parameter 'd'"},
   $.parameterValues = {SUBSEQUENCE_LENGTH, D_PARAM},
   $.similarityfn   = EditSimilarity,
   $.similarity_name = "Edit",
   $.check_equality_fn = check_equality_64
);
