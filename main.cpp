#include <cstdlib>
#include <iostream>
#include "version.h"
#include "toml.hpp"
#include "Hashlib.h"
#include "TestGlobals.h"
#include "smhasher3Random.h"
#include "LSHGlobals.h"

#include "LSHCollision.h"
#include "ApproxNearestNeighbour.h"

static bool g_testAll;
static bool g_testLSHCollision;
static bool g_testLSHApproxNearestNeighbour;

struct TestOpts {
    bool &       var;
    bool         defaultvalue;
    const char * name;
};

static TestOpts g_testopts[] = {
    { g_testLSHCollision,   false,   "LSHCollision" },
    { g_testLSHApproxNearestNeighbour,   false,    "LSHApproxNearestNeighbour" },
};

struct BioLSHasherConfig {
    uint32_t Mutation_Model;
    uint32_t Mutation_Expression_Type;

    uint32_t Sequence_Length;
    
    uint32_t Coll_Num_Hashes_To_And_Sweep_Start;
    uint32_t Coll_Num_Hashes_To_And_Sweep_End;
    uint32_t Coll_Num_Hashes_To_Or_Sweep_Start;
    uint32_t Coll_Num_Hashes_To_Or_Sweep_End;

    uint32_t Total_Sequences_In_Database_To_Search_From;
    uint32_t Total_Query_Sequences_To_Search;
    uint32_t Number_of_Runs_to_get_Average;

    uint32_t SimSearch_Num_Hashes_To_And_Sweep_Start;
    uint32_t SimSearch_Num_Hashes_To_And_Sweep_End;
    uint32_t SimSearch_Num_Hashes_To_Or_Sweep_Start;
    uint32_t SimSearch_Num_Hashes_To_Or_Sweep_End;

    double Similarity_Threshold;

    bool are_Bases_Drawn_From_Uniform_Distribution;
    std::vector<double> Categorical_Distribution_Probabilities;

    uint32_t very_Slow_NAggCases;
    uint32_t very_Slow_NSeq;
    uint32_t very_Slow_NHashes;

    uint32_t Slow_NAggCases;
    uint32_t Slow_NSeq;
    uint32_t Slow_NHashes;

    uint32_t NAggCases;
    uint32_t NSeq;
    uint32_t NHashes;

};

void load_config(const std::string& filepath) {
    try {
        // Parse the file into a TOML table
        toml::table config = toml::parse_file(filepath);


        BioLSHasherConfig configuration;
        std::vector<double> extracted_categorical_probs;
        
        if (toml::array* arr = config["Common_Configuration"]["Categorical_Distribution_Probabilities"].as_array()) {
            extracted_categorical_probs.reserve(arr->size());
            for (toml::node& node : *arr) {
                // value_or(0.0) forces the compiler to deduce double, avoiding -Wnarrowing
                extracted_categorical_probs.push_back(node.value_or(0.25)); 
            }
            
            // checking if sum = 1.0
            double sum = std::accumulate(extracted_categorical_probs.begin(), extracted_categorical_probs.end(), 0.0);
            if (std::abs(sum - 1.0) > 1e-6) {
                throw std::domain_error("Probability distribution do not sum to 1.0");
            }
        } else {
            throw std::runtime_error("Configuration error: 'Probabilities' array is missing or malformed.");
        }
        
        configuration.Mutation_Model = config["Mutation_Model_Configuration"]["Mutation_Model"].value_or(1u);
        configuration.Mutation_Expression_Type = config["Mutation_Model_Configuration"]["Mutation_Expression_Type"].value_or(0u);
        
        configuration.Coll_Num_Hashes_To_And_Sweep_Start = config["Collision_Test_Configuration"]["Coll_Num_Hashes_To_And_Sweep_Start"].value_or(1u);
        configuration.Coll_Num_Hashes_To_And_Sweep_End = config["Collision_Test_Configuration"]["Coll_Num_Hashes_To_And_Sweep_End"].value_or(1u);
        configuration.Coll_Num_Hashes_To_Or_Sweep_Start = config["Collision_Test_Configuration"]["Coll_Num_Hashes_To_Or_Sweep_Start"].value_or(1u);
        configuration.Coll_Num_Hashes_To_Or_Sweep_End = config["Collision_Test_Configuration"]["Coll_Num_Hashes_To_Or_Sweep_End"].value_or(1u);
        
        configuration.very_Slow_NSeq = config["Collision_Test_Configuration"]["very_Slow_NSeq"].value_or(1000u);
        configuration.very_Slow_NHashes = config["Collision_Test_Configuration"]["very_Slow_NHashes"].value_or(1000u);
        configuration.Slow_NSeq = config["Collision_Test_Configuration"]["Slow_NSeq"].value_or(3000u);
        configuration.Slow_NHashes = config["Collision_Test_Configuration"]["Slow_NHashes"].value_or(1500u);
        configuration.NSeq = config["Collision_Test_Configuration"]["NSeq"].value_or(5000u);
        configuration.NHashes = config["Collision_Test_Configuration"]["NHashes"].value_or(2000u);
        
        configuration.Total_Sequences_In_Database_To_Search_From = config["Similarity_Search_Test_Configuration"]["Total_Sequences_In_Database_To_Search_From"].value_or(50000u);
        configuration.Total_Query_Sequences_To_Search = config["Similarity_Search_Test_Configuration"]["Total_Query_Sequences_To_Search"].value_or(2000u);
        configuration.Number_of_Runs_to_get_Average = config["Similarity_Search_Test_Configuration"]["Number_of_Runs_to_get_Average"].value_or(5u);
        configuration.SimSearch_Num_Hashes_To_And_Sweep_Start = config["Similarity_Search_Test_Configuration"]["SimSearch_Num_Hashes_To_And_Sweep_Start"].value_or(1u);
        configuration.SimSearch_Num_Hashes_To_And_Sweep_End = config["Similarity_Search_Test_Configuration"]["SimSearch_Num_Hashes_To_And_Sweep_End"].value_or(1u);
        configuration.SimSearch_Num_Hashes_To_Or_Sweep_Start = config["Similarity_Search_Test_Configuration"]["SimSearch_Num_Hashes_To_Or_Sweep_Start"].value_or(1u);
        configuration.SimSearch_Num_Hashes_To_Or_Sweep_End = config["Similarity_Search_Test_Configuration"]["SimSearch_Num_Hashes_To_Or_Sweep_End"].value_or(1u);
        configuration.Similarity_Threshold = config["Similarity_Search_Test_Configuration"]["Similarity_Threshold"].value_or(0.95);
        configuration.are_Bases_Drawn_From_Uniform_Distribution = config["Common_Configuration"]["are_Bases_Drawn_From_Uniform_Distribution"].value_or(true);
        configuration.Categorical_Distribution_Probabilities = extracted_categorical_probs;

        configuration.Sequence_Length = config["Common_Configuration"]["Sequence_Length"].value_or(50u);

        configuration.very_Slow_NAggCases = config["Common_Configuration"]["very_Slow_NAggCases"].value_or(200000u);
        configuration.Slow_NAggCases = config["Common_Configuration"]["Slow_NAggCases"].value_or(200000u);
        configuration.NAggCases = config["Common_Configuration"]["NAggCases"].value_or(200000u);
        
        g_mutation_model = configuration.Mutation_Model;
        g_mutation_expression_type = configuration.Mutation_Expression_Type;
        g_SequenceLength = configuration.Sequence_Length;

        g_start_B = configuration.Coll_Num_Hashes_To_And_Sweep_Start;
        g_start_R = configuration.Coll_Num_Hashes_To_Or_Sweep_Start;
        g_MAX_B = configuration.Coll_Num_Hashes_To_And_Sweep_End;
        g_MAX_R = configuration.Coll_Num_Hashes_To_Or_Sweep_End;

        g_Nseq_in_Database = configuration.Total_Sequences_In_Database_To_Search_From;
        g_numQueriesForApproxNNTest = configuration.Total_Query_Sequences_To_Search;
        g_avgRunsForApproxNN = configuration.Number_of_Runs_to_get_Average;
        g_ANN_start_B = configuration.SimSearch_Num_Hashes_To_And_Sweep_Start;
        g_ANN_MAX_B = configuration.SimSearch_Num_Hashes_To_And_Sweep_End;
        g_ANN_start_R = configuration.SimSearch_Num_Hashes_To_Or_Sweep_Start;
        g_ANN_MAX_R = configuration.SimSearch_Num_Hashes_To_Or_Sweep_End;
        g_simThresholdForApproxNNTest = configuration.Similarity_Threshold;
        g_areBasesDrawnFromUniformDistribution = configuration.are_Bases_Drawn_From_Uniform_Distribution;
        // TODO:aDD THE CATEGORICAL DISTRIBUTION PROBABILITIES TO THE GLOBAL VARIABLE

        g_verySlowNAggCases = configuration.very_Slow_NAggCases;
        g_verySlowNSeq = configuration.very_Slow_NSeq;
        g_verySlowNHashes = configuration.very_Slow_NHashes;

        g_SlowNAggCases = configuration.Slow_NAggCases;
        g_SlowNSeq = configuration.Slow_NSeq;
        g_SlowNHashes = configuration.Slow_NHashes;

        g_NAggCases = configuration.NAggCases;
        g_NSeq = configuration.NSeq;
        g_NHashes = configuration.NHashes;

    } catch (const toml::parse_error& err) {
        std::cerr << "Configuration parsing failed at " 
                  << err.source().begin << ":\n" << err.description() << "\n";
        throw; 
    }
}

static void set_default_tests( bool enable ) {
    for (size_t i = 0; i < sizeof(g_testopts) / sizeof(TestOpts); i++) {
        if (enable) {
            g_testopts[i].var = g_testopts[i].defaultvalue;
        } else if (g_testopts[i].defaultvalue) {
            g_testopts[i].var = false;
        }
    }
}

static void parse_tests( const char * str, bool enable_tests ) {
    while (*str != '\0') {
        size_t       len;
        const char * p = strchr(str, ',');
        if (p == NULL) {
            len = strlen(str);
        } else {
            len = p - str;
        }

        struct TestOpts * found = NULL;
        bool foundmultiple      = false;
        for (size_t i = 0; i < sizeof(g_testopts) / sizeof(TestOpts); i++) {
            const char * testname = g_testopts[i].name;
            // Allow the user to specify test names by case-agnostic
            // unique prefix.
            if (strncasecmp(str, testname, len) == 0) {
                if (found != NULL) {
                    foundmultiple = true;
                }
                found = &g_testopts[i];
                if (testname[len] == '\0') {
                    // Exact match found, don't bother looking further, and
                    // don't error out.
                    foundmultiple = false;
                    break;
                }
            }
        }
        if (foundmultiple) {
            printf("Ambiguous test name: --%stest=%*s\n", enable_tests ? "" : "no", (int)len, str);
            goto error;
        }
        if (found == NULL) {
            printf("Invalid option: --%stest=%*s\n", enable_tests ? "" : "no", (int)len, str);
            goto error;
        }

        // printf("%sabling test %s\n", enable_tests ? "en" : "dis", testname);
        found->var = enable_tests;

        // If "All" tests are being enabled or disabled, then adjust the individual
        // test variables as instructed. Otherwise, if a material "All" test (not
        // just a speed-testing test) is being specifically disabled, then don't
        // consider "All" tests as being run.
        if (&found->var == &g_testAll) {
            set_default_tests(enable_tests);
        } else if (!enable_tests && found->defaultvalue) {
            g_testAll = false;
        }

        if (p == NULL) {
            break;
        }
        str += len + 1;
    }

    return;

  error:
    printf("Valid tests: --test=%s", g_testopts[0].name);
    for (size_t i = 1; i < sizeof(g_testopts) / sizeof(TestOpts); i++) {
        printf(",%s", g_testopts[i].name);
    }
    printf(" \n");
    exit(1);
}

void usage() {
    printf("----------------------------------------------------------------------");
    printf("\n Usage: BioLSHasher [--test=<testname>[,...]] [--verbose] [--ncpu=N]\n"
           "                    [<hashname>]\n"
           "\n"
           "        BioLSHasher [--list]|[--listnames]|[--tests]|[--version]|[--configinit]\n"
           "\n"
           " Hashnames can be supplied using any case letters.\n");
    
    printf("----------------------------------------------------------------------");
}
//-----------------------------------------------------------------------------


static void examplehash_initialisation() {
    static bool initialized = false;
    if (!initialized) {   
        //---- Variables used across tests which can be defined by user ---//
        g_verySlowNAggCases = 50000;
        g_verySlowNSeq = 1000;
        g_verySlowNHashes = 700;

        g_SlowNAggCases = 50000;
        g_SlowNSeq = 1000;
        g_SlowNHashes = 700;

        g_NAggCases = 50000;
        g_NSeq = 1000;
        g_NHashes = 700;

        g_SequenceLength = 50;

        g_start_B = 1;
        g_start_R = 1;
        g_MAX_B = 2;
        g_MAX_R = 3;

        //---------------------------------------------------------------------------------------
        
        g_Nseq_in_Database = 50000; // Number of sequences in the reference database for the Approx Nearest Neighbour test. Adjust as needed.
        g_numQueriesForApproxNNTest = 1000; // Number of query sequences to generate for the Approx Nearest Neighbour test. Adjust as needed.

        g_avgRunsForApproxNN = 2;

        g_cValuesApproxNNTest = {1}; //{0.5, 0.6, 0.7, 0.8, 0.9, 0.95}; // c-ANN approximation factors to sweep. Each c < 1 defines the boundary as c * target_sim_low.

        g_ANN_start_B = 1; // Starting value of b (hashes per table) for the Approx Nearest Neighbour test. Adjust as needed.
        g_ANN_start_R = 1; // Starting value of r (number of tables) for the Approx Nearest Neighbour test. Adjust as needed.
        g_ANN_MAX_B = 2; // Maximum value of b (hashes per table) to test in the Approx Nearest Neighbour test. Adjust as needed.
        g_ANN_MAX_R = 3; // Maximum value of r (number of tables) to test in the Approx Nearest Neighbour test. Adjust as needed.

        g_simThresholdForApproxNNTest = 0.95; // Similarity threshold for Approx Nearest Neighbour test. Adjust as needed.

        //---------------------------------------------------------------------------------------

        // Global variables for runtime communication
        // g_bincount_full = 1000;

        initialized = true;
    }
}

//-----------------------------------------------------------------------------

template <typename hashtype>
static bool test( const HashInfo * hInfo, const flags_t flags ) {
    bool result  = true;

    if(hInfo->name == std::string("exampleHash")) {
        printf("Running example hash. The parameters are set for faster computation");
        examplehash_initialisation();
    }

    if (g_testAll) {
        printf("-------------------------------------------------------------------------------\n");
    }

    if (!hInfo->Init()) {
        printf("Hash initialization failed! Cannot continue.\n");
        exit(1);
    }

    //-----------------------------------------------------------------------------
    // Sanity tests
    FILE * outfile;
    if (g_testAll) {
        outfile = stdout;
    } else {
        outfile = stderr;
    }

    fprintf(outfile, "--- Testing %s :\"%s\" ---", hInfo->name, hInfo->desc);
        
    fprintf(outfile, "\n\n");

    //-----------------------------------------------------------------------------
    // Test for Collision Curve for LSH family of hashes

    if (g_testLSHCollision) {
        result &= LSHCollisionTest<hashtype>(hInfo, flags);
        if (!result) { exit(-1); }   // Failure case.
    }

    //-----------------------------------------------------------------------------
    // Test for Nearest Neighbour search for LSH family of hashes

    if (g_testLSHApproxNearestNeighbour) {
        result &= LSHApproxNearestNeighbourTest<hashtype>(hInfo, flags);
        if (!result) { exit(-1); }   // Failure case.
    }



    return result;
}

//-----------------------------------------------------------------------------

static bool testHash( const char * name, const flags_t flags ) {
    const HashInfo * hInfo;

    if ((hInfo = findHash(name)) == NULL) {
        printf("Invalid hash '%s' specified\n", name);
        return false;
    }

    // printf("Testing hash %s\n", name);
    // printf("%u\n",hInfo->bits);

    // If you extend these statements by adding a new bitcount/type, you
    // need to adjust HASHTYPELIST in util/Instantiate.h also.
    if (hInfo->bits == 32) {
        return test<Blob<32>>(hInfo, flags);
    }
    if (hInfo->bits == 64) {
        return test<Blob<64>>(hInfo, flags);
    }
    if (hInfo->bits == 128) {
        return test<Blob<128>>(hInfo, flags);
    }
    if (hInfo->bits == 160) {
        return test<Blob<160>>(hInfo, flags);
    }
    if (hInfo->bits == 224) {
        return test<Blob<224>>(hInfo, flags);
    }
    if (hInfo->bits == 256) {
        return test<Blob<256>>(hInfo, flags);
    }

    printf("Invalid hash bit width %d for hash '%s'", hInfo->bits, hInfo->name);

    return false;
}

//-----------------------------------------------------------------------------

// uint32_t generate_single_random32(uint64_t seed) {
// 	uint32_t KeyBytes = 4; // 32 bits
// 	Rand r( seed, KeyBytes );
// 	RandSeq rs = r.get_seq(SEQ_DIST_1, KeyBytes);
// 	uint32_t random_value;
// 	rs.write((void*)&random_value, 0, 1);
// 	return random_value;
// }

int main( int argc, const char ** argv ){

    // /*----------------------------------------------*/
    // uint32_t my_random = generate_single_random32(1);
    // printf("Generated random number: 0x%08x\n", my_random);

    const char * hashToTest  = "";

    std::string config_path = "../config.toml";
    bool config_loaded = false;

    if (argc < 2) {
        printf("\n Error: No test hash given on command line\n");
        usage();
    }

    flags_t flags = FLAG_REPORT_PROGRESS;

    for (int argnb = 1; argnb < argc; argnb++) {

        const char * const arg = argv[argnb];

        if (strncmp(arg, "--", 2) == 0) {

            if (strncmp(arg, "--help", 6) == 0) {
                usage();
                exit(0);
            }
            if (strncmp(arg, "--list", 6) == 0) {
                listHashes(false);
                exit(0);
            }
            if (strncmp(arg, "--listnames", 11) == 0) {
                listHashes(true);
                exit(0);
            }
            if (strncmp(arg, "--tests", 7) == 0) {
                printf("Valid tests:\n");
                for (auto & g_testopt : g_testopts) {
                    printf("  %s\n", g_testopt.name);
                }
                exit(0);
            }
            if (strncmp(arg, "--version", 9) == 0) {
                printf("BioLSHasher %s\n", VERSION);
                exit(0);
            }
            if (strncmp(arg, "--verbose", 9) == 0) {
                flags |= FLAG_REPORT_VERBOSE;
                flags |= FLAG_REPORT_MORESTATS;
                flags |= FLAG_REPORT_DIAGRAMS;
                continue;
            }
            if (strncmp(arg, "--ncpu=", 7) == 0) {
                #if defined(HAVE_THREADS)
                    errno = 0;
                    char *   endptr;
                    long int Ncpu = strtol(&arg[7], &endptr, 0);
                    if ((errno != 0) || (arg[7] == '\0') || (*endptr != '\0') || (Ncpu < 1)) {
                        printf("Error parsing cpu number \"%s\"\n", &arg[7]);
                        exit(1);
                    }
                    // if (Ncpu > 32) { //In our machine lets not bound this.
                    //     printf("WARNING: limiting to 32 threads\n");
                    //     Ncpu = 32;
                    // }
                    g_NCPU = Ncpu;
                    continue;
                #else
                    printf("WARNING: compiled without threads; ignoring --ncpu\n");
                    continue;
                #endif
            }
            if (strncmp(arg, "--test=", 7) == 0) {
                // If a list of tests is given, only test those
                g_testAll = false;
                parse_tests(&arg[7], true);
                continue;
            }
            if (strncmp(arg, "--notest=", 9) == 0) {
                parse_tests(&arg[9], false);
                continue;
            }
            if (strncmp(arg, "--configinit", 12) == 0) {
                // This arg is used for initialising the configuration file with default values.
                //TODO: 
                continue;
            }
            if (strncmp(arg, "--configfile=", 13) == 0) {
                // This arg is used for initialising the configuration file with default values.
                config_path = std::string(&arg[13]);
                try {
                    load_config(config_path);
                    config_loaded = true;
                } catch (const std::exception& e) {
                    std::cerr << "Fatal error during configuration initialization: " << e.what() << "\n";
                    return 1;
                }
                continue;
            }
            // invalid command
            printf("Invalid command \n");
            usage();
            exit(1);
        }

        // Not a command ? => interpreted as hash name
        hashToTest = arg;
    }

    // Check that at least one test was selected
    bool anyTestSelected = false;
    for (size_t i = 0; i < sizeof(g_testopts) / sizeof(TestOpts); i++) {
        if (g_testopts[i].var) {
            anyTestSelected = true;
            break;
        }
    }
    if (!anyTestSelected) {
        printf("\n Error: No test specified. Use --test=<testname> to select a test.\n");
        printf(" Valid test names: ");
        for (size_t i = 0; i < sizeof(g_testopts) / sizeof(TestOpts); i++) {
            if (i > 0) printf(", ");
            printf("%s", g_testopts[i].name);
        }
        printf("\n");
        exit(1);
    }

    if (hashToTest[0] == '\0') {
        printf("Error: No hash name specified.\n");
        usage();
        exit(1);
    }

    if(config_loaded) {
        printf("Configuration loaded successfully from %s\n", config_path.c_str());
    } else {
        printf("Using the configuration file in project root. To load from a custom file, use --configfile=<path_to_config> \n\n");
        try {
            // 1. Load the configuration
            load_config(config_path);
        } catch (const std::exception& e) {
            std::cerr << "Fatal error during configuration initialization: " << e.what() << "\n";
            return 1;
        }
    }

    bool result = testHash(hashToTest, flags);


    return 0;
}