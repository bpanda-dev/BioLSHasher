//
// Created by Bikram Kumar Panda on 07/04/26.
//

#ifndef BIOLSHASHER_HASHINFO_H

#define BIOLSHASHER_HASHINFO_H

#include <set>
#include <cstdlib>
#include <variant>
#include <ostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <cinttypes>


typedef uint64_t  seed_t; 

#if defined(_MSC_VER)
# define FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
# define FORCE_INLINE inline __attribute__((always_inline))
#else
# define FORCE_INLINE inline
#endif


inline std::ostream& operator<<(std::ostream& os, const std::variant<int, float, double>& v) {
    std::visit([&os](const auto& val) { os << val; }, v);
    return os;
}


//--------------------------------------------------------------------------------------------------------------------//
#define HASH_FLAGS                       \
    FLAG_EXPAND(HASH_FLOATING_POINT)     \
    FLAG_EXPAND(HASH_LOOKUP_TABLE)       \
    FLAG_EXPAND(HASH_LOCALITY_SENSITIVE)    \


#define IMPL_FLAGS                          \
    FLAG_EXPAND(IMPL_SEED_WITH_HINT)        \
    FLAG_EXPAND(IMPL_SLOW)                  \
    FLAG_EXPAND(IMPL_VERY_SLOW)             \
    FLAG_EXPAND(IMPL_MULTIPLY)         \
    FLAG_EXPAND(IMPL_MULTIPLY_64_64)   \
    FLAG_EXPAND(IMPL_MULTIPLY_64_128)  \
    FLAG_EXPAND(IMPL_LICENSE_PUBLIC_DOMAIN) \
    FLAG_EXPAND(IMPL_LICENSE_BOOST)         \
    FLAG_EXPAND(IMPL_LICENSE_BSD)           \
    FLAG_EXPAND(IMPL_LICENSE_MIT)           \
    FLAG_EXPAND(IMPL_LICENSE_APACHE2)       \
    FLAG_EXPAND(IMPL_LICENSE_ZLIB)          \
    FLAG_EXPAND(IMPL_LICENSE_GPL3)

#define FLAG_EXPAND(name) FLAG_ENUM_ ## name,
typedef enum {
    HASH_FLAGS
} hashflag_enum_t;
typedef enum {
    IMPL_FLAGS
} implflag_enum_t;
#undef FLAG_EXPAND

#define FLAG_EXPAND(name) FLAG_ ## name = (1ULL << FLAG_ENUM_ ## name),
typedef enum : uint64_t {
    HASH_FLAGS
} HashFlags;
typedef enum : uint64_t {
    IMPL_FLAGS
} ImplFlags;
#undef FLAG_EXPAND
//--------------------------------------------------------------------------------------------------------------------//

class HashInfo;
typedef bool       (* HashInitFn)();

typedef seed_t     (* HashSeedfixFn)( const HashInfo * hinfo, const seed_t seed );
typedef uintptr_t  (* HashSeedFn)( const seed_t seed ); //Important for seeding.

typedef void       (* HashFn)( const void * in, const size_t len, const seed_t seed, void * out );
typedef double     (* SimilarityFn)(const std::string& seq1, const std::string& seq2, const uint32_t in1_len, const uint32_t in2_len);

typedef bool       (* CheckEqualityFn)(void* inp1, void* inp2);

class HashInfo {
    friend class HashFamilyInfo;

  protected:
    static const char * _fixup_name( const char * in );

  public:
    const char *      name;
    const char *      family;
    const char *      desc;
    uint64_t          hash_flags{};
    uint64_t          impl_flags{};
    uint32_t          bits{};
    HashInitFn        initfn;
    HashSeedfixFn     seedfixfn;
    HashSeedFn        seedfn;
    HashFn            hashfn;
    const char *      similarity_name;
    SimilarityFn      similarityfn;
    CheckEqualityFn   check_equality_fn;


    enum fixupseed : size_t {
        SEED_ALLOWFIX = 0,  // Seed via a SeedfixFn, if the hash has one
        SEED_FORCED   = 1   // Seed using the given seed, always
    };

    std::vector<std::string>    parameterNames;
    std::vector<std::string>    parameterDescriptions;
    std::vector<std::variant<int, float, double>> parameterValues;

    HashInfo( const char * n, const char * f ) :
        name( _fixup_name( n ) ), family( f ), desc( "" ),
        initfn( nullptr ), seedfixfn( nullptr ), seedfn( nullptr ),
        hashfn( nullptr ), similarity_name(""), similarityfn(nullptr), check_equality_fn(nullptr) {}

    ~HashInfo() {
        free(const_cast<char *>(name));
    }

    [[nodiscard]] FORCE_INLINE HashFn hashFn() const {  //REM: Remember Here i removed hashfn endianness.   //{ return _is_native(endian) ? hashfn_native : hashfn_bswap;
        return hashfn;
    }

    [[nodiscard]] FORCE_INLINE SimilarityFn distFn() const {
        return similarityfn;
    }

    [[nodiscard]] FORCE_INLINE CheckEqualityFn equalityFn() const {
        return check_equality_fn;
    }

    [[nodiscard]] FORCE_INLINE bool Init() const {
        if (initfn != nullptr) {
            return initfn();
        }
        return true;
    }

    FORCE_INLINE seed_t Seed( seed_t seed, enum fixupseed fixup = SEED_ALLOWFIX, uint64_t hint = 0 ) const {
        if (seedfixfn != NULL) {
            if (impl_flags & FLAG_IMPL_SEED_WITH_HINT) {
                seedfixfn(NULL, hint);
            }
            else if (fixup == SEED_ALLOWFIX) {
                seed = seedfixfn(this, seed);
            }
        }
        if (seedfn != NULL) {
            seed_t newseed = (seed_t)seedfn(seed);
            if (newseed != 0) {
                seed = newseed;
            }
        }
        return seed;
    }

    FORCE_INLINE seed_t getFixedSeed( seed_t seed ) const {
        if (seedfixfn != NULL) {
            seed = (seed_t)seedfixfn(this, seed);
        }
        return seed;
    }

    [[nodiscard]] FORCE_INLINE bool isSlow() const {
        return !!(impl_flags & (FLAG_IMPL_SLOW | FLAG_IMPL_VERY_SLOW));
    }

    [[nodiscard]] FORCE_INLINE bool isVerySlow() const {
        return !!(impl_flags & FLAG_IMPL_VERY_SLOW);
    }

    // [[nodiscard]] FORCE_INLINE bool onlyShortSequenceLength() const {
    //     return !!(impl_flags & FLAG_IMPL_SMALL_SEQUENCE_LENGTH);
    // }

    [[nodiscard]] FORCE_INLINE bool isLocalSensitive() const {
        return !!(hash_flags & FLAG_HASH_LOCALITY_SENSITIVE);
    }

    // [[nodiscard]] FORCE_INLINE bool hasTokenisationProperty() const {
    //     return !!(hash_flags & FLAG_HASH_TOKENISATION_PROPERTY);
    // }
    //
    // [[nodiscard]] FORCE_INLINE bool hasHammingSimilarity() const {
    //     return !!(hash_flags & FLAG_HASH_HAMMING_SIMILARITY);
    // }
    //
    // [[nodiscard]] FORCE_INLINE bool hasJaccardSimilarity() const {
    //     return !!(hash_flags & FLAG_HASH_JACCARD_SIMILARITY);
    // }
    //
    // [[nodiscard]] FORCE_INLINE bool hasAngularSimilarity() const {
    //     return !!(hash_flags & FLAG_HASH_ANGULAR_SIMILARITY);
    // }
    //
    // [[nodiscard]] FORCE_INLINE bool hasCosineSimilarity() const {
    //     return !!(hash_flags & FLAG_HASH_COSINE_SIMILARITY);
    // }
    // [[nodiscard]] FORCE_INLINE bool hasEditSimilarity() const {
    //     return !!(hash_flags & FLAG_HASH_EDIT_SIMILARITY);
    // }
    // [[nodiscard]] FORCE_INLINE bool hasUniverseVectorOptimisation() const {
    //     return !!(hash_flags & FLAG_HASH_UNIVERSE_VECTOR_OPTIMISATION);
    // }

    // Method to register parameters
    void registerParameters(const std::vector<std::string>& names,
                            const std::vector<std::string>& descriptions,
                            const std::vector<std::variant<int, float, double>>& defaults) {
        if ((names.size() != descriptions.size()) || (names.size() != defaults.size()) || (descriptions.size() != defaults.size())) {
            throw std::runtime_error("Parameter registration mismatch: names, descriptions, and defaults must have the same size.");
        }
        parameterNames = names;
        parameterDescriptions = descriptions;
        parameterValues = defaults;
    }

    // Method to print all parameters
    void printParameters(std::ostream& out) const {
        // Print parameter names
        out << ":4.1: ";
        for (size_t i = 0; i < parameterNames.size(); ++i) {
            out << parameterNames[i];
            if (i != parameterNames.size() - 1) {
                out << ",";
            }
        }
        out << "\n";

        // Print parameter descriptions
        out << ":4.2: ";
        for (size_t i = 0; i < parameterDescriptions.size(); ++i) {
            out << parameterDescriptions[i];
            if (i != parameterDescriptions.size() - 1) {
                out << ",";
            }
        }
        out << "\n";

        // Print parameter default values
        out << ":4.3: ";
        for (size_t i = 0; i < parameterValues.size(); ++i) {
            out << parameterValues[i];
            if (i != parameterValues.size() - 1) {
                out << ",";
            }
        }
        out << "\n";
    }

}; // class HashInfo


class HashFamilyInfo {
public:
    const char * name;
    const char * src_url;
    enum SrcStatus : uint32_t {
        SRC_UNKNOWN,
        SRC_FROZEN,    // Very unlikely to change
        SRC_STABLEISH, // Fairly unlikely to change
        SRC_ACTIVE,    // Likely to change
    }  src_status;

    HashFamilyInfo( const char * n ) :
        name( _fixup_name( n )),
        src_url( NULL ), src_status( SRC_UNKNOWN ) {}

private:
    static const char * _fixup_name( const char * in );
}; // class HashFamilyInfo


#endif //BIOLSHASHER_HASHINFO_H
