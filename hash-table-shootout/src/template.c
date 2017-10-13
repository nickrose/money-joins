#include <sys/time.h>
#include <sys/types.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <getopt.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

// TODO When generating random values to insert in the map there is no check
// to see if duplicate random values are generated. Could improve that (but the probability is so so
// low and the impact nearly null that it's not really worth it).

static const std::array<char, 62> ALPHANUMERIC_CHARS = {
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'
};

/**
 * SMALL_STRING_SIZE should be small enough so that there is no heap allocation when a std::string is created.
 */
static const std::size_t SMALL_STRING_SIZE = 15;
static const std::size_t STRING_SIZE = 50;

static const std::int64_t SEED = 0;
static std::mt19937_64 generator(SEED);


std::size_t get_memory_usage_bytes() {
    std::ifstream file("/proc/self/statm");
    
    std::size_t memory;
    file >> memory; // Ignore first
    file >> memory;
    
    return memory * getpagesize();
}

std::string get_random_alphanum_string(std::size_t size) {
    std::uniform_int_distribution<std::size_t> rd_uniform(0, ALPHANUMERIC_CHARS.size() - 1);
    
    std::string str(size, '\0');
    for(std::size_t i = 0; i < size; i++) {
        str[i] = ALPHANUMERIC_CHARS[rd_uniform(generator)];
    }
    
    return str;
}

/**
 * Generate a vector of ints from [0, nb_ints) and shuffle it
 */
std::vector<std::int64_t> get_dense_ints(std::size_t num_ints) {
    std::vector<std::int64_t> results(num_ints);
    std::iota(results.begin(), results.end(), 0);
    std::shuffle(results.begin(), results.end(), generator);
    return results;
}

/**
 * Generate vector of random ints from [min, max]
 */
std::vector<std::int64_t> get_sparse_ints(std::size_t nb_ints, 
                                          std::int64_t min = 0, 
                                          std::int64_t max = std::numeric_limits<std::int64_t>::max()) 
{
    std::uniform_int_distribution<std::int64_t> rd_uniform(min, max);
    
    std::vector<std::int64_t> random_ints(nb_ints);
    for(std::size_t i = 0; i < random_ints.size(); i++) {
        random_ints[i] = rd_uniform(generator);
    }
    
    return random_ints;
}


std::vector<std::string> get_random_alphanum_strings(std::size_t nb_strings, std::size_t string_size) {
    std::vector<std::string> random_strings(nb_strings);
    for(std::size_t i = 0; i < random_strings.size(); i++) {
        random_strings[i] = get_random_alphanum_string(string_size);
    }
    
    return random_strings;
}

class measurements {
public:    
    measurements(uint64_t keys): m_keys(keys),
                                 m_memory_usage_bytes_start(get_memory_usage_bytes()),
                                 m_chrono_start(std::chrono::high_resolution_clock::now()) { }
    
    ~measurements() {
        const auto chrono_end = std::chrono::high_resolution_clock::now();
        const std::size_t memory_usage_bytes_end = get_memory_usage_bytes();
        
        const double nb_seconds = std::chrono::duration<double>(chrono_end - m_chrono_start).count();
        // On reads or delete the used bytes could be less than initially.
        const std::size_t used_memory_bytes = (memory_usage_bytes_end > m_memory_usage_bytes_start)?
                                                    memory_usage_bytes_end - m_memory_usage_bytes_start:0;

        const auto mtps = (m_keys / 1000.0 / 1000.0) / nb_seconds;

        std::cout << std::fixed << std::setprecision(2) 
                  << mtps << " mtps "
                  << nb_seconds << " sec "
                  << (used_memory_bytes / 1024.0 / 1024.0) << " MB ";
    }
    
private: 
    uint64_t m_keys;   
    std::size_t m_memory_usage_bytes_start;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_chrono_start;
};

struct arguments {
    std::string bench;
    std::int64_t num_keys;
    std::int64_t probe_sf;
    std::int64_t num_runs;
    uint32_t key_size;
    uint32_t val_size;
    bool is_static;
    bool verbose;
};

void print_help(char *progname) {
    printf("Usage: %s [options]\n", progname);

    printf("\
    Benchmark : insert_dense, read_dense, insert_sparse, read_sparse etc.      \n\
       -b --bench=<name>    Run the given benchmark                            \n\
                                                                               \n\
       -n --nkeys=<N>       Number of keys to insert <N> [2^24]                \n\
       -p --probe_sf=<SF>   Scale factor to use when generating probe keys [1] \n\
       -r --nruns=<R>       Number of times to run the benchmark [1]           \n\
       -k --keysize=<K>     Size of key in bytes                               \n\
       -v --valsize<V>      Size of value in bytes                             \n\
       --static             Use static table (i.e., reserve space beforehand)  \n\
                                                                               \n\
    Other flags:                                                               \n\
       -h --help            Show this message                                  \n\
       --verbose            Display log information                            \n\
    \n");
}

void parse_arguments(int argc, char **argv, arguments &args) {
    static int verbose_flag = 0;
    static int static_flag = 0;
    static struct option long_options[] = {
        // These set flags
        {"verbose",   no_argument, &verbose_flag, 1},
        {"static",    no_argument, &static_flag, 1},
        // These have arguments
        {"bench",     required_argument, 0, 'b'},
        {"nkeys",     required_argument, 0, 'n'},
        {"probe_sf",  required_argument, 0, 'p'},
        {"nruns",     required_argument, 0, 'r'},  
        {"keysize",   required_argument, 0, 'k'},  
        {"valsize",   required_argument, 0, 'v'},        
        {"help",      required_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    // getopt_long stores the option index here
    int option_index = 0;
 
    while (1) {
        int c = getopt_long(argc, argv, "b:n:p:r:h", long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: {
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0) break;
                printf ("option %s", long_options[option_index].name);
                if (optarg) printf (" with arg %s", optarg);
                printf ("\n");
                break;
            }
            case 'b': {
                args.bench = std::string{optarg};
                break;
            }
            case 'n': {
                args.num_keys = std::atoi(optarg);
                break;
            }
            case 'p': {
                args.probe_sf = std::atoi(optarg);
                break;
            }
            case 'r': {
                args.num_runs = std::atoi(optarg);
                break;   
            }
            case 'k': {
                args.key_size = std::atoi(optarg);
                break;
            }
            case 'v': {
                args.val_size = std::atoi(optarg);
                break;
            }
            case 'h': 
            case '?': {
                print_help(argv[0]);
                exit(EXIT_SUCCESS);
                break;
            }
            default: break;
        }
    }
    args.verbose = verbose_flag;
    args.is_static = static_flag;
}


int main(int argc, char ** argv) {
    // Setup default arguments
    arguments args;

    args.bench = "";
    args.num_keys = 1 << 24;
    args.probe_sf = 10;
    args.num_runs = 5;
    args.key_size = 8;
    args.val_size = 8;
    args.is_static = true;
    args.verbose = false;

    // Parse everything
    parse_arguments(argc, argv, args);

    if (args.key_size % 8 != 0 || args.key_size > 256 || 
        args.val_size % 8 != 0 || args.val_size > 256) {
        printf("Key/value size must be multiple of 4 in range: [4,256]. Provided: %u/%u\n", 
               args.key_size, args.val_size);
        exit(1);
    }
    
    // Do things
    const std::int64_t num_keys = args.num_keys;
    const std::int64_t num_probe_keys = args.num_keys * args.probe_sf;
    const std::string test_type = args.bench;
    const std::int64_t value = 1;

    if (args.verbose) {
        std::cout << std::fixed << std::setprecision(2)
                  << num_keys << " " << sizeof(KeyType) << "/" << sizeof(ValueType) << " "
                  << ((sizeof(KeyType)+sizeof(ValueType))*num_keys)/1024.0/1024.0 << " MB "
                  << args.probe_sf << "x (" << num_probe_keys << ") "
                  << (args.is_static ? "static" : "dynamic") << " "
#ifdef MURMUR3
                  << "Murmur3" << std::endl;
#elif defined(CLHASH)
                  << "CLHash" << std::endl;
#elif defined(XXHASH)
                  << "xxHash" << std::endl;
#elif defined(FARMHASH)
                  << "FarmHash" << std::endl;
#elif defined(MULTHASH)
                  << "Mult" << std::endl;
#elif defined(CRCHASH)
                  << "CRC32" << std::endl;
#elif defined(SIMPLEHASH)
                  << "XOR" << std::endl;
#else
                  << "std::hash" << std::endl;
#endif
    }


    SETUP


    /**
     * Integers
     */
    if(test_type == "insert_dense") {
        const std::vector<std::int64_t> keys = get_dense_ints(num_keys);
        
        measurements m(num_keys);
        if (args.is_static) {
            RESERVE_INT(num_keys);
        }
        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_INT_INTO_HASH(keys[i], value);
        }
    }

    else if(test_type == "probe_dense") {
        const std::vector<std::int64_t> insert_keys = get_dense_ints(num_keys);         
        if (args.is_static) {
            RESERVE_INT(num_keys);
        }
        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_INT_INTO_HASH(insert_keys[i], value);
        }
        
        // Generate probe keys
        std::vector<std::int64_t> probe_keys;
        for (uint32_t i = 0; i < args.probe_sf; i++) {
            probe_keys.insert(probe_keys.end(), insert_keys.begin(), insert_keys.end());
        }
        if (args.verbose) {
            std::cout << "Generated " << probe_keys.size() << " probe keys\n";
        }
        // Shuffle probe keys
        std::shuffle(probe_keys.begin(), probe_keys.end(), generator);
        
        measurements m(num_probe_keys);
        for(std::int64_t i = 0; i < num_probe_keys; i++) {
            FIND_INT_EXISTING_FROM_HASH(probe_keys[i]);
        }
    }

    else if(test_type == "insert_sparse") {
        const std::vector<std::int64_t> keys = get_sparse_ints(num_keys);
        
        
        measurements m(num_keys);
        if (args.is_static) {
            RESERVE_INT(num_keys);
        }
        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_INT_INTO_HASH(keys[i], value);
        }
    }

    else if(test_type == "probe_sparse") {
        const std::vector<std::int64_t> insert_keys = get_sparse_ints(num_keys);
        if (args.is_static) {
            RESERVE_INT(num_keys);
        }
        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_INT_INTO_HASH(insert_keys[i], value);
        }
        
        // Generate probe keys
        std::vector<std::int64_t> probe_keys;
        for (uint32_t i = 0; i < args.probe_sf; i++) {
            probe_keys.insert(probe_keys.end(), insert_keys.begin(), insert_keys.end());
        }
        if (args.verbose) {
            std::cout << "Generated " << probe_keys.size() << " probe keys\n";
        }
        // Shuffle probe keys
        std::shuffle(probe_keys.begin(), probe_keys.end(), generator);
        
        measurements m(num_probe_keys);
        for(std::int64_t i = 0; i < num_probe_keys; i++) {
            FIND_INT_EXISTING_FROM_HASH(probe_keys[i]);
        }
    }

    else if(test_type == "read_miss_random_full") {
        const std::vector<std::int64_t> keys_insert = get_sparse_ints(num_keys, 0, std::numeric_limits<std::int64_t>::max());
        const std::vector<std::int64_t> keys_read = get_sparse_ints(num_keys, std::numeric_limits<std::int64_t>::min(), -3);
        
        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_INT_INTO_HASH(keys_insert[i], value);
        }
        
        
        measurements m(num_keys);
        for(std::int64_t i = 0; i < num_keys; i++) {
            FIND_INT_MISSING_FROM_HASH(keys_read[i]);
        }
    }   
    
    /**
     * Small strings
     */
    else if(test_type == "insert_small_string") {
        const std::vector<std::string> keys = get_random_alphanum_strings(num_keys, SMALL_STRING_SIZE);
        
        
        measurements m(num_keys);
        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_STR_INTO_HASH(keys[i], value);
        }
    }

    else if(test_type == "insert_small_string_reserve") {
        const std::vector<std::string> keys = get_random_alphanum_strings(num_keys, SMALL_STRING_SIZE);
        
        
        measurements m(num_keys);
        RESERVE_STR(num_keys);
        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_STR_INTO_HASH(keys[i], value);
        }
    }

    else if(test_type == "read_small_string") {
        std::vector<std::string> keys = get_random_alphanum_strings(num_keys, SMALL_STRING_SIZE);
        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_STR_INTO_HASH(keys[i], value);
        }
        
        std::shuffle(keys.begin(), keys.end(), generator);
        
        
        measurements m(num_keys);
        for(std::int64_t i = 0; i < num_keys; i++) {
            FIND_STR_EXISTING_FROM_HASH(keys[i]);
        }
    }

    else if(test_type == "read_miss_small_string") {
        const std::vector<std::string> keys_insert = get_random_alphanum_strings(num_keys, SMALL_STRING_SIZE);
        const std::vector<std::string> keys_read = get_random_alphanum_strings(num_keys, SMALL_STRING_SIZE);

        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_STR_INTO_HASH(keys_insert[i], value);
        }
        
        
        measurements m(num_keys);
        for(std::int64_t i = 0; i < num_keys; i++) {
            FIND_STR_MISSING_FROM_HASH(keys_read[i]);
        }
    }    
    
    /**
     * Strings
     */
    else if(test_type == "insert_string") {
        const std::vector<std::string> keys = get_random_alphanum_strings(num_keys, STRING_SIZE);
        
        
        measurements m(num_keys);
        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_STR_INTO_HASH(keys[i], value);
        }
    }
    
    else if(test_type == "insert_string_reserve") {
        const std::vector<std::string> keys = get_random_alphanum_strings(num_keys, STRING_SIZE);
        
        
        measurements m(num_keys);
        RESERVE_STR(num_keys);
        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_STR_INTO_HASH(keys[i], value);
        }
    }

    else if(test_type == "read_string") {
        std::vector<std::string> keys = get_random_alphanum_strings(num_keys, STRING_SIZE); 
        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_STR_INTO_HASH(keys[i], value);
        }
        
        std::shuffle(keys.begin(), keys.end(), generator);   
        
        
        measurements m(num_keys);
        for(std::int64_t i = 0; i < num_keys; i++) {
            FIND_STR_EXISTING_FROM_HASH(keys[i]);
        }
    }

    else if(test_type == "read_miss_string") {
        const std::vector<std::string> keys_insert = get_random_alphanum_strings(num_keys, STRING_SIZE);
        const std::vector<std::string> keys_read = get_random_alphanum_strings(num_keys, STRING_SIZE);

        for(std::int64_t i = 0; i < num_keys; i++) {
            INSERT_STR_INTO_HASH(keys_insert[i], value);
        }
        
        
        measurements m(num_keys);
        for(std::int64_t i = 0; i < num_keys; i++) {
            FIND_STR_MISSING_FROM_HASH(keys_read[i]);
        }
    }
    
    else {
        std::cout << "Unknown test type: " << test_type << "." << std::endl;
        std::exit(1);
    }
    
    
    const float load_factor = std::max(LOAD_FACTOR(hash), LOAD_FACTOR(str_hash));
    std::cout << std::fixed << std::setprecision(2) << load_factor << " LF" << std::endl;
}
