import sys, os, subprocess, signal

sys.stdout  = open('/home/jungmins/ht_result2/farmhash-15-15-dev5.txt', 'w', 0)

def best_of(runs, prog, args):
    best_mtps = 0
    best_secs = 0
    best_mem = 0
    best_load_factor = 0
    for i in range(runs):
        output = subprocess.check_output([prog] + args)
        words = output.strip().split()
        
        mtps = words[0]
        secs = words[2]
        mem = words[4]
        load_factor = words[6]

        if i == 0 or mtps > best_mtps:
            best_mtps = mtps
            best_secs = secs
            best_mem = mem
            best_load_factor = load_factor

    return (best_mtps, best_secs, best_mem, best_load_factor)

benchmarks = [
    'insert_small_string',
    'read_small_string',
    'insert_string',
    'read_string',
]

progs = [
    './build/emilib_hash_map',
    './build/google_dense_hash_map_mlf_0_9',
    './build/ska_flat_hash_map_mlf_0_9',
    './build/ska_flat_hash_map_power_of_two_mlf_0_9',
    './build/tsl_hopscotch_map',
    './build/tsl_robin_map_mlf_0_9',
    './build/libcuckoo_map',
]

if __name__ == '__main__':
    jump = 1680000
    for bench in benchmarks:
        print("===================== {0} =====================".format(bench.upper()))
        print("====> # Keys: {0}, Jumping: {1}".format((1<<24)+1000, jump))
        for prog in progs:
            nkeys = (1 << 24) + 1000
            last_lf = 0
            print("--> {0}".format(prog))
            while True:
                (mtps, secs, mem, lf) = best_of(5, prog, ['-n', str(nkeys), '-b', bench, '--static'])
                if last_lf != 0 and lf < last_lf:
                    break
                print("{0}: {1} mtps, {2} sec, {3} MB".format(lf, mtps, secs, mem))
                last_lf = lf
                nkeys += jump
