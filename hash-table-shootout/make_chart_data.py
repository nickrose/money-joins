import sys, json
from collections import OrderedDict

lines = [ line.strip() for line in sys.stdin if line.strip() ]

by_benchtype = {}

for line in lines:
    benchtype, nkeys, program, load_factor, nbytes, runtime = line.split(',')
    nkeys = int(nkeys)
    nbytes = int(nbytes)
    runtime = float(runtime)
    load_factor = float(load_factor)

    by_benchtype.setdefault("%s_runtime" % benchtype, {}).setdefault(program, []).append([nkeys, runtime, load_factor])
    if benchtype in ('insert_random_full', 'insert_small_string', 'insert_string',
                     'insert_random_full_reserve', 'insert_small_string_reserve', 'insert_string_reserve'):
        by_benchtype.setdefault("%s_memory"  % benchtype, {}).setdefault(program, []).append([nkeys, nbytes, load_factor])

proper_names = OrderedDict([
    ('std_unordered_map', 'std::unordered_map'),
    ('google_dense_hash_map', 'google::dense_hash_map'),
    ('qt_qhash', 'QHash'),
    ('spp_sparse_hash_map', 'spp::sparse_hash_map'),
    ('tsl_hopscotch_map', 'tsl::hopscotch_map'),
    ('tsl_robin_map', 'tsl::robin_map'),
    ('tsl_hopscotch_map_store_hash', 'tsl::hopscotch_map (with StoreHash)'),
    ('tsl_robin_map_store_hash', 'tsl::robin_map (with StoreHash)'),
    ('tsl_hopscotch_map_mlf_0_5', 'tsl::hopscotch_map (0.5 mlf)'),
    ('tsl_robin_map_mlf_0_9', 'tsl::robin_map (0.9 mlf)'),
    ('tsl_ordered_map', 'tsl::ordered_map'),
    ('tsl_robin_pg_map', 'tsl::robin_pg_map'),
    ('ska_flat_hash_map', 'ska::flat_hash_map'),
    ('google_dense_hash_map_mlf_0_9', 'google::dense_hash_map (0.9 mlf)'),
    ('ska_flat_hash_map_power_of_two', 'ska::flat_hash_map (power of two)'),
    ('google_sparse_hash_map', 'google::sparse_hash_map'),
    ('boost_unordered_map', 'boost::unordered_map'),
    ('emilib_hash_map', 'emilib::HashMap'),
    ('tsl_array_map', 'tsl::array_map'),
    ('tsl_array_map_mlf_1_0', 'tsl::array_map (1.0 mlf)'),
])

# do them in the desired order to make the legend not overlap the chart data
# too much
program_slugs = [
    'std_unordered_map',
    'google_dense_hash_map',
    'qt_qhash',
    'spp_sparse_hash_map',
    'tsl_hopscotch_map',
    'tsl_robin_map',
    'tsl_hopscotch_map_store_hash',
    'tsl_robin_map_store_hash',
    'tsl_hopscotch_map_mlf_0_5',
    'tsl_robin_map_mlf_0_9',
    'tsl_ordered_map',
    'tsl_robin_pg_map',
    'ska_flat_hash_map',
    'google_dense_hash_map_mlf_0_9',
    'ska_flat_hash_map_power_of_two',
    'google_sparse_hash_map',
    'boost_unordered_map',
    'emilib_hash_map',
    'tsl_array_map',
    'tsl_array_map_mlf_1_0'
]

# hashmap which will be shown (checkbox checked)
default_programs_show = [
    'std_unordered_map',
    'google_dense_hash_map',
    'qt_qhash',
    'spp_sparse_hash_map',
    'tsl_hopscotch_map',
    'tsl_robin_map',
    'tsl_hopscotch_map_store_hash',
    'tsl_robin_map_store_hash']

chart_data = {}

for i, (benchtype, programs) in enumerate(by_benchtype.items()):
    chart_data[benchtype] = []
    for j, program in enumerate(program_slugs):
        if program not in programs:
            continue
        
        data = programs[program]
        chart_data[benchtype].append({
            'program': program,
            'label': proper_names[program],
            'data': [],
        })

        for k, (nkeys, value, load_factor) in enumerate(data):
            chart_data[benchtype][-1]['data'].append([nkeys, value, load_factor])

json_text = json.dumps(chart_data)
json_text = json_text.replace("}], ", "}], \n")
print('chart_data = ' + json_text + ';')
print('\nprograms = ' + json.dumps(proper_names, indent=1) + ';')
print('\ndefault_programs_show = ' + str(default_programs_show) + ';')
