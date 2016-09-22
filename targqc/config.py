from os.path import abspath, dirname, pardir, join

code_base_path = abspath(join(dirname(abspath(__file__)), pardir))

# Defaults (can be altered here or via command line arguments):
padding = 200
depth_thresholds = [1, 5, 10, 20, 50, 100, 250, 500, 1000, 5000, 10000, 50000]
downsample_fraction = 0.05
downsample_pairs_num = 5e5
genome = 'hg19'
dedup = True

reuse_intermediate = False
is_debug = False
threads = 1
reannotate = False  # reannotate BED even if the number of columns is 4 or higher