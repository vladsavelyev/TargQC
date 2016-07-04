from os.path import abspath, dirname, pardir, join

code_base_path = abspath(join(dirname(abspath(__file__)), pardir))

padding = 200
depth_thresholds = [1, 5, 10, 25, 50, 100, 500, 1000, 5000, 10000, 50000]
downsample_to = 5e5
genome = None
features_bed_fpath = None
cds_bed_fpath = None
fai_fpath = None
bwa_prefix = None
dedup = True

reuse_intermediate = False
debug = False

parallel_cfg = None