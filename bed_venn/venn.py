#!/usr/bin/env python
from __future__ import print_function
import subprocess
import json

from ngs_utils.reporting.reporting import write_static_html_report
from os.path import dirname, basename, join, splitext, isfile, abspath, pardir


def call(cmdl):
    print(cmdl)
    subprocess.call(cmdl, shell=True)


def check_output(cmdl):
    print(cmdl)
    return subprocess.check_output(cmdl, shell=True)


def bedsize(bed):
    size = check_output("cat " + bed + " | awk -F'\\t' 'BEGIN{ SUM=0 }{ SUM+=$3-$2 }END{ print SUM }'")
    print('Size of ' + basename(bed) + ': ' + size)
    return int(size)


total_calls = 0


def intersect_pair(work_dirpath, bed1, bed2):
    bed1, bed2 = sorted([bed1, bed2])
    output_fpath = join(work_dirpath, splitext(basename(bed1))[0] + '__' + basename(bed2))
    if not isfile(output_fpath):
        print('intersect_pair: ' + splitext(basename(bed1))[0] + ' and ' + splitext(basename(bed2))[0])
        call('/usr/local/bin/bedtools intersect -a ' + bed1 + ' -b ' + bed2 + ' > ' + output_fpath)
        # global total_calls
        # total_calls += 1
    return output_fpath


def calc_set_intersection(work_dirpath, sorted_beds_subset, intersection_bed_by_subset):
    if len(sorted_beds_subset) == 1:
        return sorted_beds_subset[0]
    if tuple(sorted_beds_subset) in intersection_bed_by_subset:
        return intersection_bed_by_subset[tuple(sorted_beds_subset)]

    intersection_bed = None
    for bed in sorted_beds_subset:
        remaining_subset = [b for b in sorted_beds_subset if b != bed]  # sorted_beds_subset.index(b) > sorted_beds_subset.index(bed)]
        if not remaining_subset:
            continue
        print('comparing ' + basename(bed) + ' and ' + str([basename(k) for k in remaining_subset]))
        subset_intersection_bed = intersection_bed_by_subset.get(tuple(sorted(remaining_subset)))
        if not subset_intersection_bed:
            subset_intersection_bed = calc_set_intersection(work_dirpath, remaining_subset, intersection_bed_by_subset)
        intersection_bed = intersect_pair(work_dirpath, subset_intersection_bed, bed)
        intersection_bed_by_subset[tuple(sorted(remaining_subset))] = subset_intersection_bed
    intersection_bed_by_subset[tuple(sorted_beds_subset)] = intersection_bed
    return intersection_bed


def save_venn_diagram_data(size_by_set, names_map):
    data = []
    for venn_set, size in size_by_set.items():
        set_info = dict()
        set_info['size'] = size
        # if isinstance(venn_set, tuple):
        set_info['sets'] = [names_map.get(n, n) for n in venn_set]
        # else:
        #     set_info['sets'] = [venn_set]
        # if isinstance(venn_set, int):
        #     set_info['label'] = label_by_set[venn_set]
        data.append(set_info)
    return json.dumps(sorted(data, key=lambda x: x['sets']))


def run(work_dirpath, bed_fpaths):
    intersection_bed_by_subset = dict()
    calc_set_intersection(work_dirpath, sorted(bed_fpaths), intersection_bed_by_subset)

    intersection_size_by_subset = dict()
    # label_by_subset = dict()
    for bed_set, intersection_bed in intersection_bed_by_subset.items():
        bed_set = tuple([splitext(basename(b))[0] for b in bed_set])
        intersection_size_by_subset[bed_set] = bedsize(intersection_bed)
        # label_by_subset[bed_set] = basename(splitext(intersection_bed)[0])
        print(str(bed_set) + ': ' + basename(intersection_bed) + ', size: ' + str(intersection_size_by_subset[bed_set]))
    return intersection_size_by_subset


def write_html(output_dir, json_txt, bed_fpaths):
    output_html = join(output_dir, 'venn.html')

    # def _get_static_file(_fname):
    #     return join(dirname(abspath(reporting.__file__)), 'static', _fname)
    write_static_html_report({
            'title': 'Venn comparison for ' + ', '.join([basename(bf) for bf in bed_fpaths]),
            'diagram_data': json_txt,
        }, output_html,
        tmpl_fpath=join(dirname(abspath(__file__)), 'venn_template.html'),
        extra_js_fpaths=['venn.js', 'd3.min.js', 'd3.tip.js', 'draw_venn_diagram.js'])

    # output_html = join(output_dir, 'venn.html')
    # with open(join(dirname(__file__), 'venn_template.html')) as tmpl_f:
    #     html = tmpl_f.read()
    # html = html.replace('{title}', 'Venn comparison for ' + ', '.join([basename(bf) for bf in bed_fpaths]))
    # html = html.replace('{diagram_data}', json_txt)
    # html = html.replace('{draw_venn_diagram.js}', open(join(dirname(__file__), 'draw_venn_diagram.js')).read())
    # html = html.replace('{venn.js}', open(join(dirname(__file__), 'venn.js')).read())
    # with open(output_html, 'w') as out_f:
    #     out_f.write(html)
    
    return output_html

