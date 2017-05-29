from ngs_utils.logger import warn

metric_names = [
    'Reference size',
    'Regions size/percentage of reference (on target)',
    'Regions size/percentage of reference (on target) %',
    'Coverage Mean',
    'Coverage Mean (on target)',
    'Coverage Standard Deviation',
    'Coverage Standard Deviation (on target)',
    'Reference size',
    'Number of reads',
    'Mapped reads',
    'Mapped reads %',
    'Unmapped reads',
    'Unmapped reads %',
    'Mapped reads (on target)',
    'Mapped reads (on target) %',
    'Mapped paired reads',
    'Mapped paired reads %',
    'Paired reads',
    'Paired reads %',
    'Duplicated reads (flagged)',
    'Duplicated reads (flagged) %',
    'Duplicated reads (flagged) (on target)',
    'Duplicated reads (flagged) (on target) %',
    'Read min length',
    'Read max length',
    'Read mean length',
    'Mean Mapping Quality (on target)',
    'Mismatches (on target)',
    'Insertions (on target)',
    'Deletions (on target)',
    'Homopolymer indels (on target)',
    'Mean Mapping Quality',
    'Mismatches',
    'Insertions',
    'Deletions',
    'Homopolymer indels',
]

ALLOWED_UNITS = ['%']


def parse_qualimap_sample_report(report_fpath):
    value_by_metric = dict()

    def __get_td_tag_contents(line):
        ## examples:
        # <td class=column1>Paired reads</td>
        # <td class=column2>80,244 / 99.89%</td>
        crop_left = line.split('>')
        if len(crop_left) < 2:
            return None
        crop_right = crop_left[1].split('<')
        return crop_right[0].strip()

    def __fill_record(metric_name, line):
        val = __get_td_tag_contents(line)
        val = val.replace(' ', '').replace(',', '')
        try:
            val = val.replace(b'\xc2\xa0', '')
        except:
            val = val.replace(b'\xc2\xa0'.decode(), '')

        if metric_name == 'Read min/max/mean length':  # special case
            for metric_infix, value in zip(['min', 'max', 'mean'], val.split('/')):
                value_by_metric['Read ' + metric_infix + ' length'] = value
        else:
            if metric_name not in metric_names:
                # warn('Qualimap metric "' + metric_name + '" is not in allowed metric_names')
                return

            num_chars = []
            unit_chars = []
            i = 0
            while i < len(val) and (val[i].isdigit() or val[i] in ['.']):
                num_chars += val[i]
                i += 1
            while i < len(val):
                unit_chars += val[i]
                i += 1
            val_num = ''.join(num_chars)
            val_unit = ''.join(unit_chars)

            if val_unit and val_unit in ALLOWED_UNITS:
                # metric.unit = val_unit
                pass
            try:
                val = int(val_num)
                if val_unit == '%':
                    val = float(val) / 100
            except ValueError:
                try:
                    val = float(val_num)
                    if val_unit == '%':
                        val /= 100
                except ValueError:  # it is a string
                    val = val_num + val_unit

            value_by_metric[metric_name] = val

            if val_unit.startswith('/'):  # for values like "80,220 / 99.86%"
                meta_val = val_unit.replace('/', '').strip()
                if '%' in meta_val:
                    try:
                        val = float(meta_val.replace('%', '')) / 100.0
                    except ValueError:
                        pass
                    else:
                        value_by_metric[metric_name + ' %'] = val

    sections = [['start',                             'Summary'],
                ['globals (on target)',               'Globals (inside of regions)'],
                ['globals',                           'Globals'],
                ['coverage (on target)',              'Coverage (inside of regions)'],
                ['coverage',                          'Coverage'],
                ['mq (on target)',                    'Mapping Quality (inside of regions)'],
                ['mq',                                'Mapping Quality'],
                ['mismatches and indels (on target)', 'Mismatches and indels (inside of regions)'],
                ['mismatches and indels',             'Mismatches and indels'],
                ['finish',                            'Coverage across reference']]  # plots are starting from this line
    on_target_stats_suffix = ' (on target)'
    coverage_stats_prefix = 'Coverage '
    with open(report_fpath) as f:
        cur_section = None
        cur_metric_name = None
        for line in f:
            if 'mapped' in line.lower():
                pass

            if 'class=table-summary' in line:
                cur_section = None
                continue

            if cur_section is None:
                for name, pattern in sections:
                    if pattern in line:
                        cur_section = name
                        break
                if cur_section is None:
                    continue

            if cur_section == 'finish':
                break

            if line.find('class=column1') != -1:
                cur_metric_name = __get_td_tag_contents(line)

                if cur_section.endswith('(on target)'):
                    cur_metric_name += on_target_stats_suffix

                if cur_section.startswith('coverage'):
                    cur_metric_name = coverage_stats_prefix + cur_metric_name

                # if not metric_storage.get_metric(cur_metric_name):  # special case for Duplication rate and Clipped reads (Qualimap v.1 and v.2 difference)
                #     if metric_storage.get_metric(cur_metric_name + on_target_stats_suffix):  # extra 'on target' metrics
                #         cur_metric_name += on_target_stats_suffix

            if cur_metric_name and line.find('class=column2') != -1:
                __fill_record(cur_metric_name, line)
                cur_metric_name = None

    return value_by_metric
