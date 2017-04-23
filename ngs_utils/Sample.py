class BaseSample:
    def __init__(self, name=None, dirpath=None, work_dir=None, bam=None, bed=None, vcf=None, genome=None,
                 targqc_dirpath=None, clinical_report_dirpath=None,
                 normal_match=None, sv_fpath=None, sv_bed=None,
                 l_fpath=None, r_fpath=None, **kwargs):
        self.name = name
        self.dirpath = dirpath
        self.work_dir = work_dir
        self.bam = bam
        self.l_fpath = l_fpath
        self.r_fpath = r_fpath
        self.is_wgs = False
        self.vcf = vcf
        self.phenotype = None
        self.gender = None
        self.genome = None
        self.var_dirpath = None
        self.normal_match = normal_match
        self.min_af = None
        self.sv_fpath = sv_fpath
        self.targqc_dirpath = targqc_dirpath
        self.clinical_html = None
        for k, v in kwargs.items():
            self.__dict__[k] = v

    def __cmp__(self, other):
        return self.key_to_sort().__cmp__(other.key_to_sort())

    def key_to_sort(self):
        parts = []

        cur_part = []
        prev_was_num = False

        for c in self.name:
            if prev_was_num == c.isdigit() and c not in ['-', '.']:  # same type of symbol, but not - or .
                cur_part.append(c)
            else:
                if cur_part:
                    part = ''.join(cur_part)
                    if prev_was_num:
                        part = int(part)
                    parts.append(part)
                    cur_part = []

                if c in ['-', '.']:
                    pass
                else:
                    if c.isdigit():
                        prev_was_num = True
                    else:
                        prev_was_num = False
                    cur_part.append(c)
        if cur_part:
            part = ''.join(cur_part)
            if prev_was_num:
                part = int(part)
            parts.append(part)

        return tuple(parts)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    @classmethod
    def load(cls, data):
        sample = cls(**data)
        sample.__dict__ = data
        return sample
