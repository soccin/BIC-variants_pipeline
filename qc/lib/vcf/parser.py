
import collections
import re
import csv
import gzip
import sys
import itertools

try:
    import pysam
except ImportError:
    pysam = None



# Metadata parsers/constants
RESERVED_INFO = {
    'AA': 'String', 'AC': 'Integer', 'AF': 'Float', 'AN': 'Integer',
    'BQ': 'Float', 'CIGAR': 'String', 'DB': 'Flag', 'DP': 'Integer',
    'END': 'Integer', 'H2': 'Flag', 'MQ': 'Float', 'MQ0': 'Integer',
    'NS': 'Integer', 'SB': 'String', 'SOMATIC': 'Flag', 'VALIDATED': 'Flag'
}

RESERVED_FORMAT = {
    'GT': 'String', 'DP': 'Integer', 'FT': 'String', 'GL': 'Float',
    'GQ': 'Float', 'HQ': 'Float'
}


_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])


class _vcf_metadata_parser(object):
    '''Parse the metadat in the header of a VCF file.'''
    def __init__(self):
        super(_vcf_metadata_parser, self).__init__()
        self.info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_pattern = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.format_pattern = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.meta_pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')

    def read_info(self, info_string):
        '''Read a meta-information INFO line.'''
        match = self.info_pattern.match(info_string)
        if not match:
            raise SyntaxError(
                "One of the INFO lines is malformed: %s" % info_string)

        num_str = match.group('number')
        try:
            if not num_str:
                num = None
            else:
                num = int(num_str)
                if num < 0:
                    sys.stderr.write("INFO[%s].Number = %s which is invalid.\n" \
                        % (match.group('id'), num_str))
        except ValueError:
            # TODO check and warn about other invalid values?
            num = num_str

        if num == 0 and match.group('type') != 'Flag':
            sys.stderr.write("INFO[%s].Number equals 0 but INFO[%s].Type is not Flag but %s\n" \
                % (match.group('id'), match.group('id'), match.group('type')))

        info = _Info(match.group('id'), num,
                     match.group('type'), match.group('desc'))

        return (match.group('id'), info)

    def read_filter(self, filter_string):
        '''Read a meta-information FILTER line.'''
        match = self.filter_pattern.match(filter_string)
        if not match:
            raise SyntaxError(
                "One of the FILTER lines is malformed: %s" % filter_string)

        filt = _Filter(match.group('id'), match.group('desc'))

        return (match.group('id'), filt)

    def read_format(self, format_string):
        '''Read a meta-information FORMAT line.'''
        match = self.format_pattern.match(format_string)
        if not match:
            raise SyntaxError(
                "One of the FORMAT lines is malformed: %s" % format_string)

        num_str = match.group('number')
        try:
            if not num_str:
                num = None
            else:
                num = int(num_str)
                if num <= 0:
                    sys.stderr.write("FORMAT[%s].Number = %s which is invalid.\n" \
                        % (match.group('id'), num_str))
        except ValueError:
            # TODO check and warn about other invalid values?
            num = num_str

        form = _Format(match.group('id'), num,
                       match.group('type'), match.group('desc'))

        return (match.group('id'), form)

    def read_meta(self, meta_string):
        match = self.meta_pattern.match(meta_string)
        return match.group('key'), match.group('val')


class _Call(object):
    """ A genotype call, a cell entry in a VCF file"""

    def __init__(self, site, sample, data):
        #: The ``_Record`` for this ``_Call``
        self.site = site
        #: The sample name
        self.sample = sample
        #: Dictionary of data from the VCF file
        self.data = data
        self.gt_nums = self.data['GT']
        #: True if the GT is not ./.
        self.called = self.gt_nums is not None

    def __repr__(self):
        return "Call(sample=%s, GT=%s, GQ=%s)" % (self.sample, self.gt_nums, self.data.get('GQ', ''))

    def __eq__(self, other):
        """ Two _Calls are equal if their _Records are equal
            and the samples and ``gt_type``s are the same
        """
        return (self.site == other.site
                and self.sample == other.sample
                and self.gt_type == other.gt_type)

    @property
    def gt_bases(self):
        '''The actual genotype alleles.
           E.g. if VCF genotype is 0/1, return A/G
        '''
        # nothing to do if no genotype call
        if self.called:
            # grab the numeric alleles of the gt string; tokenize by phasing
            phase_char = "/" if not self.phased else "|"
            (a1, a2) = self.gt_nums.split(phase_char)
            # lookup and return the actual DNA alleles
            try:
                return self.site.alleles[int(a1)] + \
                       phase_char + \
                       self.site.alleles[int(a2)]
            except:
                sys.stderr.write("Allele number not found in list of alleles\n")
        else:
            return None

    @property
    def gt_type(self):
        '''The type of genotype.
           hom_ref  = 0
           het      = 1
           hom_alt  = 2  (we don;t track _which+ ALT)
           uncalled = None
        '''
        # extract the numeric alleles of the gt string
        if self.called:
            # grab the numeric alleles of the gt string; tokenize by phasing
            (a1, a2) = self.gt_nums.split("/") \
                if not self.phased else self.gt_nums.split("|")
            if a1 == a2:
                if a1 == "0": return 0
                else: return 2
            else: return 1
        else: return None

    @property
    def phased(self):
        '''A boolean indicating whether or not
           the genotype is phased for this sample
        '''
        return self.data['GT'] is not None and self.data['GT'].find("|") >= 0

    def __getitem__(self, key):
        """ Lookup value, backwards compatibility """
        return self.data[key]

    @property
    def is_variant(self):
        """ Return True if not a reference call """
        if not self.called:
            return None
        return self.gt_type != 0

    @property
    def is_het(self):
        """ Return True for heterozygous calls """
        if not self.called:
            return None
        return self.gt_type == 1


class _Record(object):
    """ A set of calls at a site.  Equivalent to a row in a VCF file.

        The standard VCF fields CHROM, POS, ID, REF, ALT, QUAL, FILTER,
        INFO and FORMAT are available as properties.

        The list of genotype calls is in the ``samples`` property.
    """

    def __init__(self, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample_indexes, samples=None):
        self.CHROM = CHROM
        self.POS = POS
        self.ID = ID
        self.REF = REF
        self.ALT = ALT
        self.QUAL = QUAL
        self.FILTER = FILTER
        self.INFO = INFO
        self.FORMAT = FORMAT
        #: 0-based start coordinate
        self.start = self.POS - 1
        #: 1-based end coordinate
        self.end = self.start + len(self.REF)
        #: list of alleles. [0] = REF, [1:] = ALTS
        self.alleles = [self.REF]
        self.alleles.extend(self.ALT)
        #: list of ``_Calls`` for each sample ordered as in source VCF
        self.samples = samples
        self._sample_indexes = sample_indexes

    def __eq__(self, other):
        """ _Records are equal if they describe the same variant (same position, alleles) """
        return (self.CHROM == other.CHROM and
                self.POS == other.POS and
                self.REF == other.REF and
                self.ALT == other.ALT)

    def __iter__(self):
        return iter(self.samples)

    def __str__(self):
        return "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s)" % self.__dict__

    def __cmp__(self, other):
        return cmp( (self.CHROM, self.POS), (other.CHROM, other.POS))

    def add_format(self, fmt):
        self.FORMAT = self.FORMAT + ':' + fmt

    def add_filter(self, flt):
        if self.FILTER is None or self.FILTER == 'PASS':
            self.FILTER = ''
        else:
            self.FILTER = self.FILTER + ';'
        self.FILTER = self.FILTER + flt

    def add_info(self, info, value=True):
        self.INFO[info] = value

    def genotype(self, name):
        """ Lookup a ``_Call`` for the sample given in ``name`` """
        return self.samples[self._sample_indexes[name]]

    @property
    def num_called(self):
        """ The number of called samples"""
        return sum(s.called for s in self.samples)

    @property
    def call_rate(self):
        """ The fraction of genotypes that were actually called. """
        return float(self.num_called) / float(len(self.samples))

    @property
    def num_hom_ref(self):
        """ The number of homozygous for ref allele genotypes"""
        return len([s for s in self.samples if s.gt_type == 0])

    @property
    def num_hom_alt(self):
        """ The number of homozygous for alt allele genotypes"""
        return len([s for s in self.samples if s.gt_type == 2])

    @property
    def num_het(self):
        """ The number of heterozygous genotypes"""
        return len([s for s in self.samples if s.gt_type == 1])

    @property
    def num_unknown(self):
        """ The number of unknown genotypes"""
        return len([s for s in self.samples if s.gt_type is None])

    @property
    def aaf(self):
        """ The allele frequency of the alternate allele.
           NOTE 1: Punt if more than one alternate allele.
           NOTE 2: Denominator calc'ed from _called_ genotypes.
        """
        # skip if more than one alternate allele. assumes bi-allelic
        if len(self.ALT) > 1:
            return None
        hom_ref = self.num_hom_ref
        het = self.num_het
        hom_alt = self.num_hom_alt
        num_chroms = float(2.0*self.num_called)
        return float(het + 2*hom_alt)/float(num_chroms)

    @property
    def nucl_diversity(self):
        """
        pi_hat (estimation of nucleotide diversity) for the site.
        This metric can be summed across multiple sites to compute regional
        nucleotide diversity estimates.  For example, pi_hat for all variants
        in a given gene.

        Derived from:
        \"Population Genetics: A Concise Guide, 2nd ed., p.45\"
          John Gillespie.
        """
        # skip if more than one alternate allele. assumes bi-allelic
        if len(self.ALT) > 1:
            return None
        p = self.aaf
        q = 1.0-p
        num_chroms = float(2.0*self.num_called)
        return float(num_chroms/(num_chroms-1.0)) * (2.0 * p * q)

    def get_hom_refs(self):
        """ The list of hom ref genotypes"""
        return [s for s in self.samples if s.gt_type == 0]

    def get_hom_alts(self):
        """ The list of hom alt genotypes"""
        return [s for s in self.samples if s.gt_type == 2]

    def get_hets(self):
        """ The list of het genotypes"""
        return [s for s in self.samples if s.gt_type == 1]

    def get_unknowns(self):
        """ The list of unknown genotypes"""
        return [s for s in self.samples if s.gt_type is None]

    @property
    def is_snp(self):
        """ Return whether or not the variant is a SNP """
        if len(self.REF) > 1: return False
        for alt in self.ALT:
            if alt not in ['A', 'C', 'G', 'T']:
                return False
        return True

    @property
    def is_indel(self):
        """ Return whether or not the variant is an INDEL """
        if len(self.REF) > 1: return True
        for alt in self.ALT:
            if alt is None:
                return True
            elif len(alt) != len(self.REF):
                return True
        return False

    @property
    def is_transition(self):
        """ Return whether or not the SNP is a transition """
        # if multiple alts, it is unclear if we have a transition
        if len(self.ALT) > 1: return False

        if self.is_snp:
            # just one alt allele
            alt_allele = self.ALT[0]
            if ((self.REF == "A" and alt_allele == "G") or
                (self.REF == "G" and alt_allele == "A") or
                (self.REF == "C" and alt_allele == "T") or
                (self.REF == "T" and alt_allele == "C")):
                return True
            else: return False
        else: return False

    @property
    def is_deletion(self):
        """ Return whether or not the INDEL is a deletion """
        # if multiple alts, it is unclear if we have a transition
        if len(self.ALT) > 1: return False

        if self.is_indel:
            # just one alt allele
            alt_allele = self.ALT[0]
            if alt_allele is None:
                return True
            if len(self.REF) > len(alt_allele):
                return True
            else: return False
        else: return False

    @property
    def var_type(self):
        """
        Return the type of variant [snp, indel, unknown]
        TO DO: support SVs
        """
        if self.is_snp:
            return "snp"
        elif self.is_indel:
            return "indel"
        else:
            return "unknown"

    @property
    def var_subtype(self):
        """
        Return the subtype of variant [ts, tv, ins, del]
        TO DO: support SV sub_types
        """
        if self.is_snp:
            if self.is_transition:
                return "ts"
            elif len(self.ALT) == 1:
                return "tv"
            else: # multiple ALT alleles.  unclear
                return "unknown"
        elif self.is_indel:
            if self.is_deletion:
                return "del"
            elif len(self.ALT) == 1:
                return "ins"
            else: # multiple ALT alleles.  unclear
                return "unknown"
        else:
            return "unknown"

    @property
    def is_monomorphic(self):
        """ Return True for reference calls """
        return len(self.ALT) == 1 and self.ALT[0] is None

class Reader(object):
    """ Reader for a VCF v 4.0 file, an iterator returning ``_Record objects`` """


    def __init__(self, fsock=None, filename=None, compressed=False, prepend_chr=False):
        """ Create a new Reader for a VCF file.

            You must specify either fsock (stream) or filename.  Gzipped streams
            or files are attempted to be recogized by the file extension, or gzipped
            can be forced with ``compressed=True``
        """
        super(VCFReader, self).__init__()

        if not (fsock or filename):
            raise Exception('You must provide at least fsock or filename')

        if filename:
            self.filename = filename
            if fsock is None:
                self.reader = file(filename)

        if fsock:
            self.reader = fsock
            if filename is None:
                if hasattr(fsock, 'name'):
                    filename = fsock.name
            self.filename = filename

        if compressed or (filename and filename.endswith('.gz')):
            self.reader = gzip.GzipFile(fileobj=self.reader)

        #: metadata fields from header
        self.metadata = None
        #: remember original order of all header ##KEYS 
        #  use the OrderedDict keys for an ordered set
        self.metadata_id_order = collections.OrderedDict() 
        #: INFO fields from header
        self.infos = None
        #: FILTER fields from header
        self.filters = None
        #: FORMAT fields from header
        self.formats = None
        self.samples = None
        self._sample_indexes = None
        self._header_lines = []
        self._tabix = None
        self._prepend_chr = prepend_chr
        self._parse_metainfo()

    def __iter__(self):
        return self

    def _parse_metainfo(self):
        '''Parse the information stored in the metainfo of the VCF.

        The end user shouldn't have to use this.  She can access the metainfo
        directly with ``self.metadata``.'''
        for attr in ('metadata', 'infos', 'filters', 'formats'):
            setattr(self, attr, collections.OrderedDict())

        parser = _vcf_metadata_parser()

        line = self.reader.next()
        while line.startswith('##'):
            self._header_lines.append(line)
            line = line.strip()

            if line.startswith('##INFO'):
                key, val = parser.read_info(line)
                self.infos[key] = val
                self.metadata_id_order["INFO"] = None

            elif line.startswith('##FILTER'):
                key, val = parser.read_filter(line)
                self.filters[key] = val
                self.metadata_id_order["FILTER"] = None

            elif line.startswith('##FORMAT'):
                key, val = parser.read_format(line)
                self.formats[key] = val
                self.metadata_id_order["FORMAT"] = None

            else:
                key, val = parser.read_meta(line.strip())
                # don't make a list of values unless we have to
                if key not in self.metadata:
                    self.metadata[key] = val
                elif type(self.metadata[key]) != type([]):
                    self.metadata[key] = [self.metadata[key], val]
                else:
                    self.metadata[key].append(val)
                self.metadata_id_order[key] = None

            line = self.reader.next()

        fields = line.rstrip().split('\t')
        self.samples = fields[9:]
        self._sample_indexes = dict([(x,i) for (i,x) in enumerate(self.samples)])

    def _map(self, func, iterable, bad='.'):
        '''``map``, but make bad values None.'''
        return [func(x) if x != bad else None
                for x in iterable]

    def _parse_info(self, info_str):
        '''Parse the INFO field of a VCF entry into a dictionary of Python
        types. If INFO.Number == 0 and INFO.Type == "Flag" the value will
        be set to True, if INFO.Number != 1 then a list is returned (i.e. 
        INFO.Number == "A", "G", or ".") if INFO.Number == 1 and there is only
        one value the value is returned (not a list), otherwise a list is returned.  
        '''
        # TODO we might want to make lists for A and G of the correct length
        # with None filled in if we have missing data (record.ALT can be [None])

        if info_str == '.':
            return collections.OrderedDict() 
            
        entries = info_str.split(';')
        retdict = collections.OrderedDict() 

        for entry in entries:
            entry = entry.split('=')
            ID = entry[0]
            try:
                entry_type = self.infos[ID].type
                entry_num = self.infos[ID].num
            except KeyError:
                entry_num = None
                try:
                    entry_type = RESERVED_INFO[ID]
                except KeyError:
                    if entry[1:]:
                        entry_type = 'String'
                    else:
                        entry_type = 'Flag'

            if (entry_num is None and len(entry) == 1) \
                or (entry_num == 0 and entry_type == 'Flag'):
                val = True
            else:
                val = entry[1].split(',')
                if entry_type == 'Integer':
                    val = self._map(int, val)
                elif entry_type == 'Float':
                    val = self._map(float, val)

                # with check for len(val) == 1 we don't lose 
                # data even if it is in violation of our spec
                if entry_num == 1 and len(val) == 1:
                    val = val[0]
                # TODO else if entry_num == 1 and len(val) > 1 warning 
                # TODO other warnings about expected length and data types 
                    
            retdict[ID] = val

        return retdict

    def _parse_samples(self, samples, samp_fmt, site):
        '''Parse a sample entry according to the format specified in the FORMAT
        column.'''
        samp_data = []# OrderedDict()
        samp_fmt = samp_fmt.split(':')

        samp_fmt_types = []
        samp_fmt_nums = []

        for fmt in samp_fmt:
            try:
                entry_type = self.formats[fmt].type
                entry_num = self.formats[fmt].num
            except KeyError:
                entry_num = None
                try:
                    entry_type = RESERVED_FORMAT[fmt]
                except KeyError:
                    entry_type = 'String'
            samp_fmt_types.append(entry_type)
            samp_fmt_nums.append(entry_num)

        for name, sample in itertools.izip(self.samples, samples):
            sampdict = self._parse_sample(sample, samp_fmt, samp_fmt_types, samp_fmt_nums)
            samp_data.append(_Call(site, name, sampdict))

        return samp_data

    def _parse_sample(self, sample, samp_fmt, samp_fmt_types, samp_fmt_nums):
        sampdict = dict([(x, None) for x in samp_fmt])

        for fmt, entry_type, entry_num, vals in itertools.izip(
                samp_fmt, samp_fmt_types, samp_fmt_nums, sample.split(':')):

            # short circuit the most common
            if vals == '.' or vals == './.':
                sampdict[fmt] = None
                continue

            # we don't need to split single entries
            if entry_num == 1 or ',' not in vals:

                if entry_type == 'Integer':
                    sampdict[fmt] = int(vals)
                elif entry_type == 'Float':
                    sampdict[fmt] = float(vals)
                else:
                    sampdict[fmt] = vals

                if entry_num != 1:
                    sampdict[fmt] = (sampdict[fmt])

                continue


            vals = vals.split(',')

            if entry_type == 'Integer':
                sampdict[fmt] = self._map(int, vals)
            elif entry_type == 'Float' or entry_type == 'Numeric':
                sampdict[fmt] = self._map(float, vals)
            else:
                sampdict[fmt] = vals


        return sampdict


    def next(self):
        '''Return the next record in the file.'''
        row = self.reader.next().split()
        chrom = row[0]
        if self._prepend_chr:
            chrom = 'chr' + chrom
        pos = int(row[1])

        if row[2] != '.':
            ID = row[2]
        else:
            ID = None

        ref = row[3]
        alt = self._map(str, row[4].split(','))

        if row[5] == '.':
            qual = None
        else:
            qual = float(row[5]) if '.' in row[5] else int(row[5])
        filt = row[6].split(';') if ';' in row[6] else row[6]
        if filt == 'PASS':
            filt = None
        info = self._parse_info(row[7])

        try:
            fmt = row[8]
        except IndexError:
            fmt = None

        record = _Record(chrom, pos, ID, ref, alt, qual, filt, info, fmt, self._sample_indexes)

        if fmt is not None:
            samples = self._parse_samples(row[9:], fmt, record)
            record.samples = samples

        return record

    def fetch(self, chrom, start, end=None):
        """ fetch records from a Tabix indexed VCF, requires pysam
            if start and end are specified, return iterator over positions
            if end not specified, return individual ``_Call`` at start or None
        """
        if not pysam:
            raise Exception('pysam not available, try "pip install pysam"?')

        if not self.filename:
            raise Exception('Please provide a filename (or a "normal" fsock)')

        if not self._tabix:
            self._tabix = pysam.Tabixfile(self.filename)

        if self._prepend_chr and chrom[:3] == 'chr':
            chrom = chrom[3:]

        # not sure why tabix needs position -1
        start = start - 1

        if end is None:
            self.reader = self._tabix.fetch(chrom, start, start+1)
            try:
                return self.next()
            except StopIteration:
                return None

        self.reader = self._tabix.fetch(chrom, start, end)
        return self


class Writer(object):
    """ VCF Writer """

    fixed_fields = "#CHROM POS ID REF ALT QUAL FILTER INFO".split()

    def __init__(self, stream, template):
        self.writer = csv.writer(stream, delimiter="\t")
        self.template = template

        # remember that after a file has been read
        # the objects can be modifed or added to
        # so add whatever might be missing
        # and if it was already in metadata_id_order 
        # position will be unchanged
        # TODO define INFO, FILTER, FORMAT somewhere
        for header_key in ["INFO", "FILTER", "FORMAT"] + template.metadata.keys():
            template.metadata_id_order[header_key] = None

        # we use the keys to preserve the overall order 
        # that was in original input file
        # e.g. metadata, INFO, FILTER, FORMAT, more metadata...
        for key in template.metadata_id_order.keys():
          if key == "INFO":
            # output all infos at once (old and new)
            for line in template.infos.values():
              stream.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % tuple(self._map(str, line)))
          elif key == "FORMAT":
            # output all formats at once (old and new)
            for line in template.formats.values():
              stream.write('##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % tuple(self._map(str, line)))
          elif key == "FILTER":
            # output all filters at once (old and new)
            for line in template.filters.values():
              stream.write('##FILTER=<ID=%s,Description="%s">\n' % tuple(self._map(str, line)))
          else:
            # only output metadata with this key
            values = template.metadata[key]
            if type(values) == type([]):
              for line in template.metadata[key]:
                stream.write('##%s=%s\n' % (key, line))
            else:
              stream.write('##%s=%s\n' % (key, values))

        self._write_header()

    def _write_header(self):
        fields = list(self.fixed_fields) # copy fixed_fields
        if self.template.samples:
            # VCF spec: if genotype data is present in the file, 
            # [the 8 fixed, mandatory columns] are followed by a FORMAT column header
            fields += ["FORMAT"]
            fields += self.template.samples
        self.writer.writerow(fields)

    def write_record(self, record):
        """ write a record to the file """
        fields = self._map(str, [record.CHROM, record.POS, record.ID, record.REF]) \
              + [self._format_alt(record.ALT), record.QUAL or '.', self._format_filter(record.FILTER),
                 self._format_info(record.INFO)]

        if record.samples: 
            # if samples aren't defined in header for some reason 
            # the header will not have FORMAT column 
            fields += [record.FORMAT]
            fields += [self._format_sample(record.FORMAT, sample)
                for sample in record.samples]

        self.writer.writerow(fields)

    def _format_alt(self, alt):
        return ','.join([x or '.' for x in alt])

    def _format_filter(self, filter):
        if not filter:
            # when read in PASS was turned to None
            # which is different from "."
            # this is risky if the user has added data
            # and doesn't realize that
            return "PASS"
        elif type(filter) == type([]):
            return ";".join(filter)
        else:
            return filter

    def _format_info(self, info):
        if not info:
            return '.'
        info_strs = []
        for x, y in info.items():
            if x in self.template.infos and self.template.infos[x].type == 'Flag':
                info_strs.append(x)
            else:
                info_strs.append("%s=%s" % (x, self._stringify(y)))
        return ';'.join(info_strs)

    def _format_sample(self, fmt, sample):
        if sample.data["GT"] is None:
            return "./."
        return ':'.join(self._stringify(sample.data[f]) for f in fmt.split(':'))

    def _stringify(self, x, none='.'):
        if type(x) == type([]):
            return ','.join(self._map(str, x, none))
        return str(x) if x is not None else none

    def _map(self, func, iterable, none='.'):
        '''``map``, but make None values none.'''
        return [func(x) if x is not None else none
                for x in iterable]

def __update_readme():
    import sys, vcf
    file('README.rst', 'w').write(vcf.__doc__)


# backwards compatibility
VCFReader = Reader
VCFWriter = Writer
