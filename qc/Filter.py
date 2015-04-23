#! /opt/bin/python

import vcf
import sys
import re
import argparse

class SNPFilter:
    """
    Filter VCF records. 

    (1) based on various parameters:
        (a) Percent called           - percentage of samples with a call 
        (b) Allele frequency         - frequency of ALT allele
        (c) Max type percent         - percentage of samples with the most frequent GT
        (d) Quality                  - QUAL

    AND/OR

    (2) based on a list of IDs given (id = [CHR]:[POS]; e.g. chr1:123456) 
    """

    def __init__(self, in_vcf, out_vcf):

        self.in_vcf = in_vcf
        self.out_vcf = out_vcf

        self.recs_by_pc = {}
        self.recs_by_af = {}
        self.recs_by_tp = {}
        self.recs_by_qual = {}

    def bin_records(self):
        """
        Bins records according to the following parameters:

        (1) Percent called           - percentage of samples with a call 
        (2) Allele frequency         - frequency of ALT allele
        (3) Max type percent         - percentage of samples with the most frequent GT
        (4) Quality                  - QUAL

        Returns nothing.
        """
        with open(self.in_vcf, 'rb') as vcf_file:
            vcf_reader = vcf.VCFReader(vcf_file)
            samples = vcf_reader.samples

            count = 0  ## count of records already examined

            for record in vcf_reader:
                count += 1
                ## keep only samples with one alternative and call rate greater than min percent called
                if len(record.ALT) == 1:

                    ## get genotypes for all samples
                    calls = []
                    for s in vcf_reader.samples:
                        if record.genotype(s).called:
                            calls = calls + [record.genotype(s)['GT']]

                    if len(set(calls)) > 1:
                        ## get percentage of samples with the most frequent GT
                        ## to make sure there is variety in hetero/homozygosity
                        max_type_percent = max(calls.count(a) for a in set(calls)) / float(len(calls))

                        id = record.CHROM + ":" + str(record.POS)

                        pc_key = round(record.call_rate,2)
                        if not pc_key in self.recs_by_pc:
                            self.recs_by_pc[pc_key] = []
                        self.recs_by_pc[pc_key].append(id)

                        af_key = round(record.INFO['AF'][0],2)
                        if not af_key in self.recs_by_af:
                            self.recs_by_af[af_key] = []
                        self.recs_by_af[af_key].append(id)

                        tp_key = round(max_type_percent,2)
                        if not tp_key in self.recs_by_tp:
                            self.recs_by_tp[tp_key] = []
                        self.recs_by_tp[tp_key].append(id)

                        qual_key = round(record.QUAL,2)
                        if not qual_key in self.recs_by_qual:
                            self.recs_by_qual[qual_key] = []
                        self.recs_by_qual[qual_key].append(id)

                if count % 10000 == 0:
                    print "Binned %d records." %count


        return

    def passing_records(self,min_pc,min_af,max_af,max_tp,min_qual=0):
        """
        Indexes record IDs by various parameters. Finds the records that meet all parameters set by finding the intersection
        of the records in the appropriate "bins".

        Returns the set of records in all "bins" that contain one or more records that meet
        all requirements.
        """

        recs_passing_pc = []
        recs_passing_af = []
        recs_passing_tp = []
        recs_passing_qual = []

        sorted_pc = sorted(self.recs_by_pc.keys(),reverse=True)
        for pc in sorted_pc:
            if pc < min_pc:
                break
            recs_passing_pc += self.recs_by_pc[pc]

        sorted_af = sorted(self.recs_by_af.keys())
        for af in sorted_af:
            if af >= min_af and af <= max_af:
                recs_passing_af += self.recs_by_af[af]

        sorted_tp = sorted(self.recs_by_tp.keys())
        for tp in sorted_tp:
            if tp > max_tp:
                break
            recs_passing_tp += self.recs_by_tp[tp]

        sorted_quals = sorted(self.recs_by_qual.keys(),reverse=True)
        for qual in sorted_quals:
            if qual < min_qual:
                break
            recs_passing_qual += self.recs_by_qual[qual]

        return set.intersection(set(recs_passing_pc),set(recs_passing_af),set(recs_passing_tp),set(recs_passing_qual))


    def write_vcf_records(self,ids=None,pick_every=1):
        """
        Filter VCF records based on a list of SNP IDs.
        """
        rec_count = 1
        with open(self.out_vcf,'w') as out:

            vcf_reader = vcf.VCFReader(open(self.in_vcf,'r'))
            vcf_writer = vcf.VCFWriter(out,vcf_reader)

            for rec in vcf_reader:
                id = rec.CHROM + ":" + str(rec.POS)
                if id in ids and rec_count % pick_every == 0:
                    vcf_writer.write_record(rec)
                rec_count += 1

        return

    def write_filtered_vcf(self,min_pc=0,min_af=0,max_af=1.0,max_tp=1.0,min_qual=0,ids=None,pick_every=1):
        """
        Filter VCF records based on parameters given. Optionally include a list of IDs to start with.
        """
        rec_count = 1
        with open(self.out_vcf,'w') as out:

            vcf_reader = vcf.VCFReader(open(self.in_vcf,'r'))
            vcf_writer = vcf.VCFWriter(out,vcf_reader)

            for rec in vcf_reader:
                id = rec.CHROM + ":" + str(rec.POS)
                if ids and not id in ids:
                    continue
                calls = []
                for s in vcf_reader.samples:
                    if rec.genotype(s).called:
                        calls += [rec.genotype(s)['GT']]
                max_type_percent = max(calls.count(a) for a in set(calls)) / float(len(calls))

                if rec.INFO['AF'][0] >= min_af and rec.INFO['AF'][0] <= max_af \
                  and max_type_percent <= max_tp \
                  and record.QUAL >= min_qual \
                  and rec_count % pick_every == 0:
                    vcf_writer.write_record(rec)
                rec_count += 1

        return

