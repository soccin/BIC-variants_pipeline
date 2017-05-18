#!/opt/common/python/python-2.7.6/bin/python
# load_stats.py STAT_TYPE PROJECT_ID INPUT_FILE
# TODO passing the queries around etc. is silly

# need mysql.connector 2.0.1 or greater because of bug: https://bugs.mysql.com/bug.php?id=73370
# >>> print mysql.connector.__version__
# 2.1.3
import mysql.connector
from mysql.connector.cursor import MySQLCursorPrepared
import sys
import os.path
import re
import datetime
import mysql_helper
import db_config
import file_helper
import log_helper
import getopt
from collections import OrderedDict

class Stat:
    def __init__(self, name, file_suffixes, method):
        self.name = name
        self.file_suffixes = file_suffixes
        self.method = method

def find_sample_id_by_name_and_barcode(pipeline_run_id, conn, sample_cursor, sample_and_barcode_query, sample, barcode, log_file):
    params = (pipeline_run_id, sample, barcode)
    #log_file.write("LOG: query = '%s', params = '%s'\n" % (sample_and_barcode_query, params))
    sample_cursor.execute(sample_and_barcode_query, params)
    row = sample_cursor.fetchone()
    if row:
       return row[0]
    raise Exception("Unknown sample %s, barcode %s for this pipeline run (DB ID = %d). Has the title file been imported?" % (sample, barcode, pipeline_run_id))

def find_sample_id_by_barcode(pipeline_run_id, conn, sample_cursor, barcode_query, barcode, log_file):
    params = (pipeline_run_id, barcode)
    #log_file.write("LOG: query = '%s', params = '%s'\n" % (barcode_query, params))
    sample_cursor.execute(barcode_query, params)
    row = sample_cursor.fetchone()
    if row:
       return row[0]
    raise Exception("Unknown barcode %s for this pipeline run (DB ID = %d). Has the title file been imported?" % (barcode, pipeline_run_id))

def find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_name, log_file):
    params = (pipeline_run_id, sample_name)
    #log_file.write("LOG: query = '%s', params = '%s'\n" % (sample_query, params))
    sample_cursor.execute(sample_query, params)
    row = sample_cursor.fetchone()
    if row:
       return row[0]
    raise Exception("Unknown sample name %s for this pipeline run (DB ID = %d). Has the title file been imported?" % (sample_name, pipeline_run_id))

def read_and_save_hs_metrics(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file):
    # e.g. /ifs/res/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_HSmetrics.txt
    # not including columns 45, 46, and 47 (1 BASED INDEX) from file, they seem to be empty for all files
    # WARNING we are relying on the query columns to be in the order of the columns in the file

    ## DMP hs metrics file contains metrics from HS,MD,CutAdapt, and some custom columns.
    ## need to parse multiple files from exome pipeline to match what is in that file
    ## This is function is a temporary hack and will be replaced once the database is updated to include
    ## tables for each individual stat source


    stat_query = """INSERT INTO sample_hs_metrics 
        (sample_id, bait_set, genome_size, bait_territory, target_territory, 
            bait_design_efficiency, total_reads, pf_reads, pf_unique_reads, pct_pf_reads,
            pct_pf_uq_reads, pf_uq_reads_aligned, pct_pf_uq_reads_aligned, pf_uq_bases_aligned, on_bait_bases,
            near_bait_bases, off_bait_bases, on_target_bases, pct_selected_bases, pct_off_bait, 
            on_bait_vs_selected, mean_bait_coverage, mean_target_coverage, pct_usable_bases_on_bait, pct_usable_bases_on_target,
            fold_enrichment, zero_cvg_targets_pct, fold_80_base_penalty, pct_target_bases_2x, pct_target_bases_10x,
            pct_target_bases_20x, pct_target_bases_30x, pct_target_bases_40x, pct_target_bases_50x, pct_target_bases_100x,
            hs_library_size, hs_penalty_10x, hs_penalty_20x, hs_penalty_30x, hs_penalty_40x,
            hs_penalty_50x, hs_penalty_100x, at_dropout, gc_dropout, library, unpaired_reads_examined, read_pairs_examined,
            unmapped_reads, unpaired_read_duplicates, read_pair_duplicates, read_pair_optical_duplicates, percent_duplication,
            estimated_library_size, both_reads_align, one_read_aligns, neither_read_aligns, read_1_total_reads, 
            read_1_trimmed_reads, per_read_1_trimmed, read_2_total_reads, read_2_trimmed_reads, per_read_2_trimmed,
            file_id) 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), 
            bait_set=%s, genome_size=%s, bait_territory=%s, target_territory=%s, 
            bait_design_efficiency=%s, total_reads=%s, pf_reads=%s, pf_unique_reads=%s, pct_pf_reads=%s,
            pct_pf_uq_reads=%s, pf_uq_reads_aligned=%s, pct_pf_uq_reads_aligned=%s, pf_uq_bases_aligned=%s, on_bait_bases=%s,
            near_bait_bases=%s, off_bait_bases=%s, on_target_bases=%s, pct_selected_bases=%s, pct_off_bait=%s, 
            on_bait_vs_selected=%s, mean_bait_coverage=%s, mean_target_coverage=%s, pct_usable_bases_on_bait=%s, pct_usable_bases_on_target=%s,
            fold_enrichment=%s, zero_cvg_targets_pct=%s, fold_80_base_penalty=%s, pct_target_bases_2x=%s, pct_target_bases_10x=%s,
            pct_target_bases_20x=%s, pct_target_bases_30x=%s, pct_target_bases_40x=%s, pct_target_bases_50x=%s, pct_target_bases_100x=%s,
            hs_library_size=%s, hs_penalty_10x=%s, hs_penalty_20x=%s, hs_penalty_30x=%s, hs_penalty_40x=%s,
            hs_penalty_50x=%s, hs_penalty_100x=%s, at_dropout=%s, gc_dropout=%s, library=%s, unpaired_reads_examined=%s, read_pairs_examined=%s,
            unmapped_reads=%s, unpaired_read_duplicates=%s, read_pair_duplicates=%s, read_pair_optical_duplicates=%s, percent_duplication=%s,
            estimated_library_size=%s, both_reads_align=%s, one_read_aligns=%s, neither_read_aligns=%s, read_1_total_reads=%s, 
            read_1_trimmed_reads=%s, per_read_1_trimmed=%s, read_2_total_reads=%s, read_2_trimmed_reads=%s, per_read_2_trimmed=%s, file_id=%s"""

    hs_file = in_file_name
    md_file = hs_file.replace('_HsMetrics.txt','_markDuplicatesMetrics.txt')
    ca_file = hs_file.replace('_HsMetrics.txt','_CutAdaptStats.txt')

    md_by_sample = {}
    with open(md_file, "r") as in_file:
        header = in_file.readline().strip("\n")
        header_columns = [x.lower() for x in header.split("\t")]
        ## [library, unpaired_reads_examined, read_pairs_examined, unmapped_reads, unpaired_read_duplicates, read_pair_duplicates, read_pair_optical_duplicates, percent_duplication, estimated_library_size]
        for line in in_file:
            cells = line.strip("\n").split("\t")
            rec = OrderedDict(zip(header_columns, cells))
            if rec['library'].endswith("__1"):
                rec['library'] = rec['library'][:-3]
            rec['both_reads_align'] = rec['read_pairs_examined']
            rec['one_read_aligns'] = rec['unpaired_reads_examined']
            rec['neither_read_aligns'] = (float(rec['unmapped_reads']) - float(rec['unpaired_reads_examined']))/2 
            sample_id = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, rec['library'], log_file)
            md_by_sample[sample_id] = rec

    ca_by_sample = {}
    with open(ca_file, "r") as in_file:
        header = in_file.readline().strip("\n")
        header_columns = [x.lower() for x in header.split("\t")]
        ## [sample, r1_processedreads, r1_trimmedreads, r1_percenttrimmedreads, r2_processedreads, r2_trimmedreads, r2_percenttrimmedreads]
        header_columns[header_columns.index('r1_processedreads')] = 'read_1_total_reads'
        header_columns[header_columns.index('r2_processedreads')] = 'read_2_total_reads'
        header_columns[header_columns.index('r1_trimmedreads')] = 'read_1_trimmed_reads'
        header_columns[header_columns.index('r2_trimmedreads')] = 'read_2_trimmed_reads'
        header_columns[header_columns.index('r1_percenttrimmedreads')] = 'per_read_1_trimmed'
        header_columns[header_columns.index('r2_percenttrimmedreads')] = 'per_read_2_trimmed'
        
        for line in in_file:
            cells = line.strip("\n").split("\t")
            rec = OrderedDict(zip(header_columns, cells))
            sample_id = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, rec['sample'], log_file)
            ca_by_sample[sample_id] = rec

    hs_by_sample = {}
    with open(in_file_name, "r") as in_file:
        # read header
        header = in_file.readline().strip("\n")
        header_columns = [x.lower() for x in header.split("\t")]
        for line in in_file:
            cells = line.strip("\n").split("\t")
            sample_name = cells[header_columns.index("sample")]
            sample_id = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_name, log_file)
            hs_by_sample[sample_id] = OrderedDict(zip(header_columns,cells)) 

    if not sorted(set(md_by_sample.keys()+ca_by_sample.keys()+hs_by_sample.keys())) == sorted(set(hs_by_sample.keys())):
        log_file.write("ERROR: Sample lists are different between *HsMetrics.txt, *CutAdaptStats.txt and *markDuplicatesMetrics.txt")
        sys.exit(-1)

    fields = ["bait_set", "genome_size", "bait_territory", "target_territory", "bait_design_efficiency", "total_reads", "pf_reads", "pf_unique_reads", "pct_pf_reads", "pct_pf_uq_reads", "pf_uq_reads_aligned", "pct_pf_uq_reads_aligned", "pf_uq_bases_aligned", "on_bait_bases", "near_bait_bases", "off_bait_bases", "on_target_bases", "pct_selected_bases", "pct_off_bait", "on_bait_vs_selected", "mean_bait_coverage", "mean_target_coverage", "pct_usable_bases_on_bait", "pct_usable_bases_on_target", "fold_enrichment", "zero_cvg_targets_pct", "fold_80_base_penalty", "pct_target_bases_2x", "pct_target_bases_10x","pct_target_bases_20x", "pct_target_bases_30x", "pct_target_bases_40x", "pct_target_bases_50x", "pct_target_bases_100x", "hs_library_size", "hs_penalty_10x", "hs_penalty_20x", "hs_penalty_30x", "hs_penalty_40x", "hs_penalty_50x", "hs_penalty_100x", "at_dropout", "gc_dropout", "library", "unpaired_reads_examined", "read_pairs_examined", "unmapped_reads", "unpaired_read_duplicates", "read_pair_duplicates", "read_pair_optical_duplicates", "percent_duplication", "estimated_library_size", "both_reads_align", "one_read_aligns", "neither_read_aligns", "read_1_total_reads", "read_1_trimmed_reads", "per_read_1_trimmed", "read_2_total_reads", "read_2_trimmed_reads", "per_read_2_trimmed"]

    for sample_id in hs_by_sample.keys():
        all_stats = ca_by_sample[sample_id]
        all_stats.update(hs_by_sample[sample_id])
        all_stats.update(md_by_sample[sample_id])
        #params = [sample_id]
        params = []
        for field in fields:
            try:
                val = all_stats[field]
                if not val:
                    val = 0   #### is this correct?
                if val == '?':
                    val = None
                params.append(val)
            except KeyError:
                print "ERROR: '%s' field missing for sample with ID %s" %(field, sample_id)
                log_file.write("ERROR: '%s' field missing for sample with ID %s" %(field, sample_id)) 
        params.append(file_id)
        params = [sample_id] + params + params
        stat_cursor.execute(stat_query, params) 

def read_and_save_insert_size_metrics(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file):
    # e.g. /ifs/res/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_insertsizemetrics.txt
    stat_query = "INSERT INTO sample_insert_size_metrics (sample_id, insert_size, value, file_id) VALUES (%s, %s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), value=%s, file_id=%s"
    with open(in_file_name, "r") as in_file:
        # read header
        header = in_file.readline().strip("\n")
        sample_names = header.split("\t")[1:]
        sample_ids = []
        for sample_name in sample_names:
            sample_id = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_name, log_file)
            sample_ids.append(sample_id)

        for line in in_file:
            line = line.strip("\n")
            cells = line.split("\t")
            values = cells[1:]
            for i in range(0, len(values)):
                value = values[i]
                sample_id = sample_ids[i]
                params = (sample_id, cells[0], value, file_id, value, file_id)
                #log_file.write("LOG: query = '%s', params = '%s'\n" % (stat_query, params))
                stat_cursor.execute(stat_query, params)


def read_and_save_fpc_results_unmatch(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file):
    # e.g. /ifs/res/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_FPCResultsUnMatch.txt
    stat_query = "INSERT INTO sample_fpc_results_unmatch (reference_sample_id, query_sample_id, fraction_discordant_alleles, file_id) VALUES (%s, %s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), fraction_discordant_alleles=%s, file_id=%s"
    with open(in_file_name, "r") as in_file:
        for line in in_file:  
            line = line.strip("\n")
            if not line or len(line) == 0:
                log_file.write("WARNING: skipping file '%s'; No unexpected matches found.\n" %(in_file_name))
                return
            sample_1, sample_2, values = line.split("\t")
            sample_1_id = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_1, log_file)
            sample_2_id = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_2, log_file)
            value_1 = None
            value_2 = None
            if "," in values:
                value_1, value_2 = values.split(",")
            else:
                value_1 = values 
            params = (sample_1_id, sample_2_id, value_1, file_id, value_1, file_id)
            #log_file.write("LOG: query = '%s', params = '%s'\n" % (stat_query, params))
            stat_cursor.execute(stat_query, params)
            if value_2:
                params = (sample_2_id, sample_1_id, value_2, file_id, value_2, file_id)
                #log_file.write("LOG: query = '%s', params = '%s'\n" % (stat_query, params))
                stat_cursor.execute(stat_query, params)
        
def read_and_save_fpc_results_unmismatch(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file):
    # e.g. /ifs/res/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_FPCResultsUnMismatch.txt
    stat_query = "INSERT INTO sample_fpc_results_unmismatch (reference_sample_id, query_sample_id, fraction_discordant_alleles, file_id) VALUES (%s, %s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), fraction_discordant_alleles=%s, file_id=%s"
    with open(in_file_name, "r") as in_file:
        for line in in_file:  
            line = line.strip("\n")
            if not line or len(line) == 0:
                log_file.write("WARNING: skipping file '%s'; No unexpected matches found.\n" %(in_file_name))
                return
            sample_1, sample_2, values = line.split("\t")
            sample_1_id = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_1, log_file)
            sample_2_id = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_2, log_file)
            value_1 = None
            value_2 = None
            if "," in values:
                value_1, value_2 = values.split(",")
            else:
                value_1 = values 
            params = (sample_1_id, sample_2_id, value_1, file_id, value_1, file_id)
            #log_file.write("LOG: query = '%s', params = '%s'\n" % (stat_query, params))
            stat_cursor.execute(stat_query, params)
            if value_2:
                params = (sample_2_id, sample_1_id, value_2, file_id, value_2, file_id)
                #log_file.write("LOG: query = '%s', params = '%s'\n" % (stat_query, params))
                stat_cursor.execute(stat_query, params)

def read_and_save_fp_summary(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file):
    # e.g. /ifs/res/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_FPsummary.txt
    stat_query = "INSERT INTO sample_fp_summary (sample_id, chr, locus, allele_1, allele_2, allele_1_count, allele_2_count, genotype, minor_allele_freq, file_id) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), allele_1=%s, allele_2=%s, allele_1_count=%s, allele_2_count=%s, genotype=%s, minor_allele_freq=%s, file_id=%s"
    with open(in_file_name, "r") as in_file:
        # read header
        header = in_file.readline().strip("\n")
        sample_headers = header.split("\t")[1:]
        sample_names = [] # needed to preserve order, can't just use barcode_to_sample_ids keys
        sample_names_to_sample_ids = {}
        # each barcode has 3 columns
        if len(sample_headers) % 3 != 0:
            raise Exception("Expected 3 columns per barcode, but found %d barcode columns in %s" % (len(sample_headers), in_file_name))
        for i in xrange(0, len(sample_headers), 3):
            # TODO consider checking other 2 columns to make sure same barcode
            sample_name = re.sub('_Counts', '', sample_headers[i])
            sample_names.append(sample_name)
            sample_names_to_sample_ids[sample_name] = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_name, log_file)

        for line in in_file:  
            line = line.strip("\n")
            cells = line.split("\t")
            target = cells[0]
            if ":" not in target:
                log_file.write("WARNING: skipping locus '%s'\n" % (target))
            else: 
                chr, locus = target.split(":")
                chr = chr.replace('chr','')
                sample_values = cells[1:]
                for i in range(0, len(sample_names)):
                    counts = sample_values[i * 3]
                    genotype = sample_values[i * 3 + 1]
                    minor_allele_freq = sample_values[i * 3 + 2] 

                    # WARNING the following will break if there not exactly 2 alleles
                    if counts.count(' ') != 1:
                        log_file.write("WARNING: skipping line in '%s' because counts = '%s' and cannot split\n" % (in_file_name, counts))
                    else:
                        allele_1_count, allele_2_count = counts.split()
                        allele_1, allele_1_count = allele_1_count.split(":")
                        allele_2, allele_2_count = allele_2_count.split(":")

                        params = (sample_names_to_sample_ids[sample_names[i]], chr, locus, allele_1, allele_2, allele_1_count, allele_2_count, genotype, minor_allele_freq, file_id, allele_1, allele_2, allele_1_count, allele_2_count, genotype, minor_allele_freq, file_id)
                        #log_file.write("LOG: query = '%s', params = '%s'\n" % (stat_query, params))
                        stat_cursor.execute(stat_query, params)


def read_and_save_fp_avg_hom(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file):
    # e.g. /ifs/res/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_FPavgHom.txt
    stat_query = "INSERT INTO sample_fp_avg_hom (sample_id, value, file_id) VALUES (%s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), value=%s, file_id=%s"
    with open(in_file_name, "r") as in_file:
        # skip header
        header = in_file.readline().strip("\n")
        for line in in_file:  
            line = line.strip("\n")
            cells = line.split("\t")
            sample_name = cells[0]
            sample_id = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_name, log_file)
            value = cells[1]
            params = (sample_id, value, file_id, value, file_id)
            #log_file.write("LOG: query = '%s', params = '%s'\n" % (stat_query, params))
            stat_cursor.execute(stat_query, params)

def read_and_save_base_qualities(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file):
    # e.g. /ifs/res/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_basequalities.txt
    stat_query = "INSERT INTO sample_base_qualities (sample_id, cycle, value, file_id) VALUES (%s, %s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), value=%s, file_id=%s"
    with open(in_file_name, "r") as in_file:
        # read header
        header = in_file.readline().strip("\n")
        sample_names = header.split("\t")[1:]
        sample_ids = []
        for sample_name in sample_names:
            sample_ids.append(find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_name, log_file))

        for line in in_file:  
            line = line.strip("\n")
            cells = line.split("\t")
            cycle = cells[0]
            sample_values = cells[1:]
            for i, value in enumerate(sample_values):
                params = (sample_ids[i], cycle, value, file_id, value, file_id)
                #log_file.write("LOG: query = '%s', params = '%s'\n" % (stat_query, params))
                stat_cursor.execute(stat_query, params)

def read_and_save_org_base_qualities(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file):
    # e.g. /ifs/res/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_orgbasequalities.txt
    stat_query = "INSERT INTO sample_org_base_qualities (sample_id, cycle, value, file_id) VALUES (%s, %s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), value=%s, file_id=%s"
    with open(in_file_name, "r") as in_file:
        # read header
        header = in_file.readline().strip("\n")
        sample_names = header.split("\t")[1:]
        sample_ids = []
        for sample_name in sample_names:
            sample_ids.append(find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_name, log_file))

        for line in in_file:  
            line = line.strip("\n")
            cells = line.split("\t")
            cycle = cells[0]
            sample_values = cells[1:]
            for i, value in enumerate(sample_values):
                params = (sample_ids[i], cycle, value, file_id, value, file_id)
                #log_file.write("LOG: query = '%s', params = '%s'\n" % (stat_query, params))
                stat_cursor.execute(stat_query, params)

def read_and_save_fpc_summary(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file):
    # e.g. /ifs/res/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_FPCsummary.txt
    stat_query = "INSERT INTO sample_fpc_summary (sample_1_id, sample_2_id, value, file_id) VALUES (%s, %s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), value=%s, file_id=%s"
    with open(in_file_name, "r") as in_file:
        # read header
        header = in_file.readline().strip("\n")
        sample_names = header.split("\t")[1:]
        #log_file.write("LOG: barcodes = %s\n" % (", ".join(barcodes)))
        sample_names_to_ids = {}
        for sample_name in sample_names:
            params = (pipeline_run_id, sample_name)
            sample_names_to_ids[sample_name] = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_name, log_file)

        for line in in_file:  
            line = line.strip("\n")
            cells = line.split("\t")
            sample_name = cells[0]
            values = cells[1:]
            for i, value in enumerate(values):
                #log_file.write("LOG: barcode sample ids %s\n" % (", ".join([ "%s = %d" % (b, id) for b, id in barcodes_to_ids.iteritems()])))
                params = (sample_names_to_ids[sample_name], sample_names_to_ids[sample_names[i]], value, file_id, value, file_id)
                #log_file.write("LOG: query = '%s', params = '%s'\n" % (stat_query, params))
                stat_cursor.execute(stat_query, params)

def read_and_save_fp_het(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file):
    # e.g. /ifs/res/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_FPhet.txt
    stat_query = "INSERT INTO sample_fp_het (sample_id, value, file_id) VALUES (%s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), value=%s, file_id=%s"
    with open(in_file_name, "r") as in_file:
        # skip header
        header = in_file.readline().strip("\n")
        for line in in_file:  
            line = line.strip("\n")
            cells = line.split("\t")
            sample_name = cells[0]
            params = (pipeline_run_id, sample_name)
            sample_id = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_name, log_file)
            value = cells[1]
            params = (sample_id, value, file_id, value, file_id)
            #log_file.write("LOG: query = '%s', params = '%s'\n" % (stat_query, params))
            stat_cursor.execute(stat_query, params)

def read_and_save_gc_bias(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file):
    # e.g. /ifs/res/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_gcbias.txt
    stat_query = "INSERT INTO sample_gc_bias (sample_id, gc_content, value, file_id) VALUES (%s, %s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), value=%s, file_id=%s"
    with open(in_file_name, "r") as in_file:
        # read header
        header = in_file.readline().strip("\n")
        #barcodes = header.split("\t")[1:]
        samples = header.split("\t")[1:]
        sample_ids = []
        for sample_name in samples:
            sample_id = find_sample_id_by_name(pipeline_run_id, conn, sample_cursor, sample_query, sample_name, log_file)
            sample_ids.append(sample_id)

        for line in in_file:
            line = line.strip("\n")
            cells = line.split("\t")
            values = cells[1:]
            for i in range(0, len(values)):
                value = values[i]
                sample_id = sample_ids[i]
                params = (sample_id, cells[0], value, file_id, value, file_id)
                #log_file.write("LOG: i=%d, query = '%s', params = '%s'\n" % (i, stat_query, params))
                stat_cursor.execute(stat_query, params)


def load(species, chipseq, paired, stat_type, project_name, pi, investigator, rerun_number, revision_number, in_file_name, conn, log_file):
    stats=dict()
    if species.lower() in ['human','hg18','hg19','hybrid','b37','grch37','xenograft']:
        stats = human_stats
        if chipseq:
            stats = human_chipseq_stats
    elif species.lower() in ['mouse','mm9','mm10']:
        stats = mouse_stats
        if chipseq:
            stats = mouse_chipseq_stats

    ## If paired add the one PE dependent stat
    if paired:
        stats.update(paired_end_stats)

    if stat_type not in stats:
        raise Exception ("Do not know sample type '%s'" % (stat_type))

    project_cursor = conn.cursor(cursor_class=MySQLCursorPrepared)
    sample_cursor = conn.cursor(cursor_class=MySQLCursorPrepared)
    file_cursor = conn.cursor(cursor_class=MySQLCursorPrepared)
    pipeline_run_file_cursor = conn.cursor(cursor_class=MySQLCursorPrepared)
    stat_cursor = conn.cursor(cursor_class=MySQLCursorPrepared)

    # TODO add list to usage statement
    project_query = "SELECT id FROM projects WHERE name=%s AND pi=%s AND investigator=%s"
    pipeline_run_query = "SELECT id FROM pipeline_runs WHERE project_id=%s AND rerun=%s"

    # file md5 should not change if file_type, file_timestamp, location and name have not changed
    file_query = "INSERT INTO files (file_type, file_timestamp, location, name, md5) VALUES (%s, %s, %s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), md5=%s"
    # nothing to update if duplicate key
    pipeline_run_file_query = "INSERT INTO pipeline_run_files (pipeline_run_id, file_id) VALUES (%s, %s) ON DUPLICATE KEY UPDATE pipeline_run_id=%s, file_id=%s"
    barcode_query = "SELECT id FROM samples WHERE sample_type='PROJECT' AND pipeline_run_id=%s AND barcode=%s"
    sample_query = "SELECT id FROM samples WHERE sample_type='PROJECT' AND pipeline_run_id=%s AND name=%s"
    sample_and_barcode_query = "SELECT id FROM samples WHERE sample_type='PROJECT' AND pipeline_run_id=%s AND name=%s AND barcode=%s"

    params = (project_name, pi, investigator)
    #log_file.write("LOG: query = '%s', params = '%s'\n" % (project_query, ", ".join(params)))
    project_cursor.execute(project_query, params)
    row = project_cursor.fetchone()
    if row:
        project_id = row[0]
    else:
        raise Exception ("Project '%s' does not exist in database. Load title file first." % (project_name))

    # get pipeline_run_id
    params = (project_id, rerun_number)
    project_cursor.execute(pipeline_run_query, params)
    row = project_cursor.fetchone()
    if row:
        pipeline_run_id = row[0]
    else:
        raise Exception ("Pipeline run for project '%s', rerun '%s' does not exist in database. Load title file first." % (project_name, rerun_number))

    # save file info
    md5 = file_helper.get_md5(in_file_name)
    params = (stat_type, file_helper.get_file_timestamp(in_file_name), os.path.dirname(in_file_name), os.path.basename(in_file_name), md5, md5)
    #log_file.write("LOG: query = '%s', params = '%s'\n" % (file_query, params))
    file_cursor.execute(file_query, params)
    file_id = mysql_helper.last_insert_id(file_cursor)
    
    # save project file link
    params = (pipeline_run_id, file_id, pipeline_run_id, file_id)
    #log_file.write("LOG: query = '%s', params = '%s'\n" % (pipeline_run_file_query, params))
    pipeline_run_file_cursor.execute(pipeline_run_file_query, params)

    # call correct function for this stat_type
    stats[stat_type].method(file_id, in_file_name, project_id, pipeline_run_id, conn, sample_cursor, stat_cursor, barcode_query, sample_and_barcode_query, sample_query, log_file)

    project_cursor.close()
    file_cursor.close()
    pipeline_run_file_cursor.close()
    sample_cursor.close()
    stat_cursor.close()

def parse_request(request_file):
    project_id = rerun_number = pi = investigator = species = None
    tumor_type = 'unknown' ## default
    pipeline = 'exome'     ## default - might change to 'variants'??
    with open(request_file,'rU') as req:
        for line in req:
            pairs = line.strip('\n').split(': ')
            if not len(pairs) == 2:
                continue
            key,val = pairs
            if key == 'ProjectID':
                project_id = val.replace('Proj_','')
            elif key == 'RunNumber':
                rerun_number = val
            elif key == 'ORIG_PI':
                pi = val.replace('@mskcc.org','')
            elif not pi and key == 'PI':
                pi = val.replace('@mskcc.org','')
            elif key == 'ORIG_Investigator':
                investigator = val.replace('@mskcc.org','')
            elif not investigator and key == 'Investigator':
                investigator = val.replace('@mskcc.org','')
            elif key == 'TumorType':
                tumor_type = val
            elif key == 'Species':
                species = val
            elif key == 'Custom species':
                custom_species = val
            #elif key == 'Pipelines':   ## exclude this for now because values in request file do not match possible values in database
            #    pipeline = val     

    # This is to handle hybrids
    if species == 'custom':
        species = custom_species

    return (project_id,rerun_number,pi,investigator,species,tumor_type,pipeline)



'''
stats = {
    # _ALL_GenotypePooledNormal.txt, _ALL_targetcoverage.txt, _ALL_tilingcoverage.txt no longer delivered to investigator
    # but still can be found in /Results/ - TODO add back and test
    "FPHet": Stat("FPHet", ["_MajorContamination.txt"], read_and_save_fp_het),
    "FPCSummary": Stat("FPCSummary", ["_DiscordantHomAlleleFractions.txt"], read_and_save_fpc_summary),
    "BaseQualities": Stat("BaseQualities", ["_post_recal_MeanQualityByCycle.txt"], read_and_save_base_qualities),
    "FPAvgHom": Stat("FPAvgHom", ["_MinorContamination.txt"], read_and_save_fp_avg_hom),
    #"ExonCoverage": Stat("ExonCoverage", ["_ALL_exoncoverage.txt"], read_and_save_exon_coverage),
    #"CanonicalExonCoverage": Stat("CanonicalExonCoverage", ["_ALL_Canonical_exoncoverage.txt"], read_and_save_canonical_exon_coverage),
    "FPSummary": Stat("FPSummary", ["_FingerprintSummary.txt"], read_and_save_fp_summary),
    "FPCResultsUnMatch": Stat("FPCResultsUnMatch", ["_UnexpectedMatches.txt"], read_and_save_fpc_results_unmatch),
    "FPCResultsUnMismatch": Stat("FPCResultsUnMismatch", ["_UnexpectedMismatches.txt"], read_and_save_fpc_results_unmismatch),
    #"GeneCoverage": Stat("GeneCoverage", ["_ALL_genecoverage.txt"], read_and_save_gene_coverage),
    #"GCBias": Stat("GCBias", ["_ALL_gcbias.txt"], read_and_save_gc_bias),
    #"GenotypePooledNormal": Stat("GenotypePooledNormal", ["_ALL_GenotypePooledNormal.txt"], read_and_save_geneotype_pooled_normal),
    "HSMetrics": Stat("HSMetrics", ["_HsMetrics.txt"], read_and_save_hs_metrics),
    #"GenotypeHotspotNormals": Stat("GenotypeHotspotNormals", ["_ALL_genotypehotspotnormals.txt"], read_and_save_geneotype_hotspot_normals),
    "InsertSizeMetrics": Stat("InsertSizeMetrics", ["_InsertSizeMetrics_Histograms.txt"], read_and_save_insert_size_metrics),
    #"IntervalCoverageLoess": Stat("IntervalCoverageLoess", ["_ALL_intervalnomapqcoverage_loess.txt", "_ALL_intervalcoverage_loess.txt"], read_and_save_interval_coverage_loess),
    #"IntervalCoverage": Stat("IntervalCoverage", ["_ALL_intervalnomapqcoverage.txt", "_ALL_intervalcoverage.txt"], read_and_save_interval_coverage),
    "OrgBaseQualities": Stat("OrgBaseQualities", ["_pre_recal_MeanQualityByCycle.txt"], read_and_save_org_base_qualities),
    #"TargetCoverage": Stat("TargetCoverage", ["_ALL_exonnomapqcoverage.txt", "_ALL_targetcoverage.txt"], read_and_save_target_coverage),
    #"TilingCoverage": Stat("TilingCoverage", ["_ALL_tilingcoverage.txt"], read_and_save_tiling_coverage),
}
'''

# This stat will be extended to the species stats based on if fastq is paired or not
paired_end_stats = {
    "InsertSizeMetrics": Stat("InsertSizeMetrics", ["_InsertSizeMetrics_Histograms.txt"], read_and_save_insert_size_metrics)
}

human_stats = {
    "FPHet": Stat("FPHet", ["_MajorContamination.txt"], read_and_save_fp_het),
    "FPCSummary": Stat("FPCSummary", ["_DiscordantHomAlleleFractions.txt"], read_and_save_fpc_summary),
    "BaseQualities": Stat("BaseQualities", ["_post_recal_MeanQualityByCycle.txt"], read_and_save_base_qualities),
    "FPAvgHom": Stat("FPAvgHom", ["_MinorContamination.txt"], read_and_save_fp_avg_hom),
    "FPSummary": Stat("FPSummary", ["_FingerprintSummary.txt"], read_and_save_fp_summary),
    "FPCResultsUnMatch": Stat("FPCResultsUnMatch", ["_UnexpectedMatches.txt"], read_and_save_fpc_results_unmatch),
    "FPCResultsUnMismatch": Stat("FPCResultsUnMismatch", ["_UnexpectedMismatches.txt"], read_and_save_fpc_results_unmismatch),
    "HSMetrics": Stat("HSMetrics", ["_HsMetrics.txt"], read_and_save_hs_metrics),
    "OrgBaseQualities": Stat("OrgBaseQualities", ["_pre_recal_MeanQualityByCycle.txt"], read_and_save_org_base_qualities),
    "GCBias": Stat("GCBias", ["_GcBiasMetrics.txt"], read_and_save_gc_bias)
}

human_chipseq_stats = {
    "FPHet": Stat("FPHet", ["_MajorContamination.txt"], read_and_save_fp_het),
    "FPCSummary": Stat("FPCSummary", ["_DiscordantHomAlleleFractions.txt"], read_and_save_fpc_summary),
    "FPAvgHom": Stat("FPAvgHom", ["_MinorContamination.txt"], read_and_save_fp_avg_hom),
    "FPSummary": Stat("FPSummary", ["_FingerprintSummary.txt"], read_and_save_fp_summary),
    "FPCResultsUnMatch": Stat("FPCResultsUnMatch", ["_UnexpectedMatches.txt"], read_and_save_fpc_results_unmatch),
    "FPCResultsUnMismatch": Stat("FPCResultsUnMismatch", ["_UnexpectedMismatches.txt"], read_and_save_fpc_results_unmismatch),
    "GCBias": Stat("GCBias", ["_GcBiasMetrics.txt"], read_and_save_gc_bias)
}

mouse_stats = {
    "BaseQualities": Stat("BaseQualities", ["_post_recal_MeanQualityByCycle.txt"], read_and_save_base_qualities),
    "HSMetrics": Stat("HSMetrics", ["_HsMetrics.txt"], read_and_save_hs_metrics),
    "OrgBaseQualities": Stat("OrgBaseQualities", ["_pre_recal_MeanQualityByCycle.txt"], read_and_save_org_base_qualities),
    "GCBias": Stat("GCBias", ["_GcBiasMetrics.txt"], read_and_save_gc_bias)
}

mouse_chipseq_stats = {
    "GCBias": Stat("GCBias", ["_GcBiasMetrics.txt"], read_and_save_gc_bias)
}

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "Usage: ./load_stats.py [-o log_file.log] [--chip] [--se] stat_type request_file revision_number infile.txt"    
        print "\tWhere stat_type:"
        for stat in sorted(stats.keys()):
            print "\t\t-", stat 
        print "\te.g. ./load_stats.py -o Proj_4773_qc_db_load.log FPHet Proj_4773_request.txt 2207 /ifs/solres/seq/faginj/knaufj/Proj_4773/FinalReport/CompiledMetrics/Proj_4773_ALL_FPhet.txt"
        sys.exit()

    log_file_name = "qc_db_scripts.%s.log" % (datetime.datetime.now().strftime("%Y%m%d%H%M%S"))
    opts, args = getopt.getopt(sys.argv[1:],"o:",["out","chip","se"])
    chipseq = False
    paired = True

    for opt, arg in opts:
        if opt in("-o","--out"):
            log_file_name = arg
        elif opt == "--chip":
            chipseq = True
        elif opt == "--se":
             paired = False

    stat_type, request_file, revision_number, in_file_name = args

    if not os.path.exists(request_file):
        print "ERROR: Can't find request file %s" %request_file
        sys.exit(-1)

    project_id, rerun_number, pi, investigator, species, tumor_type, pipeline = parse_request(request_file)

    if not project_id or not rerun_number or not pi or not investigator:
        print "ERROR: Required info missing from request file"
        sys.exit(-1)


    conn = None
    with open(log_file_name, "w") as log_file:
        try:
            conn = mysql.connector.connect(**db_config.params)
    
            in_file_name = os.path.abspath(in_file_name)
            load(species, chipseq, paired, stat_type, project_id, pi, investigator, rerun_number, revision_number, in_file_name, conn, log_file)
    
            conn.commit()
            conn.close() 
        except Exception as e:
            if conn:
                conn.close() 
            #raise
            log_helper.report_error_and_exit(log_file, str(e))
