#!/opt/common/python/python-2.7.6/bin/python
# load_project.py project_id project_dir
#
# e.g. ./load_project.py 4773 /ifs/res/seq/faginj/knaufj/Proj_4773/
#
# This takes the project_dir and loads the sample title file and stats files into MySQL.
# TODO should this be one big commit?
# This requires all stats files.
# TODO consider not running load data if we have the file stored already (based on md5).
# TODO consider reformatting files and using load data infile.

import sys
import os.path
import mysql.connector
import datetime
import db_config
import load_exome_samples
import load_exome_stats
import log_helper
import getopt

def get_title_file(project_name, project_dir):
    # try two places, but perhaps neither will be correct
    # if delivered directory then look here
    file_name = os.path.join(project_dir, "results/metrics/Proj_%s_sample_list.txt" % (project_name))
    if not os.path.exists(file_name):
        file_name = os.path.join(project_dir, "metrics/Proj_%s_sample_list.txt" % (project_name))
    return file_name

def get_stat_file(project_name, project_dir, stat):
    # pick the first file of the list of possible names that exists
    for file_suffix in stat.file_suffixes:
        file_name = os.path.join(project_dir, "results/metrics/Proj_%s%s" % (project_name, file_suffix)) 
        if os.path.exists(file_name):
            return file_name
        file_name = os.path.join(project_dir, "metrics/Proj_%s%s" % (project_name, file_suffix))
        if os.path.exists(file_name):
            return file_name
        file_name = os.path.join(project_dir, "results/metrics/fingerprint/Proj_%s%s" %(project_name, file_suffix))        
        if os.path.exists(file_name):
            return file_name
        else:
            file_name = os.path.join(project_dir, "metrics/fingerprint/Proj_%s%s" % (project_name, file_suffix))
            if os.path.exists(file_name):
                return file_name
    return None

def validate_files(project_name, project_dir, stats, log_file):
    if not os.path.exists(project_dir):
        log_helper.report_error_and_exit(log_file, "Project directory '%s' does not exist"  % (project_dir))

    title_file = get_title_file(project_name, project_dir)
    if not os.path.exists(title_file):
        log_helper.report_error_and_exit(log_file, "Sample title file '%s' does not exist"  % (title_file))
    for stat in stats:
        file = get_stat_file(project_name, project_dir, stat)
        if not file: # get_stat_file checks that it exists
            log_helper.report_error_and_exit(log_file, "%s '%s[%s]' does not exist"  % (stat.name, os.path.join(project_dir, "metrics/Proj_" + project_name), "|".join(stat.file_suffixes)))

    ## check for MarkDuplicates, CutAdapt file separately - TEMPORARILY
    ## TO DO: add these to load_stats.stats once database is updated
    if not os.path.exists(os.path.join(project_dir, "results/metrics/Proj_%s_markDuplicatesMetrics.txt" %(project_name))):
        if not os.path.exists(os.path.join(project_dir, "metrics/Proj_%s_markDuplicatesMetrics.txt" %(project_name))):
            log_helper.report_error_and_exit(log_file, "%s '%s[%s]' does not exist" %(MarkDuplicates, os.path.join(project_dir, "metrics/Proj_" + project_name), "_markDuplicatesMetrics.txt"))
    if not os.path.exists(os.path.join(project_dir, "results/metrics/Proj_%s_CutAdaptStats.txt" %(project_name))):
        if not os.path.exists(os.path.join(project_dir, "metrics/Proj_%s_CutAdaptStats.txt" %(project_name))):            
            log_helper.report_error_and_exit(log_file, "%s '%s[%s]' does not exist" %(CutAdaptStats, os.path.join(project_dir, "metrics/Proj_" + project_name), "_CutAdaptStats.txt")) 
    #####


def load(species, chipseq, paired, project_name, rerun_number, pi, investigator, tumor_type, pipeline, revision_number, project_dir, conn, log_file):
    stats = list()

    if species.lower() in ['human','hg18','hg19','hybrid','b37','grch37','xenograft']:
        stats = load_exome_stats.human_stats.values()
        if chipseq:
            stats = load_exome_stats.human_chipseq_stats.values()
    elif species.lower() in ['mouse','mm9','mm10']:
        stats = load_exome_stats.mouse_stats.values()
        if chipseq:
            stats = load_exome_stats.mouse_chipseq_stats.values()
    else:
        log_helper.report_error_and_exit(log_file, "species is not recognized: '%s'" % (species))

    ## If paired add the one PE dependent stat
    if paired:
        stats.extend(load_exome_stats.paired_end_stats.values())

    validate_files(project_name, project_dir, stats, log_file)
    # now load data 
    log_file.write("LOG: Loading title file\n")
    load_exome_samples.load(project_name, rerun_number, pi, investigator, 'unknown', pipeline, revision_number, get_title_file(project_name, project_dir), conn, log_file)
    for stat in stats:
        file_name = get_stat_file(project_name, project_dir, stat)
        log_file.write("LOG: Loading stat '%s' from '%s'\n" % (stat.name, file_name))
        load_exome_stats.load(species, chipseq, paired, stat.name, project_name, pi, investigator, rerun_number, revision_number, file_name, conn, log_file)


if __name__ == "__main__":
    #if len(sys.argv) < 8:
        #print "Usage: ./load_project.py [-o log_file.log] project_id rerun_number pi investigator tumor_type pipeline revision_number project_dir"    
        #print "\t./load_project.py -o Proj_4773_qc_db_load.log 4773 0 faginj knaufj THCA dmp 2207 /ifs/solres/seq/faginj/knaufj/Proj_4773/"
        #sys.exit()

    if not len(sys.argv) in [4,5,6,7,8]:  # with or without '-o log.txt'
        print "Usage: ./load_exome_project.py [-o log_file.log] [--chip] [--se] request_file revision_number project_dir"
        print "\t./load_exome_project.py -o Proj_4773_qc_db_load.log Proj_4773_request.txt 2207 /ifs/solres/seq/faginj/knaufj/Proj_4773/"
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

    #project_id, rerun_number, pi, investigator, tumor_type, pipeline, revision_number, project_dir = args
    request_file, revision_number, project_dir = args
    if not os.path.exists(request_file):
        print "ERROR: Can't find request file %s" %request_file
        sys.exit(-1)

    project_id, rerun_number, pi, investigator, species, tumor_type, pipeline = load_exome_stats.parse_request(request_file)

    if not project_id or not rerun_number or not pi or not investigator:
        print "ERROR: Required info missing from request file"
        sys.exit(-1)

    #'''
    print "log_file_name: ",log_file_name
    print "project_id: ",project_id
    print "rerun_number: ",rerun_number
    print "pi: ",pi
    print "investigator: ",investigator
    print "pipeline: ",pipeline
    print "revision_number: ",revision_number
    print "project_dir: ", project_dir 
    print "tumor_type: ",tumor_type
    print "species: ",species
    #'''
    #sys.exit()

    conn = None 
    with open(log_file_name, "w") as log_file:
        try:
            log_file.write("LOG: starting at %s\n" % (str(datetime.datetime.now())))
            conn = mysql.connector.connect(**db_config.params)

            project_dir = os.path.abspath(project_dir)
            load(species, chipseq, paired, project_id, rerun_number, pi, investigator, tumor_type, pipeline, revision_number, project_dir, conn, log_file)

            conn.commit()
            conn.close() 
            log_file.write("LOG: finished at %s\n" % (str(datetime.datetime.now())))
        except Exception as e:
            if conn:
                conn.close() 
            log_helper.report_error_and_exit(log_file, str(e))
