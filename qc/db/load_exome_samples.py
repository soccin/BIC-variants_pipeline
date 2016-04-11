#!/opt/common/python/python-2.7.6/bin/python
# load_samples.py Proj_*_title.txt
#
# This takes the Proj_*_title.txt file and loads it into MySQL.

import mysql.connector
from mysql.connector.cursor import MySQLCursorPrepared
import sys
import os.path
import mysql_helper
import db_config
import file_helper
import log_helper
import date_helper
import getopt

def load(project_name, rerun_number, pi, investigator, tumor_type, pipeline, revision_number, in_file_name, conn, log_file):
    project_cursor = conn.cursor(cursor_class=MySQLCursorPrepared)
    sample_cursor = conn.cursor(cursor_class=MySQLCursorPrepared)
    file_cursor = conn.cursor(cursor_class=MySQLCursorPrepared)
    project_file_cursor = conn.cursor(cursor_class=MySQLCursorPrepared)

    project_query = "INSERT INTO projects (name, pi, investigator, tumor_type, institution) VALUES (%s, %s, %s, %s, 'cmo') ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), tumor_type=%s, institution='cmo'"
    pipeline_run_query = "INSERT INTO pipeline_runs (project_id, rerun, pipeline, revision_number) VALUES (%s, %s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), revision_number=%s"
    # TODO probably should make the md5 the key - think about the timestamp though
    file_query = "INSERT INTO files (file_type, file_timestamp, location, name, md5) VALUES ('Title', %s, %s, %s, %s) ON DUPLICATE KEY UPDATE id=LAST_INSERT_ID(id), md5=%s"
    project_file_query = "INSERT INTO project_files (project_id, file_id) VALUES (%s, %s) ON DUPLICATE KEY UPDATE project_id=%s, file_id=%s"
    # should not use ON DUPLICATE KEY UPDATE because we have multiple unique indexes
    sample_query = "INSERT INTO samples (name, sample_type, pipeline_run_id) VALUES (%s, 'PROJECT', %s)"

    # save file info
    md5 = file_helper.get_md5(in_file_name)
    params = (file_helper.get_file_timestamp(in_file_name), os.path.dirname(in_file_name), os.path.basename(in_file_name), md5, md5)
    log_file.write("LOG: query = '%s', params = '%s'\n" % (file_query, params))
    file_cursor.execute(file_query, params)
    # if file_cursor.rowcount == 0, nothing changed
    # if file_cursor.rowcount == 1, insert happened 
    # if file_cursor.rowcount == 2, update happened - did md5 change?!
    file_cursor_rowcount = file_cursor.rowcount
    log_file.write("LOG: affected-rows value = %d\n" % (file_cursor_rowcount))
    file_id = mysql_helper.last_insert_id(file_cursor)

    #if file_cursor_rowcount == 1: 
    #    with open(in_file_name, "r") as in_file:
    
    params = (project_name, pi, investigator, tumor_type, tumor_type)
    log_file.write("LOG: query = '%s', params = '%s'\n" % (project_query, ", ".join(params)))
    project_cursor.execute(project_query, params)
    project_id = mysql_helper.last_insert_id(project_cursor)
    log_file.write("LOG: project_id = %s\n" % (project_id))
    log_file.write("LOG: project_cursor_rowcount = %s\n" % (project_cursor.rowcount))


    # now add pipeline_run
    params = (project_id, rerun_number, pipeline, revision_number, revision_number)
    project_cursor.execute(pipeline_run_query, params) 
    pipeline_run_id = mysql_helper.last_insert_id(project_cursor)
    log_file.write("LOG: project_cursor_rowcount = %s\n" % (project_cursor.rowcount))

    # save project file link
    params = (project_id, file_id, project_id, file_id)
    log_file.write("LOG: query = '%s', params = '%s'\n" % (project_file_query, params))
    project_file_cursor.execute(project_file_query, params)
 
    ## insert samples
    with open(in_file_name, "r") as in_file: 
        for line in in_file:
            sample_name = line.strip()
            params = (sample_name, pipeline_run_id)
            log_file.write("LOG: query = '%s', params = '%s'\n" % (sample_query, params))
            #try:
            sample_cursor.execute(sample_query, params)
            sample_id = mysql_helper.last_insert_id(sample_cursor)
            #except mysql.connector.Error as err:
            #    if err.errno == errorcode.ER_DUP_ENTRY:
            #   
            #    else:
            #        raise 
            #except:
            #    raise 
    #elif file_cursor_rowcount == 2:
        # this should never happen, but watch out for it
    #    log_file.write("WARNING: looks like md5 changed for file %s but file name and timestamp didn't\n" % (in_file_name))
    #elif file_cursor_rowcount == 0:
        # title file already in database?
    #    log_file.write("WARNING: could not insert title file %s. Does it already exist in database?\n" % (in_file_name))


    project_cursor.close()
    file_cursor.close()
    project_file_cursor.close()
    sample_cursor.close()

if __name__ == "__main__":
    if len(sys.argv) < 8:
        print "Usage: ./load_exome_samples.py [-o log_file.log] project_id rerun_number pi investigator tumor_type pipeline revision_number title_file"    
        print "\t./load_exome_samples.py -o Proj_4773_qc_db_load.log Proj_4773 0 faginj knaufj TUMOR_TYPE exome 2207 /ifs/solres/seq/faginj/knaufj/Proj_4773/Proj_4773_title.txt"
        sys.exit()

    log_file_name = "qc_db_scripts.%s.log" % (date_helper.now().strftime("%Y%m%d%H%M%S"))
    opts, args = getopt.getopt(sys.argv[1:],"o:",["out"])
   
    for opt, arg in opts:
        if opt in("-o","--out"):
            log_file_name = arg 
 

    project_id, rerun_number, pi, investigator, tumor_type, pipeline, revision_number, title_file = args

    conn = None
    with open(log_file_name, "w") as log_file:
        try:
            conn = mysql.connector.connect(**db_config.params)
    
            title_file = os.path.abspath(title_file)
            load(project_id, rerun_number, pi, investigator, tumor_type, pipeline, revision_number, title_file, conn, log_file)
    
            conn.commit()
            conn.close() 
        except Exception as e:
            if conn:
                conn.close() 
            log_helper.report_error_and_exit(log_file, str(e))
