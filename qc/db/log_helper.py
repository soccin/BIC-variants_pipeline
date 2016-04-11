import sys 

def report_error_and_exit(out_file, error):
    out_file.write("ERROR: %s" % (error, ))
    print >>sys.stderr, "ERROR: %s" % (error, )
    sys.exit(1)
