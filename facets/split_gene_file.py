# on phos: /opt/common/CentOS_6/python/python-2.7.8/bin/python
#
# python split_gene_file.py -g /home/wilson/FACETS/Proj_4553_C/Pass2/Proj_4553_C___HISENS_GeneCalls.txt -d . 
#
# Writes to out directory data_CNA.txt, data_log2CNA.txt, and data_ASCNA.txt
# NOTE: pipeline will have to modify/add meta_*.txt files for each.

import getopt
import sys
import os.path
import csv
import re

tumor_to_normal_samples = {}

# just returns tumor
def find_tumor_normal_pair(combined_sample):
  for tumor, normal in tumor_to_normal_samples.iteritems():
    if combined_sample.endswith(tumor) and combined_sample.startswith(normal):
      return tumor 
  return None

def run(gene_filename, pairing_filename, output_dir, verbose):
  with open(pairing_filename, 'r') as pairing_file:
    pairing_reader = csv.reader(pairing_file, delimiter='\t')
    for row in pairing_reader:
      tumor_to_normal_samples[row[1]] = row[0]

  samples = {}
  # read in gene file and parse sample ids
  with open(gene_filename, 'r') as gene_file:
    gene_reader = csv.reader(gene_file, delimiter='\t')
    gene_header = gene_reader.next()
    #gene_expected_header = ["Tumor_Sample_Barcode", "Hugo_Symbol", "tcn", "lcn", "chr", "seg.start", "seg.end", "frac_elev_major_cn", "Nprobes", "WGD", "mcn", "FACETS_CNA", "FACETS_CALL"]
    gene_expected_header = ["Tumor_Sample_Barcode", "Hugo_Symbol", "tcn", "lcn", "cf", "tcn.em", "lcn.em", "cf.em", "chr", "seg.start", "seg.end", "frac_elev_major_cn", "Nprobes", "WGD", "mcn", "FACETS_CNA", "FACETS_CALL"]
  
    if gene_header != gene_expected_header:
      print >>sys.stderr, "ERROR: expected header '%s' in '%s'" % (",".join(gene_expected_header), gene_filename)
      sys.exit(2)
 
    all_genes = set([])
 
    for row in gene_reader:
      # WAS sample = row[0].split("__")[-2]
      # NOW /ifs/../s_PL_1_H3G1_S64_s_PL_2_H3G2_S65_hisens.cncf.txt
      # /ifs/../s_C_001117_N001_d1_s_C_001117_M001_d1_hisens.cncf.txt
      sample_pair = row[0].split("/")[-1].split(".")[0].replace("_hisens", "", 1).replace("_purity", "", 1)
      sample = find_tumor_normal_pair(sample_pair)
      if not sample:
        print >>sys.stderr, "ERROR: could not find sample pair '%s' in '%s'" % (sample_pair, pairing_filename)
        sys.exit(2)
      gene = row[1]
      all_genes.add(gene)
      if not sample in samples:
        samples[sample] = {}
      if gene in samples[sample]:
        print >>sys.stderr, "ERROR: '%s' '%s' is duplicated in '%s', not sure what to do" % (sample, gene, gene_filename)
        sys.exit(2)
      samples[sample][gene] = row

  cna_filename = os.path.join(output_dir, "data_CNA_facets.txt")
  #log2cna_filename = os.path.join(output_dir, "data_log2CNA_facets.txt")
  ascna_filename = os.path.join(output_dir, "data_ASCNA_facets.txt")
 
  sorted_sample_names = sorted(samples.keys())
 
  with open(cna_filename, 'w') as cna_file:
    #with open(log2cna_filename, 'w') as log2cna_file:
    with open(ascna_filename, 'w') as ascna_file:
      cna_writer = csv.writer(cna_file, delimiter='\t', lineterminator='\n')
      #log2cna_writer = csv.writer(log2cna_file, delimiter='\t', lineterminator='\n')
      ascna_writer = csv.writer(ascna_file, delimiter='\t', lineterminator='\n')

      # spit out headers
      header = [ "Hugo_Symbol" ] + [ "%s" % (sample) for sample in sorted_sample_names ]
      cna_writer.writerow(header)
      #log2cna_writer.writerow(header)
      ascna_writer.writerow(header)
 
      # spit out data
      # this assumes we have data!
      # WAS 0: Tumor_Sample_Barcode, 1: Hugo_Symbol, 2: tcn, 3: lcn, 4: chr, 5: seg.start, 6: seg.end, 7: frac_elev_major_cn, 8: Nprobes, 9: WGD, 10: mcn, 11: FACETS_CNA, 12: FACETS_CALL
      # NOW 0: Tumor_Sample_Barcode, 1: Hugo_Symbol, 2: tcn, 3: lcn, 4: cf, 5: tcn.em, 6: lcn.em, 7: cf.em, 8: chr, 9: seg.start, 10: seg.end, 11: frac_elev_major_cn, 12: Nprobes, 13: WGD, 14: mcn, 15: FACETS_CNA, 16: FACETS_CALL
      FACETS_CNA_INDX = 15
      WGD_INDX = 9
      MCN_INDX = 10
      LCN_INDX = 3
      FACETS_CALL_INDX = 16
      for gene in sorted(all_genes):
        cna_writer.writerow([ gene ] + [ samples[sample][gene][FACETS_CNA_INDX] if gene in samples[sample] else "NA" for sample in sorted_sample_names])
        #log2cna_writer.writerow([ gene ] + [ samples[sample][gene][13] for sample in sorted_sample_names])
        ascna_writer.writerow([ gene ] + [ "%s;%s;%s;%s;%s" % (samples[sample][gene][WGD_INDX], samples[sample][gene][MCN_INDX], samples[sample][gene][LCN_INDX], samples[sample][gene][FACETS_CNA_INDX], samples[sample][gene][FACETS_CALL_INDX]) if gene in samples[sample] else ["NA", "NA", "NA", "NA", "NA"] for sample in sorted_sample_names ])

def usage():
  print "python split_gene_file.py --verbose --gene Proj_[PROJECT_ID]___[TYPE]_GeneCalls.txt --pairing [PROJECT_ID]_sample_pairing.txt --dir DIR" 
  print "    e.g. /opt/common/CentOS_6/python/python-2.7.8/bin/python split_gene_file.py --gene /home/wilson/FACETS/Proj_4553_C/Pass2/Proj_4553_C___HISENS_GeneCalls.txt --pairing /ifs/solres/seq/solitd/solitd/Proj_93017_B/r_001/Proj_93017_B_sample_pairing.txt --dir /home/wilson/merges/bic-mskcc/test/mskcc/wilson/4553_c/"
  print "Where:"
  print "    -g gene filename, required"
  print "    -p sample pairing filename, required (tumor second column)"
  print "    -d output directory, required"
  print "    -v for verbose logging (optional)"
  print "Parses file and writes 3 portal gene matrix files: data_CNA.txt, data_log2CNA.txt, and data_ASCNA.txt"

def main():
  try:
    opts, args = getopt.getopt(sys.argv[1:], "g:p:d:vh", ["gene=", "pairing=", "dir=", "verbose", "help"])
  except getopt.GetoptError as err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

  verbose = False
  gene_filename = None
  pairing_filename = None
  output_dir = None
  for o, a in opts:
    if o in ("-v", "--verbose"):
      verbose = True
    elif o in ("-g", "--gene"):
      gene_filename = a
    elif o in ("-p", "--pairing"):
      pairing_filename = a
    elif o in ("-d", "--dir"):
      output_dir = a
    elif o in ("-h", "--help"):
      usage()
      sys.exit(0)  
    else:
      assert False, "unhandled option '%s'" % (o)

  if not gene_filename or not os.path.isfile(gene_filename):
    print >>sys.stderr, "ERROR: gene_filename '%s' is required" % (gene_filename)
    usage()
    sys.exit(2)

  if not pairing_filename or not os.path.isfile(pairing_filename):
    print >>sys.stderr, "ERROR: pairing_filename '%s' is required" % (pairing_filename)
    usage()
    sys.exit(2)

  if not output_dir or not os.path.isdir(output_dir):
    print >>sys.stderr, "ERROR: output_dir '%s' is required" % (output_dir)
    usage()
    sys.exit(2)

  run(gene_filename, pairing_filename, output_dir, verbose)

if __name__ == "__main__":
    main()
