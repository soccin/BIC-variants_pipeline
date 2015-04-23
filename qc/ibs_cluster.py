#!/opt/bin/python

import sys

################################
## Cluster samples based on their likelihood of being paired
## as indicated by the mean and variance of their IBS scores
################################


################################
## Usage: %prog in_file max_adj_p_val
##
## in file format: tab-delimited file with 6 columns:
##
## (1) Sample1
## (2) Sample2
## (3) Total number of assays
## (4) mean IBS score
## (5) variance of IBS scores
##
################################

def create_clusters(infile,min_mn,max_var):

    clusters = []
    used_samples = []

    with open(infile,'rb') as f:
        for line in f:
            if 'Sample' not in line:
                s1,s2,total,mn,var,sd = line.strip().split("\t")
                if mn > min_mn and var < max_var:
                    ## check to see if either s1 or s2 has already been encountered
                    ## and therefore already belongs to a cluster.
                    ## if so, find that cluster and add the other sample to it
                    if s1 in used_samples or s2 in used_samples:
                        idxs = [] ## indexes of clusters to combine, if necessary
                        for i in range(len(clusters)):
                            if s1 in clusters[i] or s2 in clusters[i]:
                                idxs.append(i)
                        ## combine multiple clusters containing at least one of the two samples
                        ## and remove the individual clusters; worry about duplication later
                        for x in idxs[1:]:
                            clusters[idxs[0]] += clusters[x] 
                        if len(idxs)>1:
                            for idx in idxs[1:]:
                                del clusters[idx]
                        ## add the sample that was NOT found yet
                        if not s1 in clusters[idxs[0]]:
                            clusters[idxs[0]].append(s1)
                        elif not s2 in clusters[idxs[0]]:
                            clusters[idxs[0]].append(s2)
                    ## if neither sample belongs to a cluster, create a new one
                    else:
                        clusters.append([s1,s2])
                else:
                    if not s1 in used_samples:
                        clusters.append([s1])
                    if not s2 in used_samples:
                        clusters.append([s2])
                ## make sure both samples are in the list of ones already added
                ## to a cluster
                if not s1 in used_samples:
                    used_samples.append(s1)
                if not s2 in used_samples:
                    used_samples.append(s2)

    return clusters

def cluster_samples(args):
    """
    Cluster samples that meet the min mean and max variance thresholds set by user. 
    Sort clusters by number of elements in descending order and write
    to file the count of samples in each cluster followed by the sample names 
    """

    if not len(args) == 4:
        print >> sys.stderr, "Usage: %prog infile outfile min_mean max_var"
        sys.exit(-1)

    infile,outfile,min_mn,max_var = args
    clusters = create_clusters(infile,min_mn,max_var)
    clusters.sort(key=len,reverse=True)
    with open(outfile,'w') as out:
        for cl in clusters:
            print >> out, "%i\t%s" %(len(set(cl)),' '.join(set(cl)))
    return

if __name__ == '__main__':
    cluster_samples(sys.argv[1:])
