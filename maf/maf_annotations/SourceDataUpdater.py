#! /opt/bin/python

from ftplib import FTP
import sys
import datetime
import re
import os
import gzip
from ftpFileUpdate import *
import shutil
from Bio import SeqIO, Entrez

class SourceDataUpdater():
    '''
    A group of functions to access various FTP sites and search for files
    newer than a specified date.
    '''
    
    def __init__(self, daysAgo, local_dir, alert_dir, update_all=False, gen_build=None):
        
        self.user = 'anonymous'
        self.passwd = 'byrne@cbio.mskcc.org'
        
        self.gen_build = gen_build
        self.daysAgo = daysAgo
        self.local_dir = local_dir
        self.alert_dir = alert_dir
        self.open_alert_files = []

        ## only build-specific human data is to be updated
        if gen_build == 'hg19' or gen_build == 'hg18':
            self.mask_alert = self.alert_dir+"/"+gen_build+"_mask_alert"
            self.mask_alert_handle = open(self.mask_alert, 'w')
            self.open_alert_files = [self.mask_alert_handle]
        
        ## only mouse data is to be updated
        elif gen_build == 'mm9':
            self.mouse_mapback_alert = self.alert_dir+"/mouse_mapback_alert"
            self.mouse_mapback_alert_handle = open(self.mouse_mapback_alert, 'w')
            self.open_alert_files = [self.mouse_mapback_alert_handle]

        ## only drosophila data is to be updated
        elif gen_build == 'dm3':
            self.drosophila_mapback_alert = self.alert_dir+"/drosophila_mapback_alert"
            self.drosophila_mapback_alert_handle = open(self.drosophila_mapback_alert, 'w')
            self.open_alert_files = [self.drosophila_mapback_alert_handle]
            
        ## universal human data is being updated
        elif update_all == True:
            self.sq_alert = self.alert_dir+"/sequenom_alert"
            self.mb_alert = self.alert_dir+"/mapback_alert"
            self.ensembl_alert = self.alert_dir+"/ensembl_alert"
            self.trembl_alert = self.alert_dir+"/trembl_alert"
            self.sprot_alert = self.alert_dir+"/swissprot_alert"
            self.mammal_alert = self.alert_dir+"/mammal_alert"
            self.other_alert = self.alert_dir+"/other_alert"
            self.human_alert = self.alert_dir+"/human_alert"
        
            self.sq_alert_handle = open(self.sq_alert,'w')
            self.mb_alert_handle = open(self.mb_alert,'w')
            self.ensembl_alert_handle = open(self.ensembl_alert,'w')
            self.trembl_alert_handle = open(self.trembl_alert,'w')
            self.sprot_alert_handle = open(self.sprot_alert,'w')
            self.mammal_alert_handle = open(self.mammal_alert,'w')
            self.other_alert_handle = open(self.other_alert,'w')
            self.human_alert_handle = open(self.human_alert,'w')

            self.open_alert_files = [self.sq_alert_handle, self.mb_alert_handle, self.ensembl_alert_handle, self.trembl_alert_handle,\
             self.sprot_alert_handle, self.mammal_alert_handle, self.human_alert_handle, self.other_alert_handle]

    def __del__(self):
        for oaf in self.open_alert_files:
            oaf.close()

    def save_alerts(self):
        for oaf in self.open_alert_files:
            oaf.close()


    ################################
    def update_local_file(self, host, remote_dir, filePattern):

        (now,month_nums) = getTodaysDate()
    
        # connect to host and get list of all files in remote directory
        (ftp,filesInRemDir) = getRemoteFileList(host, self.user, self.passwd, remote_dir)
    
        # use regular expressions to get a ls -al of files of interest
        filesToCheckDates = searchForFiles(filesInRemDir, filePattern)
    
        # determine which files have been modified within the last X days
        filesToDownload = listNewFiles(filesToCheckDates, month_nums, now, self.daysAgo)

        # download new files
        if len(filesToDownload)>0:
            downloadFiles(ftp, filesToDownload, remote_dir, self.local_dir)
            # unzip downloaded file
            for file in filesToDownload:
                if '.gz' in file:
                    unzipFile(file,self.local_dir) 
            ftp.close()
            # return names of downloaded files
            return filesToDownload 
    
        ftp.close()
        return False

    
    def attempt_genbank_download(self, local_directory, acc, uid):
        
        ## Request efetch() handle
        print "  Attempting to download " + acc
        print "  Requesting efetch() handle for UID:",uid
        handle2 = Entrez.efetch(db="nucleotide", id=uid, rettype="gbwithparts")
        
        ## Write record to local file
        file_name = acc+"_"+uid
        out_file = open(local_directory + "/" + file_name, 'w')
        print "  Opened and writing to file: ",file_name
        print >>out_file, handle2.read()
        out_file.close()
        
        ## Check for "//" in file (indicates that record has not been truncated)
        print "  Closed file: ",file_name
        print "  Checking completeness of file: ",file_name
        out_file = open(local_directory + "/" + file_name, 'rU')
        file_complete = False
        for line in out_file:
            if line[0:2] == "//":
                file_complete = True
        
        return file_complete
        
    
    
    ####### GENBANK CHROMOSOME FILES ######
    def update_chromosomes(self, build, chromosome_dir):
                
        if build == "hg18":
            accessions = ("NC_000001.9","NC_000002.10","NC_000003.10",\
                    "NC_000004.10","NC_000005.8","NC_000006.10",\
                    "NC_000007.12","NC_000008.9","NC_000009.10",\
                    "NC_000010.9","NC_000011.8","NC_000012.10",\
                    "NC_000013.9","NC_000014.7","NC_000015.8",\
                    "NC_000016.8","NC_000017.9","NC_000018.8",\
                    "NC_000019.8","NC_000020.9","NC_000021.7",\
                    "NC_000022.9","NC_000023.9","NC_000024.8",\
                    "NC_012920")
        
        elif build == "hg19":
            accessions = ("NC_000001","NC_000002","NC_000003",\
                    "NC_000004","NC_000005","NC_000006",\
                    "NC_000007","NC_000008","NC_000009",\
                    "NC_000010","NC_000011","NC_000012",\
                    "NC_000013","NC_000014","NC_000015",\
                    "NC_000016","NC_000017","NC_000018",\
                    "NC_000019","NC_000020","NC_000021",\
                    "NC_000022","NC_000023","NC_000024",\
                    "NC_012920")

        elif build == "mm9":
            accessions = ("NC_000067", "NC_000068","NC_000069",\
                    "NC_000070","NC_000071","NC_000072",\
					"NC_000073","NC_000074","NC_000075",\
                    "NC_000076","NC_000077","NC_000078",\
					"NC_000079","NC_000080","NC_000081",\
                    "NC_000082","NC_000083","NC_000084",\
                    "NC_000085","NC_000086","NC_000087",\
                    "NC_005089")	 
		
        elif build == "dm3":
            accessions = ("NT_033779.4","NW_001848855.1","NT_033778.3",\
                    "NW_001848856.1","NT_037436.3","NW_001848857.1",\
                    "NT_033777.2","NW_001848858.1","NC_004353.3",\
                    "NS_000188.1","NC_004354.3","NW_001848859.1",\
                    "NW_001848860.1","NC_001709.1")
                    
        failures = []
        Entrez.email = 'byrne@cbio.mskcc.org'
        print "Attempting to retrieve:",accessions
        for acc in accessions:
            print "PERFORMING ESEARCH() FOR ACCESSION: ", acc
            handle = Entrez.esearch(db="nucleotide", term=acc)
            record = Entrez.read(handle)
            if "IdList" in record:
                print "  Retrieved the following UIDs:", record["IdList"]
                for uid in record["IdList"]:
                    file_complete = False                    
                    while not file_complete:
                        file_complete = self.attempt_genbank_download(chromosome_dir, acc, uid)
                    print "  Accession",acc,"seems to be correctly terminated."
                    out_file.close()
            else:
                print "No UIDs retrieved for accession: ",acc
        print "\n"
        return

    ########## COSMIC ############
    def update_cosmic(self):

        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('cosmicmutantexport_v\d+_\d+',re.IGNORECASE)]
        host = 'ftp.sanger.ac.uk'
        remote_dir = 'pub/CGP/cosmic/data_export'

        # 2) CHECK FOR AND DOWNLOAD NEW FILES
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)

        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME  
        if files_downloaded:
            file = files_downloaded[0]
            copy = self.local_dir + "/newCosmic"
        
            try:
                shutil.copy(self.local_dir + "/" + file, copy)
                os.remove(self.local_dir + "/" + file)

                # 4) STORE ALERT IN APPROPRIATE FILE(S)
                msg = file +" was updated."
                self.add_msg_to_mb_alert(msg)
                self.add_msg_to_sq_alert(msg)

            except IOError:
                print >> sys.stderr, "  ERROR: COULD NOT RENAME "+file
        else:
            print >> sys.stderr, "  Current local COSMIC is up to date."

        return


    ######## REFSEQ HUMAN PROTEIN GPFF #########
    def update_refseq_human_protein_gpff(self):

        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('human\.protein\.gpff',re.IGNORECASE)]
        host = 'ftp.ncbi.nih.gov'
        remote_dir = 'refseq/H_sapiens/mRNA_Prot'

        # 2) CHECK FOR AND DOWNLOAD NEW FILES
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)

        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME  
        # NA

        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            file = files_downloaded[0]
            msg1 = "\t" + file + " was updated."
            self.add_msg_to_mb_alert(msg1)
            msg2 = "\tMapback blast database was updated (triggered by update of " + file + ")."
            self.add_msg_to_sprot_alert(msg2)
            self.add_msg_to_trembl_alert(msg2)
            self.add_msg_to_ensembl_alert(msg2)
            self.add_msg_to_human_alert(msg2)
            self.add_msg_to_mammal_alert(msg2)
            self.add_msg_to_other_alert(msg2)
        else:
            print >> sys.stderr, "  Current local refseq human protein file (gpff) is up to date."

        return


    ########## REFSEQ MOUSE PROTEIN GPFF ########
    def update_refseq_mouse_protein_gpff(self):
    
        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('mouse\.protein\.gpff',re.IGNORECASE)]
        host = 'ftp.ncbi.nih.gov'
        remote_dir = 'refseq/M_musculus/mRNA_Prot'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILES
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)
    
        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA
    
        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            file = files_downloaded[0]
            msg = "\t" + file + " was updated."
            self.add_msg_to_mouse_alert(msg)
        else:
            print >> sys.stderr, "  Current local refseq mouse protein file (gpff) is up to date."
        
        return

    ########### REFSEQ DROSOPHILA PROTEIN #########
    def update_refseq_drosophila_protein_gpff(self):
    
        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('invertebrate\d+\.protein\.gpff\.gz', re.IGNORECASE)]
        host = 'ftp.ncbi.nih.gov'
        remote_dir = 'refseq/release/invertebrate'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILES
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)
    
        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA
    
        # 4) STORE ALERT IN APPROPRIATE FILE(S) 
        if files_downloaded:
            file = files_downloaded[0]
            msg = "\t" + file + " was updated."
            self.add_msg_to_dm_alert(msg)
        else:
            print >> sys.stderr, "  Current local refseq drosophila protein file (gpff) is up to date."
        
        return



    ########## VERTEBRATE MAMMAL PROTEIN ########
    def update_mammal_protein_fasta(self):
    
        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('vertebrate_mammalian.*\.protein\.faa\.gz', re.IGNORECASE)]
        host = 'ftp.ncbi.nih.gov'
        remote_dir = 'refseq/release/vertebrate_mammalian'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)
    
        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA
    
        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            for file in files_downloaded:
                msg = "\t" + file + " was updated."
                self.add_msg_to_mammal_alert(msg)
        else:
            print >> sys.stderr, "  Current local vertebrate mammal fasta files are up to date."
        
        return


    ########## VERTEBRATE OTHER PROTEIN ########
    def update_other_protein_fasta(self):
    
        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('vertebrate_other.*\.protein\.faa\.gz',re.IGNORECASE)]
        host = 'ftp.ncbi.nih.gov'
        remote_dir = 'refseq/release/vertebrate_other'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)
    
        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA
    
        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            for file in files_downloaded:
                msg = "\t" + file + " was updated."
                self.add_msg_to_other_alert(msg)
        else:
            print >> sys.stderr, "  Current local vertebrate other fasta files are up to date."
        
        return


    ########### REFSEQ HUMAN PROTEIN FASTA #######
    def update_refseq_human_protein_fasta(self):

        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('human\.protein\.faa\.gz',re.IGNORECASE)]
        host = 'ftp.ncbi.nih.gov'
        remote_dir = 'refseq/H_sapiens/mRNA_Prot'

        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)

        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA
    
        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            file = files_downloaded[0]
            msg = "\t" + file + " was updated."
            self.add_msg_to_human_alert(msg)
        else:
            print >> sys.stderr, "  Current local refseq human protein file (fasta) is up to date."

        return


    ######### ENSEMBL HUMAN PROTEIN #######
    def update_ensembl_human_protein(self):
    
        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('Homo_sapiens\.GRCh37\.57\.pep\.all\.fa\.gz',re.IGNORECASE)]
        host = 'ftp.ensembl.org'
        remote_dir = 'pub/current_fasta/homo_sapiens/pep'

        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)
    
        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA
    
        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            file = files_downloaded[0]
            msg = "\t" + file + " was updated."
            self.add_msg_to_ensembl_alert(msg)
        else:
            print >> sys.stderr, "  Current local ensembl protein data is up to date."

        return


    ########## TREMBL HUMAN PROTEIN ##########
    def update_trembl_human_protein(self):
    
        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('uniprot_trembl_human\.dat\.gz',re.IGNORECASE)]
        host = 'ftp.expasy.org'
        remote_dir = 'databases/uniprot/current_release/knowledgebase/taxonomic_divisions'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)
    
        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA
    
        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            file = files_downloaded[0]
            msg = "\t" + file + " was updated."
            self.add_msg_to_trembl_alert(msg)
        else:
            print >> sys.stderr, "  Current local trembl protein data is up to date."

        return


    ########### SWISSPROT HUMAN PROTEIN #########
    def update_sprot_human_protein(self):
    
        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('uniprot_sprot_human\.dat\.gz',re.IGNORECASE)]
        host = 'ftp.expasy.org'
        remote_dir = 'databases/uniprot/current_release/knowledgebase/taxonomic_divisions'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)
    
        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA
    
        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            file = files_downloaded[0]
            msg = "\t" + file + " was updated."
            self.add_msg_to_sprot_alert(msg)
        else:
            print >> sys.stderr, "  Current local swissprot protein data is up to date."

        return


    ########### HG18 SNP ##########
    def update_hg18_snp(self):

        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('snp\d+\.txt\.gz',re.IGNORECASE)]
        host = 'hgdownload.cse.ucsc.edu'
        remote_dir = 'goldenPath/hg18/database'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)

        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA

        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            file = files_downloaded[0]
            msg = "\t" + file + " was updated."
            self.add_msg_to_mask_alert(msg)
        else:
            print >> sys.stderr, "  Current local hg18 SNP data is up to date."

        return


    ############ HG19 SNP ###########
    def update_hg19_snp(self):

        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('snp\d+\.txt\.gz',re.IGNORECASE)]
        host = 'hgdownload.cse.ucsc.edu'
        remote_dir = 'goldenPath/hg19/database'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)

        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA

        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            file = files_downloaded[0]
            msg = "\t" + file + " was updated."
            self.add_msg_to_mask_alert(msg)
        else:
            print >> sys.stderr, "  Current local hg19 SNP data is up to date."

        return


    ######### HG18 REPEAT MASKER ########
    def update_hg18_rmsk(self):

        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('chr\d+_rmskRM\d+\.txt\.gz', re.IGNORECASE)]
        host = 'hgdownload.cse.ucsc.edu'
        remote_dir = 'goldenPath/hg18/database'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)

        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA

        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            for file in files_downloaded:
                msg = "\t" + file + " was updated."
                self.add_msg_to_mask_alert(msg)
        else:
            print >> sys.stderr, "  Current local hg18 repeatMask data is up to date."

        return


    ########## HG19 REPEAT MASKER ########
    def update_hg19_rmsk(self):

        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('rmsk\.txt\.gz', re.IGNORECASE)]
        host = 'hgdownload.cse.ucsc.edu'
        remote_dir = 'goldenPath/hg19/database'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)

        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA

        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            for file in files_downloaded:   
                msg = "\t" + file + " was updated."
                self.add_msg_to_mask_alert(msg)
        else:
            print >> sys.stderr, "  Current local hg19 repeatMask data is up to date."

        return


    ########## HG18 SIMPLE REPEAT #########
    def update_hg18_simple_repeat(self):

        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('simpleRepeat\.txt\.gz', re.IGNORECASE)]
        host = 'hgdownload.cse.ucsc.edu'
        remote_dir = 'goldenPath/hg18/database'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)

        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA

        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            for file in files_downloaded:
                msg = "\t" + file + " was updated."
                self.add_msg_to_mask_alert(msg)
        else:
            print >> sys.stderr, "  Current local hg18 simple repeat data is up to date."

        return


    ########## HG19 SIMPLE REPEAT ########
    def update_hg19_simple_repeat():

        # 1) SET UP FTP INFO
        file_pattern_to_find = [re.compile('simpleRepeat\.txt\.gz', re.IGNORECASE)]
        host = 'hgdownload.cse.ucsc.edu'
        remote_dir = 'goldenPath/hg19/database'
    
        # 2) CHECK FOR AND DOWNLOAD NEW FILE(S)
        files_downloaded = self.update_local_file(host, remote_dir, file_pattern_to_find)

        # 3) RENAME FILE TO STANDARD CURRENT FILE NAME
        # NA

        # 4) STORE ALERT IN APPROPRIATE FILE(S)
        if files_downloaded:
            for file in files_downloaded:
                msg = "\t" + file + " was updated."
                self.add_msg_to_mask_alert(msg)
        else:
            print >> sys.stderr, "  Current local hg19 simple repeat data is up to date."

        return


    def add_msg_to_sq_alert(self, msg):
        print >> self.sq_alert_handle, msg
        return

    def add_msg_to_mb_alert(self, msg):
        print >> self.mb_alert_handle, msg
        return
    
    def add_msg_to_trembl_alert(self, msg):
        print >> self.trembl_alert_handle, msg
        return
    
    def add_msg_to_ensembl_alert(self, msg):
        print >> self.ensembl_alert_handle, msg
        return

    def add_msg_to_sprot_alert(self, msg):
        print >> self.sprot_alert_handle, msg
        return

    def add_msg_to_human_alert(self, msg):
        print >> self.human_alert_handle, msg
        return

    def add_msg_to_mammal_alert(self, msg):
        print >> self.mammal_alert_handle, msg
        return
    
    def add_msg_to_other_alert(self, msg):
        print >> self.other_alert_handle, msg
        return

    def add_msg_to_mask_alert(self, msg):
        print >> self.mask_alert_handle, msg
        return

    def add_msg_to_mouse_alert(self, msg):
        print >> self.mouse_mapback_alert_handle, msg
        return

    def add_msg_to_dm_alert(self, msg):
        print >> self.drosophila_mapback_alert_handle, msg
        return
