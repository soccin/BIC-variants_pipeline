## December 2012
## C.Byrne

## Download all publicly available MAFs from TCGA
## Create a manifest file containing local file name, URL from 
## which file was downloaded, date last modified and md5

import urllib2
import hashlib
import re
from HTMLParser import HTMLParser

class MyHTMLParser(HTMLParser):
    '''
    Custom HTML parser stores all files, directories, last modified timestamp and 
    file size listed on a TCGA HTML page
    '''
    def reset(self):
        HTMLParser.reset(self)
        self.curdir = ''
        self.subdirs = {}
        self.files = {}
        self.lastElement = ''
    def handle_data(self, data):
        data = data.strip()
        if 'Index of' in data:
            self.curdir = data.replace('Index of ','')
        elif '/' in data:
            self.subdirs[data] = []
        elif '.txt' in data or '.maf' in data:
            self.files[data] = []
        elif re.search('\d{2}-[A-Za-z]{3}-\d{4}', data):
            date = data.split()[0]
            time = data.split()[1]
            size = data.split()[2]
            if self.files.has_key(self.lastElement):
                self.files[self.lastElement] = [date,time,size]
            elif self.subdirs.has_key(self.lastElement):
                self.subdirs[self.lastElement] = [date,time,size]
        self.lastElement = data

def getHTML(url):
    '''
    Get HTML string from given url
    '''
    usock = urllib2.urlopen(url)
    html = usock.read().split("\n")
    usock.close()
    return html

def md5Checksum(filePath):
    '''
    Get md5 of given file
    '''
    fh = open(filePath, 'rb')
    m = hashlib.md5()
    while True:
        data = fh.read(8192)
        if not data:
            break
        m.update(data)
    return m.hexdigest()


def downloadFile(localdir,fileURL,name=None):
    '''
    Download remote file to local directory; in order to avoid
    overwriting a file with the same name, check list of names
    already used. If duplicates are found, name new file according
    to its remote parent directory AND file name
    '''

    global localmafs

    if not name:
        name = fileURL[fileURL.rfind('/')+1:]
        if name in localmafs: 
            name = fileURL[fileURL.rfind('/',0,fileURL.rfind('/'))+1:].replace('/','_')
    localmafs.append(name)        
    
    with open(localdir+name,'w') as out:
        usock = urllib2.urlopen(fileURL)
        ## remove trailing newline
        fl = usock.read().rsplit('\n',1)
        fl = ''.join(fl)
        print>>out,fl
        usock.close()

    return name
    

def getMAFlist(url):

    global mafList
    global manifests

    ## add trailing slash for consistency if need be
    if not url[-1] == "/":
        url += "/"

    ## extract & differentiate between files & directories listed
    ## on the current HTML page
    html = getHTML(url)
    for line in html:
        elements = line.split('/>W+<')
        for element in elements:
            p.feed(element)

    ## find any maf files among the list of files extracted
    for fl in p.files.keys():
        if '.maf' in fl:
            mafList.append([url+fl]+p.files[fl])  
            manifests.append(url+'MANIFEST.txt')
 
    ## drill down into any subdirectories and continue to look for mafs
    if len(p.subdirs)>0:
        for subdir in p.subdirs.keys():
            p.reset()
            newurl = url + subdir
            getMAFlist(newurl)            
    return 


localdir = '/ifs/data/BIC/tcga/maf/'
logfile = localdir + "MANIFEST.txt"
dom = "https://tcga-data.nci.nih.gov"
dir = "/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
url = dom + dir

mafList = []  ## a list of lists; each inner list will contain: 
              ## 1) the url from which the file will be downloaded
              ## 2) date of last modification
              ## 3) time of last modification
              ## 4) file size
manifests = []  ## a list of urls to MANIFEST.txt files which
                ## contain md5s for mafs 
localmafs = []  ## some maf files on tcga site have the same name, but
                ## live in different directories; in this case, the local
                ## maf is named by its parent directory and file name.
                ## this list stores all used file names

p = MyHTMLParser()

## populate mafList[]
getMAFlist(url)

with open(logfile,'w') as log:

    print>>log, '\t'.join(['Local File','Remote File Location','Last Modified','md5'])

    ## download each maf and check the md5
    for x in range(len(mafList)):
        maf,date,time,size = mafList[x]
        timestamp = ' '.join([date,time])

        ## download file and store local file name for reporting in log
        localFileName = downloadFile(localdir,maf)

        ## get md5 of local file
        md5 = md5Checksum(localdir+localFileName)
        man_md5 = ''
        try:
            man = getHTML(manifests[x])
            for line in man:
                ## maf file (as named by tcga, not this script) should
                ## be listed with its md5 in the associated MANIFEST.txt file
                if maf[maf.rfind('/')+1:] in line:
                    man_md5 = line.split()[0]
                    if man_md5 == md5:
                        print "SUCCESSFULLY DOWNLOADED: %S" % localFileName
                        print>>log, '\t'.join([localFileName,maf,timestamp,md5])
                    else:
                        print "WARNING: Possibly incomplete or corrupt file: %s" % localFileName
                        print>>log, '\t'.join([localFileName,maf,timestamp,md5+' (DOES NOT MATCH MD5 OF REMOTE FILE!!)'])
            if not man_md5 or len(man_md5)==0:
                print "\nlooking for md5 for %s" %maf
                print "no md5 found in %s" % manifests[x]
                print "\n"       
        except:
            man_md5 = "no manifest file found"
            print>>log, '\t'.join([localFileName,maf,timestamp,man_md5])

