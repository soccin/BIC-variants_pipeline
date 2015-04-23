#!/opt/bin/python

from ftplib import FTP
import sys
import datetime
import re
import os
import gzip

def getTodaysDate():
	now = datetime.date.today()
	print>>sys.stderr, "Today's date:",now
	now = datetime.datetime(now.year, now.month, now.day)

	months=['', 'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
	month_nums={}
	x=0
	for month in months:
        	month_nums[month] = x
        	x+=1

	return (now,month_nums)

def getRemoteFileList(host,user,passwd,remote_dir):
	
    files = []
    ftp = None
    try:
        ftp=FTP(host,user,passwd)
        if not ftp==None:
            print >>sys.stderr,"Connected and Logged in to "+host
        else:
            print >>sys.stderr,"CONNECTION TO "+host+" FAILED"
            sys.exit[-1]
        ftp.dir(remote_dir, files.append)
        if not files==None:
            print >>sys.stderr, "  Got file list for "+remote_dir
        else:
            print >>sys.stderr, "  COULD NOT GET FILE LIST FOR "+remote_dir
            sys.exit[-1]
    except Exception, err:
        print>>sys.stderr, "  ERROR: %s\n" % (str(err))
        
    return (ftp, files) 

def searchForFiles(completeFileList,regExes):
	filesToCheckDate = []
	for file in completeFileList:
		for fileOfInterest in regExes:
			fileFound = fileOfInterest.search(file)
			if fileFound:
				filesToCheckDate.append(file)
				break
	return filesToCheckDate

def listNewFiles(filesToCheckDates,month_nums,now,daysAgo):
	newFiles = []
	for file in filesToCheckDates:
		info = file.split()
		fileMM = month_nums[info[5]]
		fileDD = int(info[6])
		if ":" in info[7]:
			fileYYYY = now.year
		else:
			fileYYYY = int(info[7])
		fileDate = datetime.datetime(fileYYYY,fileMM,fileDD)
		fileName = info[8]
		diff = now - fileDate
		if diff.days<daysAgo:
			newFiles.append(fileName)	
	return newFiles

def downloadFiles(ftp,filesToDownload,remote_dir,local_dir):
	os.chdir(local_dir)
	for fileName in filesToDownload:
		print>>sys.stderr, "  Downloading "+fileName+"..."
		fullPath = remote_dir+"/"+fileName
		local_file = fileName
		binaryFile = ftp.retrbinary("RETR %s" % fullPath, open(local_file,'w').write)
		if not binaryFile:
			ftp.retrlines("RETR " + filename, lambda s, w=local_file.write: w(s+"\n"))
		print>>sys.stderr, "  Download complete."
	return

def unzipFile(fileName,local_dir):
	os.chdir(local_dir)
	if ".gz" in fileName:
		#print >>sys.stderr,"\t\tUnzipping "+fileName+"..."
		try:
			f = gzip.open(fileName, 'rb')
			file_content = f.read()
		except OverflowError:
			print>>sys.stderr, "  ERROR: COULD NOT UNZIP "+fileName
			return
		unzipped_fileName = fileName.replace('.gz','')
		unzipped_handle = open(unzipped_fileName,'wb')
		print >>unzipped_handle, file_content
		unzipped_handle.close()
		print>>sys.stderr, "  Unzipped "+fileName
		os.remove(fileName)
	return

def main():
	daysAgo = int(sys.argv[1])
	(now,month_nums) = getTodaysDate()
	(ftp,filesInRemDir) = getRemoteFileList('hgdownload.cse.ucsc.edu','anonymous','byrne@cbio.mskcc.org','goldenPath/hg18/database')

	snpFiletoFind = re.compile('snp\d+\.txt\.gz', re.IGNORECASE)
        repeatFiletoFind = re.compile('(simpleRepeat.txt.gz)', re.IGNORECASE)
        rmskFiletoFind = re.compile('(rmsk\.txt\.gz)', re.IGNORECASE)

	#regExes = [snpFiletoFind,repeatFiletoFind,rmskFiletoFind]
	regExes = [repeatFiletoFind]
	filesToCheckDates = searchForFiles(filesInRemDir,regExes)
	filesToDownload = listNewFiles(filesToCheckDates,month_nums,now,daysAgo)
	if len(filesToDownload)>0:
		downloadFiles(ftp,filesToDownload,'goldenPath/hg18/database','/home/byrne/Mapback_Sequenom/HUMAN/hg18/SOURCE_DATA')
	else:
		print>>sys.stderr,"\t\tNO NEW FILES TO DOWNLOAD"
	ftp.close()

if __name__=='__main__':
	main()
