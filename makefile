all:
	g++ TrimProblematicQuality.cpp -o TrimProblematicQuality
	g++ SyncFastq.cpp -o SyncFastq

clean:
	rm TrimProblematicQuality
	rm SyncFastq