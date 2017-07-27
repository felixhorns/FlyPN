import sys
import fnmatch
import os

dir = sys.argv[1] # directory containing samples (*.fastq.gz)

# Find all *.fastq.gz files within directory
matches = []
for root, dirnames, filenames in os.walk(dir):
    for filename in fnmatch.filter(filenames, '*.fastq.gz'):
        if "_R1_" in filename:
            matches.append(root) # directory of match
            # matches.append(os.path.join(root, filename)) # full path to match
            
# Print directory containing each file
for x in sorted(matches):
    print x

