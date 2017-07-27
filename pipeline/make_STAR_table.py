import sys
import fnmatch
import os

dir = sys.argv[1] # directory containing STAR log output files (Log.final.out)
outfile = sys.argv[2] # output file

# Find all Log.final.out files within directory
matches = []
for root, dirnames, filenames in os.walk(dir):
    for filename in fnmatch.filter(filenames, 'Log.final.out'):
        matches.append(root+"/"+filename) # full path of match
            
# Print directory containing each file
print len(matches), "files"

# Parse statistics from each file

samples = []
num_input_reads = []
num_uniquely_mapped_reads = []
num_reads_multiple_loci = []
num_reads_too_many_loci = []
percent_reads_unmapped_too_many_mismatches = []
percent_reads_unmapped_too_short = []

for match in matches:

    sample = match.split("/")[-3]
    samples.append(sample)

    with open(match) as f:

        for line in f:

            if "Number of input reads" in line:
                x = line.rstrip().split()[-1]
                num_input_reads.append(x)

            if "Uniquely mapped reads number" in line:
                x = line.rstrip().split()[-1]
                num_uniquely_mapped_reads.append(x)

            if "Number of reads mapped to multiple loci" in line:
                x = line.rstrip().split()[-1]
                num_reads_multiple_loci.append(x)

            if "Number of reads mapped to too many loci" in line:
                x = line.rstrip().split()[-1]
                num_reads_too_many_loci.append(x)

            if "% of reads unmapped: too many mismatches" in line:
                x = line.rstrip().split()[-1]
                percent_reads_unmapped_too_many_mismatches.append(x)

            if "% of reads unmapped: too short" in line:
                x = line.rstrip().split()[-1]
                percent_reads_unmapped_too_short.append(x)
                
# Print output

features = ["input",
            "uniquely_mapped",
            "multiple_loci",
            "too_many_loci",
            "percent_unmapped_too_many_mismatches",
            "percent_unmapped_too_short"]

stats = [num_input_reads, num_uniquely_mapped_reads, num_reads_multiple_loci,
         num_reads_too_many_loci, percent_reads_unmapped_too_many_mismatches, percent_reads_unmapped_too_short]

with open(outfile, 'w') as out:

    header = ["name"] + features # header
    out.write("\t".join(header) + "\n")
    
    for i in xrange(len(samples)):
        sample = samples[i]
        my_stats = [stats[x][i] for x in xrange(len(stats))]
        line = "\t".join([sample] + my_stats)
        out.write(line + "\n")
