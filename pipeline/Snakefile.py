include: "Snakefile_utils.py"

##### Paths
MY_HOME=                      '/local10G/rfhorns/FlyBrain/rnaseq/'
workdir:                      MY_HOME+'/log'

RESOURCES=                    '/local10G/rfhorns/resources'
ANACONDA=		      RESOURCES+'/anaconda2'
STAR=                         RESOURCES+'/STAR/bin/Linux_x86_64/STAR'
HTSEQ_COUNT=                  RESOURCES+'/HTSeq-0.7.1/HTSeq/scripts/count.py'

##### Parameters
REFS=                         MY_HOME+'/resources'
GENOMEDIR=                    REFS+'/STAR/dmel-all-r6.10-ERCC-Transgenes/'
REFERENCE_ANNOTATION=         REFS+'/STAR/dmel-all-r6.10-ERCC-Transgenes/dmel-all-r6.10-ERCC-Transgenes.gtf'

SEEDFILE=                     config['seedfile']  # seedfile corresponds to 1 directory per line
SCRATCH=                      config['scratch']  # LOCAL_SATA or LOCAL_SSD for /local10G, LOCAL_SCRATCH for sherlock

# Load samples
SEEDS = []
with open(SEEDFILE) as f:
    for line in f:
        SEEDS.append(line.strip())

##### Rules

rule all:
  input: expand("{dir}/done", dir=SEEDS)
  params: name='all', partition='general', mem='1024'

rule zcat_R1:
  """ Concatenate fastq files """
  input: get_all_fastq_gzs_R1
  output: temp('{dir}/R1.fastq')
  params: name='zcat', partition="general", mem="5300"
  shell: 'zcat {input} > {output[0]}'

rule zcat_R2:
  """ Concatenate fastq files """
  input: get_all_fastq_gzs_R2
  output: temp('{dir}/R2.fastq')
  params: name='zcat', partition="general", mem="5300"
  shell: 'zcat {input} > {output[0]}'
  
rule star:
  """ Map reads to genome using STAR  """
  input:  rules.zcat_R1.output, rules.zcat_R2.output
  output: '{dir}/STAR_output/Aligned.sortedByCoord.out.bam'
  params: name='star', partition='general', mem='64000'
  threads: 12
  run:
      wdir = os.path.dirname(str(output[0])) + '/'
      shell("{STAR} "
            "--genomeDir {GENOMEDIR} "
            "--readFilesIn {input[0]} {input[1]} "
            "--outSAMstrandField intronMotif "
            "--runThreadN {threads} "
            "--outFileNamePrefix {wdir} "
            "--outSAMtype BAM SortedByCoordinate "
            "--outSAMattributes NH HI AS NM MD "
            "--outReadsUnmapped Fastx "
            "--clip3pAdapterSeq CTGTCTCTTATACACATCT "
            "--outFilterType BySJout "
            "--outFilterMultimapNmax 20 "
            "--outFilterScoreMinOverLread 0.4 "
            "--outFilterMatchNminOverLread 0.4 "
            "--outFilterMismatchNmax 999 "
            "--outFilterMismatchNoverLmax 0.04 "
            "--alignIntronMin 20 "
            "--alignIntronMax 1000000 "
            "--alignMatesGapMax 1000000 "
            "--alignSJoverhangMin 8 "
            "--alignSJDBoverhangMin 1 ")

rule htseq:
  """ Count reads mapping to features using htseq """
  input:  rules.star.output
  output: '{dir}/htseq_output/htseq.tab'
  params: name='htseq', partition='general', mem='5300'
  shell: "source {ANACONDA}/bin/activate {ANACONDA} && "
         "python {HTSEQ_COUNT} -s no -r pos -f bam -m intersection-strict "
         "{input} {REFERENCE_ANNOTATION} > {output}"

rule clean:
  input: rules.htseq.output
  output: "{dir}/done"
  params: name='all', partition='general', mem='1024'
  shell: 'touch {output[0]}'

