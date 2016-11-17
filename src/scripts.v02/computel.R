#######################################################################################
######      Configuration options                                                 ######
#######################################################################################
###directories


scripts.dir="scripts"
bowtie.build.path="./bin/bowtie2-2.1.0-win/bowtie2-build.exe"
bowtie.align.path="./bin/bowtie2-2.1.0-win/bowtie2-align.exe"
samtools.path="./bin/samtools-win/samtools.exe"
picard.samtofastq.jar="/bin/SamToFastq.jar"

###input_reads
read.length=76
quals = "--phred33" ##default: --phred33, alternatives: --phred64, --solexa-quals

#treat reads as single-end (single set to T) or paired-end (single set to F)
single=T
# note: paired-end reads can effectively be treated as single-end, thus single is 'T' by default
# note: compressed files can only be supplied as single-end reads

## single-end reads
# if files.with.prefix is 'T', files should be specified by filename prefix;
# if it is 'F' - by absolute paths
files.with.prefix =  T

# single-end reads specified by filename prefix (files.with.prefix is set to T)
fastq.prefix="tel_reads"
fastq.dir="examples"

#single-end reads specified by file paths (files.with.prefix is set to F)
fastq="examples/tel_reads.fq,../examples/tel_reads1.fq     ../examples/tel_reads2.fq"

# file compression type (can be 'F' (not compressed), "gz" or "bz2")
# works only for single-end reads!

file.compression = 'F'

# if paired-end read alignment is prefered, 'single' should be set to 'F', and read pairs specified below:
# paired-end reads file paths (set single to F)
fastq1="examples/tel_reads1.fq"
fastq2="examples/tel_reads2.fq"

###base_coverage_calculation_options
compute.base.cov=F # the base coverage is given or should be computed

base.cov=5.4
base.index.pathtoprefix="examples/base.index/base_index"

estimate.base.cov =T  #if files are compressed and the base coverage needs to be estimated during unzipping, set "estimate.base.cov" to T
genome.length=3101804739

###output_options

output.dir='examples/output'

###algorithm_options

pattern='TTAGGG'
num.haploid.chr=23
min.seed=12
mode.local=F
###system_options

num.proc=3
ignore.err=F

###additional_options

quals="--phred33" #default: --phred33, alternatives: --phred64, --solexa-quals
ignore.err = T


################################################################
####     assemble options in config.table for validation    ####
################################################################

config.table = NULL
config.table['scripts.dir'] = scripts.dir
config.table['bowtie.build.path'] = bowtie.build.path
config.table['bowtie.align.path'] = bowtie.align.path
config.table['samtools.path'] = samtools.path
config.table['picard.samtofastq.jar'] = picard.samtofastq.jar

config.table['fastq1'] = fastq1
config.table['fastq2'] = fastq2
config.table['single'] = single
config.table['fastq'] = fastq
config.table['files.with.prefix'] = files.with.prefix
config.table['fastq.dir'] = fastq.dir
config.table['fastq.prefix']=fastq.prefix
config.table['file.compression']=file.compression
config.table['read.length'] = read.length

config.table['pattern'] = pattern
config.table['num.haploid.chr'] = num.haploid.chr
config.table['min.seed'] = min.seed
config.table['mode.local'] = mode.local

config.table['compute.base.cov'] = compute.base.cov
config.table['estimate.base.cov'] = estimate.base.cov
config.table['genome.length'] = genome.length
config.table['base.cov'] = base.cov
config.table['base.index.pathtoprefix'] = base.index.pathtoprefix

config.table['output.dir'] = output.dir
config.table['num.proc'] = num.proc
config.table['quals']=quals
config.table['ignore.err']=ignore.err

config.table = as.matrix(config.table)

validate.R = file.path(scripts.dir, "validate.options.R")
if (!file.exists(validate.R))
  validate.R = "validate.options.R"

if (!file.exists(validate.R)){
  stop("validate.options.R not found.\n Provide scripts.dir containing the script.")
} else {
  config.set = T
  source(validate.R)
}

dir.create(output.dir, showWarnings=F)

if (!config.set){
  stop("configuration not set successfully. Scripts will not execute.\n")
} else {
  source(pipeline.R)
}











