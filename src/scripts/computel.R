#######################################################################################
######      Configuration options                                                 ######
#######################################################################################
###directories


scripts.dir="scripts"
bowtie.build.path="../windows/bowtie2-2.1.0-win/bowtie2-build.exe"
bowtie.align.path="../windows/bowtie2-2.1.0-win/bowtie2-align.exe"
samtools.path="../windows/samtools-win/samtools.exe"
picard.samtofastq.jar="../SamToFastq.jar"

###input_reads

fastq1="../examples/tel_reads1.fq"
fastq2="../examples/tel_reads2.fq"
fastq="../examples/tel_reads.fq"
single=F
read.length=rl=76

###algorithm_options

pattern='TTAGGG'
num.haploid.chr=23
min.seed=12
mode.local=F

###base_coverage_calculation_options

compute.base.cov=T
base.cov=5.4
base.index.pathtoprefix="../examples/base.index/base_index"

###output_options

output.dir="../examples/output"

###system_options

num.proc=3

###additional_options

quals="--phred33" #default: --phred33, alternatives: --phred64, --solexa-quals
ignore.err=T


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
config.table['fastq'] = fastq
config.table['single'] = single
config.table['read.length'] = read.length

config.table['pattern'] = pattern
config.table['num.haploid.chr'] = num.haploid.chr
config.table['min.seed'] = min.seed
config.table['mode.local'] = mode.local

config.table['compute.base.cov'] = compute.base.cov
config.table['base.cov'] = base.cov
config.table['base.index.pathtoprefix'] = base.index.pathtoprefix

config.table['output.dir'] = output.dir
config.table['num.proc'] = num.proc
config.table['quals']=quals
config.table['ingore.err']=F

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











