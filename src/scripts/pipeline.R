# checking and installing dependencies

if (!"seqinr" %in% installed.packages()) {  
  install.packages("seqinr", repos = "http://cran.rstudio.com/")
}
if (!"psych" %in% installed.packages()) 
  install.packages("psych", repos = "http://cran.rstudio.com/")

library("seqinr")
source(file.path(scripts.dir, "functions.R"))

# Formatting variables
#######################################################################################
#######################################################################################

tel.index.prefix = "tel_index"
tel.index.fasta = "tel_index.fasta"
tel.index.dir = file.path(output.dir, "index")
dir.create(tel.index.dir, showWarnings=F)

base.index.prefix="base_index"
tel.genome.fasta = "telomeric.genome.fasta"

align.dir = file.path(output.dir, "align")
dir.create(align.dir, showWarnings=F)

if(compute.base.cov){
  base.dir = file.path(output.dir, "base_align")
  dir.create(base.dir, showWarnings=F)
  unmapped.fq = file.path(base.dir, "unmapped.fq.gz")
}

# Build telomeric index                                                
#######################################################################################
#######################################################################################
cat("\nBuilding telomeric index in the directory", tel.index.dir, "\n")
success = build.tel.index(pattern, rl, min.seed, tel.index.dir, tel.index.prefix, 
                          bowtie.build.path, ignore.err)
if (!success){
  stop("Problem building telomeric index")
} else {
  
  # Align reads to telomeric index
  #######################################################################################
  #######################################################################################
  align.args = list()
  align.args$bowtie.align.path = bowtie.align.path
  align.args$x = file.path(tel.index.dir, tel.index.prefix)
  align.args$S = file.path(align.dir, "tel.align.sam")
  
  l = floor(rl/3)
  if (l < 6)
    l = 6
  if (l > 22)
    l = 22
  
  no.unal = ifelse(compute.base.cov, paste("--no-unal --un-gz", unmapped.fq), "--no-unal")
  
  align.args$additional.options = paste("-p", num.proc,quals, 
                                        no.unal, "--n-ceil", rl-min.seed, 
                                        "-D 20 -R 3 -N 1 -L", l, "-i S,1,0.5")
  align.args$mode = ifelse(mode.local, "--local", "--end-to-end")
  
  if(single){
    align.args$U = paste(fastqs,collapse=",")
    reads = as.character(paste('-U', align.args$U))  
  } else {
    align.args$m1 = fastq1
    align.args$m2 = fastq2
    reads = as.character(paste('-1', align.args$m1,'-2', align.args$m2))
  }
  align.args$ignore.err = ignore.err
  if(compressed){
    align.args$file.compression = file.compression
    align.args$estimate.base.cov = estimate.base.cov    
    success = do.call(bowtie.align.compressed, args=align.args)     
  } else {
    success = do.call(bowtie.align, args=align.args)  
  }
  if(estimate.base.cov){
      length.file = file.path(output.dir, "temp")
      fq.length = as.numeric(read.table(length.file, header=F))
      base.cov = fq.length/4*rl/genome.length
      cat("\nesimated base coverage: ", base.cov, "\n")
  }
  if (!success){
    stop("Problem aligning reads to telomeric index\n")
  } else {
    

    #If compute.base.cov=T, specify unmapped reads
    ###########################################################

    reads.mapped = align.args$S
    if(compute.base.cov)
	base.fastq = unmapped.fq
    
    # Compute telomeric alignment coverage
    #######################################################################################
    #######################################################################################
    cat("Computing telomeric index coverage")
    tel.coverage = coverage.from.sam(samtools.path, sam.file= reads.mapped, bam.file.name="tel.align.bam", dir=align.dir)
    if(!file.exists(tel.coverage)){
    	stop("\n****************\nError: could not compute coverage: ", tel.coverage, "\n****************\n")
    }
    
    # Compute base coverage
    #######################################################################################
    #######################################################################################
    
    if(compute.base.cov){
      cat("Computing base coverage\n")
      
      #Calculate base coverage from alignment if compute.base.cov = T
      #################################################################
      cat("Aligning unmapped reads to reference\n")
      #Align unmapped reads to reference
      #######################################################################################
      base.align.args = list()
      base.align.args$bowtie.align.path = bowtie.align.path
      base.align.args$x = file.path(base.index.pathtoprefix)
      base.align.args$S = file.path(base.dir, "base.align.sam")  
      
      base.align.args$additional.options = paste("-p", num.proc, quals)
      base.align.args$mode = ifelse(mode.local, "--local", "--end-to-end")
      
      if(single){
        base.align.args$U = base.fastq
        reads = as.character(paste('-U', base.align.args$U))  
      } else {
        base.align.args$m1 = base.fastq1
        base.align.args$m2 = base.fastq2
        reads = as.character(paste('-1', base.align.args$m1,'-2', base.align.args$m2))
      }
      
      cat("Performing alginment with command: \n\n")
      success = do.call(bowtie.align, args=base.align.args)
      if (!success){
        stop("Problem aligning reads to base index")
      } else {
        
        cat("Computing base index coverage\n")
        coverage.file = coverage.from.sam(samtools.path, sam.file= base.align.args$S, bam.file.name="base.align.bam", dir=base.dir)
        
        base.cov = base.coverage(coverage.file)
      } 
    }



	############################################## 				Compute telomere length			##########################################


    cat("estimating telomere length \n")
    tel.length.out = file.path(output.dir, "tel.length.txt")
    tel.length  = get.tel.length(tel.coverage, fastqs,
								rl, pl= nchar(pattern),
								base.cov, num.haploid.chr, 
								genome.length, min.seed,
								tel.length.out)
    
	
	########################## 				Compute telomere variants			##########################

    cat("\n\ncomputing telomeric variant composition\n")
    tel.var.out = file.path(output.dir, "tel.variants.txt")
    tel.var  = count.variants(file = reads.mapped, pattern = pattern,  
                              tel.var.out, qual.threshold = qualt)
    
  
	########################## 				Remove tel.align.sam file			##########################

	for(f in list.files(output.dir, "*.sam", recursive = T, full.names = T)){
	  file.remove(f)
	}

	########################## 				Succcess			##########################
	
	
	cat("\n**** Success ****\n")
	cat("\nComputel successfully completed the calculations")
	cat("\nThe output is stored at ", tel.length.out, "\n")
	if(tel.var){
		cat("\nThe telomeric variant composition is written to ", tel.var.out, "\n")
	}

	
	
    
    if(is.nan(tel.length) || tel.length == 0) 
      cat("\n\nWarning:\tNo telomeric reads have been identified!")
    if(base.cov < 0.1)
      cat("\nWarning:\tThe estimated base coverage of ", base.cov, "is too low! The results may not be accurate for base coverages of < 0.1.\nIf this is from the Computel test run don't worry :) Just check if the telomere length you got is equal to 10683991\n\n")
	  
	  
  }
}
