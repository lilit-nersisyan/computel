OS.name = Sys.info()["sysname"]
if (!"seqinr" %in% installed.packages()) {
  
  #   if (tolower(OS.name) != "windows") {
  #     lib.paths = .libPaths()
  #     user.dir.ind = which(regexpr("/home/", lib.paths, fixed=T)==1)
  #     if(length(user.dir.ind) > 0 ) {
  #       install.lib = lib.paths[ind]
  #     } else {
  #       
  #     }
  #     
  #   }
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

# Build telomeric index                                                
#######################################################################################
#######################################################################################
cat("\nBuilding telomeric index in the directory", tel.index.dir, "\n")
success = build.tel.index(pattern, rl, min.seed=12, tel.index.dir, tel.index.prefix, 
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
  
  no.unal = ifelse(compute.base.cov,"", "--no-unal")
  
  align.args$additional.options = paste("-p", num.proc,quals, 
                                        no.unal, "--n-ceil", rl-min.seed, 
                                        "-D 20 -R 3 -N 1 -L", l, "-i S,1,0.5")
  align.args$mode = ifelse(mode.local, "--local", "--end-to-end")
  
  if(single){
    align.args$U = fastq
    reads = as.character(paste('-U', align.args$U))  
  } else {
    align.args$m1 = fastq1
    align.args$m2 = fastq2
    reads = as.character(paste('-1', align.args$m1,'-2', align.args$m2))
  }
  align.args$ignore.err = ignore.err
  success = do.call(bowtie.align, args=align.args)  
  if (!success){
    stop("Problem aligning reads to telomeric index\n")
  } else {
    
    #If compute.base.cov=T, separate mapped from unmapped reads
    ###########################################################
    if(!compute.base.cov){  
      reads.mapped = align.args$S     
    } else {
      cat("Separating mapped from unmapped reads\n")
      sam.file = align.args$S
      reads.mapped = file.path(align.dir, "tel.align.mapped.sam")  
      reads.unmapped = file.path(align.dir, "tel.align.unmapped.sam")  
      
      command.unmapped = paste(samtools.path, "view -Sh -f 4", sam.file,'>', reads.unmapped)
      command.mapped = paste(samtools.path, "view -Sh -F 4", sam.file,'>', reads.mapped)
      call.cmd(command.unmapped)
      call.cmd(command.mapped)
      
      cat("Converting unmapped reads to fastq with command\n")  
      samtofastq.command = paste("java -jar", picard.samtofastq.jar,
                                 "VALIDATION_STRINGENCY=LENIENT", 
                                 paste("INPUT=", reads.unmapped, sep=""))
      base.dir = file.path(output.dir, "base")
      dir.create(base.dir, showWarnings=F)
      file.prefix = paste(base.dir, "/reads.unmapped", sep="")
      if (!single){
        base.fastq1 = paste(file.prefix, "_1.fastq", sep="")
        base.fastq2 = paste(file.prefix, "_2.fastq", sep="")
        base.fastq.unpaired = paste(file.prefix, "_unpaired.fastq", sep="")
        
        base.fastq1.opt = paste("FASTQ=", base.fastq1, sep="") 
        base.fastq2.opt = paste("SECOND_END_FASTQ=", base.fastq2, sep="")  
        base.fastq.unpaired.opt = paste("UNPAIRED_FASTQ=", base.fastq.unpaired, sep="")
        samtofastq.command = paste(samtofastq.command, base.fastq1.opt, 
                                   base.fastq2.opt, base.fastq.unpaired.opt)    
      }
      else {
        base.fastq = paste(file.prefix, ".fastq", sep="")
        base.fastq.opt = paste("FASTQ=", base.fastq, sep="")  
        samtofastq.command = paste(samtofastq.command, base.fastq.opt, sep=" ")
      }
      print(samtofastq.command)
      call.cmd(samtofastq.command)
    }
    
    # Compute telomeric alignment coverage
    #######################################################################################
    #######################################################################################
    cat("Computing telomeric index coverage")
    tel.coverage = coverage.from.sam(samtools.path, sam.file= reads.mapped, bam.file.name="tel.align.bam", dir=align.dir)
    
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
    # Compute telomere length
    #######################################################################################
    #######################################################################################
    cat("estimating telomere length \n")
    tel.length.out = file.path(output.dir, "tel.length.xls")
    tel.length  = tel.length(tel.coverage, rl, pl= nchar(pattern), 
                             base.cov, num.haploid.chr, tel.length.out)
    
    cat("Telomere computation successfully complete.\n")
    cat("Output stored at ", tel.length.out, "\n")
    
  }
}