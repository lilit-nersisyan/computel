
call.cmd <- function(command){
  OS.name = tolower(Sys.info()["sysname"])
  if (OS.name == "windows")
    shell(command, wait = T)
  else 
    system(paste("bash -c", "'", command, "'"), wait = T)
}

build.tel.index <- function(pattern, read.length, min.seed, 
                            tel.index.dir, tel.index.prefix, 
                            bowtie.build.path, ignore.err = F){
  pattern = string2vector(pattern)
  pattern = complement(pattern)
  pattern.length = read.length+length(pattern)-1;
  reps = floor(pattern.length/length(pattern))
  if(pattern.length-reps*length(pattern) > 0)
    left.rem = pattern[1:(pattern.length-reps*length(pattern))] else
      left.rem=NULL
  tel.index = c(rep(pattern,reps), left.rem)
  tel.index = c(tel.index, rep("N", read.length-min.seed)) 
  
  tel.index.fasta = paste(tel.index.prefix, ".fasta",sep="")  
  
  write.fasta(sequences=tel.index, names=c(">fwd: DNA"), 
              file=file.path(tel.index.dir, tel.index.fasta))
  this = getwd()
  setwd(tel.index.dir)
  command = paste(bowtie.build.path, '--quiet', tel.index.fasta,  tel.index.prefix)
  print(command)
  call.out = call.cmd(command)
  setwd(this)
  success = file.exists(paste(tel.index.dir,"/", tel.index.prefix, ".1.bt2",sep=""))
  if (success){
    if (!ignore.err){
      success = (call.out == 0)
    }
  }
  return(success)
  
}

complement <- function(pattern) {
  pattern = rev(toupper(pattern))
  letters =      c("A", "T", "G", "C", "W", "S", "M", "K", "R", "Y", "N");
  comp.letters = c("T", "A", "C", "G", "S", "W", "K", "M", "Y", "R", "N");
  comp.pattern = pattern;
  for (i in 1:length(letters)){
    ind = which(pattern == letters[i]);
    if(length(ind)>0)
      comp.pattern[ind] <-comp.letters[i];
  }
  return(comp.pattern);
}

tel.length <- function(coverage.file, rl, pl,
                       base.cov, num.haploid.chr, out.file){
  
  out.list = list()
  out.list$coverage.file = coverage.file
  out.list$read.length = rl
  out.list$pattern.length = pl
  out.list$base.cov = base.cov
  out.list$num.haploid.chr = num.haploid.chr
  out.list$tel.length = 0  
    
 
  
  
  table.tel = read.table(file=coverage.file,col.names=c("position","orientation", "depth"));   
  tl = length(table.tel$depth)
  range = c(1:(rl+pl-1))
  
  if(length(table.tel[,1]) < length(range))
    depth.tel=table.tel$depth/2/base.cov
  else
    depth.tel = table.tel$depth[range]/2/base.cov
  
  mean.length.rl.pl = mean(depth.tel)*(rl+pl-1)
  tel.length = mean.length.rl.pl/num.haploid.chr
  
  out.list$tel.length = tel.length
  out.table = as.matrix(out.list)
  
  cat("Computing telomere length from the followin parameters:\n")
  cat("coverage.file= ", coverage.file,"\n")
  cat("rl= ", rl,"\n")
  cat("pl= ", pl,"\n")
  cat("base.cov= ", base.cov,"\n")
  cat("num.haploid.chr= ", num.haploid.chr,"\n")
  cat("tel.length=",  tel.length, "\n")
  write.table(out.table, file = out.file, sep="\t")
  
  return(tel.length) 
}

base.coverage <- function(base.coverage.file){
  coverage.table = read.table(file=coverage.file,col.names=c("position","orientation", "depth"));   
  depth = coverage.table$depth
  mean.cov = mean(depth)
  return(mean.cov)
}

string2vector <- function(pattern){
  pattern = toupper(pattern)  
  if(length(pattern)==1) 
    #     pattern = substring(pattern, seq(1,nchar(pattern),1), seq(1,nchar(pattern),1))
    pattern = strsplit(pattern,"")[[1]];
  return(pattern)
}


bowtie.align <- function(bowtie.align.path, x, m1=NA, m2=NA, U=NA, S = "align.sam", 
                         mode="--end-to-end", 
                         additional.options = "", ignore.err = F) {
  args = ls();
  index = as.character(paste('-x', x))    
  
  if (!'U' %in% unlist(args) || is.na(U))
    reads = as.character(paste('-1', m1,'-2', m2))
  else  
    reads = as.character(paste('-U', U));
  
  
  samout = as.character(paste('-S', S))
  if (file.exists(S))
    file.remove(S)
    
  options = as.character(paste('-q', mode, '--quiet', additional.options))
  cat("\nPerforming alignment with command:\n")
  command = as.character(paste(bowtie.align.path, options,index, reads, samout));   
  cat(command)
  call.out = call.cmd(command)
  success = (file.exists(S) && file.info(S)$size > 0)
  if (success){
    if (!ignore.err){
      success = (call.out == 0)
    }
  }
  return(success)
}
#unzippes the files to stdout and pipes them to bowtie stdin as input. 
#command template: 
#cat fastq1 fastq2 ... | bunzip2(gunzip) - | tee >(wc -l > length) |
# bowtie2-align -p additional.options -x index -U - -S samout

bowtie.align.compressed <- function(bowtie.align.path, x, U=NA, S = "align.sam", 
                         mode="--end-to-end", 
                         additional.options = "", ignore.err = F, 
                         file.compression, estimate.base.cov) {
  args = ls();
  index = as.character(paste('-x', x))    
  
  
  reads = as.character(paste('-U', U));
  U.gz = gsub(",", replacement=" ", U)
  
  
  samout = as.character(paste('-S', S))
  if (file.exists(S))
    file.remove(S)
    
  options = as.character(paste('-q', mode, '--quiet', additional.options))
  if(file.compression == "gz")
    unzipper = "gunzip - "
  else 
    unzipper = "bunzip2 - "
  cat("\nPerforming alignment with command:\n")
  if(estimate.base.cov){    
    command = as.character(paste("cat", U.gz, "|", unzipper, "|", 
                                 "tee >(wc -l > length) |", bowtie.align.path, 
                                 options, index, "-U - ", samout))
    
    cat(command)
    call.out = call.cmd(command)
    success = (file.exists(S) && file.info(S)$size > 0)
    if (success){
      if (!ignore.err){
        success = (call.out == 0)
      }
    }
    if(success){
      success = (file.exists("length") && file.info("length")$size > 0)
      if(!success){
        cat("could not find file length for base coverage estimation \n")
      }
    }
  } else {
    command = as.character(paste("cat", U.gz, "|", unzipper, "|", 
                                 bowtie.align.path, 
                                 options, index, "-U - ", samout))
    
    cat(command)
    call.out = call.cmd(command)
    success = (file.exists(S) && file.info(S)$size > 0)
    if (success){
      if (!ignore.err){
        success = (call.out == 0)
      }
    }
  }
  return(success)
}

coverage.from.sam <- function(samtools.path, sam.file, bam.file.name, dir){
  #Convert sam to bam
  ###################
  #bam.file = file.path(dir, "tel.align.bam")
  bam.file = file.path(dir, bam.file.name)  
  samtobam = paste(samtools.path, "view -bSh", sam.file, ">", bam.file)
  cat("\nConverting sam file", sam.file,"to bam file", bam.file, "...\n")
  call.cmd(samtobam)
  
  #Sort bam
  ###################
  sorted.file = file.path(dir, paste(bam.file.name, "_sorted",sep=""))
  sort = paste(samtools.path, "sort", bam.file, sorted.file)
  cat("\nSorting bam file", bam.file, "to", sorted.file, "\n")
  call.cmd(sort)
  
  #Calculate depth
  ###################
  sorted.file = paste(sorted.file, ".bam", sep="")
  coverage.file = file.path(dir, paste(bam.file.name, ".coverage.txt", sep=""))
  depth = paste(samtools.path, "depth", sorted.file, ">", coverage.file)
  cat("\nCalculating depth of coverage from bam file", sorted.file, "to", coverage.file, "\n")
  call.cmd(depth)
  return(coverage.file)
}
