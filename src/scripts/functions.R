
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

tel.length <- function(coverage.file, fastqs,
                       rl, pl,
                       base.cov, num.haploid.chr, 
                       genome.length, min.seed,
                       out.file){
  
  out.list = list()
  out.list$coverage.file = coverage.file
  out.list$reads = paste(fastqs, collapse = ", ")
  out.list$read.length = rl
  out.list$pattern.length = pl
  out.list$base.cov = base.cov
  out.list$num.haploid.chr = num.haploid.chr
  out.list$tel.length = 0  
  out.list$genome.length = genome.length
  out.list$min.seed = min.seed  
 
  
  
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
  cat("genome.length=",  genome.length, "\n")
  
  write.table(out.table, file = out.file, sep="\t", col.names = F, quote = F)
  
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
 # the paired end implementation should be remove, so will not bother with it anymore
  length.file = file.path(output.dir, "length")
  U.cat = gsub(",", replacement=" ", U)
  command = as.character(paste("cat", U.cat, "|",  
                                 "tee >(wc -l >", length.file,") |", bowtie.align.path, 
                                 options, index, "-U - ", samout))
#  command = as.character(paste(bowtie.align.path, options,index, reads, samout));   
  cat(command)
  call.out = call.cmd(command)
  success = (file.exists(S) && file.info(S)$size > 0)
  if (success){
    if (!ignore.err){
      success = (call.out == 0)
    }
    if(success){
      success = (file.exists(length.file) && file.info(length.file)$size > 0)
      if(!success){
        cat("could not find file ", length.file , " for base coverage estimation \n")
      }
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
  
  length.file = file.path(output.dir, "length")
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
                                 "tee >(wc -l >", length.file,") |", bowtie.align.path, 
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
      success = (file.exists(length.file) && file.info(length.file)$size > 0)
      if(!success){
        cat("could not find file ", length.file , " for base coverage estimation \n")
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
  cat("\nConverting sam file", sam.file,"to bam file", bam.file, "\n")
  cat("\ncommand: ", samtobam, "\n")
  call.cmd(samtobam)
  if(!file.exists(bam.file)){
  	return("sam to bam conversion failed")
  }
  
  #Sort bam
  ###################
  sorted.file = file.path(dir, paste(bam.file.name, "_sorted",sep=""))
  sort = paste(samtools.path, "sort", bam.file, "-o", sorted.file)
  cat("\nSorting bam file", bam.file, "to", sorted.file, "\n")
  cat("\ncommand: ", sort, "\n")
  call.cmd(sort)
  
  #Calculate depth
  ###################
  if(!file.exists(sorted.file)){
  	return("bam sorting failed")
  }

  sorted.file = paste(sorted.file, sep="")
  coverage.file = file.path(dir, paste(bam.file.name, ".coverage.txt", sep=""))
  depth = paste(samtools.path, "depth -d 1000000000", sorted.file, ">", coverage.file)
  cat("\nCalculating depth of coverage from bam file", sorted.file, "to", coverage.file, "\n")
  cat("\ncommand: ", depth, "\n")
  call.cmd(depth)
  return(coverage.file)
}



############ telomeric repeat variant functions ##############


rev.comp <- function(seq){
  chars = strsplit(seq, split = "")[[1]]
  comp.chars = chars 
  comp.chars[which(chars == "A")] = "T"
  comp.chars[which(chars == "T")] = "A"
  comp.chars[which(chars == "G")] = "C"
  comp.chars[which(chars == "C")] = "G"
  
  rev.comp.seq = paste0(comp.chars[length(comp.chars):1], collapse = "")
  return(rev.comp.seq)
}

is.telomeric <- function(token, normal.pattern){
  permutted = paste0(normal.pattern, substr(normal.pattern, 1, nchar(normal.pattern)-1), collapse = "")
  return(grepl(token, permutted))
}

add.one.variant <- function(variant,n, variants, tokenqual = NULL, meanQual = NULL){
  # cat(variant, "\n")
  if(variant %in% names(variants))
    variants[variant] = variants[variant] + n
  else
    variants[variant] = n
  return(variants)
}


process.cigar <-function(cigar){
  
  d.ind  = c()
  if(grepl("D", cigar)){
    d.ind = gregexpr("D", cigar)[[1]]
  }
  
  i.ind = c()
  if(grepl("I", cigar)){
    i.ind = gregexpr("I", cigar)[[1]]
  }
  
  m.ind = c()
  if(grepl("M", cigar)){
    m.ind = gregexpr("M", cigar)[[1]]
  }
  breakpoints = c(m.ind, d.ind, i.ind)
  breakpoints = sort(breakpoints, decreasing = F)
  
  cigar.string = c()
  for(b in 1:length(breakpoints)){
    if(b == 1){
      numeral = as.numeric(substr(cigar, 1, breakpoints[b]-1))
      if(breakpoints[b] %in% m.ind)
        c = 'M'
      else if (breakpoints[b] %in% d.ind)
        c = 'D'
      else 
        c = 'I'
      cigar.string = rep(c, numeral)
    } else {
      numeral = as.numeric(substr(cigar, breakpoints[b-1]+1, breakpoints[b]-1))
      if(breakpoints[b] %in% m.ind)
        c = 'M'
      else if (breakpoints[b] %in% d.ind)
        c = 'D'
      else 
        c = 'I'
      cigar.string = c(cigar.string, rep(c,numeral))
    }
  }
  
  return(cigar.string)
}

process.line <- function(line, variants, silent = T, normal.pattern){
  sam = strsplit(line, split="\t", fixed = T)[[1]]
  seq = as.character(sam[10])
  qual = utf8ToInt(as.character(sam[11]))
  cigar = sam[6]
  pos = as.numeric(sam[4])
  
  if(!silent){
    cat(seq, "\n")
    cat(qual, "\n")
    cat(cigar, "\n")
  }
  
  cigar.string = process.cigar(cigar)
  
  pn = nchar(normal.pattern)
  tel.vec = strsplit(normal.pattern, split = "")[[1]]
  qual.threshold = mean(qual)
  
  
  if(pos > 1){
    offset = (pos-1) %% pn
    if(offset != 0){
      prefix = substr(normal.pattern, 1, offset)
      
      cigar.string = c(rep("M", offset), cigar.string)
      qual = c(rep(max(qual), offset), qual)
      seq = paste0(prefix, seq, collapse = "")
    }
  }
  seq.cursor = 1
  cigar.cursor = 1
  while(seq.cursor + pn -1 <= nchar(seq)){
    base.count = 1
    add =  T
    var = ""
    
    while(base.count <= pn){
      if(seq.cursor > nchar(seq))
        break()
      base = substr(seq, seq.cursor, seq.cursor)
      tel.base = tel.vec[base.count]
      base.cigar = cigar.string[cigar.cursor]
      if(base.cigar == "M"){
        if(base != tel.base){
          base.qual = qual[seq.cursor]
          if(base.qual < qual.threshold)
            add = F
        }
        base.count = base.count + 1
        seq.cursor = seq.cursor + 1
        cigar.cursor = cigar.cursor + 1
        var = paste0(var, base)
      } else if (base.cigar == "I"){
        base.count = base.count
        seq.cursor = seq.cursor + 1
        cigar.cursor = cigar.cursor + 1
        base.qual = qual[seq.cursor]
        if(base.qual < qual.threshold)
          add = F 
        var = paste0(var, base)
      } else{ # base.cigar = "D" 
        base.count = base.count + 1
        seq.cursor = seq.cursor
        cigar.cursor = cigar.cursor + 1
        var = var 
      }
    }
    if(add){
      if(!silent)
        cat("Added:\t", var, "\n")
      variants = add.one.variant(var, 1, variants)
    } else {
      if(!silent)
        cat("Disqualified:\t", var, "\n")
    }
  }
  
  return(variants)
}



variant.counter <- function(file, pattern){
  
  #cat(file, "\n")
  con = file(file, open = 'r')
  
  normal.pattern = rev.comp(pattern)
  variants = c()
  
  
  while(TRUE){
    line = readLines(con, n = 1)
    if(length(line) == 0)
      break()
    if(substr(line,1,1) != "@")
      break()
  }

  if(length(line) != 0) {

  count = 1
  while(TRUE){
    # cat(count, "\n")
    variants = process.line(line, variants, T, normal.pattern)
    line = readLines(con, n = 1)
    if(length(line) == 0 )
      break()
    count = count + 1
  }
  
  }
  close(con)
  
  return(variants)
}

count.variants <- function(file, normal.pattern, out.file){
  variants = variant.counter(file, normal.pattern)
  
  variants = sort(variants, decreasing = T)
  variants.perc = (variants/sum(variants))*100
  
  trv.rel = 100*variants.perc[-1]/sum(variants.perc[-1])
  
  
  
  trv.colnames = c("pattern", "abs num", "% of all patterns", "% of variants")
  trv.mat = matrix(nrow = length(variants), ncol = length(trv.colnames))
  trv.mat[,1] = unlist(lapply(names(variants), rev.comp))
  trv.mat[,2] = variants
  trv.mat[,3] = variants.perc
  trv.mat[,4] = c("0", trv.rel)
  colnames(trv.mat) = trv.colnames
  write.table(trv.mat, out.file, quote = F, sep = "\t", row.names = F)
  if(file.exists(out.file))
    return(TRUE)
  else
    return(FALSE)
}


