set.property <- function(matrix, prop.name, default.value){
  result = tryCatch({
    matrix[prop.name,1]
  }, warning = function(w) {
    cat("property \"", prop.name, "\" not found. Default value assigned: ", default.value, "\n")
    default.value
  }, error = function(e) {
    cat("property \"", prop.name, "\" not found. Default value assigned: ", default.value, "\n")
    default.value
  })
  return(result)
}

set.property.double <- function(matrix, prop.name, default.value){
  result = set.property(matrix, prop.name, default.value)
  result = tryCatch({
    as.double(result)
  }, warning = function(w) {
    cat("Warning: Illegal value for", prop.name, ": ", result,
        "Should be double. Default value assigned: ", default.value, "\n")
    default.value
  }, error = function(e) {
    cat("Error: Illegal value for", prop.name, ": ", result,
        "Should be double. Default value assigned: ", default.value, "\n")
    default.value
  })
}

set.property.integer <- function(matrix, prop.name, default.value){
  result = set.property(matrix, prop.name, default.value)
  result = tryCatch({
    as.integer(result)
  }, warning = function(w) {
    cat("Warning: Illegal value for", prop.name, ": ", result,
        "Should be integer. Default value assigned: ", default.value, "\n")
    default.value
  }, error = function(e) {
    cat("Error: Illegal value for", prop.name, ": ", result,
        "Should be integer. Default value assigned: ", default.value, "\n")
    default.value
  })
}
set.property.logical <- function(matrix, prop.name, default.value){
  result = set.property(matrix, prop.name, default.value)
  result = tryCatch({
    as.logical(result)
  }, warning = function(w) {
    cat("Warning: Illegal value for", prop.name, ": ", result,
        "Should be logical. Default value assigned: ", default.value, "\n")
    default.value
  }, error = function(e) {
    cat("Error: Illegal value for", prop.name, ": ", result,
        "Should be logical. Default value assigned: ", default.value, "\n")
    default.value
  })
}

set.property.sequence <- function(matrix, prop.name, default.value){
  result = set.property(matrix, prop.name, default.value)
  result = tryCatch({
    as.character(result)
    result = toupper(result)
    dna.chars = c ("A", "G", "C","T")
    chars = strsplit(result,"")[[1]]
    is.dna = all(chars %in% dna.chars)
    if (!is.dna)
      cat("Illegal value for ", prop.name, ": ", result,
          "Should contain DNA sequence characters ", dna.chars, ". Defualt value assigned: ",
          default.value, "\n")
    result
  }, warning = function(w) {
    cat("Warning: Illegal value for", prop.name, ": ", result,
        "Should be numeric. Default value assigned: ", default.value, "\n")
    default.value
  }, error = function(e) {
    cat("Error: Illegal value for", prop.name, ": ", result,
        "Should be numeric. Default value assigned: ", default.value, "\n")
    default.value
  })
}

set.property.executable <- function(matrix, prop.name, default.value){
  result = tryCatch({
    matrix[prop.name,1]
  }, warning = function(w) {
    cat("property \"", prop.name, "\" not found. Default value assigned: ", default.value, "\n")
    default.value
  }, error= function(e) {
    cat("property \"", prop.name, "\" not found. Default value assigned: ", default.value, "\n")
    default.value
  })
  if (length(grep("./", result, fixed = T)) > 0
      || length(grep(".\\", result, fixed = T)) > 0){
    wd = getwd()
    setwd("..")
    parent = getwd()
    setwd(wd)
  }

  if (length(grep(pattern="../", result, fixed = T)) > 0){
    result = sub(pattern="..", replacement=parent, result, fixed = T)
  } else if (length(grep(pattern="./", result, fixed = T)) > 0){
    result = sub(pattern=".", replacement=wd, result, fixed = T)
  } else   if (length(grep(pattern="..\\", result, fixed = T)) > 0){
    result = sub(pattern="..", replacement=parent, result, fixed = T)
  } else if (length(grep(pattern=".\\", result, fixed = T)) > 0){
    result = sub(pattern=".", replacement=wd, result, fixed = T)
  }

  return(result)
}

defaults = list(
  rl = 100,
  pattern = "TTAGGG",
  base.cov = 1,
  num.haploid.chr = 23,
  output.dir = "output",
  compute.base.cov = F,
  estimate.base.cov = F,
  genome.length = 3101804739,
  mode.local = F,
  single = T,
  files.with.prefix = F,
  file.compression = F,
  min.seed = 12,
  num.proc = 3,
  scripts.dir = "./",
  bowtie.dir = "../bowtie2-2.1.0/",
  bowtie.build.path = "..bowtie2-2.1.0/bowtie2-build.exe",
  bowtie.align.path = "../bowtie2-2.1.0/bowtie2-align.exe",
  samtools.path = "../samtools-0.1.19/samtools.exe",
  quals = "--phred33",
  ignore.err = F
)
rl = set.property.integer(config.table, "read.length", defaults$rl)
pattern = set.property.sequence(config.table, "pattern", defaults$pattern)

num.haploid.chr = set.property.integer(config.table, "num.haploid.chr", defaults$num.haploid.chr)
compute.base.cov = set.property.logical(config.table, "compute.base.cov",  defaults$compute.base.cov)

mode.local = set.property.logical(config.table, "mode.local", defaults$mode.local)
single = set.property.logical(config.table, "single", defaults$single)
files.with.prefix = set.property.logical(config.table, "files.with.prefix", defaults$files.with.prefix)
file.compression = set.property(config.table, "file.compression", defaults$file.compression)
compressed = F
OS.name = tolower(Sys.info()["sysname"])
if (OS.name == "windows") {
  if(file.compression != "F"){
    config.set = F
    cat("file.compression should be set to F for windows systems\n")
  }
} else {
  if(file.compression != "F"){
    if(file.compression == "gz"){
      compressed = T
    } else if (file.compression == "bz2"){
      compressed = T
    } else {
      compressed = F
      cat("file.compressed should be either gz or bz2\n")
      config.set = F
    }
  }
  if(compressed){
    if(!single){
      cat("compressed fastq files should be given with \"single\" parameter set to T")
      config.set = F
    }
  }
}
estimate.base.cov = set.property.logical(config.table, "estimate.base.cov", defaults$estimate.base.cov)
if(estimate.base.cov){
#  if(!compressed){
#    cat("Warning: base coverage is estimated only from gzip or bzip2 compressed files. Default base coverage", defaults$base.cov, #"assigned.\n")
#    estimate.base.cov = F
#    base.cov = defaults$base.cov
#  }
  if(compute.base.cov){
    cat("Error: compute.base.cov and estimate.base.cov cannot be true at the same time\n")
    config.set = F
  }
}
genome.length = set.property.double(config.table, "genome.length", defaults$genome.length)
if(genome.length < rl){
  cat("\nError: genome length (lgenome) cannot be less than the read length\n")
  config.set = F
}

min.seed = set.property.integer(config.table, "min.seed", defaults$min.seed)
if(min.seed > rl){
  cat("\nError: min.seed cannot exceed the read length\n")
  config.set = F
}
if(min.seed < defaults$min.seed){
  cat("\nError: min.seed cannot be less than 12\n")
  config.set = F
}
num.proc = set.property.integer(config.table, "num.proc", defaults$num.proc)

quals = set.property(config.table, "quals", defaults$quals)
if (quals != "--phred33"){
  if (quals != "--phred64"){
    if(quals != "--solexa-quals")
      quals = defaults$quals
  }
}

ignore.err = set.property.logical(config.table, 'ignore.err', defaults$ignore.err)

if(single){
  if (!files.with.prefix){
    fastq = tryCatch({
      config.table["fastq",1]
    }, warning = function(w) {
      cat("\"fastq\" parameter not specified for single mode alignment\n")
      config.set = F
    }, error = function(e) {
      cat("\"fastq\" parameter not specified for single mode alignment\n")
      config.set = F
    })
    names(fastq) <- NULL

    fastqs = unlist(strsplit(fastq, split='[ ,]' , fixed=F))
    fastqs = fastqs[!duplicated(fastqs)]
    empty.fastqs = c()
    for (i in 1:length(fastqs)){
      if (identical(fastqs[i], "")) {
        empty.fastqs = c(empty.fastqs,i)
      } else if (!file.exists(fastqs[i])){
        cat("fastq file ", fastqs[i],
            " does not exist. Please, specify a valid file for single mode alignment.\n")
        config.set = F
      }
    }
    if(length(empty.fastqs) > 0)
      fastqs=fastqs[-empty.fastqs]
  } else {
    fastq.prefix = tryCatch({
      config.table["fastq.prefix",1]
    }, warning = function(w) {
      cat("\"fastq.prefix\" parameter not specified for single mode alignment\n")
      config.set = F
    }, error = function(e) {
      cat("\"fastq.prefix\" parameter not specified for single mode alignment\n")
      config.set = F
    })
    fastq.dir = tryCatch({
      config.table["fastq.dir",1]
    }, warning = function(w) {
      cat("\"fastq.dir\" parameter not specified for single mode alignment\n")
      config.set = F
    }, error = function(e) {
      cat("\"fastq.dir\" parameter not specified for single mode alignment\n")
      config.set = F
    })
    if(compressed){
      suffix = file.compression
      zip.pattern = paste(fastq.prefix, ".*","\\.", suffix,"$", sep="")
      fastq.files = list.files(path=fastq.dir, zip.pattern)
    } else {
      suffix = ""
      fastq.files = list.files(path=fastq.dir, paste(fastq.prefix, ".*\\.f.*q$", sep=""))
    }


    fastqs = file.path(fastq.dir, fastq.files)
    if (length(list.files) == 0)
      cat("no files with pattern ", zip.pattern, " were present in the directory ", fastq.dir)
  }

  #   check for compression type
  if(compressed){
    suffix = file.compression
    for (fastq in fastqs){
      sout = system(paste("file",fastq), intern = T)
      if(suffix == "gz")
        grout = grep("gzip compressed", sout)
      else if(suffix == "bz2")
        grout = grep("bzip2 compressed", sout)
    }
    if(length(grout) == 0){
#     zip.pattern = paste(".*\\.", suffix, "$",sep="")
#     if(length(grep(zip.pattern, fastqs)) < length(fastqs))
      cat("fastq file ", fastq, " was not ", suffix, " compressed", "\n")
      config.set = F
    }
  }
} else {
  fastq1 = tryCatch({
    config.table["fastq1",1]
  }, warning = function(w) {
    cat("\"fastq1\" parameter not specified for paired mode alignment\n")
    config.set = F
  }, error = function(e) {
    cat("\"fastq1\" parameter not specified for paired mode alignment\n")
    config.set = F
  })

  if (!file.exists(fastq1)){
    cat("fastq1 file ", fastq1,
        " does not exist. Please, specify a valid file for paired mode alignment.\n")
    config.set = F
  }

  fastq2 = tryCatch({
    config.table["fastq2",1]
  }, warning = function(w) {
    cat("\"fastq2\" parameter not specified for paired mode alignment\n")
    config.set = F
  }, error = function(e) {
    cat("\"fastq2\" parameter not specified for paired mode alignment\n")
    config.set = F
  })

  if (!file.exists(fastq1)){
    cat("fastq2 file ", fastq1,
        " does not exist. Please, specify a valid file for paired mode alignment.\n")
    config.set = F
  }
}


base.index.exists <- function(prefix, suffix){
  filename = paste(prefix, suffix, sep="")
  if (!file.exists(filename)){
    cat("base index ", filename, " could not be found\n")
    return(F)
  } else {
    return(T)
  }
}

if (compute.base.cov){
  base.index.pathtoprefix = set.property(config.table, "base.index.pathtoprefix", "n/a")
  if(base.index.pathtoprefix == "n/a"){
    cat("compute.base.cov was equal to T, but base.index.pathtoprefix was not specified.")
    base.index.set = F
  } else{
    base.index.set = F
    if (base.index.exists(base.index.pathtoprefix, ".1.bt2"))
      if (base.index.exists(base.index.pathtoprefix, ".2.bt2"))
        if (base.index.exists(base.index.pathtoprefix, ".3.bt2"))
          if (base.index.exists(base.index.pathtoprefix, ".4.bt2"))
            if (base.index.exists(base.index.pathtoprefix, ".rev.1.bt2"))
              if (base.index.exists(base.index.pathtoprefix, ".rev.2.bt2"))
                base.index.set = T
  }
  if (!base.index.set)
    config.set = F
  picard.samtofastq.jar = set.property(config.table, "picard.samtofastq.jar", "n/a")
  if(picard.samtofastq.jar == "n/a"){
    cat("compute.base.cov was equal to T, but base.index.pathtoprefix was not specified.")
    config.set = F
  } else{
    if(!file.exists(picard.samtofastq.jar)){
      cat("Specified picard.samtofastq.jar file ", picard.samtofastq.jar, "does not exist\n")
      config.set = F
    }
  }
} else if(estimate.base.cov){
	base.cov = defaults$base.cov
} else {
  base.cov = set.property.double(config.table, "base.cov", defaults$base.cov)
}

output.dir = set.property(config.table, "output.dir", defaults$output.dir)
if(is.na(file.info(output.dir)$isdir || !file.info(output.dir)$isdir)){
  dc = dir.create(output.dir, showWarnings=F)
  if (!dc){
    cat("Warning: could not create output directory", output.dir,
        "results will be kept in default directory", defaults$output.dir, "\n")
    output.dir = defaults$output.dir
    dir.create(output.dir, showWarnings=F)
  }
}

scripts.dir = set.property(config.table, "scripts.dir", defaults$scripts.dir)
scripts.dir.set = F
#if ( file.info(scripts.dir)$isdir){
if (file.exists(file.path(scripts.dir, "pipeline.R")))
  if (file.exists(file.path(scripts.dir, "functions.R"))){
    scripts.dir.set = T
  }
#}
if (!scripts.dir.set){
  cat("scripts.dir does not exist or does not contain the required scripts\n")
  config.set = F
}

bowtie.build.path = set.property.executable(config.table, "bowtie.build.path", defaults$bowtie.build.path)
if (!file.exists(bowtie.build.path)){
  cat("bowtie.build executable could not be found at ", bowtie.build.path, "\n")
  config.set = F
}

bowtie.align.path = set.property.executable(config.table, "bowtie.align.path", defaults$bowtie.align.path)
if (!file.exists(bowtie.align.path)){
  cat("bowtie.align executable not be found at ", bowtie.align.path, "\n")
  config.set = F
}

samtools.path = set.property.executable(config.table, "samtools.path", defaults$samtools.path)
if (!file.exists(samtools.path)){
  cat("samtools.path", samtools.path, "does not exist\n")
  config.set = F
}


if (config.set) {

  cat("\nConfiguration successful\n")
  cat("Computel will execute with the following parameters:\n")
  cat("scripts.dir =", scripts.dir, "\n")
  cat("bowtie.build.path =", bowtie.build.path, "\n")
  cat("bowtie.align.path =", bowtie.align.path, "\n")
  cat("samtools.path =", samtools.path, "\n")
  if(compute.base.cov)
    cat("picard.samtofastq.jar =", picard.samtofastq.jar, "\n")
  cat("single =", single, "\n")
  if (single) {
    cat("fastqs =", fastqs, "\n")
  } else {
    cat("fastq1 =", fastq1, "\n")
    cat("fastq2 =", fastq2, "\n")
  }
  if(compressed){
    cat("file.compression =", file.compression, "\n")
  }
  cat("read.length = ", rl,"\n")
  cat("pattern = ", pattern,"\n")
  cat("num.haploid.chr = ", num.haploid.chr,"\n")
  cat ("min.seed =", min.seed, "\n")
  cat ("mode.local =", mode.local, "\n")
  cat("compute.base.cov =", compute.base.cov, "\n")
  if (compute.base.cov){
    cat("base.index.pathtoprefix =", base.index.pathtoprefix, "\n")
  } else if(estimate.base.cov){
    cat("base.cov =  will be estimated from the zipped files\n")
  } else {
    cat("base.cov = ", base.cov,"\n")
  }

  cat("output.dir= ", output.dir,"\n")
  cat ("num.proc=", num.proc, "\n")
  cat("quals=", quals, "\n")
  cat("ingore.err=", ignore.err, "\n")
  pipeline.R = file.path(scripts.dir, "pipeline.R")

}
