args <- commandArgs(TRUE)

config.file = args[1]

if (length(args) < 1){  
  cat("\nUsage: \n\n Rscript computel.cmd.R [path.to.config.file]\n\n")
  stop("No config file provided")
}
  
if (file.exists(config.file))
  config.set = T else stop("No valid config file provided")



#Read config
config.table = as.matrix(read.table(file=config.file, header=F, sep="\t",
                                    blank.lines.skip=T, as.is=TRUE, row.names=1))


prop.name = 'scripts.dir'
default.value = './'
scripts.dir = tryCatch({
    config.table[prop.name,1]
  }, warning = function(w) {
    default.value
  }, error = function(e) {
    default.value
  })
validate.R = file.path(scripts.dir, "validate.options.R")
if (!file.exists(validate.R)){
  validate.R = "validate.options.R"   
}

if (!file.exists(validate.R)){
  stop(paste("the script", validate.R, "was not found.\n Provide a valid configuration for scripts.dir, containing the script."))
} else {
  source(validate.R)
}
  
if (!config.set){
  stop("configuration not set successfully. Scripts will not execute.\n")
} else {
  source(pipeline.R)
}

