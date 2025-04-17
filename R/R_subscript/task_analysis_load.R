dir <- dirname(sys.frame(1)$ofile)

f <- list.files(paste0(dir,"/task_analysis_subscripts"))

sapply(paste0(dir,"/task_analysis_subscripts/",f),source)
