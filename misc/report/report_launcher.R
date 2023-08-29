#!/usr/bin/env Rscript

rm(list=ls())

#Rscript report_launcher.R  $proj_dir $proj_name_pref $data_dir $metadata_dir $setup $library &> ovcar3.report.stderrout


library(knitr)
args <- commandArgs(TRUE)

if (length(args) < 6) stop("Not all args are set; required: projdir proj.name.prefix data.dir mageck.dir metadata.dir setup library")

proj.dir <- args[1]
proj.name.pref <- args[2]
data.dir <- args[3]
mageck.dir <- args[4]
metadata.dir <- args[5]
setup <- args[6] ### cpossible values: all 0rm 0rm_noAlt
library <- args[7]


dir.create(proj.dir, recursive = TRUE)

wrk.dir=file.path(proj.dir,paste(proj.name.pref,"report",sep="."))
dir.create(wrk.dir, recursive = TRUE)

print(proj.dir)
print(proj.name.pref)
print(wrk.dir)
print(library)
print("Data used for report:")
print(data.dir)
print(mageck.dir)
print(metadata.dir)
print(paste("Setup",setup,sep=" "))

print(getwd())

rmarkdown::render('6257_screen_report_template_v0.2.Rmd', output_file = file.path(wrk.dir,paste("MAGeCK_report",proj.name.pref,setup,'html', sep=".")))


