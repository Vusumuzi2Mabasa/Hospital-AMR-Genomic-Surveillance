# Installation is isolated to the project library and does not change other projects.
dir.create('renv/library',recursive=TRUE,showWarnings=FALSE)
.libPaths(c(normalizePath('renv/library'),.libPaths()))
if(!requireNamespace('renv',quietly=TRUE))install.packages('renv',repos='https://cloud.r-project.org')
renv::install(c('ggplot2@4.0.0','ape@5.8-1','jsonlite@2.0.0','digest@0.6.37'))
