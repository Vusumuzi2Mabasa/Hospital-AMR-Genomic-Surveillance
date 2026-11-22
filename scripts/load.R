# Prefer existing compatible local libraries; no package installation occurs here.
library_candidates<-c(Sys.getenv('HOSPITAL_AMR_R_LIBRARY'),'renv/library','D:/AKI-Proteogenomics-Demo/AKI-Proteomics-R/renv/library/windows/R-4.5/x86_64-w64-mingw32')
.libPaths(c(library_candidates[dir.exists(library_candidates)],.libPaths()))
for(p in c('jsonlite','ggplot2','ape','digest'))if(!requireNamespace(p,quietly=TRUE))stop('Missing R package: ',p,'. See environment.yml / scripts/restore.R.')
for(f in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(f)
