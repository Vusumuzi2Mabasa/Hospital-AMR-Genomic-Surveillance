# Current epidemiological report; original builder preserved in the archive.
exe<-file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
stopifnot(system2(exe,c('--vanilla','epidemiology/analysis.R','--report'))==0)
