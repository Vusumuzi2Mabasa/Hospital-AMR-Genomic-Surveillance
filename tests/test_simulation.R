source('scripts/load.R');cfg<-jsonlite::read_json('config/study.json',simplifyVector=TRUE)
a<-simulate_population(123,cfg);b<-simulate_population(123,cfg)
stopifnot(identical(a,b),nrow(a$metadata)==150,ncol(a$sequences)==6000,!anyDuplicated(a$metadata$isolate),all(diag(a$distance)==0),identical(a$distance,t(a$distance)))
manual<-sum(a$sequences[1,]!=a$sequences[2,]);stopifnot(manual==a$distance[1,2])
for(s in cfg$scenarios){
  d<-observe_population(a,s,321,cfg);stopifnot(nrow(d)==choose(150,2),sum(d$truth)==3*choose(12,2),all(d$shared>=0&d$shared<=6000),!anyNA(d$epi))
  g<-pair_metrics(d,5,'Genomic only');e<-pair_metrics(d,5,'Genomic + epidemiology')
  stopifnot(e$TP<=g$TP,e$FP<=g$FP,e$sensitivity_all<=g$sensitivity_all,g$TP+g$FN==sum(d$truth))
  d$truth<-rev(d$truth);gg<-pair_metrics(d,5,'Genomic only');stopifnot(gg$TP+gg$FP==g$TP+g$FP)
  d$qc[]<-FALSE;z<-pair_metrics(d,5,'Genomic only');stopifnot(z$TP==0,z$FP==0,z$sensitivity_all==0,is.na(z$precision))
}
d<-observe_population(a,'baseline',321,cfg);lo<-pair_metrics(d,2,'Genomic only');hi<-pair_metrics(d,10,'Genomic only');stopifnot(lo$TP<=hi$TP,lo$FP<=hi$FP)
# Close genomes must not alone be interpreted as a planted outbreak.
stopifnot(any(d$observed<=5&!d$truth))
cat('Determinism, dimensions, exact distances, truth separation, QC abstention and all four scenarios passed.\n')
