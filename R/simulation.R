simulate_population <- function(seed,cfg) {
  set.seed(seed);n<-cfg$isolates;L<-cfg$core_sites
  ancestor<-sample(0:3,L,replace=TRUE)
  founders<-replicate(12,{s<-ancestor;j<-sample.int(L,60);s[j]<-(s[j]+sample(1:3,60,replace=TRUE))%%4;s})
  lineage<-sample(1:12,n,replace=TRUE);outbreak<-integer(n);outbreak[1:36]<-rep(1:3,each=12);lineage[1:36]<-outbreak[1:36]
  # Six close but epidemiologically unlinked comparators challenge genomic-only rules.
  lineage[37:42]<-rep(1:3,each=2)
  seqs<-t(vapply(lineage,function(g){s<-founders[,g];j<-sample.int(L,rpois(1,2));if(length(j))s[j]<-(s[j]+sample(1:3,length(j),replace=TRUE))%%4;as.integer(s)},integer(L)))
  hospital<-sample(c('Hospital A','Hospital B','Hospital C'),n,replace=TRUE)
  ward<-sample(c('Medical','Surgical','ICU'),n,replace=TRUE);day<-sample.int(365,n,replace=TRUE)-1
  for(k in 1:3){j<-which(outbreak==k);hospital[j]<-c('Hospital A','Hospital B','Hospital C')[k];ward[j]<-c('ICU','Medical','Surgical')[k];day[j]<-c(50,160,280)[k]+sample(0:20,length(j),replace=TRUE)}
  day[37:42]<-c(330,335,20,30,80,85)
  age<-round(pmax(18,pmin(90,rnorm(n,55,17))));severity<-rbinom(n,1,.3)
  gene<-rbinom(n,1,plogis(-1+.65*(hospital=='Hospital B')+.9*(outbreak>0)))
  gene_name<-ifelse(gene==0,'None',sample(c('blaKPC-2','blaNDM-1'),n,replace=TRUE))
  phenotype<-rbinom(n,1,ifelse(gene==1,.9,.08))
  death<-rbinom(n,1,plogis(-2.5+.5*gene+.8*severity+.015*(age-55)))
  meta<-data.frame(isolate=sprintf('SIM%03d',seq_len(n)),hospital,ward,day,date=as.Date('2025-01-01')+day,month=format(as.Date('2025-01-01')+day,'%Y-%m'),age,severity,gene,gene_name,phenotype,death,outbreak,lineage)
  # All 6,000 sites are known in truth; these are toy core haplotypes, not K. pneumoniae genomes.
  distance<-matrix(0L,n,n)
  for(i in seq_len(n-1))for(j in (i+1):n)distance[i,j]<-distance[j,i]<-sum(seqs[i,]!=seqs[j,])
  list(metadata=meta,sequences=seqs,distance=distance)
}
observe_population<-function(pop,scenario,seed,cfg) {
  set.seed(seed);m<-pop$metadata;n<-nrow(m);pairs<-which(upper.tri(pop$distance),arr.ind=TRUE)
  quality<-if(scenario=='poor_quality')runif(n,.87,.99) else rep(.995,n)
  q<-quality[pairs[,1]]*quality[pairs[,2]]
  shared<-rbinom(nrow(pairs),cfg$core_sites,q)
  true_d<-pop$distance[pairs]
  observed<-rbinom(length(true_d),true_d,q)+rpois(length(true_d),2*cfg$core_sites*q*if(scenario=='poor_quality')2e-4 else 1e-6)
  sampled<-if(scenario=='incomplete_sampling')rbinom(n,1,.6)==1 else rep(TRUE,n)
  if(scenario=='missing_metadata')m$ward[rbinom(n,1,.25)==1]<-NA_character_
  truth<-m$outbreak[pairs[,1]]>0&m$outbreak[pairs[,1]]==m$outbreak[pairs[,2]]
  epi<-m$hospital[pairs[,1]]==m$hospital[pairs[,2]]&!is.na(m$ward[pairs[,1]])&!is.na(m$ward[pairs[,2]])&m$ward[pairs[,1]]==m$ward[pairs[,2]]&abs(m$day[pairs[,1]]-m$day[pairs[,2]])<=cfg$maximum_days
  epi[is.na(epi)]<-FALSE
  data.frame(i=pairs[,1],j=pairs[,2],truth,observed,shared,available=sampled[pairs[,1]]&sampled[pairs[,2]],qc=shared/cfg$core_sites>=cfg$minimum_callable_fraction,epi)
}
pair_metrics<-function(d,k,method) {
  eligible<-d$available&d$qc
  called<-eligible&d$observed<=k
  if(method=='Genomic + epidemiology')called<-called&d$epi
  tp<-sum(called&d$truth);fp<-sum(called&!d$truth);fn<-sum(!called&d$truth)
  available_true<-sum(d$truth&d$available)
  data.frame(threshold=k,method,TP=tp,FP=fp,FN=fn,precision=if(tp+fp)tp/(tp+fp) else NA_real_,sensitivity_all=tp/sum(d$truth),sensitivity_sampled=if(available_true)tp/available_true else NA_real_,false_positive_rate=fp/sum(!d$truth),available_pairs=sum(d$available),qc_eligible_pairs=sum(eligible),true_pairs=sum(d$truth),sampled_true_pairs=available_true)
}
run_simulations<-function(cfg){
  rows<-list();example<-NULL
  for(r in seq_len(cfg$replicates)){
    pop<-simulate_population(cfg$seed+r,cfg)
    if(r==1)example<-pop
    for(s in cfg$scenarios){d<-observe_population(pop,s,cfg$seed+10000+r,cfg);for(k in cfg$thresholds)for(method in c('Genomic only','Genomic + epidemiology'))rows[[length(rows)+1]]<-cbind(replicate=r,scenario=s,pair_metrics(d,k,method))}
  }
  list(performance=do.call(rbind,rows),example=example)
}
