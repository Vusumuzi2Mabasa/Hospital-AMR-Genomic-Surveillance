analyse_example<-function(pop,cfg){
  m<-pop$metadata
  monthly<-aggregate(rep(1,nrow(m)),m[,c('hospital','month')],sum);names(monthly)[3]<-'isolates';csv(monthly,'results/monthly_counts.csv')
  prevalence<-do.call(rbind,lapply(split(m,m$hospital),function(d){ci<-wilson(sum(d$phenotype),nrow(d));data.frame(hospital=d$hospital[1],n=nrow(d),resistant=sum(d$phenotype),proportion=mean(d$phenotype),lower=ci[1],upper=ci[2])}));csv(prevalence,'results/resistance_proportions.csv')
  fit<-glm(death~gene+age+severity+hospital,data=m,family=binomial())
  b<-coef(summary(fit));or<-data.frame(term=rownames(b),OR=exp(b[,1]),lower=exp(b[,1]-qnorm(.975)*b[,2]),upper=exp(b[,1]+qnorm(.975)*b[,2]),p=b[,4]);csv(or,'results/mortality_associations.csv')
  ft<-fisher.test(table(m$hospital,m$phenotype))
  ct<-table(factor(m$gene,levels=0:1),factor(m$phenotype,levels=0:1))
  csv(as.data.frame(ct),'results/genotype_phenotype_counts.csv')
  jsonlite::write_json(list(hospital_resistance_fisher_p=ft$p.value,mortality_model_converged=fit$converged,deaths=sum(m$death),n=nrow(m),interpretation='Exploratory associations in one synthetic cohort; no causal or clinical claim'),'results/example_statistics.json',pretty=TRUE,auto_unbox=TRUE)
  tree<-ape::nj(as.dist(pop$distance));tree$tip.label<-m$isolate;ape::write.tree(tree,'results/full_synthetic_tree.nwk')
  ids<-m$isolate[c(1:6,13:18,25:30,37:42,43:54)]
  small<-ape::ladderize(ape::keep.tip(tree,ids));ape::write.tree(small,'results/display_synthetic_tree.nwk')
  observed<-observe_population(pop,'baseline',cfg$seed+10001,cfg)
  edges<-observed[observed$qc&observed$available&observed$epi&observed$observed<=cfg$primary_threshold,c('i','j')]
  membership<-seq_len(nrow(m));if(nrow(edges))for(k in seq_len(nrow(edges))){a<-membership[edges$i[k]];b<-membership[edges$j[k]];membership[membership==b]<-a}
  tab<-table(membership);m$candidate_component<-ifelse(tab[as.character(membership)]>1,paste0('C',match(membership,unique(membership))),'Singleton')
  csv(m[,c('isolate','hospital','ward','date','outbreak','candidate_component')],'results/candidate_cluster_membership.csv')
  # Paired method contrasts use the same population, observation and threshold.
}
make_figures<-function(out,summary,cfg){
  library(ggplot2)
  p<-out$performance
  z<-summary[summary$threshold==cfg$primary_threshold&summary$metric%in%c('precision','sensitivity_all'),]
  z$scenario<-factor(labels[z$scenario],levels=labels);z$metric<-ifelse(z$metric=='precision','Precision of candidate links','Sensitivity across all true links')
  colours<-c('Genomic only'='#b75b25','Genomic + epidemiology'='#236c83')
  g<-ggplot(z,aes(scenario,mean,color=method))+geom_errorbar(aes(ymin=pmax(0,lower),ymax=pmin(1,upper)),position=position_dodge(.45),width=.15)+geom_point(position=position_dodge(.45),size=2.6)+facet_wrap(~metric)+scale_color_manual(values=colours)+scale_y_continuous(limits=c(0,1))+labs(title='Accuracy of candidate outbreak links',subtitle='100 paired repetitions; toy cutoff of five SNPs and at least 90% jointly callable core sites',x=NULL,y='Mean and 95% Monte Carlo interval',color=NULL,caption='All records and haplotypes are synthetic. Unobserved or QC-failed true links count as missed.\nIntervals describe simulation precision, not uncertainty in real hospital performance.')+theme_study()+theme(axis.text.x=element_text(angle=25,hjust=1))
  save_plot(g,'link_accuracy',11,6)
  z<-summary[summary$metric%in%c('precision','sensitivity_all'),];z$scenario_label<-labels[z$scenario];z$metric<-ifelse(z$metric=='precision','Precision','Sensitivity across all true links')
  g<-ggplot(z,aes(threshold,mean,color=method,linetype=metric))+geom_line()+geom_point()+facet_wrap(~scenario_label,ncol=2)+scale_color_manual(values=colours)+scale_x_continuous(breaks=cfg$thresholds)+scale_y_continuous(limits=c(0,1))+labs(title='Prespecified SNP threshold sensitivity',subtitle='The same thresholds are evaluated in every scenario; none is selected after seeing truth',x='Toy core SNP cutoff',y='Mean performance',color=NULL,linetype=NULL,caption='The 6,000-site model does not define a clinical K. pneumoniae cutoff. Candidate links do not prove transmission.')+theme_study()
  save_plot(g,'threshold_sensitivity',10,7)
  m<-out$example$metadata;m$label<-ifelse(m$outbreak>0,paste('Planted cluster',m$outbreak),'Background')
  g<-ggplot(m,aes(as.Date(date),interaction(hospital,ward,sep=' / '),color=label))+geom_point(position=position_jitter(height=.12,width=0,seed=1),size=2,alpha=.85)+scale_color_manual(values=c(Background='grey65','Planted cluster 1'='#2166ac','Planted cluster 2'='#b35806','Planted cluster 3'='#6b4c9a'))+labs(title='Synthetic hospital and ward surveillance timeline',subtitle='Prespecified first replicate: 150 isolates, three hospitals and three planted clusters',x='Collection date',y=NULL,color=NULL,caption='One isolate per synthetic person. These are isolate counts, not infection incidence or evidence from NICD.')+theme_study()
  save_plot(g,'surveillance_timeline',11,6)
  z<-read.csv('results/resistance_proportions.csv')
  g<-ggplot(z,aes(hospital,proportion))+geom_errorbar(aes(ymin=lower,ymax=upper),width=.12,color='#236c83')+geom_point(size=3,color='#236c83')+scale_y_continuous(limits=c(0,1),labels=function(x)paste0(round(x*100),'%'))+geom_text(aes(label=paste0(resistant,'/',n)),vjust=-1.2)+labs(title='Simulated carbapenem non-susceptibility among sampled isolates',subtitle='Prespecified first replicate; Wilson 95% binomial intervals',x=NULL,y='Proportion of sampled isolates',caption='Binary phenotypes are generated probabilities, not measured MICs or EUCAST classifications.\nHospital differences reflect assumptions and isolate sampling; patient-day denominators are unavailable.')+theme_study()
  save_plot(g,'resistance_proportions',8,5)
  tree<-ape::read.tree('results/display_synthetic_tree.nwk');nt<-length(tree$tip.label);xx<-ape::node.depth.edgelength(tree);yy<-rep(NA_real_,length(xx));yy[seq_len(nt)]<-seq_len(nt)
  ypos<-function(node){if(is.na(yy[node]))yy[node]<<-mean(vapply(tree$edge[tree$edge[,1]==node,2],ypos,numeric(1)));yy[node]};invisible(ypos(nt+1))
  ed<-data.frame(x=xx[tree$edge[,1]],xend=xx[tree$edge[,2]],y=yy[tree$edge[,2]],yend=yy[tree$edge[,2]])
  vertical<-do.call(rbind,lapply(unique(tree$edge[,1]),function(node){j<-tree$edge[tree$edge[,1]==node,2];data.frame(x=xx[node],y=min(yy[j]),yend=max(yy[j]))}))
  tips<-data.frame(x=xx[seq_len(nt)],y=yy[seq_len(nt)],isolate=tree$tip.label);tips$cluster<-m$label[match(tips$isolate,m$isolate)]
  g<-ggplot()+geom_segment(data=ed,aes(x=x,xend=xend,y=y,yend=yend),linewidth=.4,color='grey40')+geom_segment(data=vertical,aes(x=x,xend=x,y=y,yend=yend),linewidth=.4,color='grey40')+geom_point(data=tips,aes(x,y,color=cluster),size=2)+geom_text(data=tips,aes(x=x+2,y=y,label=isolate),hjust=0,size=2.5)+scale_color_manual(values=c(Background='grey65','Planted cluster 1'='#2166ac','Planted cluster 2'='#b35806','Planted cluster 3'='#6b4c9a'))+scale_x_continuous(expand=expansion(mult=c(.02,.2)))+scale_y_continuous(breaks=NULL)+labs(title='Genomic similarity and epidemiological linkage answer different questions',subtitle='36 prespecified tips pruned from the 150-isolate neighbour-joining tree',x='SNP-distance branch length across the synthetic core',y=NULL,color=NULL,caption='Close background relatives are deliberately included. Arbitrary display root; no dates, direction of spread or branch support inferred.\nThis tree uses complete simulated truth. A separate 12-sample read-processing test assesses observed sequences.')+theme_study()+theme(panel.grid=element_blank())
  save_plot(g,'synthetic_phylogeny',10,10)
  ids<-c(1:36,37:42);d<-out$example$distance[ids,ids];grid<-expand.grid(i=seq_along(ids),j=seq_along(ids));grid$distance<-as.vector(d)
  g<-ggplot(grid,aes(i,j,fill=distance))+geom_tile()+scale_fill_gradientn(colours=c('#f7fbff','#6baed6','#08306b'),trans='sqrt',breaks=c(0,25,75,125),guide=guide_colourbar(barwidth=grid::unit(8,'cm')))+coord_equal()+scale_x_continuous(breaks=c(6.5,18.5,30.5,39.5),labels=c('Cluster 1','Cluster 2','Cluster 3','Close background'))+scale_y_continuous(breaks=c(6.5,18.5,30.5,39.5),labels=c('Cluster 1','Cluster 2','Cluster 3','Close background'))+labs(title='True core SNP distances and close background relatives',subtitle='Three planted clusters and six close background relatives',x=NULL,y=NULL,fill='SNPs',caption='Small SNP distances can occur outside planted outbreaks.\nSquare-root colour scale; all sequences are synthetic.')+theme_study()+theme(panel.grid=element_blank())
  save_plot(g,'snp_heatmap',10,8)
  paired<-merge(p[p$method=='Genomic only',],p[p$method=='Genomic + epidemiology',],by=c('replicate','scenario','threshold'),suffixes=c('_genomic','_integrated'))
  pd<-do.call(rbind,lapply(c('precision','sensitivity_all'),function(metric)data.frame(replicate=paired$replicate,scenario=paired$scenario,threshold=paired$threshold,metric,value=paired[[paste0(metric,'_integrated')]]-paired[[paste0(metric,'_genomic')]])))
  csv(pd,'results/paired_method_differences.csv');ps<-mc(pd,c('scenario','threshold','metric'));csv(ps,'results/paired_method_summary.csv')
  z<-ps[ps$threshold==cfg$primary_threshold,];z$scenario<-labels[z$scenario];z$metric<-ifelse(z$metric=='precision','Precision gain','Sensitivity difference')
  g<-ggplot(z,aes(mean,scenario))+geom_vline(xintercept=0,linetype=2,color='grey60')+geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=.12,color='#236c83')+geom_point(size=2.5,color='#236c83')+facet_wrap(~metric,scales='free_x')+labs(title='The benefit and cost of adding epidemiological evidence',subtitle='Within-replicate differences at the prespecified five-SNP cutoff',x='Integrated minus genomic-only performance',y=NULL,caption='Means and 95% t Monte Carlo intervals over 100 paired populations. Missing metadata can improve specificity while losing true links.')+theme_study()
  save_plot(g,'paired_comparison',11,5)
}
