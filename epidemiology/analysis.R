source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'helpers.R'),encoding='UTF-8')
m<-read.csv('data/synthetic_isolates.csv');stopifnot(nrow(m)==150,!anyDuplicated(m$isolate),all(m$phenotype %in% 0:1),all(m$death %in% 0:1),!anyNA(m))
m$hospital<-factor(m$hospital);m$ward<-factor(m$ward);m$date<-as.Date(m$date)
checks$people<-150;checks$deaths<-sum(m$death);checks$resistant_phenotypes<-sum(m$phenotype);checks$raw_data_regenerated<-FALSE
write.csv(data.frame(variable=names(m),missing=vapply(m,function(x)sum(is.na(x)),integer(1))),file.path(ep,'missingness.csv'),row.names=FALSE)
base<-m[c('hospital','age','ward','severity','gene','phenotype','death')]
for(v in c('severity','gene','phenotype','death'))base[[v]]<-factor(base[[v]],0:1,c('Absent','Present'))
tab1<-gtsummary::tbl_summary(base,by=hospital,missing='no',statistic=gtsummary::all_continuous()~'{mean} ({sd})',label=list(age~'Age in years',ward~'Ward category',severity~'Severe presentation',gene~'AMR gene carriage',phenotype~'Generated resistant phenotype',death~'Generated death')) |> gtsummary::add_overall() |> gtsummary::bold_labels()
save_table('population',tab1,'Table 1 Characteristics of the 150 synthetic isolate associated people')
res<-read.csv('results/resistance_proportions.csv');save_table('hospital_resistance',data.frame(Hospital=res$hospital,N=res$n,Resistant=res$resistant,`Percent resistant (95% CI)`=fmtci(100*res$proportion,100*res$lower,100*res$upper,1),check.names=FALSE),'Table 2 Generated resistance among sampled isolates by hospital')
ct<-as.data.frame.matrix(table(factor(m$gene,0:1),factor(m$phenotype,0:1)));save_table('concordance',data.frame(`AMR gene carriage`=c('Absent','Present'),`Susceptible phenotype`=ct[,1],`Resistant phenotype`=ct[,2],check.names=FALSE),'Table 3 Gene carriage and generated phenotypic resistance')
fit<-glm(phenotype~gene+hospital,data=m,family=binomial());stopifnot(fit$converged,fit$rank==length(coef(fit)),!separation(fit))
or<-tidy_hc3(fit,exponentiate=TRUE);write.csv(or,file.path(ep,'resistance_associations.csv'),row.names=FALSE)
gt<-gtsummary::tbl_regression(fit,exponentiate=TRUE,tidy_fun=tidy_hc3,label=list(gene~'AMR gene carriage present vs absent',hospital~'Hospital'),pvalue_fun=gtsummary::label_style_pvalue(digits=3)) |> gtsummary::bold_labels()
save_table('associations',gt,'Table 4 Hospital adjusted resistance associations with HC3 confidence intervals')
death<-read.csv('results/mortality_associations.csv');death_refit<-glm(death~gene+age+severity+hospital,data=m,family=binomial());stopifnot(max(abs(exp(coef(death_refit))-death$OR))<1e-8)
save_table('mortality',data.frame(Term=death$term[-1],`Adjusted OR (95% CI)`=fmtci(death$OR[-1],death$lower[-1],death$upper[-1]),P=fmtp(death$p[-1]),check.names=FALSE),'Table 5 Preserved exploratory mortality model with conventional Wald intervals')
infl<-data.frame(isolate=m$isolate,cooks=cooks.distance(fit),leverage=hatvalues(fit));write.csv(infl,file.path(ep,'influence.csv'),row.names=FALSE)
checks$resistance_model<-list(converged=TRUE,separation=FALSE,rank=fit$rank,influential_above_4n=sum(infl$cooks>4/nrow(m)),all_records_retained=TRUE);checks$legacy_mortality_OR_reproduced<-TRUE
pub<-jsonlite::fromJSON('case_study/results/descriptive_summary.json')
save_table('public',data.frame(Measure=c('Public outbreak assemblies','ST258 assignments','KPC family positive assemblies','Common callable chromosome bases','Variable common core sites','Pairwise SNP median','Pairwise SNP range'),Value=c(pub$n,pub$ST258,pub$KPC_positive,pub$common_core_bases,pub$variable_sites,pub$pairwise_snp_median,paste(pub$pairwise_snp_range,collapse=' to '))),'Table 6 Public NIH outbreak subset genomic findings')
month<-as.data.frame(table(month=format(m$date,'%Y-%m'),hospital=m$hospital));month$month<-as.Date(paste0(month$month,'-01'))
f1<-save_plot('epidemic_curve',ggplot(month,aes(month,Freq,fill=hospital))+geom_col(position='stack')+scale_fill_manual(values=c('#167D9A','#B56528','#667C54'))+labs(x='Collection month in 2025',y='Sampled isolate records',fill='Hospital',title='Distribution of sampled isolates over time',subtitle='150 synthetic records; counts are not infection incidence'),8,4)
f2<-save_plot('hospital_resistance',ggplot(res,aes(hospital,proportion))+geom_pointrange(aes(ymin=lower,ymax=upper),colour='#167D9A')+scale_y_continuous(labels=scales::label_percent(),limits=c(0,1))+labs(x=NULL,y='Generated resistant phenotype proportion',title='Resistance distribution across synthetic hospitals',subtitle='Wilson 95% intervals; selected isolates rather than an at risk population'),8,4)
z<-or[or$term!='(Intercept)',];z$label<-c(gene='AMR gene carriage',`hospitalHospital B`='Hospital B vs A',`hospitalHospital C`='Hospital C vs A')[z$term]
f3<-save_plot('adjusted_associations',ggplot(z,aes(estimate,reorder(label,estimate)))+geom_vline(xintercept=1,linetype=2)+geom_errorbar(aes(xmin=conf.low,xmax=conf.high),orientation='y',width=.15,colour='#167D9A')+geom_point(size=2.5,colour='#167D9A')+scale_x_log10()+labs(x='Adjusted odds ratio with HC3 95% CI',y=NULL,title='Gene carriage and hospital associations with resistance',subtitle='All 150 synthetic records retained'),8,3.8)
h<-or[or$term=='gene',];mort<-death[death$term=='gene',]
body<-paste0(r"(## Abstract

**Background.** Hospital surveillance asks how resistance is distributed across sampled patients and whether genetic determinants align with phenotypes and clinical outcomes. **Aim.** To describe resistance, genomic relatedness and exploratory clinical associations in a synthetic hospital cohort, with a separate descriptive reanalysis of public outbreak assemblies. **Methods.** R analysed one saved cohort of 150 people with one isolate each. We summarised hospital and time distributions, cross-tabulated gene carriage and resistance, and fitted hospital-adjusted logistic regression. A published 12-assembly NIH Klebsiella pneumoniae subset supplied a separate measured-genome case study.

)",sprintf('**Results.** There were %d generated resistant phenotypes and %d generated deaths. The hospital-adjusted gene-carriage OR for resistance was %s. All 12 public assemblies were ST258 and carried KPC-family genes. **Conclusion.** The synthetic cohort illustrates heterogeneous resistance and gene-phenotype associations; the selected public assemblies describe a closely related resistant lineage. Neither analysis establishes hospital infection incidence or direct transmission.\n\n',sum(m$phenotype),sum(m$death),fmtci(h$estimate,h$conf.low,h$conf.high)),r"(## Introduction

Resistance surveillance combines the distribution of laboratory findings with the setting and timing of patient care. Genomic relatedness can support an outbreak investigation, but similar genomes alone do not identify who infected whom. The present analysis centres on resistance distributions and clinical associations. The earlier candidate-link benchmarks and sequencing execution tests remain supplementary technical evidence, preserved without changing their results.

### Aim and objectives

The primary aim was to characterise generated phenotypic resistance among sampled hospital isolates and its association with AMR gene carriage after hospital adjustment. Objectives were to describe the cohort by hospital and sampling time; compare gene carriage with phenotype; examine the previously fitted mortality association as exploratory; and describe lineage, resistance determinants and genomic distances in the separate public outbreak subset.

## Methods

### Study populations and preprocessing

The primary analytical dataset is the saved first synthetic population: 150 distinct people, three fictional hospitals and three ward categories, with one isolate per person during 2025. All participant, gene, phenotype and mortality values are simulated. We verified unique identifiers, finite complete analysis variables and valid binary coding. No record or outcome was regenerated, imputed or excluded. Three planted clusters share hospital, ward and collection windows; random core haplotypes approximate genomic relatedness and are not reconstructed K. pneumoniae genomes.

The separate public component consists of 12 accession-versioned NIH outbreak assemblies, a convenience subset of the original publication. It contains measured bacterial sequence data but does not supply an independent patient cohort, measured phenotypic susceptibility results or complete patient contact histories. It is not pooled with the simulated cohort. Assembly QC, exact seven-locus MLST, AMRFinderPlus, masked reference projection and full common-core maximum-likelihood phylogenetics were completed previously; their inputs, versions and original outputs remain unchanged.

### Outcomes and estimands

The primary outcome is the binary generated resistant phenotype, coded 1 versus 0. Gene carriage is a binary simulated marker that receives blaKPC-2 or blaNDM-1 labels when present; those labels do not constitute biological sequencing measurements in this cohort. The secondary clinical outcome is generated death, with no fixed follow-up duration supplied. It is therefore analysed as a binary outcome rather than a survival endpoint. Hospital proportions describe selected isolate records. Patient admissions and patient-day denominators are unavailable, so infection incidence and hospital-wide prevalence cannot be estimated.

### Statistical analysis

Descriptive statistics use counts, proportions, age mean and SD, and Wilson 95% intervals for resistance proportions. With resistant outcome R, gene indicator G and hospital indicators B and C, the primary model was

$$\operatorname{logit}\Pr(R_i=1)=\alpha+\beta_GG_i+\beta_BB_i+\beta_CC_i.$$

Odds ratios were exp(beta), with HC3 sandwich Wald 95% confidence intervals. The model assumes an adequate conditional logit relationship and independent people conditional on the included covariates. HC3 addresses leverage-dependent variance estimation; it does not remove unmeasured confounding or validate the synthetic mechanism. Convergence, design rank, separation and influence were checked, with all records retained. No stepwise selection was used.

The preserved mortality model adjusted for gene carriage, age, severity and hospital. Its original conventional Wald intervals are labelled separately from HC3 intervals. Only 21 deaths support that exploratory model, so estimates may be unstable. We independently refitted it to verify the saved ORs without replacing the original result file. All association tests are exploratory, not confirmatory biomarker validation. Public-genome results are descriptive; related isolates and their pairwise SNP distances are not independent population observations.

## Results

### Population and time distribution

The cohort table describes the denominator for each hospital comparison. Differences are properties of the generated sample and were not used to select model covariates.

)",T('population'),r"(

The epidemic curve shows when isolate records were sampled and how their hospital distribution changed over the year. Peaks are partly imposed by the planted cluster windows. The plot does not estimate a disease rate or establish hospital transmission.

)",F(f1,'Figure 1 Sampling distribution by hospital and calendar month'),r"(

### Resistance and gene phenotype associations

Hospital-specific resistant proportions vary in this selected synthetic cohort. Wilson intervals describe binomial sampling uncertainty under independent observations; they do not account for a representative hospital surveillance design or residual cluster dependence.

)",T('hospital_resistance'),'\n\n',F(f2,'Figure 2 Resistance proportions among the sampled synthetic isolates'),r"(

The cross-tabulation distinguishes gene carriage from the generated phenotype. Discordance is possible because phenotype was generated probabilistically. A detected resistance gene is therefore not presented as a measured susceptibility test.

)",T('concordance'),'\n\n',sprintf('After hospital adjustment, gene carriage was associated with higher generated resistance odds: OR %s, p %s. This association reflects the simulation mechanism and does not validate a clinical test.\n\n',fmtci(h$estimate,h$conf.low,h$conf.high),fmtp(h$p.value)),T('associations'),'\n\n',F(f3,'Figure 3 Hospital adjusted associations with generated resistance'),r"(

### Exploratory mortality association

)",sprintf('The previously fitted mortality association for gene carriage was OR %s, p %s. With 21 deaths and several adjustment parameters, this result warrants particular caution. It describes a generated association and does not establish that resistance caused death.\n\n',fmtci(mort$OR,mort$lower,mort$upper),fmtp(mort$p)),T('mortality'),r"(

### Public outbreak genome findings

All 12 selected public assemblies belonged to ST258 and carried blaKPC-3. The common callable chromosome contained 42 variable sites, with pairwise distances from 0 to 25 SNPs. These results support the presence of a closely related resistant lineage in this published outbreak subset. They do not identify individual transmission events or imply that the subset represents circulating organisms across hospitals.

)",T('public'),'\n\n',F('case_study/figures/public_phylogeny.png','Figure 4 Maximum likelihood phylogeny of the public NIH outbreak subset'),r"(

## Discussion

The epidemiological findings concern resistance distribution, gene-phenotype association and the limits of interpreting genomic relatedness. In the synthetic cohort, gene carriage and hospital context describe different aspects of the observed resistance pattern. The separate public analysis supplies measured genomic evidence of an established outbreak lineage. Its concordance with the original publication is supporting context, rather than a new validation cohort.

The main limitations are generated phenotypes, a selected sample with no infection-free comparison population, limited mortality events, simplified genomic ancestry and absent patient-day denominators. Hospital differences may reflect case mix and ascertainment. A phylogeny alone cannot establish acquisition location or direction of spread. No phenotype prediction from the public assemblies, clinical threshold recommendation or causal effect is claimed.

## Conclusion

This study presents an applied hospital resistance analysis with explicit denominators, adjusted associations and genomic context. The numerical associations are entirely synthetic; the public assemblies provide a separate descriptive outbreak reanalysis. The earlier accuracy benchmarks remain available as supplementary methods results.

## Reproducibility and supporting analyses

Run `Rscript --vanilla epidemiology/analysis.R --report` from the repository root after installing the packages listed in `epidemiology/R_packages.csv` and Quarto. This command reuses saved data, creates tables and figures, checks original-file hashes and renders the Word report. Sequencing and original simulation workflows remain available through `main.nf`, `case_study/main.nf` and `RUN_WORKFLOW.txt`. Their computational requirements and results are documented in the archived manuscript. The analysis assumes no NICD affiliation or endorsement.

The archived README retains the complete original candidate-link benchmark, sequencing validation, database versions and publication references. Public genome provenance is retained in `case_study/provenance`. Data, estimates and figures that predate this refocus are protected by `epidemiology/original_results_manifest.json`.

## References

1. Snitkin ES and colleagues. Tracking a hospital outbreak of carbapenem resistant Klebsiella pneumoniae with whole genome sequencing. Science Translational Medicine 2012;4:148ra116. [Original outbreak study](https://doi.org/10.1126/scitranslmed.3004129).
2. NICD Centre for Healthcare Associated Infections Antimicrobial Resistance and Mycoses. [Surveillance context](https://www.nicd.ac.za/centres/centre-for-healthcare-associated-infections-antimicrobial-resistance-and-mycoses/).
)")
write_study('Resistance distribution and genomic relatedness in hospital surveillance','Synthetic hospital cohort with a separate measured public genome case study',body,'Hospital_AMR_Report')
