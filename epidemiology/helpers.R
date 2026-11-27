if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
required<-c('ggplot2','gtsummary','flextable','officer','jsonlite','dplyr','sandwich','lpSolve','digest','ape','knitr','rmarkdown')
missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Missing R packages: ',paste(missing,collapse=', '))
library(ggplot2)
entry<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
root<-normalizePath(file.path(dirname(entry),'..'),winslash='/',mustWork=TRUE);setwd(root)
ep<-file.path(root,'epidemiology');tabdir<-file.path(ep,'tables');figdir<-file.path(ep,'figures')
for(p in c(ep,tabdir,figdir))dir.create(p,recursive=TRUE,showWarnings=FALSE)
table_registry<-list();figure_registry<-list();checks<-list()
fmtp<-function(p)ifelse(is.na(p),'',ifelse(p<.001,'<0.001',sprintf('%.3f',p)))
fmtci<-function(e,l,u,d=2)sprintf(paste0('%.',d,'f (%.',d,'f, %.',d,'f)'),e,l,u)
wilson<-function(k,n){p<-k/n;z<-qnorm(.975);den<-1+z*z/n;centre<-(p+z*z/(2*n))/den;half<-z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den;cbind(lower=centre-half,upper=centre+half)}
tidy_hc3<-function(x,conf.int=TRUE,conf.level=.95,exponentiate=FALSE,...){b<-coef(x);se<-sqrt(diag(sandwich::vcovHC(x,type='HC3')));q<-qnorm((1+conf.level)/2);z<-data.frame(term=names(b),estimate=b,std.error=se,statistic=b/se,p.value=2*pnorm(-abs(b/se)),conf.low=b-q*se,conf.high=b+q*se,row.names=NULL);if(exponentiate)z[c('estimate','conf.low','conf.high')]<-lapply(z[c('estimate','conf.low','conf.high')],exp);z}
separation<-function(m){X<-model.matrix(m);y<-model.response(model.frame(m));B<-cbind(X*(2*y-1),-X*(2*y-1));lpSolve::lp('min',rep(0,ncol(B)),rbind(B,colSums(B)),rep('>=',nrow(B)+1),c(rep(0,nrow(B)),1))$status==0}
save_table<-function(name,data,caption){
 z<-if(inherits(data,'gtsummary'))gtsummary::as_tibble(data,col_labels=TRUE) else data
 write.csv(z,file.path(tabdir,paste0(name,'.csv')),row.names=FALSE)
 ft<-if(inherits(data,'gtsummary'))gtsummary::as_flex_table(data) else flextable::flextable(z)
 ft<-flextable::font(ft,fontname='Calibri',part='all');ft<-flextable::fontsize(ft,size=9,part='all');ft<-flextable::padding(ft,padding=3,part='all');ft<-flextable::bg(ft,bg='#E6EEF2',part='header');ft<-flextable::color(ft,color='black',part='all');ft<-flextable::autofit(ft);ft<-flextable::set_table_properties(ft,layout='autofit',width=1)
 border<-officer::fp_border(color='#D9D9D9',width=.5);ft<-flextable::border_outer(ft,border=border);ft<-flextable::border_inner_h(ft,border=border);ft<-flextable::border_inner_v(ft,border=border);ft<-flextable::paginate(ft,init=TRUE,hdr_ftr=TRUE)
 saveRDS(ft,file.path(tabdir,paste0(name,'.rds')))
 cell<-function(x){x<-as.character(x);x[is.na(x)]<-'';x<-gsub('\n','<br>',x,fixed=TRUE);gsub('|','&#124;',x,fixed=TRUE)}
 z[]<-lapply(z,cell);md<-c(paste0('| ',paste(cell(names(z)),collapse=' | '),' |'),paste0('| ',paste(rep('---',ncol(z)),collapse=' | '),' |'),apply(z,1,function(r)paste0('| ',paste(r,collapse=' | '),' |')))
 writeLines(md,file.path(tabdir,paste0(name,'.md')),useBytes=TRUE);table_registry[[name]]<<-list(caption=caption,md=paste(md,collapse='\n'));invisible(data)
}
theme_set(theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),plot.subtitle=element_text(size=10)))
save_plot<-function(name,p,width=8,height=4.5){ggsave(file.path(figdir,paste0(name,'.png')),p,width=width,height=height,dpi=300,bg='white');ggsave(file.path(figdir,paste0(name,'.pdf')),p,width=width,height=height,device=cairo_pdf);paste0('epidemiology/figures/',name,'.png')}
T<-function(name)paste0('@@TABLE_',name,'@@')
F<-function(path,caption){name<-paste0('figure',length(figure_registry)+1);figure_registry[[name]]<<-list(path=path,caption=caption);paste0('@@FIG_',name,'@@')}
validate_sources<-function(){
 manifest<-jsonlite::fromJSON(file.path(ep,'original_results_manifest.json'),simplifyVector=FALSE)
 checked<-vapply(names(manifest$files),function(n){v<-manifest$files[[n]];hash<-digest::digest(file=file.path(root,n),algo='sha256');if(hash==v$sha256)return(TRUE);if(is.null(v$canonical_sha256))return(FALSE);bytes<-readBin(file.path(root,n),'raw',n=file.info(file.path(root,n))$size);digest::digest(charToRaw(gsub('\r\n','\n',rawToChar(bytes),fixed=TRUE)),algo='sha256',serialize=FALSE)==v$canonical_sha256},logical(1))
 stopifnot(all(checked));length(checked)
}
write_study<-function(title,status,body,report_name,qmd_name=report_name,word_path=paste0('report/',report_name,'.docx')){
 readme<-body;qmd<-body
 for(n in names(table_registry)){
  x<-table_registry[[n]];readme<-gsub(T(n),paste0('**',x$caption,'**\n\n',x$md),readme,fixed=TRUE)
  chunk<-paste0('::: {custom-style="Table Caption"}\n',x$caption,'\n:::\n\n```{r}\n#| echo: false\nreadRDS("../epidemiology/tables/',n,'.rds")\n```')
  qmd<-gsub(T(n),chunk,qmd,fixed=TRUE)
 }
 for(n in names(figure_registry)){
  x<-figure_registry[[n]];stopifnot(file.exists(x$path));marker<-paste0('@@FIG_',n,'@@')
  readme<-gsub(marker,paste0('![',x$caption,'](',x$path,')'),readme,fixed=TRUE)
  qmd<-gsub(marker,paste0('![',x$caption,'](../',x$path,'){width=6.3in}'),qmd,fixed=TRUE)
 }
 stopifnot(!grepl('@@',readme),!grepl('@@',qmd))
 intro<-paste0('# ',title,'\n\n**Vusumuzi Mabasa | 27 November 2026**\n\n**',status,'**\n\n[Word report](',word_path,') | [PDF report](',sub('\\.docx$','.pdf',word_path),') | [R analysis](epidemiology/analysis.R) | [Original analyses and manuscript](archive/pre-epidemiology-refocus/README.md)\n\n')
 writeLines(paste0(intro,readme),'README.md',useBytes=TRUE)
 yaml<-paste0('---\ntitle: "',title,'"\nauthor: "Vusumuzi Mabasa"\ndate: "27 November 2026"\nsubtitle: "',status,'"\nformat:\n  docx:\n    reference-doc: epidemiology_reference.docx\n    toc: false\nexecute:\n  echo: false\n  warning: false\n  message: false\n---\n\n```{r}\nlibrary(flextable)\n```\n\n')
 qmd_path<-file.path(root,'report',paste0(qmd_name,'.qmd'));writeLines(paste0(yaml,qmd),qmd_path,useBytes=TRUE)
 preserved<-validate_sources()
 checks$original_scientific_files_preserved<-preserved;checks$tables<-length(table_registry);checks$figures<-length(figure_registry);checks$synthetic_and_public_data_distinguished<-TRUE
 jsonlite::write_json(c(list(status='passed'),checks),file.path(ep,'validation.json'),pretty=TRUE,auto_unbox=TRUE)
 write.csv(data.frame(package=required,version=vapply(required,function(n)as.character(packageVersion(n)),character(1))),file.path(ep,'R_packages.csv'),row.names=FALSE)
 capture.output(sessionInfo(),file=file.path(ep,'R_sessionInfo.txt'))
 if('--report' %in% commandArgs(TRUE)){
  quarto<-Sys.which('quarto');if(quarto=='')quarto<-'C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe'
  Sys.setenv(R_ENVIRON_USER=if(.Platform$OS.type=='windows')'NUL' else '/dev/null',R_PROFILE_USER=if(.Platform$OS.type=='windows')'NUL' else '/dev/null',QUARTO_R=file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript'))
  stopifnot(system2(quarto,c('render',shQuote(qmd_path)))==0)
  rendered<-sub('\\.qmd$','.docx',qmd_path);if(normalizePath(rendered,winslash='/')!=normalizePath(word_path,winslash='/',mustWork=FALSE))stopifnot(file.copy(rendered,word_path,overwrite=TRUE))
 }
 message('Epidemiological report complete; ',preserved,' original scientific files preserved')
}
