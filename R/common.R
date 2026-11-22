options(stringsAsFactors=FALSE)
csv <- function(x,path) write.csv(x,path,row.names=FALSE,na='')
wilson <- function(k,n) {
  if(!n)return(c(NA_real_,NA_real_))
  z<-qnorm(.975);p<-k/n;a<-(p+z*z/(2*n))/(1+z*z/n);b<-z*sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n);c(a-b,a+b)
}
mc <- function(d,groups,value='value') {
  out<-lapply(split(d,interaction(d[,groups,drop=FALSE],drop=TRUE)),function(z){v<-z[[value]];v<-v[is.finite(v)];n<-length(v);a<-if(n)mean(v) else NA_real_;se<-if(n>1)sd(v)/sqrt(n) else NA_real_;cbind(z[1,groups,drop=FALSE],n=n,mean=a,lower=a-qt(.975,max(1,n-1))*se,upper=a+qt(.975,max(1,n-1))*se)})
  do.call(rbind,out)
}
labels<-c(baseline='Good quality',poor_quality='Poor sequence quality',incomplete_sampling='60% isolate sampling',missing_metadata='25% missing ward data')
theme_study<-function()ggplot2::theme_minimal(base_size=12)+ggplot2::theme(plot.title=ggplot2::element_text(face='bold',size=15),panel.grid.minor=ggplot2::element_blank(),legend.position='bottom',plot.caption=ggplot2::element_text(hjust=0,size=10),plot.margin=ggplot2::margin(12,14,12,12))
save_plot<-function(p,name,w=10,h=6){ggplot2::ggsave(paste0('figures/',name,'.png'),p,width=w,height=h,dpi=300,bg='white');ggplot2::ggsave(paste0('figures/',name,'.pdf'),p,width=w,height=h,device='pdf')}
