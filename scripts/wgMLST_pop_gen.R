# Table 2 - Summary of genetic variations in AMOVA, estimated percent ФGT statistic 
# of population differentiation and p-values (in parentheses) based on a 999-permutation test
# between food pairs

# Amo1 is for the upper part of the table;
# Mt is for the lower part of the table;
# Both sub-tables indicate that large genetic variation within foods relative to that between foods


library(cluster)
library(adegenet)
library(genetics)
library(hierfstat)
library(pegas)
library(poppr)
library(tidyverse)
library(glmmBUGS)
library(xtable)
library(mmod)

# isolates with confirmed food source
lm_dat <- readRDS("C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/data/isolates_original_plus_new_dec_1_2021.rds") 

### cgMLST loci
cgmlst_loci <- read.csv("C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/data/cgMLST_loci.csv") %>% names

sel_rep_iso<-function(dat,ht){
  # genes<-dat %>% select(starts_with('LMO'))
  genes<-dat %>% select(all_of(cgmlst_loci))  ###### select core genes ######
  genes<-data.frame(sapply(genes,as.factor))
  epi<-dat %>% select(-starts_with('LMO'))
  
  suppressWarnings(cl<-daisy(genes,metric='gower'))
  hc<-hclust(cl,method='complete') 
  hc$labels<-epi$SRR_ID
  
  ct<-cutree(hc,h=ht)
  
  ff<-as.character(dat$SourceSite)
  ifsac<-as.character(epi$food)
  ct.mem<-lapply(unique(ct),function(x){
    w<-which(match(ct,x)==1)
    nam<-names(ct)[w]
    ss<-ff[!is.na(match(dat$SRR_ID, nam))]
    ic<-ifsac[match(nam,epi$SRR_ID)]
    if(length(nam)>1) {
      sel.nam<-names(table(ic))
      if(length(sel.nam)>1) {
        id<-sapply(sel.nam,function(y) {
          y.id<-nam[ic==y]
          set.seed(34)
          sample(y.id,1) 
        })}
      else {set.seed(23)
        id<-sample(nam[which(ic%in%sel.nam)],1)}
    } else id<-nam
    list(mem.id=nam,sourcesite=ss,ifsac=ic,sel.id=id)
  })
  sapply(ct.mem,function(x) length(x$mem.id))
  
  cluster.id<-unlist(sapply(ct.mem,function(x) x$sel.id))
  
  return(cluster.id)
}



################################################################
#  Deduplicate isolates based on core gene before imputation
################################################################

ht<-0.004 ##### proportional difference within clusters defined by ht
si<-sel_rep_iso(lm_dat,ht) ### select representative isolates


temp <-lm_dat %>% 
  filter(SRR_ID %in% si) %>% 
  select('food', starts_with("LMO")) %>% 
  mutate_if(is.integer, coalesce, 0L) %>% # integer "LMOxxxxx" as integer (52.59%) performs similar to "LMOxxxxx" as factor (51.85%)
  mutate(across(starts_with("LMO"), ~as.factor(as.character(.x)))) %>%
  as.data.frame


wgmlst=temp %>% select(starts_with('LMO'))
fd=levels(temp$food)
n.f=length(table(temp$food))

df.gid=df2genind(wgmlst,ploidy = 1)

#### calculate amova by pegas package ##########
strata(df.gid)=data.frame(temp$food)
diff_st=diff_stats(df.gid)
dd=dist(df.gid)
stra=strata(df.gid)
amo=pegas::amova(dd ~ temp.food, data = stra, nperm = 999)
amo
amo1=poppr.amova(df.gid,~temp.food,nperm=999)
amo1

amo0=ade4::amova(data.frame(df.gid@tab),dd,stra)
########## calculate pairwise amova ####################

fd.gene=as.data.frame(cbind(temp$food,wgmlst));names(fd.gene)[1]='food'
amo=lapply(1:(n.f-1),function(i) lapply((i+1):n.f,function(j) {
  dat=fd.gene[fd.gene$food %in% c(fd[i],fd[j]),]
  strata.df=data.frame(food=dat$food)
  df.g=df2genind(dat[,!names(dat)=='food'],ploidy = 1)
  strata(df.g)=strata.df
  g=poppr.amova(df.g,~food,nperm=999)
  p=randtest(g,nrepet = 999)
  
  #d=dist(df.g,method = 'binary')
  #s=strata(df.g)
  #g1=pegas::amova(d ~ food, data = s, nperm = 99)
  
  list(fd_pair=paste(fd[i],fd[j],sep='-'),amova=g,pvalue=p$pvalue)
  }))

phi=lapply(amo,function(p) unlist(lapply(p, function(g) round(g$amova$statphi*100,1))))
pvalue=lapply(amo,function(x) unlist(lapply(x, function(g) g$pvalue)))
mt=data.frame(matrix(rep(NA,length(fd)^2),ncol=length(fd)))
lapply(1:length(phi),function(p){
  mt[p,(length(fd)-length(phi[[p]])+1):length(fd)]<<-paste(paste(phi[[p]],pvalue[[p]],sep='('),')',sep='')})
names(mt)=fd
row.names(mt)=fd

print(xtable(mt),'C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/amova_tab.html',type='html',include.rownames=F)

print(xtable(cbind(amo1$componentsofcovariance,amo1$results$Df),digits = 1),'\\\\cdc.gov\\private\\M318\\vhg8\\whole genomo sequence\\data\\new data\\amova_tab0.html',type='html',include.rownames=T)

(fd=sapply(amo,function(p) unlist(lapply(p, function(g) g$fd_pair))))
(pvalues=sapply(amo,function(p) unlist(lapply(p, function(g) g$pvalue))))

# save(amo,amo1,df.gid, fd, pvalues, file='C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/pop genetic.RData')
load(file='C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/pop genetic.RData')

tab=data.frame(df.gid@tab)
##################### calculate MSE for each food ###################################
library(ICSNP)
sapply(unique(stra[,1]),function(f){
  temp=tab[stra$temp.food==f,]
  d=daisy(temp,metric = 'gower')
  d2=d^2
})

d=dist(as.matrix(tab))
d=quasieuclid(d)
counts=data.frame(sapply(fd,function(f){
  d=tab[temp$food %in% f,]
  sapply(d,function(x) sum(x,na.rm=T))
}))
amo1=ade4:::amova(counts,d)

#df.gid@pop=temp$food
#pair.pd=pairwise_Gst_Nei(df.gid, linearized = FALSE)
#pair.pd

#s.meir=Phi_st_Meirmans(df.gid)
#s.meir
#bs <- chao_bootstrap(df.gid, nreps = 100)
#summarise_bootstrap(bs, D_Jost)     # for Nei's Gst
