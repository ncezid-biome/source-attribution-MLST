
# sel_rep_iso() works well under R-3.6.1 but gives fatal errors under R-4.1.1

library(randomForestSRC) 
library(dendextend) 
library(cluster) 
library(adegenet) 
library(pegas) 
library(poppr) 
library(mmod) 
library(multiROC) 
library(phytools) 
library(data.table)
library(tidyverse)


# isolates with confirmed food source
lm_dat <- readRDS("data-raw/lm_wgs_manuscript/datasets/isolates_original_plus_new_dec_1_2021.rds") 
 
### cgMLST loci
cgmlst_loci <- read.csv("data-raw/lm_wgs_manuscript/datasets/cgMLST_loci.csv") %>% names


#################
#  Functions
#################

sel_rep_iso<-function(dat,ht){
  genes<-dat %>% select(all_of(cgmlst_loci))  ###### select core genes ######
  genes<-data.frame(sapply(genes,as.factor))
  epi<-dat %>% select(-starts_with('LMO')) # non-gene variables
  
  suppressWarnings(cl<-daisy(genes,metric='gower'))  # works well under R-3.6.1 but gives fatal errors under R-4.1.1
  hc<-hclust(cl,method='complete') 
  hc$labels<-epi$SRR_ID
  
  ct<-cutree(hc,h=ht)  # assign cluster number (1 to 320) to each isolate (SSR_ID) 
 
  ff<-as.character(dat$SourceSite)
  ifsac<-as.character(epi$food)
  ct.mem<-lapply(unique(ct),function(x){
    w<-which(match(ct,x)==1)  # position of isolates (SSR_ID) in a cluster
    nam<-names(ct)[w]         # SSR_ID of isolates in a cluster
    ss<-ff[!is.na(match(dat$SRR_ID, nam))] # SourceSites of isolates in a cluter
    ic<-ifsac[match(nam,epi$SRR_ID)]    # IFSAC food categories of isolates in a cluter
    if(length(nam)>1) {       ### if there are more than one SSR_ID in a cluster
      sel.nam<-names(table(ic))   # unique IFSAC food categories in a cluster
      if(length(sel.nam)>1) {       # if there is more than one IFSAC food category in a cluster, one SSR_ID was randomly selected from each food category
        id<-sapply(sel.nam,function(y) {  
          y.id<-nam[ic==y]
          set.seed(34)
          sample(y.id,1) 
        })}
      else {set.seed(23)            # if there is only one food category in a cluster, one SSR_ID was randomly selected from that cluster
        id<-sample(nam[which(ic%in%sel.nam)],1)}
    } else id<-nam            ### if there is only one SSR_ID in a cluster, that SSR_ID is selected
    list(mem.id=nam,sourcesite=ss,ifsac=ic,sel.id=id)
  })
  sapply(ct.mem,function(x) length(x$mem.id))  # number of isolates in each cluster

  cluster.id<-unlist(sapply(ct.mem,function(x) x$sel.id))   # SSR_ID of randomly selected isolates
  
  return(cluster.id)
}


sen_spe=function(pred,obs){
  tab=xtabs(~obs+pred)
  (correct.rate=sum(diag(tab))/sum(tab))
  mar1=apply(tab,1,sum)/sum(tab)
  mar2=apply(tab,2,sum)/sum(tab)
  exp.agg=sum(mar1*mar2)
  kappa.est=(correct.rate-exp.agg)/(1-exp.agg)
  kappa.est
  
  (t.pos=diag(tab))
  (f.neg=sapply(1:length(levels(obs)),function(x) sum(tab[x,])-tab[x,x]))
  (sen=round(t.pos/(t.pos+f.neg),2))
  
  (f.pos=sapply(1:length(levels(obs)),function(x) sum(t(tab)[x,])-t(tab)[x,x]))
  (t.neg=sapply(1:length(levels(obs)),function(x) sum(tab)-sum(tab[x,])))
  (spe=round(t.neg/(t.neg+f.pos),2))
  o.spe=sum(t.neg)/sum(t.neg+f.pos)
  
  over_accu=(sum(t.pos)+sum(t.neg))/(sum(t.pos)+sum(t.neg)+sum(f.pos)+sum(f.neg))
  prev=table(train.df$food)/nrow(train.df)
  dat=data.frame(tab,apply(tab,1,sum))
  return(list(dat=dat,se=c(sen,correct.rate),sp=c(spe,o.spe),accurate_rate=sen*prev+spe*(1-prev),over_accu=over_accu,kappa=kappa.est))
  
}


seq.rfsrc=function(dat,imp,my.seq,ntree){
  
  res=foreach(n=my.seq,.packages='randomForestSRC') %dopar% {
    set.seed(23)
    m=rfsrc(food~.,dat[,names(dat) %in% c('food',as.character(imp[1:n]))],ntree=1000)
    tab=table(m$yvar,m$class.oob)
    t.pos=diag(tab)
    f.neg=sapply(1:length(unique(dat$food)),function(x) sum(tab[x,])-tab[x,x])
    sen=c(round(t.pos/(t.pos+f.neg),2))
    over_sen=round(sum(t.pos)/sum(t.pos+f.neg),2)
    list(food=c(dimnames(tab)[[1]],'overall'),se=c(sen,over_sen))
  }
  
  dat=data.table(food=c(sapply(res,function(x) unlist(x$food))),
                 se=c(sapply(res,function(x) unlist(x$se))),loci=rep(my.seq,each=length(res[[1]]$food)))

  return(list(res=dat))
}


plot_panel_pred_prob_ind=function(pred,train=NULL){
  sp=split(pred,pred$pred.fd)
  g=lapply(1:length(sp),function(d){ #### reorder obs based on pred prob of THE food
    prob=unlist(tapply(sp[[d]]$prob,sp[[d]]$obs,function(x) rep(x[d],length(x))))  # replicate the largest probability for a given isolate X times (X is the number of food categories)
                                                                                   # As a result, in the data set for pred.fd == "fruit", the prob for fruit was replicated X times
    dat=sp[[d]][order(prob,decreasing=T),]
  })
  pred1=do.call(rbind,g)
  pred1$obs=factor(pred1$obs,levels=unique(pred1$obs)) #### keep the correct order 
  
  pred1$color=ifelse(pred1$pred.fd=='dairy','',ifelse(pred1$pred.fd=='fruit','grey',
                                                      ifelse(pred1$pred.fd=='meat','yellow',ifelse(pred1$pred.fd=='sea food','green','purple'))))
  
  pred1$color.obs=ifelse(pred1$obs.fd=='dairy','red',ifelse(pred1$obs.fd=='fruit','grey',
                                                            ifelse(pred1$obs.fd=='meat','yellow',ifelse(pred1$obs.fd=='sea food','green','purple'))))
  
  pred1$food=ifelse(pred1$food=='sea food','seafood',ifelse(pred1$food=='veg','vegetable',as.character(pred1$food)))
  
  pred1$pred.fd=ifelse(pred1$pred.fd=='sea food','seafood',ifelse(pred1$pred.fd=='veg','vegetable',as.character(pred1$pred.fd)))
  
  
  my.col=c('purple','orange','red','blue','green')
  p=ggplot(pred1,aes(obs,prob))+geom_bar(aes(fill=food),stat='identity',width = 1)+ylab('predicted probability of foods')+
    scale_fill_manual(values = my.col)+xlab('isolates')+
    theme(axis.line.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+theme(text = element_text(size=16),axis.line.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=10),
                                              axis.ticks.x=element_blank(),legend.position = c(0.8,0.3),legend.text = element_text(size=14),
                                              legend.title = element_text(size=14),strip.text.x = element_text(size = 15))+facet_wrap(~ pred.fd, scales='free_x')
  
  if(is.null(train)) p=p+geom_rug(aes(obs), sides="b")
  else p=p+geom_rug(aes(obs), subset(pred1,obs.fd==pred.fd), sides="b")
  return(list(p=p))
}


################################################################
#  Deduplicate isolates based on core gene before imputation
################################################################

ht<-0.004 ##### proportional difference within clusters defined by ht

si<-sel_rep_iso(lm_dat, ht) ### select representative isolates


###################################################################
#  importance of genes from random forest model based on all genes
###################################################################

train.df.all <- lm_dat %>% 
  filter(SRR_ID %in% si) %>%
  dplyr::select('food', starts_with("LMO")) %>% 
  mutate_if(is.integer, coalesce, 0L) %>% 
  mutate(across(starts_with("LMO"), ~as.factor(as.character(.x)))) %>%
  as.data.frame 

set.seed(23)
rf.all <- rfsrc(food ~., train.df.all, importance = T) 


################################################################################
#   Final model - Random forest model based on core gene + top 100 ranked genes
################################################################################

rf <- rf.all
c.imp=rf$importance[,1]
rk.genes=names(sort(c.imp,decreasing = T)) # top ranked genes

train.df.core_top100 <-lm_dat %>% 
  filter(SRR_ID %in% si) %>% 
  select('food', all_of(c(cgmlst_loci, rk.genes[1:100]))) %>% 
  mutate_if(is.integer, coalesce, 0L) %>% 
  mutate(across(starts_with("LMO"), ~as.factor(as.character(.x)))) %>%
  as.data.frame 

set.seed(23)
rf.core_top100 <- rfsrc(food ~., train.df.core_top100, importance = F) 


###################################################
#  Table 1 Summary of genetic variations in AMOVA
###################################################

temp <- train.df.all
wgmlst=temp %>% select(starts_with('LMO'))
fd=levels(temp$food)
n.f=length(table(temp$food))

# Table 1 - Genetic variance (σ)of food isolates and their contribution (%) to the total variance
df.gid=df2genind(wgmlst,ploidy = 1)
strata(df.gid)=data.frame(temp$food)
amo1=poppr.amova(df.gid,~temp.food,nperm=999)
cbind(amo1$componentsofcovariance,amo1$results$Df)

# Table 1 - Genetic differentiation between food pairs
set.seed(432)
fd.gene=as.data.frame(cbind(temp$food,wgmlst));names(fd.gene)[1]='food'
amo=lapply(1:(n.f-1),function(i) lapply((i+1):n.f,function(j) {
  dat=fd.gene[fd.gene$food %in% c(fd[i],fd[j]),]
  strata.df=data.frame(food=dat$food)
  df.g=df2genind(dat[,!names(dat)=='food'],ploidy = 1)
  strata(df.g)=strata.df
  g=poppr.amova(df.g,~food,nperm=999)
  p=randtest(g,nrepet = 999)
  
  list(fd_pair=paste(fd[i],fd[j],sep='-'),amova=g,pvalue=p$pvalue)
}))

phi=lapply(amo,function(p) unlist(lapply(p, function(g) round(g$amova$statphi*100,1))))
pvalue=lapply(amo,function(x) unlist(lapply(x, function(g) g$pvalue)))
mt=data.frame(matrix(rep(NA,length(fd)^2),ncol=length(fd)))
lapply(1:length(phi),function(p){
  mt[p,(length(fd)-length(phi[[p]])+1):length(fd)]<<-paste(paste(phi[[p]],pvalue[[p]],sep='('),')',sep='')})
names(mt)=fd
row.names(mt)=fd
mt


###########################################
#   Table 2 Summary of predictive metrics
###########################################

rf <- rf.core_top100 
train.df <- train.df.core_top100

cond=data.frame(rf$predicted.oob)
pd=c("_pred_rf", "_pred_rf", "_pred_rf", "_pred_rf", "_pred_rf")
names(cond)=paste(names(cond),pd,sep='')
res=data.frame(to.matrix(train.df$food,levels(train.df$food)))
names(res)=paste(names(res),'_true',sep='')
dat=cbind(res,cond)

m=multi_roc(dat) # Area under ROC curve

se=sen_spe(rf$class.oob,train.df$food)$se # sensitivity
sp=sen_spe(rf$class.oob,train.df$food)$sp # specificity

# Table 2: sensitivity, specificity, and area under ROC curve
dd=data.frame(food=names(se),se=round(se,2),sp=round(sp,2),AUC=round(as.numeric(unlist(m$AUC)),2)[1:length(se)])
dd=dd[-dim(dd)[1],]
dd

# Cohen’s kappa 
sen_spe(rf$class.oob,train.df$food)$kappa


#####################
#    Table 3
#####################
top_n <- c(10, 20, 30, 50, 70, 100, 200, 500, 1000, 1748) 
n_top_in_core <- lapply(top_n, function(i) rk.genes[1:i] %in% cgmlst_loci %>% sum) %>% do.call(rbind, .)
data.frame(top_n = top_n, n_top_in_core = n_top_in_core) %>% 
  mutate(row_pct = round(n_top_in_core / top_n * 100, 0))


##########################################################################################################
#      Figure 1 prediction accuracy based on different missing value thresholds and cutoffs for cluster
##########################################################################################################

######## selection criteria #########################
# this section examines the prediction accuracy of random forest models built from dataset created 
# under different cluster criteria (h) and different proportions of missing values (m)

h=c(0,0.001,0.004,0.01) ## row selection, h is equivalent to ht above
m<-c(0,0.01,0.05,0.1,0.3,1) #### column selection, m is equivalent to miss_thres above

hh=vector('list',length(h)*length(m))
hh[[1]]=c(0,.01)
hh[[2]]=c(0,.05)
hh[[3]]=c(0,.1)
hh[[4]]=c(0,1)
hh[[5]]=c(.001,.01)
hh[[6]]=c(.001,.05)
hh[[7]]=c(.001,.1)
hh[[8]]=c(.001,1)
hh[[9]]=c(.004,.01)
hh[[10]]=c(.004,.05)
hh[[11]]=c(.004,.1)
hh[[12]]=c(.004,1)
hh[[13]]=c(.01,.01)
hh[[14]]=c(.01,.05)
hh[[15]]=c(.01,.1)
hh[[16]]=c(.01,1)
hh[[17]]=c(0,0)
hh[[18]]=c(.001,0)
hh[[19]]=c(.004,0)
hh[[20]]=c(.01,0)
hh[[21]]=c(0,0.3)
hh[[22]]=c(.001,0.3)
hh[[23]]=c(.004,0.3)
hh[[24]]=c(.01,0.3)

library(doParallel)
library(foreach)

ncores <- detectCores() - 1 
cl=makeCluster(ncores)
registerDoParallel(cl) #number of cores on the machine

res.se=foreach(n=1:length(hh),.packages=c('tidyverse','cluster','randomForestSRC')) %dopar% {
  
  if(hh[[n]][1]==0) ssi=as.character(lm_dat$SRR_ID) else ssi=sel_rep_iso(lm_dat,hh[[n]][1]) ### select representative isolates at different levels of ht
  dat<-lm_dat %>% filter(SRR_ID %in% ssi) %>% select('food',starts_with('LMO'))
  dat$food<-factor(dat$food) 
  dat <- dat[,apply(dat,2,function(x) sum(is.na(x))/length(x))<= hh[[n]][2]]   ### excluding loci with proportion of missingness > miss_thres
  dat <- dat %>% 
    mutate_if(is.integer, coalesce, 0L) %>% 
    mutate(across(starts_with("LMO"), ~as.factor(as.character(.x)))) %>%
    as.data.frame
  set.seed(78)
  m=rfsrc(food~.,dat,ntree=1000) 
  tab=table(m$yvar,m$class.oob) # m$yvar contains values for food in dat, wheares m$class.oob contains values for food predicted by the random forest model for the purpose of caluclating OOB
  accu=c(diag(tab)/table(m$yvar),sum(diag(tab))/length(m$yvar))
  # alternative code for the two rows above: accu = 1-m$err.rate[m$ntree,]
  return(list(model=m,matrix=tab,accu=accu,rep=hh[[n]][1],dsize=length(ssi),dat=dat, miss=hh[[n]][2]))
}


df=data.table(do.call(rbind,lapply(res.se,function(x) x$matrix))) # combine the 24 m$yvar - m$class.oob to a df
df$row_sele=unlist(lapply(res.se,function(x) rep(x$rep,length(unique(lm_dat$food))))) # h
df$column_sele=unlist(lapply(res.se,function(x) rep(x$miss,length(unique(lm_dat$food))))) # m
df$n.genes=unlist(lapply(res.se,function(x) rep(dim(x$dat)[2] - 1,length(unique(lm_dat$food))))) # -1 is to exclude "food" from counting
df$n.iso=unlist(lapply(res.se,function(x) rep(dim(x$dat)[1],length(unique(lm_dat$food)))))
df$n.iso1=unlist(lapply(res.se,function(x) rep((x$dsize),length(unique(lm_dat$food))))) # exactly the same as df$n.iso
df$accuracy=100*unlist(lapply(res.se,function(x) x$accu[-length(x$accu)])) # accuracy for each food under each h and each m

df.lg=(gather(df,index,value,1:length(unique(lm_dat$food))))  # transform df to long form
df.lg$food=rep(levels(lm_dat$food),nrow(df))
df.lg$row_sele=factor(df.lg$row_sele)
df.lg$column_sele=factor(df.lg$column_sele)

# Figure

plot.col=c("purple", "orange", "red", "blue", "green")  

df %>%
  select(row_sele, column_sele, n.genes, n.iso) %>%
  distinct %>%
  arrange(column_sele, row_sele) -> df1  

lapply(unique(df1$column_sele), function(pct_missing) df1 %>% 
         filter(column_sele == pct_missing) %>% 
         `$`(n.genes) %>% 
         paste0(collapse = ", ") %>%
         paste0("(No. of genes: ", ., ")")) -> n_genes_by_pctNA  

  
strip.txt=list('0'='without missing alleles','0.01'='missingness <= 1%','0.05'='missingness <= 5%',
               '0.1'='missingness <= 10%','0.3'='missingness <= 30%','1'='all wgMLST genes included')  
strip.txt=lapply(1:length(strip.txt),function(x) paste(strip.txt[[x]],n_genes_by_pctNA[x],sep='\n'))  
labeller.fun=function(variable,value){
  return(strip.txt[value])
} 
############ modify x labels #########################
xt=rev(names(table(df.lg$n.iso))) 

x.txt=paste(paste(paste(as.numeric(as.character(unique(df.lg$row_sele)))*100,'%',sep=''),xt,sep='\n(n='),')',sep='') 

lt1=c("dashed", "dotted", "dotdash", "longdash", "twodash", "F1", "4C88C488", "12345678")   
lt=lt1[1:length(levels(lm_dat$food))]   


png(filename= "fig1.png", width = 8*400, height =8*300, units = "px", res=400,pointsize = 15, bg = "white", restoreConsole = TRUE)

ggplot(df.lg,aes(row_sele,accuracy,group=food,color=food,linetype=food))+geom_line(size=1)+geom_point(aes(shape = food), size = 2)+
  scale_linetype_manual(values=lt)+
  facet_wrap(~column_sele,labeller = labeller.fun)+
  scale_color_manual(values = plot.col)+
  scale_shape_manual(values = c(3, 6, 15, 9, 17)) + 
  scale_y_continuous(breaks = seq(0, 100, 20), 
                     labels = as.character(seq(0, 100, 20))) +
  scale_x_discrete(labels=x.txt)+
  xlab('cutoff values of proportional pairwise difference\n (number of isolates selected)')+
  ylab('prediction accuracy (%)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'top', 
        legend.title = element_blank())
dev.off()


############################################################
#    Figure 2 predictability over ranked loci 
############################################################

library(data.table)
library(doParallel)
library(foreach)

ncores <- detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl) #number of cores on the machine
train.imp <- train.df.all
my.seq=c(dim(train.imp)[2]-1,seq(1600,1000,-300),seq(900,100,-200),seq(90,50,-20),seq(45,10,-5))
my.seq=c(my.seq,c(8,5,3))
rk.genes=names(sort(c.imp,decreasing = T)) # top ranked genes
rfsrc.seq=seq.rfsrc(train.imp,rk.genes,my.seq,ntree=1000) 

dat=rfsrc.seq$res
ltype=c("dotdash", "dotted", "solid", "dashed", "twodash", "longdash")
my.col=c("red", "orange", "black",   "purple", "green", "blue")

dat$food=factor(dat$food,levels=levels(reorder(dat$food,-dat$se)))
ltype=ifelse(levels(as.factor(dat$food))=='overall','solid','longdash')


# figure 
png(filename= "fig2.png", width = 8*400, height =8*300, units = "px", res=400,pointsize = 15, bg = "white", restoreConsole = TRUE)

ggplot(dat,aes(loci,se,group=food,color=food))+geom_line(aes(linetype=food),size=1)+ylab('estimated sensitivity')+
  geom_point(aes(shape = food), size=2)+
  coord_trans(x = "log")+scale_x_continuous(breaks=c(3,5,10,20,30, 50, 70,100,200,500,1000,2000, 4804)) + 
  xlab('number of the ranked wgMLST loci')+
  scale_linetype_manual(values=ltype)+
  scale_shape_manual(values = c(15,6, 16, 3, 17, 9)) + 
  scale_color_manual(values=my.col)+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = sprintf("%4.2f", seq(0, 1, 0.2))) +
  theme(legend.position='right',
        text=element_text(size = 15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        panel.grid.minor.x = element_blank())+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1), 
         linetype=guide_legend(keywidth = 3, keyheight = 1),
         colour=guide_legend(keywidth = 3, keyheight = 1))

dev.off()


#################################################################
#    Figure 3 panel plot of conditional probability of isolates 
#################################################################

rf <- rf.core_top100
train.df <- train.df.core_top100


pred=data.frame(rf$predicted.oob)   # predicted probability for each food category
names(pred)=levels(train.df$food)
pred$pred.fd=rf$class.oob           # food predicted with the largest probability

pred$obs=row.names(train.df)        # row numbers in train.df
pred$obs.fd=train.df$food           # observed food category

pred.l=data.table(gather(pred,food,prob,1:(dim(pred)[2]-3)))  # transfer wide-form data to long-form in which each isolates have X rows, 
                                                              # the food column in the long-form stores the X food categories in the wide-form,
                                                              # the prob column in the long-form stores the predicted probabilities for the X food categories in the wide-form,
                                                              # Other columns in the wide-form (pred.fd, obs. obs.fd) were replicated to X rows in the long-form

pp=plot_panel_pred_prob_ind(pred.l,train=TRUE) ### train indicate train data and correct prediction is rugged


png(filename= "fig3.png", width = 8*400, height =8*300, units = "px", res=400,pointsize = 15, bg = "white", restoreConsole = TRUE)

pp$p

dev.off()





