
# final code for Weidong's Listeria manuscript as of September 2020
# Table 2 is produced using wgMLST_pop_gen.R
# sel_rep_iso() works well under R-3.6.1 but gives fetal errors under R-4.1.2

library(tidyverse)
library(randomForestSRC)
library(dendextend)
library(party)
library(ape) 
library(cluster)
library(data.table)

data_folder <- "../data"
results_folder <- "../results/"
plot_folder <- "../plots"

# isolates with confirmed food source
lm_dat <- readRDS(paste0(data_folder, "/isolates_original_plus_new_dec_1_2021.rds"))
 
### cgMLST loci
cgmlst_loci <- read.csv(paste0(data_folder, "/cgMLST_loci.csv")) %>% names


# rm(list = ls())

sel_rep_iso<-function(dat,ht){
  #genes<-dat %>% select(starts_with("LMO")) # select all available genes
  genes<-dat %>% select(all_of(cgmlst_loci))  ###### select core genes ######
  genes<-data.frame(sapply(genes,as.factor))
  epi<-dat %>% select(-starts_with('LMO')) # non-gene variables
  
  suppressWarnings(cl<-daisy(genes,metric='gower'))  # works well under R-3.6.1 but gives fetal errors under R-4.1.2
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


################################################################
#  Deduplicate isolates based on core gene before imputation
################################################################

ht<-0.004 ##### proportional difference within clusters defined by ht

si<-sel_rep_iso(lm_dat, ht) ### select representative isolates


###################################################################
#  importance of genes from random forest model based on all genes
###################################################################
train.df.all <-lm_dat %>% 
  filter(SRR_ID %in% si) %>% 
  select('food', starts_with("LMO")) %>% 
  mutate_if(is.integer, coalesce, 0L) %>% # integer "LMOxxxxx" as integer (52.59%) performs similar to "LMOxxxxx" as factor (51.85%)
  mutate(across(starts_with("LMO"), ~as.factor(as.character(.x)))) %>%
  as.data.frame # rfsrc() doesn't work with a tibble

set.seed(23)
rf.all <- rfsrc(food ~., train.df.all, importance = T) # Weidong used importance = T for sensitivity analysis using ranked importance, ZC edited to importance = F because only the full model will be used for source prediction

# saveRDS(rf, paste0(results_folder, "rf_importance.rds"))
# rf.all <- readRDS(paste0(results_folder, "rf_importance.rds"))
# rf_importance.rds was then copied and pasted to git

################################################################################
#   Final model - Random forest model based on core gene + top 100 ranked genes
################################################################################

# rf.all <- readRDS("data-raw/lm_wgs_manuscript/datasets/rf_importance.rds")
rf <- rf.all
c.imp=rf$importance[,1]
rk.genes=names(sort(c.imp,decreasing = T)) # [ZC: top ranked genes]

train.df.core_top100 <-lm_dat %>% 
  filter(SRR_ID %in% si) %>% 
  select('food', all_of(c(cgmlst_loci, rk.genes[1:100]))) %>% 
  mutate_if(is.integer, coalesce, 0L) %>% 
  mutate(across(starts_with("LMO"), ~as.factor(as.character(.x)))) %>%
  as.data.frame 

set.seed(23)
rf.core_top100 <- rfsrc(food ~., train.df.core_top100, importance = F) 


##########################################
#  Table 1 Summary of PulseNet isolates
##########################################
table(lm_dat$food)
table(train.df$food)




#### no missingness of wgMLST #############
train.complete=train.df[,sapply(train.df,function(x) sum(is.na(x))==0)]
set.seed(23)
rf.complete=rfsrc(food~.,train.complete,ntree=1000,importance = T)
rf
rf.complete
#train.df<-data.frame(sapply(train.df,factor))


### party package ###########
#crf<-cforest(food~.,train.df,control =cforest_control(ntree=1000,mtry=sqrt(dim(train.df)[2]),minbucket=3))
#pred=predict(crf,OOB=T)
#table(train.df$food,pred)
###########################################

#rf<-randomForest(food~.,train.df,mtry=sqrt(dim(train.df)[2]),ntree=1000)






#### sensitivity analysis of ranked loci on bootstrap samples ####################################

dt=lapply(1:length(res.se),function(x) {
  d=data.frame(loci=names(sort(res.se[[x]]$model$importance[,1],decreasing = T)),
  vimp=sort(res.se[[x]]$model$importance[,1],decreasing = T))
  d$row_sele=hh[[x]][1]
  d$column_sele=hh[[x]][2]
  d})
dt=data.table(do.call(rbind,dt))
tp=30
dt.s=dt[,.((loci[1:tp])),by=c('row_sele','column_sele')]
dt.s$V1=factor(dt.s$V1)
dt.s=data.frame(sort(table(dt.s$V1),decreasing = T))
#ggplot(dt.s[dt.s$Freq>2,],aes(Var1,Freq,color=Var1))+geom_point()+theme(legend.position='none',axis.text.x = element_text(angle = 90, vjust = 0.5))+
#  xlab("")+ylab('frequency')+ggtitle('Frequency of top ranked wgMLST in 24 different datasets')
ggplot(dt.s[dt.s$Freq > 2,], aes(Var1, Freq, color = Var1)) +
  geom_point() + theme(legend.position = 'none', axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("") + ylab('frequency') + ggtitle('Frequency of top ranked wgMLST in 24 different datasets') +
  ggsave(paste0(plot_folder, "/frequency_plot.png"))

##################################################################################################

# source('\\\\cdc.gov\\private\\M318\\vhg8\\r scripts\\wgs\\wgMLST funs.R')
# source("//cdc.gov/project/CCID_NCZVED_DFBMD_EDEB/Analytics/____ARCHIVE/Weidong/r scripts\\wgs\\wgMLST funs.R")
#c.imp=varimp(rf,OOB=T)
#p.imp=names(sort(c.imp,decreasing = T))[1:100]
#d.imp=train.df[,names(train.df) %in% c('food',p.imp)]
rf <- rf.all
train.df <- train.df.all
c.imp=rf$importance[,1]
p.imp=names(sort(c.imp,decreasing = T))[1:100]
d.imp=train.df[,names(train.df) %in% c('food',p.imp)]

#rf.imp<-cforest(food~.,d.imp,control =cforest_control(ntree=1000,mtry=sqrt(dim(d.imp)[2]),minbucket=3))


#save(train.df,rf,rf1,file='\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\rf.RData')
#load('\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\rf.RData')

# load("//cdc.gov/project/CCID_NCZVED_DFBMD_EDEB/Analytics/____ARCHIVE/Weidong/manuscripts docs/wgs/EID/revision/rf.RData")

# rf=rf70    # [object 'rf70' not found]






###########################################
#   Table 3 Summary of predictive metrics
###########################################
rf <- rf.core_top100 
train.df <- train.df.core_top100

library(multiROC)
library(data.table)
library(phytools)
#cond=data.frame(do.call(rbind,rf@cond_distr_response()))
#cond=data.frame(predict(rf,type='prob'))
cond=data.frame(rf$predicted.oob)
pd=c("_pred_rf", "_pred_rf", "_pred_rf", "_pred_rf", "_pred_rf")
names(cond)=paste(names(cond),pd,sep='')
res=data.frame(to.matrix(train.df$food,levels(train.df$food)))
#res=data.frame(rf@responses@transformations)
names(res)=paste(names(res),'_true',sep='')
dat=cbind(res,cond)

m=multi_roc(dat)

####### estimate prediction uncertainty using multinomial distributio of conditional probability##

# source('\\\\cdc.gov\\private\\M318\\vhg8\\r scripts\\wgs\\wgMLST funs.R')
# source("//cdc.gov/project/CCID_NCZVED_DFBMD_EDEB/Analytics/____ARCHIVE/Weidong/r scripts\\wgs\\wgMLST funs.R")
#temp=do.call(rbind,rf@cond_distr_response())
#temp1=data.frame(obs=1:nrow(train.df),pred.fd=levels(train.df$food)[rf@predict_response()],temp)

temp=cond
temp1=data.frame(obs=1:nrow(train.df),pred.fd=rf$class.oob,temp)

temp2=temp1 %>% gather('food','prob',-c(1:2))
#temp2$food=sapply(temp2$food,function(x) strsplit(x,split = '\\.')[[1]][2])
temp2=data.table(temp2[order(temp2$obs),])
conf=cal_conf_attr(n.rep=500,temp2,train=TRUE)
conf$attribution[,2:5]=round(conf$attribution[,2:5],2)
out=conf$attribution
out$report=paste(paste(paste(out$mean,out$lower,sep='('),out$upper,sep='-'),')',sep='')
out$food=levels(train.df$food)
out=out[order(out$food),]
#conf$conf95_se
#conf$conf95_kappa     # [ZC: kappa and 95% CI]
se=sen_spe(rf$class.oob,train.df$food)$se
sp=sen_spe(rf$class.oob,train.df$food)$sp

dd=data.frame(food=names(se),se=round(se,2),sp=round(sp,2),AUC=round(as.numeric(unlist(m$AUC)),2)[1:length(se)])
dd=dd[-dim(dd)[1],]

library(xtable)

#print(xtable(dd),'\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\Epidemiology\\se_sp_auc.html',type='html',include.rownames = F)
print(xtable(dd),'C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/se_sp_auc.html',type='html',include.rownames = F)

sen_spe(rf$class.oob,train.df$food)






###### Figure 1 fan phylo tree ############################################################
library(ggtree)
library(ggdendro)

# dat<-read.csv("\\\\cdc.gov\\project\\CCID_NCZVED_DFBMD_EDEB\\Analytics\\Weidong\\lm model\\suppl.csv")
# use suppl.csv because it contains WGS_id
# dat<-read.csv("//cdc.gov/project/CCID_NCZVED_DFBMD_EDEB/Analytics/____ARCHIVE/Weidong/Whole Genomo Sequence/data/eid data/suppl.csv")
dat <- lm_dat
train.df <- train.df.all

sel.fd=levels(train.df$food)
# plot.col=c("blue", "brown", "green", "pink", "red") # original colors
plot.col=c("purple", "orange", "red", "blue", "green")
plot.id='food'
plot.dat=train.df.all %>% bind_cols(lm_dat %>% filter(SRR_ID %in% si) %>% select(-names(train.df.all)))

p=vis.food.iso.clust(plot.dat,sel.fd,ylim=c(0,1),draw_ori = NULL,plot.id,plot.col)

# [wgs/EID/fig1_fan_tree.png is figure 1 on the manuscript]
# none of the three figures below is the finalized
# win.metafile("C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/fig1_fan_plot.wmf", 
#              width = 11.5, height = 7)
# pdf("C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/fig1_fan_plot_cgMLST1000.pdf", width = 8, height = 6)
png(filename= "C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/fig1_fan_plot_dot2.png",
    width = 8*800, height =8*800, units = "px", res=400,pointsize = 15, bg = "white", restoreConsole = TRUE)

p$fan.tree_dot

dev.off()


png(filename= "C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/fig1_fan_plot_label2.png",
    width = 8*800, height =8*800, units = "px", res=400,pointsize = 15, bg = "white", restoreConsole = TRUE)

p$fan.tree_label

dev.off()



##########################################################################################################
#      Figure 2 prediction accuracy based on different missing value thresholds and cutoffs for cluster
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
  dat$food<-factor(dat$food) # identical to train.df when ht<-0.004 and miss_thres<-0.01
  dat <- dat[,apply(dat,2,function(x) sum(is.na(x))/length(x))<= hh[[n]][2]]   ### excluding loci with proportion of missingness > miss_thres
  dat <- dat %>% 
    mutate_if(is.integer, coalesce, 0L) %>% 
    mutate(across(starts_with("LMO"), ~as.factor(as.character(.x)))) %>%
    as.data.frame
  set.seed(78)
  m=rfsrc(food~.,dat,ntree=1000) # impute missing values and build random forest model
  tab=table(m$yvar,m$class.oob) # m$yvar contains values for food in dat, wheares m$class.oob contains values for food predicted by the random forest model for the purpose of caluclating OOB
  accu=c(diag(tab)/table(m$yvar),sum(diag(tab))/length(m$yvar))
  # alternative code for the two rows above: accu = 1-m$err.rate[m$ntree,]
  return(list(model=m,matrix=tab,accu=accu,rep=hh[[n]][1],dsize=length(ssi),dat=dat, miss=hh[[n]][2]))
  #df<-grp_alle_impu(dat,hh[[n]][2])
  #names(df)[1]<-'food'
  #df<-data.frame(sapply(df,factor))
  #set.seed(434)
  #rf<-randomForest(food~.,df,mtry=sqrt(dim(df)[2]),ntree=1000)
  #list(matrix=rf$confusion,rep=hh[[n]][1],dsize=length(ssi),df=df,miss=hh[[n]][2])
}

# save(lm_dat,train.df,train.imp,rf,res.se,file='\\\\cdc.gov\\project\\CCID_NCZVED_DFBMD_EDEB\\Analytics\\Weidong\\lm model\\sen_1.RData') 
# load('//cdc.gov/project/CCID_NCZVED_DFBMD_EDEB/Analytics/____ARCHIVE/Weidong/Whole Genomo Sequence/Lm WGS/data/sen_1.RData') 

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

# Figure 2 in the manuscript dated 8_26_2020
ggplot(df.lg,aes(food,accuracy,group=as.factor(row_sele),color=as.factor(row_sele)))+
  geom_line()+facet_wrap(~column_sele)

plot.col=c("purple", "orange", "red", "blue", "green")
###### modify panel strip text
rn <- df %>% 
  select(row_sele, column_sele, n.genes) %>% 
  distinct %>%
  group_by(column_sele, n.genes) %>% 
  summarise(n = n()) %>%
  filter(n >= 2)
# rn=df[,.(list(range(n.genes))),by=column_sele][order(column_sele)][,V1]
rn=sapply(rn$n.genes,function(x){
  if (rn %>% filter(n.genes == x) %>% `$`(n) %>% `==`(4)) {
    paste0("(", x, ")*")
  } else {
    paste0("(", x, ")**")
  }
})

df %>%
  select(row_sele, column_sele, n.genes, n.iso) %>%
  distinct %>%
  arrange(column_sele, row_sele) -> df1

lapply(unique(df1$column_sele), function(pct_missing) df1 %>% 
         filter(column_sele == pct_missing) %>% 
         `$`(n.genes) %>% 
         paste0(collapse = ", ") %>%
         paste0("(No. of genes: ", ., ")")) -> n_genes_by_pctNA

  
# rn[length(rn)]=as.character(dim(lm_dat %>% select(starts_with('LM')))[2])
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


# win.metafile("\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\epidemiology\\sensitivity1.wmf",width = 11.5, 
#              height = 7)
#pdf("\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\epidemiology\\sensitivity1.pdf", width = 11.5, height = 7)
png(filename= "C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/fig2_sensitivity2.png", width = 8*400, height =8*300, units = "px", res=400,pointsize = 15, bg = "white", restoreConsole = TRUE)

ggplot(df.lg,aes(row_sele,accuracy,group=food,color=food,linetype=food))+geom_line(size=1)+geom_point(aes(shape = food), size = 2)+
  # geom_vline(data=subset(df.lg,column_sele==1),aes(xintercept=3),linetype=2,color='red')+
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

# overall prediction accurary by col_sele and row_sele
# lapply(1:length(hh), function(i) data.frame(row_sele = hh[[i]][1], 
#                                             col_sele = hh[[i]][2], 
#                                             over_accu = res.se[[i]]$accu[[6]])) %>% 
#   bind_rows() %>% arrange(col_sele, row_sele) %>% 
#   write.xlsx("C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/overall_accu.xlsx")





####### Figure 3 predictability over ranked loci #########################################
library(data.table)
library(doParallel)
library(foreach)

ncores <- detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl) #number of cores on the machine
train.imp <- train.df.all
my.seq=c(dim(train.imp)[2]-1,seq(1600,1000,-300),seq(900,100,-200),seq(90,50,-20),seq(45,10,-5))
my.seq=c(my.seq,c(8,5,3))
##cf.seq=seq.cforest(train.imp,names(sort(c.imp,decreasing = T)),my.seq,ntree=1000) ## party package
rk.genes=names(sort(c.imp,decreasing = T)) # [ZC: top ranked genes]
rfsrc.seq=seq.rfsrc(train.imp,rk.genes,my.seq,ntree=1000) ## randomforestsrc package
# save(train.df,rf,rk.genes, rfsrc.seq, file='\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\rf.RData')
# load('\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\rf.RData')

dat=rfsrc.seq$res
ltype=c("dotdash", "dotted", "solid", "dashed", "twodash", "longdash")
my.col=c("red", "orange", "black",   "purple", "green", "blue")

dat$food=factor(dat$food,levels=levels(reorder(dat$food,-dat$se)))
ltype=ifelse(levels(as.factor(dat$food))=='overall','solid','longdash')

# gsave("\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\fig2_pred_accu_loci.eps", width = 10.5, height = 7)


# [figure 3 on the manuscript]
png(filename= "C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/fig3_pred_accu_loci2.png", width = 8*400, height =8*300, units = "px", res=400,pointsize = 15, bg = "white", restoreConsole = TRUE)

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



############# Figure 4 panel plot of conditional probability of isolates ########################
# [the colors of meat and vegetables are inconsistent with those on the manuscript]

#pred=data.frame(do.call(rbind,rf@cond_distr_response()))
rf <- rf.core_top100
train.df <- train.df.core_top100


pred=data.frame(rf$predicted.oob)   # [ZC: predicted probability for each food category]
names(pred)=levels(train.df$food)
pred$pred.fd=rf$class.oob           # [ZC: food predicted with the largest probability]

pred$obs=row.names(train.df)        # [ZC: row numbers in train.df]
pred$obs.fd=train.df$food           # [ZC: observed food category]

pred.l=data.table(gather(pred,food,prob,1:(dim(pred)[2]-3)))  # [ZC: transfer wide-form data to long-form in which each isolates have X rows, 
                                                              #     [the food column in the long-form stores the X food categories in the wide-form,
                                                              #     [the prob column in the long-form stores the predicted probabilities for the X food categories in the wide-form,
                                                              #     [Other columns in the wide-form (pred.fd, obs. obs.fd) were replicated to X rows in the long-form]

pp=plot_panel_pred_prob_ind(pred.l,train=TRUE) ### train indicate train data and correct prediction is rugged

# win.metafile("\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\fig3_pred_prob_ind.wmf", width = 10.5, height = 7)
# pdf("\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\fig3_pred_prob_ind.pdf", width = 10.5, height = 7)

png(filename= "C:/Users/nyv5/OneDrive - CDC/Weidong/Lm WGS/results/fig4_pred_prob2.png", 
    width = 8*400, height =8*300, units = "px", res=400,pointsize = 15, bg = "white", restoreConsole = TRUE)

pp$p

dev.off()






##### tanglegram supplement figure ################################

library(phytools)
library(dendextend)
library(ape)
library(tidyverse)
library(DECIPHER)

eid_iso=read.delim('\\\\cdc.gov\\project\\CCID_NCZVED_DFBMD_EDEB\\Analytics\\Weidong\\LM model\\ncbi_access_dat.csv',sep=',',header=T,stringsAsFactors = F)

my.col=c('blue','brown','green','pink','red')

s801=read.tree('\\\\cdc.gov\\project\\CCID_NCZVED_DFBMD_EDEB\\Analytics\\Weidong\\LM model\\weeded_variants_fasttree_801.dnd')
m801=read.tree('\\\\cdc.gov\\project\\CCID_NCZVED_DFBMD_EDEB\\Analytics\\Weidong\\LM model\\DendroExport_801.dnd')

r_s801=(midpoint.root(s801))
r_m801=(midpoint.root(m801))

r_s801b <- compute.brlen(r_s801)
r_m801b <- compute.brlen(r_m801)

dend_s801=chronos(r_s801b)
dend_m801=chronos(r_m801b)

dend_s801c <- multi2di(dend_s801)

conn_l_col2=eid_iso$food[match(dend_s801c$tip.label,eid_iso$WGS_id )]
no_sel=match(dend_s801c$tip.label,eid_iso$WGS_id[eid_iso$selected=='no'] )

col.s2=my.col[as.factor(conn_l_col2)]
col.s3=ifelse(!is.na(no_sel),alpha('grey',0.85),col.s2)
col.s4=ifelse(col.s3=='grey', alpha(muted(col.s2, l = 75), 0.85),col.s3)

dendl801 <- dendextend::untangle(as.dendrogram(dend_s801c),
                                 as.dendrogram(dend_m801),
                                 method = "step2side")

save(dend_s801,dend_m801,dendl801, file="\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\dendl801.RData")

load("\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\dendl801.RData")

dendl801 %>% set("branches_lwd", 1) %>%
  set("labels_col", "white") -> dendl801

#win.metafile("~/Downloads/801_tanglegram.wmf", width = 10.5, height = 7)
#pdf("~/Downloads/801_tanglegram.pdf", width = 10.5, height = 7)
win.metafile("\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\suppl_tanglegram.wmf", width = 10.5, height = 7)
pdf("\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\suppl_tanglegram.pdf", width = 10.5, height = 7)

tanglegram(dendl801,
           main_left='hqSNP tree',
           main_right='wgMLST tree',
           lab.cex=0.1,margin_inner=0.3, lwd=1,
           color_lines=col.s4,
           highlight_distinct_edges = FALSE)
dev.off()
