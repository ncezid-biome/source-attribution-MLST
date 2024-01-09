fan_plot_source_human=function(dat){
  library(phytools)
  #dat=dat[sample(1:nrow(dat),200,replace=F),]
  dat1=dat %>% select(starts_with('LMO'))
  row.names(dat1)=1:dim(dat1)[1]
  cl1=daisy(dat1,metric='gower')
  hc1=hclust(cl1,method = 'complete')
  mytree1=(as.phylo(hc1))
  mytree1$tip.label=1:length(mytree1$tip.label)
  #x=dat$food
  #mytree1$edge.length[mytree1$edge.length==0]<-max(nodeHeights(mytree1))*1e-6
  #m<-ace(x,mytree1,model="ER",type="discrete",marginal = F)
  #m1<-ace(x,mytree1,model="ER",type="discrete",marginal = T)
  
  my.col=c('blue','brown','red','pink','green') ## for dairy, fruit, meat, seafood, veg
  
  d1.tipcol=setNames(sapply(dat$food,function(x) my.col[which(levels(dat1$food)==x)]),dat1$food)
  col=setNames(my.col,levels(dat1$food))
  
  png(filename='\\\\cdc.gov\\private\\M318\\vhg8\\Whole Genomo Sequence\\Lm WGS\\prediction of human isolates\\fan_plot_so_human.png', 
      width = 8*450, height =8*300, units = "px", res=200,pointsize = 15, bg = "white", restoreConsole = TRUE)
  
  
  plot(mytree1,tip.color=d1.tipcol,show.tip.label=F,type='fan')
  nodelabels(node=1:nrow(dat1),cex=0.15,pie=to.matrix(dat$food,levels(dat$food)),piecol = col)
  so=which(dat$method=='source')
  tiplabels(tip=so,pie=to.matrix(dat$food[so],levels(dat$food)),
            piecol = col,cex=.15,offset = 0.025)
  
  legend("topright",levels(dat$food),fill=my.col,cex=1.5)
  dev.off()
}

fit_ace_model=function(dat){
  library(geiger)

  cl=daisy(dat,metric='gower')
  hc=hclust(cl,method = 'complete')
  mytree=(as.phylo(hc))
  x=dat$food
  names(x)=rownames(dat)
  mytree$edge.length[mytree$edge.length==0]<-max(nodeHeights(mytree))*1e-6
  m=fitDiscrete(mytree,dat=x,model='ARD')
  m1=fitDiscrete(mytree,dat=x,model='ER')
}
plot_selected_node=function(dat,hum.id,lv.node,xlim=NULL){
  dat1=rbind(dat,imp.human[hum.id,])
  row.names(dat1)=1:dim(dat1)[1]
  cl1=daisy(dat1,metric='gower')
  hc1=hclust(cl1,method = 'complete')
  mytree1=(as.phylo(hc1))
  mytree1$tip.label=1:length(mytree1$tip.label)
  x=dat1$food
  mytree1$edge.length[mytree1$edge.length==0]<-max(nodeHeights(mytree1))*1e-6
  m<-ace(x,mytree1,model="ER",type="discrete",marginal = F)
  m1<-ace(x,mytree1,model="ER",type="discrete",marginal = T)
  
  my.col=c('blue','brown','red','pink','green') ## for dairy, fruit, meat, seafood, veg
  
  d1.tipcol=setNames(sapply(dat1$food,function(x) my.col[which(levels(dat1$food)==x)]),dat1$food)
  col=setNames(my.col,levels(dat1$food))
  plot(mytree1,tip.color=d1.tipcol,label.offset=0.007,show.tip.label=F)
  tiplabels(pie=to.matrix(dat1$food,levels(dat1$food)),piecol = col,cex=.3)
  nodelabels(cex=0.5)
  nodelabels(pie=m$lik.anc,piecol=col,cex=0.5)
  
  np1=nodepath(mytree1)
  (nd1=np1[length(np1)][[1]])
  s.node=nd1[length(nd1)-(1+lv.node)]
  #s.tr=extract.clade(mytree1,s.node)
  #### my.tips needs some work####### 
  
  my.tips=prop.part(mytree1)[[s.node - length(mytree1$tip.label)]]
  s.tr=drop.tip(mytree1,mytree1$tip.label[!mytree1$tip.label %in% my.tips],subtree = T)
  
  lst=sapply(np1,function(x) if(length(which(x==s.node)>0)) x[(length(x)-which(x==s.node)-1):(length(x)-1)])
  inv.nodes=sort(unique(unlist(lst)))
  corr.ori.nodes=inv.nodes-Ntip(mytree1)
  s.like=m$lik.anc[corr.ori.nodes,]
  
  #s.col=ifelse(s.tr$tip.label==nrow(dat1),'red','black')
  
  t.col=sapply(dat1$food[my.tips],function(f) my.col[levels(dat1$food)==f])
  pm=to.matrix(dat1$food[my.tips],levels(dat1$food))
  my.t.lab=ifelse(s.tr$tip.label %in% mytree1$tip.label,ifelse(s.tr$tip.label==dim(dat1)[1],'human',
          as.character(x[as.numeric(s.tr$tip.label)])),s.tr$tip.label)
  
  png(filename= paste('\\\\cdc.gov\\private\\M318\\vhg8\\Whole Genomo Sequence\\Lm WGS\\prediction of human isolates\\sel_hm_ace',paste(hum.id,paste(lv.node,'.png',sep=''),sep='_')), 
      width = 8*430, height =8*300, units = "px", res=350,pointsize = 15, bg = "white", restoreConsole = TRUE)
  
  plot(s.tr,tip.col=t.col,show.tip.label=F,cex=.6,x.lim=xlim)
  tiplabels(text=my.t.lab,col=t.col,frame='none',bg=t.col,cex=0.7)
  nodelabels(pie=s.like,piecol = col,cex=.5)
  add.arrow(tip=which(s.tr$tip.label==nrow(dat1)),col='orange')
  
  my.usr=par()$usr
  add.simmap.legend(colors=col,prompt=FALSE,x=1*par()$usr[1],
                    y=11.1*my.usr[3],fsize=0.5)
  
  dev.off()
}

ace_ard=function(dat,hm.dat){
  ran=lapply(1:nrow(hm.dat),function(n) {
    print(n)
    #temp=sapply(1:50,function(n) {
    dat1=rbind(dat,hm.dat[n,])
    cl1=daisy(dat1,metric='gower')
    hc1=hclust(cl1,method = 'complete')
    mytree1=(as.phylo(hc1))
    np1=nodepath(mytree1)
    (nd1=np1[length(np1)][[1]])
    (nd1=nd1[length(nd1)-1]-Ntip(mytree1))
    x=dat1$food
    mytree1$edge.length[mytree1$edge.length==0]<-max(nodeHeights(mytree1))*1e-6
    
    lapply(levels(as.factor(dat$food)),function(f){
      x[length(x)]=f
      m<-ace(x,mytree1,model="ARD",type="discrete",marginal = F)
      prob=as.numeric(m$lik.anc[nd1,])
      prob[prob<0]=0
      list(prob=prob,assign=f,pred=levels(dat$food)[which.max(prob)],
           simp=diversity(prob,index='simpson'))
    })
  })
}

ace_discrete=function(dat,hm.dat){
  #lab.col=data.frame(sapply(lab.col,function(x) as.factor(x)))
  dat=data.frame(sapply(dat,as.factor))
  cl=daisy(dat,metric='gower')
  hc=hclust(cl,method = 'complete')
  mytree=(as.phylo(hc))
  x=factor(dat$food)
  mytree$edge.length[mytree$edge.length==0]<-max(nodeHeights(mytree))*1e-6
  fitER<-ace(x,mytree,model="ER",type="discrete",marginal = F) ## ER equal rate
  f.ARD=ace(x,mytree,model='ARD',type='discrete',marginal=F) ### ARD all rates different
  post.pred=fitER$lik.anc
  np=nodepath(mytree)
  
  temp=sapply(1:nrow(hm.dat),function(n) {
    print(n)
    #temp=sapply(1:50,function(n) {
    dat1=rbind(dat,hm.dat[n,])
    cl1=daisy(dat1,metric='gower')
    hc1=hclust(cl1,method = 'complete')
    mytree1=(as.phylo(hc1))
    np1=nodepath(mytree1)
    (nd1=np1[length(np1)][[1]])
    (nd1=nd1[length(nd1)-1])
    id=which(sapply(np1[-length(np1)],function(n) any(n==nd1)))
    id.set=np[id]
    min.n=min(sapply(id.set,length))
    
    if (length(id)==1){
      ori.nd=id.set[[1]][length(id.set[[1]])-1]-Ntip(mytree)
    } else {
      m=do.call(cbind,id.set)
      lg=apply(m,1,function(x) unique(x))
      ori.nd=lg[which.min(sapply(lg, function(x) length(x)==1))-1][[1]]-Ntip(mytree)
    }
    post.pred[ori.nd,]
    #list(like=post.pred[ori.nd,],pred.fd=names(which.max(post.pred[ori.nd,])),simpson=diversity(post.pred[ori.nd,],index='simpson'))
  })
    
  ran=lapply(1:nrow(hm.dat),function(n) {
    print(n)
    #temp=sapply(1:50,function(n) {
    dat1=rbind(dat,hm.dat[n,])
    cl1=daisy(dat1,metric='gower')
    hc1=hclust(cl1,method = 'complete')
    mytree1=(as.phylo(hc1))
    np1=nodepath(mytree1)
    (nd1=np1[length(np1)][[1]])
    (nd1=nd1[length(nd1)-1]-Ntip(mytree1))
    x=dat1$food
    mytree1$edge.length[mytree1$edge.length==0]<-max(nodeHeights(mytree1))*1e-6
    
    lapply(levels(as.factor(dat$food)),function(f){
      x[length(x)]=f
      m<-ace(x,mytree1,model="ER",type="discrete",marginal = F)
      list(prob=m$lik.anc[nd1,],assign=f,pred=names(which.max(m$lik.anc[nd1,])),
           simp=diversity(m$lik.anc[nd1,],index='simpson'))
    })
  })
    #plot(mytree)
    #nodelabels(node=1:mytree$Nnode+Ntip(mytree))
    #plot(mytree1,'fan')
    #nodelabels(node=1:mytree1$Nnode+Ntip(mytree1))
    #list(iso=n,pred.fd=names(which.max( post.pred[g,])),prob=post.pred[g,],simpson=diversity(post.pred[g,],index='simpson'))

  #plot(mytree,show.tip.label=FALSE)
  
  #cols<-setNames(my.col,levels(dat$food))
  #plot(mytree1, 'fan', no.margin=TRUE,show.tip.label=FALSE, edge.width=1)
  #nodelabels(node=1:mytree$Nnode+Ntip(mytree),pie=fitER$lik.anc,piecol=cols,cex=0.8)
  #tiplabels(col=cols,cex=0.8)
  #tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.8)
  #add.simmap.legend(colors=cols,prompt=F,x=-.7*par()$usr[1],y=0.5*par()$usr[3])
  #return(list(assign=dat$food[nrow(dat)],pred.fd=names(which.max(post.pred[nd,])),post.pred=post.pred[nd,]))
  return(list(nnb=temp,ran_assg=ran))
  }

fan_panel=function(dat,my.col){
  ## this function plots hc using labels and colors defined in lab.col db
  #lab.col=data.frame(sapply(lab.col,function(x) as.factor(x)))
  dat= data.frame(sapply(dat,as.factor))
  dat=dat %>% select(c('food','method'),starts_with('LMO'))
  cl=daisy(dat,metric='gower')
  hc=hclust(cl,method = 'complete')
  
  ape.tree=(as.phylo(hc))
  ape.tree$tip.label=paste(dat$food,ape.tree$tip.label,sep = "_")
  grpinfo=split(ape.tree$tip.label,gsub('_\\w+','',ape.tree$tip.label))
  ape.tree=groupOTU(ape.tree,grpinfo)
  
  p.fan=ggtree(ape.tree,aes(color=group),layout='circular')
  #Modify (tip) labels 
  #This could be easily done via the %<+% operator to attach the modified version of the labels and than use geom_tiplab to display the modified version.
  d=data.frame(id=ape.tree$tip.label,label2=dat$food,method=ifelse(dat$method=='human',1,0))
  p.fan=p.fan %<+% d+geom_tiplab(size=1,aes(label=label2,angle=angle))+
    scale_color_manual(values=my.col)+theme(legend.position = 'right')
  
  p0=ggtree(ape.tree,aes(color=group))+geom_tippoint(size=1)+scale_color_manual(values=my.col)+theme(legend.position = 'right')
  
  p1=facet_plot(p0,panel='',d,geom=geom_point,aes(x=method),size=0.8)

  gt = ggplot_gtable(ggplot_build(p1))
  #gtable_show_layout(gt) # will show you the layout - very handy function
  gt # see plot layout in table format
  gt$layout$l[grep('panel-2-1', gt$layout$name)] # you want to find the column specific to panel-2
  gt$widths[7] = 0.02*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
  return(list(fan=p.fan,panel=gt))
}

tanglegram.tree=function(dat){  ##### forced comparative study of SNP vs wgMLST trees
  library(phytools)
  library(dendextend)
  library(ape)
  library(tidyverse)
  library(DECIPHER)
  
  eid_iso=read.delim('\\\\cdc.gov\\project\\CCID_NCZVED_DFBMD_EDEB\\Analytics\\Weidong\\LM model\\ncbi_access_dat.csv',sep=',',header=T,stringsAsFactors = F)
  s325=read.tree("\\\\cdc.gov\\project\\CCID_NCZVED_DFBMD_EDEB\\Analytics\\Weidong\\LM model\\weeded_variants_variants_fasttree_325.dnd")
  m325=read.tree("\\\\cdc.gov\\project\\CCID_NCZVED_DFBMD_EDEB\\Analytics\\Weidong\\LM model\\DendroExport_325.dnd")
  
  dend_s325=chronos(s325)
  dend_m325=chronos(m325)
  
  r_s325=(midpoint.root(dend_s325))
  r_m325=(midpoint.root(dend_m325))
  
  my.col=c('blue','brown','green','pink','red')
  
  conn_l_col=eid_iso$food[match(r_s325$tip.label,eid_iso$WGS_id )]
  col.s=my.col[as.factor(conn_l_col)]
  
  dendl <- dendextend::untangle(as.dendrogram(r_s325),
                                as.dendrogram(r_m325),
                                method = "step2side")
  
  dendl %>% set("branches_lwd", 1) %>%
    set("labels_col", "white") -> dendl
  
  win.metafile("~/Downloads/325_tanglegram.wmf", width = 10.5, height = 7)
  #pdf("~/Downloads/325_tanglegram.pdf", width = 10.5, height = 7)
  tanglegram(dendl,
             main_left='hqSNP tree',
             main_right='wgMLST tree',
             lab.cex=0.1,
             color_lines=col.s,
             highlight_distinct_edges = FALSE)
  dev.off()
  entanglement(dendl)
  cor.dendlist(dendl, method = "cophenetic")
  cor.dendlist(dendl, method = "baker")
  
  
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
  col.s3=ifelse(!is.na(no_sel),alpha('grey',0.25),col.s2)
  col.s4=ifelse(col.s3=='grey', alpha(muted(col.s2, l = 75), 0.25),col.s3)
  
  dendl801 <- dendextend::untangle(as.dendrogram(dend_s801c),
                                   as.dendrogram(dend_m801),
                                   method = "step2side")
  
  save(dendl801, file="\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\dendl801.RData")

  
  dendl801 %>% set("branches_lwd", 1) %>%
    set("labels_col", "white") -> dendl801
  
  #win.metafile("~/Downloads/801_tanglegram.wmf", width = 10.5, height = 7)
  #pdf("~/Downloads/801_tanglegram.pdf", width = 10.5, height = 7)
  png(filename= '\\\\cdc.gov\\private\\M318\\vhg8\\manuscripts docs\\wgs\\EID\\revision\\tangle.png', 
      width = 8*300, height =8*400, units = "px", res=400,pointsize = 15, bg = "white", restoreConsole = TRUE)
  
  tanglegram(dendl801,
             main_left='hqSNP tree',
             main_right='wgMLST tree',
             lab.cex=0.1,margin_inner=0.3, lwd=1,
             color_lines=col.s3,
             highlight_distinct_edges = FALSE)
  dev.off()
}

grp.lasso=function(dat){
  library(msgl)
  library(mlogit)
  library(caret)
  load('\\\\cdc.gov\\private\\M318\\vhg8\\r scripts\\wgs\\working.RData')
  dat=train.df[,sapply(train.df,function(x) length(levels(x))>1)]
  d=model.matrix(food~.,dat)[,-1] ### no intercept
  g=sapply(dat[,-1],function(x) length(unique(x))-1)
  groups=rep(1:length(g),g)
  gw=rep(0.1,length(g))
  #m.msgl=msgl(x=d,classes=dat$food,alpha=0.5,grouping=groups,lambda=0.2,standardize = F)
  
  cl <- makeCluster(5)
  registerDoParallel(cl)
  
  # Do 4-fold cross validation on a lambda sequence of length 100.
  # The sequence is decreasing from the data derived lambda.max to 0.2*lambda.max
  fit.cv <- msgl::cv(x=d,classes=dat$food, grouping = groups,fold = 5, alpha=0.5,lambda = 0.01, use_parallel = TRUE,standardize = F)
  
  stopCluster(cl)
  
  # Print information about models
  # and cross validation errors (estimated expected generalization error)
  fit.cv
  b.m=msgl(x=d,classes=dat$food,alpha=0.5,grouping=groups,lambda=0.029,standardize = F)
  
  
  save(fit.cv,file='\\\\cdc.gov\\private\\M318\\vhg8\\r scripts\\wgs\\grp.lasso.RData')
  load('\\\\cdc.gov\\private\\M318\\vhg8\\r scripts\\wgs\\grp.lasso.RData')
  
  dummy.m=dummyVars(food~.,dat)
  trn=predict(dummy.m, dat)
  trn.df=data.frame(trn)
  trn.df=d
  ids=sample(1:nrow(dat),150,replace=F)
  trn.model=train(y=dat$food,x=trn.df,method='multinom',trControl = trainControl(method='cv',number=3,allowParallel = TRUE))
}

comp_two_predions=function(rf.pred.prob,xgb.pred.prob){
 tab=table(rf.pred.prob$pred.fd,xgb.pred.prob$pred.fd)
 dat=data.frame(rbind(rf.pred.prob,xgb.pred.prob))
 dat=dat[,names(dat)!='pred.fd']
 dat$method=rep(c('rf','xgb'),each=nrow(rf.pred.prob))
 (dat %>% gather(food,1:(dim(dat)[2]-1)))
}

##### XGBOOST model #####################################
xgb_fun=function(dat,imp.human){
library(xgboost)
library(Matrix)
  
del.v=sapply(dat,function(x)length(levels(x)))<2
dat=dat[,!names(dat) %in% names(del.v)[del.v]]
xgb_tr=xgb.DMatrix(data=sparse.model.matrix(food~.,dat)[,-1],label=as.numeric(dat$food)-1)

bst.cv=xgb.cv(data=xgb_tr,max_depth=3,nthread=20,nrounds=10, nfold=5,params=list(objective='multi:softprob','num_class'=length(unique(dat$food))),prediction=T)

bst=xgboost(data=xgb_tr,max_depth=3,nthread=20,nrounds=10,params=list(objective='multi:softprob',
    'num_class'=length(unique(dat$food))),prediction=T)
xgb_tst=xgb.DMatrix(data=sparse.model.matrix(food~.,imp.human[,!names(imp.human) %in% names(del.v)[del.v]])[,-1])
pred_h=predict(bst,newdata=xgb_tst,reshape = T)
return(list(cv=bst.cv,pred_human=pred_h))
}
###########################################################################
plot_perc_inno_loci=function(){
  inno.id=data.frame(alle=unlist(inno_alle[-1]))
  ggplot(inno.id,aes(alle))+geom_histogram()
  
  novel_alle_by_loci=data.frame(loci=names(h.dat)[-1],n.novel=sapply(inno_alle[-1],length))
  novel_alle_by_loci$loci=factor(novel_alle_by_loci$loci,levels=levels(reorder(novel_alle_by_loci$loci,novel_alle_by_loci$n.novel)))
  
  quntl=quantile(novel_alle_by_loci$n.novel,probs = c(.20,.5,.75,.9))
  novel_alle_by_loci$col=sapply(novel_alle_by_loci$n.novel,function(v) plot.col[min(which(v<quntl))])
  novel_alle_by_loci$col[is.na(novel_alle_by_loci$col)]=plot.col[1+length(quntl)]
  novel_alle_by_loci$col=factor(novel_alle_by_loci$col,levels=plot.col)
  legd=c('<20th','20th-50th','50th-75th','75th-90th','>90th')
  ggplot(novel_alle_by_loci,aes(loci,n.novel,colour=col))+geom_point()+ylab('Number of novel alleles')+
    scale_color_manual(name='percentile',values = plot.col,labels=legd)+ggtitle('percentile of loci with novel alleles')+
    theme(axis.text.x = element_blank(),axis.line.y = element_line(color="black"),
          axis.line.x = element_line(color="black",size=0.5))+xlab('ranked loci by increaseing numbers of novel alleles')
  
  #ggscatterstats(novel_alle_by_loci,y=n.novel,x=1:(ncol(f.dat)-1),marginal = T,margins='y')
  ### fill innovative alleles as NA
  #row.names(fil_hum_imp)=row.names(h.dat)
}

plot_perc_inno_iso=function(){
  plot.col=c('red','blue','green','salmon','purple')
  perc_inno_h_iso=data.frame(iso=row.names(h.dat),perc_inno=apply(fil_hum_imp,1,function(r) sum(is.na(r))/length(r)))
  perc_inno_h_iso$iso=factor(perc_inno_h_iso$iso,levels=levels(reorder(perc_inno_h_iso$iso,perc_inno_h_iso$perc_inno)))
  
  quntl=quantile(perc_inno_h_iso$perc_inno,probs = c(.20,.5,.75,.9))
  perc_inno_h_iso$col=sapply(perc_inno_h_iso$perc_inno,function(v) plot.col[min(which(v<quntl))])
  
  perc_inno_h_iso$col[is.na(perc_inno_h_iso$col)]=plot.col[1+length(quntl)]
  perc_inno_h_iso$col=factor(perc_inno_h_iso$col,levels=plot.col)
  legd=c('<20th','20th-50th','50th-75th','75th-90th','>90th')
  p=ggplot(perc_inno_h_iso,aes(iso,perc_inno,colour=col))+geom_point()+theme(axis.ticks.x = element_blank())+
    scale_color_manual(name='percentile',values = plot.col,labels=legd)+
    ylab('proportional missing allele per isolate')+ggtitle('percentile of isolates with novel alleles')+
    theme(axis.text.x = element_blank(),axis.line.y = element_line(color="black"),
          axis.line.x = element_line(color="black",size=0.5))+xlab('ranked human isolates by increaseing numbers of novel alleles')
  
  ###### impute innovative alleles (as NA) by randomly assign food id as one of food categories
  
  #fil_hum_imp$food=sample(levels(f.dat$food),size=nrow(fil_hum_imp),replace=T,prob=table(f.dat$food)/sum(table(f.dat$food)))
  return(list(perc=perc_inno_h_iso,p=p))
}

prep_dat=function(){ ### re-categorize some foods for analysis
  ############### 
  ################################################################################################
  rem_cls=c('environ','complex food','NA')
  ana.dat=merged[!merged$food %in% rem_cls &!is.na(merged$food),] ### include human isolates for predicability
  ana.dat$food=factor(ana.dat$food)
  ana.dat=ana.dat[ana.dat$SourceState!='Western Cape',]
  ana.dat$SourceState=factor(ana.dat$SourceState)
  
  ############# recategorize foods #####################
  meat=c('beef','chicken','pork','turkey')
  #poultry=c('chicken','turkey')
  #produce=c('vegetable','fruit')
  ana.dat$food=ifelse(ana.dat$food %in% meat,'meat',
                      ifelse(ana.dat$food %in% 'cheese', 'dairy',as.character(ana.dat$food)))
  #ana.dat$food=ifelse(ana.dat$food %in% produce,'produce',as.character(ana.dat$food))
  ana.dat$food=ifelse(ana.dat$food =='ice cream','dairy',as.character(ana.dat$food))
  
  #sum(ana.dat$wgs_Id %in% eid_iso$WGS_id)
  ##### summarizing missing calls
  miss.plot(ana.dat)
  
  #### output for verify food assignment #########################################
  #ana.dat %>% select(- starts_with('LMO')) %>% write.table('\\\\cdc.gov\\private\\M318\\vhg8\\food_assig.csv',sep=',', row.names=F)
  ############################################################################################
  
  ################### illustration of clusters #########################################
  library(ggtree)
  sel.fd=c('fruit','vegetable')
  plot.col=c('red','green')
  plot.id='Key' ## lables to be plotted
  plot.id='Key'
  plot.id='food'
  plot.dat=ana.dat
  v=vis.food.iso.clust(plot.dat,sel.fd,ylim=c(0,1),draw_ori = NULL,plot.id,plot.col)
  ######################################################################################
  #####################################################################################################
  
  ####### estimate pairwise differences of isolates within PulseNet outbreak clusters ##########################################
  #p=est.pair.dist.clust(merged)
  png(filename= '\\\\cdc.gov\\private\\M318\\vhg8\\whole genomo sequence\\data\\new data\\pairwise difference outbreak cluster.png', 
      width = 8*500, height =8*300, units = "px", res=350,pointsize = 15, bg = "white", restoreConsole = TRUE)
  #p
  dev.off()
  ############# imputation missing calls ###########################################
  ####################################################################################
  #ana.dat=ana.dat %>% select(-starts_with('LMO'),which(names(ana.dat) %in% core.id))
  fd.dat=ana.dat %>% filter(food!='human')
  fd.dat$food=factor(fd.dat$food)
  
  
  ###### number of isolates from solved outbreaks ############
  tab=xtabs(~food+source,fd.dat)
  sum(tab[,dimnames(tab)$source=='human'])
  dim(fd.dat)[1]-sum(tab[,dimnames(tab)$source=='human'])
  return(fd.dat)
}

gen_tab_ind_pred_prob_orpb=function(pred){
  # produce sheet for ORPB verification #####################
  pred.db=pred[,unique(pred.fd),obs]
  names(pred.db)[2]='predicted food'
  pred.db$pred.prob=pred[,max(prob),obs][,V1]
  
  pred.db$wgs.id=ana.dat$WGS_id.x[row.names(ana.dat) %in% pred.db$obs]
  pred.db$miss=round(perc_inno_h_iso$perc_inno[perc_inno_h_iso$iso %in% sel.id]*100,2)
  pred.db=pred.db[,c(4,2,3,5)]
  tab=head(pred.db)
  print(tab)
  print(xtable(tab,digits = 2),'\\\\cdc.gov\\private\\M318\\vhg8\\tab1.html',type='HTML',include.rownames = F)
  
  #########################################################################
}

gen_sum_tab_attr=function(pred,sel.id,h.dat){
  ### summary of predicted foods### 
  pred[,unique(pred.fd),obs][,table(V1)]
  db=data.frame(table(train.df$food)/nrow(train.df))
  names(db)=c('food','training data')
  db$human=pred[,unique(pred.fd),obs][,table(V1)]/nrow(h.pred)
  
  #exclude duplicated human isolates
  ht=0.004 ##### proportion difference within clusters defined by ht, height parameter in cutree
  genes=data.frame(sapply(h.dat[row.names(h.dat) %in% sel.id,grep('LM',names(h.dat))],as.factor))
  row.names(genes)=sel.id
  cl=daisy(genes,metric='gower')
  hc=hclust(cl,method='complete') # what kind of linkage
  
  ct=cutree(hc,h=c(ht))# assign clade memberships of individual isolates
  length(unique(ct)) ### number of clades
  ct.mem=sapply(unique(ct),function(x){
    w=which(match(ct,x)==1) ### find individuals belong to cluster x 
    nam=names(ct)[w] #### find names of individuals belong to x
    if(length(nam)>1) {
      set.seed(45)
      id=sample(nam,1) 
    }else id=nam
  })
  
  temp=pred[pred$obs %in% ct.mem,unique(pred.fd),obs] ### summarize selected representative isolates
  nrow(temp) #### number of representative isolates
  db$human_exc_dup=temp[,table(V1)]/nrow(temp)
  print(db)
  print(xtable(db,digits = 2),'\\\\cdc.gov\\private\\M318\\vhg8\\Whole Genomo Sequence\\Lm WGS\\prediction of human isolates\\summary_tab.html',type='HTML',include.rownames = F)
  
}

cal_conf_attr=function(n.rep,pred,train=NULL){
  pred=data.table(pred)
  #### calculate intervals of uncertainty of predicted class ####################################
  #oob.accu=cf.seq$res[which.max(cf.seq$res$overall)] ### adjust estimation by oob accuracy of training data
  #oob.accu=unlist(oob.accu)[match(levels(train.df$food),names(oob.accu))]
  oob.accu=rep(1,length(unique(train.df$food)))

  pred=pred[order(as.numeric(pred$obs))]
  temp=pred[,list(perm.ind=c(rmultinom(n.rep,1,prob = prob/sum(prob) )),
                  id=rep(c(1:n.rep),each=length(unique(train.df$food)))),obs]
  temp[,perm.class:= ifelse(perm.ind==1,unique(pred$food)[which(perm.ind==1)],NA),by=c('obs','id')]
  
  ####### calculate accuracy uncertainty of bootstrap samples ############
  if(!is.null(train)){
    g=temp[,.(obs=unique(obs),pred=perm.class[!is.na(perm.class)]),id]
    sp=split(g,g$id)
    gg=lapply(sp,function(x){
      s=sen_spe(x$pred,train.df$food)
      #print(diag(table(x$pred,train.df$food))/table(train.df$food))
      list(accu=c(s$accurate_rate,s$over_accu),se=s$se,sp=s$sp,kappa=s$kappa)
    })
    my.f=function(x) { c(mean(x),quantile(x,probs=c(.025,.975)))}
    accu=apply(do.call(rbind,lapply(gg,function(x) x$accu)),2,function(y) my.f(y))
    se=apply(do.call(rbind,lapply(gg,function(x) x$se)),2,function(y) my.f(y))
    sp=apply(do.call(rbind,lapply(gg,function(x) x$sp)),2,function(y) my.f(y))
    k=apply(do.call(rbind,lapply(gg,function(x) x$kappa)),2,function(y) my.f(y))
  } 
  #######################################################################################################
  
  #temp[,agree:=ifelse(perm.class==train.df$food[obs],1,0)]
  temp1=temp[,.(attr=c(table(perm.class)/sum(table(perm.class))),food=unique(pred$food)),id]
  #temp1[,quantile(attr,probs=c(.025,.5,.975)),food]
  temp2=temp1[,.(attr=c(mean(attr),sd=sd(attr),quantile(attr,probs=c(.025,.975))),nam=c('mean','sd','lower','upper')),food]
 # over=temp[,.(over=c(sum(agree,na.rm=T)/nrow(train.df))),id]
  
  temp2=temp2 %>% spread(key=nam,value=attr) 
  temp2$food=factor(temp2$food,levels=levels(reorder(temp2$food,-temp2$mean)))
  if(!is.null(train))
    (list(attribution=temp2,conf95_accuracy=accu,conf95_se=se,conf95_sp=sp,conf95_kappa=k))
  else (list(attribution=temp2))
}

plot_panel_pred_prob_ind=function(pred,train=NULL){
  ###########################################################     
  sp=split(pred,pred$pred.fd)
  g=lapply(1:length(sp),function(d){ #### reorder obs based on pred prob of THE food
    prob=unlist(tapply(sp[[d]]$prob,sp[[d]]$obs,function(x) rep(x[d],length(x))))  # [ZC: replicate the largest probability for a given isolate X times (X is the number of food categories). 
                                                                                   #     [As a result, in the data set for pred.fd == "fruit", the prob for fruit was replicated X times
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
  
  
  ##### plot predicted probablity of foods for each isolate
  #png(filename= '\\\\cdc.gov\\private\\M318\\vhg8\\Whole Genomo Sequence\\Lm WGS\\prediction of human isolates\\pred_prob.png', 
   #   width = 8*420, height =8*300, units = "px", res=400,pointsize = 15, bg = "white", restoreConsole = TRUE)
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

plot_fantree_pred_fd=function(pred,sel.fd,imp.human){
  ##################### plot clusters of predicted foods #########################################
  source('\\\\cdc.gov\\private\\M318\\vhg8\\r scripts\\wgs\\wgMLST funs.R')
  imp.human$food=pred[,food[which.max(prob)],obs][,V1] ### change the randomly assigned foods to predicted one
  
  
  plot.col=c('brown','blue','red','green','pink')[1:length(sel.fd)]
  plot.id='food'
  
  #f=1:length(sel.fd)
  #lapply(1:length(sel.fd),function(f) {
  #temp=imp.human[imp.human$food==sel.fd[f],]
  temp=imp.human
  temp1=temp[,grep('LM',names(temp))]
  #fd=sel.fd[f]
  #p.col=plot.id[f]
  fd=sel.fd
  p.col=plot.id
  
  lab.col=data.frame(id=row.names(temp1), plot.lab=as.character(temp$food),
                     food=temp$food,col=plot.col[match(temp$food,unique(temp$food))])
  # hc_plot(temp1,lab.col,draw_ori = TRUE)
  hh=hc_plot(temp1,lab.col,ylim=c(0,1),draw_ori =NULL)
  
  
  
  hh$fan.tree
  
  #dev.off()
  return(hh$fan.tree)

}

plot_pred_prob_fd=function(ctrees){
######## plot conditional probabilities of foods of individual isolates #################
pred.prob=predict(ctrees,OOB=T,type='prob')
pred=data.frame(prob=unlist(pred.prob),obs=rep(rep(1:nrow(train.df),each=length(unique(train.df$food)))),
                food=rep(levels(train.df$food),nrow(train.df)),
                pred.fd=rep(predict(ctrees,OOB = T),each=length(unique(train.df$food))),
                obs.fd=rep(train.df$food,each=length(unique(train.df$food))))

sp=split(pred,pred$pred.fd)
g=lapply(1:length(sp),function(d){ #### reorder obs based on pred prob of THE food
  prob=unlist(tapply(sp[[d]]$prob,sp[[d]]$obs,function(x) rep(x[d],length(x))))
  dat=sp[[d]][order(prob,decreasing=T),]
})
pred1=do.call(rbind,g)
pred1$obs=factor(pred1$obs,levels=unique(pred1$obs)) #### keep the correct order 
pred1$color=ifelse(pred1$pred.fd=='dairy','red',ifelse(pred1$pred.fd=='fruit','grey',
                                                       ifelse(pred1$pred.fd=='meat','yellow',ifelse(pred1$pred.fd=='sea food','green','purple'))))

pred1$color.obs=ifelse(pred1$obs.fd=='dairy','red',ifelse(pred1$obs.fd=='fruit','grey',
                                                          ifelse(pred1$obs.fd=='meat','yellow',ifelse(pred1$obs.fd=='sea food','green','purple'))))
pred1$correct=ifelse(pred1$pred.fd==pred1$obs.fd,'black','red')

pred1$food=ifelse(pred1$food=='sea food','seafood',ifelse(pred1$food=='veg','vegetable',as.character(pred1$food)))

pred1$pred.fd=ifelse(pred1$pred.fd=='sea food','seafood',ifelse(pred1$pred.fd=='veg','vegetable',as.character(pred1$pred.fd)))



tick.col=tapply(pred1$color,as.factor(pred1$obs),function(x) unique(x))
tick.col.obs=tapply(pred1$color.obs,as.factor(pred1$obs),function(x) unique(x))
tick.col.correct=tapply(pred1$correct,as.factor(pred1$obs),function(x) unique(x))
g=unlist(sapply(table(predict(ctrees,OOB=T)),function(x) 1:x))
rug.df=data.frame(x=g,indicator=tick.col.correct,food=tapply(pred1$pred.fd,as.factor(pred1$obs),function(x) unique(x)))
rug.df$y=ifelse(rug.df$indicator=='red',NA,0.1)
rug.df$color=ifelse(rug.df$indicator=='black','red',NA)
#rug.df$color=ifelse(is.na(rug.df$color),NA,ifelse(rug.df$food=='fruit','grey',ifelse(rug.df$food=='meat','yellow',
 #         ifelse(rug.df$food=='seafood','green',ifelse(rug.df$food=='vegetable','purple','red')))))
#tick.col=tick.col[order(as.numeric(names(tick.col)))]
ggplot(pred1,aes(obs,prob))+geom_bar(aes(fill=food),stat='identity',width = 1)+
  theme(axis.line.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_rug(aes(obs), sides="b")+facet_wrap(~ pred.fd, scales='free_x')
  #geom_rug(aes(obs), subset(sam.df, obs.fd==pred.fd), sides="b")+



pred1 %>%
  mutate_at(vars(food, pred.fd, obs.fd), funs(recode(., `sea food` = "seafood"))) %>% 
  mutate_at(vars(food, pred.fd, obs.fd), funs(recode(., `veg` = "vegetable"))) ->
  pred1

p=ggplot(pred1,aes(obs,prob))+geom_bar(aes(fill=food),stat='identity',width = 1)+
  scale_fill_manual(values = c('red','grey','yellow','green','purple'))+
  xlab('isolates')+ylab('predicted probability')+
  theme(text = element_text(size=16),axis.line.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=10),
        axis.ticks.x=element_blank(),legend.position = c(0.8,0.3),legend.text = element_text(size=14),
        legend.title = element_text(size=14),strip.text.x = element_text(size = 15))+
  geom_rug(aes(obs), subset(pred1, obs.fd==pred.fd), sides="b")+
  facet_wrap(~ pred.fd, scales='free_x')
#theme(axis.text.x=element_text(color=tick.col,angle=90))+facet_grid(~pred.fd,scale='free',space='free_x')

return(p)
}

parallel.seq.forest=function(dat,imp,my.seq,ntree){

  res=foreach(n=1:length(my.seq),.packages='randomForest') %dopar% {
    d=dat[,names(dat) %in% c('food',as.character(imp[1:my.seq[n]]))]
    m=randomForest(y=d$food,x=d[,names(d)!='food'],ntree=ntree)
    list(1-m$confusion[,which(dimnames(m$confusion)[2][[1]]=='class.error')],overall=sum(diag(m$confusion))/sum(m$confusion))
  }

  dat=data.table(food=names(unlist(res)),accu=unlist(res),loci=rep(my.seq,each=(length(unique(dat$food))+1)))
  
  dat$food=factor(dat$food,levels=levels(reorder(dat$food,-dat$accu)))

  dat.w=(spread(dat,food,accu))
  
  ltype=ifelse(levels(dat$food)=='overall','solid','longdash')
  mycol=c( "#B79F00", "#00BA38","#F8766D","#619CFF", "#F564E3","#00BFC4") #### ggplot default colours by ggplot_build(p)$data[[1]]colour
  
  p=ggplot(dat,aes(loci,accu,group=food,color=food))+geom_line(aes(linetype=food),size=1)+ylab('prediction accuracy')+
    coord_trans(x = "log")+scale_x_continuous(breaks=c(5,10,20,30, 50, 70,100,200,500,1000,2000))+xlab('number of the ranked wgMLST loci')+
    scale_linetype_manual(values=ltype)+scale_color_manual(values=mycol)+
    theme(legend.position='right',
          text=element_text(size = 15),axis.text = element_text(size=10),legend.text = element_text(size=14),
          legend.title = element_text(size=14))+guides(fill = guide_legend(keywidth = 1, keyheight = 1), linetype=guide_legend(keywidth = 3, keyheight = 1),
                                                       colour=guide_legend(keywidth = 3, keyheight = 1))
  
  return(list(res=dat.w,p=p))
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
  
  #levels(dat$food)[c(4,6)]=c("vegetable", "seafood")
  #dat.w=(spread(dat,food,accu))
  
  return(list(res=dat))
}

seq.cforest=function(dat,imp,my.seq,ntree){
  library(doParallel)
  library(foreach)

  cl=makeCluster(22)
  registerDoParallel(cl) #number of cores on the machine
  
  
  res=foreach(n=my.seq,.packages='party') %dopar% {
    m=cforest(food~.,dat[,names(dat) %in% c('food',as.character(imp[1:n]))],
              control=cforest_control(ntree=ntree,mtry=sqrt(n),
              minbucket=3))
    print(n)
    tab=table(dat$food,predict(m,OOB=T))
    #pred=diag(tab)/table(dat$food)
    t.pos=diag(tab)
    f.neg=sapply(1:length(unique(dat$food)),function(x) sum(tab[x,])-tab[x,x])
    sen=c(round(t.pos/(t.pos+f.neg),2))
    over_sen=round(sum(t.pos)/sum(t.pos+f.neg),2)
    
    f.pos=sapply(1:length(unique(dat$food)),function(x) sum(t(tab)[x,])-t(tab)[x,x])
    t.neg=sapply(1:length(unique(dat$food)),function(x) sum(tab)-sum(tab[x,]))
    spe=c(round(t.neg/(t.neg+f.pos),2))
    over_spe=round(sum(t.neg)/sum(t.neg+f.pos),2)
    #over_accu=(sum(t.pos)+sum(t.neg))/(sum(t.pos)+sum(t.neg)+sum(f.pos)+sum(f.neg))
    #prev=table(dat$food)/nrow(dat)
    #dat=data.frame(rbind(tab,apply(tab,1,sum),sen,spe))
    #list(food=c(dimnames(tab)[[1]],'overall'),pred_accu=c(sen*prev+spe*(1-prev),over_accu),
         #se=c(sen,over_sen),sp=c(spe,over_spe))
    list(food=c(dimnames(tab)[[1]],'overall'),se=c(sen,over_sen),sp=c(spe,over_spe))
  }
  
  dat=data.table(food=c(sapply(res,function(x) unlist(x$food))),
                 se=c(sapply(res,function(x) unlist(x$se))),sp=c(sapply(res,function(x) unlist(x$sp))),loci=rep(my.seq,each=length(res[[1]]$food)))
  
  #levels(dat$food)[c(4,6)]=c("vegetable", "seafood")
  #dat.w=(spread(dat,food,accu))
  
  ltype=ifelse(levels(dat$food)=='overall','solid','longdash')
  my.col=c('blue','green','black','red','brown','pink')
  
  dat$food=factor(dat$food,levels=levels(reorder(dat$food,-dat$se)))
   
  p.se=ggplot(dat,aes(loci,se,group=food,color=food))+geom_line(aes(linetype=food),size=1)+ylab('estimated sensitivity')+
    coord_trans(x = "log")+scale_x_continuous(breaks=c(5,10,20,30, 50, 70,100,200,500,1000,2000))+xlab('number of the ranked wgMLST loci')+
    scale_linetype_manual(values=ltype)+scale_color_manual(values=my.col)+
    theme(legend.position='right',
          text=element_text(size = 15),axis.text = element_text(size=10),legend.text = element_text(size=14),
          legend.title = element_text(size=14))+guides(fill = guide_legend(keywidth = 1, keyheight = 1), linetype=guide_legend(keywidth = 3, keyheight = 1),
                                                       colour=guide_legend(keywidth = 3, keyheight = 1))

  return(list(res=dat,p.se=p.se))
}

for_EDLB_only=function(si){
  ###### temporary data for EDLB #####################
  df=ana.dat[,names(ana.dat) %in% c('Key','Key','Outbreak','food')]
  df$Outbreak=ifelse(df$Outbreak=='NA','not_cluster',df$Outbreak)
  df$steven=ifelse(df$Key %in% dup$Key,'rep',NA)
  df$auto=ifelse(df$Key %in% si,'rep',NA)
  df=df[order(df$Outbreak),]
  sp=split(df,df$Outbreak)
  df$test=unlist(sapply(sp,function(x){
    if(sum(is.na(x$steven))<nrow(x) & x$Outbreak!='not_cluster') x$steven else x$auto
    
  }))
  si.steven=df$Key[!is.na(df$test)]
  
  temp=df %>%  filter(Key %in% si.steven) %>% select(starts_with('LMO'))
  temp <- temp[,apply(temp,2,function(x)sum(is.na(x))/length(x))<.3]
  
  gg=sapply(temp,function(x) { #### note this module only works with numerical variables
    tab=table(x)
    if(length(tab)>53){
      low_f=as.numeric(names(tab)[tab<5])
      x[unlist(x) %in% low_f]=9999
    }
    if(length(table(x))>53){
      tab=table(x)
      temp=as.numeric(names(tab)[tab<4])
      x[unlist(x) %in% temp]=9998
    }
    if(length(table(x))>53){
      tab=table(x)
      temp=as.numeric(names(tab)[tab<5])
      x[unlist(x) %in% temp]=9997
    }
    if(length(table(x))>53){
      tab=table(x)
      temp=as.numeric(names(tab)[tab<6])
      x[unlist(x) %in% temp]=9996
    }
    x
  })
  x.st=data.frame(apply(gg,2,as.factor))
  x.st=cbind(x.st,sapply(dat[dat$Key %in% si.steven, names(dat) %in% conf],as.factor))
  y.st=as.factor(dat$food[dat$Key %in% si.steven])
  
  set.seed(3434)
  steven.imp=rfImpute(y=y.st,x=x.st)
  names(steven.imp)[1]='food'
  
  rf.st=randomForest(y=steven.imp$food,x=steven.imp[,-1],importance=T,ntree=1000)
  
  df=data.frame(imp=randomForest::importance(rf)[,'MeanDecreaseAccuracy'])
  df$id=rownames(df)
  df$method='autoselection'
  df=df[order(df$imp,decreasing = T),]
  imp.auto=df$id
  df=df[1:60,]
  df = df %>% filter(!id %in% conf)
  
  df.st=data.frame(imp=randomForest::importance(rf.st)[,'MeanDecreaseAccuracy'])
  df.st$id=rownames(df.st)
  df.st$method='steven'
  df.st=df.st[order(df.st$imp,decreasing = T),]
  imp.st=df.st$id
  df.st=df.st[1:60,]
  df.st = df.st %>% filter(!id %in% conf)
  
  comb=rbind(df,df.st)
  comb$id=factor(comb$id,levels=levels(reorder(comb$id,(comb$imp))))
  
  ggplot(comb,aes(id,imp,colour=method))+geom_point()+coord_flip()
  
  temp=randomForest(y=steven.imp$food,x=steven.imp[,names(steven.imp) %in% imp.st[1:500]],ntree=1000)
  temp0=randomForest(y=rf.dat$food,x=rf.dat[,names(rf.dat) %in% imp.auto[1:500]],ntree=1000)
  
  pred=data.frame(accuracy=1-c(temp0$confusion[,6],temp$confusion[,6]),food=rep(
    row.names(temp$confusion),2),method=rep(c('autoselection','steven'),each=5))
  ggplot(pred,aes(food,accuracy,colour=method))+geom_point()
}

grp_alle_impu=function(dat,conf){
  #F_dat <-dat %>%  filter(Key %in% si) %>% select(starts_with('LMO'))
  
  #df=data.frame(alle=sapply(cs_dat,function(x) length(table(x))))
  #df$gene=ifelse(rownames(df)%in% mlst.id,'MLST',ifelse(rownames(df)%in% core.id,'core','non-core'))
  #df$gene=factor(df$gene,levels=c('MLST','core','non-core'))
  #png(filename= '\\\\cdc.gov\\private\\M318\\vhg8\\whole genomo sequence\\data\\new data\\num_alle_call per locus.png', width = 8*400, 
    #  height =8*400, units = "px", res=400,pointsize = 15, bg = "white", restoreConsole = TRUE)
  #ggplot(df,aes(1:dim(df)[1],alle,color=gene))+geom_point()+xlab('')+ylab('number of allelles')+ylim(0,80)+
    #scale_color_manual(values=c('green',"red",'blue'))+ ggtitle('Number of allelles before grouping')
  #dev.off()
  
  ### combine rare alleles
  #temp=dat %>%  filter(Key %in% si) %>% select(starts_with('LMO'))
  temp=dat %>%  dplyr::select(starts_with('LMO'))
  gg=data.frame(sapply(temp,function(x) { #### note this module only works with numerical variables
    tab=table(x)
    if(length(tab)>53){
      or=order(table(x),decreasing=T)
      low_f=as.numeric(names(tab)[or[-c(1:52)]])
      x[unlist(x) %in% low_f]=9999
    }
    
    x
  }))
  
  #df=data.frame(alle=apply(gg,2,function(x) length(levels(as.factor(x)))))
  #df$gene=ifelse(rownames(df)%in% mlst.id,'MLST',ifelse(rownames(df)%in% core.id,'core','non-core'))
  #df$gene=factor(df$gene,levels=c('MLST','core','non-core'))
  #ggplot(df,aes(1:dim(df)[1],alle,color=gene))+geom_point()+xlab('')+ylab('number of allelles')+ylim(0,80)+
    #scale_color_manual(values=c('green',"red",'blue'))+ ggtitle('Number of allelles after grouping')
  
  x=data.frame(sapply(gg,as.factor))
  #x=cbind(x,sapply(dat[,names(dat) %in% conf],as.factor))
  y=as.factor(dat$food)
  
  set.seed(3434)
  core.imp=rfImpute(y=y,x=x)

  ################################# for EDLB ####################################
  
  ##### imputation before selection of representative isolates to maximize information in duplicates for imputation
  
  
  return(core.imp)
}

dissim_range=function(dat){
  temp <-dat %>%  filter(Key %in% si) 
  temp=data.frame(sapply(temp,as.factor))
  sp=split(temp,temp$food)
  b=sapply(sp,function(x) {
    xx=x %>% dplyr::select(starts_with('LMO'))
    cl=daisy(x,metric='gower')
    c(mean(cl),quantile(cl,probs=c(0,0.5,1)))
  })
  tb=data.frame(t(b))
  names(tb)=c('mean','minimum','median','maximum')
  return(tb)
}

est.pair.dist.clust=function(dat){
  library(StatMatch)
  tab=table(dat$Outbreak)
  ot=names(tab)[tab>1]
  ot=ot[ot!='NA']
  temp=dat %>% filter(Outbreak %in% ot)
  temp$Outbreak=factor(temp$Outbreak)
  temp=split(temp,temp$Outbreak)
  g=sapply(1:length(temp),function(x) {
    print(x)
    xx=select(temp[[x]],starts_with('LMO')) 
    c(gower.dist(xx))
  })
  gg=data.frame(out=rep(names(temp),sapply(temp,nrow)^2),value=unlist(g)*100)
  gg$out=factor(gg$out,levels=levels(reorder(gg$out,gg$value)))
  p=ggplot(gg,aes(out,value,colour=out))+geom_point()+geom_hline(yintercept = 0.4,linetype='dashed')+
    scale_y_log10() +coord_flip()+ylab('percent difference (log)')+xlab('Outbreak ID')+
    theme(legend.position = 'none',axis.text.x = element_text(size=5,angle = 90,vjust=0.5, hjust = 1))
  return(p)
}

miss.plot=function(dat){
  wg.mlst=dat[,grep('LMO',names(dat))]
  miss=data.frame(sapply(wg.mlst,function(x) sum(is.na(x))/length(x)*100))
  names(miss)='missing'
  miss$gene=ifelse(rownames(miss)%in% mlst.id,'MLST',ifelse(rownames(miss)%in% core.id,'core','non-core'))
  miss$gene=factor(miss$gene,levels=c('MLST','core','non-core'))
  
  #png(filename= '\\\\cdc.gov\\private\\M318\\vhg8\\whole genomo sequence\\data\\new data\\overall_missing.png', 
     # width = 8*400, height =8*300, units = "px", res=300,pointsize = 15, bg = "white", restoreConsole = TRUE)
  
  p=ggplot(miss,aes(1:dim(miss)[1],missing,color=gene))+geom_point()+xlab('')+
    scale_color_manual(values=c('green',"red",'blue'))+ ggtitle('percentage of missing allelle calls')
  #dev.off()
  return(p)
  #################################
}

vis.food.iso.clust=function(dat,sel.fd,ylim,draw_ori=NULL,plot.id,plot.col){
  temp=dat[dat$food %in% sel.fd,]
  temp$food=factor(temp$food)
  temp1=temp[,grep('LM',names(temp))]
  lab.col=data.frame(id=as.character(temp$SRR_ID), plot.lab=if(plot.id=='Key') as.character(temp$Key) 
          else as.character(temp$food),food=temp$food,col=plot.col[match(temp$food,levels(temp$food))])
  # hc_plot(temp1,lab.col,draw_ori = TRUE)
  hh=hc_plot(temp1,lab.col,ylim=ylim,draw_ori =NULL)
  #hh$gd
  #png(filename= '\\\\cdc.gov\\private\\M318\\vhg8\\whole genomo sequence\\data\\new data\\ice cream cluster.png', 
      #width = 8*500, height =8*300, units = "px", res=350,pointsize = 15, bg = "white", restoreConsole = TRUE)
  #plot(hh$h,main=sel.fd,sub='',xlab='',ylab='proportional difference in wgMLST loci',cex=0.9)
  #rect.hclust(hh$hc,h=0.004,,border=2:length(table(cutree(hh$hc,h=0.005))))
  
  return(list(gd=hh$gd,fan.tree_dot=hh$fan.tree_dot, fan.tree_label=hh$fan.tree_label, legend_coord = hh$legend_coord))
  #dev.off()
}

hc_plot=function(dat,lab.col,ylim,draw_ori=NULL){
  ## this function plots hc using labels and colors defined in lab.col db
  #lab.col=data.frame(sapply(lab.col,function(x) as.factor(x)))
  dat=data.frame(sapply(dat,as.factor))
  cl=daisy(dat,metric='gower')
  hc=hclust(cl,method = 'complete')
  
  hc$labels=as.character(lab.col$id)
  
  dend <- as.dendrogram (hc) 
  
  set.seed(13)
  col.dat=data.frame(nam=names(table(lab.col$col)),freq=as.numeric(table(lab.col$col)),
                     color=sample(colors(),length(table(lab.col$col)),replace=F),stringsAsFactors = F)
  #color=c('red','blue','blueviolet','darkgoldenrod','green','brown','cornflowerblue','aquamarine4',
  
  
  dendr <- dendro_data(hc, type="rectangle") 
  
  if(!is.null(draw_ori)){
    dendr$segments$y=dendr$segments$y*dim(dat)[2]
    dendr$segments$yend=dendr$segments$yend*dim(dat)[2]
  } 
  
  
  lab=dendr$labels
  
  
  g=left_join(lab,lab.col,by=c('label'='id'))
  #gg=inner_join(g,col.dat,by=c('col'='nam'))
  
  gd=ggplot() + geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
    geom_point(aes(x=1:dim(g)[1],y=0),color='blue')+
    geom_text(data=dendr$labels, aes(x=x, y=y, label=g$plot.lab,hjust=0, color=g$col), size=2) +
    coord_flip(ylim=ylim)+scale_y_reverse(expand=c(0.2, 0))+
    scale_color_discrete(labels=g$food[match(levels(as.factor(plot.col)),g$col)])+
    guides(color=guide_legend("source of isolates"))+ylab('proportion of difference in wgMLST')+xlab('')+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),panel.grid=element_blank(),panel.background=element_rect(fill="white"))
  #scale_color_manual(name="",values=as.character(col.dat$color),labels=as.character(col.dat$nam))+theme(legend.position="bottom")+
  #theme(axis.line.y=element_blank(),
  #axis.ticks.y=element_blank(),
  #axis.text.y=element_blank(),
  #axis.title.y=element_blank(), axis.title.x=element_blank(),
  #panel.background=element_rect(fill="lightgrey"),
  #theme(panel.grid=element_blank())
  
  
  ape.tree=(as.phylo(hc))
  df=data.frame(id=ape.tree$tip.label)
  df1=left_join(df,g,by=c('id'='label'))
  #df2=data.frame(id=ape.tree$tip.label,food=df1$food,val=rnorm(1196,sd=3))
  #ape.tree$tip.label=as.character(df1$food)
  #ape.tree$tip.color=as.character(df1$col)
  df2=data.frame(id=ape.tree$tip.label,food=df1$id)
  
  p0=ggtree(ape.tree)+geom_tippoint(color=df1$col)
  p1=facet_plot(p0,panel='host',data=df2,geom=geom_point,aes(x=-2),size=0.7)
  #+scale_color_manual(values=my.col)
  library(grid)
  library(gtable)
  
  gt = ggplot_gtable(ggplot_build(p1))
  gtable_show_layout(gt) # will show you the layout - very handy function
  gt # see plot layout in table format
  gt$layout$l[grep('panel-2-1', gt$layout$name)] # you want to find the column specific to panel-2
  gt$widths[7] = 0.01*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
  grid.draw(gt) # plot with grid draw
  
  plot(ape.tree,tip.color=as.character(df1$col),cex=0.2) ### plot.phylo function in ape package
  axisPhylo(side = 1, root.time = NULL, backward = TRUE)
  legend(x='topleft',cex=0.7,title='source of isolate',title.col='black',
         legend=levels(df1$food),text.col=as.character(df1$col[match(levels(df1$food),df1$food)]))
  phylo.tree=recordPlot()
  
  ape.tree$tip.label <- as.character(lab.col$food)
  par(mar = c(0, 0, 0, 2), xpd=TRUE)
  plot(ape.tree, type = "fan",tip.color=as.character(df1$col),cex=1.2, show.tip.label = FALSE)
  tiplabels(pch = 19, col = as.character(df1$col[match(levels(df1$food),df1$food)]), cex = 1)
  # legend(x='topright',cex=1.5, title='source of isolate',title.col='black', 
  #               legend=levels(df1$food),text.col=as.character(df1$col[match(levels(df1$food),df1$food)]),
  #               inset = c(-0.2, 0.1), y.intersp = 0.4, x.intersp = 0.1, text.width = 0)
  leg <- legend(x='topright',cex=1.8, 
         legend=levels(df1$food),text.col=as.character(df1$col[match(levels(df1$food),df1$food)]),
         pch = 19, col = as.character(df1$col[match(levels(df1$food),df1$food)]), 
         x.intersp = 0.2, y.intersp = 0.35, plot = FALSE)
  leftx <- leg$rect$left + 0.2
  rightx <- leftx + leg$rect$w * 0.52
  topy <- leg$rect$top * 0.75
  bottomy <- topy - leg$rect$h * 0.9
  legend(x = c(leftx, rightx), y = c(topy, bottomy), cex=1.8,  bty = "n",
         legend=levels(df1$food),text.col=as.character(df1$col[match(levels(df1$food),df1$food)]),
         pch = 19, col = as.character(df1$col[match(levels(df1$food),df1$food)]), 
         x.intersp = 0.2, y.intersp = 0.35)
  text(x = leftx + 0.01, y = topy - 0.02, "Source of isolate", cex = 1.8, pos = 4)
  
  fan.tree_dot=recordPlot()
  
  
  plot(ape.tree, type = "fan",tip.color=as.character(df1$col),cex=1.2)

  leg <- legend(x='topright',cex=1.8,  
                legend=levels(df1$food),text.col=as.character(df1$col[match(levels(df1$food),df1$food)]),
                xjust = 0, x.intersp = 0, y.intersp = 0.35, plot = FALSE)
  leftx <- leg$rect$left + 0.2
  rightx <- leftx + leg$rect$w * 0.52
  topy <- leg$rect$top * 0.75
  bottomy <- topy - leg$rect$h * 0.9
  legend(x = c(leftx, rightx), y = c(topy, bottomy), cex=1.8, bty = "n", 
         legend=levels(df1$food),text.col=as.character(df1$col[match(levels(df1$food),df1$food)]), 
         xjust = 0, x.intersp = 0.2, y.intersp = 0.35)
  text(x = leftx + 0.02, y = topy - 0.02, "Source of isolate", cex = 1.8, pos = 4)
  fan.tree_label=recordPlot()
  
  return(list(gd=gd,dend=dend,hc=hc,fan.tree_dot=fan.tree_dot, fan.tree_label=fan.tree_label, phylo.tree=phylo.tree, legend_coord = leg))
}


comp_two_sel=function(si){
  
 ###########################################################################  
  
  temp=ana.dat[ana.dat$Outbreak %in% unique(dup$Outbreak),] ## check with steven's selection
  temp=ana.dat
  temp$my.sel=ifelse(temp$Key %in% si,'red','grey')
  #temp$steven.sel=ifelse(temp$Key %in% si,'green','grey')
  temp$steven.sel=ifelse(temp$Key %in% dup$Key,'green','purple')
  temp$out.id=ifelse(temp$Outbreak!='NA',temp$Outbreak,temp$food)
  
  fd=names(table(temp$food))
  sapply(fd,function(f) {
    print(f)
    temp.m=temp[temp$food %in% f,]
    temp1=temp.m[,grep('LM',names(temp.m))]
    
    ### note two colours for two selection methods
    lab.col=data.frame(lab=as.character(temp.m$Key),col1=temp.m$steven.sel,col2=temp.m$my.sel,
                       lab1=temp.m$out.id)
    
    #lab.col=data.frame(lab=as.character(temp.m$Key),col1=temp.m$steven.sel,col2=temp.m$my.sel)
    
    hh=comparison(temp1,lab.col,ylim=c(0,0.02),draw_ori=NULL)
    
    p=hh$gd+theme(legend.position = 'right',axis.text = element_text(size=14),axis.text.y=element_blank())+labs(colour=f)+geom_hline(yintercept = ht,linetype='dashed')+
      ylab('proportional difference in wgMLST')+xlab('')
      
    p_iden=hh$gd_iden+theme(legend.position = 'right',axis.text = element_text(size=14),axis.text.y=element_blank())+labs(colour=f)+geom_hline(yintercept = ht,linetype='dashed')+
      ylab('proportional difference in wgMLST')+xlab('')
      
    png(filename= paste('\\\\cdc.gov\\private\\M318\\vhg8\\whole genomo sequence\\data\\new data\\auto_sel_comp_',paste(f,'.png',sep=''),sep=''),, width = 8*400, 
        height =8*500, units = "px", res=500,pointsize = 15, bg = "white", restoreConsole = TRUE)
    print(p)
    dev.off()
    
    png(filename= paste('\\\\cdc.gov\\private\\M318\\vhg8\\whole genomo sequence\\data\\new data\\sel_iso_',paste(f,'.png',sep=''),sep=''),, width = 8*400, 
        height =8*500, units = "px", res=500,pointsize = 15, bg = "white", restoreConsole = TRUE)
    print(p_iden)
    dev.off()
    
    
  })
  
  sel_comparison=function(dat,lab.col,ylim=NULL,draw_ori=NULL){
    
    dat=data.frame(sapply(dat,as.factor))
    cl=daisy(dat,metric='gower')
    hc=hclust(cl,method='complete')
    
    hc$labels=as.character(lab.col$lab)
    
    dend <- as.dendrogram (hc) 
    
    set.seed(13)
    col.dat=data.frame(nam=names(table(lab.col$col2)),freq=as.numeric(table(lab.col$col2)),
                       color=sample(colors(),length(table(lab.col$col2)),replace=F),stringsAsFactors = F)
    
    dendr <- dendro_data(hc, type="rectangle") 
    
    if(!is.null(draw_ori)){
      dendr$segments$y=dendr$segments$y*dim(dat)[2]
      dendr$segments$yend=dendr$segments$yend*dim(dat)[2]
    } 
    
    
    lab=dendr$labels
    
    g=left_join(lab,lab.col,by=c('label'='lab'))
    gg=inner_join(g,col.dat,by=c('col2'='nam'))
    
    
    gd= ggplot() + geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + geom_point(aes(x=1:dim(gg)[1],y=0),color=gg$col2)+
      geom_text(data=gg, aes(x=x, y=y, label=label,hjust=-0.1), color=gg$col2, size=2.5) +
      coord_flip(ylim=ylim)+scale_y_reverse(expand=c(0.2, 0))+labs(color='food')+theme(panel.grid=element_blank())
    
    gd_iden=ggplot() + geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
      geom_point(aes(x=1:dim(gg)[1],y=0),color=gg$col2)+
      geom_text(data=dendr$labels, aes(x=x, y=y, label=label,hjust=0), color=gg$col2, size=2)+
      coord_flip(ylim=ylim)+scale_y_reverse(expand=c(0.2, 0))+labs(color='food')+theme(panel.grid=element_blank())
    
    return(list(gd=gd,gd_iden=gd_iden))
  }
  
  
  #################################################################################
}
hc_prediction=function(dat,ht){
  genes=dat %>% select(starts_with('LM'))
  genes=data.frame(sapply(genes,as.factor))
  epi=dat %>% select(-starts_with('LM'))
  
  cl=daisy(genes,metric='gower')
  hc=hclust(cl,method='complete') # what kind of linkage
  hc$labels=epi$Key

  ct=cutree(hc,h=c(ht))# Cut interpretation for single linkage: for each point Xi, there is another point Xj in

  #plot(hc)
  #rect.hclust(hc,h=ht)
  
  ff=as.character(dat$food)
  t.ct=table(ct)
  m.ct=names(t.ct)[t.ct>1]
  tot.clust.n=length(m.ct)
  tot.clust.m=sum(t.ct[names(t.ct) %in% m.ct])
  ct.accu=lapply(m.ct,function(x){
    w=which(match(ct,x)==1) 
    nam=names(ct)[w]
    ss=ff[!is.na(match(epi$Key, nam))] ### observed food category
    if(length(nam)>1) {
      pred.food=names(which.max(table(ss)))
      list(fd=pred.food,accu=sum(match(ss,pred.food),na.rm=T),nn=length(nam))
    }})
  
  uni.fd=unique(sapply(ct.accu,function(x) x$fd))
  
  gp.fd=sapply(uni.fd,function(x) which(sapply(ct.accu,function(y) y$fd)==x))
  
  gp.accu=sapply(gp.fd,function(x) sapply(x,function(y) ct.accu[[y]]$accu))
  gp.nn=sapply(gp.fd,function(x) sapply(x,function(y) ct.accu[[y]]$nn))
  
  gp.sum=lapply(1:length(gp.fd),function(x) c(sum(gp.accu[[x]]),sum(gp.nn[[x]]),length(gp.nn[[x]])))
  
  gp.pred=data.frame(food=names(gp.fd),accuracy=sapply(gp.sum,function(x)x[1]/x[2]*100),
        n.clust=sapply(gp.sum,function(x)x[3]),n.ind=sapply(gp.sum,function(x)x[2]),stringsAsFactors = F)
  gp.pred[(1+dim(gp.pred)[1]),1]='overall'
  gp.pred[(dim(gp.pred)[1]),2]=sum(sapply(gp.sum,function(x)x[[1]]))/sum(sapply(gp.sum,function(x)x[[2]]))*100
  gp.pred[(dim(gp.pred)[1]),3]=tot.clust.n
  gp.pred[(dim(gp.pred)[1]),4]=tot.clust.m
  gp.pred=gp.pred[order(gp.pred$food),]
  
  return(gp.pred)

}


chk_ind_cls=function(sel_id){
  ####### HC individual clusters ###########################
  dt=dat[dat$Key %in% sel_id,grep('LM',names(dat))]
  dt1=data.frame(lab=sel_id,col=as.character(dat$SourceSite[dat$Key %in% sel_id]))
  h=hc_plot(dt,dt1)
  return(h)
}
sim_bhi=function(dat,fd,imp_index,start,n_step1,n_step2,ht,n_rep){
  
  #fc <- tapply(rownames(dat),fd, c)
  #fc <- annotationListToMatrix(fc, rownames(dat))
  
  #### bhi index for ranked loci ###########
  r_g=list()
  rr_rep=1
  print('imp_loci')
  for(i in 1:rr_rep){
    set.seed(i)

    r_g[[i]]=sapply(seq(start_n,dim(dat)[2],n_step1),function(x) {
      print(x)
      temp=dat[,names(dat)%in% imp_index[1:x]]
      hc=hclust(daisy(temp,metric='gower'),method='complete')
      cluster <- cutree(hc,h=ht)
      #c(x,BHI(cluster, fc))
      bhi1=rem.single.cluster.hc(dat,fd,cluster,ht)
      c(x,bhi1)
    })
  }
  
  bhi=data.frame(method='ranked',t(Reduce('+',r_g)/rr_rep),stringsAsFactors = F)
  names(bhi)[2:3]=c('loci','est')
  
  ##### for random selection of loci ############
  g=list()
  
  for(i in 1:n_rep){
    print(i)
    g[[i]]=sapply(seq(start_n,dim(dat)[2],n_step2),function(x) {
      
      set.seed(i)
      temp1=dat[,sample(1:dim(dat)[2],x,replace=F)]
      hc1=hclust(daisy(temp1,metric='gower'),method='complete')
      cluster1 <- cutree(hc1,h=ht)
      #c(x,BHI(cluster1, fc))
      bhi2=rem.single.cluster.hc(dat,fd,cluster1,ht)
      c(x,bhi2)
    })
  }
  
  ss= sapply(g,function(x) x[2,])
  cl=data.frame(method='random',loci=g[[1]][1,],t(apply(ss,1,function(x) c(mean(x),quantile(x, probs=c(.025,.975))))),
                stringsAsFactors = F)
  names(cl)[3:5]=c('est','Low95','upper95')
  
  
  return(list(bhi=bhi,cl=cl))
  
}

cal_BHI=function(dat,fd,imp_index,start,n_step,ht){
  ### calculate biological homogeneity index
  #n_step=20
  remain=1:dim(dat)[2]
  in_sample=0
  fc <- tapply(rownames(dat),fd, c)
  fc <- annotationListToMatrix(fc, rownames(dat))
  g=sapply(seq(start,dim(dat)[2],n_step),function(x) {
    temp=dat[,names(dat)%in% imp_index[1:x]]
    hc=hclust(daisy(temp,metric='gower'),method='complete')
    cluster <- cutree(hc,h=ht)
    ### remove clusters with only one isolate ###
    bhi1=rem.single.cluster.hc(dat,fd,cluster,ht)
    
    sam=sample(remain,n_step,replace=F)
    in_sample=c(in_sample,sam)
    remain=remain[-sam]
    #c(length(table(in_sample)),length(remain))
    
    temp1=dat[,sample(1:dim(dat)[2],x,replace=F)]
    hc1=hclust(daisy(temp1,metric='gower'),method='complete')
    cluster1 <- cutree(hc1,h=ht)
    bhi2=rem.single.cluster.hc(dat,fd,cluster1,ht)
    c(x,bhi1,bhi2)
    #c(x,BHI(cluster, fc),BHI(cluster1, fc))
  })
  return(g)
}

rem.single.cluster.hc=function(dat,fd,cluster,ht){
  w1=table(cluster)
  w1=names(w1)[w1==1]
  if(length(w1)>0){
  excl=sapply(w1,function(x) which(match(cluster,x)==1))
  new_dat=dat[-excl,]
  fd1=fd[-excl]
  } else {
    new_dat=dat
    fd1=fd
    }
  fc <- tapply(rownames(new_dat),fd1, c)
  fc <- annotationListToMatrix(fc, rownames(new_dat))
  hc=hclust(daisy(new_dat,metric='gower'),method='complete')
  nc <- cutree(hc,h=ht)
  return(BHI(nc, fc))
}




mymosaic=function(){
  library(vcd)
  eqal=data.frame(source=unique(as.character(merged$source)),SourceState='ZA')
  temp=merged[,c('source','SourceState')]
  temp=rbind(temp,eqal)
  temp$method='original'
  
  temp1=mod.dat[,c('source','SourceState')]
  temp1=rbind(temp1,eqal)
  temp1$method='rem_rep'
  
  temp0=rbind(temp,temp1)
  temp0=temp0[temp0$source!='complex food',]
  temp0$source=factor(temp0$source)
  
  png(filename= '\\\\cdc.gov\\private\\M318\\vhg8\\whole genomo sequence\\mosaicplot.png', width = 8*400, height =8*400, units = "px", res=300,pointsize = 15, bg = "white", restoreConsole = TRUE)
  mosaic(~SourceState+source|method,temp0,rot_labels=c(90,0,90,90))
  dev.off()
}

hmap=function(){
  #### food wide variability #############################################
  
  temp=data.frame(sapply(core,as.numeric))
  food.st <- temp[order(food$source),]
  fn=c(19,27,13,15,2,16,25,36,3,23,35,11,14,21,1,4)[1:length(unique(food))]
  fn=as.numeric(food$source[order(food$source)])
  #fn=as.numeric(as.factor(food[order(food)]))
  fd <- matrix(rep(fn,50),ncol=50)
  fd <- cbind(fd,as.matrix(food.st))
  max_l <- max(apply(fd,2,function(x)length(table(x))))
  my.col <- colors()[2:max_l]
  set.seed(1690)
  my.col <- colors()[sample(1:length(colors()),max_l,replace=F)]
  
  png(filename= '\\\\cdc.gov\\private\\M318\\vhg8\\whole genomo sequence\\var_food.png', width = 8*400, height =8*400, units = "px", res=200,pointsize = 15, bg = "white", restoreConsole = TRUE)
  #heatmap.2(fd,dendrogram='none',scale='none',col=redblue(35),main='variability between food')
  
  
  heatmap(fd,Rowv=NA,Colv=NA,labCol=NA,labRow=NA,scale='none',col=my.col)
  
  #heatmap(fd,Rowv=NA,Colv=NA,labCol=NA,labRow=NA,scale='none',,col=my.col,main='variability between food')
  table(food$source)
  
  dev.off()
}

sen_spe=function(pred,obs){
  #pred<-predict(tree1,newdata=comp.ue,type='class')
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

intro_var=function(){
  ### examine intra-variability of loci in food sources ###########################
  sp=split(merged,merged$source)
  lapply(sp,function(d) {
    d=d[,grep('LMO_',names(d))]
    d=data.frame(sapply(d,factor))
    tb=table(apply(d,2,function(x) length(table(x))))
    dp=duplicated(d)
    list(tb=tb,dp=dp)
  })
  ### for all loci compared, 0=missing, 1=identical, 2=two alleles, 3...

}
