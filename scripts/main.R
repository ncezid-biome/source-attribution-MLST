# Revised code by Weidong Gu December 8, 2023
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(hierfstat))
suppressPackageStartupMessages(library(poppr))
suppressPackageStartupMessages(library(mmod))
suppressPackageStartupMessages(library(randomForestSRC))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(multiROC))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(geiger))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(msgl))
suppressPackageStartupMessages(library(mlogit))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(xgboost))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(StatMatch))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gtable))
suppressPackageStartupMessages(library(vcd))

base_dir <- "//wsl.localhost/Ubuntu-18.04/home/gzu2/src/source-attribution-MLST/";
data_folder <- paste0(base_dir,"data")
results_folder <- paste0(base_dir,"results")
plot_folder <- paste0(base_dir,"plots")
script_folder <- paste0(base_dir, "scripts")

#print("DEBUG: changed filter from starts_with to starts_with('LMO003') to reduce it to 42 loci")
loci_start_with <- "LMO"

source(paste0(script_folder,"/wgMLST_funs_update.R"))

#lm_dat <- readRDS(paste0(data_folder,"/isolates_original_plus_new_dec_1_2021.rds"))
lm_dat <- read.csv(paste0(data_folder,"/isolates_original_plus_new_dec_1_2021.csv.gz"),header = TRUE)
cgmlst_loci <- read.csv(paste0(data_folder, "/cgMLST_loci.csv")) %>% names

### ht defines the threshold of the proportional difference within which isolates were treated
#### as originated from the same outbreaks or the collection from the same facilities

ht <- 0.004

si <- sel_rep_iso(lm_dat, ht) ### select representative isolates


###################################################################
#  importance of genes from random forest model based on all genes
###################################################################
train.df.all <- lm_dat %>%
  filter(SRR_ID %in% si) %>%
  #select("food", starts_with(loci_start_with)) %>%
  select("food", starts_with(loci_start_with)) %>%
  mutate_if(is.integer, coalesce, 0L) %>% # integer "LMOxxxxx" as integer (52.59%) performs similar to "LMOxxxxx" as factor (51.85%)
  mutate(across(starts_with(loci_start_with), ~ as.factor(as.character(.x)))) %>%
  as.data.frame() # rfsrc() doesn't work with a tibble

set.seed(23)
rf.all <- rfsrc(food ~ ., train.df.all, importance = T) #

rf <- rf.all
c.imp <- rf$importance[, 1]
rk.genes <- names(sort(c.imp, decreasing = T)) # [ZC: top ranked genes]

train.df.core_top100 <- lm_dat %>%
  filter(SRR_ID %in% si) %>%
  select("food", all_of(c(cgmlst_loci, rk.genes[1:100]))) %>%
  mutate_if(is.integer, coalesce, 0L) %>%
  mutate(across(starts_with(loci_start_with), ~ as.factor(as.character(.x)))) %>%
  as.data.frame()


rf.core_top100 <- rfsrc(food ~ ., train.df.core_top100, importance = F)
rf <- rf.core_top100
train.df <- train.df.core_top100

set.seed(23)
rf.core_top100 <- rfsrc(food ~ ., train.df, importance = F)

###### summary of PulseNet isolates and training dataset #############
table(lm_dat$food)
table(train.df$food) ### note change
################################################################

rf <- rf.all
train.df <- train.df.all
c.imp <- rf$importance[, 1]
p.imp <- names(sort(c.imp, decreasing = T))[1:100]
d.imp <- train.df[, names(train.df) %in% c("food", p.imp)]

############################################################################################
################# Calculation of AMOVA model ###############################################
############################################################################################

temp <- lm_dat %>%
  filter(SRR_ID %in% si) %>%
  select("food", starts_with(loci_start_with)) %>%
  mutate_if(is.integer, coalesce, 0L) %>%
  mutate(across(starts_with(loci_start_with), ~ as.factor(as.character(.x)))) %>%
  as.data.frame()

wgmlst <- temp %>% select(starts_with(loci_start_with))
fd <- levels(temp$food)
n.f <- length(table(temp$food))

fd.gene <- as.data.frame(cbind(temp$food, wgmlst))
names(fd.gene)[1] <- "food"

#### Note calculation amo is time consuming .... ######################

amo <- lapply(1:(n.f - 1), function(i) {
  lapply((i + 1):n.f, function(j) {
    print(paste(fd[i], fd[j], sep = "-"))
    dat <- fd.gene[fd.gene$food %in% c(fd[i], fd[j]), ]
    strata.df <- data.frame(food = dat$food)
    df.g <- df2genind(dat[, !names(dat) == "food"], ploidy = 1)
    strata(df.g) <- strata.df
    g <-suppressMessages(poppr.amova(df.g, ~food, nperm = 499))
    p <- randtest(g, nrepet = 499)
    list(fd_pair = paste(fd[i], fd[j], sep = "-"), amova = g, pvalue = p$pvalue)
  })
})

# save(amo, file = "amo.RData")
# load("amo.RData")

phi <- lapply(amo, function(p) unlist(lapply(p, function(g) round(g$amova$statphi * 100, 1))))
pvalue <- lapply(amo, function(x) unlist(lapply(x, function(g) g$pvalue)))
mt <- data.frame(matrix(rep(NA, length(fd)^2), ncol = length(fd)))
lapply(1:length(phi), function(p) {
  mt[p, (length(fd) - length(phi[[p]]) + 1):length(fd)] <<- paste(paste(phi[[p]], pvalue[[p]], sep = "("), ")", sep = "")
})
names(mt) <- fd
row.names(mt) <- fd

mt %>% gt()

###########################################
#   Table 2 Summary of predictive metrics
###########################################

cond <- data.frame(rf$predicted.oob)
pd <- c("_pred_rf", "_pred_rf", "_pred_rf", "_pred_rf", "_pred_rf")
names(cond) <- paste(names(cond), pd, sep = "")
res <- data.frame(to.matrix(train.df$food, levels(train.df$food)))
# res=data.frame(rf@responses@transformations)
names(res) <- paste(names(res), "_true", sep = "")
dat <- cbind(res, cond)

m <- multi_roc(dat)

####### estimate prediction uncertainty using multinomial distributio of conditional probability##

temp <- cond
temp1 <- data.frame(obs = 1:nrow(train.df), pred.fd = rf$class.oob, temp)

temp2 <- temp1 %>% gather("food", "prob", -c(1:2))
# temp2$food=sapply(temp2$food,function(x) strsplit(x,split = '\\.')[[1]][2])
temp2 <- data.table(temp2[order(temp2$obs), ])
conf <- cal_conf_attr(n.rep = 500, temp2, train = TRUE)
conf$attribution[, 2:5] <- round(conf$attribution[, 2:5], 2)
out <- conf$attribution
out$report <- paste(paste(paste(out$mean, out$lower, sep = "("), out$upper, sep = "-"), ")", sep = "")
out$food <- levels(train.df$food)
out <- out[order(out$food), ]

se <- sen_spe(rf$class.oob, train.df$food)$se
sp <- sen_spe(rf$class.oob, train.df$food)$sp

dd <- data.frame(food = names(se), se = round(se, 2), sp = round(sp, 2), AUC = round(as.numeric(unlist(m$AUC)), 2)[1:length(se)])
dd <- dd[-dim(dd)[1], ]

##### Table 2 #######################
dd %>% gt()


##########################################################################################################
#      Figure 1 prediction accuracy based on different missing value thresholds and cutoffs for cluster
##########################################################################################################

######## sensitivity analysis of cluster definition and handling of loci with missing alleles #########################
# based on the prediction accuracy of random forest models built from the training datasets created
# under different cluster criteria (ht) and proportions of missing threshold

### the algorithm uses parallel computering based on multi-core CUP and it may take long time for the systems with a small number of cores 

ht <- c(0, 0.001, 0.004, 0.01) ## row selection, h is equivalent to ht above
missing_threshold <- c(0, 0.01, 0.05, 0.1, 0.3, 1) #### column selection, loci would be excluded if missing prop is greater than the threshold value

hh <- c(sapply(ht,function(x) lapply(missing_threshold,function(y) c(x,y))))

ncores <- ifelse(detectCores() - 1 < 1, 1, detectCores() - 1)
cl <- makeCluster(ncores)
registerDoParallel(cl) # number of cores on the machine

### Note large numbers of trees may take longer time to complete ########################
n_tree <- 300
res.se <- foreach(n = 1:length(hh), .packages = c("dplyr", "cluster", "randomForestSRC")) %dopar% {
  if (hh[[n]][1] == 0) ssi <- as.character(lm_dat$SRR_ID) else ssi <- sel_rep_iso(lm_dat, hh[[n]][1]) ### select representative isolates at different levels of ht
  dat <- lm_dat %>%
    filter(SRR_ID %in% ssi) %>%
    select("food", starts_with(loci_start_with))
  dat$food <- factor(dat$food) # identical to train.df when ht<-0.004 and miss_thres<-0.01
  dat <- dat[, apply(dat, 2, function(x) sum(is.na(x)) / length(x)) <= hh[[n]][2]] ### excluding loci with proportion of missingness > miss_thres
  dat <- dat %>%
    mutate_if(is.integer, coalesce, 0L) %>%
    mutate(across(starts_with(loci_start_with), ~ as.factor(as.character(.x)))) %>%
    as.data.frame()
  set.seed(78)
  m <- rfsrc(food ~ ., dat, ntree = n_tree, importance = TRUE) # impute missing values and build random forest model
  tab <- table(m$yvar, m$class.oob) # m$yvar contains values for food in dat, wheares m$class.oob contains values for food predicted by the random forest model for the purpose of caluclating OOB
  accu <- c(diag(tab) / table(m$yvar), sum(diag(tab)) / length(m$yvar))

  return(list(model = m, matrix = tab, accu = accu, rep = hh[[n]][1], dsize = length(ssi), dat = dat, miss = hh[[n]][2]))
}
stopCluster(cl)

# save(res.se, file = "sensitivity.RData")
# load("sensitivity.RData")

df <- data.table(do.call(rbind, lapply(res.se, function(x) x$matrix))) # combine the 24 m$yvar - m$class.oob to a df
df$row_sele <- unlist(lapply(res.se, function(x) rep(x$rep, length(unique(lm_dat$food))))) # h
df$column_sele <- unlist(lapply(res.se, function(x) rep(x$miss, length(unique(lm_dat$food))))) # m
df$n.genes <- unlist(lapply(res.se, function(x) rep(dim(x$dat)[2] - 1, length(unique(lm_dat$food))))) # -1 is to exclude "food" from counting
df$n.iso <- unlist(lapply(res.se, function(x) rep(dim(x$dat)[1], length(unique(lm_dat$food)))))
df$n.iso1 <- unlist(lapply(res.se, function(x) rep((x$dsize), length(unique(lm_dat$food))))) # exactly the same as df$n.iso
df$accuracy <- 100 * unlist(lapply(res.se, function(x) x$accu[-length(x$accu)])) # accuracy for each food under each h and each m

df.lg <- (gather(df, index, value, 1:length(unique(lm_dat$food)))) # transform df to long form
df.lg$food <- rep(levels(lm_dat$food), nrow(df))
df.lg$row_sele <- factor(df.lg$row_sele)
df.lg$column_sele <- factor(df.lg$column_sele)

plot.col <- c("purple", "orange", "red", "blue", "green")
###### modify panel strip text
rn <- df %>%
  select(row_sele, column_sele, n.genes) %>%
  distinct() %>%
  group_by(column_sele, n.genes) %>%
  summarise(n = n()) %>%
  filter(n >= 2)

rn <- sapply(rn$n.genes, function(x) {
  if (rn %>% filter(n.genes == x) %>% `$`(n) %>% `==`(4)) {
    paste0("(", x, ")*")
  } else {
    paste0("(", x, ")**")
  }
})

df %>%
  select(row_sele, column_sele, n.genes, n.iso) %>%
  distinct() %>%
  arrange(column_sele, row_sele) -> df1

lapply(unique(df1$column_sele), function(pct_missing) {
  df1 %>%
    filter(column_sele == pct_missing) %>%
    `$`(n.genes) %>%
    paste0(collapse = ", ") %>%
    paste0("(No. of genes: ", ., ")")
}) -> n_genes_by_pctNA

strip.txt <- list(
  "0" = "without missing alleles", "0.01" = "missingness <= 1%", "0.05" = "missingness <= 5%",
  "0.1" = "missingness <= 10%", "0.3" = "missingness <= 30%", "1" = "all wgMLST genes included"
)
strip.txt <- lapply(1:length(strip.txt), function(x) paste(strip.txt[[x]], n_genes_by_pctNA[x], sep = "\n"))
labeller.fun <- function(variable, value) {
  return(strip.txt[value])
}

xt <- rev(names(table(df.lg$n.iso)))

x.txt <- paste(paste(paste(as.numeric(as.character(unique(df.lg$row_sele))) * 100, "%", sep = ""), xt, sep = "\n(n="), ")", sep = "")

lt1 <- c("dashed", "dotted", "dotdash", "longdash", "twodash", "F1", "4C88C488", "12345678")
lt <- lt1[1:length(levels(lm_dat$food))]

#########  Figure 1 ####################################################################
ggplot(df.lg, aes(row_sele, accuracy, group = food, color = food, linetype = food)) +
  geom_line(size = 1) +
  geom_point(aes(shape = food), size = 2) +
  scale_linetype_manual(values = lt) +
  facet_wrap(~column_sele, labeller = labeller.fun) +
  scale_color_manual(values = plot.col) +
  scale_shape_manual(values = c(3, 6, 15, 9, 17)) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    labels = as.character(seq(0, 100, 20))
  ) +
  scale_x_discrete(labels = x.txt) +
  xlab("cutoff values of proportional pairwise difference\n (number of isolates selected)") +
  ylab("prediction accuracy (%)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

#ggsave(filename = "fig1_sensitivity2.png", width = 15, height = 15, units = "in")


####### Figure 2 predictability over ranked loci #########################################
train.imp <- train.df.all
my.seq <- c(dim(train.imp)[2] - 1, seq(1600, 1000, -300), seq(900, 100, -200), seq(90, 50, -20), seq(45, 10, -5))
my.seq <- c(my.seq, c(8, 5, 3))

rk.genes <- names(sort(c.imp, decreasing = T)) # [ZC: top ranked genes]

cl <- makeCluster(ncores)
registerDoParallel(cl)
rfsrc.seq <- seq.rfsrc(train.imp, rk.genes, my.seq, ntree = 500)
stopCluster(cl)

# save(rfsrc.seq, file = "rfsrc_seq.RData")
# load("rfsrc_seq.RData")

dat <- rfsrc.seq$res
ltype <- c("dotdash", "dotted", "solid", "dashed", "twodash", "longdash")
my.col <- c("red", "orange", "black", "purple", "green", "blue")

dat$food <- factor(dat$food, levels = levels(reorder(dat$food, -dat$se)))
ltype <- ifelse(levels(as.factor(dat$food)) == "overall", "solid", "longdash")

########### Figure 2 ########################################################

ggplot(dat, aes(loci, se, group = food, color = food)) +
  geom_line(aes(linetype = food), size = 1) +
  ylab("estimated sensitivity") +
  geom_point(aes(shape = food), size = 2) +
  coord_trans(x = "log") +
  scale_x_continuous(breaks = c(3, 5, 10, 20, 30, 50, 70, 100, 200, 500, 1000, 2000, 4804)) +
  xlab("number of the ranked wgMLST loci") +
  scale_linetype_manual(values = ltype) +
  scale_shape_manual(values = c(15, 6, 16, 3, 17, 9)) +
  scale_color_manual(values = my.col) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = sprintf("%4.2f", seq(0, 1, 0.2))) +
  theme(
    legend.position = "right",
    text = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid.minor.x = element_blank()
  ) +
  guides(
    fill = guide_legend(keywidth = 1, keyheight = 1),
    linetype = guide_legend(keywidth = 3, keyheight = 1),
    colour = guide_legend(keywidth = 3, keyheight = 1)
  )

#ggsave(filename = "fig2_pred_accu_loci.png", width = 15, height = 15, units = "in")


############# Figure 3 panel plot of conditional probability of isolates ########################
# [the colors of meat and vegetables are inconsistent with those on the manuscript]
rf <- rf.core_top100
train.df <- train.df.core_top100
pred <- data.frame(rf$predicted.oob) # [ZC: predicted probability for each food category]
names(pred) <- levels(train.df$food)
pred$pred.fd <- rf$class.oob # [ZC: food predicted with the largest probability]
pred$obs <- row.names(train.df) # [ZC: row numbers in train.df]
pred$obs.fd <- train.df$food # [ZC: observed food category]
pred.l <- data.table(gather(pred, food, prob, 1:(dim(pred)[2] - 3)))
pp <- plot_panel_pred_prob_ind(pred.l, train = TRUE) ### train indicate train data and correct prediction is rugged

###### Figure 3 ###################################################################
pp$p
#ggsave(filename = "fig3_pred_prob2.png", width = 15, height = 15, units = "in")

