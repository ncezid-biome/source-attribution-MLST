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

lm_dat <- readRDS("isolates_original_plus_new_dec_1_2021.rds")
### cgMLST loci



### ht defines the threshold of the proportional difference within which isolates were treated
#### as originated from the same outbreaks or the collection from the same facilities
ht <- 0.004
si <- sel_rep_iso(lm_dat, ht) ### select representative isolates

train.df.all <- lm_dat %>%
  select("food", starts_with("LMO")) %>%

rf.all <- rfsrc(food ~ ., train.df.all, importance = T) #
rf <- rf.all
c.imp <- rf$importance[, 1]
rk.genes <- names(sort(c.imp, decreasing = T)) # [ZC: top ranked genes]
  filter(SRR_ID %in% si) %>%
  select("food", all_of(c(cgmlst_loci, rk.genes[1:100]))) %>%
  mutate_if(is.integer, coalesce, 0L) %>%
  mutate(across(starts_with("LMO"), ~ as.factor(as.character(.x)))) %>%
  as.data.frame()

train.df <- train.df.core_top100


p.imp <- names(sort(c.imp, decreasing = T))[1:100]

############################################################################################
temp <- lm_dat %>%
  as.data.frame()
wgmlst <- temp %>% select(starts_with("LMO"))
fd <- levels(temp$food)
n.f <- length(table(temp$food))

#### Note calculation amo is time consuming .... ######################
    strata(df.g) <- strata.df
    g <-suppressMessages(poppr.amova(df.g, ~food, nperm = 499))
    p <- randtest(g, nrepet = 499)
    list(fd_pair = paste(fd[i], fd[j], sep = "-"), amova = g, pvalue = p$pvalue)
})


phi <- lapply(amo, function(p) unlist(lapply(p, function(g) round(g$amova$statphi * 100, 1))))
pvalue <- lapply(amo, function(x) unlist(lapply(x, function(g) g$pvalue)))
mt <- data.frame(matrix(rep(NA, length(fd)^2), ncol = length(fd)))
lapply(1:length(phi), function(p) {
  mt[p, (length(fd) - length(phi[[p]]) + 1):length(fd)] <<- paste(paste(phi[[p]], pvalue[[p]], sep = "("), ")", sep = "")
names(mt) <- fd
row.names(mt) <- fd


cond <- data.frame(rf$predicted.oob)
pd <- c("_pred_rf", "_pred_rf", "_pred_rf", "_pred_rf", "_pred_rf")
names(cond) <- paste(names(cond), pd, sep = "")
res <- data.frame(to.matrix(train.df$food, levels(train.df$food)))
names(res) <- paste(names(res), "_true", sep = "")



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


##### Table 2 #######################
dd %>% gt()

##########################################################################################################

# under different cluster criteria (ht) and proportions of missing threshold
### the algorithm uses parallel computering based on multi-core CUP and it may take long time for the systems with a small number of cores 
missing_threshold <- c(0, 0.01, 0.05, 0.1, 0.3, 1) #### column selection, loci would be excluded if missing prop is greater than the threshold value

cl <- makeCluster(ncores)
registerDoParallel(cl) # number of cores on the machine
n_tree <- 300
oportion of missingness > miss_thres
  dat <- dat %>%
    mutate_if(is.integer, coalesce, 0L) %>%
    mutate(across(starts_with("LMO"), ~ as.factor(as.character(.x)))) %>%
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


