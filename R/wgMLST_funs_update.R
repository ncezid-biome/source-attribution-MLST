#' @title sel_rep_iso
#'
#' @param dat 
#' @param ht 
#' @param cgmlst_loci 
#'
#' @return cluster.id
#' @export
#'
#' @examples TODO
sel_rep_iso <- function(dat, ht, cgmlst_loci) {
  genes <- dat[, names(dat) %in% cgmlst_loci] ###### select core genes ######
  genes <- data.frame(lapply(genes, function(x) factor((x))))
  epi <- dat %>% select(-starts_with("LMO")) # non-gene variables

  (cl <- daisy(genes, metric = "gower")) # works well under R-3.6.1 but gives fetal errors under R-4.1.2
  hc <- hclust(cl, method = "complete")
  hc$labels <- epi$SRR_ID

  ct <- cutree(hc, h = ht) # assign cluster number (1 to 320) to each isolate (SSR_ID)

  ff <- as.character(dat$SourceSite)
  ifsac <- as.character(epi$food)
  ct.mem <- lapply(unique(ct), function(x) {
    w <- which(match(ct, x) == 1) # position of isolates (SSR_ID) in a cluster
    nam <- names(ct)[w] # SSR_ID of isolates in a cluster
    ss <- ff[!is.na(match(dat$SRR_ID, nam))] # SourceSites of isolates in a cluter
    ic <- ifsac[match(nam, epi$SRR_ID)] # IFSAC food categories of isolates in a cluter
    if (length(nam) > 1) { ### if there are more than one SSR_ID in a cluster
      sel.nam <- names(table(ic)) # unique IFSAC food categories in a cluster
      if (length(sel.nam) > 1) { # if there is more than one IFSAC food category in a cluster, one SSR_ID was randomly selected from each food category
        id <- sapply(sel.nam, function(y) {
          y.id <- nam[ic == y]
          set.seed(34)
          sample(y.id, 1)
        })
      } else {
        set.seed(23) # if there is only one food category in a cluster, one SSR_ID was randomly selected from that cluster
        id <- sample(nam[which(ic %in% sel.nam)], 1)
      }
    } else {
      id <- nam
    } ### if there is only one SSR_ID in a cluster, that SSR_ID is selected
    list(mem.id = nam, sourcesite = ss, ifsac = ic, sel.id = id)
  })
  sapply(ct.mem, function(x) length(x$mem.id)) # number of isolates in each cluster

  cluster.id <- unlist(sapply(ct.mem, function(x) x$sel.id)) # SSR_ID of randomly selected isolates

  return(cluster.id)
}


#' @title cal_conf_attr 
#'
#' @param n.rep 
#' @param pred 
#' @param train 
#'
#' @return I have no idea
#' @import data.table
#' @export
#'
#' @examples TODO
cal_conf_attr <- function(n.rep, pred, train = NULL) {
  pred <- data.table(pred)
  #### calculate intervals of uncertainty of predicted class ####################################
  # oob.accu=cf.seq$res[which.max(cf.seq$res$overall)] ### adjust estimation by oob accuracy of training data
  # oob.accu=unlist(oob.accu)[match(levels(train.df$food),names(oob.accu))]
  oob.accu <- rep(1, length(unique(train.df$food)))

  pred <- pred[order(as.numeric(pred$obs))]
  temp <- pred[, list(
    perm.ind = c(rmultinom(n.rep, 1, prob = prob / sum(prob))),
    id = rep(c(1:n.rep), each = length(unique(train.df$food)))
  ), obs]
  temp[, perm.class := ifelse(perm.ind == 1, unique(pred$food)[which(perm.ind == 1)], NA), by = c("obs", "id")]

  ####### calculate accuracy uncertainty of bootstrap samples ############
  if (!is.null(train)) {
    g <- temp[, .(obs = unique(obs), pred = perm.class[!is.na(perm.class)]), id]
    sp <- split(g, g$id)
    gg <- lapply(sp, function(x) {
      s <- sen_spe(x$pred, train.df$food)
      # print(diag(table(x$pred,train.df$food))/table(train.df$food))
      list(accu = c(s$accurate_rate, s$over_accu), se = s$se, sp = s$sp, kappa = s$kappa)
    })
    my.f <- function(x) {
      c(mean(x), quantile(x, probs = c(.025, .975)))
    }
    accu <- apply(do.call(rbind, lapply(gg, function(x) x$accu)), 2, function(y) my.f(y))
    se <- apply(do.call(rbind, lapply(gg, function(x) x$se)), 2, function(y) my.f(y))
    sp <- apply(do.call(rbind, lapply(gg, function(x) x$sp)), 2, function(y) my.f(y))
    k <- apply(do.call(rbind, lapply(gg, function(x) x$kappa)), 2, function(y) my.f(y))
  }
  #######################################################################################################

  # temp[,agree:=ifelse(perm.class==train.df$food[obs],1,0)]
  temp1 <- temp[, .(attr = c(table(perm.class) / sum(table(perm.class))), food = unique(pred$food)), id]
  # temp1[,quantile(attr,probs=c(.025,.5,.975)),food]
  temp2 <- temp1[, .(attr = c(mean(attr), sd = sd(attr), quantile(attr, probs = c(.025, .975))), nam = c("mean", "sd", "lower", "upper")), food]
  # over=temp[,.(over=c(sum(agree,na.rm=T)/nrow(train.df))),id]

  temp2 <- temp2 %>% spread(key = nam, value = attr)
  temp2$food <- factor(temp2$food, levels = levels(reorder(temp2$food, -temp2$mean)))
  if (!is.null(train)) {
    (list(attribution = temp2, conf95_accuracy = accu, conf95_se = se, conf95_sp = sp, conf95_kappa = k))
  } else {
    (list(attribution = temp2))
  }
}


#' @title plot_panel_pred_prob_ind
#'
#' @param pred 
#' @param train 
#'
#' @return list(p=p)
#' @export
#'
#' @examples TODO
plot_panel_pred_prob_ind <- function(pred, train = NULL) {
  sp <- split(pred, pred$pred.fd)
  g <- lapply(1:length(sp), function(d) { #### reorder obs based on pred prob of THE food
    prob <- unlist(tapply(sp[[d]]$prob, sp[[d]]$obs, function(x) rep(x[d], length(x)))) # [ZC: replicate the largest probability for a given isolate X times (X is the number of food categories).
    #     [As a result, in the data set for pred.fd == "fruit", the prob for fruit was replicated X times
    dat <- sp[[d]][order(prob, decreasing = T), ]
  })
  pred1 <- do.call(rbind, g)
  pred1$obs <- factor(pred1$obs, levels = unique(pred1$obs)) #### keep the correct order

  pred1$color <- ifelse(pred1$pred.fd == "dairy", "", ifelse(pred1$pred.fd == "fruit", "grey",
    ifelse(pred1$pred.fd == "meat", "yellow", ifelse(pred1$pred.fd == "sea food", "green", "purple"))
  ))

  pred1$color.obs <- ifelse(pred1$obs.fd == "dairy", "red", ifelse(pred1$obs.fd == "fruit", "grey",
    ifelse(pred1$obs.fd == "meat", "yellow", ifelse(pred1$obs.fd == "sea food", "green", "purple"))
  ))

  pred1$food <- ifelse(pred1$food == "sea food", "seafood", ifelse(pred1$food == "veg", "vegetable", as.character(pred1$food)))

  pred1$pred.fd <- ifelse(pred1$pred.fd == "sea food", "seafood", ifelse(pred1$pred.fd == "veg", "vegetable", as.character(pred1$pred.fd)))

  my.col <- c("purple", "orange", "red", "blue", "green")
  p <- ggplot(pred1, aes(obs, prob)) +
    geom_bar(aes(fill = food), stat = "identity", width = 1) +
    ylab("predicted probability of foods") +
    scale_fill_manual(values = my.col) +
    xlab("isolates") +
    theme(
      axis.line.x = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    theme(
      text = element_text(size = 16), axis.line.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
      axis.ticks.x = element_blank(), legend.position = c(0.8, 0.3), legend.text = element_text(size = 14),
      legend.title = element_text(size = 14), strip.text.x = element_text(size = 15)
    ) +
    facet_wrap(~pred.fd, scales = "free_x")

  if (is.null(train)) {
    p <- p + geom_rug(aes(obs), sides = "b")
  } else {
    p <- p + geom_rug(aes(obs), subset(pred1, obs.fd == pred.fd), sides = "b")
  }
  return(list(p = p))
}


#' @title seq.rfsrc
#'
#' @param dat 
#' @param imp 
#' @param my.seq 
#' @param ntree 
#'
#' @return list(res = dat)
#' @import foreach
#' @export
#'
#' @examples TODO
seq.rfsrc <- function(dat, imp, my.seq, ntree) {
  res <- foreach(n = my.seq, .packages = "randomForestSRC") %dopar% {
    set.seed(23)
    m <- rfsrc(food ~ ., dat[, names(dat) %in% c("food", as.character(imp[1:n]))], ntree = 1000)
    tab <- table(m$yvar, m$class.oob)
    t.pos <- diag(tab)
    f.neg <- sapply(1:length(unique(dat$food)), function(x) sum(tab[x, ]) - tab[x, x])
    sen <- c(round(t.pos / (t.pos + f.neg), 2))
    over_sen <- round(sum(t.pos) / sum(t.pos + f.neg), 2)
    list(food = c(dimnames(tab)[[1]], "overall"), se = c(sen, over_sen))
  }

  dat <- data.table(
    food = c(sapply(res, function(x) unlist(x$food))),
    se = c(sapply(res, function(x) unlist(x$se))), loci = rep(my.seq, each = length(res[[1]]$food))
  )

  return(list(res = dat))
}

#' @title sen_spe
#'
#' @param pred 
#' @param obs 
#'
#' @return large list
#' @export
#'
#' @examples TODO
sen_spe <- function(pred, obs) {
  # pred<-predict(tree1,newdata=comp.ue,type='class')
  tab <- xtabs(~ obs + pred)
  (correct.rate <- sum(diag(tab)) / sum(tab))
  mar1 <- apply(tab, 1, sum) / sum(tab)
  mar2 <- apply(tab, 2, sum) / sum(tab)
  exp.agg <- sum(mar1 * mar2)
  kappa.est <- (correct.rate - exp.agg) / (1 - exp.agg)
  kappa.est

  (t.pos <- diag(tab))
  (f.neg <- sapply(1:length(levels(obs)), function(x) sum(tab[x, ]) - tab[x, x]))
  (sen <- round(t.pos / (t.pos + f.neg), 2))

  (f.pos <- sapply(1:length(levels(obs)), function(x) sum(t(tab)[x, ]) - t(tab)[x, x]))
  (t.neg <- sapply(1:length(levels(obs)), function(x) sum(tab) - sum(tab[x, ])))
  (spe <- round(t.neg / (t.neg + f.pos), 2))
  o.spe <- sum(t.neg) / sum(t.neg + f.pos)

  over_accu <- (sum(t.pos) + sum(t.neg)) / (sum(t.pos) + sum(t.neg) + sum(f.pos) + sum(f.neg))
  prev <- table(train.df$food) / nrow(train.df)
  dat <- data.frame(tab, apply(tab, 1, sum))
  return(list(dat = dat, se = c(sen, correct.rate), sp = c(spe, o.spe), accurate_rate = sen * prev + spe * (1 - prev), over_accu = over_accu, kappa = kappa.est))
}

