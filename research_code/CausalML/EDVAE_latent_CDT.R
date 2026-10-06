data=read.csv("C:/Users/ASUS/Desktop/cg working/Variational autoencoder/vae_compare/data4000.csv")
data=data[,-1]
names(data)
X=data[,c(5:19)]
Xy=data[,c(1,3,5:19)]
edvae_latent=read.csv("C:/Users/ASUS/Desktop/cg working/Variational autoencoder/vae_compare/edvae_latentout.csv")
edvae_latent=edvae_latent[,-1]
names(edvae_latent)
dim(edvae_latent)

#graphics.off()
par(mfrow=c(1,1))
cor(edvae_latent)
heatmap(cor(edvae_latent))

cor(X,edvae_latent)
heatmap(cor(X,edvae_latent))
cor(Xy,edvae_latent)
heatmap(cor(Xy,edvae_latent))

edvae_latent1=edvae_latent[data$Treat==1,]
edvae_latent0=edvae_latent[data$Treat==0,]
X1=X[data$Treat==1,]
X0=X[data$Treat==0,]
cor(X1,edvae_latent1)
heatmap(cor(X1,edvae_latent1))
cor(X0,edvae_latent0)
heatmap(cor(X0,edvae_latent0))


library(ggplot2)
library(reshape2)

latentcor=melt(cor(edvae_latent))
Xlatentcor=melt(cor(X,edvae_latent))
Xlatentcory=melt(cor(Xy,edvae_latent))
cor_over=cor(Xy,edvae_latent)
cor_over[abs(cor(Xy,edvae_latent))<0.1]=0
Xlatentcory_over=melt(cor_over)
ggp1=ggplot(latentcor,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  labs(x="latent variable",y="latent variavle")+
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint = 0,limits=c(-1,1))
ggp1

ggp2=ggplot(Xlatentcor,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  labs(x="covariates",y="latent variavle")+
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint = 0,limits=c(-1,1))
ggp2

ggp2=ggplot(Xlatentcory,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  labs(x="covariates",y="latent variavle")+
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint = 0)
ggp2

Xlatentcor_abs=Xlatentcory
Xlatentcor_abs$value=abs(Xlatentcor_abs$value)
ggp2.5=ggplot(Xlatentcor_abs,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  labs(x="covariates",y="latent variavle")+
  scale_fill_gradient(low="white",high="red",limits=c(0,1))
ggp2.5

ggp2.6=ggplot(Xlatentcory_over,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  labs(x="covariates",y="latent variavle")+
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint = 0)
ggp2.6



library(energy)
library(minerva)
#dCov, dCor
datadcov=data[,c(1,3,5:19)]
dcov_matrix=matrix(0,dim(edvae_latent)[2],dim(datadcov)[2])
dcor_matrix=matrix(0,dim(edvae_latent)[2],dim(datadcov)[2])
for(i in c(1:dim(datadcov)[2])){
  for(j in c(1:dim(edvae_latent)[2])){
    a=dcov.test(datadcov[,i],edvae_latent[,j])
    dcov_matrix[j,i]=a$estimates[1]
    dcor_matrix[j,i]=a$estimates[2]
    print(c(j,i,a$estimates[1],a$estimates[2]))
  }
}
colnames(dcov_matrix)=c("outcome","treatment","x1","x2","x3","x4","x5","x6","x7","x8",
                        "x9","x10","x11","x12","x13","x14","x15")
row.names(dcov_matrix)=c("zt1","zt2","zc1","zc2","zc3","zy1","zy2","zy3","zy4","zy5")
#row.names(dcov_matrix)=c("zt1","zt2","zt3","zt4","zt5","zt6","zt7","zt8","zt9",
#                         "zc1","zc2","zc3","zc4","zc5","zc6","zy1","zy2")
colnames(dcor_matrix)=c("outcome","treatment","x1","x2","x3","x4","x5","x6","x7","x8",
                        "x9","x10","x11","x12","x13","x14","x15")
row.names(dcor_matrix)=c("zt1","zt2","zc1","zc2","zc3","zy1","zy2","zy3","zy4","zy5")
#row.names(dcor_matrix)=c("zt1","zt2","zt3","zt4","zt5","zt6","zt7","zt8","zt9",
#                         "zc1","zc2","zc3","zc4","zc5","zc6","zy1","zy2")
latentdcov=melt(dcov_matrix)
latentdcor=melt(dcor_matrix)
ggp_dcov=ggplot(latentdcov,aes(x=Var2,y=Var1,fill=value))+
  geom_tile()+
  labs(x="covariates",y="latent variavle")+
  scale_fill_gradient(low="white",high="red",limits=c(0,1.5))
ggp_dcov
ggp_dcor=ggplot(latentdcor,aes(x=Var2,y=Var1,fill=value))+
  geom_tile()+
  labs(x="covariates",y="latent variavle")+
  scale_fill_gradient(low="white",high="red",limits=c(0,1))
ggp_dcor

#MIC, TIC
dataMIC=data[,c(1,3,5:19)]
MIC_matrix=matrix(0,dim(edvae_latent)[2],dim(dataMIC)[2])
TIC_matrix=matrix(0,dim(edvae_latent)[2],dim(dataMIC)[2])
for(i in c(1:dim(dataMIC)[2])){
  for(j in c(1:dim(edvae_latent)[2])){
    a=mine(dataMIC[,i],edvae_latent[,j])
    MIC_matrix[j,i]=a$MIC
    TIC_matrix[j,i]=a$TIC
    print(c(j,i,a$MIC,a$TIC))
  }
}
colnames(MIC_matrix)=c("outcome","treatment","x1","x2","x3","x4","x5","x6","x7","x8",
                        "x9","x10","x11","x12","x13","x14","x15")
row.names(MIC_matrix)=c("zt1","zt2","zc1","zc2","zc3","zy1","zy2","zy3","zy4","zy5")
colnames(TIC_matrix)=c("outcome","treatment","x1","x2","x3","x4","x5","x6","x7","x8",
                        "x9","x10","x11","x12","x13","x14","x15")
row.names(TIC_matrix)=c("zt1","zt2","zc1","zc2","zc3","zy1","zy2","zy3","zy4","zy5")
latentMIC=melt(MIC_matrix)
latentTIC=melt(TIC_matrix)
ggp_MIC=ggplot(latentMIC,aes(x=Var2,y=Var1,fill=value))+
  geom_tile()+
  labs(x="covariates",y="latent variavle")+
  scale_fill_gradient(low="white",high="red")
ggp_MIC
ggp_TIC=ggplot(latentTIC,aes(x=Var2,y=Var1,fill=value))+
  geom_tile()+
  labs(x="covariates",y="latent variavle")+
  scale_fill_gradient(low="white",high="red")
ggp_TIC

print(cor(Xy,edvae_latent)[16,6])
print(dcor_matrix[6,16])
print(MIC_matrix[6,16])
print(TIC_matrix[6,16])

print(cor(Xy,edvae_latent)[3,6])
print(dcor_matrix[6,3])
print(MIC_matrix[6,3])
print(TIC_matrix[6,3])


#plot covariate-----------------------------------------------------------------
edvae_out=read.csv("C:/Users/ASUS/Desktop/cg working/Variational autoencoder/vae_compare/edvae_to_metalearner_output.csv")
edvae_out = edvae_out[,-1]
hist(edvae_out$propensity_hat,breaks=40)
plot(edvae_out$edvae_ite, edvae_out$ite)
abline(0,1,col="red")
library(MASS)
library(ggplot2)
X <- edvae_out[,4:18]
W <- edvae_out[,3]
Y <- edvae_out[,1]
e.hat <- edvae_out$propensity_hat
IPW <- ifelse(W == 1, 1 / e.hat, 1 / (1 - e.hat))
plot.df <- data.frame(X,
                      W = as.factor(W),
                      IPW = IPW)

ggplot(plot.df, aes(x = X5, fill = W)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30)
ggplot(plot.df, aes(x = X5, weight = IPW, fill = W)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30)
ggplot(plot.df, aes(x = X8, fill = W)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30)
ggplot(plot.df, aes(x = X8, weight = IPW, fill = W)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30)

#covariates balance
#X=covariates, W=treatment, phat=propensity score hat
love_table=function(X,W,phat){ 
  X_ad=X*matrix(W/phat+(1-W)/(1-phat),dim(X)[1],dim(X)[2])
  mycaculate=matrix(0,length(names(X)),2)
  for(i in c(1:length(names(X)))){
    mycaculate[i,1]=abs(mean(X[,i][W==1])-mean(X[,i][W==0]))/sqrt(var(X[,i][W==1])+var(X[,i][W==0]))
    mycaculate[i,2]=abs(sum(X_ad[,i][W==1])/sum(W/phat)-sum(X_ad[,i][W==0])/sum((1-W)/(1-phat)))/
      sqrt(cov.wt(as.matrix(X[,i]), (W/phat))$cov + cov.wt(as.matrix(X[,i]), ((1-W)/(1-phat)))$cov)
  }
  mycaculate_temp=mycaculate
  mycaculate_temp=matrix(c(t(mycaculate_temp)),dim(mycaculate_temp)[1]*2,1)
  mycaculate_temp=as.data.frame(mycaculate_temp)
  mycaculate_temp=cbind(c(t(matrix(names(X),length(names(X)),2))),c("Before adjustment","After adjustment"),mycaculate_temp)
  colnames(mycaculate_temp)=c("covariate_name","Cohort","value")
  od=order(mycaculate_temp$value[mycaculate_temp$Cohort=="Before adjustment"])
  c(t(matrix(c(2*od,2*od-1),length(od),2)))
  mycaculate_temp=mycaculate_temp[c(t(matrix(c(2*od-1,2*od),length(od),2))),]
  return(mycaculate_temp)
}
love_plot=function(X,lovetable){
  ggplot2::ggplot(lovetable,
                  ggplot2::aes(
                    x = value,
                    y = c(t(matrix(c(1:length(names(X))),length(names(X)),2))),
                    colour = Cohort
                  )
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(
      ggplot2::aes(group = Cohort),
      orientation = "y"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 1) +
    ggplot2::geom_vline(xintercept = 0.1, linetype = 2) +
    ggplot2::xlab("Absolute standardized mean difference") +
    ggplot2::ylab("")+
    scale_y_continuous(breaks=c(t(matrix(c(1:length(names(X))),length(names(X)),2))),
                       labels = lovetable$covariate_name)
}
lovetable=love_table(X,W,e.hat)
lovetable
loveplot=love_plot(X,lovetable)
loveplot

#-------------------------------------------------------------------------------


#CDT----------------------------------------------------------------------------
library(rpart)
library(partykit)
library(ggparty)
library(Rcpp)
library(glmnet)
library(stringr)
library(geepack)
library(grf)

#Distillation tree function
get_rpart_paths <- function(rpart_fit) {
  leaf_node_ids <- rpart_fit$frame |>
    tibble::rownames_to_column("id") |>
    dplyr::filter(var == "<leaf>") |>
    dplyr::pull("id") |>
    as.numeric()
  subgroups <- rpart::path.rpart(rpart_fit, leaf_node_ids, print.it = FALSE) |>
    purrr::map(~ setdiff(.x, "root")) |>
    purrr::compact()
  return(subgroups)
}
get_rpart_tree_info <- function(rpart_fit, digits = getOption("digits")) {
  out <- NULL
  splits <- rpart_fit$splits
  if (!is.null(splits) && isTRUE(nrow(splits) > 0)) {
    ff <- rpart_fit$frame
    is.leaf <- ff$var == "<leaf>"
    n <- nrow(splits)
    nn <- ff$ncompete + ff$nsurrogate + !is.leaf
    ix <- cumsum(c(1L, nn))
    ix_prim <- unlist(
      mapply(ix, ix + c(ff$ncompete, 0), FUN = seq, SIMPLIFY = F)
    )
    type <- rep.int("surrogate", n)
    type[ix_prim[ix_prim <= n]] <- "primary"
    type[ix[ix <= n]] <- "main"
    left <- character(nrow(splits))
    side <- splits[, 2L]
    for (i in seq_along(left)) {
      left[i] <- if (side[i] == -1L)
        paste("<", format(signif(splits[i, 4L], digits)))
      else if (side[i] == 1L)
        paste(">=", format(signif(splits[i, 4L], digits)))
      else {
        catside <- rpart_fit$csplit[splits[i, 4L], 1:side[i]]
        paste(c("L", "-", "R")[catside], collapse = "", sep = "")
      }
    }
    nodeids <- rep(as.integer(row.names(ff)), times = nn)
    out <- cbind(
      data.frame(
        var = rownames(splits),
        type = type,
        node = nodeids,
        ix = rep(seq_len(nrow(ff)), nn),
        depth = trunc(log(nodeids, base = 2)) + 1,
        left = left
      ),
      as.data.frame(splits, row.names = F)
    ) |>
      dplyr::filter(type == "main") |>
      dplyr::rename(thr = index) |>
      dplyr::mutate(
        cat_thr = purrr::pmap_chr(
          list(v = var, l = left),
          function(v, l) {
            if (grepl("[RL]", l)) {
              l_split <- strsplit(l, "")[[1]]
              R_idx <- which(l_split == "R")[1]
              L_idx <- which(l_split == "L")[1]
              if (R_idx < L_idx) {
                th <- levels(X[[v]])[L_idx]
              } else {
                th <- levels(X[[v]])[R_idx]
              }
            } else {
              return(NA)
            }
          }
        )
      )
  }
  return(out)
}
get_party_paths <- function(party_fit) {
  purrr::map(
    .list.rules.party(party_fit),
    ~ tibble::tibble(
      subgroup = list(stringr::str_split(.x, " & ")[[1]])
    )
  ) |>
    dplyr::bind_rows(.id = "leaf_id") |>
    dplyr::mutate(
      leaf_id = as.numeric(leaf_id)
    )
}
get_party_node_depths <- function(party_fit, return_features = FALSE) {
  printed_tree <- capture.output(party_fit)
  id_counter <- 1
  depths <- rep(NA, length(party_fit))
  names(depths) <- 1:length(party_fit)
  features <- rep(NA, length(party_fit))
  names(features) <- 1:length(party_fit)
  for (idx in seq_along(printed_tree)) {
    if (grepl("\\[[0-9]+\\]", printed_tree[idx])) {
      depths[id_counter] <- stringr::str_count(printed_tree[idx], "\\|")
      if (return_features) {
        features[id_counter] <- stringr::str_extract(
          printed_tree[idx],
          "(?<=\\]).*?(?=<|>)"
        ) |>
          stringr::str_trim()
      }
      id_counter <- id_counter + 1
    }
  }
  if (return_features) {
    return(
      tibble::tibble(
        depth = depths,
        feature = features
      )
    )
  } else {
    return(depths)
  }
}
.list.rules.party <- function(x, i = NULL, ...) {
  if (is.null(i)) {
    i <- partykit::nodeids(x, terminal = TRUE)
  }
  if (length(i) > 1) {
    ret <- sapply(i, .list.rules.party, x = x)
    names(ret) <- if (is.character(i)) i else names(x)[i]
    return(ret)
  }
  if (is.character(i) && !is.null(names(x))) {
    i <- which(names(x) %in% i)
  }
  stopifnot(length(i) == 1 & is.numeric(i))
  stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  dat <- partykit::data_party(x, i)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    fit <- dat[, findx:ncol(dat), drop = FALSE]
    dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0) {
      dat <- x$data
    }
  } else {
    fit <- NULL
    dat <- x$data
  }
  
  rule <- c()
  
  recFun <- function(node) {
    if (partykit::id_node(node) == i) return(NULL)
    kid <- sapply(partykit::kids_node(node), partykit::id_node)
    whichkid <- max(which(kid <= i))
    split <- partykit::split_node(node)
    ivar <- partykit::varid_split(split)
    svar <- names(dat)[ivar]
    index <- partykit::index_split(split)
    if (is.factor(dat[, svar])) {
      if (is.null(index)) {
        index <- ((1:nlevels(dat[, svar])) > breaks_split(split)) + 1
      }
      slevels <- levels(dat[, svar])[index == whichkid]
      srule <- paste(
        svar, " %in% c(\"",
        paste(slevels, collapse = "\", \"", sep = ""), "\")",
        sep = ""
      )
    } else {
      if (is.null(index)) index <- 1:length(kid)
      breaks <- cbind(
        c(-Inf, partykit::breaks_split(split)),
        c(partykit::breaks_split(split), Inf)
      )
      sbreak <- breaks[index == whichkid,]
      right <- partykit::right_split(split)
      srule <- c()
      if (is.finite(sbreak[1])) {
        srule <- c(srule, paste(svar, ifelse(right, ">", ">="), sbreak[1]))
      }
      if (is.finite(sbreak[2])) {
        srule <- c(srule, paste(svar, ifelse(right, "<=", "<"), sbreak[2]))
      }
      srule <- paste(srule, collapse = " & ")
    }
    rule <<- c(rule, srule)
    return(recFun(node[[whichkid]]))
  }
  node <- recFun(partykit::node_party(x))
  paste(rule, collapse = " & ")
}
student_rpart <- function(X, y, method = "anova", rpart_control = NULL,
                          prune = c("none", "min", "1se"), fit_only = FALSE) {
  y=c(as.matrix(y))
  prune <- match.arg(prune)
  df <- data.frame(X, y)
  
  # if tauhat is constant, return NULL model (no subgroups)
  if (length(unique(y)) == 1) {
    if (fit_only) {
      out <- NULL
    } else {
      out <- list(
        fit = NULL,
        tree_info = NULL,
        subgroups = list(),
        predictions = rep(unique(y), nrow(df))
      )
    }
  } else {
    if (is.null(rpart_control)) {
      fit <- rpart::rpart(
        y ~ ., data = df, method = method
      )
    } else {
      fit <- rpart::rpart(
        y ~ ., data = df, method = method, control = rpart_control
      )
    }
    
    # pruning
    if (prune != "none") {
      best_cp <- as.data.frame(fit$cptable) |>
        dplyr::filter(xerror == min(xerror, na.rm = TRUE)) |>
        dplyr::slice(1)
      if (prune == "min") {
        fit <- rpart::prune(fit, cp = best_cp$CP)
      } else if (prune == "1se") {
        best1se_cp <- as.data.frame(fit$cptable) |>
          dplyr::filter(xerror <= (best_cp$xerror + best_cp$xstd)) |>
          dplyr::filter(nsplit == min(nsplit, na.rm = TRUE)) |>
          dplyr::slice(1)
        fit <- rpart::prune(fit, cp = best1se_cp$CP)
      }
    }
    
    if (fit_only) {
      out <- fit
    } else {
      subgroups <- get_rpart_paths(fit)
      tree_info <- get_rpart_tree_info(fit)
      predictions <- predict(fit)
      out <- list(
        fit = fit,
        tree_info = tree_info,
        subgroups = subgroups,
        predictions = predictions
      )
    }
  }
  return(out)
}
student_model <- function(X,y,prune){
  student_rpart(X,y,method = "anova", rpart_control = rpart::rpart.control(minbucket = 30),
                prune, fit_only = FALSE)
}
estimate_group_cates_Y <- function(fit, X, Y, Z) {
  Z=c(as.matrix(Z))
  Y=c(as.matrix(Y))
  if (!is.null(fit)) {
    if ("rpart" %in% class(fit)) {
      fit <- partykit::as.party(fit)
    }
    leaf_ids <- tryCatch(
      predict(fit, data.frame(X), type = 'node'),
      error = function(e) as.numeric(as.factor(predict(fit, data.frame(X))))
    )
  } else {
    leaf_ids <- NULL
  }
  group_cates <- tibble::tibble(
    Z = Z,
    Y = Y,
    leaf_id = leaf_ids
  ) |>
    dplyr::group_by(dplyr::across(tidyselect::any_of("leaf_id"))) |>
    dplyr::summarise(
      estimate = mean(Y[Z == 1]) - mean(Y[Z == 0]),
      variance = 1 / sum(Z == 1) * var(Y[Z == 1]) +
        1 / sum(Z == 0) * var(Y[Z == 0]),
      .var1 = var(Y[Z == 1]),
      .var0 = var(Y[Z == 0]),
      .n1 = sum(Z == 1),
      .n0 = sum(Z == 0),
      .sample_idxs = list(dplyr::cur_group_rows()),
      .groups = "drop"
    )
  if ("party" %in% class(fit)) {
    group_cates <- dplyr::left_join(
      get_party_paths(fit),
      group_cates,
      by = "leaf_id"
    )
  }
  return(group_cates)
}
estimate_group_cates_pseudo_DR <- function(fit, X, Y, Z) {
  Z=c(as.matrix(Z))
  Y=c(as.matrix(Y))
  if (!is.null(fit)) {
    if ("rpart" %in% class(fit)) {
      fit <- partykit::as.party(fit)
    }
    leaf_ids <- tryCatch(
      predict(fit, data.frame(X), type = 'node'),
      error = function(e) as.numeric(as.factor(predict(fit, data.frame(X))))
    )
  } else {
    leaf_ids <- NULL
  }
  group_cates <- tibble::tibble(
    Z = Z,
    Y = Y,
    leaf_id = leaf_ids
  ) |>
    dplyr::group_by(dplyr::across(tidyselect::any_of("leaf_id"))) |>
    dplyr::summarise(
      estimate = mean(Y),
      .n1 = sum(Z == 1),
      .n0 = sum(Z == 0),
      .sample_idxs = list(dplyr::cur_group_rows()),
      .groups = "drop"
    )
  if ("party" %in% class(fit)) {
    group_cates <- dplyr::left_join(
      get_party_paths(fit),
      group_cates,
      by = "leaf_id"
    )
  }
  return(group_cates)
}
estimate_group_cates_pseudo_R <- function(fit, X, Y, Z, W) {
  Z=c(as.matrix(Z))
  Y=c(as.matrix(Y))
  W=c(as.matrix(W))
  if (!is.null(fit)) {
    if ("rpart" %in% class(fit)) {
      fit <- partykit::as.party(fit)
    }
    leaf_ids <- tryCatch(
      predict(fit, data.frame(X), type = 'node'),
      error = function(e) as.numeric(as.factor(predict(fit, data.frame(X))))
    )
  } else {
    leaf_ids <- NULL
  }
  group_cates <- tibble::tibble(
    Z = Z,
    Y = Y,
    W = W,
    leaf_id = leaf_ids
  ) |>
    dplyr::group_by(dplyr::across(tidyselect::any_of("leaf_id"))) |>
    dplyr::summarise(
      estimate = weighted.mean(Y,W),
      .n1 = sum(Z == 1),
      .n0 = sum(Z == 0),
      .sample_idxs = list(dplyr::cur_group_rows()),
      .groups = "drop"
    )
  if ("party" %in% class(fit)) {
    group_cates <- dplyr::left_join(
      get_party_paths(fit),
      group_cates,
      by = "leaf_id"
    )
  }
  return(group_cates)
}
sourceCpp("C:/Users/ASUS/Desktop/cg working/CausalML/causaltree/stability.cpp")
evaluate_subgroup_stability <- function(estimator, fit, X, y, Z = NULL,rpart_control = NULL,
                                        prune=prune,B = 100,max_depth = NULL) {
  y=c(as.matrix(y))
  if (!("rpart" %in% class(fit))) {
    warning(
      "fit is not an rpart object. ",
      "Stability diagnostics have only been implemented for the rpart student model. ",
      "Skipping stability diagnostics."
    )
    return(NULL)
  } else if ((B == 0) || is.null(fit)) {
    return(NULL)
  }
  
  fit_orig <- partykit::as.party(fit)
  node_depths_orig <- get_party_node_depths(fit_orig)
  leaf_ids_orig <- predict(fit_orig, data.frame(X), type = "node")
  if (is.null(max_depth)) {
    max_depth <- max(max(node_depths_orig), 4)
  }
  
  # modify rpart controls so that the tree is forced to make a split when possible
  rpart_control[["minsplit"]] <-  2
  rpart_control[["minbucket"]] <- 1
  rpart_control[["cp"]] <- 0
  rpart_control[["maxdepth"]] <- max_depth
  estimator <- purrr::partial(estimator,method = "anova",prune=prune,
                              rpart_control = rpart_control)
  
  bootstrap_out <- purrr::map(
    1:(2 * B),
    function(b) {
      bootstrap_idx <- sample(1:nrow(X), size = nrow(X), replace = TRUE)
      X_b <- X[bootstrap_idx,]
      y_b <- y[bootstrap_idx]
      if (is.null(Z)) {
        fit_b <- estimator(X = X_b, y = y_b, fit_only = TRUE)
      } else {
        Z_b <- Z[bootstrap_idx]
        fit_b <- estimator(X = X_b, Y = y_b, Z = Z_b)
      }
      if (!is.null(fit_b)) {
        fit_b <- partykit::as.party(fit_b)
        node_depths_b <- get_party_node_depths(fit_b)
        return(
          list(
            "fit" = fit_b,
            "node_depths" = node_depths_b
          )
        )
      } else {
        return(NULL)
      }
    }
  ) |>
    purrr::compact()
  
  bootstrap_fits <- purrr::map(bootstrap_out, "fit")
  node_depths <- purrr::map(bootstrap_out, "node_depths")
  
  Js <- list()
  preds_mean <- list()
  preds_var <- list()
  for (n_depth in 1:max_depth) {
    # if (any(node_depths_orig > n_depth)) {
    #   fit_orig_pruned <- partykit::nodeprune(
    #     fit_orig, ids = names(node_depths_orig)[node_depths_orig == n_depth]
    #   )
    # } else {
    #   fit_orig_pruned <- fit_orig
    # }
    # node_depths_orig_pruned <- get_party_node_depths(fit_orig_pruned)
    # leaf_ids_orig <- predict(fit_orig_pruned, data.frame(X), type = "node")
    
    bootstrap_leaf_ids <- purrr::map2(
      bootstrap_fits, node_depths,
      function(fit_b, node_depths_b) {
        if (any(node_depths_b > n_depth)) {
          fit_b_pruned <- partykit::nodeprune(
            fit_b, ids = names(node_depths_b)[node_depths_b == n_depth]
          )
        } else {
          fit_b_pruned <- fit_b
        }
        leaf_ids_b <- predict(fit_b_pruned, data.frame(X), type = "node")
        return(leaf_ids_b)
      }
    )
    
    bootstrap_leaf_preds <- purrr::map2(
      bootstrap_fits, node_depths,
      function(fit_b, node_depths_b) {
        if (any(node_depths_b > n_depth)) {
          fit_b_pruned <- partykit::nodeprune(
            fit_b, ids = names(node_depths_b)[node_depths_b == n_depth]
          )
        } else {
          fit_b_pruned <- fit_b
        }
        leaf_preds <- predict(fit_b_pruned, data.frame(X))
        return(leaf_preds)
      }
    )
    
    J <- purrr::map_dbl(
      1:floor(length(bootstrap_fits) / 2),
      ~ jaccardSSI(
        as.numeric(as.factor(bootstrap_leaf_ids[[.x * 2 - 1]])) - 1,
        as.numeric(as.factor(bootstrap_leaf_ids[[.x * 2]])) - 1
      )
    )
    Js[[n_depth]] <- J
    
    preds_mean[[n_depth]] <- do.call(cbind, bootstrap_leaf_preds) |>
      rowMeans()
    preds_var[[n_depth]] <- do.call(cbind, bootstrap_leaf_preds) |>
      apply(1, var)
  }
  
  feature_dist <- purrr::map(
    bootstrap_fits, ~ get_party_node_depths(.x, return_features = TRUE)
  ) |>
    dplyr::bind_rows(.id = "bootstrap_idx") |>
    dplyr::filter(
      !is.na(feature)
    ) |>
    dplyr::group_by(depth, feature) |>
    dplyr::summarise(
      freq = dplyr::n()
    ) |>
    dplyr::ungroup()
  
  out <- list(
    "jaccard_mean" = sapply(Js, mean),
    "jaccard_distribution" = Js,
    "feature_distribution" = feature_dist,
    "bootstrap_predictions_mean" = preds_mean,
    "bootstrap_predictions_var" = preds_var,
    "leaf_ids" = leaf_ids_orig
  )
  return(out)
}
plot_cdt <- function(cdt, show_digits = 3) {
  
  party_obj <- partykit::as.party(cdt$student_fit$fit)
  
  plt <- ggparty::ggparty(party_obj) +
    ggparty::geom_edge() +
    ggparty::geom_edge_label(
      ggplot2::aes(
        label = substr(breaks_label, start = 1, stop = 12 + show_digits)
      ),
      fill = "#BFBFBF",
      size = 2,
    ) +
    ggparty::geom_node_label(
      ggplot2::aes(label = splitvar),
      ids = "inner",
      fill = "white",
      size = 4,
      label.size = 0.1
    )
  
  ## === Subgroup ATE 资料 ===
  subgroup_ates <- data.frame(id = plt$data$id) |>
    dplyr::left_join(cdt$estimate, by = c("id" = "leaf_id")) |>
    dplyr::mutate(
      label = sprintf("Subgroup ATE\n= %.3f", estimate),
      sign  = estimate > 0,
      fill_col = dplyr::if_else(
        sign,
        "#F4D6D0",
        "#DDEBDC"
      )
    )
  
  ## 合并回 ggparty data
  plt$data <- plt$data |>
    dplyr::left_join(subgroup_ates, by = "id")
  
  ## === Terminal node label（背景变色，文字黑）===
  plt <- plt +
    ggparty::geom_node_label(
      ggplot2::aes(
        label = label,
        fill  = fill_col
      ),
      ids = "terminal",
      colour = "black",   # 文字维持黑色
      label.size = 0.4
    ) +
    ggplot2::scale_fill_identity()
  
  return(plt)
}
plot_jaccard <- function(...) {
  dots_ls <- rlang::dots_list(...)
  
  default_names <- paste0("Model", 1:length(dots_ls))
  if (is.null(names(dots_ls))) {
    names(dots_ls) <- default_names
  } else {
    names(dots_ls)[names(dots_ls) == ""] <- default_names[names(dots_ls) == ""]
  }
  
  ssi_df <- purrr::map(
    dots_ls,
    function(stability_diagnostics) {
      tibble::tibble(
        `Tree Depth` = 1:length(stability_diagnostics$jaccard_mean),
        `Jaccard SSI` = stability_diagnostics$jaccard_mean
      )
    }
  ) |>
    dplyr::bind_rows(.id = "Teacher Model")
  
  plt <- ggplot2::ggplot(ssi_df) +
    ggplot2::aes(x = `Tree Depth`, y = `Jaccard SSI`, color = `Teacher Model`) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::theme_classic()
  return(plt)
}


#EDVAE
X_train=edvae_out[,c(4:18)]
tauhat=as.data.frame(edvae_out$edvae_ite)
colnames(tauhat)="tauhat"
treat=as.data.frame(edvae_out$Treat)
colnames(treat)="treat"
Y_train=edvae_out$outcome
DR_train=edvae_out$pseudo.outcome_DR
R_train=edvae_out$pseudo.outcome_R
R_train_wt=edvae_out$pseudo.outcome_R_w
true_ite=edvae_out$ite

ground_true=F
have_pseudo="Y"
if(have_pseudo=="DR"){
  outcome_train=DR_train
  if(ground_true){
    outcome_train=true_ite
  }
  estimate_group_cates=estimate_group_cates_pseudo_DR
}else if(have_pseudo=="R"){
  outcome_train=R_train
  estimate_group_cates=estimate_group_cates_pseudo_R
}else{
  outcome_train=Y_train
  estimate_group_cates=estimate_group_cates_Y
}

student_fit_out=student_model(X_train,tauhat,prune="min")
if(have_pseudo=="R"){
  group_cates=estimate_group_cates(student_fit_out$fit,X_train,outcome_train,treat,R_train_wt)
}else{
  group_cates=estimate_group_cates(student_fit_out$fit,X_train,outcome_train,treat)
}
out <- list(
  estimate = group_cates,
  student_fit = student_fit_out,
  teacher_predictions = c(as.matrix(tauhat))
)
set.seed(85619)
stability_out <- evaluate_subgroup_stability(
  estimator = student_rpart,
  fit = student_fit_out$fit,
  X = X_train,
  y = tauhat,
  Z = NULL,
  prune="min",
  rpart_control = NULL,
  B = 100,
  max_depth = NULL
)
plot_cdt(out)
group_cates
stability_out$jaccard_mean
stability_out$feature_distribution

mean(edvae_out$ite)
mean(edvae_out$edvae_ite)
mean(edvae_out$pseudo.outcome_DR)
weighted.mean(edvae_out$pseudo.outcome_R,edvae_out$pseudo.outcome_R_w)


HTE_data = cbind(edvae_out$outcome,
                 edvae_out$Treat,
                 edvae_out$ite,
                 edvae_out$edvae_ite,
                 edvae_out$pseudo.outcome_DR,
                 edvae_out$pseudo.outcome_R,
                 edvae_out$pseudo.outcome_R_w,
                 out$student_fit$fit$where)
HTE_data = as.data.frame(HTE_data)
colnames(HTE_data)=c("outcome","Treat","ite","edvae_ite","DR_pseudo","R_pseudo","R_pseudo_wt","group")
HTE_data$group = as.factor(HTE_data$group)

HTE_est=glm(ite~group,data=HTE_data)
GATE_true=HTE_est$coefficients
GATE_true[-1]=GATE_true[-1]+GATE_true[1]

HTE_est=glm(R_pseudo~group,data=HTE_data,weights=R_pseudo_wt)
GATE=HTE_est$coefficients
HTE_coef=HTE_est$coefficients
GATE[-1]=GATE[-1]+GATE[1]
HTE_coef
GATE
GATE_true

#bootstrap for CI
library(dplyr)
B=10000
DR_coef=matrix(0,B,length(unique(out$student_fit$fit$where))-1)
R_coef=matrix(0,B,length(unique(out$student_fit$fit$where))-1)
DR_GATE=matrix(0,B,length(unique(out$student_fit$fit$where)))
R_GATE=matrix(0,B,length(unique(out$student_fit$fit$where)))
RC_GATE=matrix(0,B,length(unique(out$student_fit$fit$where)))
set.seed(85619)
for(b in c(1:B)){
  id=sample(c(1:dim(HTE_data)[1]),size=dim(HTE_data)[1],replace = T)
  HTE_data_b = HTE_data[id,]
  
  HTE_est=glm(DR_pseudo~group,data=HTE_data_b)
  DR_coef_b=HTE_est$coefficients
  DR_GATE_b=HTE_est$coefficients
  DR_GATE_b[-1]=DR_GATE_b[-1]+DR_GATE_b[1]
  
  HTE_est=glm(R_pseudo~group,data=HTE_data_b,weights=R_pseudo_wt)
  R_coef_b=HTE_est$coefficients
  R_GATE_b=HTE_est$coefficients
  R_GATE_b[-1]=R_GATE_b[-1]+R_GATE_b[1]
  
  data1=HTE_data_b[HTE_data_b$Treat==1,c(1,7)]
  data1=data1%>%
    group_by(group)%>%
    summarise(mean_y=mean(outcome))
  data0=HTE_data_b[HTE_data_b$Treat==0,c(1,7)]
  data0=data0%>%
    group_by(group)%>%
    summarise(mean_y=mean(outcome))
  if(length(data1$mean_y)==15 & length(data0$mean_y)==15){
    RC_GATE_b=data1$mean_y-data0$mean_y
  }else{
    RC_GATE_b=rep(NA,15)
  }
  
  DR_coef[b,]=DR_coef_b[-1]
  R_coef[b,]=R_coef_b[-1]
  DR_GATE[b,]=DR_GATE_b
  R_GATE[b,]=R_GATE_b
  RC_GATE[b,]=RC_GATE_b
  
  if((b/1000-floor(b/1000))==0){print(c("B =",b))}
}
sqrt(diag(var(DR_GATE)))
sqrt(diag(var(R_GATE)))
sqrt(diag(var(RC_GATE)))
colMeans(DR_GATE)
colMeans(R_GATE)
colMeans(RC_GATE)
as.numeric(GATE_true)
j=15
hist(R_GATE[,j],breaks=50,freq=F)
lines(seq(min(R_GATE[,j]),max(R_GATE[,j]),by=0.001),
      dnorm(seq(min(R_GATE[,j]),max(R_GATE[,j]),by=0.001),
            colMeans(R_GATE)[j],sqrt(diag(var(R_GATE)))[j]))


for(j in c(1:dim(DR_GATE)[2])){
  print(as.numeric(round(quantile(DR_GATE[,j],c(0.025,0.975)),4)))
}
for(j in c(1:dim(R_GATE)[2])){
  print(as.numeric(round(quantile(R_GATE[,j],c(0.025,0.975)),4)))
}
a=na.omit(RC_GATE)
for(j in c(1:dim(a)[2])){
  print(as.numeric(round(quantile(a[,j],c(0.025,0.975)),4)))
}


j=14
hist(R_coef[,j],breaks=50,freq=F)
lines(seq(min(R_coef[,j]),max(R_coef[,j]),by=0.001),
      dnorm(seq(min(R_coef[,j]),max(R_coef[,j]),by=0.001),
            colMeans(R_coef)[j],sqrt(diag(var(R_coef)))[j]))

HTE_data = cbind(edvae_out$Composite,
                 edvae_out$fluvac,
                 edvae_out$edvae_ite,
                 edvae_out$pseudo.outcome_DR,
                 edvae_out$pseudo.outcome_R,
                 edvae_out$pseudo.outcome_R_w,
                 out$student_fit$fit$where)
HTE_data = as.data.frame(HTE_data)
colnames(HTE_data)=c("outcome","Treat","edvae_ite","DR_pseudo","R_pseudo","R_pseudo_wt","group")
HTE_data$group = as.factor(HTE_data$group)

HTE_est=glm(DR_pseudo~group,data=HTE_data)
HTE_coef=HTE_est$coefficients
GATE=HTE_est$coefficients
GATE[-1]=GATE[-1]+GATE[1]
DR_HTE_coef=matrix((GATE-GATE[15])[-15],1,14)
DR_HTE_test = DR_HTE_coef%*%solve(var(DR_GATE)[-15,-15])%*%t(DR_HTE_coef)
DR_HTE_pvalue = 1-pchisq(DR_HTE_test,14)
DR_HTE_pvalue

HTE_est=glm(R_pseudo~group,data=HTE_data,weights=R_pseudo_wt)
HTE_coef=HTE_est$coefficients
GATE=HTE_est$coefficients
GATE[-1]=GATE[-1]+GATE[1]
R_HTE_coef=matrix((GATE-GATE[15])[-15],1,14)
R_HTE_test = R_HTE_coef%*%solve(var(R_GATE)[-15,-15])%*%t(R_HTE_coef)
R_HTE_pvalue = 1-pchisq(R_HTE_test,14)
R_HTE_pvalue

