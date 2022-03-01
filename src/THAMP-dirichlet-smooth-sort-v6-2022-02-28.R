rm(list=ls())
library(readxl)
library(tidyverse)
library(ggplot2)
library(spatstat) 
library(sparr)
library(RColorBrewer)
library(ggrepel)   
library(ggnewscale)
library(Cairo)
library(showtext)
library(optparse)
library(car)
library(xlsx)
library(stringr)
library(pROC)

mkdir<- function(p){
  if (! dir.exists(p)) {
    dir.create(p,recursive = TRUE)
  }
}



#command line parameters
option_list <- list(
  make_option(c('-i','--indir'), type = 'character', help = 'path of the results from PCA, eg.2d_full_data'),
  make_option(c('-d','--data'), type = 'character', help = 'original all feature data'),
  make_option(c('-o','--outdir'), type = 'character', help = 'output dir')
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

root<- "THAMP-data-and-code/src"
setwd(root)
opt$indir <- '../output/2d_full_data/'
opt$outdir <- '../output/'
opt$data <- '../data/data-2022-02-28.xlsx'




#-----------------------------------------------------------part1: read in the data, set the traits, karyotype and the number of columns where the mutation is located
if (is.null(opt$indir) | is.null(opt$outdir) | is.null(opt$data)){
  print_help(opt_parser)
  stop('missing parameters, please check carefully !', call.=FALSE)
}

indir <- opt$indir
outdir <- opt$outdir
mkdir(opt$outdir)
###start####
t1 <- Sys.time()
###Data reading and tiding####
datapath<- '../output/axis_rotation/myeloid'
savepath <- opt$outdir

cligene.data <- as.data.frame(read_excel(opt$data,sheet = 1))
MYE <- c('AML', 'MDS', 'MDS/MPN', 'MPN')
cligene.data.MYE <- cligene.data[cligene.data$Diagnosis %in% MYE, ]
LYM <- c('ALL', 'CLL', 'MM')
cligene.data.LYM <- cligene.data[cligene.data$Diagnosis %in% LYM, ]
tra.num.start <- 5
tra.num.end <- 96
kar.num.start <- 97
kar.num.end <- 168
mut.num.start <- 169
mut.num.end <- 342




#----------------------------------------------------------part2:  axis rotation for results from PCA----------
#-----axis rotation----------
axis_roation_path<- '../output/axis_rotation'
mkdir(axis_roation_path)
center_projection<- function(dd,key, type){
  print(type)
  if (type=="MYE") {
    group<- c("MPN","AML")
    g1<- "MPN"
    g2<- "AML"
  }else if (type=="LYM") {
    group<- c("MM","ALL")
    g1<- "MM"
    g2<- "ALL"
  }else{
    print("wrong type")
  }
  
  if (type=="MYE") {
    path=paste0(axis_roation_path,"/myeloid")
  }else if (type=="LYM") {
    path= paste0(axis_roation_path,"/lymphoid")
  }
  if (! dir.exists(path)){
    dir.create(path)
  }
  
  #
  get_center<- function(df){
    Z = array(0, dim = c(length(group),2))
    rownames(Z)<- group
    for (i in 1:length(group)) {
      dty<- df[df$Diagnosis ==  group[i],]
      ##calculate centers，or we can use the MASS::cov.trob funciton if we have a preview of the data distribution
      center<- cov.wt(dty[c('X','Y')])$center
      Z[i,]<- as.vector(center)
    }
    return(Z)
  }
  Z<- get_center(dd)
  
  x1 = Z[g1,1]; y1 = Z[g1,2]
  x2 = Z[g2,1];  y2 = Z[g2,2]
  k = (y2-y1)/ (x2-x1)
  b = y1-k*x1
  
  theta<- atan(k)
  if (x1>x2) {
    theta<- theta+pi
  }
  
  X_new<- c()
  Y_new<- c()
  for (i in 1: dim(dd)[1]) {
    x<- dd$X[i]
    y<- dd$Y[i]
    xt<- x*cos(theta) + y*sin(theta)
    yt<- y*cos(theta) - x*sin(theta)
    X_new<- c(X_new, xt)
    Y_new<- c(Y_new, yt)
  }
  
  dd_new<- data.frame(X=X_new, Y=Y_new, Diagnosis = dd$Diagnosis, ID=dd$ID)
  
  pdf(file = paste0(path,"/axis-rotation-by-center-location_",key,"_.pdf"),width = 10, height = 8)
  p<- ggplot(data = dd, aes(x=X, y=Y, col=Diagnosis))+geom_point()+theme_bw()+
    stat_ellipse(type = "norm",level = 0.95)+scale_color_manual(values = rainbow(4))
  print(p)
  
  p<- ggplot(data = dd_new, aes(x=X, y=Y, col=Diagnosis))+geom_point()+theme_bw()+
    stat_ellipse(type = "norm",level = 0.95)+scale_color_manual(values = rainbow(4))
  print(p)
  
  dataEllipse(dd$X, dd$Y,levels = 0.7, groups = as.factor(dd$Diagnosis))
  abline(b,k, col="green", lwd=2)
  
  dev.off()
  write.xlsx(dd_new, file = paste0(path,"/axis-rotation-by-center-location_",key,"_.xlsx"))
  return(x1>x2)
}


for (f in list.files(indir)) {
  print(f)
  key = str_split(f,"_")[[1]][5]
  df<- as.data.frame(read_excel(path = paste0(indir,"/",f), sheet = 1))
  dd<- df[c("X","Y","Diagnosis","ID")]
  
  type="LYM"
  dd1<- dd[dd$Diagnosis %in% c('ALL','CLL','MM'),]
  loccom<- center_projection(dd1,key=key, type=type)
  
  type="MYE"
  dd2<- dd[dd$Diagnosis %in% c('AML','MDS','MPN','MDS/MPN'),]
  loccom<- center_projection(dd2,key=key, type=type)
}

print("------------------------ part 2: axis rotation done !!-------------")



###--------------------------------------------part 3: dirichlet cluster--------------------------------------------------
#Perform Dirichlet clustering to generate Figure 2a and Figure 2b.
###organize data###
karyotype <- c(kar.num.start:kar.num.end)
mutation <- c(mut.num.start:mut.num.end)

mutation_or_karyotype<-karyotype
Mutation <- cligene.data[, c(1:5,mutation_or_karyotype)]

for (j in 6:ncol(Mutation)) {
  Mutation[,j] <- as.numeric(Mutation[,j])
}

genes <- colnames(Mutation)[6:ncol(Mutation)]
diseases <- unique(Mutation$Diagnosis)

diagnosis.n <- array(NA, length(diseases))
for (m in 1:length(diseases)) {
  diagnosis.n[m] <- sum(Mutation$Diagnosis == diseases[m])
}

diagnosis.pos <- array(NA, dim=c(length(genes),length(diseases)))
diagnosis.tot <- array(NA, dim=c(length(genes),length(diseases)))
for (j in 1:length(genes)) {
  for (m in 1:length(diseases)) {
    diagnosis.pos[j, m] <- sum(Mutation$Diagnosis == diseases[m] &
                                 Mutation[, which(genes==genes[j])+5] == 1,
                               na.rm = TRUE)
    diagnosis.tot[j, m] <- sum(Mutation$Diagnosis == diseases[m] &
                                 !is.na(Mutation[, which(genes==genes[j])+5]),
                               na.rm = TRUE)
  }
}
colnames(diagnosis.pos) <- diseases
rownames(diagnosis.pos) <- genes
colnames(diagnosis.tot) <- diseases
rownames(diagnosis.tot) <- genes


###dirichlet process clustering###

### helper functions ###
rs <- (0:.5e3)/1e3*2

binomial.likelihood.scan <- function(n, pos) {
  prob <- c()
  for (t in 1:length(rs)) {
    r <- rs[t]
    prob <- c(prob, r^pos * (1-r)^(n-pos))
  }
  prob
}

sample.beta <- function(n, pos) {
  log.prob <- c()
  for (t in 1:length(rs)) {
    r <- rs[t]
    if (pos==0) {
      log.prob <- c(log.prob, log(1-r)*(n-pos))
    }
    if (pos==n) {
      log.prob <- c(log.prob, log(r)*pos)
    }
    if (pos!=0 & pos!=n) {
      log.prob <- c(log.prob, log(r)*pos + log(1-r)*(n-pos))
    }
  }
  log.prob <- log.prob - max(log.prob)
  prob <- exp(log.prob)
  
  sample(rs, 1, prob=prob/sum(prob))
}

naive.log.likelihood <- function(diagnosis.tot.j, diagnosis.pos.j) {
  log.likelihood <- 0
  for (k in 1:length(diseases)) {
    log.likelihood <- log.likelihood +
      log(mean(binomial.likelihood.scan(diagnosis.tot.j[k], diagnosis.pos.j[k])))
  }
  log.likelihood
}

relabel.clusters <- function(cluster.labels, cluster.betas) {
  cluster.labels.original <- cluster.labels
  cluster.labels <- cluster.labels.original - 
    min(cluster.labels.original) + max(cluster.labels.original) + 1
  clusters <- sort(unique(cluster.labels), decreasing=FALSE)
  for (k in 1:length(clusters)) {
    cluster.labels[cluster.labels==clusters[k]] <- k
  }
  
  cluster.betas.original <- cluster.betas
  cluster.betas <- list()
  for (k in 1:length(clusters)) {
    old.label <- cluster.labels.original[cluster.labels==k][1]
    cluster.betas[[k]] <- cluster.betas.original[[old.label]]
  }
  
  list(cluster.labels=cluster.labels, cluster.betas=cluster.betas)
}

###
burn.in <- 2500
posterior.samples <- 2500

cluster.labels <- 1:length(genes)
cluster.betas <- list()
for (j in 1:length(genes)) {
  cluster.betas[[j]] <- array(NA, length(diseases)) 
  for (m in 1:length(diseases)) {
    cluster.betas[[j]][m] <- sample.beta(diagnosis.tot[j,m], diagnosis.pos[j,m])
  }
}

###
alpha <- 1e-4

naive.logs <- array(NA, length(genes))
for (j in 1:length(genes)) {
  naive.logs[j] <- naive.log.likelihood(diagnosis.tot[j,], diagnosis.pos[j,]) + log(alpha)
}

for (t in 1:(burn.in+posterior.samples)) {
  cat('\n', t, ':')
  for (j in sample(length(genes), length(genes), replace=FALSE)) {
    
    naive.log <- naive.logs[j]
    existing.clusters.log <- array(0, max(cluster.labels))
    for (k in 1:max(cluster.labels)) {
      beta <- cluster.betas[[k]]
      n <- diagnosis.n
      pos <- diagnosis.pos[j,]
      existing.clusters.log[k] <- 
        sum(log(beta^pos)) + sum(log((1-beta)^(n-pos))) + 
        log(sum(cluster.labels[-j] == k))
    }
    
    offset <- max(c(naive.log, existing.clusters.log))
    native.log <- naive.log - offset
    existing.clusters.log <- existing.clusters.log - offset
    
    prob <- exp(c(existing.clusters.log, native.log))
    if (sum(cluster.labels==cluster.labels[j])==1) {
      prob[cluster.labels[j]] <- 0
    }
    prob <- prob / sum(prob)
    label.j <- sample(1:(max(cluster.labels)+1), 1, prob=prob)
    
    if (label.j > max(cluster.labels)) {
      cluster.labels[j] <- label.j
      cluster.betas[[label.j]] <- array(NA, length(diseases)) 
      for (m in 1:length(diseases)) {
        cluster.betas[[label.j]][m] <- 
          sample.beta(diagnosis.tot[j,m], diagnosis.pos[j,m])
      }
    } else {
      cluster.labels[j] <- label.j
      cluster.betas[[label.j]] <- array(NA, length(diseases)) 
      for (m in 1:length(diseases)) {
        I <- cluster.labels==label.j
        cluster.betas[[label.j]][m] <- 
          sample.beta(diagnosis.tot[j,m]*sum(I), sum(diagnosis.pos[I,m]))
      }
    }
    
    relabel.out <- relabel.clusters(cluster.labels, cluster.betas)
    cluster.labels <- relabel.out$cluster.labels
    cluster.betas <- relabel.out$cluster.betas
    
    if (j %% 10 == 0) cat('.')
  }
  
  cat(': ', max(cluster.labels), '\n')
  cat(' ', sort(table(cluster.labels),decreasing=T)[1:50], '\n')
}

###
betas <- array(NA, dim=c(max(cluster.labels), length(diseases)))
for (k in 1:max(cluster.labels)) {
  betas[k,] <- cluster.betas[[k]]
}
rownames(betas) <- paste('gene cluster', 1:max(cluster.labels))
colnames(betas) <- diseases

if ( !dir.exists(paste0(savepath, '/dirichlet cluster/'))) {
  dir.create(paste0(savepath, '/dirichlet cluster/'), recursive = TRUE)
}

write.csv(betas, file=paste0(savepath, '/dirichlet cluster/kt-betas-2021-11-21.csv'))

###
for (label.j in 1:max(cluster.labels)) {
  cluster.betas[[label.j]] <- array(NA, length(diseases)) 
  for (m in 1:length(diseases)) {
    I <- cluster.labels==label.j
    cluster.betas[[label.j]][m] <- 
      sample.beta(diagnosis.tot[j,m]*sum(I), sum(diagnosis.pos[I,m]))
  }
}

###sample from posterior distribution###

cluster.labels.log <- c()
cluster.betas.log <- list()

for (t in 1:500) {
  cat('\n', t, ':')
  for (j in sample(length(genes), length(genes), replace=FALSE)) {
    naive.log <- naive.log.likelihood(diagnosis.n, diagnosis.pos[j,]) + log(alpha)
    existing.clusters.log <- array(0, max(cluster.labels))
    for (k in 1:max(cluster.labels)) {
      beta <- cluster.betas[[k]]
      n <- diagnosis.n
      pos <- diagnosis.pos[j,]
      existing.clusters.log[k] <- 
        sum(log(beta^pos)) + sum(log((1-beta)^(n-pos))) + 
        log(sum(cluster.labels[-j] == k))
    }
    
    offset <- max(c(naive.log, existing.clusters.log))
    native.log <- naive.log - offset
    existing.clusters.log <- existing.clusters.log - offset
    
    prob <- exp(c(existing.clusters.log, native.log))
    if (sum(cluster.labels==cluster.labels[j])==1) {
      prob[cluster.labels[j]] <- 0
    }
    prob <- prob / sum(prob)
    label.j <- sample(1:(max(cluster.labels)+1), 1, prob=prob)
    
    if (label.j > max(cluster.labels)) {
      cluster.labels[j] <- label.j
      cluster.betas[[label.j]] <- array(NA, length(diseases)) 
      for (m in 1:length(diseases)) {
        cluster.betas[[label.j]][m] <- 
          sample.beta(diagnosis.n[m], diagnosis.pos[j,m])
      }
    } else {
      cluster.labels[j] <- label.j
      cluster.betas[[label.j]] <- array(NA, length(diseases)) 
      for (m in 1:length(diseases)) {
        I <- cluster.labels==label.j
        cluster.betas[[label.j]][m] <- 
          sample.beta(diagnosis.n[m]*sum(I), sum(diagnosis.pos[I,m]))
      }
    }
    
    relabel.out <- relabel.clusters(cluster.labels, cluster.betas)
    cluster.labels <- relabel.out$cluster.labels
    cluster.betas <- relabel.out$cluster.betas
    
    if (j %% 10 == 0) cat('.')
  }
  
  cat(' ', max(cluster.labels))
  
  if (t %% 1 == 0) {
    betas <- array(NA, dim=c(max(cluster.labels), length(diseases)))
    for (k in 1:max(cluster.labels)) {
      betas[k,] <- cluster.betas[[k]]
    }
    
    cluster.labels.log <- cbind(cluster.labels.log, cluster.labels)
    cluster.betas.log[[t]] <- betas
  }
}

save(diseases, genes, cluster.labels.log, cluster.betas.log,
     file = paste0(savepath, '/dirichlet cluster/posterior-2021-11-21.RData'))

####
top1 <- array(NA, length(genes))
perc1 <- array(NA, length(genes))
top2 <- array(NA, length(genes))
perc2 <- array(NA, length(genes))
top3 <- array(NA, length(genes))
perc3 <- array(NA, length(genes))

for (j in 1:length(genes)) {
  tab <- sort(table(cluster.labels.log[j,]), decreasing=TRUE)
  if (length(tab)>=1) {
    top1[j] <- names(tab)[1]
    perc1[j] <- tab[1]/500
  }
  if (length(tab)>=2) {
    top2[j] <- names(tab)[2]
    perc2[j] <- tab[2]/500
  }
  if (length(tab)>=3) {
    top3[j] <- names(tab)[3]
    perc3[j] <- tab[3]/500
  }
}

write.csv(data.frame(gene=genes,
                     top1.cluster=top1,
                     perc1=perc1,
                     top2.cluster=top2,
                     perc2=perc2,
                     top3.cluster=top3,
                     perc3=perc3,
                     diagnosis.pos/diagnosis.tot),
          file= paste0(savepath, '/dirichlet cluster/label-sample-2021-11-21.csv'))

####
mean.betas <- array(NA, dim=c(max(cluster.labels), length(diseases)))
sd.betas <- array(NA, dim=c(max(cluster.labels), length(diseases)))
for (k in 1:max(cluster.labels)) {
  for (m in 1:length(diseases)) {
    temp <- c()
    for (t in 1:length(cluster.betas.log)) {
      temp[t] <- cluster.betas.log[[t]][k,m]
    }
    mean.betas[k,m] <- mean(temp)
    sd.betas[k,m] <- sd(temp)
  }
}

rownames(mean.betas) <- paste('gene cluster', 1:max(cluster.labels))
colnames(mean.betas) <- diseases

rownames(sd.betas) <- paste('gene cluster', 1:max(cluster.labels))
colnames(sd.betas) <- diseases

write.csv(mean.betas, file = paste0(savepath, '/dirichlet cluster/mean-betas-2021-11-21.csv'))
write.csv(sd.betas, file = paste0(savepath, '/dirichlet cluster/sd-betas-2021-11-21.csv'))

####
for (i in 1:max(cluster.labels)) {
  I <- top1==i & perc1>0.0
  cat('cluster', i, ": ")
  genes.subset <- genes[I]
  for (j in 1:length(genes.subset)) {
    cat(genes.subset[j], " ")
  }
  cat('\n')
}
I <- perc1<=0.0
genes.subset <- genes[I]
cat("Not confidently clustered: ")
for (j in 1:length(genes.subset)) {
  cat(genes.subset[j], " ")
}
cat('\n')


print("-----------------derichlet  cluster done -------------")


####------------------------------------------------part 4: functions for kernel estimation-------------------------####
#*******************************************************************************
#Judgment error type
is.error <- function(x) inherits(x, 'try-error')
#Digital format function
num.format <- function(x) {
  xx <- c()
  for (i in 1:length(x)) {
    if ((x[i] - floor(x[i])) < 1e-3) {
      xx[i] <- paste(x[i], '.0', sep='')
    } else {
      xx[i] <- as.character(x[i])
    }
  }
  xx
}
#kde funtion
kde <- function(x = x1, y = y1, lamda = NULL, adapt = T, data = clinical.single){
  dx <- data$X - x
  dy <- data$Y - y
  if (adapt) {
    sum(data$V * exp(-(dx ^ 2 + dy ^ 2)/(2 * lamda ^ 2)) / lamda ^ 2)/
      (5 + sum(exp(-(dx ^ 2 + dy ^ 2) / (2 * lamda ^ 2)) / lamda ^ 2))
  }else{
    sum(data$V * exp(-(dx ^ 2 + dy ^ 2) / (2 * lamda ^ 2)) / lamda ^ 2)
  }
}
#Clinical trait kde function
clinicalkde <- function(colnum, data = clinical.data, pl = T, adapt=T, point.pl=T, 
                        color = colorRampPalette(c('white', 'Orange3'))(50)){
  name.ind <- names(data)[colnum]
  cat(paste0(name.ind, '\n'))
  
  clinical.single <- data[, c('X', 'Y', name.ind)]
  names(clinical.single) <- c('X', 'Y', 'V')
  clinical.single <- clinical.single[clinical.single$V != '0', ]
  
  #Find global bandwidth h0
  sigma.x <- sd(clinical.single[clinical.single$V != '0', ]$X)
  sigma.v <- sd(clinical.single[clinical.single$V != '0', ]$V)
  sigma <- (sigma.x * sigma.v) ^ (1 / 2)
  n <- nrow(clinical.single[clinical.single$V != '0', ])
  RK <- 1 / (sqrt(2 * pi))
  h0 <- 3 * (RK / (35 * n)) ^ (1 / 5) * sigma #Maximum smoothness principle
  
  if (adapt) {
    #Calculate variable bandwidth #DOI:10.1002/sim.7577
    fi.data <- clinical.single
    fi.data$V <- ifelse(fi.data$V != 0, 1, 0) 
    fi <- 0                                                         #Using a fixed bandwidth, calculate the initial kernel density estimate for each point 
    for (i in 1:nrow(clinical.single)) {
      fi[i] <- kde(clinical.single$X[i], clinical.single$Y[i], h0, adapt = F, data = fi.data)
    }
    y <- exp(mean(log(1 / sqrt(fi))))                               #Find the geometric mean of the initial kernel density estimates for all points
    sqfi <- pmin(1 / sqrt(fi), 10 * y)                              #In order to prevent the bandwidth from being too large due to too small fi, a limit of 10*y is given
    hi <- h0 * sqfi / y                                             #Calculate variable bandwidth
  }else{
    hi <- h0
  }
  
  #Calculate the kernel density estimation matrix
  xseq <- seq(1, 100, by = 1)
  yseq <- seq(1, 100, by = 1)
  kde.data <- matrix(0, nrow = length(xseq), ncol = length(yseq))
  for (i in 1:length(xseq)) {                   
    for (j in 1:length(yseq)) {
      kde.data[i, j]<-kde(xseq[i], yseq[j], lamda = hi, adapt = T, data = clinical.single)  
    }
  }
  
  max.ind <- which(kde.data[, 50] == max(kde.data[, 50]), arr.ind = T)
  max.x <- max.ind[1]
  core.x <- sum(kde.data[, 50] %*% seq(1, 100)) / sum(kde.data[, 50])
  
  if (pl) {
    jpeg(filename = paste0(savepath, gsub('/|%', '', name.ind), '.jpg'))
    opar <- par(pin = c(4, 4))
    plot(1:100, kde.data[, 50], type = 'l',
         xlim = c(0, 100), ylim = c(0, 20),
         ann = T, xaxt = 'n',
         xlab = '', ylab = '',
         xaxs='i', yaxs='i',
         main = name.ind, cex.main = 3)
    par(new = T) 
    plot(meanx, 0, type = 'p', pch = 17,
         xlim = c(0, 100), ylim = c(0, 20),
         ann = F, xaxt = 'n', yaxt = 'n',
         xlab = '', ylab = '',
         xaxs = 'i', yaxs = 'i') 
    box()
    grid()
    axis(1, at = seq(0, 100, 10), seq(0, 100, 10))
    par(opar)
    dev.off()
  }
  if (point.pl) {
    ggplot(data=subset(clinical.single, V != '0'), aes(x = X, y = Y, colour = V)) +
      geom_point(alpha = 0.6) +
      scale_color_gradientn(colors = c('gold', 'orange', 'darkred')) +
      scale_x_continuous(limits = c(0, 100),breaks = seq(0, 100, 10),expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 100),breaks = seq(0, 100, 10),expand = c(0, 0)) +
      ggtitle(name.ind)+
      theme_bw()
    ggsave(filename = paste0(savepath, 'zpoint', gsub('/|%', '', name.ind), '.jpg'), 
           width = 6, height = 5, dpi = 600)
  }
  
  cat('Done.\n')
  c(max.x, core.x, kde.data[, 50])
} 
#Gene kde function
genekde <- function(colnum, data = clinical.data, epsilon = 0, color = colorRampPalette(c('white', 'blue'))(50), pl = T){
  
  name.ind <- names(data)[colnum]
  cat(paste0(name.ind, '\n'))
  
  #Select mutant genes
  data.mut <- data[data[, colnum] == 1, c('X','Y')]
  data.nmut <- data[data[, colnum] == 0, c('X','Y')]
  data.all <- data[, c('X', 'Y')]
  
  #Generate ppp objects based on map boundaries and mutant gene coordinates
  points.li <- list(structure(list(x.n = data.mut$X, y.n = data.mut$Y), class = 'data.frame', .Names = c('x.n', 'y.n')), 
                    structure(list(x.n = data.nmut$X, y.n = data.nmut$Y), class = 'data.frame', .Names = c('x.n', 'y.n')),
                    structure(list(x.n = data.all$X, y.n = data.all$Y),  class = 'data.frame', .Names = c('x.n', 'y.n'))) 
  data.a <- mapply(function(x, y) {ppp(x$x.n, x$y.n, window = y)}, x = points.li, y = bound.li, SIMPLIFY = FALSE) 
  
  h0 <- OS(data.a[[3]], nstar = 'npoints')
  
  gene.rr <- risk(data.a[[1]], data.a[[3]], h0 = h0, adapt = F, tolerate = F,log = F,
                  hp=c(OS(data.a[[3]]) / 2, OS(data.a[[3]]) / 2), pilot.symmetry = 'none', epsilon = 1, davies.baddeley = 0.05)
  
  gene.f <- nrow(data.mut) * gene.rr$f$z$v
  gene.g <- nrow(data.all) * gene.rr$g$z$v
  gene  <- gene.f / ((gene.g + epsilon * max(gene.g)))
  
  core.x <- sum(gene[65, ]%*%seq(1, 128))/sum(gene[65, ]) 
  
  if (pl) {
    jpeg(filename = paste0(datapath, gsub('/', '', name.ind), '.jpg'))
    opar <- par(pin = c(4, 4))
    plot(1:128, gene[65, ], type = 'l',
         xlim = c(0, 128), ylim = c(0, 0.2),
         ann = T, xaxt = 'n',
         xlab = '', ylab = '',
         xaxs = 'i', yaxs='i',
         main = name.ind, cex.main = 3)
    par(new = T) 
    plot(meanx, 0, type = 'p', pch = 17,
         xlim = c(0, 128), ylim = c(0, 0.2),
         ann = F, xaxt = 'n',yaxt = 'n',
         xlab = '',ylab = '',
         xaxs = 'i', yaxs = 'i')  
    box()
    grid()
    axis(1, at = seq(0, 128, 12.8), seq(0, 100, 10))
    par(opar)
    dev.off()
  }
  cat('Done.\n')
  
  c(core.x / 1.28, gene[65, ])
}
#Disease kde function
diseasekde <- function(diseasename, data = clinical.data, epsilon = 0, color = colorRampPalette(c('white', 'blue'))(50), pl = T){
  cat(paste0(diseasename, '\n'))
  
  data.mut <- data[data$Diagnosis == diseasename, c('X', 'Y')]
  data.nmut <- data[data$Diagnosis != diseasename, c('X', 'Y')]
  data.all <- data[, c('X', 'Y')]
  
  points.li <- list(structure(list(x.n = data.mut$X, y.n = data.mut$Y), class = 'data.frame', .Names = c('x.n', 'y.n')), 
                    structure(list(x.n = data.nmut$X, y.n = data.nmut$Y),  class = 'data.frame', .Names = c('x.n', 'y.n')),
                    structure(list(x.n = data.all$X, y.n = data.all$Y),  class = 'data.frame', .Names = c('x.n', 'y.n'))) 
  data.a <- mapply(function(x, y) {ppp(x$x.n, x$y.n, window = y)}, x = points.li, y = bound.li, SIMPLIFY = FALSE) 
  
  h0 <- OS(data.a[[3]], nstar = 'npoints')
  
  disease.rr <- risk(data.a[[1]], data.a[[3]], h0 = h0, adapt = F, tolerate = F,log = F,
                     hp=c(OS(data.a[[3]])/2, OS(data.a[[3]])/2), pilot.symmetry = 'none', epsilon = 1, davies.baddeley = 0.05)
  
  disease.f <- nrow(data.mut) * disease.rr$f$z$v
  disease.g <- nrow(data.all) * disease.rr$g$z$v
  disease  <- disease.f / ((disease.g + epsilon * max(disease.g)))
  
  core.x <- sum(disease[65, ] %*% seq(1, 128)) / sum(disease[65, ]) 
  
  if (pl) {
    jpeg(filename = paste0(datapath, gsub('/', '', diseasename), '.jpg'))
    opar <- par(pin = c(4, 4))
    plot(1:128, disease[65, ], type = 'l',
         xlim = c(0, 128), ylim = c(0, 0.2),
         ann = T, xaxt = 'n',
         xlab = '',ylab = '',
         xaxs = 'i',yaxs = 'i',
         main = diseasename, cex.main = 3)
    par(new = T) 
    plot(meanx, 0, type = 'p', pch = 17,
         xlim = c(0, 128),ylim = c(0, 0.2),
         ann = F, xaxt = 'n',yaxt = 'n',
         xlab = '',ylab = '',
         xaxs = 'i',yaxs = 'i')  
    box()
    grid()
    axis(1, at = seq(0, 128, 12.8), seq(0, 100, 10))
    par(opar)
    dev.off()
  }
  cat('Done.\n')
  
  c(core.x / 1.28, disease[65,])
}
#Disease 2D kde function
diseasekde2d <- function(diseasename, data = clinical.data, epsilon = 0, color = colorRampPalette(c("white", "blue"))(50), pl=T){
  cat(paste0(diseasename, '\n'))
  
  data.mut <- data[data$Diagnosis == diseasename, c('X', 'Y')]
  data.nmut <- data[data$Diagnosis != diseasename, c('X', 'Y')]
  data.all <- data[, c('X', 'Y')]
  
  points.li <- list(structure(list(x.n = data.mut$X, y.n = data.mut$Y), class = "data.frame", .Names = c("x.n", "y.n")), 
                    structure(list(x.n = data.nmut$X, y.n = data.nmut$Y),  class = "data.frame", .Names = c("x.n", "y.n")),
                    structure(list(x.n = data.all$X, y.n = data.all$Y),  class = "data.frame", .Names = c("x.n", "y.n"))) 
  data.a <- mapply(function(x, y) {ppp(x$x.n, x$y.n, window = y)}, x = points.li, y = bound.li, SIMPLIFY = FALSE) 
  
  h0 <- OS(data.a[[3]], nstar = "npoints")
  
  disease.rr <- risk(data.a[[1]], data.a[[3]], h0 = h0, adapt = F, tolerate = F,log = F,
                     hp=c(OS(data.a[[3]]) / 2, OS(data.a[[3]]) / 2), pilot.symmetry = "none", epsilon = 1, davies.baddeley = 0.05)
  
  disease.f <- nrow(data.mut) * disease.rr$f$z$v
  disease.g <- nrow(data.all) * disease.rr$g$z$v
  disease  <- disease.f / ((disease.g + epsilon * max(disease.g)))
  
  max.ind <- which(disease == max(disease),arr.ind = T)
  
  max.x <- max.ind[2]
  max.y <- max.ind[1]
  
  mean.x <- colMeans(data.mut,na.rm = T)[1]
  mean.y <- colMeans(data.mut,na.rm = T)[2]
  
  
  core.x <- sum(disease %*% seq(1, 128)) / sum(disease) 
  core.y <- sum(t(disease) %*% seq(1, 128)) / sum(disease)
  
  if (pl) {
    jpeg(filename = paste0(savepath, '/MYE/2d/disease-plot/', gsub('/', '', diseasename), '.jpg'))
    disease[1, 1] <- 1
    opar <- par(usr = c(0, 100, 0, 100), pin = c(4, 4))
    image(t(disease), col=color, ann = F, xaxt = "n", yaxt = "n")
    par(new = T) 
    plot(max.x / 1.28, max.y / 1.28, type = 'p', pch = 20, xlim = c(0, 100), ylim = c(0, 100)) 
    box()
    grid()
    title(main = diseasename)
    axis(1, at = seq(0, 100, 10) / 100, seq(0, 100, 10))
    axis(2, at = seq(0, 100, 10) / 100, seq(0, 100, 10))
    par(opar)
    dev.off()
  }
  cat('Done.\n')
  
  write.csv(disease.f,paste0(savepath, '/MYE/2d/f-', gsub('/', '', diseasename), '.csv'))
  write.csv(disease, paste0(savepath, '/MYE/2d/', gsub('/', '', diseasename), '.csv'))
  
  c(mean.x, mean.y, max.x / 1.28, max.y / 1.28, core.x / 1.28, core.y / 1.28, disease[order(disease, decreasing = TRUE)[2]], nrow(data.mut) / nrow(data))
}
#Gene 2D kde function
genekde2d <- function(colnum,data = clinical.data, epsilon = 0, color = colorRampPalette(c("white", "blue"))(50), pl = T){
  name.ind <- names(data)[colnum]
  cat(paste0(name.ind, '\n'))
  
  data.mut <- data[data[, colnum] == 1, c('X', 'Y')]
  data.nmut <- data[data[, colnum] == 0, c('X', 'Y')]
  data.all <- data[, c('X', 'Y')]
  
  points.li <- list(structure(list(x.n = data.mut$X, y.n = data.mut$Y), class = "data.frame", .Names = c("x.n", "y.n")), 
                    structure(list(x.n = data.nmut$X, y.n = data.nmut$Y),  class = "data.frame", .Names = c("x.n", "y.n")),
                    structure(list(x.n = data.all$X, y.n = data.all$Y),  class = "data.frame", .Names = c("x.n", "y.n"))) 
  data.a <- mapply(function(x, y) {ppp(x$x.n, x$y.n, window = y)}, x = points.li, y = bound.li, SIMPLIFY = FALSE) 
  
  h0 <- OS(data.a[[3]], nstar = "npoints")
  
  gene.rr <- risk(data.a[[1]], data.a[[3]], h0 = h0, adapt = F, tolerate = F,log = F,
                  hp=c (OS(data.a[[3]]) / 2, OS(data.a[[3]]) / 2), pilot.symmetry = "none", epsilon = 1, davies.baddeley = 0.05)
  
  gene.f <- nrow(data.mut) * gene.rr$f$z$v
  gene.g <- nrow(data.all) * gene.rr$g$z$v
  #epsilon <- 2
  gene  <- gene.f / ((gene.g + epsilon * max(gene.g)))
  
  max.ind <- which(gene == max(gene), arr.ind = T)
  max.x <- max.ind[2]
  max.y <- max.ind[1]
  
  mean.x <- colMeans(data.mut, na.rm = T)[1]
  mean.y <- colMeans(data.mut, na.rm = T)[2]
  
  core.x <- sum(gene %*% seq(1, 128)) / sum(gene) 
  core.y <- sum(t(gene) %*% seq(1, 128)) / sum(gene)
  
  if (pl) {
    jpeg(filename = paste0(savepath, '/MYE/2d/gene-plot/', gsub('/', '', name.ind), '.jpg'))
    gene[1, 1] <- 0.2
    opar <- par(usr = c(0, 100, 0, 100), pin = c(4, 4))
    image(t(gene), col = color, ann = F, xaxt = "n", yaxt = "n")
    par(new = T) 
    plot(max.x / 1.28, max.y / 1.28, type = 'p', pch = 20, xlim = c(0, 100),ylim=c(0, 100))  
    box()
    grid()
    title(main = name.ind)
    axis(1, at = seq(0, 100, 10) / 100, seq(0, 100, 10))
    axis(2, at = seq(0, 100, 10) / 100, seq(0, 100, 10))
    par(opar)
    dev.off()
  }
  cat('Done.\n')
  
  c(mean.x, mean.y, max.x / 1.28, max.y / 1.28, core.x / 1.28, core.y / 1.28, gene[order(gene, decreasing = TRUE)[2]], nrow(data.mut) / nrow(data))
}
print("-------------------define the functions for calculating kernal density done !-----------------")


###2D####
####------------------------------------------------part 5: 2d map for Fig3,ab -------------------------####
#Generate the two-dimensional projection map of the disease (Figure 3a) and the barycentric projection map of the two-dimensional kernel density estimation of the gene (Figure 3b); the results are stored in the 2d folder.
###data reading and tiding###
#cluster.data <- read_xlsx('./data/rawdata/encode_100_65_bs_0_2d.xlsx', sheet = 1)
cluster.data <- read_xlsx('../output/2d_full_data/encode_100_65_bs_0_2d.xlsx', sheet = 1)

filenumber <- '2d'
names(cluster.data)[names(cluster.data) == "Diagnosis"] <- "Diagnosis"

xmin <- floor(range(cluster.data$X)[1] / 10) * 10
xmax <- ceiling(range(cluster.data$X)[2] / 10) * 10
ymin <- floor(range(cluster.data$Y)[1] / 10) * 10
ymax <- ceiling(range(cluster.data$Y)[2] / 10) * 10
cluster.data$X <- 100 * (cluster.data$X - xmin) / (xmax - xmin)
cluster.data$Y <- 100 * (cluster.data$Y - ymin) / (ymax - ymin)

clinical.data <- merge(cluster.data[, c('X', 'Y', 'ID')], cligene.data, by = 'ID', all.x = T)

ref.li <- list(structure(list(x.ref = c(0, 0, 100, 100), y.ref = c(100, 0, 0, 100)),
                         class = "data.frame", row.names = c(NA, -4L), .Names = c("x.ref", "y.ref")),
               structure(list(x.ref = c(0, 0, 100, 100), y.ref = c(100, 0, 0, 100)),
                         class = "data.frame", row.names = c(NA, -4L), .Names = c("x.ref", "y.ref")),
               structure(list(x.ref = c(0, 0, 100, 100), y.ref = c(100, 0, 0, 100)),
                         class = "data.frame", row.names = c(NA, -4L), .Names = c("x.ref", "y.ref")))
bound.li <- lapply(ref.li, function(x) {owin(poly = list(x = x$x.ref, y = x$y.ref))})

if (! dir.exists(paste0(savepath, '/MYE/2d/disease-plot/'))) {
  dir.create(paste0(savepath, '/MYE/2d/disease-plot/'), recursive = TRUE)
}

###disease###
disease <- c('AML', 'MDS', 'MDS/MPN', 'MPN', 'ALL', 'CLL', 'MM', 'CML')
disease2d.data <- sapply(disease, function(x) try(diseasekde2d(x, epsilon = 0.1, pl = F)))
colnames(disease2d.data) <- disease

if (! dir.exists(paste0(savepath, '/MYE/2d/csv_data/'))) {
  dir.create(paste0(savepath, '/MYE/2d/csv_data/'), recursive = TRUE)
}
write.csv(disease2d.data,paste0(savepath, '/MYE/2d/csv_data/', 'MYE_disease2D_', filenumber, '.csv'))
failed.disease <- sapply(disease2d.data, is.error)
print(disease[failed.disease])
if (sum(failed.disease) > 0) {
  disease2d.data <- as.data.frame(disease2d.data[!failed.disease])
  names(disease2d.data) <- disease[!failed.disease]
  write.csv(disease2d.data, paste0(savepath, '/MYE/2d/csv_data/', 'MYE_disease2D_', filenumber, '.csv'))
}

getPalette.AML <- colorRampPalette(c("transparent", "blue"))
getPalette.MDS <- colorRampPalette(c("transparent", "red"))
getPalette.MPN <- colorRampPalette(c("transparent", "green"))
getPalette.MDSMPN <- colorRampPalette(c("transparent", "black"))

getPalette.MM <- colorRampPalette(c("transparent", "goldenrod"))
getPalette.ALL <- colorRampPalette(c("transparent", "cyan"))
getPalette.CLL <- colorRampPalette(c("transparent", "darkorchid1"))
getPalette.CML <- colorRampPalette(c("transparent", "magenta"))

gg.AML <- ggplot()+
  geom_density2d_filled(data = clinical.data[clinical.data$Diagnosis == 'AML', ],
                        mapping = aes(X, Y), alpha = 0.5, contour_var = "ndensity") + theme_bw() +
  scale_fill_manual(values = getPalette.AML(10)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(c(0, 100)) + xlim(c(0, 100)) +   #The x and y axis boundaries are adjusted according to the graph                                         
  ggtitle(label = 'AML')
gg.AML
ggsave(gg.AML, file = paste0(savepath, '/MYE/2d/disease-plot/AML.pdf'), width = 10, height = 10)

gg.MDS <- ggplot()+
  geom_density2d_filled(data = clinical.data[clinical.data$Diagnosis == 'MDS', ],
                        mapping = aes(X, Y), alpha = 0.5, contour_var = "ndensity") + theme_bw() +
  scale_fill_manual(values = getPalette.MDS(10)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(c(0, 100)) + xlim(c(0, 100)) +
  ggtitle(label = 'MDS')
gg.MDS
ggsave(gg.MDS, file = paste0(savepath, '/MYE/2d/disease-plot/MDS.pdf'), width = 10, height = 10)

gg.MPN <- ggplot() +
  geom_density2d_filled(data = clinical.data[clinical.data$Diagnosis == 'MPN', ],
                        mapping = aes(X, Y), alpha = 0.5, contour_var = "ndensity") + theme_bw() +
  scale_fill_manual(values = getPalette.MPN(10)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylim(c(0, 100)) + xlim(c(0, 100)) +
  ggtitle(label = 'MPN')
gg.MPN
ggsave(gg.MPN, file = paste0(savepath, '/MYE/2d/disease-plot/MPN.pdf'), width = 10, height = 10)

gg.MDSMPN <- ggplot()+
  geom_density2d_filled(data = clinical.data[clinical.data$Diagnosis == 'MDS/MPN',],
                        mapping = aes(X, Y), alpha = 0.5, contour_var = "ndensity") + theme_bw() +
  scale_fill_manual(values = getPalette.MDSMPN(10)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylim(c(0, 100)) + xlim(c(0, 100)) +
  ggtitle(label = 'MDS/MPN')
gg.MDSMPN
ggsave(gg.MDSMPN, file = paste0(savepath, '/MYE/2d/disease-plot/MDSMPN.pdf'), width = 10, height = 10)

gg.MM <- ggplot() +
  geom_density2d_filled(data = clinical.data[clinical.data$Diagnosis == 'MM', ],
                        mapping = aes(X, Y), alpha = 0.5, contour_var = "ndensity") + theme_bw() +
  scale_fill_manual(values = getPalette.MM(10)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(c(0, 100)) + xlim(c(0, 100)) +
  ggtitle(label = 'MM')
gg.MM
ggsave(gg.MM, file = paste0(savepath, '/MYE/2d/disease-plot/MM.pdf'), width = 10, height = 10)

gg.ALL <- ggplot() +
  geom_density2d_filled(data = clinical.data[clinical.data$Diagnosis == 'ALL', ],
                        mapping = aes(X, Y), alpha = 0.5, contour_var = "ndensity") + theme_bw() +
  scale_fill_manual(values = getPalette.ALL(10)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(c(0, 100)) + xlim(c(0, 100)) +
  ggtitle(label = 'ALL')
gg.ALL
ggsave(gg.ALL, file = paste0(savepath, '/MYE/2d/disease-plot/ALL.pdf'), width = 10, height = 10)

gg.CLL <- ggplot() +
  geom_density2d_filled(data = clinical.data[clinical.data$Diagnosis == 'CLL', ],
                        mapping = aes(X, Y), alpha = 0.5, contour_var = "ndensity") + theme_bw() +
  scale_fill_manual(values = getPalette.CLL(10)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(c(0, 100)) + xlim(c(0, 100)) +
  ggtitle(label = 'CLL')
gg.CLL
ggsave(gg.CLL, file = paste0(savepath, '/MYE/2d/disease-plot/CLL.pdf'), width = 10, height = 10)

gg.CML <- ggplot() +
  geom_density2d_filled(data = clinical.data[clinical.data$Diagnosis == 'CML', ],
                        mapping = aes(X, Y), alpha = 0.5, contour_var = "ndensity") + theme_bw() +
  scale_fill_manual(values = getPalette.CML(10)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(c(0, 100)) + xlim(c(0, 100)) +
  ggtitle(label = 'CML')
gg.CML
ggsave(gg.CML, file = paste0(savepath, '/MYE/2d/disease-plot/CML.pdf'), width = 10, height = 10)

###gene###
clinical.data.MYE <- merge(cluster.data[, c('X', 'Y', 'ID')], cligene.data.MYE, by = 'ID', all.y = T)

clinical.data.LYM <- merge(cluster.data[,c('X', 'Y', 'ID')], cligene.data.LYM, by = 'ID', all.y = T)

if (! dir.exists(paste0(savepath, '/MYE/2d/gene-plot/'))) {
  dir.create(paste0(savepath, '/MYE/2d/gene-plot/'), recursive = TRUE)
}

#Myeloid 
gene2d.data.MYE <- sapply(seq((kar.num.start + 2), (mut.num.end + 2)), function(x) try(genekde2d(x, data = clinical.data.MYE, epsilon = 2, pl = F)))
colnames(gene2d.data.MYE) <- names(clinical.data.MYE)[(kar.num.start + 2):(mut.num.end + 2)]
write.csv(gene2d.data.MYE, paste0(savepath, '/MYE/2d/csv_data/', 'MYE_gene2D_', filenumber, '.csv'))
failed <- sapply(gene2d.data.MYE, is.error)
print(colnames(clinical.data.MYE)[(kar.num.start + 2):(mut.num.end + 2)][failed])
if (sum(failed) > 0) {
  gene2d.data.MYE <- as.data.frame(gene2d.data.MYE[!failed])
  names(gene2d.data.MYE) <- colnames(clinical.data.MYE)[!failed]
  write.csv(gene2d.data.MYE, paste0(savepath, '/MYE/2d/csv_data/', 'MYE_gene2D_', filenumber, '.csv'))
}

#Lymphoid
gene2d.data.LYM <- sapply(seq((kar.num.start + 2), (mut.num.end + 2)), function(x) try(genekde2d(x, data = clinical.data.LYM, epsilon = 2, pl = F)))
colnames(gene2d.data.LYM) <- names(clinical.data.LYM)[(kar.num.start + 2):(mut.num.end + 2)]
write.csv(gene2d.data.LYM, paste0(savepath, '/MYE/2d/csv_data/', 'LYM_gene2D_', filenumber, '.csv'))
failed <- sapply(gene2d.data.LYM, is.error)
print(colnames(clinical.data.LYM)[(kar.num.start + 2):(mut.num.end + 2)][failed])
if (sum(failed) > 0){
  gene2d.data.LYM <- as.data.frame(gene2d.data.LYM[!failed])
  names(gene2d.data.LYM) <- colnames(clinical.data.LYM)[!failed]
  write.csv(gene2d.data.LYM, paste0(savepath, '/MYE/2d/csv_data/', 'LYM_gene2D_', filenumber, '.csv'))
}

ind.MYE <- read.csv(paste0(savepath, '/MYE/2d/csv_data/', 'MYE_gene2D_', filenumber, '.csv'))[, -1]
colnames(ind.MYE) <- colnames(gene2d.data.MYE)
ind.LYM <- read.csv(paste0(savepath, '/MYE/2d/csv_data/', 'LYM_gene2D_', filenumber, '.csv'))[, -1]
colnames(ind.LYM) <- colnames(gene2d.data.LYM)

ind.data.MYE <- as.data.frame(t(ind.MYE[, !(colnames(ind.MYE) %in% c("MRD"))])) 
names(ind.data.MYE) <- c('x1', 'y1', 'x2', 'y2', 'x', 'y', 'value', 'radio')

ind.data.MYE$type <- c(rep(1, (kar.num.end - kar.num.start + 1)), rep(2, nrow(ind.data.MYE) - (kar.num.end - kar.num.start + 1)))
write.csv(ind.data.MYE, paste0(savepath, '/MYE/2d/csv_data/inddata_MYE.csv')) 
radio  <- ind.data.MYE$radio
radio[!is.na(ind.data.MYE$radio) & ind.data.MYE$radio < 0.01] <- (radio[!is.na(ind.data.MYE$radio) & ind.data.MYE$radio < 0.01] / 0.01) ^ 1.5 * 0.01 
ind.data.MYE$radio <- radio

plot.data.MYE <- ind.data.MYE
I.MYE <- order(ind.data.MYE$radio,decreasing = T)[1:25]
write.csv(plot.data.MYE,paste0(savepath, '/MYE/2d/csv_data/', 'MYE_plot2D_', filenumber, '.csv'))

ind.data.LYM <- as.data.frame(t(ind.LYM[, !(colnames(ind.MYE) %in% c("MRD"))])) 
names(ind.data.LYM) <- c('x1', 'y1', 'x2', 'y2', 'x', 'y', 'value', 'radio')
ind.data.LYM$type <- c(rep(1, (kar.num.end - kar.num.start + 1)), rep(2, nrow(ind.data.LYM) - (kar.num.end - kar.num.start + 1)))
write.csv(ind.data.LYM, paste0(savepath, '/MYE/2d/csv_data/inddata_LYM.csv')) 
radio  <- ind.data.LYM$radio
radio[!is.na(ind.data.LYM$radio) & ind.data.LYM$radio < 0.01] <- (radio[!is.na(ind.data.LYM$radio) & ind.data.LYM$radio < 0.01] / 0.01) ^ 1.5 * 0.01
ind.data.LYM$radio <- radio

plot.data.LYM <- ind.data.LYM
I.LYM <- order(ind.data.LYM$radio,decreasing = T)[1:15]
write.csv(plot.data.LYM, paste0(savepath, '/MYE/2d/csv_data/', 'LYM_plot2D_', filenumber, '.csv'))

set.seed(01)
gg.gene <- ggplot() +
  geom_point(data = plot.data.LYM, 
             mapping = aes(x, y, 
                           color = factor(type, labels = c('karyotype', 'mutation')),
                           size = radio), alpha = 0.6) +
  guides(color = guide_legend(''), size = guide_legend("frequency")) + theme_bw() +
  xlab(label = "PC1") + ylab(label = 'PC2') + ggtitle('Center of gravity projection') +
  scale_color_manual(values = c("orange", "darkgreen")) +
  ylim(c(0, 100)) + xlim(c(0, 100)) +
  scale_size(limits = c(0, max(c(plot.data.LYM$radio, plot.data.MYE$radio))+0.01), range = c(0.01, 8)) + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

gg.gene <- gg.gene+new_scale('color') +
  geom_point(data = plot.data.MYE, 
             mapping = aes(x, y, 
                           color = factor(type, labels = c('karyotype', 'mutation')),
                           size = radio),alpha = 0.6) +
  scale_color_manual(values = c("red", "blue"))
gg.gene
ggsave(gg.gene, file = paste0(savepath, '/MYE/2d/gene-plot/gene.pdf'), width = 8, height = 8)

gg.genetext <- ggplot() +
  geom_point(data = plot.data.LYM, 
             mapping = aes(x, y, 
                           color = factor(type, labels = c('karyotype', 'mutation')),
                           size = radio), alpha = 0.6) +
  guides(color = guide_legend(''), size = guide_legend("frequency")) + theme_bw() +
  xlab(label = "PC1") + ylab(label = 'PC2') + ggtitle('Center of gravity projection') +
  scale_color_manual(values = c("orange", "darkgreen")) +
  ylim(c(0, 100)) + xlim(c(0, 100)) +
  geom_label_repel(aes(plot.data.LYM[I.LYM, ]$x, plot.data.LYM[I.LYM, ]$y,
                       label = rownames(plot.data.LYM[I.LYM, ])), size = 6, max.overlaps = Inf, fontface = "bold", color = "black",
                   box.padding = unit(2, "lines"), point.padding = unit(0, "lines"),
                   segment.colour = "grey50") + 
  scale_size(limits = c(0, max(c(plot.data.LYM$radio, plot.data.MYE$radio))+0.01), range = c(0.01, 8)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

gg.genetext <- gg.genetext + new_scale('color')+
  geom_point(data = plot.data.MYE, 
             mapping = aes(x, y, 
                           color = factor(type, labels = c('karyotype', 'mutation')),
                           size = radio), alpha = 0.6) +
  geom_label_repel(aes(plot.data.MYE[I.MYE, ]$x, plot.data.MYE[I.MYE, ]$y, label = rownames(plot.data.MYE[I.MYE, ])), size = 6, max.overlaps = Inf,
                   fontface = "bold", color = "black",
                   box.padding = unit(2, "lines"), point.padding = unit(0, "lines"),
                   segment.colour = "grey50") +
  scale_color_manual(values = c("red", "blue"))
gg.genetext
ggsave(gg.genetext, file = paste0(savepath, '/MYE/2d/gene-plot/genetext.pdf'), width = 8, height = 8)

print("---------------------2d map done !--------------")


####-------------------------------------------------------part 6: 100 boostrap data to generate a sorted heat map  -------------------------####
#*******************************************************************************
###100 boostrap####

file_list <- list.files(path = datapath, pattern = '*.xlsx$', recursive = TRUE,full.names = TRUE)
for(j in 1:length(file_list)){  
  cluster.data <- read_xlsx(file_list[j], sheet = 1)
  filenumber <- sub('.*_([0-9]+)_.*', '\\1', file_list[j], perl = TRUE) 
  clinical.data <- merge(cluster.data[, c('X','Y','ID')], cligene.data.MYE, by = 'ID')
  
  xmin <- floor(range(clinical.data$X)[1] / 10) * 10
  xmax <- ceiling(range(clinical.data$X)[2] / 10) * 10
  ymin <- floor(range(clinical.data$Y)[1] / 10) * 10
  ymax <- ceiling(range(clinical.data$Y)[2] / 10) * 10
  
  clinical.data$X <- 100 * (clinical.data$X - xmin) / (xmax - xmin)
  clinical.data$Y <- 100 * (clinical.data$Y - ymin) / (ymax - ymin) * 0 + 50 
  
  ref.li <- list(structure(list(x.ref = c(0, 0, 100, 100), y.ref = c(100, 0, 0, 100)),
                           class = 'data.frame', row.names = c(NA, -4L), .Names = c('x.ref', 'y.ref')),
                 structure(list(x.ref = c(0, 0, 100, 100), y.ref = c(100, 0, 0, 100)),
                           class = 'data.frame', row.names = c(NA, -4L), .Names = c('x.ref', 'y.ref')),
                 structure(list(x.ref = c(0, 0, 100, 100), y.ref = c(100, 0, 0, 100)),
                           class = 'data.frame', row.names = c(NA, -4L), .Names = c('x.ref', 'y.ref')))
  
  bound.li <- lapply(ref.li, function(x) {owin(poly = list(x = x$x.ref, y = x$y.ref))})
  
  ###disease###
  disease.data <- sapply(MYE, function(x) try(diseasekde(x, epsilon = 0.1, pl = F)))
  colnames(disease.data) <- MYE
  if ( !dir.exists(paste0(savepath, '/MYE/csv_data/'))) {
    dir.create(paste0(savepath, '/MYE/csv_data/'), recursive = TRUE)
  }
  
  write.csv(disease.data, paste0(savepath, '/MYE/csv_data/', 'MYE_disease1D_', filenumber, '.csv'))
  failed.disease <- sapply(disease.data, is.error)
  print(MYE[failed.disease])
  if (sum(failed.disease)>0) {
    disease.data <- as.data.frame(disease.data[!failed.disease])
    names(disease.data) <- MYE[!failed.disease]
    write.csv(disease.data, paste0(savepath, '/MYE/csv_data/', 'MYE_disease1D_', filenumber, '.csv'))
  }
  
  ###disease sort####
  
  x <- seq(0.5, 100, 0.78125)
  df <- disease.data[-1, ]
  for (i in 1:ncol(df)) {
    df[, i] <- as.numeric(df[, i])
  }
  
  if (! dir.exists(paste0(savepath, '/MYE/disease/'))) {
    dir.create(paste0(savepath, '/MYE/disease/'), recursive = TRUE)
  }
  
  jpeg(filename = paste0(savepath,'/MYE/disease/', 'mye disease sort', filenumber, '.jpg'),
       width = 700 * 3, height = 3 * 600, res = 72 * 3)
  
  heatmap(t(as.matrix(df[, c(4, 3, 2, 1)])), Rowv = NA, Colv = NA, keep.dendro = TRUE, 
          labRow = rev(MYE), labCol = num.format(round(x * 10) / 10), xlab = 'Myeloid axis')
  dev.off()
  
  ###trait###
  clinical1d.data<-sapply(seq((tra.num.start + 2), (tra.num.end + 2)),function(x) try(clinicalkde(x, adapt = T, point.pl = F, pl = F)))
  colnames(clinical1d.data)<- names(clinical.data)[(tra.num.start + 2):(tra.num.end + 2)]
  write.csv(clinical1d.data, paste0(savepath, '/MYE/csv_data/', 'MYE_clinical1D_', filenumber, '.csv'))
  failed.all <- sapply(clinical1d.data, is.error)
  print(names(clinical.data)[(tra.num.start + 2):(tra.num.end + 2)][failed.all])
  if (sum(failed.all) > 0) {
    clinical1d.data <- as.data.frame(clinical1d.data[!failed.all])
    names(clinical1d.data) <- colnames(clinical1d.data)[!failed.all]
    write.csv(clinical1d.data, paste0(savepath, '/MYE/csv_data/', 'MYE_clinical1D_', filenumber, '.csv'))
  }
  
  ###trait order###  
  center.of.mass <- as.numeric(clinical1d.data[1, ])
  traits <- colnames(clinical1d.data)
  x <- seq(0.5, 99.5)
  df <- clinical1d.data[-1:-2, ]
  for (i in 1:ncol(df)) {
    df[, i] <- as.numeric(df[, i])
  }
  
  df[is.na(df)] <- 0
  I <- colSums(df) != 0
  df <- df[, I]
  traits <- traits[I]
  center.of.mass <- center.of.mass[I]
  
  for (i in 1:length(traits)) {
    df[,i] <- df[,i] / mean(df[, i])
  }
  df[df > 5] <- 5
  
  peak.of.mass <- array(NA, length(traits))
  peak.val <- array(NA, length(traits))
  for (i in 1:length(traits)) {
    peak.of.mass[i] <- x[which.max(df[, i])]
    peak.val[i] <- max(df[, i])
  }
  
  pos.score <- center.of.mass
  sort.order <- sort(pos.score[peak.val > 0], index.return = TRUE)$ix 
  
  if (!dir.exists(paste0(savepath, '/MYE/clinical/'))) {
    dir.create(paste0(savepath, '/MYE/clinical/'), recursive = TRUE)
  }
  
  jpeg(filename = paste0(savepath, '/MYE/clinical/', 'mye clinical sort', filenumber, '.jpg'),
       width = 600 * 3,height = 3 * 600,res = 72 * 3)
  
  heatmap(t(as.matrix(df[, peak.val > 0][, sort.order])), Rowv = NA, Colv = NA, keep.dendro = TRUE, 
          labRow = traits[peak.val > 0][sort.order], labCol = num.format(round(x * 10) / 10), xlab = 'Myeloid axis',
          margins = c(5, 8),cexRow = 0.6 )
  dev.off()
  
  traits.order  <- traits[sort(pos.score, index.return = TRUE)$ix]
  clinicalpeak.order <- peak.of.mass[sort(pos.score, index.return = TRUE)$ix]
  clinicalpeak.data <- cbind(traits.order, clinicalpeak.order)
  write.csv(clinicalpeak.data, paste0(savepath, '/MYE/csv_data/', 'MYE_clinicalpeak_', filenumber, '.csv'))
  
  ###gene###
  
  gene1d.data <- sapply(seq((kar.num.start + 2), (mut.num.end + 2)), function(x) try(genekde(x, epsilon = 0.1, pl = F)))
  colnames(gene1d.data) <- names(clinical.data)[(kar.num.start + 2):(mut.num.end + 2)]
  write.csv(gene1d.data, paste0(savepath, '/MYE/csv_data/', 'MYE_gene1D_', filenumber, '.csv'))
  failed <- sapply(gene1d.data, is.error)
  print(colnames(clinical1d.data)[(kar.num.start + 2):(mut.num.end + 2)][failed])
  if (sum(failed) > 0) {
    gene1d.data <- as.data.frame(gene1d.data[!failed])
    names(gene1d.data) <- colnames(clinical1d.data)[!failed]
    write.csv(ind.data, paste0(savepath, '/MYE/csv_data/', 'MYE_gene1D_', filenumber, '.csv'))
  }
  
  ###gene order###
  center.of.mass <- as.numeric(gene1d.data[1, ])
  aberrations <- colnames(gene1d.data)
  x <- seq(0.5, 100, 0.78125)
  df <- gene1d.data[-1, ]
  for (i in 1:ncol(df)) {
    df[, i] <- as.numeric(df[, i])
  }
  
  I <- colSums(is.na(df)) == 0
  df <- df[, I]
  aberrations <- aberrations[I]
  center.of.mass <- center.of.mass[I]
  
  peak.of.mass <- array(NA, length(aberrations))
  peak.val <- array(NA, length(aberrations))
  for (i in 1:length(aberrations)) {
    peak.of.mass[i] <- x[which.max(df[, i])]
    peak.val[i] <- max(df[, i])
  }
  
  pos.score <- center.of.mass
  sort.order <- sort(pos.score[peak.val > 0.015], index.return = TRUE)$ix 
  
  if (! dir.exists(paste0(savepath, '/MYE/gene/'))) {
    dir.create(paste0(savepath, '/MYE/gene/'), recursive = TRUE)
  }
  
  jpeg(filename = paste0(savepath, '/MYE/gene/', 'mye gene sort', filenumber, '.jpg'),
       width = 600 * 3, height = 3 * 600, res = 72 * 3)
  
  heatmap(t(as.matrix(df[, peak.val > 0.015][, sort.order])), Rowv = NA, Colv = NA, keep.dendro = TRUE, 
          labRow = aberrations[peak.val > 0.015][sort.order], labCol = num.format(round(x * 10) / 10), xlab = 'Myeloid axis',
          margins = c(5, 8),cexRow = 0.7)
  dev.off()
  
  aberrations.order  <- aberrations[sort(pos.score, index.return = TRUE)$ix]
  genepeak.order <- peak.of.mass[sort(pos.score, index.return = TRUE)$ix]
  genepeak.data <- cbind(aberrations.order, genepeak.order)
  write.csv(genepeak.data, paste0(savepath, '/MYE/csv_data/', 'MYE_genepeak_', filenumber, '.csv'))
}
print("------------------------------100 boostrap data done !!!!!---------")





####-------------------------------------------------------part 7: Calculate the average of the sort  -------------------------####

#*******************************************************************************
###Calculate the average of the sort####
###ID###
file_list <- list.files(path = datapath, pattern = '*.xlsx$', recursive = TRUE, full.names = TRUE)
disease.id <- cligene.data[cligene.data$Diagnosis %in% MYE, c('ID', 'Diagnosis')]
idrank <- data.frame(array(NA, dim = c(nrow(disease.id), length(file_list))))
idrank <- cbind(disease.id, idrank)

for(j in 1:length(file_list)){
  cluster.data <- read_xlsx(file_list[j], sheet = 1)
  cluster.data$rank <- rank(cluster.data$X)
  cluster.data <- cluster.data[!duplicated(cluster.data$ID), ]
  idrank1 <- merge(disease.id, cluster.data[, c('ID', 'rank')], by = 'ID', all.x = T)
  idrank[,j+2] <- idrank1$rank
}

if (! dir.exists(paste0(savepath, '/MYE/mean_test/'))) {
  dir.create(paste0(savepath, '/MYE/mean_test/'), recursive = TRUE)
}
write.csv(idrank,paste0(savepath, '/MYE/mean_test/', 'MYE_idrank.csv'))

idorderdata <- read.csv(paste0(savepath, '/MYE/mean_test/', 'MYE_idrank.csv'))
sdrank <- apply(idorderdata[, -1:-3], 1, sd, na.rm = TRUE)
medianrank <- apply(idorderdata[, -1:-3], 1, median, na.rm = TRUE)
meanrank <- rowMeans(idorderdata[, -1:-3], na.rm = T)
idrankdata <- cbind(idorderdata[, 1:3], meanrank, medianrank, sdrank)
write.csv(idrankdata, paste0(savepath, '/MYE/mean_test/', 'MYE_id_rank.csv'))

#Calculate the stability of sorting
idorderdata <- as.data.frame(idorderdata)
iddata <- idorderdata[, -1]

for (k in 1:nrow(iddata)) {
  iddata[k, ][is.na(iddata[k, ])] <- mean(as.numeric(iddata[k, -1:-2]), na.rm = T)
}

idkendall.data <- matrix(NA,nrow = ncol(iddata)-2, ncol = ncol(iddata) - 2)
for (i in 1:(ncol(iddata) - 2)){
  for (j in 1:(ncol(iddata) - 2)) {
    idkendall.data[i, j] <- cor(iddata[, 2 + i], iddata[, 2 + j], method = 'kendall')
  }
}
mean(idkendall.data[idkendall.data < 1]) 
quantile(idkendall.data[idkendall.data < 1],c(0.025, 0.975))

###clinical trait###

file_list <- list.files(path = paste0(savepath, '/MYE/csv_data/'), pattern = '*clinicalpeak.*.csv$', recursive = TRUE, full.names = TRUE)  #获得csv文件列表
order.data <- data.frame(clinicalname = names(cligene.data.MYE)[tra.num.start:tra.num.end])
for(j in 1:length(file_list)){
  cluster.data <- read.csv(file_list[j])
  clinicalname.data <- cluster.data[, c('X', 'traits.order')]
  order.data <- merge(order.data, clinicalname.data, by.x = 'clinicalname', by.y = 'traits.order', all = T)
}
write.csv(order.data, paste0(savepath, '/MYE/mean_test/', 'MYE_clinicalorderdata.csv'))


clinicalorderdata <- read.csv(paste0(savepath, '/MYE/mean_test/', 'MYE_clinicalorderdata.csv'))
I <- rowSums(!is.na(clinicalorderdata[, -1:-2])) != 0

trait.number <- colSums(cligene.data.MYE[, clinicalorderdata$clinicalname[I]] != 0, na.rm = T)
AMLtrait.number <- colSums(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'AML', clinicalorderdata$clinicalname[I]] != 0, na.rm = T)
MDStrait.number <- colSums(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'MDS', clinicalorderdata$clinicalname[I]] != 0, na.rm = T)
MPNtrait.number <- colSums(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'MPN', clinicalorderdata$clinicalname[I]] != 0, na.rm = T)
MDSMPNtrait.number <- colSums(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'MDS/MPN', clinicalorderdata$clinicalname[I]] != 0, na.rm = T)

trait.sd <- sapply(cligene.data.MYE[, clinicalorderdata$clinicalname[I]], sd)
AMLtrait.sd <- sapply(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'AML', clinicalorderdata$clinicalname[I]], sd, na.rm = T)
MDStrait.sd <- sapply(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'MDS', clinicalorderdata$clinicalname[I]], sd, na.rm = T)
MPNtrait.sd <- sapply(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'MPN', clinicalorderdata$clinicalname[I]], sd, na.rm = T)
MDSMPNtrait.sd <- sapply(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'MDS/MPN', clinicalorderdata$clinicalname[I]], sd, na.rm = T)

trait.mean <- colMeans(cligene.data.MYE[,clinicalorderdata$clinicalname[I]], na.rm = T)
AMLtrait.mean <- colMeans(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'AML', clinicalorderdata$clinicalname[I]], na.rm = T)
MDStrait.mean <- colMeans(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'MDS', clinicalorderdata$clinicalname[I]], na.rm = T)
MPNtrait.mean <- colMeans(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'MPN', clinicalorderdata$clinicalname[I]], na.rm = T)
MDSMPNtrait.mean <- colMeans(cligene.data.MYE[cligene.data.MYE$Diagnosis %in% 'MDS/MPN', clinicalorderdata$clinicalname[I]], na.rm = T)

trait.cv <- trait.sd / trait.mean
AMLtrait.cv <- AMLtrait.sd / AMLtrait.mean
MDStrait.cv <- MDStrait.sd / MDStrait.mean
MPNtrait.cv <- MPNtrait.sd / MPNtrait.mean
MDSMPNtrait.cv <- MDSMPNtrait.sd / MDSMPNtrait.mean

rankclinical <- clinicalorderdata[I, ] 
sdrank <- apply(rankclinical[, -1:-2], 1, sd, na.rm=T)
medianrank <- apply(rankclinical[, -1:-2], 1, median, na.rm=T)
meanrank <- rowMeans(rankclinical[, -1:-2], na.rm = T)
clinicalrankdata <- cbind(clinicalname=rankclinical$clinicalname, meanrank, medianrank, sdrank,
                          trait.number, AMLtrait.number, MDStrait.number, MPNtrait.number, MDSMPNtrait.number,
                          trait.sd, AMLtrait.sd, MDStrait.sd, MPNtrait.sd, MDSMPNtrait.sd,
                          trait.mean, AMLtrait.mean, MDStrait.mean, MPNtrait.mean, MDSMPNtrait.mean,
                          trait.cv, AMLtrait.cv, MDStrait.cv, MPNtrait.cv, MDSMPNtrait.cv)
write.csv(clinicalrankdata, paste0(savepath, '/MYE/mean_test/', 'MYE_clinical_rank.csv'))


clinicalrankdata <- as.data.frame(clinicalrankdata)
clinicalrankdata.sc <- clinicalrankdata[as.numeric(clinicalrankdata$sdrank) <= quantile(as.numeric(clinicalrankdata$sdrank), 0.5) & 
                                          clinicalrankdata[, 'trait.cv'] > 0.5&                                                     
                                          clinicalrankdata$clinicalname %in% names(cligene.data.MYE)[tra.num.start:tra.num.end], ]
I <- rowSums(is.na(clinicalorderdata[, -1:-2])) == 0 & clinicalorderdata$clinicalname %in% clinicalrankdata.sc$clinicalname
clinicaldata <- clinicalorderdata[I, ]
clikendall.data <- matrix(NA, nrow = ncol(clinicaldata) - 2, ncol = ncol(clinicaldata) - 2)
for (i in 1:(ncol(clinicaldata) - 2)){
  for (j in 1:(ncol(clinicaldata) - 2)) {
    clikendall.data[i, j]<-cor(clinicaldata[, 2 + i], clinicaldata[, 2 + j], method = 'kendall')
  }
}
mean(clikendall.data[clikendall.data < 1]) 
quantile(clikendall.data[clikendall.data < 1], c(0.025, 0.975))

###gene###

file_list <- list.files(path = paste0(savepath, '/MYE/csv_data/'), pattern = '*genepeak.*.csv$', recursive = TRUE, full.names = TRUE)  
order.data <- data.frame(genename = names(cligene.data.MYE)[kar.num.start:mut.num.end])
for(j in 1:length(file_list)){
  cluster.data <- read.csv(file_list[j])
  genename.data <- cluster.data[, c('X', 'aberrations.order')]
  order.data <- merge(order.data, genename.data, by.x = 'genename', by.y = 'aberrations.order', all = T)
}
write.csv(order.data, paste0(savepath, '/MYE/mean_test/', 'MYE_geneorderdata.csv'))

geneorderdata <- read.csv(paste0(savepath, '/MYE/mean_test/', 'MYE_geneorderdata.csv'))
geneordio <- colSums(cligene.data[cligene.data$Diagnosis %in% MYE, kar.num.start:mut.num.end], na.rm = T) > 0 

AMLmutationnumber <- colSums(cligene.data[cligene.data$Diagnosis %in% 'AML', kar.num.start:mut.num.end], na.rm = T)[geneordio]
MDSmutationnumber <- colSums(cligene.data[cligene.data$Diagnosis %in% 'MDS', kar.num.start:mut.num.end], na.rm = T)[geneordio]
MPNmutationnumber <- colSums(cligene.data[cligene.data$Diagnosis %in% 'MPN', kar.num.start:mut.num.end], na.rm = T)[geneordio]

mutationnumber <- colSums(cligene.data[cligene.data$Diagnosis %in% MYE, kar.num.start:mut.num.end], na.rm = T)[geneordio]
geneorderdata1 <- geneorderdata[geneorderdata$genename %in% names(geneordio[geneordio == T]), ]
rankgene <- sapply(geneorderdata1[, -1:-2], rank)
rankgene <- cbind(geneorderdata1$genename, rankgene)
I <- rowSums(is.na(rankgene[, -1])) == 0
genedata <- as.data.frame(rankgene[I, ])
genedata.num <- as.data.frame(sapply(genedata[, -1], as.numeric))
sdrank <- apply(genedata.num, 1, sd, na.rm = T)
medianrank <- apply(genedata.num, 1, median, na.rm = T)
meanrank <- rowMeans(genedata.num, na.rm = T)
generankdata <- cbind(genename = geneorderdata1$genename, meanrank, medianrank, sdrank,
                      mutationnumber = mutationnumber[geneorderdata1$genename],
                      AMLmutationnumber = AMLmutationnumber[geneorderdata1$genename],
                      MDSmutationnumber = MDSmutationnumber[geneorderdata1$genename],
                      MPNmutationnumber = MPNmutationnumber[geneorderdata1$genename])
write.csv(generankdata, paste0(savepath, '/MYE/mean_test/', 'MYE_gene_rank.csv'))

generankdata <- as.data.frame(generankdata)
generankdata.sc <- generankdata[generankdata$genename %in% names(cligene.data)[mut.num.start:mut.num.end], ]
generankdata.sc <- generankdata.sc[as.numeric(generankdata.sc$sdrank) <= quantile(as.numeric(generankdata.sc$sdrank), 1) &
                                     as.numeric(generankdata.sc$mutationnumber) > nrow(cligene.data.MYE) * 0.01, ]
genedata.sc <- genedata[genedata$V1 %in% generankdata.sc$genename, ]
genedata.sc.num <- as.data.frame(sapply(genedata.sc[, -1], as.numeric))
kendall.data <- matrix(NA, nrow = ncol(genedata.sc.num), ncol = ncol(genedata.sc.num))
for (i in 1:(ncol(genedata.sc.num))){
  for (j in 1:(ncol(genedata.sc.num))) {
    kendall.data[i, j] <- cor(genedata.sc.num[, i], genedata.sc.num[, j], method = 'kendall')
  }
}
mean(kendall.data[kendall.data < 1])                      
quantile(kendall.data[kendall.data < 1], c(0.025, 0.975)) 

print("--------------------------Calculate the average of the sort  done!!------------------")



####-------------------------------------------------------part 8: generate a ranking heat map of one-dimensional kernel density estimates of stable genes and features for all patients  -------------------------####

#*******************************************************************************
###All patient results####
header.names <- colnames(cligene.data)
header.names <- as.data.frame(header.names)

if (!dir.exists(paste0(savepath, '/MYE/all/'))) {
  dir.create(paste0(savepath, '/MYE/all/'), recursive = TRUE)
}

##List of stable genes and traits
gene_ranks.raw <- read.csv(paste0(savepath, '/MYE/mean_test/', 'MYE_gene_rank.csv'))
gene_ranks.raw <- arrange(gene_ranks.raw, medianrank)
gene_ranks.raw <- gene_ranks.raw[gene_ranks.raw$genename %in% header.names[mut.num.start:mut.num.end, ], ] 
write.csv(gene_ranks.raw, paste0(savepath, "/MYE/all/MYE_gene_rankraw.csv"))

trait_ranks.raw <- read.csv(paste0(savepath, '/MYE/mean_test/', 'MYE_clinical_rank.csv'))
trait_ranks.raw <- arrange(trait_ranks.raw, medianrank)

gene_ranks <- gene_ranks.raw[gene_ranks.raw$sdrank <= quantile(gene_ranks.raw$sdrank, 1), ]
gene_ranks <- gene_ranks[gene_ranks[, 'mutationnumber'] > nrow(cligene.data.MYE) * 0.01, ] 
trait_ranks <- trait_ranks.raw[trait_ranks.raw$sdrank <= quantile(trait_ranks.raw$sdrank, 0.5) & 
                                 trait_ranks.raw$clinicalname %in% header.names[tra.num.start:tra.num.end, ], ]
trait_ranks <- trait_ranks[trait_ranks[, 'trait.cv'] > 0.5 & !is.na(trait_ranks[, 'trait.cv']), ] 

write.csv(gene_ranks, paste0(savepath, "/MYE/all/MYE_gene_rank.csv"))
write.csv(trait_ranks, paste0(savepath, "/MYE/all/MYE_trait_rank.csv"))

stable.tg <- list(genename = gene_ranks$genename, traitname = trait_ranks$clinicalname)

cluster.data.raw <- read_xlsx('../output/axis_rotation/myeloid/axis-rotation-by-center-location_0_.xlsx', sheet = 1)

cluster.data <- cluster.data.raw
filenumber <- 'allpatient' 
MYE <- c('AML', 'MDS', 'MDS/MPN', 'MPN')
cluster.data <- cluster.data[cluster.data$Diagnosis %in% MYE, ]

xmin <- floor(range(cluster.data$X, na.rm = T)[1] / 10) * 10
xmax <- ceiling(range(cluster.data$X, na.rm = T)[2] / 10) * 10
ymin <- floor(range(cluster.data$Y, na.rm = T)[1] / 10) * 10
ymax <- ceiling(range(cluster.data$Y, na.rm = T)[2] / 10) * 10

cluster.data$X <- 100 * (cluster.data$X - xmin) / (xmax - xmin)
cluster.data$Y <- 100 * (cluster.data$Y - ymin) / (ymax - ymin)

clinical.data <- merge(cluster.data[, c('X', 'Y', 'ID')], cligene.data.MYE, by = 'ID')
clinical.data$Y <- clinical.data$Y * 0 + 50 

ref.li <- list(structure(list(x.ref = c(0, 0, 100, 100), y.ref = c(100, 0, 0, 100)),
                         class = "data.frame", row.names = c(NA, -4L), .Names = c("x.ref", "y.ref")),
               structure(list(x.ref = c(0, 0, 100, 100), y.ref = c(100, 0, 0, 100)),
                         class = "data.frame", row.names = c(NA, -4L), .Names = c("x.ref", "y.ref")),
               structure(list(x.ref = c(0, 0, 100, 100), y.ref = c(100, 0, 0, 100)),
                         class = "data.frame", row.names = c(NA, -4L), .Names = c("x.ref", "y.ref")))

bound.li <- lapply(ref.li, function(x){owin(poly = list(x = x$x.ref, y = x$y.ref))})

###disease###
disease.data <- sapply(MYE, function(x) try(diseasekde(x, epsilon = 0.1, pl = F)))
colnames(disease.data) <- MYE
if (!dir.exists(paste0(savepath, '/MYE/all/csv_data/'))) {
  dir.create(paste0(savepath, '/MYE/all/csv_data/'), recursive = TRUE)
}
write.csv(disease.data,paste0(savepath, '/MYE/all/csv_data/', 'MYE_disease1D_', filenumber, '.csv'))
failed.disease <- sapply(disease.data, is.error)
print(MYE[failed.disease])
if (sum(failed.disease) > 0) {
  disease.data <- as.data.frame(disease.data[!failed.disease])
  names(disease.data) <- MYE[!failed.disease]
  write.csv(disease.data, paste0(savepath, '/MYE/all/csv_data/', 'MYE_disease1D_', filenumber, '.csv'))
}

###disease order###

x <- seq(0.5, 100, 0.78125)
df <- disease.data[-1, ]
for (i in 1:ncol(df)) {
  df[,i] <- as.numeric(df[, i])
}

if (! dir.exists(paste0(savepath, '/MYE/all/disease/'))) {
  dir.create(paste0(savepath, '/MYE/all/disease/'), recursive = TRUE)
}

CairoPDF(file = paste0(savepath, '/MYE/all/disease/', 'mye disease order', filenumber, '.pdf'), width = 10, height = 10)
showtext_begin()
heatmap(t(as.matrix(df[, c(4, 3, 2, 1)])), Rowv = NA, Colv = NA, keep.dendro = TRUE, 
        labRow = rev(MYE), labCol = num.format(round(x * 10) / 10), xlab = 'Myeloid axis')
showtext_end()
dev.off()

write.csv(t(as.matrix(df[, c(4, 3, 2, 1)])), paste0(savepath, '/MYE/all/csv_data/', 'MYE_diseaseheat_', filenumber, '.csv'))

###clinical trait###
clinical1d.data <- sapply(seq((tra.num.start + 2), (tra.num.end + 2)), function(x) try(clinicalkde(x, adapt = T, point.pl = F, pl = F)))
colnames(clinical1d.data) <- names(clinical.data)[(tra.num.start + 2):(tra.num.end + 2)]
write.csv(clinical1d.data, paste0(savepath, '/MYE/all/csv_data/', 'MYE_clinical1D_', filenumber, '.csv'))
failed.all <- sapply(clinical1d.data, is.error)
print(names(clinical.data)[(tra.num.start + 2):(tra.num.end + 2)][failed.all])
if (sum(failed.all) > 0) {
  clinical1d.data <- as.data.frame(clinical1d.data[!failed.all])
  names(clinical1d.data) <- colnames(clinical1d.data)[!failed.all]
  write.csv(clinical1d.data, paste0(savepath, '/MYE/all/csv_data/', 'MYE_clinical1D_', filenumber, '.csv'))
}

###clinical trait order###  
center.of.mass <- as.numeric(clinical1d.data[1, ])
traits <- colnames(clinical1d.data)
x <- seq(0.5, 99.5)
df <- clinical1d.data[-1:-2, ]
for (i in 1:ncol(df)) {
  df[,i] <- as.numeric(df[, i])
}

df[is.na(df)] <- 0
I <- colSums(df) != 0 & colnames(df)%in%stable.tg$traitname #Pick stable traits
df <- df[, I]
traits <- traits[I]
center.of.mass <- center.of.mass[I]

for (i in 1:length(traits)) {
  df[,i] <- df[, i] / mean(df[, i])
}
df[df > 5] <- 5

peak.of.mass <- array(NA, length(traits))
peak.val <- array(NA, length(traits))
for (i in 1:length(traits)) {
  peak.of.mass[i] <- x[which.max(df[, i])]
  peak.val[i] <- max(df[, i])
}

pos.score <- center.of.mass
sort.order <- sort(pos.score[peak.val > 0], index.return = TRUE)$ix 

if (! dir.exists(paste0(savepath, '/MYE/all/clinical/'))) {
  dir.create(paste0(savepath, '/MYE/all/clinical/'), recursive = TRUE)
}

CairoPDF(file = paste0(savepath, '/MYE/all/clinical/', 'mye clinical order-', filenumber, '.pdf'), width = 10, height = 10)
showtext_begin()
heatmap(t(as.matrix(df[, peak.val > 0][, sort.order])), Rowv = NA, Colv = NA, keep.dendro = TRUE, 
        labRow = traits[peak.val > 0][sort.order], labCol = num.format(round(x * 10) / 10), xlab = 'Myeloid axis',
        margins = c(5, 8), cexRow = 0.6)
showtext_end()
dev.off()

write.csv(t(as.matrix(df[, peak.val > 0][, sort.order])), paste0(savepath, '/MYE/all/csv_data/', 'MYE_clinicalheat_', filenumber, '.csv'))

traits.order  <- traits[sort(pos.score, index.return = TRUE)$ix]
clinicalpeak.order <- peak.of.mass[sort(pos.score, index.return = TRUE)$ix]
clinicalpeak.data <- cbind(traits.order, clinicalpeak.order)
write.csv(clinicalpeak.data, paste0(savepath, '/MYE/all/csv_data/', 'MYE_clinicalpeak-', filenumber, '.csv'))

###gene### 

gene1d.data <- sapply(seq((kar.num.start + 2), (mut.num.end + 2)), function(x) try(genekde(x, epsilon = 2, pl = F)))
colnames(gene1d.data) <- names(clinical.data)[(kar.num.start + 2):(mut.num.end + 2)]
write.csv(gene1d.data, paste0(savepath, '/MYE/all/csv_data/', 'MYE_gene1D_', filenumber, '.csv'))
failed <- sapply(gene1d.data, is.error)
print(colnames(gene1d.data)[(kar.num.start + 2):(mut.num.end + 2)][failed])
if (sum(failed) > 0) {
  gene1d.data <- as.data.frame(gene1d.data[!failed])
  names(gene1d.data) <- colnames(clinical1d.data)[!failed]
  write.csv(gene1d.data, paste0(savepath, '/MYE/all/csv_data/', 'MYE_gene1D_', filenumber, '.csv'))
}

###gene order###
center.of.mass <- as.numeric(gene1d.data[1, ])
aberrations <- colnames(gene1d.data)
x <- seq(0.5, 100, 0.78125)
df <- gene1d.data[-1, ]
for (i in 1:ncol(df)) {
  df[, i] <- as.numeric(df[, i])
}

I <- colSums(is.na(df)) == 0 & colnames(df) %in% stable.tg$genename #Pick stable genes
df <- df[, I]
aberrations <- aberrations[I]
center.of.mass <- center.of.mass[I]

for (i in 1:length(aberrations)) {
  df[, i] <- df[, i] / mean(df[, i])
}
df[df>5] <- 5

peak.of.mass <- array(NA, length(aberrations))
peak.val <- array(NA, length(aberrations))
for (i in 1:length(aberrations)) {
  peak.of.mass[i] <- x[which.max(df[, i])]
  peak.val[i] <- max(df[, i])
}

pos.score <- center.of.mass
sort.order <- sort(pos.score[peak.val > 0.015], index.return = TRUE)$ix 

if (! dir.exists(paste0(savepath, '/MYE/all/gene/'))) {
  dir.create(paste0(savepath, '/MYE/all/gene/'), recursive = TRUE)
}

CairoPDF(file = paste0(savepath, '/MYE/all/gene/', 'mye gene order-', filenumber, '.pdf'), width = 10, height = 10)
showtext_begin()
heatmap(t(as.matrix(df[, peak.val > 0.015][, sort.order])), Rowv = NA, Colv = NA, keep.dendro = TRUE, 
        labRow = aberrations[peak.val > 0.015][sort.order], labCol = num.format(round(x * 10) / 10), xlab = 'Myeloid axis',
        margins = c(5, 8), cexRow = 0.7)
showtext_end()
dev.off()

write.csv(t(as.matrix(df[, peak.val > 0.015][, sort.order])), paste0(savepath, '/MYE/all/csv_data/', 'MYE_geneheat_', filenumber, '.csv'))

aberrations.order  <- aberrations[sort(pos.score, index.return = TRUE)$ix]
genepeak.order <- peak.of.mass[sort(pos.score, index.return = TRUE)$ix]
center.order <- center.of.mass[sort(pos.score, index.return = TRUE)$ix]
genepeak.data <- cbind(aberrations.order, genepeak.order, center.order)
write.csv(genepeak.data, paste0(savepath, '/MYE/all/csv_data/', 'MYE_genepeak_', filenumber, '.csv'))
print("----------------------------All patient results done!-----------")




####-------------------------------------------------------part 9: Group testing of gene and traits , to generate Figure 5a,5b -------------------------####
#*******************************************************************************

###Group testing of gene and traits####
cligene.data1 <- apply(cligene.data[, -1:-(tra.num.start - 1)], 2, as.numeric)
cligene.data <- cbind(cligene.data[, 1:(tra.num.start - 1)], cligene.data1)
cligene.data.MYE1 <- apply(cligene.data.MYE[, -1:-(tra.num.start - 1)], 2, as.numeric)
cligene.data.MYE <- cbind(cligene.data.MYE[, 1:(tra.num.start - 1)], cligene.data.MYE1)
cligene.data.LYM1 <- apply(cligene.data.LYM[, -1:-(tra.num.start - 1)], 2, as.numeric)
cligene.data.LYM <- cbind(cligene.data.LYM[, 1:(tra.num.start - 1)], cligene.data.LYM1)

header.names <- colnames(cligene.data.MYE)
header.names <- as.data.frame(header.names)

gene_ranks.raw <- read.csv(paste0(savepath, '/MYE/mean_test/', "MYE_gene_rank.csv"))
trait_ranks.raw <- read.csv(paste0(savepath, '/MYE/mean_test/', "MYE_clinical_rank.csv"))
trait_ranks.raw <- arrange(trait_ranks.raw, medianrank)


genes.MR <- c('DNMT3A', 'TET2', 'CREBBP', 'NOTCH1', 'SETD2', 'TP53')
genes.R <- c('KRAS', 'NRAS', 'CEBPA-dm', 'CEBPA-sm', 'FLT3-ITD', 'GATA2', 'IDH2', 'KIT', 'NPM1', 'WT1',
             'CSF3R', 'FLT3-Other', 'FLT3-TKD', 'IKZF1', 'RAD21')
genes.ML <- c('CALR', 'IDH1', 'ARID1A', 'ARID2', 'ASXL1', 'ATM', 'BCOR', 'BCORL1', 'BRAF', 
              'CBL', 'CUX1', 'DDX41', 'DNM2', 'EP300', 'ETV6', 'EZH2', 'FAT1', 
              'GNAS', 'JAK3', 'KMT2A', 'KMT2B', 'KMT2D', 'MPL', 'MYC', 'MYOM2', 
              'NF1', 'PHF6', 'PIK3CD', 'PRKDC', 'PRPF8', 'PTPN11', 'RELN', 
              'RUNX1', 'SETBP1', 'SF3B1', 'SH2B3', 'SMC1A', 'SMC3', 'SRSF2', 'STAG2', 'U2AF1')
genes.L <- c('JAK2')
genes <- c(genes.R, genes.MR, genes.ML, genes.L)

gene.cluster.data.raw <- read.csv(paste0(savepath, '/dirichlet cluster_mutation/', "label-sample-2021-11-21.csv"))
kar.cluster.data.raw <- read.csv(paste0(savepath, '/dirichlet cluster-karyotype/', "label-sample-2021-11-21.csv"))
gene.cluster.group <- data.frame(cluster = seq(1, 12, 1), cluster_group = sprintf("G%02d",c(2, 4, 6, 5, 9, 12, 1, 10, 7, 11, 8, 3)))
kar.cluster.group <- data.frame(cluster = seq(1, 7, 1), cluster_group = sprintf("K%02d",c(5, 7, 4, 6, 3, 1, 2)))
gene.group.data1 <- merge(gene.cluster.data.raw, gene.cluster.group, by.x = 'top1.cluster', by.y = 'cluster')
gene.group.data <- gene.group.data1[gene.group.data1$gene %in% genes, ]
kar.group.data <- merge(kar.cluster.data.raw, kar.cluster.group,  by.x = 'top1.cluster', by.y = 'cluster')

lambda.sd <- 0.5 
trait_ranks <- trait_ranks.raw[trait_ranks.raw$sdrank <= quantile(trait_ranks.raw$sdrank, lambda.sd) & 
                                 trait_ranks.raw$clinicalname %in% header.names[tra.num.start:tra.num.end, ], ]
trait_ranks <- trait_ranks[trait_ranks[, 'trait.cv']>0.5 & !is.na(trait_ranks[,'trait.cv']), ]
traits <- trait_ranks$clinicalname
gene.group <- c('G07', 'G09', 'G01', 'G06', 'G02', 'G08', 'G04', 'G12', 'G05', 'G03')
kar.group <- rev(c('K03', 'K06', 'K01', 'K05', 'K04', 'K07', 'K02'))

###gene group test###
L.VS.ML <- wilcox.test(gene_ranks.raw[gene_ranks.raw$genename %in% genes.L, ]$medianrank, 
                       gene_ranks.raw[gene_ranks.raw$genename %in% genes.ML, ]$medianrank)$p.value
ML.VS.MR <- wilcox.test(gene_ranks.raw[gene_ranks.raw$genename %in% genes.ML, ]$medianrank, 
                        gene_ranks.raw[gene_ranks.raw$genename %in% genes.MR, ]$medianrank)$p.value
MR.VS.R <- wilcox.test(gene_ranks.raw[gene_ranks.raw$genename %in% genes.MR, ]$medianrank, 
                       gene_ranks.raw[gene_ranks.raw$genename %in% genes.R, ]$medianrank)$p.value
group.median.rank.p <- data.frame(L.VS.ML = L.VS.ML, ML.VS.MR = ML.VS.MR, MR.VS.R = MR.VS.R)

sum.data <- data.frame(ID = cligene.data.MYE$ID, Diagnosis = cligene.data.MYE$Diagnosis,
                       R = rowSums(as.data.frame(cligene.data.MYE[, genes.R]), na.rm = T),
                       MR = rowSums(as.data.frame(cligene.data.MYE[, genes.MR]), na.rm = T),
                       ML = rowSums(as.data.frame(cligene.data.MYE[, genes.ML]), na.rm = T),
                       L = rowSums(as.data.frame(cligene.data.MYE[, genes.L]), na.rm = T),
                       other = rowSums(as.data.frame(cligene.data.MYE[, header.names$header.names[mut.num.start:mut.num.end][!header.names$header.names[mut.num.start:mut.num.end] %in% genes]]), na.rm = T))
if (! dir.exists(paste0(savepath, '/MYE/cor/'))) {
  dir.create(paste0(savepath, '/MYE/cor/'), recursive = TRUE)
}
write.csv(group.median.rank.p, paste0(savepath, '/MYE/cor/group_median_rank_p.csv'))
write.csv(sum.data, paste0(savepath, '/MYE/cor/sum_data.csv'))

###traits mutation MYE###
cor.gene.trait.data <- array(NA, dim = c(length(traits), length(gene.group)))
for (i in 1:length(traits)) {
  trait <- traits[i]
  trait.data <- cligene.data.MYE[, header.names == trait]
  for (j in 1:length(gene.group)) {
    gene <- gene.group.data[gene.group.data$cluster_group == gene.group[j], ]$gene
    gene.data <- rowSums(as.data.frame(cligene.data.MYE[, gene]), na.rm = T)
    cor.num <- try(cor(trait.data, gene.data))
    
    if(!is.error(cor.num)){
      if (!is.nan(cor.num)) {
        cor.gene.trait.data[i, j] <- cor.num
      }else{
        cor.gene.trait.data[i, j] <- 0
      }
    }else{
      cor.gene.trait.data[i, j] <- 0
    }
  }
  
}

colnames(cor.gene.trait.data) <- gene.group
rownames(cor.gene.trait.data) <- trait_ranks$clinicalname
cor.gene.trait.data[is.na(cor.gene.trait.data)] <- 0
write.csv(t(as.matrix(cor.gene.trait.data[, 1:length(gene.group)])), paste0(savepath, '/MYE/cor/cor_gene_trait.csv'))

my_palette <- colorRampPalette(c("green", "gray", "red"))
jpeg(filename = paste0(savepath, '/MYE/cor/cor_gene_trait.jpg'),
    width = 900 * 3, height = 3 * 600, res = 96 * 3)
gplots::heatmap.2(t(as.matrix(cor.gene.trait.data[, 1:length(gene.group)])), Rowv = NA, Colv = NA,
                  dendrogram = "none", density.info = "none",
                  labRow = gene.group[1:length(gene.group)], labCol = trait_ranks[, 2], trace = "none",
                  col = my_palette, keysize = 1, margins = c(18, 5), 
                  main = 'gene mutations vs. phenotypic features')
dev.off()

###by disease###
disease1 <- c('AML', 'MDS')
disease <- c()
for (k in 1:length(disease1)) {
  I.disease <- cligene.data.MYE$Diagnosis == disease1[k]
  disease[k] <- sub('/', '', disease1[k])
  cor.gene.trait.data <- array(NA, dim = c(length(traits), length(gene.group)))
  for (i in 1:length(traits)) {
    trait <- traits[i]
    trait.data <- cligene.data.MYE[I.disease, header.names == trait]
    for (j in 1:length(gene.group)) {
      gene <- gene.group.data[gene.group.data$cluster_group == gene.group[j], ]$gene
      gene.data <- rowSums(as.data.frame(cligene.data.MYE[I.disease, gene]), na.rm = T)
      cor.num <- try(cor(trait.data, gene.data))
      
      if(!is.error(cor.num)){
        if (!is.nan(cor.num)) {
          cor.gene.trait.data[i, j] <- cor.num
        }else{
          cor.gene.trait.data[i, j] <- 0
        }
      }else{
        cor.gene.trait.data[i, j] <- 0
      }
    }
  }
  
  colnames(cor.gene.trait.data) <- gene.group
  rownames(cor.gene.trait.data) <- trait_ranks$clinicalname
  cor.gene.trait.data[is.na(cor.gene.trait.data)] <- 0
  write.csv(t(as.matrix(cor.gene.trait.data[, 1:length(gene.group)])), paste0(savepath, '/MYE/cor/', disease[k], '_cor_gene_trait.csv'))
  
  if (disease[k] == 'AML') {
    row.split <- 20
    col.split <- 3
  }else{
    row.split <- 15
    col.split <- 8
  }
  
  r.p <- wilcox.test(cor.gene.trait.data[(row.split + 1):nrow(cor.gene.trait.data), 1:col.split], 
                     cor.gene.trait.data[(row.split + 1):nrow(cor.gene.trait.data), (col.split + 1):ncol(cor.gene.trait.data)])$p.value
  l.p <- wilcox.test(cor.gene.trait.data[1:row.split, 1:col.split],
                     cor.gene.trait.data[1:row.split, (col.split + 1):ncol(cor.gene.trait.data)])$p.value
  u.p <- wilcox.test(cor.gene.trait.data[1:row.split, 1:col.split],
                     cor.gene.trait.data[(row.split + 1):nrow(cor.gene.trait.data), 1:col.split])$p.value
  d.p <- wilcox.test(cor.gene.trait.data[1:row.split, (col.split + 1):ncol(cor.gene.trait.data)],
                     cor.gene.trait.data[(row.split + 1):nrow(cor.gene.trait.data), (col.split + 1):ncol(cor.gene.trait.data)])$p.value
  cortest <- data.frame(r.p = r.p, l.p = l.p, u.p = u.p, d.p = d.p)
  write.csv(cortest, paste0(savepath, '/MYE/cor/', disease[k],'_cor_pvalue.csv'))
  
  my_palette <- colorRampPalette(c("green", "gray", "red"))
  jpeg(filename = paste0(savepath, '/MYE/cor/', disease[k], '_cor_gene_trait.jpg'),
       width = 900 * 3, height = 3 * 600, res = 96 * 3)
  gplots::heatmap.2(t(as.matrix(cor.gene.trait.data[, 1:length(gene.group)])), Rowv = NA, Colv = NA,
                    dendrogram = "none", density.info = "none",
                    labRow = gene.group[1:length(gene.group)], labCol = trait_ranks[, 2], trace = "none",
                    col = my_palette, keysize = 1, margins = c(18, 5), 
                    main = paste0(disease[k], '\ngene mutations vs. phenotypic features'))
  dev.off()
}

###karyotype mutation MYE + LYM###
cor.gene.kar.data <- array(NA, dim = c(length(kar.group), length(gene.group)))
for (i in 1:length(kar.group)) {
  kar <- kar.group.data[kar.group.data$cluster_group == kar.group[i], ]$gene
  kar.data <- rowSums(as.data.frame(cligene.data[, kar]), na.rm = T)
  for (j in 1:length(gene.group)) {
    gene <- gene.group.data[gene.group.data$cluster_group == gene.group[j], ]$gene
    gene.data <- rowSums(as.data.frame(cligene.data[, gene]), na.rm = T)
    cor.num <- try(cor(kar.data, gene.data))
    
    if(!is.error(cor.num)){
      if (!is.nan(cor.num)) {
        cor.gene.kar.data[i, j] <- cor.num
      }else{
        cor.gene.kar.data[i, j] <- 0
      }
    }else{
      cor.gene.kar.data[i, j] <- 0
    }
  }
}

colnames(cor.gene.kar.data) <- gene.group
rownames(cor.gene.kar.data) <- kar.group
write.csv(t(as.matrix(cor.gene.kar.data[, 1:length(gene.group)])), paste0(savepath, '/MYE/cor/cor_kar_gene_total.csv'))

my_palette <- colorRampPalette(c("green", "gray", "red"))
jpeg(filename = paste0(savepath, '/MYE/cor/cor_kar_gene_total.jpg'),
    width = 900 * 3, height = 3 * 600, res = 96 * 3)
gplots::heatmap.2(t(as.matrix(cor.gene.kar.data[, 1:length(gene.group)])), Rowv = NA, Colv = NA,
                  dendrogram = "none", density.info = "none",
                  labRow = gene.group[1:length(gene.group)], labCol = kar.group, trace = "none",
                  col = my_palette, keysize = 1, margins = c(5, 5), 
                  main = 'gene mutations vs. cytogenetics')
dev.off()



print("-----------------------------Group testing of gene and traits ----------------")





####-------------------------------------------------------part 10: AML MDS response  -------------------------####
#*******************************************************************************
############## AML response ################ 

AML.data_subset <- cligene.data[cligene.data$Diagnosis == 'AML'&
                                     cligene.data$`M3 (for AML)` !=1 &
                                     !is.na(cligene.data$`MRD (for AML)`) &
                                     cligene.data$Age > 16 &
                                     cligene.data$treatment == "Pure conventional treatment",]


I <-   AML.data_subset$SRSF2 == 0 &
  AML.data_subset$SF3B1 == 0 &
  AML.data_subset$U2AF1 == 0 &
  AML.data_subset$ZRSR2 == 0 &
  AML.data_subset$ASXL1 == 0 &
  AML.data_subset$EZH2 == 0 &
  AML.data_subset$BCOR == 0 &
  AML.data_subset$STAG2 == 0 &
  AML.data_subset$`AML type`!='AML-MRC'

mrd.thres <- .1

AML.data_subset$mrd.class <- as.numeric(AML.data_subset$`MRD (for AML)` <= mrd.thres)
AML.data_subset$num.MR <-  rowSums(AML.data_subset[,genes.MR])
AML.data_subset$num.R <-  rowSums(AML.data_subset[,genes.R])
AML.data_subset$num.ML <-  rowSums(AML.data_subset[,genes.ML])
AML.data_subset$num.L <-  rowSums(as.data.frame(AML.data_subset[,genes.L]))
AML.data_subset$num.others <-  rowSums(AML.data_subset[,169:342])-
  rowSums(AML.data_subset[,genes.MR])-rowSums(AML.data_subset[,genes.R])-
  rowSums(AML.data_subset[,genes.ML])-rowSums(as.data.frame(AML.data_subset[,genes.L]))


AML.data_subset$MR.hi <- as.numeric(AML.data_subset$num.MR>quantile(AML.data_subset$num.MR, .6))
AML.data_subset$R.hi <- as.numeric(AML.data_subset$num.R>quantile(AML.data_subset$num.R, .6))
AML.data_subset$ML.hi <- as.numeric(AML.data_subset$num.ML>quantile(AML.data_subset$num.ML, .6))
AML.data_subset$L.hi <- as.numeric(AML.data_subset$num.L>quantile(AML.data_subset$num.L, .6))
AML.data_subset$others.hi <- as.numeric(AML.data_subset$num.others>quantile(AML.data_subset$num.others, .6))
AML.data_subset$MRR.hi<- 0
AML.data_subset$MRR.hi[AML.data_subset$num.MR+AML.data_subset$num.R > quantile(AML.data_subset$num.MR+AML.data_subset$num.R, .6)]<- 1
AML.data_subset$MLL.hi <- as.numeric(AML.data_subset$num.ML+AML.data_subset$num.L>quantile(AML.data_subset$num.ML+AML.data_subset$num.L, .6))
AML.data_subset$num.common <- rowSums(AML.data_subset[, c(169:342)[colSums(AML.data_subset[,169:342])>quantile(colSums(AML.data_subset[,169:342]), .75)]])
AML.data_subset$common.hi <- as.numeric(AML.data_subset$num.common>quantile(AML.data_subset$num.common, .6)) 

AML.data_subset$ELN.MRRhi <- paste(AML.data_subset$`2017 ELN risk score (for AML)`, AML.data_subset$MRR.hi, sep='')
AML.data_subset$ELN.group <- as.numeric(AML.data_subset$`2017 ELN risk score (for AML)`=='Favorable')*1 +
  as.numeric(AML.data_subset$`2017 ELN risk score (for AML)`=='Intermediate')*2 +
  as.numeric(AML.data_subset$`2017 ELN risk score (for AML)`=='Adverse')*3

AML.data_subset$ELN.enhanced.group <- as.numeric(AML.data_subset$ELN.MRRhi =='Favorable0')*1 +
  as.numeric(AML.data_subset$ELN.MRRhi=='Favorable1')*2 +
  as.numeric(AML.data_subset$ELN.MRRhi %in% c('Intermediate0'))*2 +
  as.numeric(AML.data_subset$ELN.MRRhi %in% c('Intermediate1','Adverse0','Adverse1'))*3


#Fig.7B


# all de novo AML
model <- glm(mrd.class ~ Age + `2017 ELN risk score (for AML)` + MRR.hi + MLL.hi + others.hi, 
             data=AML.data_subset, family='binomial')
summary(model)

# "bona fide" all de novo AML
model <- glm(mrd.class ~ Age + `2017 ELN risk score (for AML)` + MRR.hi + MLL.hi + others.hi, 
             data=AML.data_subset[I,], family='binomial')
summary(model)


#Fig.7D

# all de novo AML
roc(AML.data_subset$mrd.class, AML.data_subset$ELN.group,direction = '>')
roc(AML.data_subset$mrd.class,AML.data_subset$ELN.enhanced.group,direction = '>')

# "bona fide" all de novo AML
roc(AML.data_subset$mrd.class[I], AML.data_subset$ELN.group[I],direction = '>')
roc(AML.data_subset$mrd.class[I],AML.data_subset$ELN.enhanced.group[I],direction = '>')




# Fig.7E 

# ELN

mat <- c()
mat <- rbind(mat, c(sum(AML.data_subset$mrd.class[AML.data_subset$ELN.group==1 & I]==1),
                    sum(AML.data_subset$mrd.class[AML.data_subset$ELN.group==1 & I]==0)))
mat <- rbind(mat, c(sum(AML.data_subset$mrd.class[AML.data_subset$ELN.group==2 & I]==1),
                    sum(AML.data_subset$mrd.class[AML.data_subset$ELN.group==2 & I]==0)))
mat <- rbind(mat, c(sum(AML.data_subset$mrd.class[AML.data_subset$ELN.group==3 & I]==1), 
                    sum(AML.data_subset$mrd.class[AML.data_subset$ELN.group==3 & I]==0)))
prop.table(mat,1)
fisher.test(mat[1:2,])
fisher.test(mat[2:3,])

# ELN enhanced

mat <- c()
mat <- rbind(mat, c(sum(AML.data_subset$mrd.class[AML.data_subset$ELN.enhanced.group==1 & I]==1),
                    sum(AML.data_subset$mrd.class[AML.data_subset$ELN.enhanced.group==1 & I]==0)))
mat <- rbind(mat, c(sum(AML.data_subset$mrd.class[AML.data_subset$ELN.enhanced.group==2 & I]==1),
                    sum(AML.data_subset$mrd.class[AML.data_subset$ELN.enhanced.group==2 & I]==0)))
mat <- rbind(mat, c(sum(AML.data_subset$mrd.class[AML.data_subset$ELN.enhanced.group==3 & I]==1), 
                    sum(AML.data_subset$mrd.class[AML.data_subset$ELN.enhanced.group==3 & I]==0)))
prop.table(mat,1)
fisher.test(mat[1:2,])
fisher.test(mat[2:3,])



############## MDS response ################ 

MDS.data <- cligene.data[cligene.data$Diagnosis == 'MDS',]

I <- MDS.data$`2016 WHO category (for MDS)` == 'EB' & 
  MDS.data$Age > 16 

MDS.data_subset <- MDS.data[I,]

MDS.data_subset$num.MR <-  rowSums(MDS.data_subset[,genes.MR])
MDS.data_subset$num.R <- rowSums(MDS.data_subset[,genes.R])
MDS.data_subset$num.ML <- rowSums(MDS.data_subset[,genes.ML])
MDS.data_subset$num.L <- rowSums(as.data.frame(MDS.data_subset[,genes.L]))
MDS.data_subset$num.others <-  rowSums(MDS.data_subset[,169:342])-
  rowSums(MDS.data_subset[,genes.MR])-rowSums(MDS.data_subset[,genes.R])-
  rowSums(MDS.data_subset[,genes.ML])-rowSums(as.data.frame(MDS.data_subset[,genes.L]))


MDS.data_subset$MR.hi <- as.numeric(MDS.data_subset$num.MR>quantile(MDS.data_subset$num.MR, .6))
MDS.data_subset$R.hi <- as.numeric(MDS.data_subset$num.R>quantile(MDS.data_subset$num.R, .6))
MDS.data_subset$ML.hi <- as.numeric(MDS.data_subset$num.ML>quantile(MDS.data_subset$num.ML, .6))
MDS.data_subset$L.hi <- as.numeric(MDS.data_subset$num.L>quantile(MDS.data_subset$num.L, .6))
MDS.data_subset$others.hi <- as.numeric(MDS.data_subset$num.others>quantile(MDS.data_subset$num.others, .6))
MDS.data_subset$MRR.hi <- as.numeric(MDS.data_subset$num.MR+MDS.data_subset$num.R>quantile(MDS.data_subset$num.MR+MDS.data_subset$num.R, .6))
MDS.data_subset$MLL.hi <- as.numeric(MDS.data_subset$num.ML+MDS.data_subset$num.L>quantile(MDS.data_subset$num.ML+MDS.data_subset$num.L, .6))
MDS.data_subset$num.common <- rowSums(MDS.data_subset[, c(169:342)[colSums(MDS.data_subset[,169:342])>quantile(colSums(MDS.data[,169:342]), .75)]]) 
MDS.data_subset$common.hi <- as.numeric(MDS.data_subset$num.common>quantile(MDS.data_subset$num.common, .6)) 

#Fig.8A 

model <- glm(`CR after 6 cycles of chemotherapy (for MDS)`~ Age + `IPSS-R score (for MDS)` + MRR.hi + MLL.hi + others.hi,
             data=MDS.data_subset, family='binomial')
summary(model)

#Fig.8B

mat <- c()
I <- MDS.data_subset$MLL.hi==1
mat <- rbind(mat, c(sum(MDS.data_subset$`CR after 6 cycles of chemotherapy (for MDS)`[I]==1, na.rm=T), sum(MDS.data_subset$`CR after 6 cycles of chemotherapy (for MDS)`[I]==0, na.rm=T)))
I <- MDS.data_subset$MLL.hi==0
mat <- rbind(mat, c(sum(MDS.data_subset$`CR after 6 cycles of chemotherapy (for MDS)`[I]==1, na.rm=T), sum(MDS.data_subset$`CR after 6 cycles of chemotherapy (for MDS)`[I]==0, na.rm=T)))

fisher.test(mat,alternative="two.sided")

#Fig.8C

roc(MDS.data_subset$`CR after 6 cycles of chemotherapy (for MDS)`[!is.na(MDS.data_subset$`CR after 6 cycles of chemotherapy (for MDS)`)],
    MDS.data_subset$`IPSS-R score (for MDS)`[!is.na(MDS.data_subset$`CR after 6 cycles of chemotherapy (for MDS)`)],
    direction = '>')
roc(MDS.data_subset$`CR after 6 cycles of chemotherapy (for MDS)`[!is.na(MDS.data_subset$`CR after 6 cycles of chemotherapy (for MDS)`)],
    MDS.data_subset$MLL.hi[!is.na(MDS.data_subset$`CR after 6 cycles of chemotherapy (for MDS)`)],
    direction = '>')


print("---------------------AML-MDS response done !!----------")
#*******************************************************************************

###end####
print('-----------------ALL FINISHED!!-------------')
print('-----All eclipse time:')
print(Sys.time() - t1)





