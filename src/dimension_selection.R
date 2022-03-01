#This script is for evaluating the PCA projections and assistant the selection of dimensions.
#Calling R packages
library(readxl)
library(car)
library(ggplot2)
library(xlsx)
library(stringr)
library(reshape2)
library(car)
library(sp)
library(readxl)
library(optparse)

#Command line parameter parsing
option_list <- list(
  make_option(c("-i","--indir"),type="character",help="2d_full_data"),
  make_option(c("-o","--outdir"),type="character",help="output dir")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

root= "THAMP-data-and-code/src/"
setwd<-(root)
opt$indir<-"../output/2d_full_data"
opt$outdir<-"../output/dim_sel"

indir <- opt$indir
outdir<- opt$outdir

#parameters check
if (is.null(opt$indir) | is.null(opt$outdir)){
  print_help(opt_parser)
  stop("missing parameters, please check carefully !!", call.=FALSE)
}
# create output directory
if (! dir.exists(opt$outdir)) {
  dir.create(opt$outdir)
}


print("----------step1 Calculate the ratio of eigenvalues---------------")
#--------step1 Calculate the ratio of eigenvalues
a = list.files(opt$indir)
vec1<- c()
vec2<- c()
x<- c()
for (i in 1:length(a)){
  key = str_split(a[i],"_")[[1]][3]
  df<- as.data.frame(read_excel(path = paste0(opt$indir,"/",a[i]),sheet = 1))
  dd<- df[c("X","Y","Diagnosis")]
  dd1<- dd[dd$Diagnosis %in% c('AML','MDS','MPN','MDS/MPN'),]
  dd2<- dd[dd$Diagnosis %in% c('ALL','CLL','MM'),]

  dd1.pca<- prcomp(dd1[1:2],scale=F)
  dd2.pca<- prcomp(dd2[1:2],scale=F)
  #ratio of eigenvalue from PC1 to PC2
  vec1<- c(vec1,dd1.pca$sdev[1]^2/dd1.pca$sdev[2]^2)
  vec2<- c(vec2,dd2.pca$sdev[1]^2/dd2.pca$sdev[2]^2)
  x<- c(x, as.numeric(key))
}
#Aggregate and output calculation result
Eigen.ratio<- data.frame(lymphoid = vec2, myeloid = vec1, dimension = x)
#result2<- melt(result, measure.vars = c("lymphoid","myeloid"))



print("-------------step2 calculate weighted average of purity-------")
#--------step2 calculate weighted average of purity
MS<- c('AML','MDS','MPN','MDS/MPN')
evaluation_reduc<- function(df){
  dd<- df[c("X","Y","Diagnosis","ID")]
  colnames(dd)<- c("X","Y","Diagnosis","ID")
  dd<- dd[dd$Diagnosis %in% c('AML','MDS','MPN','MDS/MPN'),]
  #get confidence ellipse for each groups
  a<- dataEllipse(dd$X, dd$Y,levels = 0.8, groups = as.factor(dd$Diagnosis))
  eva<- function(ex, g){
    typ_in <- dd[ex==1,]$Diagnosis
    pur.rate<- sum((typ_in == g))/length(typ_in)
    cov.rate<- sum((typ_in == g))/sum(dd$Diagnosis ==g)
    return(list("pur"=pur.rate, "cov"=cov.rate))
  }
  #calculate the purity and coverage of each group
  rate.AML<- eva(point.in.polygon(dd$X, dd$Y, a$AML[,"x"], a$AML[,"y"]), "AML")
  rate.MDS<- eva(point.in.polygon(dd$X, dd$Y, a$MDS[,"x"], a$MDS[,"y"]), "MDS")
  rate.MPN<- eva(point.in.polygon(dd$X, dd$Y, a$MPN[,"x"], a$MPN[,"y"]), "MPN")
  rate.MDS_MPN<- eva(point.in.polygon(dd$X, dd$Y, a$`MDS/MPN`[,"x"], a$`MDS/MPN`[,"y"]), "MDS/MPN")
  result<- list("Diagnosis"=c('AML','MDS','MPN','MDS/MPN'),
                "pur"= c(rate.AML$pur, rate.MDS$pur, rate.MPN$pur, rate.MDS_MPN$pur),
                "cov"= c(rate.AML$cov, rate.MDS$cov, rate.MPN$cov, rate.MDS_MPN$cov))
  return(result)
}


#Merging and Sorting result
flist<- list.files(opt$indir)
all.pur<- matrix(data = 1:length(flist)*5,nrow = length(flist), ncol = 5)
all.cov<-  matrix(data = 1:length(flist)*5,nrow = length(flist), ncol = 5)
for (i in 1:length(flist) ) {
  df<- as.data.frame(read_excel(path = paste0(opt$indir,"/",flist[i]), sheet = 1))
  e<- evaluation_reduc(df)
  key = str_split(flist[i],"_")[[1]][3]
  all.pur[i,] <- c(key,e$pur)
  all.cov[i,] <- c(key,e$cov)
}


print("------------step3 plot and output results--------")
#Export to local files
colnames(all.pur)<- c("dim",'AML','MDS','MPN','MDS/MPN')
colnames(all.cov)<- c("dim",'AML','MDS','MPN','MDS/MPN')
# write.csv(as.data.frame(all.pur), file = paste0(opt$outdir,"/all.purity_of_all_dimensions.csv"))
# write.csv(as.data.frame(all.cov), file = paste0(opt$outdir,"/all.coverage_of_all_dimensions.csv"))

#plot scatter diagram of dimension evaluation.
#Then we could select the proper dimension in the upper right area
dd<- df[df$Diagnosis %in% c('AML','MDS','MPN','MDS/MPN'),]
all.pur <- as.data.frame(all.pur)
stat <- table(dd$Diagnosis)
#head(all.pur)

all.pur$weighted_average<- (as.numeric(all.pur$AML) * stat["AML"] + as.numeric(all.pur$MDS) * stat["MDS"] + as.numeric(all.pur$MPN)*stat["MPN"] + as.numeric(all.pur$`MDS/MPN`)*stat["MDS/MPN"])/sum(stat)
all.pur$eigenvalue_rate <- Eigen.ratio$myeloid
#plot and select
png(paste0(opt$outdir,"/dimension_evaluation_",Sys.Date(),".png"),width = 800,height = 600,res=100)

plot(all.pur$weighted_average, all.pur$eigenvalue_rate, type="n", xlab="Weighted average purity",ylab="Ratio of eigenvalues")
text(all.pur$weighted_average, all.pur$eigenvalue_rate, labels = all.pur$dim)
dev.off()
write.xlsx(all.pur, file = paste0(opt$outdir,"/dimension_evaluation_",Sys.Date(),".xlsx"))
print("----------dimension selection done-------------")
