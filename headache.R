## Danni Shi
## R codes for Section 3.2 of the Master's Thesis
## 'Behavioral Carry-over Effects & Power Considerations'
## dannis5@uw.edu | dshidanni@gmail.com

rm(list=ls())
# devtools::install_github("mbannick/RobinCar")
library(RobinCar)

# setwd("C:/Users/danni/Dropbox/Thesis")
dat <- read.table("headache.txt",header=TRUE)
# head(dat)

##### Sorting the data #####
# ABBA
AB.dat <- dat[dat$TrmtSeq==2,] # AB
n.AB <- nrow(AB.dat)/2
BA.dat <- dat[dat$TrmtSeq==1,] # BA
n.BA <- nrow(BA.dat)/2
AB.A <- AB.dat$Response[AB.dat$Treatment=="A"]
AB.B <- AB.dat$Response[AB.dat$Treatment=="B"]
BA.B <- BA.dat$Response[BA.dat$Treatment=="B"]
BA.A <- BA.dat$Response[BA.dat$Treatment=="A"]

ABBA.ID <- c(AB.dat$ID[seq(from=1,to=nrow(AB.dat),by=2)],
             BA.dat$ID[seq(from=1,to=nrow(BA.dat),by=2)])
ABBA.Center <- c(AB.dat$Center[seq(from=1,to=nrow(AB.dat),by=2)],
                 BA.dat$Center[seq(from=1,to=nrow(BA.dat),by=2)])
ABBA.trmtseq <- c(AB.dat$TrmtSeq[seq(from=1,to=nrow(AB.dat),by=2)],
                  BA.dat$TrmtSeq[seq(from=1,to=nrow(BA.dat),by=2)])
ABBA.trmtseq <- as.numeric(ABBA.trmtseq==2) # 1:AB; 0:BA
ABBA.Y1 <- c(AB.dat$Response[seq(from=1,to=nrow(AB.dat),by=2)],
             BA.dat$Response[seq(from=1,to=nrow(BA.dat),by=2)])
ABBA.Y2 <- c(AB.dat$Response[seq(from=2,to=nrow(AB.dat),by=2)],
             BA.dat$Response[seq(from=2,to=nrow(BA.dat),by=2)])
ABBA.deltaY <- ABBA.Y1-ABBA.Y2 # Y1-Y2
ABBA.dat <- data.frame(ID=ABBA.ID, center=ABBA.Center, trmtseq=ABBA.trmtseq,
                       Y1=ABBA.Y1, Y2=ABBA.Y2, deltaY=ABBA.deltaY)
ABBA.dat$trmt <- c(rep("AB",n.AB),rep("BA",n.BA))



# BPPB
BP.dat <- dat[dat$TrmtSeq==4,] # BP
PB.dat <- dat[dat$TrmtSeq==3,] # PB
PB.P <- with(PB.dat, Response[Treatment=="P"])
PB.B <- with(PB.dat, Response[Treatment=="B"])
n.PB <- nrow(PB.dat)/2
BP.P <- with(BP.dat, Response[Treatment=="P"])
BP.B <- with(BP.dat, Response[Treatment=="B"])
n.BP <- nrow(BP.dat)/2

BPPB.ID <- c(BP.dat$ID[seq(from=1,to=nrow(BP.dat),by=2)],
             PB.dat$ID[seq(from=1,to=nrow(PB.dat),by=2)])
BPPB.Center <- c(BP.dat$Center[seq(from=1,to=nrow(BP.dat),by=2)],
                 PB.dat$Center[seq(from=1,to=nrow(PB.dat),by=2)])
BPPB.trmtseq <- c(BP.dat$TrmtSeq[seq(from=1,to=nrow(BP.dat),by=2)],
                  PB.dat$TrmtSeq[seq(from=1,to=nrow(PB.dat),by=2)])
BPPB.trmtseq <- as.numeric(BPPB.trmtseq==4) # 1:BP; 0:PB
BPPB.Y1 <- c(BP.dat$Response[seq(from=1,to=nrow(BP.dat),by=2)],
             PB.dat$Response[seq(from=1,to=nrow(PB.dat),by=2)])
BPPB.Y2 <- c(BP.dat$Response[seq(from=2,to=nrow(BP.dat),by=2)],
             PB.dat$Response[seq(from=2,to=nrow(PB.dat),by=2)])
BPPB.deltaY <- BPPB.Y1-BPPB.Y2 # Y1-Y2
BPPB.dat <- data.frame(ID=BPPB.ID, center=BPPB.Center, trmtseq=BPPB.trmtseq,
                       Y1=BPPB.Y1, Y2=BPPB.Y2, deltaY=BPPB.deltaY)



# APPA
AP.dat <- dat[dat$TrmtSeq==6,] # AP
PA.dat <- dat[dat$TrmtSeq==5,] # PA
PA.P <- with(PA.dat, Response[Treatment=="P"])
PA.A <- with(PA.dat, Response[Treatment=="A"])
n.PA <- nrow(PA.dat)/2
AP.P <- with(AP.dat, Response[Treatment=="P"])
AP.A <- with(AP.dat, Response[Treatment=="A"])
n.AP <- nrow(AP.dat)/2

APPA.ID <- c(AP.dat$ID[seq(from=1,to=nrow(AP.dat),by=2)],
             PA.dat$ID[seq(from=1,to=nrow(PA.dat),by=2)])
APPA.Center <- c(AP.dat$Center[seq(from=1,to=nrow(AP.dat),by=2)],
                 PA.dat$Center[seq(from=1,to=nrow(PA.dat),by=2)])
APPA.trmtseq <- c(AP.dat$TrmtSeq[seq(from=1,to=nrow(AP.dat),by=2)],
                  PA.dat$TrmtSeq[seq(from=1,to=nrow(PA.dat),by=2)])
APPA.trmtseq <- as.numeric(APPA.trmtseq==6) # 1:AP; 0:PA
APPA.Y1 <- c(AP.dat$Response[seq(from=1,to=nrow(AP.dat),by=2)],
             PA.dat$Response[seq(from=1,to=nrow(PA.dat),by=2)])
APPA.Y2 <- c(AP.dat$Response[seq(from=2,to=nrow(AP.dat),by=2)],
             PA.dat$Response[seq(from=2,to=nrow(PA.dat),by=2)])
APPA.deltaY <- APPA.Y1-APPA.Y2 # Y1-Y2
APPA.dat <- data.frame(ID=APPA.ID, center=APPA.Center, trmtseq=APPA.trmtseq,
                       Y1=APPA.Y1, Y2=APPA.Y2, deltaY=APPA.deltaY)


##### Table One #####
library(tableone)
ABBA.dat <- ABBA.dat[,c(1:6)]
ABBA.dat$study <- "ABBA"
APPA.dat$study <- "APPA"
BPPB.dat$study <- "BPPB"
dat.all <- rbind(ABBA.dat,APPA.dat,BPPB.dat)
dat.all$center <- as.factor(dat.all$center)
dat.all$trmtseq <- as.factor(dat.all$trmtseq)
CreateTableOne(data=dat.all,
               vars=c("center"),
               strata=c("study","trmtseq"),
               includeNA=TRUE,addOverall=TRUE,test=FALSE)



##### Descriptive statistics #####
headache.descriptive <- 
  data.frame(type=c("AB","BA","BP","PB","AP","PA"),
             n=c(n.AB,n.BA,n.BP,n.PB,n.AP,n.PA),
             period.1.mean=c(mean(AB.A),mean(BA.B),mean(BP.B),mean(PB.P),
                             mean(AP.A),mean(PA.P)),
             period.1.sd=c(sd(AB.A)/sqrt(n.AB),sd(BA.B)/sqrt(n.BA),
                           sd(BP.B)/sqrt(n.BP),sd(PB.P)/sqrt(n.PB),
                           sd(AP.A)/sqrt(n.AP),sd(PA.P)/sqrt(n.PA)),
             period.2.mean=c(mean(AB.B),mean(BA.A),mean(BP.P),mean(PB.B),
                             mean(AP.P),mean(PA.A)),
             period.2.sd=c(sd(AB.B)/sqrt(n.AB),sd(BA.A)/sqrt(n.BA),
                           sd(BP.P)/sqrt(n.BP),sd(PB.B)/sqrt(n.PB),
                           sd(AP.P)/sqrt(n.AP),sd(PA.A)/sqrt(n.PA)),
             diff.mean=c(mean(AB.A-AB.B),mean(BA.A-BA.B),
                         mean(BP.B-BP.P),mean(PB.B-PB.P),
                         mean(AP.A-AP.P),mean(PA.A-PA.P)),
             diff.sd=c(sd(AB.A-AB.B)/sqrt(n.AB),sd(BA.A-BA.B)/sqrt(n.BA),
                       sd(BP.B-BP.P)/sqrt(n.BP),sd(PB.B-PB.P)/sqrt(n.PB),
                       sd(AP.A-AP.P)/sqrt(n.AP),sd(PA.A-PA.P)/sqrt(n.PA)),
             rho=c(cor(AB.A,AB.B),cor(BA.B,BA.A),cor(BP.B,BP.P),
                   cor(PB.P,PB.B),cor(AP.A,AP.P),cor(PA.P,PA.A)))
headache.descriptive[,c(3:9)] <- round(headache.descriptive[,c(3:9)],3)



##### Carry-over effects of A & B #####
# compare treatment effect from AB & PB sequences
mean(AB.B)+c(-1,1)*qnorm(.975)*sd(AB.B)/sqrt(length(AB.B))
mean(PB.B)+c(-1,1)*qnorm(.975)*sd(PB.B)/sqrt(length(PB.B))

# compare treatment effect from BA & PA sequences
mean(BA.A)+c(-1,1)*qnorm(.975)*sd(BA.A)/sqrt(length(BA.A))
mean(PA.A)+c(-1,1)*qnorm(.975)*sd(PA.A)/sqrt(length(PA.A))




##### 4 estimators #####
library(RobinCar)
# parallel
temp.pr <- robincar_linear(df=ABBA.dat, response_col="Y1", treat_col="trmtseq",
                           car_scheme="simple", adj_method="ANOVA", contrast_h="diff")$contrast$result
pr.est <- temp.pr$estimate # 0.6157199 
pr.se <- temp.pr$se # 0.4555336
pr.est+c(-1,1)*qnorm(.975)*pr.se # -0.2771096  1.5085494

# paralle, adjusted
ABBA.dat$center <- as.factor(ABBA.dat$center)
temp.pradj <- robincar_linear(df=ABBA.dat, response="Y1", treat_col="trmtseq",
                              covariate_cols=c("center"),
                              car_scheme="simple", adj_method="ANHECOVA",
                              contrast_h="diff")$contrast$result
pradj.est <- temp.pradj$estimate # 0.6266414
pradj.se <- temp.pradj$se # 0.4475047
pradj.est+c(-1,1)*qnorm(.975)*pradj.se # -0.2504516  1.5037344

# crossover, unadjusted
temp.cr <- robincar_linear(df=ABBA.dat, response_col="deltaY", treat_col="trmtseq",
                           car_scheme="simple", adj_method="ANOVA", contrast_h="diff")$contrast$result
cr.est <- .5*temp.cr$estimate # 1.127141 
cr.se <- .5*temp.cr$se # 0.2732061
cr.est+c(-1,1)*qnorm(.975)*cr.se # 1.718807 2.789755

# crossover, adjusted
temp.cradj <- robincar_linear(df=ABBA.dat,response_col="deltaY",treat_col="trmtseq",
                              covariate_cols=c("center"),
                              car_scheme="simple",adj_method="ANHECOVA",contrast_h="diff")$contrast$result
cradj.est <- temp.cradj$estimate # 1.13938 
cradj.se <- .5*temp.cradj$se # 0.2645545
cradj.est+c(-1,1)*qnorm(.975)*cradj.se # 1.760243 2.797278
