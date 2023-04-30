## Danni Shi
## R codes for Section 3.1 of the Master's Thesis
## 'Behavioral Carry-over Effects & Power Considerations'
## dannis5@uw.edu | dshidanni@gmail.com

rm(list=ls())
# devtools::install_github("mbannick/RobinCar")
library(RobinCar)
library(ggplot2)
library(reshape2)


## Set up 
N <- 2500
alpha <- 0.025

thetas <- seq(0,0.5,by=0.05)
lambdas <- c(-0.1,0,0.1,0.3)
bs <- c(0,1/3)
tau <- -0.05



##### formula powers #####
formula.power.pr <- function(n,theta) {
  sigma.pr <- sqrt(13)
  return(pnorm(qnorm(alpha)+sqrt(n)*theta/sigma.pr))
}
formula.power.pr(550,0.5) # 0.9018654
n <- 550
# set n=550 so that power_pr is around 90% when theta=0.5


formula.power.pradj <- function(n,theta) {
  sigma.pradj <- 2
  return(pnorm(qnorm(alpha)+sqrt(n)*theta/sigma.pradj))
}


formula.power.cr <- function(n,theta,lambda,b) {
  sigma.cr <- sqrt((1-b)^2+2)
  return(pnorm(qnorm(alpha)+sqrt(n)*(theta-lambda)/sigma.cr))
}
# power_cr,unadj is the only one that can be affected by b.


formula.power.cradj <- function(n,theta,lambda) {
  sigma.cradj <- sqrt(2)
  return(pnorm(qnorm(alpha)+sqrt(n)*(theta-lambda)/sigma.cradj))
}

# generate and store formula powers
form.result <- data.frame(theta=rep(thetas,2*4),
                          b=rep(bs,each=length(thetas)*4),
                          lambda=rep(rep(lambdas,each=length(thetas)),2),
                          pr=NA,pradj=NA,cr=NA,cradj=NA,type="formula")
form.result$pr <- rep(formula.power.pr(n,thetas),8)
form.result$pradj <- rep(formula.power.pradj(n,thetas),8)

form.cr.res <- c()
for (i in 1:2) {
  b.tmp <- bs[i]
  for (j in 1:4) {
    lambda.tmp <- lambdas[j]
    form.cr.res <- c(form.cr.res, 
                     formula.power.cr(n,thetas,lambda.tmp,b.tmp))
  }
}
form.result$cr <- form.cr.res

form.cradj.res <- c()
for (i in 1:4) {
  lambda.tmp <- lambdas[j]
  form.cradj.res <- c(form.cradj.res,
                      formula.power.cradj(n,thetas,lambda.tmp))
}
form.result$cradj <- form.cradj.res
write.csv(form.result,"formu_dat.csv")




##### data generating & simulation power #####
simu.power <- function(n,N,thetas,lambda,tau,b) {
  simu.res <- matrix(data=0, ncol=8, nrow=length(thetas))
  colnames(simu.res) <- c("theta","b","lambda","pr","pradj","cr","cradj","type")
  simu.res <- as.data.frame(simu.res)
  simu.res$b <- b
  simu.res$lambda <-lambda
  simu.res$type <- "simulation"
  
  for (i in 1:length(thetas)) {
    theta <- thetas[i]
    simu.res[i,1] <- theta
    
    for (j in 1:N) {
      X1 <- rnorm(n,0,1)
      X2 <- rnorm(n,0,1)
      X3 <- rbinom(n,1,.5)
      
      Ai <- rbinom(n,1,.5)
      Yi1_0 <- X1+X2+X3+rnorm(n,0,1)
      Yi1_1 <- theta+X1+X2+X3+rnorm(n,0,1)
      Yi2_10 <- tau+lambda+X1+b*X2+X3+rnorm(n,0,1)
      Yi2_01 <- tau+theta-lambda+X1+b*X2+X3+rnorm(n,0,1)
      
      dat.simu <- data.frame(Ai=Ai, X1=X1, X2=X2, X3=X3)
      dat.simu$Y1 <- Yi1_1*Ai+Yi1_0*(1-Ai)
      dat.simu$Y2 <- Yi2_10*Ai+Yi2_01*(1-Ai)
      dat.simu$delta.Y <- with(dat.simu, Y1-Y2)
      
      # parallel
      temp.pr <- robincar_linear(df=dat.simu, response_col="Y1", treat_col="Ai",
                                 car_scheme="simple", adj_method="ANOVA", contrast_h="diff")$contrast$result
      simu.res[i,4] <- simu.res[i,4]+1*(temp.pr$estimate/temp.pr$se>qnorm(1-alpha))
      
      # parallel, adjusted
      temp.pradj <- robincar_linear(df=dat.simu, response="Y1", treat_col="Ai",
                                    covariate_cols=c("X1","X2","X3"),
                                    car_scheme="simple", adj_method="ANHECOVA",
                                    contrast_h="diff")$contrast$result
      simu.res[i,5] <- simu.res[i,5]+1*(temp.pradj$estimate/temp.pradj$se>qnorm(1-alpha))
      
      # crossover, unadjusted
      temp.cr <- robincar_linear(df=dat.simu, response_col="delta.Y", treat_col="Ai",
                                 car_scheme="simple", adj_method="ANOVA", contrast_h="diff")$contrast$result
      simu.res[i,6] <- simu.res[i,6]+1*(temp.cr$estimate/temp.cr$se>qnorm(1-alpha))
      
      # crossover, adjusted
      temp.cradj <- robincar_linear(df=dat.simu,response_col="delta.Y",treat_col="Ai",
                                    covariate_cols=c("X1","X2","X3"),
                                    car_scheme="simple",adj_method="ANHECOVA",contrast_h="diff")$contrast$result
      simu.res[i,7] <- simu.res[i,7]+1*(temp.cradj$estimate/temp.cradj$se>qnorm(1-alpha))
    }
  }
  
  simu.res[,4:7] <- simu.res[,4:7]/N
  return(simu.res)
}

simu.result <- data.frame(theta=NA,b=NA,lambda=NA,pr=NA,pradj=NA,cr=NA,cradj=NA,type="simulation")
for (i in 1:2) {
  b.tmp <- bs[i]
  for (j in 1:4) {
    res.tmp <- simu.power(n,N,thetas,lambdas[j],tau,b.tmp)
    simu.result <- rbind(simu.result,res.tmp)
  }
}
simu.result <- simu.result[-1,]
write.csv(simu.result, "simu_dat.csv")




##### Data analysis #####
# setwd("C:/Users/danni/Dropbox/Thesis")
formu.dat <- read.csv("formu_dat.csv")
formu.dat <- formu.dat[,-1]
simu.dat <- read.csv("simu_dat.csv")
simu.dat <- simu.dat[,-1]

# expanded data frame for formula powers
df.long <- melt(formu.dat[,-8],id.vars=c("b","lambda","theta"))
df.long$b <- as.factor(df.long$b)
df.long$lambda <- as.factor(df.long$lambda)
names(df.long)[4] <- "Test"
df.subset <- subset(df.long,Test=="cr")
df.subset <- rbind(df.subset,subset(df.long,Test=="pr" & b==0))
df.subset <- rbind(df.subset,subset(df.long,Test=="pradj" & b==0))
df.subset <- rbind(df.subset,subset(df.long,Test=="cradj" & b==0))

levels(df.subset$b) <- c("0","1/3")
levels(df.subset$Test) <- c("Parallel","Parallel (adj)",
                            "Crossover (unadj)","Crossover (adj)")
levels(df.subset$lambda) <- c(expression(paste(lambda," = - 0.1")),
                              expression(paste(lambda," = 0")), 
                              expression(paste(lambda," = 0.1")),
                              expression(paste(lambda," = 0.3")))


# expanded data frame for simulation powers
df.long.sim <- melt(simu.dat[,-8],id.vars=c("b","lambda","theta"))
df.long.sim$b <- as.factor(df.long.sim$b)
df.long.sim$lambda <- as.factor(df.long.sim$lambda)
names(df.long.sim)[4] <- "Test"
df.long.sim$type <- "Simulation"
levels(df.long.sim$b) <- c("0","1/3")
levels(df.long.sim$Test) <- c("Parallel","Parallel (adj)",
                              "Crossover (unadj)","Crossover (adj)")
levels(df.long.sim$lambda) <- c(expression(paste(lambda," = - 0.1")),
                                expression(paste(lambda," = 0")), 
                                expression(paste(lambda," = 0.1")),
                                expression(paste(lambda," = 0.3")))


# formula powers when b=0
df.subset0 <- subset(df.subset,b==0)
df.subset0$type <- "Formula"
# simulation powers when b=0
df.long.sim0 <- subset(df.long.sim,b=="0")
# combine both
df.long.both0 <- rbind(df.long.sim0,df.subset0)
names(df.long.both0)[6] <- "Type"
levels(df.long$b) <- c("0","1/3")


# formula powers when b=1/3
df.subset13 <- subset(df.long,Test=="cr")
df.subset13 <- subset(df.subset13,b=="1/3")
df.subset13 <- rbind(df.subset13,subset(df.long,Test=="pr" & b=="1/3"))
df.subset13 <- rbind(df.subset13,subset(df.long,Test=="pradj" & b=="1/3"))
df.subset13 <- rbind(df.subset13,subset(df.long,Test=="cradj" & b=="1/3"))
df.subset13$Type <- "Formula"
levels(df.subset13$Test) <- c("Parallel","Parallel (adj)",
                              "Crossover (unadj)","Crossover (adj)")
levels(df.subset13$lambda) <- c(expression(paste(lambda," = - 0.1")),
                                expression(paste(lambda," = 0")), 
                                expression(paste(lambda," = 0.1")),
                                expression(paste(lambda," = 0.3")))
# simulation powers when b=1/3
df.long.sim13 <- subset(df.long.sim,b=="1/3")
names(df.long.sim13)[6] <- "Type"
# combine both
df.long.both13 <- rbind(df.long.sim13,df.subset13)
df.long.both13 <- df.long.both13[,-1]
df.long.both13$b <- "1/3"


## ggplot for b=0
ggplot(df.long.both0, aes(theta, value))+geom_line(aes(color=Type,linetype=Test),linewidth=1)+
  theme_bw()+theme(plot.title=element_text(hjust=0.5,size=12))+
  ggtitle("(a) b = 0")+
  xlab(expression(theta))+ ylab("power")+
  scale_linetype_manual(values=c("dotted","dashed","solid","dotdash"),
                        labels=c(expression(T[pr]),expression(T["pr,adj"]),
                                 expression(T[cr]),expression(T["cr,adj"])))+
  scale_color_manual(values=c("black","plum3"))+
  facet_grid(.~lambda,labeller = label_parsed)+
  scale_y_continuous(breaks=c(0.00,0.025,0.2,0.4,0.6,0.8,1.00))+ 
  theme(legend.position="bottom", legend.text=element_text(size=12))



## ggplot for b=1/3
ggplot(df.long.both13, aes(theta, value))+geom_line(aes(color=Type,linetype=Test),linewidth=1)+
  theme_bw()+theme(plot.title=element_text(hjust=0.5,size=12))+
  ggtitle("(b) b = 1/3")+
  xlab(expression(theta))+ ylab("power")+
  scale_linetype_manual(values=c("dotted","dashed","solid","dotdash"),
                        labels=c(expression(T[pr]),expression(T["pr,adj"]),
                                 expression(T[cr]),expression(T["cr,adj"])))+
  scale_color_manual(values=c("black","plum3"))+
  facet_grid(.~lambda,labeller = label_parsed)+
  scale_y_continuous(breaks=c(0.00,0.025,0.2,0.4,0.6,0.8,1.00))+ 
  theme(legend.position="bottom", legend.text=element_text(size=12))
