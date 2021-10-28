########################################################

#Estimação do erro tipo I do teste de Hotelling.

########################################################
rm(list=ls())## Limpar a Memória do R

########## simulação via Monte Carlo ##############

#### Pacotes Utilizados

library(mvtnorm)
library(ggplot2)
library(ICSNP)
library(xtable)

#==============================================================

#######    Estimação para alpha = 0.01 ##############
#n= 10
#m= 1000
#==============================================================

set.seed(1987)
n<-10
alpha<- .01 # fazendo alpha a 1%
mu<- rep(0,3) #vetor de média
variancia <- length(mu)
sigma <-diag(variancia)
p<- length(mu)
m<- 1000 #Número de simulação Monte Carlo
pvalue<- numeric(m)#armazenamento de p-valores
tempo.inicio<- Sys.time()

require(mvtnorm)


for(j in 1:m){
  x<-rmvnorm(n, mu, sigma) # gerando valores aleatórios seguindo dist normal bivariada
  testet2 <- n%*%t((apply(x,2,mean)-mu))%*%solve(sigma)%*%(apply(x,2,mean)-mu) # teste t^2 de hotelling
  pvalue[j]<- testet2
}

p.10.1000.1hat<- mean(pvalue<(((n-1)*p)/(n-p))*qf(alpha,p,n-p))
se.10.1000.1hat<-sqrt(p.10.1000.1hat*(1-p.10.1000.1hat)/m)

tempo.fim<- Sys.time()
tempo.exec<-tempo.fim - tempo.inicio;tempo.exec
#==============================================================

#para n=50 e m= 1000

#==============================================================

n<-50
alpha<- .01 # fazendo alpha a 1%
mu<- rep(0,3) #vetor de médias
variancia <- length(mu)
sigma <-diag(variancia)
p<- length(mu)
m<- 1000 #Número de simulação Monte Carlo
pvalue<- numeric(m)#armazenamento de p-valores
tempo.inicio<- Sys.time()

require(mvtnorm)

for(j in 1:m){
  x<-rmvnorm(n, mu, sigma) # gerando valores aleatórios seguindo dist normal bivariada
  testet2 <- n%*%t((apply(x,2,mean)-mu))%*%solve(sigma)%*%(apply(x,2,mean)-mu) # teste t^2 de hotelling
  pvalue[j]<- testet2
}

p.50.1000.1hat<- mean(pvalue<(((n-1)*p)/(n-p))*qf(alpha,p,n-p))
se.50.1000.1hat<-sqrt(p.50.1000.1hat*(1-p.50.1000.1hat)/m)

tempo.fim<- Sys.time()
tempo.exec<-tempo.fim - tempo.inicio;tempo.exec

#==============================================================