########------practical 1: least square estimation and testing-------############
#__Q1:-
y=c(56.7,65.9,75.7,64.2,129.2,90.7)
y
x=matrix(c(1,1,1,1,1,1,1,0,0,1,1,1,1,1,0,0,1,1,0,1,1,0,1,1,0,0,1,1,1,0),byrow=FALSE,nrow=6,ncol=5)
x
#estimate weight of object using least sq
bhat=solve(t(x)%*%x)%*%t(x)%*%y
bhat
#for obtain the estimate of disp matrix
v=var(y)
disp=v*solve(t(x)%*%x)
disp
#for checking estamablity i.e lamda_prime*h =labmda_prime   then estimable
H=solve(t(x)%*%x)%*%(t(x)%*%x)
H
h=round(H,2)
h
#(i):
lamda=c(0,1,1,1,0)
ans=t(lamda)%*%h
ans
#hence given function is estimable
#(ii):-
lamda1=c(0,1,1,-2,0)
ans1=t(lamda1)%*%h
ans1
#hence given function is estimable
#(iii):-
lamda2=c(0,1,1,0,-1)
ans2=t(lamda2)%*%h
ans2
#hence given function is estimable

##--Q2-----


























####--------------------practical 2 Two way annova----------------------######
##Q1:-------
y=c(580,568,570,1080,1087,1392,1380,1386,550,1070,1035,1000,1328,1312,575,599,1045,1066,867,904,889)
y
n=length(y);n
Y=matrix(c(1718,2167,4158,550,3105,2640,1174,2111,2660),nrow=3,ncol=3,byrow=T)   ##sum of each block column wise
Y                                                                      #i.e 580+568+570=1718

N=matrix(c(3,2,3,1,3,2,2,2,3),nrow=3,ncol=3,byrow=T)
N
A=rowSums(Y)
A
B=colSums(Y)
B
r=rowSums(N)
r
k=colSums(N)
k
da=diag(r);da
db=diag(k);db
tss=sum(y^2)-(sum(Y)^2/n)   #Y...=sum(Y)    #tss=sum(yijk^2)-((Y..)^2/n)  
tss
CT=sum(Y)^2/n    #for more convinience
c=da-(N%*%solve(db)%*%t(N))
c
Q=A-(N%*%solve(db)%*%B)
Q
library(MASS)
alpha=ginv(c)%*%Q
alpha
ssa_unadj=sum(A^2/r)-CT
ssa_unadj
ssb_unadj=sum(B^2/k)-CT
ssb_unadj
ssa_adj=t(Q)%*%alpha
ssa_adj
sse=tss-ssa_adj-ssb_unadj
sse
ssb_adj=tss-sse-ssa_unadj
ssb_adj

#genrally blocks are column 
p=3;q=3
ftab=qf(0.05,(p-1),(n-p-q+1))
ftab
F1=(ssa_adj/(p-1))/(sse/(n-p-q+1))
F1   
#  fcal > ftab hence Reject h0 at 5%
F2=(ssb_adj/(q-1))/(sse/(n-p-q+1))
F2          
#fcal>ftab hence Definitely reject h0 at 5%

##-----> Q.2:-

y=c(580,568,570,1080,1087,1085,1392,1380,1386,550,530,579,1070,1035,1000,1328,1312,1299,546,575,579,1045,1053,1066,867,904,889)
n=length(y)
n
Y=matrix(c(1718,3252,4158,1659,3105,3939,1700,3164,2660),nrow=3,ncol=3,byrow = T);Y
N=matrix(c(3,3,3,3,3,3,3,3,3),nrow=3,ncol=3)
N
A=rowSums(Y);A
B=colSums(Y);B
r=rowSums(N);r
k=colSums(N);k
Da=diag(r);Da
Db=diag(k);Db

C=Da-(N%*%solve(Db)%*%t(N));C
Q=A-(N%*%solve(Db)%*%B);Q

CT=sum(Y)^2/n;CT
TSS=sum(y^2)-CT ; TSS

SSA_unadj=sum(A^2/r)-CT 
SSA_unadj
SSB_unadj=sum(B^2/k)-CT
SSB_unadj

library(MASS)
alpha =ginv(C)%*%Q
alpha

SSA_adj=t(Q)%*%alpha
SSA_adj
SSE=TSS-SSA_adj-SSB_unadj;SSE
SSB_adj=TSS-SSE-SSA_unadj
SSB_adj

p=q=3
F1=(SSA_adj/(p-1))/(SSE/(n-p-q+1))
F1
F2=(SSB_adj/(q-1))/(SSE/(n-p-q+1))
F2
ftab=qf(0.05,(p-1),(n-p-q+1))
ftab
###here both are facl>ftab hence reject H0 '[][=]========at 5%

####--------Practical 3: stochastic process 1--------####
#Q1_----

tpm = matrix(c(1/3,1/3,1/3,0,1/3,2/3,3/4,1/4,0),ncol=3,nrow=3,byrow=TRUE);tpm
states = c(1,2,3)
intdist = c(1/3,1/3,1/3)
p1 = c(tpm[1,]);p1
p2 = c(tpm[2,]);p2
p3 = c(tpm[3,]);p3

X = c() 
X[1] = sample(states,1)
for(i in 1:500){
  if(X[i]==1){
    X[i+1] = sample(states,1,prob=p1) 
  }
  else if(X[i]==2){
    X[i+1] = sample(states,1,prob=p2)  
  }
  else if(X[i]==3){
    X[i+1] = sample(states,1,prob=p3)   
  }
}
x
#est_tpm=markovchainFit(data=X)


n = c(500,1000,5000,10000)

for(j in 1:length(n)){
  X = c() 
  X[1] = sample(states,1)
  for(i in 1:n[j]){
    if(X[i]==1){
      X[i+1] = sample(states,1,prob=p1) 
    }
    else if(X[i]==2){
      X[i+1] = sample(states,1,prob=p2)  
    }
    else if(X[i]==3){
      X[i+1] = sample(states,1,prob=p3)   
    }
  }
  
  c1=c2=c3=0
  sum11=sum12=sum13=sum21=sum22=sum23=sum31=sum32=sum33=0
  
  for(i in 1:n[j]){
    if(X[i]==1){
      c1=c1+1
      if(X[i+1]==1){
        sum11=sum11+1
      }
      if(X[i+1]==2){
        sum12=sum12+1
      }
      if(X[i+1]==3){
        sum13=sum13+1
      }
    }
    if(X[i]==2){
      c2=c2+1
      if(X[i+1]==1){
        sum21=sum21+1
      }
      if(X[i+1]==2){
        sum22=sum22+1
      }
      if(X[i+1]==3){
        sum23=sum23+1
      }
    }
    if(X[i]==3){
      c3=c3+1
      if(X[i+1]==1){
        sum31=sum31+1
      }
      if(X[i+1]==2){
        sum32=sum32+1
      }
      if(X[i+1]==3){
        sum33=sum33+1
      }
    }
  }
  est_tpm_count = matrix(c(sum11,sum12,sum13,sum21,sum22,sum23,sum31,sum32,sum33),nrow=3,ncol=3,byrow=TRUE)
  count=c(c1,c2,c3)
  est_tpm=matrix(nrow=3,ncol=3)
  for(i in 1:3){
    est_tpm[i,] = est_tpm_count[i,]/count[i]
  }
  print(round(est_tpm,2))
}
round(tpm,2)
est_tpm
#############################
#que 2:

tpm2 = matrix(c(1/2,1/2,0,0,3/4,1/4,0,0,0,0,1/4,3/4,0,0,1/2,1/2),ncol=4,nrow=4,byrow=TRUE)
states = c(1,2,3,4)
intdist = c(1/4,1/4,1/4,1/4)
p1 = c(tpm2[1,])
p2 = c(tpm2[2,])
p3 = c(tpm2[3,])
p4 = c(tpm2[4,])

n = c(500,1000,5000,10000)

for(j in 1:length(n)){
  X = c() 
  X[1] = sample(states,1)
  for(i in 1:n[j]){
    if(X[i]==1){
      X[i+1] = sample(states,1,prob=p1) 
    }
    else if(X[i]==2){
      X[i+1] = sample(states,1,prob=p2)  
    }
    else if(X[i]==3){
      X[i+1] = sample(states,1,prob=p3)   
    }
    else if(X[i]==4){
      X[i+1] = sample(states,1,prob=p4)   
    }
  }
  
  c1=c2=c3=c4=0
  sum11=sum12=sum13=sum14=sum21=sum22=sum23=sum24=sum31=sum32=sum33=sum34=sum41=sum42=sum43=sum44=0
  
  for(i in 1:n[j]){
    if(X[i]==1){
      c1=c1+1
      if(X[i+1]==1){
        sum11=sum11+1
      }
      if(X[i+1]==2){
        sum12=sum12+1
      }
      if(X[i+1]==3){
        sum13=sum13+1
      }
      if(X[i+1]==4){
        sum14=sum14+1
      }
    }
    if(X[i]==2){
      c2=c2+1
      if(X[i+1]==1){
        sum21=sum21+1
      }
      if(X[i+1]==2){
        sum22=sum22+1
      }
      if(X[i+1]==3){
        sum23=sum23+1
      }
      if(X[i+1]==4){
        sum24=sum24+1
      }
    }
    if(X[i]==3){
      c3=c3+1
      if(X[i+1]==1){
        sum31=sum31+1
      }
      if(X[i+1]==2){
        sum32=sum32+1
      }
      if(X[i+1]==3){
        sum33=sum33+1
      }
      if(X[i+1]==4){
        sum34=sum34+1
      }
    }
    if(X[i]==4){
      c4=c4+1
      if(X[i+1]==1){
        sum41=sum41+1
      }
      if(X[i+1]==2){
        sum42=sum42+1
      }
      if(X[i+1]==3){
        sum43=sum43+1
      }
      if(X[i+1]==4){
        sum44=sum44+1
      }
    }
  }
  est_tpm_count = matrix(c(sum11,sum12,sum13,sum14,sum21,sum22,sum23,sum24,sum31,sum32,sum33,sum34,sum41,sum42,sum43,sum44),nrow=4,ncol=4,byrow=TRUE)
  count=c(c1,c2,c3,c4)
  est_tpm=matrix(nrow=4,ncol=4)
  for(i in 1:4){
    est_tpm[i,] = est_tpm_count[i,]/count[i]
  }
  print(round(est_tpm,2))
}


####-----prac 4:-   branching process --------######
#Q1:---------

ex=sum(0.2+2*0.15+3*0.1+4*0.05)
ex2=sum(0.2+4*0.15+9*0.1+16*0.05)
var=ex2-(ex^2)
var;ex

z = c()
z[1] = 1
p = c(0.5,0.2,0.15,0.1,0.05);p
s = c(0:4);s     #0:4
z[2]=sample(s,z[1],prob=p)

for(i in 2:5){
  print(sample(x=s,size=z[i],prob=p,replace=TRUE))
  z[i+1] = sum(sample(x=s,size=z[i],prob=p,replace=TRUE)) 
}
z
sum(p*s)
(sum(p*s))^(length(z))

f= function(s){
  return(0.5+0.2*s+0.15*s^2+0.1*s^3+0.05*s^4)
}
s = seq(from=0,to=1,by=0.001)
fs = f(s)
plot(s,fs,type='l')
abline(0,1)


#######-----practiacl 4;; BIBD  ------######
###Q1:---------
p=5    #no of tratment
q=5    #no of block
r=4    #no of element each row or no of tratment in row or replication
k=4    #no tratment in each block
lamda=(k-1)*r/(p-1)
lamda
n=r*p
n
y= matrix(c(17,0,12,13,14,14,14,10,13,0,0,12,9,12,13,11,13,0,12,11,12,11,8,0,10),nrow=5,ncol=5,byrow=T);y
nij=ifelse(y!= 0,1,0);nij
yidot=rowSums(y)
ydotj=colSums(y)
ct=sum(y)^2/n
tss=sum(y^2)-ct
tss
ssb=(1/k)*sum(ydotj^2)-ct
ssb
Q=c()
for(i in 1:p){
  Q[i]=yidot[i] - (1/k)* nij[i,]%*%ydotj
}
Q
sst_adj = (k/(lambda*p)) * sum (Q^2) ; sst_adj
sse= tss - sst_adj - ssb ;ssb
mst = SST_adjust /(r-1) ;mst
mse= sse/ (n-p-q+1) ; mse
fcal = mst / mse ;fcal
ftab=qf(0.05,(r-1),(n-p-q+1))
ftab
