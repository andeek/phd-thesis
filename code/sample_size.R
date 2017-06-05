

get.quantile<-function(data,p){
## p is the quantile of interest
## just use the default quantile function in R

a<-as.vector(quantile(data,p))
return(a)
}



### estimates the variance/spectral density at 0 needed
### in large-sample variance formula
### of the p sample quantile


spectral.est<-function(data,p,M){
## p is the quantile of interest
## M is the bandwidth in flat-top kernel

ep<-get.quantile(data,p)
y<-data*0
y[data<ep]<-1

res<-acf(y, lag.max = M, type = "covariance", plot = FALSE)[[1]] 
res<-as.vector(res)

t<-seq(0,M,1)/M
weight<-2*(1-t)
weight[t<.5]<-1

a<-2*res*weight
a[1]<-a[1]/2
return(sum(a))
}


## estimate the bandwidth M
## for spectral density estimation

spectral.bandwidth<-function(data,p){
## p is the quantile of interest
## M is the bandwidth in flat-top kernel
## n>5 needed


n<-length(data)
K<-5
bound<-2*sqrt(log(n)/n)


ep<-get.quantile(data,p)
y<-data*0
y[data<ep]<-1

### compute acf & remove sample acf at 0
res<-acf(y, lag.max=n-1,type = "correlation", plot = FALSE)[[1]] 
res<-abs(as.vector(res))
res<-res[-1]
 
m<-0
a<-2*bound



while(a>bound){
m<-m+1
a<-max(res[(m+1):(m+K)]) 
if(m>(n-K-1)){a<-2*bound}
}

 
M<-2*m
return(M)
}




## estimate the bandwidth H
## for pdf estimation

pdf.bandwidth<-function(data){
## p is the quantile of interest
## M is the bandwidth in flat-top kernel
## n>5 needed


n<-length(data)
K<-5
bound<-2*sqrt(log(n)/n)


f<-function(t){
a1<-data*cos(seq(1,n,1)*(t+m))
a2<-data*sin(seq(1,n,1)*(t+m))
a1<-sum(a1)
a2<-sum(a2)
a<-sqrt(a1^2+a2^2)/n
return(a)
}
 
m<-0
b<-2*bound



while(b>bound){
m<-m+1

f<-function(t){
a1<-data*cos(seq(1,n,1)*(t+m))
a2<-data*sin(seq(1,n,1)*(t+m))
a1<-sum(a1)
a2<-sum(a2)
a<-sqrt(a1^2+a2^2)/n
return(a)
}

b<-optimize(f, interval = c(0,K), maximum = TRUE)$objective
 if(m>(n-K-1)){b<-2*bound}
}


H<-2*m
return(H)
}



pdf.hat<-function(data,x,H){
### kernel estimate of pdf at x
### using bandwidth H

t<-x-data
n<-length(data)
t<-t[t!=0]
a4<-(n-length(t))*(3*H/2-3*H^2/4)/n



a1<- sin(H*t)/t
a1<-sum(a1)*(2-2*H)/n


a2<- sin(H*t/2)/t
a2<-sum(a2)*(H-1)/n

a3<-cos(H*t)-cos(H*t/2)
a3<-a3/(t^2)
a3<-sum(a3)*(-2)/n

return((a1+a2+a3+a4)/(pi))
}


#### estimate the sample size
#### for precision d and confidence 
#### level C (as %) with quantile p

S.hat<-function(data,p,d,C){


ep<-get.quantile(data,p)

H<-pdf.bandwidth(data)

f.hat<-pdf.hat(data,ep,H)
  
M<-spectral.bandwidth(data,p)
V<-spectral.est(data,p,M)

q<-qnorm((1-C)/2)

s<- ceiling(q^2*V/(ep*d*f.hat)^2)
return(s)
}




