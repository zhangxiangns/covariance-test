### Population N(mu,sigma)
### The function is to test whether sigma is equal to sigma0. ###
### input:
### mydata--a matrix of n*p
### sigma0--a matrix of p*p
### note:even though that the methods in books test respectively the sigma0 where the
###      sigma0 is equal to Ip or some others, the function integrates both of the methods
###      together.
cov.test.constants=function(mydata,sigma0){
  n = nrow(mydata)
  p = ncol(mydata)
  A = var(mydata)*(n-1)
  lamda = exp(sum(diag(-1/2*A%*%solve(sigma0))))*det(A%*%solve(sigma0))^(n/2)*(exp(1)/n)^(1/2*n*p)
  statistic = -2*log(lamda)
  p_value = 1-pchisq(statistic,p*(p+1)/2)
  output = list()
  output[['p_value']] = p_value
  output[['statistics']] = statistic
  return(output) 
}

### Population N(mu1,sigma1),N(mu2,sigma2),...,N(muk,sigmak)
### H0:sigma1=sigma2=...=sigmak H1:sigma1,sigma2,...,sigmak are not the same.
### input:mydata
### mydata is a data frame such as:
####### X1   X2   X3   X4   group
####### x11  x12  x13  x14  1
####### x21  x22  x23  x24  1
####### ...  ...  ...  ...  ...
####### xn1  xn2  xn3  xn4  k
### Attention: the classification column needs to be 'group'!
### note:the function is applied in function analysis_variance.
cov.test.multi = function(mydata){
  groups = unique(mydata$group)
  k = length(groups)
  n = rep(0,k)
  N = nrow(mydata)
  p = ncol(mydata)-1
  AA = list()
  A = 0
  for (i in 1:k){
    label = groups[i]
    subdata = subset(mydata,group==label,-group)
    subdata = as.matrix(subdata)
    n[i] = nrow(subdata)
    AA[[i]] = var(subdata)*(n[i]-1)
    A = A + AA[[i]]
  }
  M = (N-k)*log(det(A/(N-k)))
  for (i in 1:k){
    M = M - (n[i]-1)*log(det(AA[[i]]/(n[i]-1)))
  }
  if (sum(n==n[1])==k){
    d = (2*p^2+3*p-1)*(k+1)/(6*(p+1)*(N-k))
  }
  if (sum(n==n[1])!=k){
    d = (2*p^2+3*p-1)/(6*(p+1)*(k-1))*(sum(1/(n-1))-1/(N-k))
  }
  statistic = (1-d)*M
  p_value = 1-pchisq(statistic,1/2*p*(p+1)*(k-1))
  output = list()
  output[['p_value']] = p_value
  output[['statistics']] = statistic
  return(output) 
}

### Population N(mu1,sigma1),N(mu2,sigma2)
### The function is to test whether sigma1 is equal to sigma2.
### input:
### mydata1--a matrix of n1*p.
### mydata2--a matrix of n2*p.
### note:this function is applied the function bi.mean.test.

cov.test.bi = function(mydata1,mydata2){
  n1 = nrow(mydata1)
  n2 = nrow(mydata2)
  A1 = var(mydata1)*(n1-1)
  A2 = var(mydata2)*(n2-1)
  p = ncol(mydata1)
  lamda = (n1+n2)^(p*(n1+n2-2)/2)*(det(A1))^((n1-1)/2)*(det(A2))^((n2-1)/2)
  lamda = lamda/((n1-1)^(p*(n1-1)/2)*(n2-1)^(p*(n2-1)/2)*(det(A1+A2))^((n1+n2-2)/2))
  d = (2*p^2+3*p-1)/(6*(p+1))*(1/(n1-1)+1/(n2-1)-1/(n1+n2-2))
  statistic = -(1-d)*log(lamda)
  p_value = 1 - pchisq(statistic,p*(p+1)/2)
  output = list()
  output[['p_value']] = p_value
  output[['statistics']] = statistic
  return(output) 
}