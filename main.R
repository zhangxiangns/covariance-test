setwd('/Users/zhangxiang/GitHub/协方差检验')
source('协方差检验.R')
data = read.csv('example3_3_2.csv')
cov.test.multi(data)


data = read.csv('Japan.csv')
data1 = subset(data,group==1,-group)
data1 = as.matrix(data1)
data2 = subset(data,group==2,-group)
data2 = as.matrix(data2)
cov.test.bi(data1,data2)


covs = var(data1)
cov.test.single(data1,covs)
