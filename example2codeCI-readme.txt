example1codeCI  has five functions: lm.model，SDB.lm，BLB.lm,  lm.boot and twostep.

glm.model

Description
The function implemrnts  DAC-EL for logistic regression

Usage
glm.model(x, y, beta0=NULL,  K=NULL, type="consecutive", family=binomial(link = "logit"), Intercept=FALSE, alph=0.05)

Arguments
x: the design matrix.
y: the response vector.
beta0:    null hypothesis of  the parameter. The default is 0.
K:  the number of blocks.
type:  a character string specifying the type of blocks to be generated.  "random"  means is randomly  splited observation into K blocks.
      "consecutive" is  splited by natural number. 
Intercept: whether to calculate the intercept for this model. If set to False, no intercept will be used in calculations.
family： a description of the error distribution and link function to be used in the model. 
alph: Significant level.  The default is 0.05.


SDB.glm

Description
The function implemrnts  SDB for logistic regression.

Usage
SDB.glm (x, y, S=NULL , b=NULL, beta0=NULL, alpha=0.05)
Arguments
x: the design matrix.
y: the response vector.
S: the number of subsets.
b: the  number of sampled subset.
beta0:    null hypothesis of  the parameter. The default is 0.
alpha: significant level.  The default is 0.05.


BLB.glm

Description
The function implemrnts  SDB for logistic regression

Usage
BLB.glm (x, y, S=NULL, r=100, b=NULL, beta0=NULL, alpha=0.05)

Arguments
x: the design matrix.
y: the response vector.
S: the number of subsets.
r: the number of sampled subset.
b: the  number of sampled subset.
beta0:    null hypothesis of  the parameter. The default is 0.
alpha: sgnificant level.  The default is 0.05.


glm.boot

Description
The function implemrnts  TB for logistic regression.

Usage
lm.boot(x, y, B=100, beta0=NULL, alpha=0.05)

Arguments
x: the design matrix.
y: the response vector.
B: the  number of replicatment.
beta0:    null hypothesis of  the parameter. The default is 0.
alpha: significant level.  The default is 0.05.



twostep

Description
The function implemrnts  the method of Wang, Zhu, and Ma (2018).  for logistic regression

Usage
twostep(X, Y, beta0=NULL, r1, r2, method=c("mvc", "mmse", "uni"),alph=0.05)


Arguments
x: the design matrix.
y: the response vector.
B: the  number of replicatment.
beta0:    null hypothesis of  the parameter. The default is 0.
r1: the subsample size of the first step.
r2: the subsample size of the second step.
methd:  the kind of the method.
alph: significant level.  The default is 0.05.

