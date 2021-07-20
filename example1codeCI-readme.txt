example1codeCI  has four functions: lm.model，SDB.lm，BLB.lm  and lm.boot. 


lm.model

Description
The function implemrnts  DAC-EL for linear model.

Usage
lm.moel(x, y, beta0=NULL,  K=NULL, type="consecutive", Intercept=FALSE, alph=0.05)

Arguments
x: the design matrix.
y: the response vector.
beta0:    null hypothesis of  the parameter. The default is 0.
K:  the number of blocks.
type:  a character string specifying the type of blocks to be generated.  "random"  means is randomly  splited observation into K blocks.
      "consecutive" is  splited by natural number. 
Intercept: whether to calculate the intercept for this model. If set to False, no intercept will be used in calculations.
alph: significant level.  The default is 0.05.


SDB.lm

Description
The function implemrnts  SDB for linear model.

Usage
SDB.lm (x, y, S=NULL , b=NULL, beta0=NULL, alpha=0.05)
Arguments
x: the design matrix.
y: the response vector.
S: the number of subsets.
b: the  number of sampled subset.
beta0:    null hypothesis of  the parameter. The default is 0.
alpha: significant level.  The default is 0.05.


BLB.lm

Description
The function implemrnts  SDB for linear model.

Usage
BLB.lm (x, y, S=NULL, r=100, b=NULL, beta0=NULL, alpha=0.05)

Arguments
x: the design matrix.
y: the response vector.
S: the number of subsets.
r: the number of sampled subset.
b: the  number of sampled subset.
beta0:    null hypothesis of  the parameter. The default is 0.
alpha: significant level.  The default is 0.05.


lm.boot

Description
The function implemrnts  TB for linear model.

Usage
lm.boot(x, y, B=100, beta0=NULL, alpha=0.05)

Arguments
x: the design matrix.
y: the response vector.
B: the  number of replicatment.
beta0:    null hypothesis of  the parameter. The default is 0.
alpha: significant level.  The default is 0.05.



