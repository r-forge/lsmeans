# Package lsmeans
# Version 1.0
#
# Calculates the lsmeans and their variances for lm, gls or lme models.
# It also contains a function for mean comparisons and for contrasting linear
# combination of parameters in the model
# The present version only works with default parametrization.
#
# Writen by
#
# Julio A. Di Rienzo dirienzo@agro.unc.edu.ar
# http://sites.google.com/site/juliodirienzo/
#
# The package contains the following main functions.
# Other functions are for internal use.
#
#     CalculateMmatrix
#     CalculateMatrixofMeansAndStandardErrors
#     MeanAndVariance
#     MeansComparison
#     contrastLC
#     RSquare
#
# Examplification of lsmeans package.
# Load the package
# Let suppose that myData is a dataframe containing the columns:
# "y", "Treatment" and "Sex". "y" represent the response, "Treatment" a factor
# with levels T1,T2,T3 and "Sex" a factor wiht levels F and M.
# Let an object named "myModel" contains the fitted lm, gls or mle model
# let suppose that the formula for this model is y~Treatment+Sex+Treatment:Sex
# let an object named "myData" contains the data included in the "data"
# argument of the lm, gls or mle functions.
#
# Mmatrix<-CalculateMmatrix(myModel,myData)
#(CalculateMmatrix calculates a matrix needed for other routines
# this function is computationally remanding. There is no need to call it
# more than once for a given model)
#
# The means and standard errors for each level of "Treatment" are requested by
#
# CalculateMatrixofMeansAndStandardErrors(myModel,myData,"Treatment",Mmatrix)
#
# The means and estandard errors for each level of the interaction Treatment by
# Sex are requested by
#
# CalculateMatrixofMeansAndStandardErrors(myModel,myData,"Treatment:sex",Mmatrix)
#
# To obtain de coefficient of linear combination need to get the matrix of means
# and estandard error use:
#
# LinearCombinationsCoefficientsToGetMean(myModel,myData,myModelTerm,Mmatrix)
#
# It works in the same maner as it does  CalculateMatrixofMeansAndStandardErrors
# but insted of returning a matrix of means and standard errors it returns
# the coefficientes of the linear combination of parametres to get those
# estimates.
#
# To obtain a single mean, its variance and the linear combination used to
# generate these estimates call MeanAndVariance(). For example if
# the mean of say level T2 of "Treatment" is required try the following
# call
#
# MeanAndVariance("TreatmentT2", Mmatrix,myModel)
#
# The function returns a vector containing the mean, the variance and a
# vector of coefficients used to calculate the mean and variance.
# For example let x be the vector of coefficients returned by MeanAndVariance,
# b the vector of estimated parameters of the model and S the
# variance-covariance matrix of b. Then mean=x'b and variance=x'Sx
# Other example:
# To obtain the mean of level T2 of factor Treatment and the level "F"
# factor sex make the following call
#
# MeanAndVariance("TreatmentT2:SexF", Mmatrix,myModel)
#
# To compare the means of levels T1 and T2 of "Treatment"  try
#
# MeansComparison("TreatmentT1,"TreatmentT2",Mmatrix,myModel,myDF)
# This function returns the p-value for the comparison. the argument myDF is
# optional. If specified, it represents the degree of freedom o the appropriate
# term for the comparison and the function performs a test like simmilar to
# a T-test. Otherwise it uses the anova function for a linear combination of
# the parameter. Usually you will not specify this argument or set it to zero
#
# To test a single or a set of linear combinations use
#
# contrastLC(myModel,myCombination)
#
# The function returns $p=p-value, $Fobs= the calculated F-value, and $numDF=
# the numerator degree of freedom of the F-test.
# The argument myCombination can be a single vector or a matrix of combinations,
# each combination corresponding to a row.
#
# To obtain the R-square of the fitted model use RSquare(myModel). If the model
# is a mixed effects model this function will retur an R-square for the fixed
# and the nested random effects (the R-square will increase with the
# consecutive random effects)

#---------------------------------------------------------------------------

SearchForASubstringInAStringWithSeparators<-function(s1,s2,separator)
 {
 Result<-FALSE
 MyListstrings<-strsplit(s2,separator)[[1]]
 myindices<-grep(s1,MyListstrings,fixed =TRUE)
 i<-0
 if (length(myindices)>0)
{
 while ((Result==FALSE)&(i<length(myindices)))
     {
      i<-i+1
      Result<-(Result|(s1==MyListstrings[myindices[i]]))
     }
 }
 Result
 }
#---------------------------------------------------------------------------

 SearchForASubstringInAListOfStringsWithSeparators<-function(s1,MyList,separator)
 {
 Result<-grep(s1,c())
 myindices<-grep(s1,MyList,fixed =TRUE)
 if (length(myindices)>0) for (i in (1:length(myindices))) if (SearchForASubstringInAStringWithSeparators(s1,MyList[myindices[i]],separator)) Result<-c(Result,myindices[i])
 Result
 }

#---------------------------------------------------------------------------

CalculateMmatrix<-function(myModel,myData)
{
 #---------------------------------------------------------------------------

 for (i in (1:ncol(myData))) if (is.factor(myData[,i]))  myData[,i]<-as.factor(as.character(myData[,i]))

 if (class(myModel)=="lme") f<-colnames(attributes(terms(myModel))$factors)
 if (class(myModel)=="gls") f<-names(myModel$parAssign)
 if (class(myModel)=="lm") f<-labels(myModel$terms)

 f<-strsplit(f, ":")
 myModelFactors<- c()
 for (i in (1:length(f)))  myModelFactors<- c(myModelFactors,f[[i]])
 myModelFactors<- unique(myModelFactors)
 myDataColnames<-colnames(myData)
 factores<-c()
 for (i in (1:length(myModelFactors)))
 {
  x<-SearchForASubstringInAListOfStringsWithSeparators(myModelFactors[i],myDataColnames,":")
 if (length(x)>0) {if (is.factor(myData[,x])==TRUE) factores<-c(factores,myDataColnames[x])}
 }

 MyList=paste(rep(factores[1],nrow(myData)),myData[,factores[1]],sep='')
 if (length(factores)>1) for (i in (2:length(factores)))
    {
    MyList=paste(MyList,paste(rep(factores[i],nrow(myData)),myData[,factores[i]],sep=''),sep=':')
    }
 MyList=unique(MyList)

 #---------------------------------------------------------------------------
 # MyList contains all the combination of levels of fixed effects in myModel
 #---------------------------------------------------------------------------
 if (class(myModel)=="lme") NameOfCoefficients<-names(myModel$coefficients$fixed)
 if (class(myModel)=="gls") NameOfCoefficients<-names(myModel$coefficients)
 if (class(myModel)=="lm")  NameOfCoefficients<-names(myModel$coefficients[complete.cases(myModel$coefficients)])

 MyCoefficientes<-matrix(rep(0,length(MyList)*length(NameOfCoefficients)),nrow=length(MyList),ncol=length(NameOfCoefficients))
 constant<-as.integer(length(grep("(Intercept)",NameOfCoefficients,fixed =TRUE))>0)

 for (j in (1:length(NameOfCoefficients)))
     {
      ss<-strsplit(NameOfCoefficients[j], ":")[[1]]
      myindices<-SearchForASubstringInAListOfStringsWithSeparators(ss,MyList,':')
      if (length(myindices)>0)
      {
      for (i in (1:length(myindices)))
      {
      x<-SearchForASubstringInAListOfStringsWithSeparators(ss[1],MyList[myindices[i]],":")
      if (length(ss)>1) for (ii in (2:length(ss))) x<-intersect(x,SearchForASubstringInAListOfStringsWithSeparators(ss[ii],MyList[myindices[i]],":"))
      if (length(x)>0) MyCoefficientes[myindices[i],j]<-1
      }
      }
     }
 if (constant>0) MyCoefficientes[,1]<-constant

 for (i in (1:length(colnames(myData))))
 if (is.numeric(myData[,i]))
 {
 g<-SearchForASubstringInAListOfStringsWithSeparators(colnames(myData)[i],NameOfCoefficients,":")
 if (length(g)>0) for (j in (1:length(g)))
        {
            sc=colnames(myData)[i]
            ss=strsplit(NameOfCoefficients[g[j]], ":")[[1]]
            ss=setdiff(ss,sc)
            if (length(ss)>1) ss=paste(ss,collapse=":")
            myindices<-c()
            if (length(ss)>0) myindices<-grep(ss,MyList,fixed =TRUE)
            mx=mean(myData[,i],na.rm = TRUE)
            if   (length(myindices)==0) MyCoefficientes[,g[j]]<-mx
            if   (length(myindices)>0)
                 for (k in (1:length(myindices))) MyCoefficientes[myindices[k],g[j]]=mx
        }
 }

 rownames(MyCoefficientes)<-MyList
 colnames(MyCoefficientes)<-NameOfCoefficients
 MyCoefficientes
 }
 #---------------------------------------------------------------------------

 linear.combination.to.estimate.a.mean<-function(tratamiento,Mmatrix,MyCoefficientes)
 {
  ss<-strsplit(tratamiento, ":")[[1]]
  x<-SearchForASubstringInAListOfStringsWithSeparators(ss[1],rownames(Mmatrix),":")
  if (length(ss)>1) for (i in (2:length(ss))) x<-intersect(x,SearchForASubstringInAListOfStringsWithSeparators(ss[i],rownames(Mmatrix),":"))

  ifelse (length(x)>0, SM<-Mmatrix[x,],SM<-c(1,rep(0,(length(MyCoefficientes)-1))))
  x<-as.matrix(rep((1/nrow(SM)),nrow(SM)))
  ifelse ((nrow(x)==0),x<-as.matrix(SM), x<-as.matrix(t(x)%*%SM))
  if (nrow(x)>1) x<-t(x)
  colnames(x)<-names(MyCoefficientes)
  x
  }

 MeanAndVariance<-function(MeanLabel,Mmatrix,myModel)
  {

  if (class(myModel)=="lme") MyCoefficientes<-myModel$coefficients$fixed
  if (class(myModel)=="gls") MyCoefficientes<-myModel$coefficients
  if (class(myModel)=="lm")  MyCoefficientes<-myModel$coefficients[complete.cases(myModel$coefficients)]

  x<-linear.combination.to.estimate.a.mean(MeanLabel,Mmatrix,MyCoefficientes)

  if (class(myModel)=="lme") Result<-c(x%*%MyCoefficientes,x%*%myModel$varFix%*%t(x),x)
  if (class(myModel)=="gls") Result<-c(x%*%MyCoefficientes,x%*%myModel$varBeta%*%t(x),x)
  if (class(myModel)=="lm") Result<-c(x%*%MyCoefficientes,x%*%vcov(myModel)%*%t(x),x)

  Result
  }

#--------------------------------------------------------------------------
 RebuildTreatmentNames<-function(myData,myModelTerm)
{
 f<-strsplit(myModelTerm, ":")[[1]]
 myindices=complete.cases(myData[,f])
 MyList=paste(rep(f[1],nrow(myData)),myData[myindices,f[1]],sep='')
 if (length(f)>1) for (i in (2:length(f)))
    {
    MyList=paste(MyList,paste(rep(f[i],nrow(myData)),myData[myindices,f[i]],sep=''),sep=':')
    }
 MyList=unique(MyList)

 MyList
}

#--------------------------------------------------------------------------
CalculateMatrixofMeansAndStandardErrors<-function(myModel,myData,myModelTerm,Mmatrix)
{

 factores<-strsplit(myModelTerm, ":")
 f<-factores[[1]]

 MyList=paste(rep(f[1],nrow(myData)),myData[,f[1]],sep='')
 if (length(f)>1) for (i in (2:length(f)))
    {
    MyList=paste(MyList,paste(rep(f[i],nrow(myData)),myData[,f[i]],sep=''),sep=':')
    }
 MyList=unique(MyList)

 MyTable<-c()
 for (i in (1:length(MyList))) MyTable<-rbind(MyTable,c(MeanAndVariance(MyList[i],Mmatrix,myModel)[1],sqrt(MeanAndVariance(MyList[i],Mmatrix,myModel)[2])))
 rownames(MyTable)<-MyList
 colnames(MyTable)<-c("Mean","S.E.")
 MyTable
}

#--------------------------------------------------------------------------
LinearCombinationsCoefficientsToGetMean<-function(myModel,myData,myModelTerm,Mmatrix)
{
 factores<-strsplit(myModelTerm, ":")
 f<-factores[[1]]

 MyList=paste(rep(f[1],nrow(myData)),myData[,f[1]],sep='')
 if (length(f)>1) for (i in (2:length(f)))
    {
    MyList=paste(MyList,paste(rep(f[i],nrow(myData)),myData[,f[i]],sep=''),sep=':')
    }
 MyList=unique(MyList)

 MyTable<-c()
 for (i in (1:length(MyList))) MyTable<-rbind(MyTable,MeanAndVariance(MyList[i],Mmatrix,myModel)[c(-1,-2)])
 rownames(MyTable)<-MyList

  if (class(myModel)=="lme") MyCoefficientes<-myModel$coefficients$fixed
  if (class(myModel)=="gls") MyCoefficientes<-myModel$coefficients
  if (class(myModel)=="lm")  MyCoefficientes<-myModel$coefficients[complete.cases(myModel$coefficients)]

 colnames(MyTable)<-names(MyCoefficientes);
 MyTable
}
#---------------------------------------------------------------------------
MeansComparison<-function(MeanLabel1,MeanLabel2,Mmatrix,myModel,myDF=0)
{
  m1<-MeanAndVariance(MeanLabel1,Mmatrix,myModel)
  m2<-MeanAndVariance(MeanLabel2,Mmatrix,myModel)
  dif<-m1[3:length(m1)]-m2[3:length(m2)]

  if (class(myModel)=="lme") MyCoefficientes<-myModel$coefficients$fixed
  if (class(myModel)=="gls") MyCoefficientes<-myModel$coefficients
  if (class(myModel)=="lm" ) MyCoefficientes<-myModel$coefficients[complete.cases(myModel$coefficients)]

  if (class(myModel)=="lme") estadisticoT<-dif%*%MyCoefficientes/sqrt(dif%*%myModel$varFix%*%dif)
  if (class(myModel)=="gls") estadisticoT<-dif%*%MyCoefficientes/sqrt(dif%*%myModel$varBeta%*%dif)
  if (class(myModel)=="lm" ) estadisticoT<-dif%*%MyCoefficientes/sqrt(dif%*%vcov(myModel)%*%dif)

  if (  myDF==0)  ifelse((m1==m2),Result<-1,Result<-anova(myModel,L=dif)$p)
  if (!(myDF==0)) ifelse((m1==m2),Result<-1,Result<-2*(1-pt(abs(estadisticoT),df=myDF)))

  Result
  }
#---------------------------------------------------------------------------

contrastar<-function(myModel,myData,myModelTerm,Mmatrix,myContrast,myDF)
{
 factores<-strsplit(myModelTerm, ":")
 f<-factores[[1]]

 MyList=paste(rep(f[1],nrow(myData)),myData[,f[1]],sep='')
 if (length(f)>1) for (i in (2:length(f)))
    {
    MyList=paste(MyList,paste(rep(f[i],nrow(myData)),myData[,f[i]],sep=''),sep=':')
    }
 MyList=unique(MyList)

 MyTable<-c()
 for (i in (1:length(MyList)))
 {
 a<-MeanAndVariance(MyList[i],Mmatrix,myModel);
 MyTable<-rbind(MyTable,a[3:length(a)])
 }

  dif<-as.matrix(t(myContrast)%*%MyTable)
  if (class(myModel)=="lme") MyCoefficientes<-myModel$coefficients$fixed
  if (class(myModel)=="gls") MyCoefficientes<-myModel$coefficients
  if (class(myModel)=="lm")  MyCoefficientes<-myModel$coefficients[complete.cases(myModel$coefficients)]

  numDF<-nrow(dif)
  if (class(myModel)=="lme") ChiStatistic<-t(dif%*%MyCoefficientes)%*%solve(dif%*%myModel$varFix%*%t(dif))%*%(dif%*%MyCoefficientes)
  if (class(myModel)=="gls") ChiStatistic<-t(dif%*%MyCoefficientes)%*%solve(dif%*%myModel$varBeta%*%t(dif))%*%(dif%*%MyCoefficientes)
  if (class(myModel)=="lm")  ChiStatistic<-t(dif%*%MyCoefficientes)%*%solve(dif%*%vcov(myModel)%*%t(dif))%*%(dif%*%MyCoefficientes)

  Fstatistic<-ChiStatistic/numDF
  if (!(denDF==0)) ifelse((Fstatistic[1]==0),p<-1,p<-(1-pf(Fstatistic[1],numDF,myDF)))
  Result<-list("F"=Fstatistic[1],"numDF"=numDF,"denDF"=myDF,"p"=p)
  Result

}
#---------------------------------------------------------------------------

contrastLC<-function(myModel,myCombination)
{
 if (!(class(myModel)=="lm")) {
  a<-anova(myModel,L=t(as.matrix(myCombination)))
  Result<-list("p"=a$p,"Fobs"=a$"F-value","numDF"=a$numDF)
  Result
  } else
  {
  myCombination=as.matrix(myCombination)
  HB=t(myCombination)%*%myModel$coefficients[complete.cases(myModel$coefficients)]
  Fobs=t(HB)%*%solve(t(myCombination)%*%vcov(myModel)%*%myCombination)%*%HB
  p=1-pf(as.numeric(Fobs),1,myModel$df)
  a=c("p"=p,"Fobs"=as.numeric(Fobs),"numDF"=ncol(myCombination))
  Result<-list("p"=a[1],"Fobs"=a[2],"numDF"=a[3])
  Result
  }
}

RSquare<-function(mymodel)

{
  mydata=eval(mymodel$call[[3]])
  y=mydata[,as.character(formula(mymodel)[[2]])]
  if (class(mymodel)=="lme")
        {
        index=as.numeric(rownames(as.data.frame(mymodel$fitted)))
        r=cor(mymodel$fitted,y[index])
        }
  if (class(mymodel)=="gls")
        {
        index=as.numeric(rownames(as.data.frame(mymodel$fitted)))
        r=(cor(mymodel$fitted,y[index]))
        }
  if (class(mymodel)=="lm")
        {
        index=as.numeric(rownames(as.data.frame(mymodel$fitted)))
        r=(cor(mymodel$fitted,y[index]))
        }
        rr=as.data.frame(t(r*r))
        colnames(rr)=paste("R2",seq(0,length(rr)-1),sep="_")
        rr
}
#--------------------------------------------------------------------------