# Package lsmeans
# Version 2.0
#
# Calculates the lsmeans and their variances for lm, gls, lme or mer models.
# It also contains a function for mean comparisons and for contrasting linear
# combination of parameters in the myModel
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
#     TheModelCoefficients
#     TheModelCoefficientsCovar
#     TheModelTerms
#
# Examplification of lsmeans package.
# Load the package
# Let suppose that myData is a dataframe containing the columns:
# "y", "Treatment" and "Sex". "y" represent the response, "Treatment" a factor
# with levels T1,T2,T3 and "Sex" a factor wiht levels F and M.
# Let an object named "myModel" contains the fitted lm, gls, lme or mer model.
# Let suppose that the formula for this myModel is y~Treatment+Sex+Treatment:Sex
#
# The first step is to calculate:
# Mmatrix<-CalculateMmatrix(myModel,myData)
#(CalculateMmatrix calculates a matrix needed for other routines
# this function is computationally demanding. There is no need to call it
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
# It works in the same way as it does  CalculateMatrixofMeansAndStandardErrors
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
# For example let result be the vector of coefficients returned by MeanAndVariance,
# b the vector of estimated parameters of the myModel and S the
# variance-covariance matrix of b. Then mean=result'b and variance=result'Sx
#
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
# myModelTerm for the comparison and the function performs a test like simmilar to
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
# To obtain the R-square of the fitted myModel use RSquare(myModel). If the myModel
# is a mixed effects myModel this function will retur an R-square for the fixed
# and the nested random effects (the R-square will increase with the
# consecutive random effects)



 TheModelCoefficients<-function(myModel)
 {
  if (class(myModel)=="lme") result<-myModel$coefficients$fixed
  if (class(myModel)=="gls") result<-myModel$coefficients
  if (class(myModel)=="lm")  result<-myModel$coefficients[complete.cases(myModel$coefficients)]
  if (class(myModel)=="glm") result<-myModel$coefficients[complete.cases(myModel$coefficients)]
  if (class(myModel)=="mer") result<-myModel@fixef[complete.cases(myModel@fixef)]
  if (class(myModel)=="glmerMod") result<-fixef(myModel)
  if (class(myModel)=="lmerMod") result<-fixef(myModel)
  result
 }

TheModelCoefficientsCovar<-function(myModel)
{
   if (class(myModel)=="lme") result<-myModel$varFix
   if (class(myModel)=="gls") result<-myModel$varBeta
   if (class(myModel)=="lm")  result<-vcov(myModel)
   if (class(myModel)=="glm") result<-summary(myModel)$cov.unscale
   if (class(myModel)=="mer") result<-as.matrix(summary(myModel)@vcov)
   if (class(myModel)=="glmerMod") result<-as.matrix(summary(myModel)$vcov)
   if (class(myModel)=="lmerMod") result<-as.matrix(summary(myModel)$vcov)
 result
}

TheModelTerms<-function(myModel)
{
 if (class(myModel)=="lme")      result<-colnames(attributes(terms(myModel))$factors)
 if (class(myModel)=="gls")      result<-names(myModel$parAssign)
 if (class(myModel)=="lm")       result<-labels(myModel$terms)
 if (class(myModel)=="glm")      result<-labels(myModel$terms)
 if (class(myModel)=="mer")      result<-labels(terms(myModel))
 if (class(myModel)=="glmerMod") result<-labels(terms(myModel))
 if (class(myModel)=="lmerMod") result<-labels(terms(myModel))
result
}


SearchForASubstringInAStringWithSeparators<-function(s1,s2,caracterseparadordentrodes2)
 {
 myResult<-FALSE
 listastrings<-strsplit(s2,caracterseparadordentrodes2)[[1]]
 indices<-grep(s1,listastrings,fixed =TRUE)
 i<-0
 if (length(indices)>0)
{
 while ((myResult==FALSE)&(i<length(indices)))
     {
      i<-i+1
      myResult<-(myResult|(s1==listastrings[indices[i]]))
     }
 }
 myResult
 }
#---------------------------------------------------------------------------

 SearchForASubstringInAListOfStringsWithSeparators<-function(s1,myList,caracterseparadordentrodes2)
 {
 result<-grep(s1,c())
 indices<-grep(s1,myList,fixed =TRUE)
 if (length(indices)>0) for (i in (1:length(indices))) if (SearchForASubstringInAStringWithSeparators(s1,myList[indices[i]],caracterseparadordentrodes2)) result<-c(result,indices[i])
 result
 }

#---------------------------------------------------------------------------

 SearchForASubstringInAListOfStrings<-function(s1,myList,caracterseparadordentrodes2)
 {
 result<-grep(s1,c())
 indices<-grep(s1,myList,fixed =TRUE)
 if (length(indices)>0) for (i in (1:length(indices))) if (s1==myList[indices[i]]) result<-c(result,indices[i])
 result
 }
#---------------------------------------------------------------------------

CalculateMmatrix<-function(myModel,myData)
{
 #---------------------------------------------------------------------------

 for (i in (1:ncol(myData))) if (is.factor(myData[,i]))  myData[,i]<-as.factor(as.character(myData[,i]))

 f=TheModelTerms(myModel)

# f<-strsplit(f, ":")
# factoresEnELmodelo<- c()
# for (i in (1:length(f)))  factoresEnELmodelo<- c(factoresEnELmodelo,f[[i]])
# factoresEnELmodelo<- unique(factoresEnELmodelo)
# colnamedatos<-colnames(myData)
# factores<-c()
# for (i in (1:length(factoresEnELmodelo)))
# {
#  result<-SearchForASubstringInAListOfStringsWithSeparators(factoresEnELmodelo[i],colnamedatos,":")
# if (length(result)>0) {if (is.factor(myData[,result])==TRUE) factores<-c(factores,colnamedatos[result])}
# }
#
# myList=paste(rep(factores[1],nrow(myData)),myData[,factores[1]],sep='')
# if (length(factores)>1) for (i in (2:length(factores)))
#    {
#    myList=paste(myList,paste(rep(factores[i],nrow(myData)),myData[,factores[i]],sep=''),sep=':')
#    }
# myList=unique(myList)

# myModelFactors<-unique(unlist(strsplit(f, ":")))
# myModelFactors<-intersect(myModelFactors,colnames(myData))
# myModelFactors=myModelFactors[unlist(lapply(myData[,myModelFactors],is.factor))]
# indices=complete.cases(myData[,myModelFactors])
# aa=as.data.frame(sapply(myModelFactors,function(result){paste(result,myData[indices,result],sep="")}))
# myList=as.character(unique(interaction(aa[,myModelFactors],sep=":")))


 indices=grep(":",f)
 if (length(indices)>0) myModelFactors<-unique(f[-indices]) else myModelFactors<-unique(f)
 myModelFactors<-intersect(myModelFactors,colnames(myData))
 myModelFactors=myModelFactors[unlist(sapply(as.data.frame(myData[,myModelFactors]),is.factor))]
 LevelsList=lapply(myModelFactors,function(result,myData) {paste(result,levels(myData[,result]),sep="")},myData)
 list1=levels(interaction(LevelsList,sep=":"))

 myModelFactors<-unique(unlist(strsplit(f, ":")))
 myModelFactors<-intersect(myModelFactors,colnames(myData))
 myModelFactors=myModelFactors[unlist(lapply(as.data.frame(myData[,myModelFactors]),is.factor))]
 indices=complete.cases(myData[,myModelFactors])
 aa=as.data.frame(sapply(myModelFactors,function(result){paste(result,myData[indices,result],sep="")}))
 list2=as.character(unique(interaction(aa[,myModelFactors],sep=":")))
 myList=union(list1,list2)


 #---------------------------------------------------------------------------
 # in myList, the list of all treatments. In case of nested factors, only the
 # the convinations for the levels observed in the data
 #---------------------------------------------------------------------------


 myCoefficientsNames=names(TheModelCoefficients(myModel))

 myCoefficients<-matrix(rep(0,length(myList)*length(myCoefficientsNames)),nrow=length(myList),ncol=length(myCoefficientsNames))
 myConstantTerm<-as.integer(length(grep("(Intercept)",myCoefficientsNames,fixed =TRUE))>0)


 for (j in (1:length(myCoefficientsNames)))
     {
      ss<-unlist(strsplit(myCoefficientsNames[j], ":"))
      indices<-SearchForASubstringInAListOfStringsWithSeparators(ss,myList,':')
      if (length(indices)>0)
      {
      for (i in (1:length(indices)))
      {
      result<-SearchForASubstringInAListOfStringsWithSeparators(ss[1],myList[indices[i]],":")
      if (length(ss)>1) for (ii in (2:length(ss))) result<-intersect(result,SearchForASubstringInAListOfStringsWithSeparators(ss[ii],myList[indices[i]],":"))
      if (length(result)>0) myCoefficients[indices[i],j]<-1
      }
      }
     }
 if (myConstantTerm>0) myCoefficients[,1]<-myConstantTerm


 for (i in (1:length(colnames(myData))))
 if (is.numeric(myData[,i]))
 {
 g<-SearchForASubstringInAListOfStringsWithSeparators(colnames(myData)[i],myCoefficientsNames,":")
 if (length(g)>0) for (j in (1:length(g)))
        {
            sc=colnames(myData)[i]
            ss=strsplit(myCoefficientsNames[g[j]], ":")[[1]]
            ss=setdiff(ss,sc)
            if (length(ss)>1) ss=paste(ss,collapse=":")
            indices<-c()
            if (length(ss)>0) indices<-grep(ss,myList,fixed =TRUE)
            mx=mean(myData[,i],na.rm = TRUE)
            if   (length(indices)==0) myCoefficients[,g[j]]<-mx
            if   (length(indices)>0)
                 for (k in (1:length(indices))) myCoefficients[indices[k],g[j]]=mx
        }
 }
 rownames(myCoefficients)<-myList
 colnames(myCoefficients)<-myCoefficientsNames
 myCoefficients
 }
 #---------------------------------------------------------------------------

 linear.combination.to.estimate.a.mean<-function(MeanLabel,Mmatrix,myCoefficients)
 {
  ss<-strsplit(MeanLabel, ":")[[1]]
  result<-SearchForASubstringInAListOfStringsWithSeparators(ss[1],rownames(Mmatrix),":")
  if (length(ss)>1) for (i in (2:length(ss))) result<-intersect(result,SearchForASubstringInAListOfStringsWithSeparators(ss[i],rownames(Mmatrix),":"))

  ifelse (length(result)>0, SM<-Mmatrix[result,],SM<-c(1,rep(0,(length(myCoefficients)-1))))
  result<-as.matrix(rep((1/nrow(SM)),nrow(SM)))
  ifelse ((nrow(result)==0),result<-as.matrix(SM), result<-as.matrix(t(result)%*%SM))
  if (nrow(result)>1) result<-t(result)
  colnames(result)<-names(myCoefficients)
  result
 }


 MeanAndVariance<-function(MeanLabel,Mmatrix,myModel)
  {
  myCoefficients=TheModelCoefficients(myModel)
  covar=TheModelCoefficientsCovar(myModel)
  result<-(linear.combination.to.estimate.a.mean(MeanLabel,Mmatrix,myCoefficients))
  myResult<-c(result%*%myCoefficients,result%*%covar%*%t(result),result)
  myResult
  }

#--------------------------------------------------------------------------
 RebuildTreatmentNames<-function(myData,myModelTerm)
 {
 f<-strsplit(myModelTerm, ":")[[1]]
 indices=complete.cases(myData[,unlist(f)])
 myList=paste(rep(f[1],nrow(as.data.frame(myData[indices,]))),myData[indices,f[1]],sep='')
 if (length(f)>1) for (i in (2:length(f)))
    {
    myList=paste(myList,paste(rep(f[i],nrow(myData[indices,])),myData[indices,f[i]],sep=''),sep=':')
    }
 myList=unique(myList)
 myList
 }

#--------------------------------------------------------------------------
CalculateMatrixofMeansAndStandardErrors<-function(myModel,myData,myModelTerm,Mmatrix,MissingCells=FALSE)
{
 factores<-strsplit(myModelTerm, ":")
 f<-factores[[1]]
 myList=paste(rep(f[1],nrow(myData)),myData[,f[1]],sep='')
 if (length(f)>1) for (i in (2:length(f)))
    {
    if (MissingCells==TRUE) myList=levels(interaction(myList,paste(rep(f[i],nrow(myData)),myData[,f[i]],sep=''),sep=':')) else myList=paste(myList,paste(rep(f[i],nrow(myData)),myData[,f[i]],sep=''),sep=':')
    }
 if (MissingCells==TRUE) myList=as.character(myList)
 myList=unique(myList)
 myTable<-c()
 for (i in (1:length(myList))) myTable<-rbind(myTable,c(MeanAndVariance(myList[i],Mmatrix,myModel)[1],sqrt(MeanAndVariance(myList[i],Mmatrix,myModel)[2])))
 rownames(myTable)<-myList
 colnames(myTable)<-c("Media","E.E.")
 myTable
}

#---------------------------------------------------------------------------
FindGreatestDFForAGivenListOfMeansForA_lme_model<-function(myModel,meanslist)
{
  elgl<-function(myMean,myTerm)
     {
        result=(attributes(gregexpr (myTerm,myMean,fixed = T)[[1]])[[1]]==nchar(myTerm))
     }
  maxgl<-function(myterm,meanslist)
      {
      sapply(meanslist,elgl,myterm)
      }
   if (class(myModel)=="lme")
  {
   myterms=names(myModel$fixDF$terms)
   m=as.matrix(sapply(myterms,maxgl,meanslist))
   max(myModel$fixDF$terms%*%t(m))
  }
}
#--------------------------------------------------------------------------
AvgSEMeansDiff_and_AvgDMS<-function(myModel,myData,myModelTerm,Mmatrix,alfa=0.05)
{
 compare<-function(listofmeans,myModel)
 {

  COVBETA=TheModelCoefficientsCovar(myModel)

   result=c()
   for (i in (1:(length(listofmeans)-1)))
      {
      m1<-MeanAndVariance(listofmeans[i],Mmatrix,myModel)
      for (j in ((i+1):length(listofmeans)))
          {
          m2<-MeanAndVariance(listofmeans[j],Mmatrix,myModel)
          dif<-m1[3:length(m1)]-m2[3:length(m2)]
          result=c(result,(dif%*%COVBETA%*%dif))
         }
      }
    result
 }

 #-----------------------------------------
 factores<-strsplit(myModelTerm, ":")
 f<-factores[[1]]

 listofmeans=paste(rep(f[1],nrow(myData)),myData[,f[1]],sep='')
 if (length(f)>1) for (i in (2:length(f)))
    {
    listofmeans=paste(listofmeans,paste(rep(f[i],nrow(myData)),myData[,f[i]],sep=''),sep=':')
    }
 listofmeans=unique(listofmeans)
 #-----------------------------------------
  if (class(myModel)=="lm")   n=myModel$df.residual
  if (class(myModel)=="glm")  n=summary(myModel)$df.residual
  if (class(myModel)=="gls")  n=myModel$dims$N-myModel$dims$p
  if (class(myModel)=="lme")  n=FindGreatestDFForAGivenListOfMeansForA_lme_model(myModel,c(listofmeans[1],listofmeans[2]))
  if (class(myModel)=="mer")  n=myModel@dims["n"]-myModel@dims["p"]
  if (class(myModel)=="glmerMod")  n=nrow(getME(myModel,"X"))-ncol(getME(myModel,"X"))
  if (class(myModel)=="lmerMod")  n=nrow(getME(myModel,"X"))-ncol(getME(myModel,"X"))


  SEM=sqrt(mean(compare(listofmeans,myModel),na.rm=T))

 result=list(AvgSEM=SEM,AvgDMS=qt((1-alfa/2),n)*SEM)
 result
}

#---------------------------------------------------------------------------
MeansComparison<-function(MeanLabel1,MeanLabel2,Mmatrix,myModel,myDF=0)
{
  m1<-MeanAndVariance(MeanLabel1,Mmatrix,myModel)
  m2<-MeanAndVariance(MeanLabel2,Mmatrix,myModel)
  dif<-m1[3:length(m1)]-m2[3:length(m2)]
  
  myCoefficients=TheModelCoefficients(myModel)
  COVBETA=TheModelCoefficientsCovar(myModel)
  SEdiff<-sqrt(dif%*%COVBETA%*%dif)

  estadisticoT<-dif%*%myCoefficients/SEdiff

  if (class(myModel)=="lm")  myDF=myModel$df.residual
  if (class(myModel)=="glm") myDF=summary(myModel)$df.residual
  if ((class(myModel)=="lme")&(myDF==0)) myDF=FindGreatestDFForAGivenListOfMeansForA_lme_model(myModel,c(MeanLabel1,MeanLabel2))
  if (class(myModel)=="gls") myDF=myModel$dims$N-myModel$dims$p
  if (class(myModel)=="mer") myDF=myModel@dims["myDF"]-myModel@dims["p"]
  if (class(myModel)=="glmerMod")  myDF=nrow(getME(myModel,"X"))-ncol(getME(myModel,"X"))
  if (class(myModel)=="lmerMod")  myDF=nrow(getME(myModel,"X"))-ncol(getME(myModel,"X"))

  if (myDF==0)    ifelse((m1==m2),myResult<-1,myResult<-anova(myModel,L=dif)$p)
  if (!(myDF==0)) ifelse((m1==m2),myResult<-1,myResult<-2*(1-pt(abs(estadisticoT),df=myDF)))

  myResult
  }
#---------------------------------------------------------------------------

Test_Contrast<-function(myModel,myData,myModelTerm,Mmatrix,myContrast,denDF)
{
 factores<-strsplit(myModelTerm, ":")
 f<-factores[[1]]
 if (c((formula(myModel)[[2]])) %in% colnames(myData))  cc=complete.cases(myData[,c(unlist(factores),as.character(formula(myModel)[[2]]))]) else cc=complete.cases(myData[,c(unlist(factores))])
 rowdatos=nrow(as.data.frame(myData[cc,]))

 myList=paste(rep(f[1],rowdatos),myData[cc,f[1]],sep='')
 if (length(f)>1) for (i in (2:length(f)))
    {
    myList=paste(myList,paste(rep(f[i],rowdatos),myData[cc,f[i]],sep=''),sep=':')
    }
 myList=unique(myList)

 myTable<-c()
 for (i in (1:length(myList)))
 {
 a<-MeanAndVariance(myList[i],Mmatrix,myModel);
 myTable<-rbind(myTable,a[3:length(a)])
 }
 myTable=matrix(as.numeric(as.character(myTable)),ncol=ncol(myTable))
 dif<-as.matrix(t(myContrast)%*%myTable)

  myCoefficients=TheModelCoefficients(myModel)
  COVBETA=TheModelCoefficientsCovar(myModel)

  numDF<-nrow(dif)

  Fstatistic<-t(dif%*%myCoefficients)%*%solve(dif%*%COVBETA%*%t(dif))%*%(dif%*%myCoefficients)
  EE<-sqrt(diag(dif%*%COVBETA%*%t(dif)))


  if (class(myModel)=="glm") denDF=summary(myModel)$df.residual
  if (class(myModel)=="mer") denDF=myModel@dims["n"]-myModel@dims["p"]
  if (class(myModel)=="glmerMod")  n=nrow(getME(myModel,"X"))-ncol(getME(myModel,"X"))
  if (class(myModel)=="lmerMod")  n=nrow(getME(myModel,"X"))-ncol(getME(myModel,"X"))

  if (("glm" %in% class(myModel))|(class(myModel)=="mer")|(class(myModel)=="glmerMod")|(class(myModel)=="lmerMod"))
  {
          p<-(1-pchisq(Fstatistic[1],numDF))
          myResult<-list("Chi"=Fstatistic[1],"numDF"=numDF,"p"=p,"contrastvalue"=dif%*%myCoefficients,"contrastSE"=EE[1])
  }
  if (!(("glm" %in% class(myModel))|(class(myModel)=="mer")|(class(myModel)=="glmerMod")|(class(myModel)=="lmerMod")))
  {
          Fstatistic<-Fstatistic/numDF# es una chi dividida sus grados de libertad
          if (!(denDF==0)) ifelse((Fstatistic[1]==0),p<-1,p<-(1-pf(Fstatistic[1],numDF,denDF)))
          myResult<-list("F"=Fstatistic[1],"numDF"=numDF,"denDF"=denDF,"p"=p,"contrastvalue"=dif%*%myCoefficients,"contrastSE"=EE[1])
  }
  myResult

}
#---------------------------------------------------------------------------
checklinearindependence<-function(mycombination)
{
  mycombination=as.matrix(mycombination)
  result=FALSE
  tryCatch(error=function(e) "error", expr=
  {
   solve(t(mycombination)%*%mycombination)
   result=TRUE
  })
  result
}
#---------------------------------------------------------------------------
contrastLC<-function(myModel,mycombination)
{

 if (class(myModel)=="gls")
 {
  mycombination=as.matrix(mycombination)
  myanova=anova(myModel,L=mycombination)
  HB=t(mycombination)%*%myModel$coefficients[complete.cases(myModel$coefficients)]
  if (ncol(mycombination)==1) EE=sqrt(as.numeric(t(mycombination)%*%vcov(myModel)%*%mycombination)) else EE=NA
  myResult<-list("p"=myanova$p,"Fobs"=myanova$F,"numDF"=ncol(mycombination),"cl"=HB,"EE"=EE)
  }


 if (class(myModel)=="lme")
 {
 myResult=NULL
 mycombination=as.matrix(mycombination)
 tryCatch(error=function(e) "error", expr=
  {
  myanova=anova(myModel,L=mycombination)

  HB=t(mycombination)%*%myModel$coefficients$fixed[complete.cases(myModel$coefficients$fixed)]
  if (ncol(mycombination)==1) EE=sqrt(as.numeric(t(mycombination)%*%vcov(myModel)%*%mycombination)) else EE=NA
  myResult<-list("p"=myanova$p,"Fobs"=myanova$F,"numDF"=ncol(mycombination),"cl"=HB,"EE"=EE)
  })

  if (is.null(myResult))
          {
          dd=summary(myModel)$dims
          denDF=as.numeric(dd$N-sum(dd$ncol))
          HB=t(mycombination)%*%myModel$coefficients$fixed[complete.cases(myModel$coefficients$fixed)]
          if (ncol(mycombination)==1) EE=sqrt(as.numeric(t(mycombination)%*%vcov(myModel)%*%mycombination)) else EE=NA
          Fobs=t(HB)%*%solve(t(mycombination)%*%vcov(myModel)%*%mycombination)%*%HB
          p=1-pf(as.numeric(Fobs),ncol(mycombination),denDF)
          a=c("p"=p,"Fobs"=as.numeric(Fobs),"numDF"=ncol(mycombination))
          myResult<-list("p"=p,"Fobs"=as.numeric(Fobs),"numDF"=ncol(mycombination),"cl"=HB,"EE"=EE)
          }
  }

 if (class(myModel)=="lm")
  {
          mycombination=as.matrix(mycombination)
          denDF=myModel$df
          HB=t(mycombination)%*%myModel$coefficients[complete.cases(myModel$coefficients)]
          Fobs=t(HB)%*%solve(t(mycombination)%*%vcov(myModel)%*%mycombination)%*%HB
          if (ncol(mycombination)==1) EE=sqrt(as.numeric(t(mycombination)%*%vcov(myModel)%*%mycombination)) else EE=NA
          Fobs:=Fobs/ncol(mycombination)
          p=1-pf(as.numeric(Fobs),ncol(mycombination),denDF)
          a=c("p"=p,"Fobs"=as.numeric(Fobs),"numDF"=ncol(mycombination))
          myResult<-list("p"=a[1],"Fobs"=a[2],"numDF"=a[3],"cl"=HB,,"EE"=EE)
  }

 if (((class(myModel)=="glm") |(class(myModel)=="mer")|(class(myModel)=="glmerMod")|(class(myModel)=="lmerMod")))
  {
  mycombination=as.matrix(mycombination)

  if (class(myModel)=="glm")
    {
     HB=t(mycombination)%*%myModel$coefficients[complete.cases(myModel$coefficients)]
     if (ncol(mycombination)==1) EE=sqrt(as.numeric(t(mycombination)%*%summary(myModel)$cov.scaled%*%mycombination)) else EE=NA
     denDF=summary(myModel)$df.residual
     Fobs=t(HB)%*%solve(t(mycombination)%*%summary(myModel)$cov.scaled%*%mycombination)%*%HB
    }

  if (class(myModel)=="mer")
    {
     HB=t(mycombination)%*%myModel@fixef[complete.cases(myModel@fixef)]
     if (ncol(mycombination)==1) EE=sqrt(as.numeric(t(mycombination)%*%as.matrix(summary(myModel)@vcov)%*%mycombination)) else EE=NA
     Fobs=t(HB)%*%solve(t(mycombination)%*%as.matrix(summary(myModel)@vcov)%*%mycombination)%*%HB
    }

  if ((class(myModel)=="glmerMod")|(class(myModel)=="lmerMod"))
    {
     HB=t(mycombination)%*%fixef(myModel)[complete.cases(fixef(myModel))]
     if (ncol(mycombination)==1) EE=sqrt(as.numeric(t(mycombination)%*%as.matrix(vcov(myModel))%*%mycombination)) else EE=NA
     Fobs=t(HB)%*%solve(t(mycombination)%*%as.matrix(vcov(myModel))%*%mycombination)%*%HB
    }
  p=1-pchisq(as.numeric(Fobs),ncol(mycombination))
  a=c("p"=p,"Chi"=as.numeric(Fobs),"numDF"=ncol(mycombination))
  myResult<-list("p"=a[1],"Fobs"=a[2],"numDF"=a[3],"cl"=HB,"EE"=EE)
  }
  myResult
}

RSquare<-function(myModel)
{

  if ((class(myModel)=="mer")|(class(myModel)=="glmerMod")|(class(myModel)=="lmerMod")) y=myModel@resp$y

  if (!((class(myModel)=="mer")|(class(myModel)=="glmerMod")|(class(myModel)=="lmerMod")))
     {
     mydata=eval(myModel$call[[3]])
     y=mydata[,as.character(formula(myModel)[[2]])]
     }
  if (class(myModel)=="lme")
        {
        index=as.numeric(rownames(as.data.frame(myModel$fitted)))
        r=cor(myModel$fitted,y[index])
        }
  if (class(myModel)=="gls")
        {
        index=as.numeric(rownames(as.data.frame(myModel$fitted)))
        r=(cor(myModel$fitted,y[index]))
        }
  if (class(myModel)=="lm")
        {
        index=as.numeric(rownames(as.data.frame(myModel$fitted)))
        r=(cor(myModel$fitted,y[index]))
        }
  if (class(myModel)=="glm")
        {
        index=as.numeric(rownames(as.data.frame(myModel$fitted)))
        r=(cor(myModel$fitted,y[index]))
        }

  if ((class(myModel)=="mer")|(class(myModel)=="glmerMod")|(class(myModel)=="lmerMod"))
        {
        index=as.numeric(rownames(as.data.frame(myModel@resp$mu)))
        r=(cor(myModel@resp$family$linkinv(getME(myModel,"X")%*%fixef(myModel)),y[index]))
        }

        rr=as.data.frame(t(r*r))
        colnames(rr)=paste("R2",seq(0,length(rr)-1),sep="_")
        rr
}
#--------------------------------------------------------------------------

UpdateFormula<-function(myModel,myModelTerm){result=paste(as.character(formula(myModel))[[2]],as.character(formula(myModel))[[1]],as.character(formula(myModel))[[3]],"-",myModelTerm,sep="");result}

ReduceModel<-function(myModel,myModelTerm)
 {
  mynewModel=update(myModel,form=UpdateFormula(myModel,myModelTerm))
  result=anova(myModel,mynewModel)
  if (class(myModel)=="mer") result=list("UpdatedModel"=mynewModel,"Test"=c("Term"=myModelTerm,"Chi-square"=abs(2*result[,4]%*%c(-1,1)),"df"=result[2,6],"p-value"=1-pchisq(abs(2*result[,4]%*%c(-1,1)),result[2,6])))
  if ((class(myModel)=="glmerMod")|(class(myModel)=="lmerMod")) result=list("UpdatedModel"=mynewModel,"Test"=c("Term"=myModelTerm,"Chi-square"=result[2,6],"df"=result[2,7],"p-value"=result[2,8]))
  result
 }


#--------------------------------------------------------------------------
dfMixedModelTypeIAnova<-function(myModel,TheTitle="Building Anova Type I Table")
   {
         f=TheModelTerms(myModel)
         pb <- winProgressBar(title=TheTitle, min = 0,  max = length(f), width = 300)

         mynewModel=myModel
         result=c()
         for (i in (length(f):1))
            {
             r=ReduceModel(mynewModel,f[i])
             mynewModel=r$UpdatedModel;
             result=rbind(r$Test,result)
             setWinProgressBar(pb, (length(f)-i+1), title = TheTitle)
            }
        Sys.sleep(1)
        close(pb)
        result=as.data.frame(result)
        result[,2]=as.numeric(as.character(result[,2]))
        result[,3]=as.numeric(as.character(result[,3]))
        result[,4]=as.numeric(as.character(result[,4]))
        result
   }

#---------------------------------------------------------------------------
 RebuildFactorLevels<-function(myData,elfactor){unique(paste(elfactor,levels(myData[,elfactor]),sep=""))}
#---------------------------------------------------------------------------