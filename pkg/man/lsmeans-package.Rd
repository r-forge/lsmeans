\name{lsmeans-package}
\alias{lsmeans-package}
\alias{lsmeans}
\docType{package}
\title{lsmeans}
\description{Calculates the lsmeans and their variances for lm, gls, mle, glm, mer, lmerMod or glmerMod models
	It also contains a function to test contrasts and to test linear combinations.
	The present version only works with default parametrization.
}
\details{
\tabular{ll}{
Package: \tab lsmeans\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2012-10-06\cr
License: \tab GPL \cr
LazyLoad: \tab yes\cr
}


The package contains the following main functions: CalculateMmatrix, CalculateMatrixofMeansAndStandardErrors, MeanAndVariance, MeansComparison and RSquare.

}
\author{
	Julio A. Di Rienzo  <email: dirienzo@agro.unc.edu.ar> and Maria del Carmen Romero <email: mariadelc.romero@gmail.com>
	
	Maintainer: 
	Julio A. Di Rienzo  <email: dirienzo@agro.unc.edu.ar> and Maria del Carmen Romero <email: mariadelc.romero@gmail.com>
 }
\references{

}
\keyword{ package }
\seealso{
	\code{\link{CalculateMmatrix}}, 
	\code{\link{CalculateMatrixofMeansAndStandardErrors}},
	\code{\link{MeanAndVariance}},
	\code{\link{contrastLC}}, 
	\code{\link{Test_Contrast}},
	\code{\link{RSquare}}

}
\examples{
library(nlme)
library(lme4)
library(lsmeans)
Oats[,"nitro"]=as.factor(Oats[,"nitro"])

myModel<-lme(yield~1+Variety+nitro+Variety:nitro,random =~1|Block/Variety,,method="REML",data=Oats)
anova(myModel)
Mmatrix=CalculateMmatrix(myModel,Oats)
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"Variety",Mmatrix,MissingCells=FALSE)
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"nitro",Mmatrix,MissingCells=FALSE)
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"Variety:nitro",Mmatrix,MissingCells=FALSE)
MeanAndVariance("VarietyVictory", Mmatrix,myModel)
MeanAndVariance("VarietyVictory:nitro0.6", Mmatrix,myModel)
Test_Contrast(myModel,myData,"Variety",Mmatrix,c(1,-1,0),10)
contrastLC(myModel,c(0,0,0,1,-1,0,0,0,0,0,0,0))

RSquare(myModel)

}
