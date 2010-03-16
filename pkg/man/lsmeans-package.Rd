\name{lsmeans-package}
\alias{lsmeans-package}
\alias{lsmeans}
\docType{package}
\title{lsmeans}
\description{Calculates the lsmeans and their variances for lm, gls or lme models.
	It also contains a function for means comparisons and for contrasting linear
	combination of parameters.
	The present version only works with default parametrization.
}
\details{
\tabular{ll}{
Package: \tab lsmeans\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-03-03\cr
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
	\code{\link{MeansComparison}},
	\code{\link{RSquare}}


}
\examples{
library(nlme)
library(lsmeans)
Oats[,"nitro"]=as.factor(Oats[,"nitro"])

myModel<-lme(yield~1+Variety+nitro+Variety:nitro,random =~1|Block/Variety,,method="REML",data=Oats)
anova(myModel)
Mmatrix=CalculateMmatrix(myModel,Oats)
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"Variety",Mmatrix)
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"nitro",Mmatrix)
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"Variety:nitro",Mmatrix)
MeanAndVariance("VarietyVictory", Mmatrix,myModel)
MeanAndVariance("VarietyVictory:nitro0.6", Mmatrix,myModel)
MeansComparison("VarietyVictory:nitro0.6","VarietyVictory:nitro0.4",Mmatrix,myModel,myDF=10)
MeansComparison("VarietyVictory:nitro0.6","VarietyMarvellous:nitro0.4",Mmatrix,myModel,myDF=45)
RSquare(myModel)

}
