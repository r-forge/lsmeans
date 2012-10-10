\name{CalculateMatrixofMeansAndStandardErrors}
\alias{CalculateMatrixofMeansAndStandardErrors}
\title{CalculateMatrixofMeansAndStandardErrors}
\description{Returns the means and their standard errors for each factor or factors combination included in formula argument of the call to lm, gls, mle,lmer or glmer functions.}
\usage{
CalculateMatrixofMeansAndStandardErrors(myModel, myData,myModelTerm, Mmatrix,MissingCells)
}
\arguments{
  \item{myModel}{
	An object of class lm, gls, mle, glm, mer, lmerMod or glmerMod
}
  \item{myData}{
	A data frame included in the data argument of the call to the lm, gls, mle, glm, lmer and glmer functions.
}
  \item{myModelTerm}{
	myModelTerm is a string identifying the term of the model for which the means and estandar errors are requested. These terms must not contain references to covariables. Only the factors involved in the model formula can be included. Let suppose that the model formula is Y~FactorA+FactorB+FactorA:FactoB +CovariableX, then admissible calls to CalculateMatrixofMeansAndStandardErrors are:

CalculateMatrixofMeansAndStandardErrors(myModel, myData,"FactorA",Mmatrix,MissingCells)
CalculateMatrixofMeansAndStandardErrors(myModel, myData,"FactorB",Mmatrix,MissingCells)
CalculateMatrixofMeansAndStandardErrors(myModel, myData,"FactorA:FactorB",Mmatrix,MissingCells)

}
  \item{Mmatrix}{

	 A matrix returned by CalculateMmatrix(myModel,myData) function.
}
  \item{MissingCells}{

	 MissingCells=FALSE.
	 If TRUE it will return, when possible, estimations for means for levels combinations not present in the
data set. 
}

}
\details{

}
\value{

     A two columns matrix containing means and their estandard error for the levels of the model term identified by myModelTerm.

}
\references{

}
\author{
	Julio A. Di Rienzo  <email: dirienzo@agro.unc.edu.ar>
       }
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{MeanAndVariance}}
}
\examples{
library(nlme)
library(lme4)
library(lsmeans)
Oats[,"nitro"]=as.factor(Oats[,"nitro"])

myModel<-lme(yield~1+Variety+nitro+Variety:nitro,random =~1|Block/Variety,,method="REML",data=Oats)
Mmatrix=CalculateMmatrix(myModel,Oats)
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"Variety",Mmatrix,MissingCells=FALSE)
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"nitro",Mmatrix,MissingCells=FALSE)
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"Variety:nitro",Mmatrix,MissingCells=FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
