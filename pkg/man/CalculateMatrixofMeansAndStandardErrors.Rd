\name{CalculateMatrixofMeansAndStandardErrors}
\alias{CalculateMatrixofMeansAndStandardErrors}
\title{CalculateMatrixofMeansAndStandardErrors}
\description{Returns the means and standard errors for each factor or factors combination included in formula argument of the call to lm, gls o mle functions.}
\usage{
CalculateMatrixofMeansAndStandardErrors(myModel, myData, myModelTerm, Mmatrix)
}
\arguments{
  \item{myModel}{
	An object of class lm, gls or mle.
}
  \item{myData}{
	A data frame included in the data argument of the lm, gls or mle functions.
}
  \item{myModelTerm}{
	myModelTerm is a string identifying the term of the model for which the means and estandar
errors are requested. These terms must not contain references to covariables and the factors involved had to be
declared in the model formula. Let suppose that the model formula is Y~FactorA+FactorB+FactorA:FactorB+CovariableX,
then the admissible calls to CalculateMatrixofMeansAndStandardErrors are:

CalculateMatrixofMeansAndStandardErrors(myModel, myData,"FactorA",Mmatrix)
CalculateMatrixofMeansAndStandardErrors(myModel, myData,"FactorB",Mmatrix)
CalculateMatrixofMeansAndStandardErrors(myModel, myData,"FactorA:FactorB",Mmatrix)

}
  \item{Mmatrix}{

	 A matrix returned by CalculateMmatrix(myModel,myData) function.
}
}
\details{

}
\value{

     A two columns matrix containing means and estandard error for the levels of the model term identified by myModelTerm.

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
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"Variety",Mmatrix)
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"nitro",Mmatrix)
CalculateMatrixofMeansAndStandardErrors(myModel,Oats,"Variety:nitro",Mmatrix)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
