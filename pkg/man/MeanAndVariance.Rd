\name{MeanAndVariance}
\alias{MeanAndVariance}
\title{
	MeanAndVariance
}
\description{
Returns a vector containing the mean, the variance and the coefficients needed 
to combine the model's parameters to obtain the required mean and variance. 
}
\usage{
MeanAndVariance(tratamiento, Mmatrix, myModel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tratamiento}{
String identifying a level of a factor or a combination of levels of interacting factors.
The factors involved had to be declared in the model formula. 
Let Suppose that the model formula is Y~FactorA+FactorB+FactorA:FactorB+CovariableX,
if FactorA has three levels A1,A2 and A3, FactorB has two level B1 and B2, examples of 
admissible calls are:

MeanAndVariance("FactorAA1", Mmatrix, myModel)
MeanAndVariance("FactorBB1", Mmatrix, myModel)
MeanAndVariance("FactorAA1:FactorBB1", Mmatrix, myModel)


}
  \item{Mmatrix}{
	A matrix returned by CalculateMmatrix(myModel,myData)
}
  \item{myModel}{
	An object of class lm, gls or mle
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector containing the mean, the variance and a the coefficients needed to 
combine the model's parameters to obtain the  required mean and variance

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
	Julio A. Di Rienzo  <email: dirienzo@agro.unc.edu.ar>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(nlme)
library(lsmeans)
Oats[,"nitro"]=as.factor(Oats[,"nitro"])

myModel<-lme(yield~1+Variety+nitro+Variety:nitro,random =~1|Block/Variety,,method="REML",data=Oats)
Mmatrix=CalculateMmatrix(myModel,Oats)
MeanAndVariance("VarietyVictory", Mmatrix,myModel)
MeanAndVariance("VarietyVictory:nitro0.6", Mmatrix,myModel)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
