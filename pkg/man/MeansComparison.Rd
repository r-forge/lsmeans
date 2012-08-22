\name{MeansComparison}
\alias{MeansComparison}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	MeansComparison
}
\description{
It tests the equality of two means. 
}
\usage{
MeansComparison(MeanLabel1, MeanLabel2, Mmatrix, myModel, myDF = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MeanLabel1}{

String identifying a level of a factor or a combination of levels of interacting factors.
The factors involved had to be declared in the model formula. 
Let suppose that the model formula is Y~FactorA+FactorB+FactorA:FactorB+CovariableX,
if FactorA has three levels A1,A2 and A3, FactorB has two level B1 and B2, examples of 
admissible calls are:

MeansComparison("FactorAA1","FactorAA3", Mmatrix, myModel,myDF)
MeansComparison("FactorBB1","FactorBB2", Mmatrix, myModel,myDF)
MeansComparison("FactorAA1:FactorBB1","FactorAA1:FactorBB2", Mmatrix, myModel,myDF)

}
  \item{MeanLabel2}{
   See MeanLabel2
}
  \item{Mmatrix}{
	A matrix returned by CalculateMmatrix(myModel,myData)
}
  \item{myModel}{
	An object of class lm, gls or mle
}
  \item{myDF}{

It represents the degree of freedom of the appropriate error term for the comparison. 
These degrees of freedom can be obtained from tha anova table.
The function that performs a test is similar to a T-test.
If myDF is not specified or set to zero it uses the anova function for a linear combination of
the parameter.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	Returns the p-value for the comparison of the means of two levels of Treatment. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
	Julio A. Di Rienzo
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
library(lme4)
library(lsmeans)
Oats[,"nitro"]=as.factor(Oats[,"nitro"])

myModel<-lme(yield~1+Variety+nitro+Variety:nitro,random =~1|Block/Variety,,method="REML",data=Oats)
anova(myModel)
Mmatrix=CalculateMmatrix(myModel,Oats)
MeansComparison("VarietyVictory:nitro0.6","VarietyVictory:nitro0.4",Mmatrix,myModel,myDF=10)
MeansComparison("VarietyVictory:nitro0.6","VarietyMarvellous:nitro0.4",Mmatrix,myModel,myDF=45)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
