\name{contrastLC}
\alias{contrastLC}
\title{
	contrastLC
}
\description{
	This function tests a single or a set of linear combinations. 
}
\usage{
contrastLC(myModel, myCombination)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myModel}{
	An object of class lm, gls or mle.
}
  \item{myCombination}{
	A vector of lengh p (p=the number of parameters of the fixed part of the fitted model) containing  a set of coefficients.
        These coefficients define a linear combination of the parameters which expected value under the null hypotesis is zero.
        It also admit a k*p matrix of coefficients, each row representing a linear combination.
}
}
\details{
Let myParameters be the vector of model parameters. Then the null hypotesis is myCombination%*%myParameters=0
}
\value{
An objetc of class named-list containing the parts 
p: p-value for the null hypotesis, Fobs: Observed F-estatistics, numDF: 1 or the number of independent rows of myCombination. 
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
myCombination<-as.vector(c(0,1,-1,0,0,0))
contrastLC(myModel, myCombination)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
