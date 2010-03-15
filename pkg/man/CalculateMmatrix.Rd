\name{CalculateMmatrix}
\alias{CalculateMmatrix}
\title{CalculateMmatrix}
\description{
Calculates a matrix needed for other routines in the lsmeans package. 
This function is computationally demanding and there is no need to call it more than once for a given model.
}
\usage{
CalculateMmatrix(myModel, myData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myModel}{
	An object of class lm, gls or mle.
}
  \item{myData}{
	A data frame included in the data argument of the lm, gls or mle functions.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 An object of class matrix.  This matrix is required as an argument of the following function of the lsmeans package:

 CalculateMatrixofMeansAndStandardErrors()
 MeanAndVariance()
 MeansComparison()

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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
