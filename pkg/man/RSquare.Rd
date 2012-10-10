\name{RSquare}
\alias{RSquare}
\title{
	RSquare
}
\description{
It returns the R-square of a fitted model. 
If the model is a mixed effects model this function will return a named lists of R-squares.
The R2_0 for the fixed part and the R2_1, R2_2,... for the random parts. 
The R-squares will increase with the each nested random effect.
}
\usage{
RSquare(myModel)
}

\arguments{
  \item{myModel}{
	An object of class lm, gls, mle, glm, mer, lmerMod or glmerMod.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A named list of variable length containing R_squares.
At least it has one element R2_0, which is the R-square of the fixed part of the model,
then for each random effect an R-square is added with names R2_1,R2_2,...
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
library(lme4)
library(lsmeans)
Oats[,"nitro"]=as.factor(Oats[,"nitro"])

myModel<-lme(yield~1+Variety+nitro+Variety:nitro,random =~1|Block/Variety,,method="REML",data=Oats)
RSquare(myModel)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
