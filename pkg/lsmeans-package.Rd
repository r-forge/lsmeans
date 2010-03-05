\name{lsmeans-package}
\alias{lsmeans-package}
\alias{lsmeans}
\docType{package}
\title{lsmeans}
\description{Calculates the lsmeans and their variances for lm, gls or lme models
	It also contains a function for mean comparisons and for contrasting linear
	combination of parameters in the model
	The present version only works with default parametrization.
}
\details{
\tabular{ll}{
Package: \tab lsmeans\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-03-03\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}


The package contains the following main functions.\cr
\tabular{rlll}{
\tab \tab \tab CalculateMmatrix(myModel,myData,"Treatment",Mmatrix)\cr
\tab \tab \tab CalculateMatrixofMeansAndStandardErrors(myModel,myData,"Treatment",Mmatrix)\cr
\tab \tab \tab MeanAndVariance("TreatmentT1", Mmatrix,myModel)\cr
\tab \tab \tab MeansComparison("TreatmentT1,"TreatmentT2",Mmatrix,myModel,myDF)\cr
\tab \tab \tab contrastLC(myModel,myCombination)\cr
\tab \tab \tab RSquare(myModel)\cr
}
Load the package

 Let suppose that myData is a dataframe containing the columns:
 "y", "Treatment" and "Sex". "y" represent the response, "Treatment" a factor
 with levels T1,T2,T3 and "Sex" a factor wiht levels F and M.
 Let an object named "myModel" contains the fitted lm, gls or mle model
 let suppose that the formula for this model is y~Treatment+Sex+Treatment:Sex
 let an object named "myData" contains the data included in the "data"
 argument of the lm, gls or mle functions.

 Mmatrix<-CalculateMmatrix(myModel,myData) (CalculateMmatrix calculates a matrix needed for other routines)

 The means and standard errors for each level of "Treatment" are requested by

 CalculateMatrixofMeansAndStandardErrors(myModel,myData,"Treatment",Mmatrix)
 The means and estandard errors for each level of the interaction Treatment by
 Sex are requested by

 CalculateMatrixofMeansAndStandardErrors(myModel,myData,"Treatment:sex",Mmatrix)

 To obtain a single mean, its variance and the linear combination used to
 generate these estimates call MeanAndVariance(). For example if
 the mean of say level T2 of "Treatment" is required try the following
 call

 MeanAndVariance("TreatmentT2", Mmatrix,myModel)

 The function returns a vector containing the mean, the variance and a
 vector of coefficients used to calculate the mean and variance.
 For example let x be the vector of coefficients returned by MeanAndVariance,
 b the vector of estimated parameters of the model and S the
 variance-covariance matrix of b. Then mean=x'b and variance=x'Sx
 Other example:

 To obtain the mean of level T2 of factor Treatment and the level "F"
 factor sex make the following call

 MeanAndVariance("TreatmentT2:SexF", Mmatrix,myModel)

 To compare the means of levels T1 and T2 of "Treatment"  try

 MeansComparison("TreatmentT1,"TreatmentT2",Mmatrix,myModel,myDF)

 This function returns the p-value for the comparison. the argument myDF is
 optional. If specified, it represents the degree of freedom o the appropriate
 term for the comparison and the function performs a test like simmilar to
 a T-test. Otherwise it uses the anova function for a linear combination of
 the parameter. Usually you will not specify this argument or set it to zero

 To test a single or a set of linear combinations use

 contrastLC(myModel,myCombination)

 The function returns $p=p-value, $Fobs= the calculated F-value, and $numDF=
 the numerator degree of freedom of the F-test.
 The argument myCombination can be a single vector or a matrix of combinations,
 each combination corresponding to a row.

 To obtain the R-square of the fitted model use 
 
 RSquare(myModel). 

 If the model is a mixed effects model this function will retur an R-square for the fixed
 and the nested random effects (the R-square will increase with the
 consecutive random effects)

}


\author{
Julio A. Di Rienzo <dirienzo@agro.unc.edu.ar>\cr
Maria del Carmen Romero <mariadelc.romero@gmail.com>\cr\cr
Maintainer: Julio A. Di Rienzo <dirienzo@agro.unc.edu.ar> and Maria del Carmen Romero <mariadelc.romero@gmail.com>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
 
}
