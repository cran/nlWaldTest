\name{nlConfint}
\alias{nlConfint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence intervals for nonlinear functions of parameters
}
\description{
Computes confidence intervals for nonlinear functions of a model parameters. Delta method is used to compute standard errors. Applicable after any model provided  estimates of parameters and their covariance matrix are available.
}
\usage{
nlConfint(obj = NULL, texts, level = 0.95, coeff = NULL,
          Vcov = NULL, df2 = NULL, x = NULL)
# Standard:
# nlConfint(obj, texts) # based on z-statistics
# nlWaldtest(obj, texts, df2 = T) # based on z-statistics

# If coef(obj) and vcov(obj) are not available
# nlWaldtest(texts = funcions, coeff = vector, Vcov = matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
model object of any class, for which \code{vcov.class(obj)} and 
\code{coef.class(obj)} methods are defined. Otherwise, both \code{coeff} and
\code{Vcov} should be inputted directly.
}
  \item{texts}{
function(s) of parameters, b[i], as string or vector of strings.
Several functions can be inputted as a string, separated by semicolon, or as a
character vector, e.g. \code{texts = "b[1]^b[2]-1; b[3]"}, or \code{texts = c("b[1]^b[2]-1", "b[3]")}; \code{b}'s should be numbered as in \code{coeff} vector.
}
  \item{level}{
confidence level, a number in (0, 1). Default is 0.95.
}
  \item{coeff}{
vector of parameter estimates. If missing, it is set for \code{coef(obj)}
when available. It allows, for example, to compute CI for functions of marginal
effects and elasticities provided their covariance matrix is inputted.
}
  \item{Vcov}{
covariance matrix of parameters. If missing, it is set to \code{coef(obj)} when 
available. If \code{coeff} and/or \code{Vcov} are inputed, theirs counterparts from \code{obj}  are superseded.
}
  \item{df2}{
defines whether CI will be computed based on z (the default method) or t statistics. To compute t-based intervals, one can use \code{df2 = T}, provided a method for \code{\link[stats]{df.residual}}  is available. Otherwise, one could input \code{df2 = n}, where \code{n} is a natural number. \code{df2} is the  df in the t statistics. If \code{df2 = T} but \code{df.residuals(obj)} doesn't exist, z-based intervals are forced, followed by a message.
}
  \item{x}{
number, or numeric vector. Provides a way to supply cumbersome coefficients
into functions,  e.g. \code{texts = "b[1]^x[1] + x[2]"},
\code{x = c(0.1234, 5.6789)} to compute CI for b[1]^0.1234 + 5.6789.
}
}

\details{The function should be applicable after (almost) any regression-type model, estimated using cross-section, time series, or panel data. If there are no methods for \code{coef(obj)} and/or \code{vcov(obj)}, \code{coeff} and \code{Vcov} arguments should be inputted directly. To realize the delta-method, the function first tries to compute analytical derivatives using \code{\link[stats]{deriv}}. If failed, it computes numerical derivatives, calling \code{\link[stats]{numericDeriv}}.
}
\value{an r by 3 matrix, where r is the number of functions in \code{texts} argument. The first column is formed of values of the functions computed at
parameters estimates. The two last columns are confidence bounds.
}
\references{
Greene, W.H. (2011). Econometric Analysis, 7th edition. Upper Saddle River, NJ: Prentice Hall

}
\author{
Oleh Komashko
}

\seealso{
\code{\link{nlWaldtest}}
}
\examples{
set.seed(13)
x1<-rnorm(30);x2<-rnorm(30);x3<-rnorm(30);y<-rnorm(30)
set.seed(NULL)
lm1a<-lm(y~x1+x2+x3)
nlConfint(lm1a, c("b[2]^3+b[3]*b[1]","b[2]"))
}
