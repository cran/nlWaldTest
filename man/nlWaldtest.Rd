\name{nlWaldtest}
\alias{nlWaldtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonlinear restriction(s) Wald test
}
\description{
Tests restriction(s) on model parameters of the form R(b)=q, where R is vector or scalar valued (non)linear function of b, the vector of model parameters, and q is numeric vector or scalar. Delta method is used for covariance matrix. Applicable after any model provided parameters estimates and their covariance matrix are available.
}
\usage{
nlWaldtest(obj = NULL, texts, rhss = NULL, coeff = NULL, 
           Vcov = NULL, df2 = NULL, x = NULL)
# Standard:
# nlWaldtest(obj, texts) # Chi square test
# nlWaldtest(obj, texts, df2 = T) # F test

# Force different covariance matrix:
# nlWaldtest(obj, texts, Vcov =  vcovHC(obj))

# If coef(obj) and vcov(obj) are not available
# nlWaldtest(texts = restrictions, coeff = vector, Vcov = matrix)

# Backward compatibility:
# nlWaldtest(obj, texts, rhss)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
model object of any class, for which \code{vcov.class(obj)} and 
\code{coef.class(obj)} methods are defined. If missing, both \code{coeff} and
\code{Vcov} should be inputted.
}
  \item{texts}{
left-side(s) of normalized restriction(s), R(b), as string or vector of strings.
Multiple restrictions can be inputted as a character vector or as a character,
separated by semicolon. Right-hand sides can be included either separated 
by "=", or substracted, e.g. \code{texts = "b[1]^b[2] = 1; b[3] = 2"}, or, the same, \code{texts = c("a[1]^a[2] - 1", "a[3] = 2")}; \code{b}'s should be 
numbered as in \code{coeff} vector.
}
  \item{rhss}{
right-side(s) of normalized restriction(s) as number or vector. Retained
mostly for backward compatibility. Set to zero(s), if missing.
}
  \item{coeff}{
vector of parameter estimates. If missing, it is set to \code{coef(obj)}
when available. It allows, for example, to test hypotheses in terms of marginal
effects and elasticities provided their covariance matrix is inputted.
}
  \item{Vcov}{
covariance matrix of parameters. If missing, it is set to \code{coef(obj)} when 
available. If \code{coeff} and/or \code{Vcov}  are inputed, theirs counterparts from \code{obj}  are superseded.
}
  \item{df2}{
defines the type of the test. By default, Chi square test is performed. To 
perfom F test one can use \code{df2 = T}, if a method for \code{\link[stats]{df.residual}}  is available. Otherwise, one could input \code{df2 = n}, where \code{n} is a natural number. \code{df2} is the denominator df in the F statistics. If \code{df2 = T} but \code{df.residuals(obj)} doesn't exist, Chi square test is forced, followed by a message.
}
  \item{x}{
number, or numeric vector. Provides a way to supply cumbersome coefficients
into restrictions,  e.g. \code{texts = "b[1]^x[1] = x[2]"},
\code{x = c(0.1234, 5.6789)} to test b[1]^0.1234 = 5.6789. Instead of "b", one can use any valid variable name excluding "x".  The "cumbersome" coefficients 
must be named only as x[i].}
}
\details{The test should be applicable after (almost) any regression-type model, estimated using cross-section, time series, or panel data. If there are no methods for \code{coef(obj)} and/or \code{vcov(obj)}, \code{coeff} and \code{Vcov} arguments should be inputted directly. To realize the delta-method, the function first tries to compute analytical derivatives using \code{\link[stats]{deriv}}. If failed, it computes numerical derivatives, calling \code{\link[stats]{numericDeriv}}.
}
\value{
an object of "htest" class.
}
\references{
Greene, W.H. (2011). Econometric Analysis, 7th edition. Upper Saddle River, NJ: Prentice Hall
}
\author{
Oleh Komashko
}

\seealso{
\code{\link{nlConfint}}
}

\examples{
set.seed(13)
x1<-rnorm(30);x2<-rnorm(30);x3<-rnorm(30);y<-rnorm(30)
set.seed(NULL)
lm1<-lm(y~x1+x2+x3)
nlConfint(lm1, "b[2]^3+b[3]*b[1];b[2]")
nlWaldtest(lm1,"a[2]^3+a[3]*a[1] = x[1]; a[2]", x = -0.07)
nlWaldtest(lm1,c("b[2]^3+b[3]*b[1]+0.07", "b[2]"))



# Reproduce example in EVievs 8 Users Guide II, pp. 149-151.

\dontrun{
require(nlme) 
nl1<-nls(log(q)~c1+c2*log(c3*(k^c4)+(1-c3)*(l^c4)),
data=CESdata,start=list(c1=-2.6,c2=1.8,c3=0.0001,c4=-6),
nls.control(maxiter = 100, tol = 1e-05,minFactor = 1/2^15))
nlWaldtest(nl1,"b[2]-1/b[4]",0)
nlWaldtest(nl1,"b[2]*b[4]",1)
}
}
