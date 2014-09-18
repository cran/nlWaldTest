\name{nlWaldtest}
\alias{nlWaldtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Nonlinear restriction(s) Wald test}
\description{Tests restriction(s) on regression parameters of the form R(b)=q, where R is vector or scalar valued (non)linear function of b, the vector of regression parameters, and q is numeric vector or scalar. Delta method is used for covariance matrix.
}
\usage{
nlWaldtest(obj, texts, rhss, Vcov = vcov(obj))
# nlWaldtest(obj, texts, rhss, vcovHC(obj))
# nlWaldtest(obj, texts, rhss)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
estimated model object of class "lm","gls","nls","Arima", or "ivreg".
}
  \item{texts}{left-side(s) of normalized restriction(s), R(b), as string or vector of strings, see Examples.
}
  \item{rhss}{right-side(s) of normalized restriction(s) as number or vector.
  }
  \item{Vcov}{covariance matrix of the model parameters. HC or HAC estimators are possible for "lm" and "ivreg" objects.
}
}
\details{Model parameters denoted as param[1],param[2],... . They are counted as they go in coef(obj). "param" can be any sequence of small latin characters, e.g. b[1],b[2],...,rumpelstiltskin[1],rumpelstiltskin[2],... . Be careful with "Arima" objects, as they bave ARMA terms at the top.
}
\value{1 by 4 data.frame with Wald F-statistics, two F-distribution parameters and p-value. Names are "F", "df1", "df2", and "Pr > F".
}
\references{
Greene, W.H. (2011). Econometric Analysis, 7th edition. Upper Saddle River, NJ: Prentice Hal
}
\author{
Oleh Komashko
}

\examples{
x1<-rnorm(20);x2<-rnorm(20);x3<-rnorm(20);y<-rnorm(20);z1<-rnorm(20);

l35<-lm(y~x1+x2+x3)
nlWaldtest(l35,"a[2]^3+a[3]*a[1]",0)

# Require {sandwich}
\dontrun{
nlWaldtest(l35,c("a[2]+a[3]^4","a[1]+a[2]^2"),c(1,0),vcovHC(l35))
}

# Require {forecast}
\dontrun{
ddd<-ts(data.frame(y,x1,x2,x3))
tl35<-tslm(y~x1+x2+x3,data=ddd)
nlWaldtest(tl35,"b[2]^3+b[3]*b[1]",0)
}

# Require {nlme}
\dontrun{
gl35<-gls(y~x1+x2+x3,corr = corARMA(p = 1, q = 0))
nlWaldtest(gl35,"ab[2]^3+ab[3]*ab[1]",0)
}

# Require {forecast}
\dontrun{
al35<-auto.arima(y,xreg=cbind(x1,x2,x3))
nlWaldtest(al35,"a[2]^3+a[3]*a[1]",0)
}

# Require {forecast}
\dontrun{
Arl35<-Arima(y,order=c(1,0,0),xreg=cbind(x1,x2,x3))
nlWaldtest(Arl35,"a[4]^3+a[5]*a[3]",0)
}

arl35<-arima(y,order=c(1,0,0),xreg=cbind(x1,x2,x3))
nlWaldtest(arl35,"a[4]^3+a[5]*a[3]",0)

# Require {AER}, {sandwich}
\dontrun{
il35<-ivreg(y~x1+x2+x3|x2+x3+z1)
nlWaldtest(il35,"a[2]^3+a[3]*a[1]",0,vcovHC(il35))
}

# Reproduce example in EVievs 8 Users Guide II, pp. 149-151.
# Require {nlme}
\dontrun{
nl1<-nls(log(q)~c1+c2*log(c3*(k^c4)+(1-c3)*(l^c4)),
data=CESdata,start=list(c1=-2.6,c2=1.8,c3=0.0001,c4=-6),
nls.control(maxiter = 100, tol = 1e-05,minFactor = 1/2^15))
nlWaldtest(nl1,"b[2]-1/b[4]",0)
nlWaldtest(nl1,"b[2]*b[4]",1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonlinear}
\keyword{restriction}
