\name{lm1}
\alias{lm1}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("lm1")}
\format{
  The format is:
List of 12
 $ coefficients : Named num [1:2] -0.357 4.034
  ..- attr(*, "names")= chr [1:2] "(Intercept)" "Bwt"
 $ residuals    : Named num [1:144] -0.711 -0.311 1.789 -0.915 -0.815 ...
  ..- attr(*, "names")= chr [1:144] "1" "2" "3" "4" ...
 $ effects      : Named num [1:144] -127.567 23.411 1.87 -0.838 -0.738 ...
  ..- attr(*, "names")= chr [1:144] "(Intercept)" "Bwt" "" "" ...
 $ rank         : int 2
 $ fitted.values: Named num [1:144] 7.71 7.71 7.71 8.11 8.11 ...
  ..- attr(*, "names")= chr [1:144] "1" "2" "3" "4" ...
 $ assign       : int [1:2] 0 1
 $ qr           :List of 5
  ..$ qr   : num [1:144, 1:2] -12 0.0833 0.0833 0.0833 0.0833 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:144] "1" "2" "3" "4" ...
  .. .. ..$ : chr [1:2] "(Intercept)" "Bwt"
  .. ..- attr(*, "assign")= int [1:2] 0 1
  ..$ qraux: num [1:2] 1.08 1.12
  ..$ pivot: int [1:2] 1 2
  ..$ tol  : num 1e-07
  ..$ rank : int 2
  ..- attr(*, "class")= chr "qr"
 $ df.residual  : int 142
 $ xlevels      : Named list()
 $ call         : language lm(formula = Hwt ~ Bwt, data = cats)
 $ terms        :Classes 'terms', 'formula'  language Hwt ~ Bwt
  .. ..- attr(*, "variables")= language list(Hwt, Bwt)
  .. ..- attr(*, "factors")= int [1:2, 1] 0 1
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:2] "Hwt" "Bwt"
  .. .. .. ..$ : chr "Bwt"
  .. ..- attr(*, "term.labels")= chr "Bwt"
  .. ..- attr(*, "order")= int 1
  .. ..- attr(*, "intercept")= int 1
  .. ..- attr(*, "response")= int 1
  .. ..- attr(*, ".Environment")=<environment: 0x10ae9e990> 
  .. ..- attr(*, "predvars")= language list(Hwt, Bwt)
  .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "numeric"
  .. .. ..- attr(*, "names")= chr [1:2] "Hwt" "Bwt"
 $ model        :'data.frame':	144 obs. of  2 variables:
  ..$ Hwt: num [1:144] 7 7.4 9.5 7.2 7.3 7.6 8.1 8.2 8.3 8.5 ...
  ..$ Bwt: num [1:144] 2 2 2 2.1 2.1 2.1 2.1 2.1 2.1 2.1 ...
  ..- attr(*, "terms")=Classes 'terms', 'formula'  language Hwt ~ Bwt
  .. .. ..- attr(*, "variables")= language list(Hwt, Bwt)
  .. .. ..- attr(*, "factors")= int [1:2, 1] 0 1
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:2] "Hwt" "Bwt"
  .. .. .. .. ..$ : chr "Bwt"
  .. .. ..- attr(*, "term.labels")= chr "Bwt"
  .. .. ..- attr(*, "order")= int 1
  .. .. ..- attr(*, "intercept")= int 1
  .. .. ..- attr(*, "response")= int 1
  .. .. ..- attr(*, ".Environment")=<environment: 0x10ae9e990> 
  .. .. ..- attr(*, "predvars")= language list(Hwt, Bwt)
  .. .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "numeric"
  .. .. .. ..- attr(*, "names")= chr [1:2] "Hwt" "Bwt"
 - attr(*, "class")= chr "lm"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(lm1)
## maybe str(lm1) ; plot(lm1) ...
}
\keyword{datasets}
