.numbersintext<-function(n,text)
{ h<-c();
  for (i in (1:n)) 
  {if (grepl(paste("[",as.character(i),"]",sep=""),text,fixed=TRUE)) h<-c(h,i)};
  h}
 .func1<-function(...){NULL}
.extractname<-function(text)
{ un<-unique(str_extract_all(text,"[a-z]+")[[1]])
  n<-length(un)
  h<-0
  for (i in (1:n)) 
  {if (grepl(paste(un[i],"[",sep=""),text,fixed=TRUE)) h<-un[i]}
  h}  

.mygrad<-function(text,values,nm="rumpelstiltskin2014")
{ if (nm=="rumpelstiltskin2014") name<-.extractname(text) else name<-nm;
  eval(parse(text=paste(".func1 <- function(",name,"){",text,"}")))  
  outp<-grad(.func1,values)
  remove(.func1)
  outp}

.myfunc<-function(text,values,nm="rumpelstiltskin2014")
{ if (nm=="rumpelstiltskin2014") name<-.extractname(text) else name<-nm;
  eval(parse(text=paste(".func1 <- function(",name,"){",text,"}")))  
  outp<-.func1(values);
  remove(.func1)
  outp}

.Grad1<-function(text,rhs,beta)
{ n<-length(beta)
  sub<-.numbersintext(n,text); 
  beta2<-rep(0,n)
  beta2[sub]<-beta[sub]
  m<-sub[length(sub)]
  grd<-.mygrad(text,beta2[1:m])
  grd2<-c(grd,rep(0,n-m))
  Rb<-.myfunc(text,beta2[1:m])-rhs
  list(grd2,Rb)} 

.Grads<-function(texts,rhss,beta)
{k<-length(rhss);
 Grads<-c()
 rhs<-c()
 for (i in 1:k) 
 {ll<-.Grad1(texts[i],rhss[i],beta);
 Grads<-rbind(Grads,ll[[1]]);
 rhs<-c(rhs,ll[[2]])} 
 list(Grads,rhs)}

.MyWald21<-function(obj,texts,rhss,Vcov=vcov(obj))
{ beta<-as.numeric(coef(obj));
  vc<-Vcov;rownames(vc)<-NULL;colnames(vc)<-NULL;
  Gradsrss<-.Grads(texts,rhss,beta);
  Grads<-Gradsrss[[1]];
  Rbs<-Gradsrss[[2]];
  rcm<-Grads%*%vc%*%t(Grads);
  ccm<-solve(rcm);
  n<-length(Rbs)
  Rbsl<-matrix(Rbs,1,n)
  Rbs2<-matrix(Rbs,n,1)
  Fst<- Rbsl%*%ccm%*%Rbs2 
  as.numeric(Fst)}

.MyWald31<-function(obj,texts,rhss,Vcov=vcov(obj))
{ beta<-as.numeric(coef(obj));
  vc<-Vcov;rownames(vc)<-NULL;colnames(vc)<-NULL;
  Gradsrss<-.Grad1(texts,rhss,beta);
  Grads<-Gradsrss[[1]];
  Rbs<-Gradsrss[[2]];
  n<-length(beta);
  gl<-matrix(Grads,1,n);
  gr<-matrix(Grads,n,1)
  rcm<-as.numeric(gl%*%vc%*%gr)
  Fst<-Rbs^2/rcm
  Fst}


nlWaldtest<-function(obj,texts,rhss,Vcov=
                       vcov(obj))
{ m<-length(rhss) 
  if (m==1) Fst<-.MyWald31(obj,texts,rhss,Vcov) else Fst<-.MyWald21(obj,texts,rhss,Vcov);
  df<-0
  if (class(obj)=="lm") df<-obj$df;
  if (class(obj)=="nls") df<-(summary(obj)$df)[2];
  if (class(obj)=="ivreg") df<-(summary(obj)$df)[2];
  if (class(obj)=="gls") df<-obj$dims$N-obj$dims$p
  if (class(obj)=="Arima") df<-length(obj$residuals)-length(obj$coef) 
    tabl<-data.frame(Fst,m,df,1-pf(Fst,m,df));
  colnames(tabl)<-c("F","df1","df2","Pr > F");
  rownames(tabl)<-""
  if (df==0) {print(paste("Be careful: has not tested for",class(obj),"yet.")); outp<-Fst} else outp<- tabl;
  # tabl
outp}