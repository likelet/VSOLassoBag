# for simulating data
options(digits = 10)

SimulateData <- function(n.sample=1000,n.feature=100,n.fun=50,n.poisson=2,n.gaussian=3,n.binom=4,
                         poisson.lambda=8,random.range=1000,gaussian.mean.max=100,
                         gaussian.sd.max=10,random.seed=NULL,out.type="regression",
                         fun.ratio.for.y=NULL,bias=0.5,add.noise=T) {

  # n.feature is the number of features we want to simulate
  # n.fun is the number of funtion
  # n.poisson and n.othe is the number of features that with the distribution of
  # poisson and other. The sum of these number must be less than n.feature and the rest
  # will be occupied by random number(uniform distribution)
  # a very important things: Dependent Var is the linear(or other connection) combination of
  # features, I set default: randomly pick up some features and combine them
  # up to now this function just support regression
  # more:
  # poisson.lambda is the parameter lambda of poisson distribution: how many events happened during given time;
  # random.range is the upper bound of creating random distribution
  # gaussian.mean.max is the upper bound of mean that is set to normal distribution
  # gaussian.sd.max is of the same meaning
  # random.seed, as name suggested
  # out.type now can just for regression data, but not for classification data
  # fun.ratio.for.y should be a combination of the probability of each function to be selected
  # bias is the bias add the linear combination of X
  # add.noise is whether to add "gaussian" noise to the output(that is y)
  # returns:
  # x,y,information dataframe

  if (is.null(random.seed)) {
    set.seed(sample(c(1:1000),1))
  } else {
    set.seed(random.seed)
  }


  # init X
  X <- matrix(0,nrow=n.sample,ncol=n.feature,dimnames = list(c(1:n.sample),c(1:n.feature)))

  n.feature.useless <- n.feature - sum(c(n.poisson,n.gaussian,n.binom))

  # first let us define some distribution
  poisson.ft <- rep("poisson",n.poisson)
  gaussian.ft <- rep("gaussian",n.gaussian)
  binom.ft <- rep("binomial",n.binom)
  useless.ft <- rep("random",n.feature.useless)
  all.ft <- c(poisson.ft,gaussian.ft,binom.ft,useless.ft)
  # try to get sample probs
  fillin.probs <- runif(sum(c(n.poisson,n.gaussian,n.binom)))
  all.ft <- sample(all.ft,size=length(all.ft),replace=F)

  # apply distribution to each feature and sample from them
  col.distribution <- function(distribution,len) {
    if (distribution == "random") {
      # randomly set min and max
      mini <- sample(c((-random.range):random.range),1)
      maxi <- mini + sample(c(1:random.range),1)
      less.precise <- runif(n=len,min=mini,max=maxi)
      out <- round(less.precise,digits = 3)
      return(out)
    } else if (distribution=="binomial") {
      probs <- runif(1,min=0,max=1)
      return(rbinom(n=len,size=1,prob=probs))  # actually in R, rbinom can even create number of multiclassification
    } else if (distribution=="poisson") {
      return(rpois(n=len,lambda = poisson.lambda))
    } else if (distribution=="gaussian") {
      rand.mean <- runif(1,min=-gaussian.mean.max,max=gaussian.mean.max)
      rand.sd <- runif(1,min=0,max=gaussian.sd.max)
      less.precise <- rnorm(n=len,mean=rand.mean,sd=rand.sd)
      out <- round(less.precise,digits = 3)
      return(out)
    }
  }
  # substitude the origin X by different distribution
  for (i in c(1:n.feature)) {
    distribution <- all.ft[i]
    len <- n.sample
    this.col <- col.distribution(distribution,len)
    X[,i] <- this.col
  }

  # second, construct y. y is the linear/sigmoid combination of selected features, which is chosen
  # according to the pick up possibility
  # (1) init a function pool;(2) define the relationship of each features;(3) linear combination of random sampled function
  # we set each function contains 5 features(randomly sampled)
  fun.name.list <- c()
  ft.name.list <- c()

  for (i in c(1:n.fun)) {
    fun.name <- paste(c("fun",i),sep="",collapse="")
    ft.select <- sample(c(1:n.feature),5)
    ft.name <- paste(ft.select,sep="",collapse = ",")
    fun.name.list <- c(fun.name.list,fun.name)
    ft.name.list <- c(ft.name.list,ft.name)
  }
  if (is.null(fun.ratio.for.y)) {
    affect.ratio.list <- runif(length(fun.name.list),min=0,max=1)  # adjustable, we may change it to meet our requirement
  } else {
    affect.ratio.list <- fun.ratio.for.y
  }

  fun.ft.df <- data.frame(FunName=fun.name.list,FtName=ft.name.list,AffectRatio=affect.ratio.list)

  # the coefficient(from 1-10) in each function (default the function is linear function)
  coef.list <- c()
  for (f in fun.name.list) {
    # coef.sample <- runif(5,min=0,max=10)  # a more complicated form
    coef.sample <- rep(1,5)  # a simpler form
    coef.info <- paste(coef.sample,sep="",collapse = ",")
    coef.list <- c(coef.list,coef.info)
  }
  fun.ft.df <- cbind(fun.ft.df,Coef=coef.list)

  # make output,the output of y is the sum of the functions that selected
  pt.y.list <- c()
  pt.func.list <- c()  # for saving the info
  for (pt in c(1:n.sample)) {  # pt means patient, it's just a representation
    pt.x <- X[pt,]
    pt.ratio <- runif(1,min=0,max=1)
    info.need <- fun.ft.df[which(fun.ft.df$AffectRatio>=pt.ratio),]
    pt.func <- paste(info.need$FunName,sep="",collapse=",")
    pt.func.list <- c(pt.func.list,pt.func)  # for saving the info
    result.list <- c()
    for (func in info.need$FunName) {
      ft.thisfunc <- as.character(info.need[which(info.need$FunName==func),"FtName"])
      coef.thisfunc <- as.character(info.need[which(info.need$FunName==func),"Coef"])
      ft.thisfunc <- strsplit(ft.thisfunc,",")[[1]]
      coef.thisfunc <- strsplit(coef.thisfunc,",")[[1]]
      ft.thisfunc <- as.numeric(ft.thisfunc)
      coef.thisfunc <- as.numeric(coef.thisfunc)
      ft.value <- as.numeric(X[pt,ft.thisfunc])
      result <- ft.value * coef.thisfunc
      result.list <- c(result.list,result)
    }
    pt.y <- sum(result.list)  # ע??С????
    pt.y.list <- c(pt.y.list,pt.y)
  }

  df.pt.use.func <- data.frame(Patient=c(1:n.sample),UseFunc=pt.func.list)
  if (out.type=="regression") {
    output.info <- list(X,pt.y.list,fun.ft.df,df.pt.use.func)
    return(output.info)
  } else {
    pt.y.list.class <- c()
    med <- median(pt.y.list)
    for (y in pt.y.list) {
      if (y >= med) {
        pt.y.list.class <- c(pt.y.list.class,1)
      } else {
        pt.y.list.class <- c(pt.y.list.class,0)
      }
    }
    output.info <- list(X,pt.y.list.class,fun.ft.df,df.pt.use.func)
    return(output.info)
  }
}

# e.g.
m<-SimulateData(n.sample = 2500,n.feature=400,n.fun=3,n.poisson = 50,n.binom = 40,n.gaussian = 10,out.type = "classification")
mx.write <- as.data.frame(m[[1]])
my.write <- as.data.frame(m[[2]])
write.table(mx.write,file="/Users/tianqinli/Code/Renlab/lassoBag/data/SimulatedXv2.txt",
            quote=F,row.names = F)
write.table(my.write,file="/Users/tianqinli/Code/Renlab/lassoBag/data/SimulatedYv2.txt",
            quote=F,col.names = F,row.names=F)

# try <- read.table("/Users/tianqinli/Code/Renlab/lassoBag/data/SimulatedX.txt",
#                   header=T)
write.table(m[[3]],file="SavedFtFunc.txt",quote = F,row.names = F)
write.table(m[[4]],file="/Users/tianqinli/Code/Renlab/lassoBag/data/SavedUseFunc.txt",quote=F,row.names = F)
