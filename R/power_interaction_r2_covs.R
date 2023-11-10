#' Creates the input to `power_interaction_r2_covs()`
#'
#' Companion function to 'power_interaction_r2_covs()'. Generates a formatted list for users to specify the analysis parameters.
#'
#' @param c.num Number of covariates in the model.
#'
#' @return A list to be used with the 'power_interaction_r2_covs()' function.
#' @export
#'
#' @examples
#' ex1 = generate.interaction.cov.input(c.num=2)
#' ex1$correlations$r.y.x1x2 = c(0.1,0.2,0.3)
generate.interaction.cov.input = function(c.num){

tot.var=4+c.num+c.num*2

cor.mat = base::matrix(data = 0,nrow = tot.var,ncol = tot.var) %>% base::as.data.frame()

base::diag(cor.mat) = 1

if(c.num>0){
var.names = c("y","x1","x2",
              base::paste("c",base::seq(1,c.num,1),sep = ""),
  "x1x2")}else{
    var.names = c("y","x1","x2",
                  #paste("c",seq(1:c.num),sep = ""),
                  "x1x2")

  }
if(c.num>0){
for(i in 1: c.num){
  var.names = c(var.names,
                base::paste("c",i,c("x1","x2"),sep="")
                )
}
}
base::colnames(cor.mat) = var.names
base::rownames(cor.mat) = var.names


ind <- base::which( base::lower.tri(cor.mat,diag=F) , arr.ind = TRUE )

cor.mat2 = base::data.frame( col = dimnames(cor.mat)[[2]][ind[,2]] ,
            row = dimnames(cor.mat)[[1]][ind[,1]] ,
            val = cor.mat[ ind ] )

cor.mat2$m1 = 0

for(i in 1: dim(cor.mat2)[1]){
  a1 = base::grep(pattern = cor.mat2$col[i],x = cor.mat2$row[i])
  a2 = base::grep(pattern = cor.mat2$row[i],x = cor.mat2$col[i])
    if(base::length(a1) > 0  | base::length(a2) > 0){  cor.mat2$m1[i] = 1}

  if(base::nchar(cor.mat2$col[i]) == 4  && base::nchar(cor.mat2$row[i]) == 4) {cor.mat2$m1[i] = 1}


}


fixed = cor.mat2[cor.mat2$m1 == 1,]
cor.mat2 = cor.mat2[cor.mat2$m1 == 0,]

cor.mat3 = base::matrix(0,nrow = dim(cor.mat2)[1],ncol = 1) %>% base::as.data.frame()
cor.mat3$V1 = base::paste("r",cor.mat2$col,cor.mat2$row,sep = ".")
cor.mat3 = t(cor.mat3) %>% base::as.data.frame()
base::colnames(cor.mat3) = cor.mat3[1,]
cor.mat3 = cor.mat3[-1,]
cor.mat3[1,] = 0
base::rownames(cor.mat3) =NULL

cor.mat3 = apply(cor.mat3, 2,base::as.numeric) %>% base::as.data.frame() %>% t()%>% base::as.data.frame() %>% base::as.list()


if(c.num>0){
  var.names2 = c("y","x1","x2",
                 base::paste("c",seq(1,c.num,1),sep = ""))}else{
                  var.names2 = c("y","x1","x2")

                }


rel.mat = base::matrix(0,nrow =base::length(var.names2),ncol = 1) %>% base::as.data.frame()
rel.mat$V1 = base::paste("rel",var.names2,sep = ".")
rel.mat = t(rel.mat) %>% base::as.data.frame()
colnames(rel.mat) = rel.mat[1,]
rel.mat = rel.mat[-1,]
rel.mat[1,] = 1
base::rownames(rel.mat) =NULL

rel.mat = apply(rel.mat, 2,base::as.numeric) %>% base::as.data.frame() %>% t()%>% base::as.data.frame() %>% as.list()

  out.mat=list()
  out.mat$correlations = cor.mat3
  out.mat$reliability = rel.mat
  return(out.mat)

}




#' Analytic interaction power analysis with covariates
#'
#' Analytic power analysis of an interaction model with covariates. Additional covariate x main effect interaction terms are additionally added.
#'
#' @param cov.input Output of 'power_interaction_r2_covs()'. Variable correlations and reliabilities are set by first modifying this list.
#' @param N Sample size. Must be a positive integer. Has no default value. Can be a single value or a vector of values.
#' @param alpha The alpha. At what p-value is the interaction deemed significant? Default is 0.05.
#' @param cl Number of clusters to use for running simulations in parallel. Default is NULL (i.e. not in parallel). Useful when running several thousand analyses at once.
#' @param detailed_results Default is FALSE. Should detailed results be reported?
#'
#' @importFrom dplyr "%>%"
#' @importFrom foreach "%dopar%"
#'
#' @return  A data frame containing the analytic power for each unique setting combination.
#' @export
#'
#' @examples
#' ex1 = generate.interaction.cov.input(c.num=2)
#' ex1$correlations$r.y.x1x2 = c(0.1,0.2,0.3)
#' power_interaction_r2_covs(cov.input = ex1,N=100)


power_interaction_r2_covs=function(cov.input, N,alpha=0.05,detailed_results = FALSE,cl=NULL){

c.num = base::length(cov.input$reliability)-3

# range checks
if(base::max(base::abs(base::expand.grid(c(cov.input$correlations)))) > 1){stop("Correlations cannot be greater than |1|")}
if(base::max((base::expand.grid(c(cov.input$reliability)))) > 1 | min((expand.grid(c(cov.input$reliability)))) <= 0){stop("Reliabilities cannot be greater than 1 or less than/equal to 0")}

cov.input$other$alpha=alpha
cov.input$other$N=N

full.dat = base::expand.grid(c(cov.input$correlations,cov.input$reliability,cov.input$other))

message(base::paste("Performing",(base::dim(full.dat)[1]) ,"analyses",sep=" "))


if(!is.null(cl)){
  clus <- parallel::makeCluster(cl)
  doParallel::registerDoParallel(clus)
  chunks<-  base::sample(base::rep(c(1:cl),length=base::dim(full.dat)[1]),replace = FALSE)
  chunks=chunks[base::order(chunks,decreasing = F)]

}else{
    chunks<-base::rep(1,dim(full.dat)[1])
    }


d=NULL
power_test<-foreach::foreach(d = 1:base::length(unique(chunks)),
                             .combine = 'rbind',
                             .packages = c('dplyr')) %dopar% {


invalid.settings=list()
power.out=list()

full.dat2 = full.dat[chunks == d,]

for(i in 1: dim(full.dat2)[1]){

    settings = full.dat2[i,]

  ####
  # create correlation data

  tot.var=4+c.num+c.num*2

  cor.mat = base::matrix(data = 0,nrow = tot.var,ncol = tot.var) %>% as.data.frame()

  base::diag(cor.mat) = 1

  if(c.num>0){
    var.names = c("y","x1","x2",
                  base::paste("c",seq(1,c.num,1),sep = ""),
                  "x1x2")}else{
                    var.names = c("y","x1","x2",
                                  #paste("c",seq(1:c.num),sep = ""),
                                  "x1x2")

                  }
  if(c.num>0){
    for(i in 1: c.num){
      var.names = c(var.names,
                    base::paste("c",i,c("x1","x2"),sep="")
      )
    }
  }
  base::colnames(cor.mat) = var.names
  base::rownames(cor.mat) = var.names


  ind <- base::which( base::lower.tri(cor.mat,diag=F) , arr.ind = TRUE )

  cor.mat2 = data.frame( col = base::dimnames(cor.mat)[[2]][ind[,2]] ,
                         row = base::dimnames(cor.mat)[[1]][ind[,1]]  )
  cor.mat2$V1 = base::paste("r",cor.mat2$col,cor.mat2$row,sep = ".")
  cor.mat2$index = base::as.numeric(rownames(cor.mat2))
  ###
  cor.settings = t(settings[,dplyr::starts_with(match = "r.",vars = base::colnames(settings))]) %>% base::as.data.frame()
  colnames(cor.settings)[1]="val"
  cor.settings$V1=base::rownames(cor.settings)

  cor.mat2=  base::merge(cor.mat2,cor.settings,all=T)
  cor.mat2=cor.mat2[base::order(cor.mat2$index,decreasing = F),]
  ###

  # create reliability data

  rel.settings = t(settings[,dplyr::starts_with(match = "rel.",vars = base::colnames(settings))]) %>% base::as.data.frame()
  colnames(rel.settings)[1]="rel"
  rel.settings$V1 =  base::as.data.frame(base::strsplit(rownames(rel.settings),"[.]"))[2,] %>% t() %>% base::unname()


 ###
  # adjust main correlations for reliability
for(j in 1: base::dim(cor.mat2)[1]){

  if(!is.na(base::match(cor.mat2$row[j],rel.settings$V1)) & !is.na(base::match(cor.mat2$col[j],rel.settings$V1))){

    j1=rel.settings$rel[base::match(cor.mat2$row[j],rel.settings$V1)]
    j2=rel.settings$rel[base::match(cor.mat2$col[j],rel.settings$V1)]

    cor.mat2$val[j] = cor.mat2$val[j]*sqrt(j1*j2)
    }
  }

###
  # compute reliability & sd of interaction terms

  if(c.num>0){
    int.vars = c("x1,x2")
    for(j in 1: c.num){
      int.vars = c(int.vars,
                   base::paste("c",j,",",c("x1","x2"),sep="")
      )
    }
  }else{ int.vars = c("x1,x2")}

  int.vars=  base::strsplit(int.vars,",") %>% base::as.data.frame() %>% t() %>% base::unname() %>% base::as.data.frame()

  colnames(int.vars) = c("a1","a2")
  int.vars$rel=NA
  int.vars$V1 = base::paste(int.vars$a1,int.vars$a2,sep="")

  int.vars$sd=NA

  for(j in 1: base::dim(int.vars)[1]){

  int.cor= cor.mat2$val[base::intersect(which(cor.mat2$col == int.vars$a1[j]), base::which(cor.mat2$row == int.vars$a2[j]))]

  if(base::length(int.cor) == 0){   int.cor= cor.mat2$val[base::intersect(which(cor.mat2$row == int.vars$a1[j]), base::which(cor.mat2$col == int.vars$a2[j]))]}

  j1=rel.settings$rel[base::match(int.vars$a1[j],rel.settings$V1)]
  j2=rel.settings$rel[base::match(int.vars$a2[j],rel.settings$V1)]

  int.vars$rel[j] = ((j1*j2)+int.cor^2)/(1 + int.cor^2)

  int.vars$sd[j] = sqrt(1 + int.cor^2)
  #int.vars$sd[j] = 1


 }

  rel.settings$sd=1

  int.vars = int.vars[,c(3:5)]
  rel.settings = base::rbind(rel.settings,int.vars)

###
   # adjust correlations btw. main and interactions for reliability
   # compute correlations between interaction terms

  for(j in 1: dim(cor.mat2)[1]){

    if(!is.na(cor.mat2$val[j])){
      if(base::nchar(cor.mat2$col[j]) == 4 | base::nchar(cor.mat2$row[j]) == 4 ){

        j1=rel.settings$rel[base::match(cor.mat2$row[j],rel.settings$V1)]
        j2=rel.settings$rel[base::match(cor.mat2$col[j],rel.settings$V1)]

        cor.mat2$val[j] = cor.mat2$val[j]*sqrt(j1*j2)

      }
    }


    if(is.na(cor.mat2$val[j])){
      a1 = base::grep(pattern = cor.mat2$col[j],x = cor.mat2$row[j])
      a2 = base::grep(pattern = cor.mat2$row[j],x = cor.mat2$col[j])
      if(base::length(a1) > 0  |base::length(a2) > 0){  cor.mat2$val[j] = 0}else{


        d1=base::paste(base::unlist(base::strsplit(cor.mat2$row[j],""))[c(1,2)],sep="",collapse="")
        d2=base::paste(base::unlist(base::strsplit(cor.mat2$row[j],""))[c(3,4)],sep="",collapse="")
        d3=base::paste(base::unlist(base::strsplit(cor.mat2$col[j],""))[c(1,2)],sep="",collapse="")
        d4=base::paste(base::unlist(base::strsplit(cor.mat2$col[j],""))[c(3,4)],sep="",collapse="")

      d1r=c(base::which(cor.mat2$col==d1),base::which(cor.mat2$row==d1))
      d2r=c(base::which(cor.mat2$col==d2),base::which(cor.mat2$row==d2))
      d3r=c(base::which(cor.mat2$col==d3),base::which(cor.mat2$row==d3))
      d4r=c(base::which(cor.mat2$col==d4),base::which(cor.mat2$row==d4))

      d12 = cor.mat2$val[base::intersect(d1r,d2r)]
      d34 = cor.mat2$val[base::intersect(d3r,d4r)]
      d13 = cor.mat2$val[base::intersect(d1r,d3r)]
      d24 = cor.mat2$val[base::intersect(d2r,d4r)]
      d14 = cor.mat2$val[base::intersect(d1r,d4r)]
      d32 = cor.mat2$val[base::intersect(d3r,d2r)]

      if(base::length(d12)>1){d12 = 1}
      if(base::length(d34)>1){d34 = 1}
      if(base::length(d13)>1){d13 = 1}
      if(base::length(d24)>1){d24 = 1}
      if(base::length(d14)>1){d14 = 1}
      if(base::length(d32)>1){d32 = 1}

      cor.mat2$val[j] =  (((d12*d34) + (d13*d24) + (d14*d32))-   (d12*d34))/ sqrt((rel.settings$sd[base::match(cor.mat2$row[j],rel.settings$V1)]^2)*(rel.settings$sd[base::match(cor.mat2$col[j],rel.settings$V1)]^2))

      }

    }
  }


  cor.mat[ind]=cor.mat2$val
  ind2=ind
 ind2[,1]=ind[,2]
 ind2[,2]=ind[,1]
 cor.mat[ind2]=cor.mat2$val

if(base::min(eigen(cor.mat,only.values = T)$values)<0){
invalid.settings[[length(invalid.settings)+1]]=settings
}else{

  rel.settings$sd1 = 1

  covmat_all<- base::diag(rel.settings$sd1) %*% base::as.matrix(cor.mat) %*% base::diag(rel.settings$sd1) # cor 2 cov

  cor.mat.no.int = base::as.matrix(cor.mat)[-base::which(rownames(cor.mat)=="x1x2"),-base::which(colnames(cor.mat)=="x1x2")]
  covmat.no.int <- base::diag(rel.settings$sd1[-base::which(rel.settings$V1 == "x1x2")]) %*% cor.mat.no.int %*% base::diag(rel.settings$sd1[-base::which(rel.settings$V1 == "x1x2")])# cor 2 cov


 ###  path tracing
  b<-covmat_all[1,-1] # the y-row
  A<-covmat_all[-1,-1] # everything but the y-row

  betas_all<-base::solve(A,b) # get the standardized regression effect sizes
  totalr2 = base::sum(betas_all*base::as.matrix(cor.mat)[1,-1])
  ##
  b<-covmat.no.int[1,-1] # the y-row
  A<-covmat.no.int[-1,-1] # everything but the y-row

  betas_noint<-base::solve(A,b) # get the standardized regression effect sizes
  nullr2 = base::sum(betas_noint*base::as.matrix(cor.mat.no.int)[1,-1])

### power

  df_denom <- (settings$N - (base::dim(rel.settings)[1]-2) )-1
  f2 <- (totalr2 -nullr2 )/(1-nullr2)
  lambda = f2*df_denom
  minusalpha<-1-settings$alpha
  Ft<-stats::qf(minusalpha, 1, df_denom)
  Power<-1-stats::pf(Ft, 1, df_denom,lambda)

  settings=cbind(Power,settings)
  colnames(settings)[1]="pwr"

  if(detailed_results == TRUE){

    covmat_all<- base::diag(rel.settings$sd) %*% base::as.matrix(cor.mat) %*% base::diag(rel.settings$sd) # cor 2 cov
    ###  path tracing
    b<-covmat_all[1,-1] # the y-row
    A<-covmat_all[-1,-1] # everything but the y-row
    betas_all<-base::solve(A,b) # get the standardized regression effect sizes

    details = c(f2,totalr2,nullr2,betas_all) %>% t() %>% base::as.data.frame()
    colnames(details)=c("f2","totalr2","nullr2",base::paste("b",rel.settings$V1[-1],sep=".") )
    settings=base::cbind(settings,details)

  }


  power.out[[base::length(power.out)+1]]=base::as.data.frame(settings)


}



  }##### end of fulldat loop

out.power=list()
out.power$power = power.out
out.power$error = invalid.settings
return(out.power)
} # end of dopar


if(!is.null(cl)){
power.out = power_test[,1]
invalid.settings = power_test[,2]

##########

power.out2=base::list()
for(X in 1:base::length(power.out)){
  power.out2[[X]] =   base::as.data.frame(base::do.call(base::rbind, power.out[[X]]))}

power.final = base::as.data.frame(base::do.call(base::rbind, power.out2))



#########

remove=0
for(k in 1:base::length(invalid.settings)){
  if(base::length(invalid.settings[[k]]) == 0){remove=c(remove,k)}
}
if(base::length(remove) >1){
  invalid.settings = invalid.settings[-remove]
  }

}else{
  power.out = power_test$power
  invalid.settings = power_test$error

  power.final = base::as.data.frame(base::do.call(base::rbind, power.out))

}



if(length(invalid.settings) > 0){

    invalid2=base::list()
  for(X in 1:base::length(invalid.settings)){
    invalid2[[X]] =   base::as.data.frame(base::do.call(base::rbind, invalid.settings[[X]]))}

  invalid.settings = base::as.data.frame(base::do.call(base::rbind, invalid2))

  warning(base::paste("N=",base::length(invalid.settings)),"power analyses could not be completed due to a correlation matrix that is not postive semi-definite.")
  warning(invalid2)
}



#power.final = as.data.frame(do.call(rbind, power.out))

if(!is.null(cl)){
  parallel::stopCluster(clus)
  foreach::registerDoSEQ()
}


return(power.final)
}

