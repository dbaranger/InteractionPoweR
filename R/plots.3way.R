#' See the correlation matrix for a 3-way interaction
#'
#' Prints or plots the correlation matrix for a 3-way interaction
#'
#' @param power.results Data frame of results from power_interaction_3way_r2().
#' @param row.num Which row to show? Can only be a single number. Default is 1.
#' @param return.plot Return a matrix (FALSE, default), or a plot (TRUE)?
#' @importFrom dplyr "%>%"

#' @return A matrix or a ggplot2 object
#' @export
#'
#' @examples
#' power_analysis = power_interaction_3way_r2(detailed_results = TRUE,N = c(1000),
#'r.x1.y = .2,r.x2.y = .3,r.x3.y = .1,r.x1x2.y = .01,r.x1x3.y = .05,r.x2x3.y = .1,
#'b.x1x2x3 = 0.1,r.x1.x2 = .1,r.x1.x3 = .1,r.x2.x3 = .1,
#'rel.x1 = 1,rel.x2 = 1,rel.x3 = 1,rel.y = 1 )
#'cor.mat.3way(power_analysis)

cor.mat.3way = function(power.results,row.num = 1,return.plot = FALSE){
  if(length(row.num)!=1){stop("row.num must be a single value")}
  cors = power.results[row.num,dplyr::starts_with(vars = colnames(power.results),match = "obs.r.")]


  cor.mat = base::matrix(data = 0,nrow = 8,ncol = 8) %>% base::as.data.frame()
  base::diag(cor.mat) = 1
  var.names = c("y","x1","x2","x3","x1x2","x1x3","x2x3","x1x2x3")
  base::colnames(cor.mat) = var.names
  base::rownames(cor.mat) = var.names
  ind <- base::which( base::lower.tri(cor.mat,diag=F) , arr.ind = TRUE )

  cor.mat[ind]=cors %>% c() %>% unlist()
  ind2=ind
  ind2[,1]=ind[,2]
  ind2[,2]=ind[,1]
  cor.mat[ind2]=cors%>% c() %>% unlist()

  if(return.plot == FALSE){out = cor.mat
  }else{

    ind = expand.grid(row=c(1:8),col=c(1:8)) %>% as.matrix()

    cor.mat2 = base::data.frame( col = dimnames(cor.mat)[[2]][ind[,2]] ,
                                 row = dimnames(cor.mat)[[1]][ind[,1]] ,
                                 val = cor.mat[ ind ] ) %>% as.data.frame()

    cor.mat2$col = factor(cor.mat2$col,levels = rev(var.names))
    cor.mat2$row = factor(cor.mat2$row,levels = (var.names))

    segment1 = data.frame(x = seq(0.5,8.5,1),xend = seq(0.5,8.5,1), y = .5, yend = 8.5)
    segment2 = data.frame(y = seq(0.5,8.5,1),yend = seq(0.5,8.5,1), x = .5, xend = 8.5)

    segment3=rbind(segment1,segment2)

    out =  ggplot2::ggplot(cor.mat2,ggplot2::aes(x = .data$row, y = .data$col, fill = .data$val,label = round(.data$val,3)))+
      ggplot2::geom_raster()+
      ggplot2::geom_segment(data = segment3,inherit.aes = F,ggplot2::aes(x = .data$x,xend = .data$xend, y = .data$y, yend = .data$yend),color="black")+
      ggplot2::scale_fill_gradient2(limits = c(-1,1),
                                    high = "#C42503FF",
                                    low = "#455BCDFF",
                                    mid = "white",
                                    midpoint = 0)+
      ggplot2::geom_text()+
      ggplot2::theme_classic()+
      ggplot2::coord_equal()+
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::labs(fill = "Correlation (r)")+
      ggplot2::theme(legend.position = "top",
            axis.title = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.text = ggplot2::element_text(color = "black"))+

      ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",ticks.colour = "black"))
  }
  return(out)
}


#' See the simple slopes for a 3-way interaction
#'
#' Prints or plots the simple slopes for a 3-way interaction
#'
#' @param power.results Data frame of results from power_interaction_3way_r2().
#' @param row.num Which row to show? Can only be a single number. Default is 1.
#' @param return.plot Return a matrix (FALSE, default), or a plot (TRUE)?
#' @importFrom dplyr "%>%"

#' @return A matrix or a ggplot2 object
#' @export
#'
#' @examples
#' power_analysis = power_interaction_3way_r2(detailed_results = TRUE,N = c(1000),
#'r.x1.y = .2,r.x2.y = .3,r.x3.y = .1,r.x1x2.y = .01,r.x1x3.y = .05,r.x2x3.y = .1,
#'b.x1x2x3 = 0.1,r.x1.x2 = .1,r.x1.x3 = .1,r.x2.x3 = .1,
#'rel.x1 = 1,rel.x2 = 1,rel.x3 = 1,rel.y = 1 )
#'simple.slopes.3way(power_analysis)

simple.slopes.3way = function(power.results,row.num = 1,return.plot = FALSE){


if(length(row.num)!=1){stop("row.num must be a single value")}
betas = power.results[row.num,dplyr::starts_with(vars = colnames(power.results),match = "obs.b.")]

row.names(betas)=NULL

slopes = expand.grid(M = c(-1,0,1),Z = c(-1,0,1),intercept1=NA,slope1=NA)

slopes = cbind(slopes,betas)

slopes$intercept1 = apply(slopes,1,function(X){
  X=as.data.frame(t(X))
  colnames(X) = colnames(slopes)
  (X$obs.b.x2*X$M) + (X$obs.b.x3*X$Z) + (X$obs.b.x2x3*X$M*X$Z)
})

slopes$slope1 = apply(slopes,1,function(X){
  X=as.data.frame(t(X))
  colnames(X) = colnames(slopes)
  X$obs.b.x1 + (X$obs.b.x1x2*X$M) + (X$obs.b.x1x3 * X$Z) + (X$obs.b.x1x2x3*X$M*X$Z)
})


if(return.plot == FALSE){

  slopes$X2 = paste(slopes$M,"SD")
  slopes$X3 = paste(slopes$Z,"SD")
  # slopes$X2 = stringr::str_replace(string = slopes$X2,pattern = "0 SD",replacement = "Mean")
  # slopes$X3 = stringr::str_replace(string = slopes$X3,pattern = "0 SD",replacement = "Mean")
  slopes$X2 = sub(x = slopes$X2,pattern = "0 SD",replacement = "Mean")
  slopes$X3 = sub(x = slopes$X3,pattern = "0 SD",replacement = "Mean")

  out = data.frame(X2 = slopes$X2,X3 = slopes$X3,intercept = slopes$intercept1,slope = slopes$slope1)

}

if(return.plot == TRUE){



  slopes$X2 = paste(slopes$M,"SD")
  slopes$X3 = paste("X3 =",slopes$Z,"SD")

  #slopes$X2 = stringr::str_replace(string = slopes$X2,pattern = "0 SD",replacement = "Mean")
  slopes$X2 = sub(x = slopes$X2,pattern = "0 SD",replacement = "Mean")
  #slopes$X3 = stringr::str_replace(string = slopes$X3,pattern = "0 SD",replacement = "Mean")
  slopes$X3 = sub(x = slopes$X3,pattern = "0 SD",replacement = "Mean")

  slopes$X2 = factor(slopes$X2,levels = c("-1 SD","Mean", "1 SD"))
  slopes$X3 = factor(slopes$X3,levels = c("X3 = -1 SD","X3 = Mean", "X3 = 1 SD"))


  out = ggplot2::ggplot(data = slopes)+

    ggplot2::geom_hline(yintercept = 0,color="grey",size=1,linetype="solid")+
    ggplot2::geom_vline(xintercept = 0,color="grey",size=1,linetype="solid")+

    ggplot2::geom_abline(ggplot2::aes(intercept = .data$intercept1,slope=.data$slope1,color = .data$X2),size=1.5,show.legend = F)+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=0,ymax=0,x=0,fill = .data$X2))+


    ggplot2::xlab(label = "X1 (SDs)")+
    ggplot2::ylab(label = "Y (SDs)")+
    ggplot2::scale_y_continuous(limits = c(-3,3),breaks = seq(-3,3,1))+
    ggplot2::scale_x_continuous(breaks = seq(-3,3,1))+
    ggplot2::coord_cartesian(xlim = c(-3,3),clip = 'off')+
    ggplot2::ggtitle(label = "Simple slopes plot")+

    ggplot2::theme_minimal()+
    ggplot2::theme(
      # plot.title =  ggplot2::element_text(size = 18),
      # axis.title =  ggplot2::element_text(size = 12),
      # axis.text = ggplot2::element_text(size = 10),
      axis.ticks = ggplot2::element_line(size=1,color="black"),
      axis.ticks.length=ggplot2::unit(0.2,"cm"),
      axis.line.x = ggplot2::element_line(color = "black",size=1),
      axis.line.y = ggplot2::element_line(color = "black",size=1),
      legend.position='top',
        # strip.text = ggplot2::element_text(size = 12),
        plot.margin = ggplot2::unit(c(1,4,1,1), "lines"),
      legend.key=ggplot2::element_rect(size=1),
      legend.key.size = ggplot2::unit(1, "lines"))+

    ggplot2::scale_color_viridis_d(option = "D",direction = -1)+
    ggplot2::scale_fill_viridis_d(option = "D",direction = -1)+

    ggplot2::guides(color=ggplot2::guide_legend(title.position = "top"  ))+

    ggplot2::facet_wrap(ncol = 3, ~ X3)
}


return(out)
}
