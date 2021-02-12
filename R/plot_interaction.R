
#' Plot interaction
#'
#' Plots a single simulated interaction data set
#'
#' @param data Output of generate_interaction().
#' @param q Simple slope quantiles. Default is 2. X2 is the moderator. Must be a positive integer > 1.
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' dataset <- generate_interaction(N = 250,r.x1.y = 0,r.x2.y = .1,r.x1x2.y = -.2,r.x1.x2 = .3)
#' plot_interaction(dataset,q=3)
#' }
plot_interaction<-function(data,q=2){

 # if(length(table(data$x2)) == 2){q<-2}
  #data$simpleslopes<-as.factor(as.numeric(gtools::quantcut(data$x2,q=q)))
  data$simpleslopes<-as.factor(as.numeric(cut(data$x2,
                              include.lowest = T,
                              breaks = stats::quantile(data$x2,probs = seq(0,1,by = 1/q)))))

  # #cols<-seq(0,1,by=1/(q+1))[-c(1, (q+2) )]
  # if(length(table(data$x1)) == 2 && length(table(data$x2)) == 2){
  #   data$x1<-as.factor(as.numeric(as.factor( data$x1 )))
  #   data$x2<-as.factor(as.numeric(as.factor( data$x2 )))
  #
  #   plot1<-ggplot2::ggplot(data = data,ggplot2::aes(x = x1,y=y,color = x2,fill=x2))+
  #     ggplot2::scale_color_viridis_d(option = c("D"))+
  #     ggplot2::scale_fill_viridis_d(option = c("D"))+
  #     ggplot2::geom_hline(yintercept = 0,color="grey",size=1,linetype="solid")+
  #
  #     ggplot2::geom_violin(alpha=0.25, position = ggplot2::position_dodge(width = .75),size=0,color=NA) +
  #     ggplot2::geom_boxplot(notch = T,  outlier.size = -1, color="black",lwd=1, alpha = 0.5,show.legend = F)+
  #
  #     ggplot2::geom_quasirandom(shape = 21,alpha=.5,dodge.width = .75,color="black",
  #                      width = .2,show.legend = F)+
  #
  #     ggplot2::theme_minimal()+
  #     ggplot2::xlab(label = "x")+
  #     # ggtitle(label = "Cross-over. rXM = 0.5")+
  #
  #     ggplot2::theme(
  #       axis.ticks = ggplot2::element_line(size=1,color="black"),
  #       axis.ticks.length=ggplot2::unit(0.2,"cm"),
  #       axis.line.x = ggplot2::element_line(color = "black",size=1),
  #       axis.line.y = ggplot2::element_line(color = "black",size=1),
  #       legend.position='top',
  #       legend.key=ggplot2::element_rect(size=.5),
  #       legend.key.size = ggplot2::unit(1, "lines"))+
  #     ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1,color="black")))+
  #     ggplot2::labs(fill="Moderator group")
  #
  #   #plot1<-ggpar(p = plot1,  legend.title = "Moderator group")
  #   return(plot1)
  #
  #
  # }else{

    plot1<-ggplot2::ggplot(data = data,ggplot2::aes(x = x1,y=y,color = simpleslopes,fill=simpleslopes))+
      ggplot2::scale_color_viridis_d(option = c("D"))+
      ggplot2::scale_fill_viridis_d(option = c("D"))+
      ggplot2::geom_hline(yintercept = 0,color="grey",size=1,linetype="solid")+
      ggplot2::geom_smooth(data = data,ggplot2::aes(x = x1,y=y),method = "stats::lm",inherit.aes = F,alpha=.5,
                  formula = y~x,se=F,size=1.5,color="black",linetype="dashed")+
      ggplot2::geom_point(alpha=0.5,shape=21,color="black")+
      ggplot2::geom_smooth(method = lm,formula = y~x,se=T,size=1.5,linetype="solid")+
      ggplot2::theme_minimal()+
      ggplot2::xlab(label = "x")+
      # ggtitle(label = "Cross-over. rXM = 0.5")+

      ggplot2::theme(
        axis.ticks = ggplot2::element_line(size=1,color="black"),
        axis.ticks.length=ggplot2::unit(0.2,"cm"),
        axis.line.x = ggplot2::element_line(color = "black",size=1),
        axis.line.y = ggplot2::element_line(color = "black",size=1),
        legend.position='top',
        legend.key=ggplot2::element_rect(size=1),
        legend.key.size = ggplot2::unit(1, "lines"))+
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1)))+
      ggplot2::labs(fill="Moderator group", color="Moderator group")

    #plot1<-ggpar(p = plot1,  legend.title = "Moderator quantile group")
    return(plot1)
  #}
}

