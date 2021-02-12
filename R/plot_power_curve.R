#' Plot power curve
#'
#' Plot the output of power_interaction().
#'
#' @param power_data Data frame of results from power_interaction(). Can accept the raw results if up to 3 parameters were varied during simulation. Any more and data should be filtered first.
#' @param x Optional, the x-axis of the plot. Default is the first variable after 'pwr'.
#' @param group Optional, grouping variable for the line color. Default is the second variable after 'pwr', if present.
#' @param facets Optional, grouping variable for plot facets. Default is the third variable after 'pwr' if present.
#' @param power_target The target power. Default is 80%.
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' #'\dontrun{
#' dataset <- generate_interaction(N = seq(100,300,by=10),r.x1.y = 0,r.x2.y = .1,r.x1x2.y = -.2,r.x1.x2 = .3)
#' power_curve(dataset)
#' }
plot_power_curve<-function(power_data,x=NULL,group=NULL,facets=NULL,power_target=.8){

  y = "pwr"

  # Default behavior is to use columns in the order they are present in the power analysis output
  if(is.null(x) && is.null(group) && is.null(facets)){
    power_data2<-power_data[,c(1: which(colnames(power_data) == "pwr"))]
    x = colnames(power_data2)[1]

    if(dim(power_data2)[2]>2){ group = colnames(power_data2)[2] }
    if(dim(power_data2)[2]>3){ facets = colnames(power_data2)[3] }
    if(dim(power_data2)[2]>4){
      print("Too many variable combinations in input data. Please select a subset of data to plot")
      stop()
    }
  }

  y_col<-which(colnames(power_data) == y)
  x_col<-which(colnames(power_data) == x)
  group_col<-which(colnames(power_data) == group)
  facets_col<-which(colnames(power_data) == facets)
  col_keep<- c(unname(cbind(x_col,y_col,group_col,facets_col)))
  power_data<-as.data.frame(power_data[,col_keep])



  if(dim(power_data)[2] == 2){

    x_name<-InteractionPoweR::name_key[match(colnames(power_data)[1],InteractionPoweR::name_key[,1]  ),2]
    y_name<-InteractionPoweR::name_key[match(colnames(power_data)[2],InteractionPoweR::name_key[,1]  ),2]
    # group_name<-InteractionPoweR::name_key[match(colnames(power_data)[3],InteractionPoweR::name_key[,1] ),2]

    power_plot<-ggplot2::ggplot(data = power_data,ggplot2::aes(x = power_data[,1],y = power_data[,2] ))+
      #scale_color_viridis_d(option = c("C"),end = .95)+
      #scale_fill_viridis_d(option = c("C"),end = .95)+
      ggplot2::geom_hline(yintercept = power_target,color="black")+
      ggplot2::geom_hline(yintercept = 0.05,color="darkgrey")+
      ggplot2::geom_hline(yintercept = 1,color="darkgrey")+
      ggplot2::geom_line(linetype="solid",size=.5,color = "darkblue")+
      ggplot2::geom_smooth(formula = y ~ x,method = 'loess',se=T,alpha=0,size=1,
                  color = "darkblue")+
      ggplot2::xlab(label = x_name)+
      ggplot2::ylab(label = y_name)+
      ggplot2::geom_point(shape=21,color="black",show.legend = F,size=2,fill="white")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.key=ggplot2::element_rect(size=1,color="black"))+
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1)))

    return(power_plot)

  }


  if(dim(power_data)[2] == 3){

    x_name<-InteractionPoweR::name_key[match(colnames(power_data)[1],InteractionPoweR::name_key[,1]  ),2]
    y_name<-InteractionPoweR::name_key[match(colnames(power_data)[2],InteractionPoweR::name_key[,1]  ),2]
    group_name<-InteractionPoweR::name_key[match(colnames(power_data)[3],InteractionPoweR::name_key[,1] ),2]

    power_plot<-ggplot2::ggplot(data = power_data,ggplot2::aes(x = power_data[,1],
                                                               y = power_data[,2],
                                             color=as.factor(power_data[,3]),
                                             fill=as.factor(power_data[,3]) ))+
      ggplot2::scale_color_viridis_d(option = c("C"),end = .95)+
      ggplot2::scale_fill_viridis_d(option = c("C"),end = .95)+
      ggplot2::geom_hline(yintercept = power_target,color="black")+
      ggplot2::geom_hline(yintercept = 0.05,color="darkgrey")+
      ggplot2::geom_hline(yintercept = 1,color="darkgrey")+
      ggplot2::geom_line(linetype="solid",size=.5,show.legend = F)+
      ggplot2::geom_smooth(formula = y ~ x,method = 'loess',se=T,alpha=0,size=1)+
      ggplot2::xlab(label = x_name)+
      ggplot2::ylab(label = y_name)+
      ggplot2::geom_point(shape=21,color="black",show.legend = F)+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.key=ggplot2::element_rect(size=1,color="black"))+
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1)))+
      ggplot2::labs(color = group_name,fill=group_name )

    return(power_plot)

  }


  if(dim(power_data)[2] == 4){

    x_name<-InteractionPoweR::name_key[match(colnames(power_data)[1],InteractionPoweR::name_key[,1]  ),2]
    y_name<-InteractionPoweR::name_key[match(colnames(power_data)[2],InteractionPoweR::name_key[,1]  ),2]
    group_name<-InteractionPoweR::name_key[match(colnames(power_data)[3],InteractionPoweR::name_key[,1] ),2]
    facets_name<-InteractionPoweR::name_key[match(colnames(power_data)[4],InteractionPoweR::name_key[,1] ),2]

    power_data[,4]<-paste(facets_name,"=", power_data[,4])

    power_plot<-ggplot2::ggplot(data = power_data,ggplot2::aes(x = power_data[,1],y = power_data[,2],
                                             color=as.factor(power_data[,3]),
                                             fill=as.factor(power_data[,3]) ))+
      ggplot2::scale_color_viridis_d(option = c("C"),end = .95)+
      ggplot2::scale_fill_viridis_d(option = c("C"),end = .95)+
      ggplot2::geom_hline(yintercept = power_target,color="black")+
      ggplot2::geom_hline(yintercept = 0.05,color="darkgrey")+
      ggplot2::geom_hline(yintercept = 1,color="darkgrey")+
      ggplot2::geom_line(linetype="solid",size=.5,show.legend = F)+
      ggplot2::geom_smooth(formula = y ~ x,method = 'loess',se=T,alpha=0,size=1)+
      ggplot2::xlab(label = x_name)+
      ggplot2::ylab(label = y_name)+
      ggplot2::geom_point(shape=21,color="black",show.legend = F)+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.key=ggplot2::element_rect(size=1,color="black"))+
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1)))+
      ggplot2::labs(color = group_name,fill=group_name )+
      ggplot2::facet_wrap(facets = colnames(power_data)[4], scales = "free_y", strip.position = "top") +
      ggplot2::theme(strip.background = ggplot2::element_blank(), strip.placement = "outside")

    return(power_plot)

  }

}
