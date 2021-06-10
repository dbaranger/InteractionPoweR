
#' Simple slope plot
#'
#' Plots the simple slope min and max estimates from power_interaction().
#'
#' @param power_data Data frame of results from power_interaction(). Can accept the raw results if up to 2 parameters were varied during simulation. Any more and data should be filtered first.
#' @param x Optional, the x-axis of the plot. Default is the first variable after 'pwr'.
#' @param facets Optional, grouping variable for plot facets. Default is the second variable after 'pwr' if present.
#' @importFrom dplyr "%>%"

#' @return A ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' dataset <- generate_interaction(N = seq(100,300,by=10),
#' r.x1.y = 0,r.x2.y = .1,r.x1x2.y = -.2,r.x1.x2 = .3)
#' simple_slope_plot(dataset)
#' }
plot_simple_slope<-function(power_data,x=NULL,facets=NULL){


  y = "pwr"

  # Default behavior is to use columns in the order they are present in the power analysis output
  if(is.null(x) &&  is.null(facets)){
    power_data2<-power_data[,c(1: which(colnames(power_data) == "pwr"))]
    x = colnames(power_data2)[1]

    if(dim(power_data2)[2]>2){ facets = colnames(power_data2)[2] }
    #if(dim(power_data2)[2]>3){ facets = colnames(power_data2)[3] }
    if(dim(power_data2)[2]>3){
      print("Too many variable combinations in input data. Please select a subset of data to plot")
      stop()
    }
  }

  x_col<-which(colnames(power_data) == x)
  #group_col<-which(colnames(power_data) == group)
  facets_col<-which(colnames(power_data) == facets)

  col_keep<- c(unname(cbind(x_col,facets_col)))
  power_data2<-as.data.frame(power_data[,col_keep])

  #bounds<-power_data[,c((dim(power_data)[2]-3): (dim(power_data)[2]) )]

  col_names = c("min.lwr", "min.upr","max.lwr","max.upr" )

  bounds <- power_data[ , base::match(col_names,colnames(power_data))]
  power_data2<-cbind(power_data2,bounds)

  if(dim(power_data2)[2] == 5){

    x_name<-InteractionPoweR::name_key[match(colnames(power_data)[1],InteractionPoweR::name_key[,1]  ),2]


    power_data2<-power_data2 %>% tidyr::pivot_longer(
      cols = c("min.lwr","min.upr","max.lwr","max.upr"),
      names_to = "bound",
      #names_sep = "[.]",
      values_to = "slopes")
    power_data2$bound2<- matrix(unlist(strsplit(power_data2$bound,split  = "[.]")),nrow   = 2)[1,]
    #power_data2$bound2<-stringr::str_split(power_data2$bound,pattern = "[.]",simplify = T)[,1] # need to redo this
    power_data2$bound3<-paste(as.matrix(power_data2[,1]),power_data2$bound2,sep="_")

    slope_plot<-ggplot2::ggplot(data = power_data2,ggplot2::aes(x =as.matrix(power_data2[,1]),
                                                                y = .data$slopes,
                                                                fill= .data$bound2,
                                                                group=.data$bound3))+
      ggplot2::scale_fill_viridis_d()+
      ggplot2::geom_hline(yintercept = 0,color="darkgrey")+
      ggplot2::geom_line(linetype="solid",size=.5,color = "black")+
      ggplot2::geom_point(shape=21,color="black",show.legend = F,size=2)+

      ggplot2::ylab(label = "Simple slope extrema IQR")+
      ggplot2::xlab(label = x_name)+
      #scale_color_viridis_d(option = c("C"),end = .5)+
      ggplot2::theme_minimal()




    return(slope_plot)

  }


  if(dim(power_data2)[2] == 6){

    x_name<-InteractionPoweR::name_key[match(colnames(power_data)[1],InteractionPoweR::name_key[,1]  ),2]
    #y_name<-InteractionPoweR::name_key[match(colnames(power_data)[2],InteractionPoweR::name_key[,1]  ),2]
    facet_name<-InteractionPoweR::name_key[match(colnames(power_data)[2],InteractionPoweR::name_key[,1] ),2]


    power_data2<-power_data2 %>% tidyr::pivot_longer(
      cols = c("min.lwr","min.upr","max.lwr","max.upr"),
      names_to = "bound",
      #names_sep = "[.]",
      values_to = "slopes")

    power_data2$bound2<- matrix(unlist(strsplit(power_data2$bound,split  = "[.]")),nrow   = 2)[1,]
    power_data2$bound3<-paste(as.matrix(power_data2[,1]),power_data2$bound2,sep="_")


    slope_plot<-

      ggplot2::ggplot(data = power_data2,ggplot2::aes(x =as.matrix(power_data2[,1]),

                                    color=as.factor(as.matrix(power_data2[,2])),
                                    fill=as.factor(as.matrix(power_data2[,2])),

                                    y = .data$slopes,
                                    group=.data$bound3))+
      ggplot2::geom_hline(yintercept = 0,color="black")+
      ggplot2::scale_color_viridis_d(option = c("C"),end = .95)+
      ggplot2::scale_fill_viridis_d(option = c("C"),end = .95)+
      ggplot2::geom_line(linetype="solid",size=.5)+
      ggplot2::geom_point(shape=21,color="black",show.legend = T,size=2)+

      ggplot2::ylab(label = "Simple slope extrema IQR")+
      ggplot2::xlab(label = x_name)+
      #scale_color_viridis_d(option = c("C"),end = .5)+
      ggplot2::theme_minimal()+
      ggplot2::facet_wrap(facets = colnames(power_data)[2], scales = "free_y",
                 strip.position = "top",labeller = ggplot2::label_both ) +
      ggplot2::theme(strip.background = ggplot2::element_blank(), strip.placement = "outside")+
      ggplot2::labs(color = facet_name,fill=facet_name )

    return(slope_plot)



  }



}
