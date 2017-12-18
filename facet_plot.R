##' plot tree associated data in an additional panel
##'
##'
##' @title facet_plot
##' @param p tree view
##' @param panel panel name for plot of input data
##' @param data data to plot by 'geom', first column should be matched with tip label of tree
##' @param geom geom function to plot the data
##' @param mapping aes mapping for 'geom'
##' @param grid_space space parameter for facet_grid
##' @param ... additional parameters for 'geom'
##' @return ggplot object
##' @export
##' @author Guangchuang Yu
facet_plot <- function(p, panel, data, geom, mapping=NULL, grid_space=NULL, ...) {
	p <- add_panel(p, panel, grid_space)
    df <- p %+>% data
    p + geom(data=df, mapping=mapping, ...)
}

##' @importFrom ggplot2 facet_grid
add_panel <- function(p, panel, grid_space) {
	df <- p$data
    if (is.null(df$.panel)) {
       df$.panel <- factor("Tree")
    }
    levels(df$.panel) %<>% c(., panel)
    p$data <- df
    if (is.null(grid_space)){
    	p + facet_grid(.~.panel, scales="free_x")
	}else{
		p + facet_grid(.~.panel, scales="free_x", space=grid_space)
	}
}

