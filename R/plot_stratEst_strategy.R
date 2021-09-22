#' Plot Method for stratEst.strategy
#' @param x An object of class \code{stratEst.strategy}.
#' @param y Argument two of the generic function.
#' @param ... Further arguments passed to or from other methods.
#' @param title String. The title of the plot.
#' @param show.legend Logical. Hide plot legend if FALSE. Default is TRUE.
#' @param show.title Logical. Hide plot title if FALSE. Default is TRUE.
#' @param node.fontsize Font-size of the plot labels.
#' @param main.fontsize Font-size of the plot title.
#' @param arrow.fontsize Font-size of the arrow labels.
#' @param legend.fontsize Font-size of the legend.
#' @param legend.width Width of the legend items.
#' @param node.width Width of the nodes.
#' @param arrowsize Size of the arrowhead.
#' @param node.penwidth Width of the nodes.
#' @param arrow.penwidth Width of the nodes.
#' @param fillcolor Vector of hex-color codes of the choices.
#' @param ranksep Separation of nodes with the same rank.
#' @param file String. A valid path followed by a file name. Should end with .svg. Default is NA and no file is written.
#' @return No return value, called to create a plot.
#' @export

plot.stratEst.strategy <- function( x, y, ... , title = NULL, show.legend = TRUE, show.title = TRUE,  node.fontsize = 25, main.fontsize = 40, arrow.fontsize = 20, legend.fontsize = 20, legend.width = 1, node.width = 1, arrowsize = 1, node.penwidth = 1, arrow.penwidth = 1, fillcolor = NULL, ranksep = 0, file = NA ){

  if( is.null(fillcolor) ){
    def.palette <- c("#B3CDE3","#DECBE4","#CCEBC5","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#FBB4AE")
  }else{
    def.palette <- fillcolor
  }
  # plot options
  legend_start_options <-     paste("fontsize =",legend.fontsize,", style=doublecircle, peripheries = 2, fontname='sans-serif', color = gray60, fontcolor = black")
  legend_choice_options <-    paste(",fontsize =",legend.fontsize,", width=",legend.width,", shape=box, style = filled, fontname='sans-serif', color = white")
  start_state_options <-      paste(",fontsize =",node.fontsize,", width=",node.width,", penwidth=",node.penwidth,", shape = doublecircle, fontname='sans-serif',style='wedged', color = gray60")
  pure_start_state_options <- paste(",fontsize =",node.fontsize,", width=",node.width,", penwidth=",node.penwidth,", shape = doublecircle, fontname='sans-serif',style='filled', color = gray60")
  state_options <-            paste(",fontsize =",node.fontsize,", width=",node.width,", shape = doublecircle, fontname='sans-serif', style='wedged', color = white")
  pure_state_options <-       paste(",fontsize =",node.fontsize,", width=",node.width,", shape = doublecircle, fontname='sans-serif', style='filled', color = white")
  arrow_options <-            paste("fontsize =",arrow.fontsize,", arrowsize =",arrowsize,", penwidth=",arrow.penwidth,", fontname='sans-serif', fontsize = 18")

  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package \"DiagrammeR\" is needed for this function to work. Please install it.",
         call. = FALSE)
  }

    # retrieve information from strategy
    strategy <- x
    strategy_name <- deparse(substitute(x))
    plot_title <- ifelse( is.null(title), strategy_name, title )
    #strategy_name <- gsub(".*\\$","",strategy_name)
    states <- rownames(strategy)
    num.states <- length(states)
    choice_indices = grepl("prob.", colnames(strategy), fixed = TRUE)
    choices = gsub('prob.', '', colnames(strategy)[choice_indices])
    input_indices <- grepl("tr(", colnames(strategy), fixed = TRUE)
    inputs = gsub('tr', '', colnames(strategy)[input_indices])
    inputs = gsub('[[:punct:]]', '', inputs)
    num.choices = length(choices)
    num.inputs = length(inputs)
    probs <- strategy[,choice_indices]
    probs_NA <- apply(is.na(probs),1,sum)> 0
    if( "tremble" %in% colnames(strategy) ){
      for( i in 1:nrow(probs)){
        if( is.na(strategy$tremble[i]) == FALSE ){
          probs[i,] = probs[i,]*(1-strategy$tremble[i]) + ((1-probs[i,])*strategy$tremble[i])/(num.choices-1)
        }
      }
    }
    pure_row <- apply(probs == 1,1,sum) == 1
    pure_row[is.na(pure_row)] = FALSE
    pure_choice <- rep(NA,length(pure_row))
    for( i in 1:length(pure_row)){
      if( pure_row[i] ){
        pure_choice[i] = c(1:num.choices)[probs[i,] == 1]
      }else{
        pure_choice[i] = 0
      }
    }
    if( sum(input_indices) > 0 ){
      transitions <- strategy[,input_indices]
    }else{
      transitions <- NULL
    }

    string <- paste(
      "digraph finite_state_machine {",
      "rankdir=LR;",
      paste("graph [ranksep='",ranksep,"'];"),
      "forcelabels=true;",
      "subgraph cluster_B{",
      "penwidth = 0;",
      "labelloc='t';",
      paste("fontsize = ",main.fontsize,";",sep=""),
      "fontname='sans-serif';",
      ifelse( show.title , paste("label = '",plot_title,"';",sep=""), ""),
      ifelse( show.legend , paste(
      "placeholder[label=' ', color = white];",
      "start[",
      legend_start_options,"];",
      paste(sapply(c(1:num.choices), function(x) paste(choices[x],"[fillcolor='",def.palette[x],"'",legend_choice_options,"];",sep="")),collapse=""),
      paste("start","->",choices[1],"[style=invis];",sep=""),
      paste(sapply(c(1:(num.choices-1)), function(x) paste(choices[x],"->",choices[x+1],"[style=invis];",collapse="")),collapse=""),
      sep=""),""),
      "subgraph cluster_B{",
      "label = ' ';",
      "penwidth = 0;",
      "fixedsize=true;",
      ifelse( pure_row[1] | probs_NA[1] ,
              paste(  states[1],"[fillcolor= '",
                      ifelse(probs_NA[1],"#999999",def.palette[pure_choice[1]]),
                      "'",
                      pure_start_state_options,"];", collapse = '')
              ,
              paste(  states[1],"[fillcolor= '",
                      paste(sapply(c(1:num.choices), function(x) paste(def.palette[x],";",as.character(probs[1,x]),ifelse(x<num.choices,":",""),sep = "")),collapse=""),
                      "'",
                      start_state_options,"];", collapse = '')
      ),
      ifelse( num.states > 1 ,
              paste(
                sapply(c(2:num.states), function(s) ifelse( pure_row[s]  | probs_NA[s],
                                                            paste(  states[s],"[fillcolor= '",
                                                                    ifelse(probs_NA[s],"#999999",def.palette[pure_choice[s]]),
                                                                    "'",
                                                                    pure_state_options,"];", collapse = '')
                                                            ,
                                                            paste(  states[s],"[fillcolor= '",
                                                                    paste(sapply(c(1:num.choices), function(x) paste(def.palette[x],";",as.character(probs[s,x]),ifelse(x<num.choices,":",""),sep = "")),collapse=""),
                                                                    "'",
                                                                    state_options,"];", collapse = '')
                                                    )


                )
              , collapse = ''),
              ""),
      ifelse( is.null(transitions) == FALSE,
              paste( sapply(c(1:num.states),  function(s)  paste(
                paste( sapply(c(1:num.states), function(s2) ifelse( sum(transitions[s,] == s2) > 0, paste(states[s],"->",states[s2],"[label='",paste(paste(inputs[transitions[s,] == s2],sep=" "),collapse=", "),"'",arrow_options,"];",sep=""), "" )), collapse = "")
                ,collapse="")), collapse=""),""),
      "}",
      "}",
      "}"
    )

    if( is.na(file) == FALSE ){
        cat(DiagrammeRsvg::export_svg(DiagrammeR::grViz(string)),file=file)
    }

    DiagrammeR::grViz(string)
}

