### Functions for AIMS analysis

## Data is published openly and available on Mendeley Data

#--------------------------------------
## Install packages if missing and load
#--------------------------------------

install_load_packages <- function(packages_to_load){
  for (pkg in packages_to_load) {
    if (!require(pkg, character.only = TRUE)){
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }  
    else {
      library(pkg, character.only = TRUE)
    }
  } 
}

#------------------------------------
## Frequencies and proportions table
#------------------------------------

prop_table <- function(var){
  tab1 <- as_tibble(table(var))
  tab2 <- as_tibble(100*round(prop.table(table(var)), 3))
  # Combine tab1 and tab2, suppressing new name messages
  tab <- suppressMessages(bind_cols(tab1, tab2[2]))
  colnames(tab) <- c("Category", "N", "Perc")
  return(tab)
}

#--------------------------------------------
## Export multiple objects to a single sheet
#--------------------------------------------

export_objects <- function(wb, sheet_title, object_list){
  addWorksheet(wb, sheet_title)
  
  curr_row <- 1
  for(i in seq_along(object_list)) {
    writeData(wb, sheet_title,
              names(object_list)[i], startCol = 1, startRow = curr_row)
    writeData(wb, sheet_title,
              object_list[[i]], startCol = 1, startRow = curr_row+1)
    curr_row <- curr_row + nrow(object_list[[i]]) + 3
  }
}

#--------------------------------------------
## Correlation matrix with significance stars
#--------------------------------------------

# x is a matrix containing the data
# method : correlation method. "pearson"" or "spearman"" is supported
# requires Hmisc package
# code from http://www.sthda.com/english/wiki/elegant-correlation-table-using-xtable-r-package

cor_stars <-function(x, method=c("pearson", "spearman")){
  
  #Compute correlation matrix
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coefficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))
  
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle of correlation matrix
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)
  
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew)
} 

#----------------------------------------
## Weighted pairwise t-tests (one sample)
#----------------------------------------

# Requires "weights" package

weighted_pw_ttests <- function(dataset, variables, weight){
  
  weighted_results <- NULL
  
  variables2 <- variables
  
  for (i in variables){
    for (j in variables2){
      if (i != j){
        ttest_result <- (wtd.t.test(dataset[[i]] - dataset[[j]], weight = weight))
        
        tmp_res <- NULL
        tmp_res$var1 <- i
        tmp_res$var2 <- j
        tmp_res$difference <- ttest_result$additional[1]
        tmp_res$std_error <- ttest_result$additional[4]
        tmp_res$t_stat <- ttest_result$coefficients[1]
        tmp_res$df <- ttest_result$coefficients[2]
        tmp_res$p_value <- ttest_result$coefficients[3]
        
        weighted_results <- bind_rows(weighted_results, tmp_res)
      }
    }
    
    if(!is.null(variables2)){
      variables2 <- variables2[-1]
    }
    
  }
  
  return(weighted_results)
  
}

#-------------------------------------------------
# Weighted pairwise t-tests (independent samples)
#-------------------------------------------------

### Survey weighted pairwise t-tests (independent samples)

# Requires "survey" package

svy_pw_ttests_ind <- function(dataset, y, x, survey_design){
  
  weighted_results <- NULL
  
  x_vars <- unique(dataset[[x]])
  x_vars2 <- x_vars
  
  for (i in x_vars){
    for (j in x_vars2){
      if (i != j){
        
        subset_design <- subset(survey_design, dataset[[x]] %in% c(i, j))
        
        ttest_result <- svyttest(as.formula(str_c(y, "~", x)), design = subset_design)
        tmp_res <- NULL
        tmp_res$var1 <- if_else(i < j, j, i)
        tmp_res$var2 <- if_else(i < j, i, j)
        tmp_res$difference <- ttest_result$estimate
        tmp_res$t_stat <- ttest_result$statistic
        tmp_res$df <- ttest_result$parameter
        tmp_res$ci_lower <- ttest_result$conf.int[1]
        tmp_res$ci_upper <- ttest_result$conf.int[2]
        tmp_res$p_value <- ttest_result$p.value
        
        weighted_results <- bind_rows(weighted_results, tmp_res)
      }
    }
    
    if(!is.null(x_vars2)){
      x_vars2 <- x_vars2[-1]
    }
    
  }
  
  names(weighted_results) <- c("Variable 1", "Variable 2", "Difference",
                               "t-statistic", "Degrees of freedom",
                               "CI lower", "CI upper", "p-value")
  
  return(weighted_results)
  
}


#---------------------------------
## Layered, interactive histogram
#---------------------------------

# Requires packages "survey" and "matrixStats",
# "tidyverse", "plotly" and "htmlwidgets"

plotly_histogram <- function(dataset, top_var, svy_design, var_names,
                             var_n, var_labs, var_titles,
                             bins = NULL, x_scale = NULL, x_labels = NULL,
                             tick_text = NULL, tick_vals = NULL,
                             filename = NULL){
  
  # Initialize plot
  layered_plot <- plot_ly(dataset)
  
  # Data for top layer
  plot_data <- as_tibble(svytable(formula = top_var, design = svy_design))
  colnames(plot_data) <- c("plot_var", "n")
  plot_data$perc <- plot_data$n/sum(plot_data$n)*100
  plot_data$plot_var <- as.numeric(plot_data$plot_var)
  
  # Create bins of data if required
  if(!is.null(bins)){
    plot_data <- plot_data %>%
      group_by(group = cut(plot_var, breaks = bins[[1]])) %>%
      summarise(n = sum(n),
                perc = sum(perc))
    plot_data <- rename(plot_data, c("plot_var" = "group"))
    plot_data$plot_var <- x_scale[[1]]
  }
  
  # Histogram
  layered_plot <- add_bars(layered_plot, 
                           x = plot_data$plot_var,
                           y = plot_data$perc,
                           visible = T,
                           marker = list(color = "#004b8c",
                                         line = list(color = "darkgray",
                                                     width = 1)),
                           hovertemplate = paste('%{y:.2f}%<extra></extra>'),
                           name = "Percent")
  
  # Mean line
  if(top_var != "~md$F1"){
    layered_plot <- add_segments(layered_plot, 
                                 x = svymean(top_var, design = svy_design, na.rm = T), 
                                 xend = svymean(top_var, design = svy_design, na.rm = T),
                                 y = 0, 
                                 yend = max(plot_data$perc),
                                 text = sqrt(svyvar(top_var, design = svy_design, na.rm = T))[1],
                                 #visible = T, 
                                 line = list(color="red",
                                             width = 3),
                                 hovertemplate = paste('Mean: %{x:.2f}',
                                                       '<br>Standard deviation: %{text:.2f}<extra></extra>'),
                                 name = "Mean (SD)")
  }
  
  # Empty plot for F1 as no mean line (categorical variable)
  if(top_var == "~md$F1"){
    layered_plot <- add_segments(layered_plot,
                                 x = 0,
                                 xend = 0,
                                 y= 0,
                                 yend = 0,
                                 visible = F)
  }
  
  # Add remaining plot layers
  for (i in var_names[2:var_n]) {
    var_nr <- (which(var_names == i))
    plot_data <- as_tibble(svytable(bquote(~.(as.name(i))), design = svy_design))
    colnames(plot_data) <- c(i, "n")
    plot_data$perc <- plot_data$n/sum(plot_data$n)*100
    plot_data[[i]] <- as.numeric(plot_data[[i]])
    
    # Create bins of data if required
    if(!is.null(bins)){
      plot_data <- plot_data %>%
        group_by(group = cut(plot_data[[i]], breaks = bins[[var_nr]])) %>%
        summarise(n = sum(n),
                  perc = sum(perc))
      plot_data <- rename(plot_data, c(i = "group"))
      plot_data[[i]] <- x_scale[[var_nr]]
    }
    
    # Histogram
    layered_plot <- add_bars(layered_plot, 
                             x = plot_data[[i]],
                             y = plot_data$perc,
                             visible = F,
                             marker = list(color = "#004b8c",
                                           line = list(color = "darkgray",
                                                       width = 1)),
                             hovertemplate = paste('%{y:.2f}%<extra></extra>'),
                             name = "Percent")
    
    # Mean line if not F2 and F3, for which we plot medians
    if(i != "F2_cat" & i != "F3_cat"){
      layered_plot <- add_segments(layered_plot, 
                                   x = svymean(md[[i]], design = svy_design, na.rm = T), 
                                   xend = svymean(md[[i]], design = svy_design, na.rm = T),
                                   y = 0, 
                                   yend = max(plot_data$perc),
                                   text = sqrt(svyvar(md[[i]], design = svy_design, na.rm = T))[1],
                                   visible = F, 
                                   line = list(color="red",
                                               width = 3),
                                   hovertemplate = paste('Mean: %{x:.2f}',
                                                         '<br>Standard deviation: %{text:.2f}<extra></extra>'),
                                   name = "Mean (SD)")
    }
    # Median line for F2 and F3
    if(i == "F2_cat" | i == "F3_cat"){
      i <- str_remove(i, "_cat")
      layered_plot <- add_segments(layered_plot, 
                                   x = weightedMedian(md[[i]], md$weight, na.rm = T), 
                                   xend = weightedMedian(md[[i]], md$weight, na.rm = T),
                                   y = 0, 
                                   yend = max(plot_data$perc),
                                   text = weightedMad(md[[i]], md$weight, na.rm=T),
                                   visible = F, 
                                   line = list(color="red",
                                               width = 3),
                                   hovertemplate = paste('Median: %{x:.2f}',
                                                         '<br>Median absolute deviation: %{text:.2f}<extra></extra>'),
                                   name = "Median (MAD)")
    }
  }
  
  # Setup button for showing traces before plot formatting
  
  # Initial button shows first two traces
  init_button <- rep(c(TRUE, FALSE), times = c(2, (var_n-1)*2))
  
  # Create a vector of Fs for remaining remaining traces
  init_button_rep <- rep(FALSE, var_n*2)
  
  # Replicate previous vector for each remaining trace
  button_config <- replicate(var_n-1, rep(init_button_rep), simplify = FALSE)
  
  # Create a list of all the buttons
  button_config <- c(list(init_button),button_config)
  
  # For the buttons after the first, recode to "TRUE" to show desired traces
  for (i in 2:var_n){
    button_config[[i]][i*2-1] <- TRUE
    button_config[[i]][i*2] <-  TRUE
  }
  
  # For the F1 and F11, hide the mean lines
  if (top_var == "~md$F1"){
    for (i in 1:var_n){
      button_config[[i]][2] <- FALSE
      button_config[[i]][4] <- FALSE
    }
  }
  
  # Plot formatting
  layered_plot <- layout(layered_plot,
                         xaxis = list(
                           ticktext = tick_text[[1]], 
                           tickvals = tick_vals[[1]],
                           showline = F),
                         title = list(text = paste(
                           strwrap(var_titles[1], width = 100),
                           collapse = "\n"),
                           font = list(size = 13), x = 0.43),
                         yaxis = list(title = "Percent", hoverformat = ".1f"),
                         xaxis = list(title = "Response", zeroline = FALSE),
                         bargap = 0,
                         showlegend = T,
                         updatemenus = list(
                           list(type = "dropdown", y = 1.4, x = 0.3,
                                buttons = lapply(1:var_n, function(i){
                                  list(method = "update",
                                       args = list(
                                         list(visible = button_config[[i]]),
                                         list(title = list(text = paste(
                                           strwrap(var_titles[i],width = 100),
                                           collapse = "\n"), 
                                           font = list(size = 13), x = 0.43),
                                           xaxis = list(ticktext = tick_text[[i]], tickvals = tick_vals[[i]],
                                                        showline = F),
                                           yaxis = list(title = "Percent", hoverformat = ".1f"),
                                           bargap = 0)),
                                       label = var_labs[i])
                                }
                                ))
                         )
  )
  
  # Save widget (suppress warnings that file already exists)
  suppressWarnings(htmlwidgets::saveWidget(as_widget(layered_plot), filename))
  
  return(layered_plot)
}

#------------------------
## Add significance stars
#------------------------

# Define a function to assign stars
get_sig_stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*","")
                )
         )
}
