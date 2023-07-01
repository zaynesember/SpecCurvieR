# FUNCTIONS TO PERFORM SPECIFICATION CURVE ANALYSIS

# Load packages
require(tidyverse)
require(combinat)
require(lfe)
require(stringr)
require(pbapply)
require(parallel)

# Pastes together x and controls to be fed into formula()
# Arguments:
#   controls = vector of strings containing control variable column names
#   x = string name of independent variable of interest's column name
# Returns: String of the independent variable and controls pasted together, 
# ready to be put in a formula object
paste_factory <- function(controls, x){
  if(T %in% str_detect(controls, x)){
    return(paste(controls, collapse=" + "))
  }
  else return(paste(x, paste(controls, collapse=" + "), sep=" + "))
}

# Removes controls already in interaction
# i.e. lm(y ~ x + c1 + x*c1) does not need the solo c1 term
# Arguments:
#   controls = vector of control variable column names as strings
#   x = string name of independent variable of interest's column name
# Returns: A vector of strings containing control variables and interaction 
# names
duplicate_remover <- function(controls, x){
  # Check for interactions
  if(T %in% str_detect(controls, "\\*")){
    # Find interaction terms
    indices <- which(T==str_detect(controls, "\\*"))
    # Find controls that are in interaction terms
    extraTerms <- str_replace(str_replace(controls[indices], 
                                          pattern=x, 
                                          replacement=""),
                              pattern="\\*",
                              replacement="")
    # Remove controls that are already present in interaction
    return(controls[!controls %in% extraTerms])
  }
  else return(controls)
}

# Builds all regression formula combinations possible with given controls
# Arguments:
#   y = string name of dependent variable column
#   x = string name of independent variable column
#   controls = vector of strings containing desired controls and interactions
#   fixedEffects = string name of variable to use for fixed effects if desired
# Returns: vector of formula objects
formula_builder <- function(y, x, controls, fixedEffects=NA, cluster=NA){
  
  # Get all combinations of controls
  powerset <- unlist(lapply(1:length(controls),
                            combinat::combn, 
                            x = controls,
                            simplify = FALSE),
                     recursive=F)
  
  # Remove duplicate controls that are already in the interaction
  powerset <- unique(sapply(X=powerset, FUN=duplicate_remover, x=x))
  
  # Build right hand side of the formulae
  if(is.na(fixedEffects)){
    if(is.na(cluster)) RHS <- unique(sapply(powerset, paste_factory, x))
    else RHS <- paste(unique(sapply(powerset, paste_factory, x)), "", 
                      cluster, sep=" | ")
  }
  else{
    if(is.na(cluster)) RHS <- paste(unique(sapply(powerset, paste_factory, x)), 
                                   fixedEffects, sep=" | ")
    else RHS <- paste(unique(sapply(powerset, paste_factory, x)), fixedEffects, 
                      cluster, sep=" | ")
  }
  # Build formulae
  formulae <- sapply(paste(y, RHS, sep=" ~ "), formula)
  
  return(formulae)
}

# TODO: UPDATE DOCUMENTATION
# Runs every possible combination of regression models using lm() or felm()
# Arguments:
#   y = string name of dependent variable column
#   x = string name of independent variable column
#   controls = vector of strings containing desired controls and interactions
#   data = dataframe object containing y, x, controls, and fixed effects
#          variables
#   fixedEffects = string adding fixed effects variables if desired
# Returns: dataframe object containing the coefficient, standard error, t-value,
#          p-value, list of terms in the regression, and significance level
sca <- function(y, x, controls, data, fixedEffects=NULL,
                cluster=NULL, returnModels=F, plot=F, colorControls=T,
                combinePlots=T, plotFits=F, progressBar=T, 
                parallel=FALSE, workers=2){
  # With parallel computing
  if(parallel){
    
    cl <- makePSOCKcluster(rep("localhost", workers))
    
    # Load needed package into each cluster
    clusterEvalQ(cl, library(lfe))
    
    # No fixed effects specified
    if(is.null(fixedEffects)){
      # Build the formulae
      formulae <- formula_builder(y=y, x=x, controls=controls)
      
      # Estimate the models with lm()
      clusterExport(cl, "formulae", envir=environment()) 
      clusterExport(cl, "data", envir=environment()) 
      
      # Show progress bar if desired
      if(progressBar){
        print.noquote(paste("Estimating", length(formulae), "models in parallel with", 
                    workers, "workers"))
        system.time(models <- pbsapply(formulae, 
                                       function(x2) summary(lm(x2, data=data)),
                                       cl=cl))
      }
      else{
        models <- parSapply(cl, formulae, function(x2) summary(lm(x2, data=data)))
      }
    }
    # Fixed effects specified
    else{
      # Build the formulae
      formulae <- formula_builder(y=y, x=x, controls=controls, 
                                       fixedEffects=fixedEffects)
      
      clusterExport(cl, "formulae", envir=environment()) 
      clusterExport(cl, "data", envir=environment())
      
      if(progressBar){
        print.noquote(paste("Estimating", length(formulae), "models in parallel with", 
                    workers, "workers"))
        system.time(models <- pbsapply(formulae, 
                                       function(x2) summary(felm(x2, data=data)),
                                       cl=cl))
      }
      else{
        models <- parSapply(cl, formulae, function(x2) summary(felm(x2, data=data)))
      }
    }
  }
  # Without parallel computing
  else{
    # No fixed effects specified
    if(is.null(fixedEffects)){
      # Build the formulae
      formulae <- formula_builder(y=y, x=x, controls=controls)
      
      if(progressBar){
        print.noquote(paste("Estimating", length(formulae), "models"))
        system.time(models <- pbsapply(formulae, 
                                       function(x2) summary(lm(x2, data=data))))
      }
      else{
        models <- sapply(formulae, function(x2) summary(lm(x2, data=data)))
      }
    }
    # Fixed effects specified
    else{
      # Build the formulae
      formulae <- formula_builder(y=y, x=x, controls=controls,
                                       fixedEffects=fixedEffects)
      
      if(progressBar){
        print.noquote(paste("Estimating", length(formulae), "models"))
        system.time(models <- pbsapply(X=formulae, 
                                       function(x2) summary(felm(x2,data=data))))
      }
      else{
        models <- sapply(X=formulae, function(x2) summary(felm(x2, data=data)))
      }
    }
  }
  
  # Garbage collection for parallel connections
  if(parallel) stopCluster(cl=cl)
  
  # If the user just wants a list of the model summary objects return that and 
  # we're done
  if(returnModels) return(models)
  
  # Extract results for IV
  vals <- apply(X=models, MARGIN=2, 
                 FUN=function(x2) list(x2$coefficients[x,],
                                       sqrt(mean(x2$residuals^2)),
                                       x2$adj.r.squared))
  # Get each value of interest across models
  coef <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 1))
  se <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 2))
  t <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 3))
  p <- unlist(lapply(lapply(vals, `[[`, 1), `[[`, 4))
  terms <- names(apply(X=models, MARGIN=2, FUN=function(x2) x2$terms[[1]]))
  RMSE <- unlist(lapply(lapply(vals, `[[`, 2), `[[`, 1))
  adjR <- unlist(lapply(lapply(vals, `[[`, 3), `[[`, 1))
  # Put into a dataframe
  retVal <- data.frame(terms, coef, se, t, p, RMSE, adjR) %>% 
    mutate(
      sig.level=case_when(
        p < .005 ~ "p < .005",
        p < .05 ~ "p < .05",
        p < .1 ~ "p < .1",
        p >= .1 ~ "p >= .1",
        T ~ NA_character_
      )) %>% 
    arrange(coef) %>% 
    mutate(index=row_number())
  
  # Build dummy columns for terms present in each model for visualization
  temp <- data.frame(matrix(ncol = length(controls), nrow = nrow(retVal)))
  
  colnames(temp) <- controls
  
  retVal <- cbind(retVal, temp)
  
  for(c in controls){
    retVal[c] <- ifelse(str_detect(retVal$terms, fixed(c)), 1, 0)
  }
  
  # Generate plots if desired
  if(plot){
    grid::grid.newpage()
    sc1 <- plotCurve(retVal)
    sc2 <- plotVars(retVal, colorControls)
    
    if(plotFits){
      sc3 <- plotRMSE(retVal)
      sc4 <- plotR2Adj(retVal)
    }
    
    if(combinePlots){
      
      if(plotFits){
        plots <- grid::grid.draw(rbind(ggplotGrob(sc1), 
                                       ggplotGrob(sc3),
                                       ggplotGrob(sc4),
                                       ggplotGrob(sc2)))
      }
      else plots <- grid::grid.draw(rbind(ggplotGrob(sc1), ggplotGrob(sc2)))
      
      
      return(plots)
    }
    else{
      if(plotFits) return(list(coefficient=sc1, controls=sc2, 
                               RMSE=sc3, R2Adj=sc4))
      else return(list(coefficient=sc1, controls=sc2))
    }
  }
  # Or just return model parameters
  else return(retVal)
}

# Takes in the output of sca() and returns a list with the dataframe and 
# labels to make a plot to visualize the controls included in each spec curve
# model
# Arguments:
#   spec_data = dataframe object with output from sca()
# Returns: list containing dataframe, controls, and control IDs
scp <- function(spec_data){
  df <- spec_data %>% 
    select(-terms, -coef, -se, -t, -p, -sig.level) %>% 
    pivot_longer(!index, names_to="control", values_to="value") %>% 
    filter(value==1) %>% 
    mutate(controlID = with(.,match(control, unique(control)))) %>% 
    select(-value)
  
  df_labels <- df %>% select(control, controlID) %>% unique()
  
  return(list(df, setNames(as.character(df_labels$control), 
                           df_labels$controlID)))
}

# TODO: Documentation and testing
plotCurve <- function(sca_data, title="", 
                         y_lab="Coefficient"){
  
  pointSize <- -.25*(ncol(sca_data)-7)+(13/4)
  
  sc <- ggplot(data=sca_data, aes(y=coef, x=index, 
                                  fill=factor(sig.level))) +
    geom_hline(yintercept = 0, color="red", linetype="dashed", linewidth=.75) +
    geom_ribbon(aes(ymin=coef-se, ymax=coef+se), alpha=.4, fill="red", color="darkred") +
    geom_point(aes(fill=as.factor(sig.level)),size=pointSize) +
    labs(title=title, x="", y=y_lab) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(vjust=-0.5),
      legend.position="none",
      legend.title=element_blank(),
      legend.key.size = unit(.4, 'cm'),
      plot.margin = unit(c(-15,1,-5,1), "points")
    )
  
  return(sc)
}

# TODO: Documentation and testing
plotVars <- function(sca_data, colorControls=F){
  scp_data <- scp(sca_data)

  markSize <- 10/length(scp_data[[2]])
  
  if(colorControls){
    sc <- ggplot(data=scp_data[[1]],
                  aes(x=index,y=factor(controlID), color=factor(controlID))
    ) +
      geom_point(shape="|", size=markSize) +
      labs(y="", x="") +
      scale_y_discrete(labels=scp_data[[2]], expand=c(.25,.25)) +
      theme_void() +
      theme(
        legend.position = "none",
        axis.text.y = element_text(size=6, hjust=0),
        axis.text.x = element_blank(),
        plot.margin = unit(c(-5,1,-5,1), "points")
      )
  }
  else{
      sc <- ggplot(data=scp_data[[1]],
                   aes(x=index,y=factor(controlID))
      ) +
        geom_point(shape="|", size=markSize) +
        labs(y="", x="") +
        scale_y_discrete(labels=scp_data[[2]], expand=c(.25,.25)) +
        theme_void() +
        theme(
          legend.position = "none",
          axis.text.y = element_text(size=6, hjust=0),
          axis.text.x = element_blank(),
          plot.margin = unit(c(-5,1,-5,1), "points")
        )
    }    
  return(sc)
}

# TODO: Documentation and testing
plotRMSE <- function(sca_data, title=""){
  
  pointSize <- -.25*(ncol(sca_data)-7)+(13/4)
  
  sc <- ggplot(data=sca_data, aes(y=RMSE, x=index)) +
    geom_point(size=pointSize) +
    labs(title=title, x="", y="RMSE") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      legend.title=element_blank(),
      legend.key.size = unit(.4, 'cm'),
      plot.margin = unit(c(-20,1,-5,1), "points")
    )
}

# TODO: Documentation and testing
plotR2Adj <- function(sca_data, title=""){
  
  pointSize <- -.25*(ncol(sca_data)-7)+(13/4)
  
  sc <- ggplot(data=sca_data, aes(y=adjR, x=index)) +
    geom_point(size=pointSize) +
    labs(title=title, x="", y=bquote('Adj. R'^2)) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      legend.title=element_blank(),
      legend.key.size = unit(.4, 'cm'),
      plot.margin = unit(c(-20,1,-5,1), "points")
    )
}

# TODO: Variance decomposition
# TODO: Add support for grouping models
# TODO: Swap lm() for glm() and ensure support for all glm() models
# TODO: Vary control tick colors
# TODO: Add facet-wrapped histograms of control coefficients
# TODO: Add other measures of model fit (AIC, BIC, Log. Likelihood)
# TODO: Add support for IV via felm
# TODO: Add support random effects
# TODO: Add support for clustered SEs
# TODO: Add support for robust SEs
# TODO: Add support for comparison between SE types
# TODO: Add plot of number of observations across models
#       Maybe add some kind of power analysis?
# TODO: Compare model estimation speed to specR