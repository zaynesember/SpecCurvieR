# FUNCTIONS TO PERFORM SPECIFICATION CURVE ANALYSIS

# Load packages
require(tidyverse)
require(combinat)
require(lfe)
require(stringr)
require(snow)

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
formula_builder <- function(y, x, controls, fixedEffects=NA){
  
  # Get all combinations of controls
  powerset <- unlist(lapply(1:length(controls),
                            combinat::combn, 
                            x = controls,
                            simplify = FALSE),
                     recursive=F)
  
  # Remove duplicate controls that are already in the interaction
  powerset <- unique(sapply(X=powerset, FUN=duplicate_remover, x=x))
  
  #Build right hand side of the formulae
  if(is.na(fixedEffects)){
    RHS <- unique(sapply(powerset, paste_factory, x))
  }
  else{
    RHS <- paste(unique(sapply(powerset, paste_factory, x)), fixedEffects, sep=" | ")
  }
  # Build formulae
  formulae <- sapply(paste(y, RHS, sep=" ~ "), formula)
  
  return(formulae)
}

# Runs every possible combination of regression models using lm() or felm()
# Arguments:
#   y = string name of dependent variable column
#   x = string name of independent variable column
#   controls = vector of strings containing desired controls and interactions
#   data = dataframe object containing y, x, controls, and fixed effects
#          variables
#   fixedEffects = string name of variable to use for fixed effects if desired
# Returns: dataframe object containing the coefficient, standard error, t-value,
#          p-value, list of terms in the regression, and significance level
sca <- function(y, x, controls, data, fixedEffects = NULL, 
                parallel=FALSE, workers=2){
  # With parallel computing
  if(parallel){
    
    # No fixed effects specified
    if(is.null(fixedEffects)){
      
      cl <- makeSOCKcluster(rep("localhost", workers))
      
      # Build the formulae
      formulae <- formula_builder(y=y, x=x, controls=controls)
      
      # Estimate the models with lm()
      clusterExport(cl, "formulae", envir=environment()) 
      clusterExport(cl, "data", envir=environment()) 
      models <- parSapply(cl, formulae, function(x2) summary(lm(x2, data=data)))
      }
    # Fixed effects specified
    else{
      formulae <- formula_builder(y=y, x=x, controls=controls, fixedEffects=fixedEffects)
      
      clusterExport(cl, "formulae", envir=environment()) 
      clusterExport(cl, "data", envir=environment())
      
      # Estimate the models with felm()
      models <- parSapply(cl, formulae, function(x2) summary(felm(x2, data=data)))
      }
    }
  # Without parallel computing
  else{
    # No fixed effects specified
    if(is.null(fixedEffects)){
      models <- sapply(formulae, function(x2) summary(lm(x2, data=data)))
    }
    # Fixed effects specified
    else{
      # Build the formulae
      formulae <- formula_builder(y=y, x=x, controls=controls, fixedEffects=fixedEffects)
      
      # Estimate the models with felm()
      models <- sapply(X=formulae, function(x2) summary(felm(x2, data=data)))      
    }
  }
  
  # Extract results for IV
  xVals <- apply(X=models, MARGIN=2, FUN=function(x2) x2$coefficients[x,])
  
  # Get each value of interest across models
  coef <- xVals[1,]
  se <- xVals[2,]
  t <- xVals[3,]
  p <- xVals[4,]
  terms <- names(apply(X=models, MARGIN=2, FUN=function(x2) x2$terms[[1]]))
  
  # Put into a dataframe
  retVal <- data.frame(terms, coef, se, t, p) %>% 
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

  return(retVal)
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
  
  return(list(df, setNames(as.character(df_labels$control), df_labels$controlID)))
}

# TODO: Documentation and testing
plotCurve <- function(sca_data, title="", 
                         y_lab="Independent variable coefficient"){
  
  pointSize <- -.25*(ncol(sca_data)-7)+(13/4)
  
  sc <- ggplot(data=sca_data, aes(y=coef, x=index, fill=factor(sig.level))) +
    geom_hline(yintercept = 0, color="red", linetype="dashed", linewidth=.75) +
    geom_ribbon(aes(ymin=coef-se, ymax=coef+se), alpha=.5) +
    geom_point(size=pointSize) +
    labs(title=title, x="", y=y_lab) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      legend.position="top",
      legend.title=element_blank(),
      legend.key.size = unit(.4, 'cm')
    )
  
  return(sc)
}

# TODO: Documentation and testing
plotVars <- function(sca_data){
  scp_data <- scp(sca_data)

  markSize <- 10/length(scp_data[[2]])
  
  sc <- ggplot(data=scp_data[[1]],
                aes(x=index,y=factor(controlID))
  ) +
    geom_point(shape="|", size=markSize) +
    labs(y="", x="") +
    scale_y_discrete(labels=scp_data[[2]], expand=c(.25,.25)) +
    theme_void() +
    theme(
      axis.text.y = element_text(size=6, hjust=0),
      axis.text.x = element_blank()
    )
  return(sc)
}

# TODO: Documentation and testing
plotSCV <- function(y, x, controls, data, fixedEffects = NULL, combine=T){
  df_sca <- sca(y=y, x=x, controls=controls, data=data, fixedEffects=fixedEffects)
  df_scp <- scp(df_sca)

  sc1 <- plotCurve(df_sca)
  sc2 <- plotVars(df_sca)
  
  if(combine){
    grid::grid.newpage()
    return(grid::grid.draw(rbind(ggplotGrob(sc1), ggplotGrob(sc2))))
  }
  else return(list(sc1, sc2))
}

# TODO: Add estimated time/updates in console
# TODO: Fix fixedEffects in parallel
