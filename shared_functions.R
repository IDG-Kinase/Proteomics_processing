LoB <- function(x)
{
  ##method1:Linett paper
  ##x: a vecotr of values
  x <- na.exclude(sort(x))
  Nb <- length(x)
  obs.idx <- Nb*0.95+0.5
  if(obs.idx%in%(1:Nb)) quan <- x[obs.idx]
  else if (obs.idx>1 & obs.idx<Nb)      
  {
    #quan <- 0.5*(x[floor(obs.idx)]+x[ceiling(obs.idx)])
    quan <- x[floor(obs.idx)]+(obs.idx-floor(obs.idx))*(x[ceiling(obs.idx)]-x[floor(obs.idx)])##linear interpolation
  }
  else if(obs.idx>Nb) quan <- x[Nb]
  else if(obs.idx<1) quan <- x[1]
  else
  {
    stop("sth wrong!")
  }
  return(quan)
}

LoB2 <- function(x)
{
  ##method2: Linett paper, remove the "0.5" in obs.idx
  ##x: a vecotr of values
  
  x <- na.exclude(sort(x))
  Nb <- length(x)
  obs.idx <- Nb*0.95##remove 0.5
  if(obs.idx%in%(1:Nb)) quan <- x[obs.idx]
  else if (obs.idx>1 & obs.idx<Nb)      
  {
    #quan <- 0.5*(x[floor(obs.idx)]+x[ceiling(obs.idx)])
    quan <- x[floor(obs.idx)]+(obs.idx-floor(obs.idx))*(x[ceiling(obs.idx)]-x[floor(obs.idx)])##linear interpolation
  }
  else if(obs.idx>Nb) quan <- x[Nb]
  else if(obs.idx<1) quan <- x[1]
  else {
    browser()
    stop("sth wrong!")
  }
  return(quan)
}

LoB3 <- function(x)
{
  ##method3: complement with negative values
  ##x: a vecotr of values
  
  x <- na.exclude(sort(x))
  Nb <- length(x)
  neg.x <- -1*x
  quan <- as.numeric(quantile(c(neg.x,x),prob=0.95))
  return(quan)
}

LoB4 <- function(x)
{
  ##method4: Armbruster paper, mean+1.645*sd of blank obs
  ##x: a vecotr of values
  quan <- as.numeric(mean(x,na.rm=T)+1.645*sd(x,na.rm=T))
  return(quan)
}

LoB.Est.method <- function(x)
{
  ###wrap for all the four methods
  est.quan1 <- LoB(x)
  est.quan2 <- LoB2(x)
  est.quan3 <- LoB3(x)
  est.quan4 <- LoB4(x)
  return(c(method1=est.quan1,method2=est.quan2,method3=est.quan3,method4=est.quan4))
  
}


my_theme <- theme_bw() +
  theme(
    axis.text.x  = element_text(angle=0, vjust=0.5, size=12),   # rotate x-axis label and define label size
    # axis.ticks.x  = element_blank(), # x-axis ticks
    axis.text.y  = element_text(angle=0, vjust=0.5, size=12),   # rotate y-axis label and define label size
    axis.title.x = element_text(colour="black", size=16),  # x-axis title size
    axis.title.y = element_text(colour="black", size=16),  # y-axis title size
    plot.title = element_text(face="bold", colour="black", size=20, hjust=0.5, vjust=0.5),   # Main title size
    
    strip.text = element_text(size=8), 
    
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(), 
    
    legend.key = element_rect(colour = NA, fill = 'transparent'), 
    legend.background = element_rect(colour = NA,  fill = "transparent"),
    # legend.position = c(0.2, 0.3),  
    legend.position = "none",  # to remove all legends
    # legend.title = element_blank(),
    legend.title = element_text(colour="black", size=16),
    legend.text = element_text(colour="black", size=14),
    legend.text.align = 0, 
    legend.box = NULL
  )


regression_eqn_r2 <- function (summary_data) { # obtain the regression equation
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(summary_data$rlm_intercept, digits = 2), 
                        b = format(summary_data$rlm_slope, digits = 2), 
                        r2 = format(summary_data$rlm_Rsquared, digits = 2)))
  as.character(as.expression(eq));
}

regression_eqn_r <- function (summary_data) { # obtain the regression equation
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)~"="~r2, 
                   list(a = format(summary_data$rlm_intercept, digits = 2), 
                        b = format(summary_data$rlm_slope, digits = 2), 
                        r2 = format( sqrt(summary_data$rlm_Rsquared) , digits = 2)))
  as.character(as.expression(eq));
}