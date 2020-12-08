AFtoDelta <- function (AF, Rst) {
  R<- AF/(1-AF)
  delta<- (R/Rst-1)*1000
  delta
}

calcPercN15<- function(deltaTarg){
  RstN <- 0.003678
  AFadd <- 0.984
  AFinit <- RstN/(1+RstN)
  Rtarg <- ((deltaTarg/1000) +1)*RstN
  AFtarg <- (Rtarg/(1+Rtarg))
  (AFinit-AFtarg)/(AFtarg-AFadd)
}

# Assumptions: Dissolved N2 gas concentration at steady state before tracer experiment
predictNGas <- function(N15NO30, NO30, deltaNO30, B0Tot, B1Tot, precTot, B0Den, B1Den, precDen, k2, n, time = seq(0, 60*24, length.out = 50)){
  kTotMu <- B0Tot + B1Tot * log(NO30) 
  kTot <- rlnorm(n = length(kTotMu), meanlog = kTotMu, sdlog = sqrt(1/precTot))
  kDenMu <- B0Den + B1Den * log(kTot) 
  kDen <- rlnorm(n = length(kDenMu), meanlog = kDenMu, sdlog = sqrt(1/precDen))
  kAssim <- kTot-kDen
  delta <- matrix(NA,n,length(time))
  NGas <- matrix(NA,n,length(time))
  for(i in 1:length(time)){
    N15 <- (kDen*N15NO30/(k2 - kDen - kAssim))*(exp(-(kDen + kAssim)*time[i])- exp(-k2*time[i])) + AFinit * kDen*NO30/k2
    NTot <- (kDen*NO30/(k2 - kDen - kAssim))*(exp(-(kDen + kAssim)*time[i])- exp(-k2*time[i])) + kDen*NO30/k2
    delta[,i] <- AFtoDelta(N15/NTot, 0.003678)
    NGas[,i] <- NTot
  }
  return(list(delta = delta, NGas = NGas, kTot = kTot, kDen = kDen, kAssim = kAssim, deltaNO30 = deltaNO30, NO30 = NO30))
}

getPI <- function(ensemble, time){
  pi <- apply(ensemble$delta,2,quantile,c(0.025,0.5,0.975)) %>%
    t() %>%
    as_tibble() %>%
    mutate(time = time/60) %>%
    rename("Lower" = "2.5%", "Median" = "50%", "Upper" = "97.5%") %>% 
    mutate(delta = ensemble$deltaNO30, NO30 = ensemble$NO30)
  return(pi)
}

plotPI <- function(PI){
  plot <- ggplot(data = PI, aes(x = time)) +
    geom_ribbon(aes(x = time, ymin = Lower, ymax = Upper), fill = "grey80") + 
    geom_line(aes(y = Median)) +
    facet_grid(rows = vars(NO30), cols = vars(delta), labeller = labeller(.rows = label_both, .cols = label_both)) +
    xlab("Time (hours)") +
    ylab("Delta 15N-NGas")+
    theme_minimal()
  return(plot)
}