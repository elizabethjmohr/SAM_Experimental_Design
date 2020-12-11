AFtoDelta <- function (AF, Rst) {
  R<- AF/(1-AF)
  delta<- (R/Rst-1)*1000
  delta
}

RtoDelta <- function (R, Rst) {
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
predictNGas <- function(N15NO30, NO30, NO3Pre, deltaNO30, B0Tot, B1Tot, precTot, B0Den, B1Den, precDen, k2, n, time = seq(0, 60*24, length.out = 50)){
  kTotMu <- B0Tot + B1Tot * log(NO30) 
  kTot <- rlnorm(n = length(kTotMu), meanlog = kTotMu, sdlog = sqrt(1/precTot))
  kDenMu <- B0Den + B1Den * log(kTot) 
  kDen <- rlnorm(n = length(kDenMu), meanlog = kDenMu, sdlog = sqrt(1/precDen))
  kAssim <- kTot-kDen
  N15Nitrate <- matrix(NA,n,length(time))
  N15NitrateFrac <- matrix(NA,n,length(time))
  NitrateTotal <- matrix(NA,n,length(time))
  N15Gas <- matrix(NA,n,length(time))
  NGasTotal <- matrix(NA,n,length(time))
  deltaN15NO3 <- matrix(NA,n,length(time))
  deltaN15NGas <- matrix(NA,n,length(time))
  AFPre <- 0.003664522
  for(i in 1:length(time)){
    N15Nitrate[,i] <- NO3Pre * (AFPre) + (N15NO30-AFPre*NO3Pre)*exp(-(kAssim + kDen) * time[i])
    N15NitrateFrac[,i] <- NO3Pre * (AFPre)/N15NO30 + (1-AFPre*NO3Pre/N15NO30)*exp(-(kAssim + kDen) * time[i])
    NitrateTotal[,i] <- NO3Pre + (NO30-NO3Pre)* exp(-(kAssim + kDen) * time[i])
    N15Gas[,i] <- AFPre*NO3Pre*kDen/k2+(N15NO30-AFPre*NO3Pre) * (kDen/(k2 - kDen - kAssim))*(exp(-(kDen + kAssim)*time[i])- exp(-k2*time[i])) 
    NGasTotal[,i] <- NO3Pre*kDen/k2 + (NO30-NO3Pre)*(kDen/(k2 - kDen - kAssim))*(exp(-(kDen + kAssim)*time[i])- exp(-k2*time[i])) 
    deltaN15NO3[,i] <- AFtoDelta(N15Nitrate[,i]/NitrateTotal[,i], 0.003678)
    deltaN15NGas[,i] <-  AFtoDelta(N15Gas[,i]/NGasTotal[,i], 0.003678)
  }
  return(list(N15Nitrate = N15Nitrate, deltaN15NO3 = deltaN15NO3, N15Gas = N15Gas, deltaN15NGas = deltaN15NGas, kTot = kTot, kDen = kDen, kAssim = kAssim, k2 = k2, deltaNO30 = deltaNO30, NO30 = NO30, NO3Pre = NO3Pre))
}

getPI <- function(ensemble, varName, time){
  pi <- apply(ensemble[[varName]],2,quantile,probs = c(0.025,0.5,0.975), na.rm = TRUE) %>%
    t() %>%
    as_tibble() %>%
    mutate(time = time/60) %>%
    rename("Lower" = "2.5%", "Median" = "50%", "Upper" = "97.5%") %>% 
    mutate(delta = ensemble$deltaNO30, NO3Pre = ensemble$NO3Pre, k2 = ensemble$k2)
  return(pi)
}

plotPI <- function(PI, xFacet = "NO3Pre", yFacet = "delta"){
  plot <- ggplot(data = PI, aes(x = time)) +
    geom_ribbon(aes(x = time, ymin = Lower, ymax = Upper), fill = "grey80") + 
    geom_line(aes(y = Median)) +
    xlab("Time (hours)") +
    theme_minimal()
  return(plot)
}