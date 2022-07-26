#COVID Air Quality Project R Code

library(sp)
library(gstat) 
library(fields)
library(classInt)
library(maps)
library(nlme) 
library(tidyverse)
library(lubridate)
library(spdplyr)

ploteqc <- function(spobj, z, breaks, ...){
  pal <- two.colors(n=length(breaks)-1, start="darkgreen", end="red", middle="yellow",
                    alpha=1.0)
  fb <- classIntervals(z, n = length(pal), 
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(spobj, col = col, ...)
  image.plot(legend.only = TRUE, zlim = range(breaks), col = pal)
}


aqidata <- read.csv("C:/Users/Chals/Downloads/monthy_aqi_ca.csv")
aqiCA <- aqidata %>% filter(StateCode == 6)
CAmonthly <- read.csv("C:/Users/Chals/Downloads/monthly")

siteData <- read.csv("C:/Users/Chals/Downloads/aqs_sites edited.csv")
siteData <- siteData %>% filter(State.Code == 6)

CAyrweather <- read.csv("C:/Users/Chals/Downloads/yearly_CA_weather.csv")
head(CAyrweather)

aqiweath <- read.csv("C:/Users/Chals/Downloads/avg_aqi_with_weather.csv")

###################
#create spatial objects
###################
joineddata <- left_join(aqidata, siteData, by = c("StateCode" = "State.Code", "CountyCode" = "County.Code", "SiteNumber" = "Site.Number"))
columns <- joineddata %>% dplyr::select(year = Year, month = ï..Month, AQI, Defining.Parameter, Defining.Site, StateCode , CountyCode, SiteNumber, Elevation, Land.Use, Location.Setting, Longitude, Latitude)

monitor_locations <- cbind(joineddata$Longitude, joineddata$Latitude)
CAData <- SpatialPointsDataFrame(monitor_locations, columns, 
                                 proj4string = CRS("+proj=longlat"),
                                 match.ID = TRUE)
plot(CAData, main = "monitors with observations")
maps::map("county", region = "california", add = TRUE)


columns <- siteData %>% dplyr::select(Longitude, Latitude, Elevation, State.Code, County.Code, Site.Number, Land.Use, Location.Setting)
monitor_locations <- cbind(siteData$Longitude, siteData$Latitude)
monitors <- SpatialPointsDataFrame(monitor_locations, columns, 
                                   proj4string = CRS("+proj=longlat"),
                                   match.ID = TRUE)


CAData <- CAData %>% mutate(Location.SettingNum = case_when(Location.Setting == "URBAN AND CENTER CITY" ~ 0,
                                                            Location.Setting == "SUBURBAN" ~ 1,
                                                            Location.Setting == "RURAL" ~ 2))
monitors <- monitors %>% mutate(Location.SettingNum = case_when(Location.Setting == "URBAN AND CENTER CITY" ~ 0,
                                                                Location.Setting == "SUBURBAN" ~ 1,
                                                                Location.Setting == "RURAL" ~ 2))
# 
# #############################
# % 
# % range(CAtemps$weather)
# % breaks = 30:95
# % ploteqc(CAtemps , CAtemps$weather, breaks, pch = 19) #plot sample 
# % maps::map("county", region = "california", add = TRUE)
# % title(main = "Average Annual Temperature, 2000-2021, degrees F")
# % 
# % range(CAwind$weather)
# % breaks = 1:11
# % ploteqc(CAwind , CAwind$weather, breaks, pch = 19) #plot sample 
# % maps::map("county", region = "california", add = TRUE)
# % title(main = "Average Annual Wind Speed, 2000-2021, Knots")
# % 
# % range(CApres$weather)
# % breaks = 775:2000
# % ploteqc(CApres , CApres$weather, breaks, pch = 19) #plot sample 
# % maps::map("county", region = "california", add = TRUE)
# % title(main = "Average Annual Barometric pressure, 2000-2021, Millibars")
# % 

######CA PM yearly averages#####


CApm <- CAData %>% filter(CAData$Defining.Parameter == "PM2.5")
range(CApm$AQI)
breaks <- 0:300
ploteqc(CApm, CApm$AQI, breaks,pch = 19, main = "CA PM2.5 AQI, 2010-2021")
maps::map("county", region = "california", add = TRUE)

CAoz <- CAData %>% filter(CAData$Defining.Parameter == "Ozone")
ploteqc(CAoz, CAoz$AQI, breaks,pch = 19)
maps::map("county", region = "california", add = TRUE, main= "CA Ozone AQI, 2010-2021")


####################################
#Likelihood estimation and kriging#
###################################

m <- nrow(monitors)
npm <- nrow(CApm)
noz <- nrow(CAoz)


years <- 2010:2021
for (i in years){
  print(i)
  datam <- CApm %>% filter(year == i)
  MLEandKriging(datam, monitors)
}


MLEandKriging <- function(CApm, monitors){
  
  
  lmpm <- lm(AQI ~ Longitude + Latitude + Elevation + Location.SettingNum, data = CAoz)
  summary(lmpm)
  
  fitted <- predict(lmpm, newdata = CAoz, na.action = na.pass)
  ehat <- CAoz$AQI - fitted #get MOM estimates
  
  
  
  #plot the fitted and the residuals
  #fitted
  jpeg(file=paste("saving_mlefit",i,".jpeg"))
  fittedbreaks <- seq(min(fitted, na.rm = TRUE), max(fitted, na.rm = TRUE), by = .025)
  ploteqc(CAoz, fitted, fittedbreaks, pch = 19)
  maps::map("county", region = "california", add = TRUE)
  title(main = paste(i, ", Fitted Mean Estimates"))
  dev.off()
  
  #residuals
  jpeg(file=paste("saving_mleres",i,".jpeg"))
  ehatbreaks <- seq(min(ehat, na.rm = TRUE),max(ehat, na.rm = TRUE), by = .025)
  ploteqc(CAoz, ehat, ehatbreaks, pch = 19)
  maps::map("county", region = "california", add = TRUE)
  title(main =  paste(i, ", Residuals"))
  dev.off()
  
  # From a subjective viewing of the fitted and residuals plots, 
  # you can see some areas of clustering. 
  # This indicates that there is spatial association 
  # and it would be best to conduct some spatial modeling. 
  
  #Nonparametric estimation of the variogram
  #using the MOM estimates from earlier
  CAoz$ehat <- ehat
  CApms <- CAoz[!is.na(ehat),] #removing NA values
  
  #view and plot the moment estimator points with bin size set by width
  vg <- variogram(ehat ~ 1, data = CApms, width = 25)
  
  #plot(vg, xlab = "Distance", ylab = "Semi-variogram estimate")
  fitvg <- fit.variogram(vg, vgm(0.5, "Exp", 100, 0.5))
  #print(fitvg)
  s2.hat <- fitvg$psill[2]
  rho.hat <- fitvg$range[2]
  tau2.hat <- fitvg$psill[1]
  #plot exponential variogram diagram
  expparavar <- plot(vg, fitvg, xlab = "Distance", ylab = "Semi-variogram estimate")
  #expparavar
  
  
  d <- rdist.earth(coordinates(CApms))
  n <- nrow(CApms)
  #cov matrix = s2*rho + t2*I
  gamma <- exp(-d/rho.hat)
  sigma <- s2.hat*gamma + tau2.hat*diag(n)
  invsigma <- solve(sigma)
  
  y <- CApms$AQI
  #create x matrix for observed locations
  x <- cbind(rep(1,n), CApms$Longitude, CApms$Latitude, CApms$Elevation, CApms$Location.SettingNum)
  
  
  beta.hat <- solve(t(x) %*% invsigma %*% x) %*% t(x) %*% invsigma %*% y
  #beta.hat
  
  #create new distance matrix with unobserved locations
  dcross <- rdist.earth(coordinates(CApms), coordinates(monitors))
  
  #cross-covariances
  gamma2 <- exp(-dcross/rho.hat)
  sigmacross <- s2.hat * gamma2
  
  #other componants
  xpred <- cbind(1, coordinates(monitors), monitors$Elevation, monitors$Location.SettingNum)
  colnames(xpred) = c("Intercept", "lon", "lat", "elevation", "location")
  
  #krigging equation for the mean estimates
  ypred <- xpred %*% beta.hat + t(sigmacross) %*%
    solve(sigma, y - x %*% beta.hat)
  #plot
  jpeg(file=paste("saving_plot1",i,".jpeg"))
  ypredbreaks <- 0:300
  ploteqc(monitors, ypred, ypredbreaks, pch = 19)
  points(coordinates(CApms), pch = 4, cex = 2)
  title(main = paste(i, ", Predicted Average PM2.5 AQI, 2010-2021"))
  maps::map("county", region = "california", add = TRUE)
  dev.off()
  
  b <- t(xpred) - t(x) %*% solve(sigma, sigmacross)
  vpred <- s2.hat -
    diag(t(sigmacross) %*% solve(sigma, sigmacross) +
           t(b) %*% solve(t(x) %*% solve(sigma, x), b))
  sepred <- sqrt(vpred)
  
  #plot the standard errors
  jpeg(file=paste("saving_plot2",i,".jpeg"))
  sepredbreaks <- seq(.5, 12, length.out = 10)
  ploteqc(monitors, sepred, sepredbreaks, pch = 19)
  points(coordinates(CApms), pch = 4, cex = 2)
  title(main  = paste(i, ", S.E. of Interpolations using Kriging"))
  maps::map("county", region = "california", add = TRUE)
  dev.off()
  print(paste(i," done"))
  save(ploteqc, monitors, CAData, CApm, CAoz, fittedbreaks, fitted, ehatbreaks, ehat, vg, expparavar, ypred, ypredbreaks, sepredbreaks, sepred, file = "someplotsb.rdata")
}


#####################################3
#Point-wise regression#
##################################3
library(broom)
library(dplyr)
library(spdplyr)
library(ggplot2)
library(tidyr)
library(purrr)

#get data
aqidata <- read.csv("C:/Users/Chals/Downloads/monthy_aqi_ca.csv")

siteData <- read.csv("C:/Users/Chals/Downloads/aqs_sites edited.csv")
siteData <- siteData %>% filter(State.Code == 6)




###################
#create spatial objects
###################
joineddata <- left_join(aqidata, siteData, by = c("StateCode" = "State.Code", "CountyCode" = "County.Code", "SiteNumber" = "Site.Number"))
columns <- joineddata %>% dplyr::select(year = Year, month = ï..Month, AQI, Defining.Parameter, Defining.Site, StateCode , CountyCode, SiteNumber, Elevation, Land.Use, Location.Setting, Longitude, Latitude)

monitor_locations <- cbind(joineddata$Longitude, joineddata$Latitude)
CAData <- SpatialPointsDataFrame(monitor_locations, columns, 
                                 proj4string = CRS("+proj=longlat"),
                                 match.ID = TRUE)
plot(CAData, main = "monitors with observations")
maps::map("county", region = "california", add = TRUE)


columns <- siteData %>% dplyr::select(Longitude, Latitude, Elevation, State.Code, County.Code, Site.Number, Land.Use, Location.Setting)
monitor_locations <- cbind(siteData$Longitude, siteData$Latitude)
monitors <- SpatialPointsDataFrame(monitor_locations, columns, 
                                   proj4string = CRS("+proj=longlat"),
                                   match.ID = TRUE)


CAData <- CAData %>% mutate(Location.SettingNum = case_when(Location.Setting == "URBAN AND CENTER CITY" ~ 0,
                                                            Location.Setting == "SUBURBAN" ~ 1,
                                                            Location.Setting == "RURAL" ~ 2))
monitors <- monitors %>% mutate(Location.SettingNum = case_when(Location.Setting == "URBAN AND CENTER CITY" ~ 0,
                                                                Location.Setting == "SUBURBAN" ~ 1,
                                                                Location.Setting == "RURAL" ~ 2))

#seperate into pm2.5 and ozone observations
CApm <- CAData %>% filter(CAData$Defining.Parameter == "PM2.5")
range(CApm$AQI)
breaks <- 0:300
ploteqc(CApm, CApm$AQI, breaks,pch = 19, main = "CA PM2.5 AQI, 2010-2021")
maps::map("county", region = "california", add = TRUE)

CAoz <- CAData %>% filter(CAData$Defining.Parameter == "Ozone")
ploteqc(CAoz, CAoz$AQI, breaks,pch = 19)
maps::map("county", region = "california", add = TRUE, main= "CA Ozone AQI, 2010-2021")


CADatadf <- as.data.frame(CAData)
CADatadf$date <- paste(CADatadf$year, CADatadf$month)

ca <- map_data("county", "california")

CApmdf <- as.data.frame(CApm)
CAozdf <- as.data.frame(CAoz)
#jpeg(file=paste("yearlyAQImaps.jpeg"))
#dev.new(width=5, height=4)

#yearly plots of observations
yearlyPMmap <- ggplot(CApmdf) +
  geom_point(aes(Longitude, Latitude, colour = AQI)) +
  geom_polygon(data=ca,aes(x=long,y=lat, group = group),
               colour='black', fill=NA, lwd=.0005) +
  scale_colour_distiller( palette="Spectral", limits = c(0,300)) +
  facet_wrap(~year, dir = "v", ncol = 4) + # facet by date    +
  labs(title="PM2.5 AQI data by year")
#dev.off()
yearlyozmap <- ggplot(CAozdf) +
  geom_point(aes(Longitude, Latitude, colour = AQI)) +
  geom_polygon(data=ca,aes(x=long,y=lat, group = group),
               colour='black', fill=NA, lwd=.0005) +
  scale_colour_distiller( palette="Spectral", limits = c(0,300)) +
  facet_wrap(~year, dir = "v", ncol = 4) + # facet by date    +
  labs( title= "Ozone AQI data by year")

save(yearlyPMmap, yearlyozmap, file = "yearlymaps.rdata")

#lm function
fit_one_pixel <- function(data)
  mod <- lm(AQI ~ 1 + year + month, data = data)

#use data up until 2019
til99 <- CApmdf %>% filter(year <= 2019)
groups <- til99 %>%
  filter(!is.na(AQI)) %>%
  group_by(Longitude, Latitude) %>%
  nest()

#group the data by monitor and run the lm for each monitor to get the model
groups <- groups %>% 
  mutate(model = map(data, fit_one_pixel)) %>%
  mutate(model_df = map(model, tidy))

lm_pars <- groups %>%
  unnest(model_df)

aqi_pred <- matrix(NA, nrow = 12, ncol = 2)
aqi_pred[,1] <- rep(2020, 12)
aqi_pred[,2] <- 1:12
colnames(aqi_pred) <- c("year","month")

#predictions function
predict_one_pixel <- function(lm, aqi_pred) {
  predict(lm, # linear model
          newdata = aqi_pred, # pred. covariates
          interval = "prediction") %>% # output intervals
    data.frame() %>% # convert to df
    mutate(se = (upr-lwr)/(2 * 1.96)) %>% # comp pred. se
    select(fit, se) # return fit & se
}

#get predictions for each monitor
groups2020 <- groups %>%
  mutate(preds = map(model,
                     predict_one_pixel,
                     aqi_pred = as.data.frame(aqi_pred))) %>%
  unnest(preds)

#plots of predictions
pred <- ggplot(groups2020) +
  geom_point(aes(Longitude, Latitude,
                 colour = fit)) +
  geom_polygon(data=ca,aes(x=long,y=lat, group = group),
               colour='grey', fill=NA, lwd=.0005) +
  scale_colour_distiller( palette="Spectral", limits = c(0,300)) + 
  theme_bw() + coord_fixed() + # fix scale and theme 
  labs(title = "Predicted Ozone AQI, 2020")

predse <- ggplot(groups2020) +
  geom_point(aes(Longitude, Latitude,
                 colour = se)) +
  geom_polygon(data=ca,aes(x=long,y=lat, group = group),
               colour='grey', fill=NA, lwd=.0005) +  scale_colour_distiller( palette="Spectral", limits = c(-50,50)) + 
  theme_bw() + coord_fixed() + # fix scale and theme 
  labs(title = "Predicted Ozone SE, 2020")


real <- ggplot(filter(CApmdf, year == 2020)) +
  geom_point(aes(Longitude, Latitude,
                 colour = AQI)) +
  facet_wrap(~year, dir = "v") + # facet by date
  geom_polygon(data=ca,aes(x=long,y=lat, group = group),
               colour='grey', fill=NA, lwd=.0005) +
  scale_colour_distiller( palette="Spectral", limits = c(0,300)) + 
  theme_bw() + coord_fixed() + # fix scale and theme 
  labs(title = "Actual Ozone AQI, 2020")

save(pred, predse, real, file = "ozplots.rdata")
# save(pred, predse, real, file = "pmplots.rdata")

#save(yearlyPMmap, yearlyozmap, pred, predse, real, file = "moreplots.rdata")

#get the differences
counts <- CApmdf %>% group_by(Longitude, Latitude) %>% summarize(n = n())
CAoz2020 <- CApmdf %>% filter(year == 2020) %>% select(year, month, Longitude, Latitude, AQI)
CAoz2020 <- CAoz2020 %>% group_by(Longitude, Latitude, year) %>% summarise(meanaqi = mean(AQI), n = n())
preds <- groups2020 %>% group_by(Longitude, Latitude) %>% summarise(meanfit = mean(fit), meanse = mean(se))
joined <- left_join(CAoz2020, preds, by = c("Longitude" = "Longitude", "Latitude" = "Latitude"))
joined <- joined %>% filter(!is.na(fit)) 
joined <- left_join(joined, counts, by = c("Longitude" = "Longitude", "Latitude" = "Latitude"))
diffs <- joined %>% mutate(diff = abs(meanaqi - meanfit))
hist(diffs$diff)
diffs <- diffs %>% filter(!is.na(diff))
#hist(diffs$diff)
pmdiffs <- ggplot(diffs) +
  geom_point(aes(Longitude, Latitude,
                 colour = diff, size = n.y, alpha = .2)) +
  geom_polygon(data=ca,aes(x=long,y=lat, group = group),
               colour='grey', fill=NA, lwd=.0005) +
  scale_colour_distiller( palette="Spectral", limits = c(0,100)) + 
  theme_bw() + coord_fixed() + # fix scale and theme 
  labs(title = "PM2.5 diffs")
# save(ozdiffs, file = "ozdiffsplot.rdata")
# save(pmdiffs, file = "pmdiffsplot.rdata")


##############################################3
#Spatio-temporal GLM#
####################################




##########################
library("ape")
library("dplyr")
library("FRK")
library("ggplot2")
library("gstat")
library("sp")
library("spacetime")
library("tidyr")
library(MASS)
library(nlme) 

CApmold <- filter(CApm, year <= 2021)
#using negative binomial instead of poisson to account for overdispersion
capm_GLM <- glm(AQI ~ (Longitude + Latitude + year)^2 + Elevation, # formula
                family = negative.binomial(1), # Poisson + log link
                data = CApmold) # data set


capm_GLM$deviance / capm_GLM$df.residual
# <1 shows no oversdispersion
capm_GLM$df.residual
capm_GLM$deviance
#significant no overdispersion
1 - pchisq(q = capm_GLM$deviance, df = capm_GLM$df.residual)


pmmonitors <- distinct(CApmold, Defining.Site, Longitude, Latitude, Elevation, year)
pred_grid <- as.data.frame(pmmonitors)

G <- auto_basis(data = pred_grid[,c("Longitude","Latitude")] %>% # Take Tmax
                  SpatialPoints(), # To sp obj
                nres = 1, # One resolutio
                type = "Gaussian")

S_pred <- eval_basis(basis = G, # basis functs
                     s = pred_grid[,c("Longitude","Latitude")] %>% # pred locs
                       as.matrix()) %>% # conv. to matrix
  as.matrix() # as matrix
colnames(S_pred) <- paste0("B", 1:ncol(S_pred)) # assign names
pred_grid <- cbind(pred_grid,S_pred) # attach to grid

aqi_preds <- predict(capm_GLM,
                     newdata = pred_grid,
                     type = "link",
                     se.fit = TRUE)

pred_grid <- pred_grid %>%
  mutate(log_cnt = aqi_preds$fit,
         se = aqi_preds$se.fit)


CApmold$residuals <- residuals(capm_GLM)

CApmolddf <- as.data.frame(CApmold)

glmres <- ggplot(CApmolddf) +
  geom_point(aes(Longitude, Latitude, colour = residuals)) +
  scale_colour_distiller(name = "residuals",   palette="Spectral") +
  geom_polygon(data=ca,aes(x=long,y=lat, group = group),
               colour='black', fill=NA, lwd=.0005) +
  facet_wrap(~year, nrow = 3) +
  labs(title = "PM2.5 Residuals")
#residuals don't look too noisy, indicated spatial correllation.

glmpreds <- ggplot(pred_grid) +
  geom_point(aes(Longitude, Latitude, colour = exp(log_cnt))) +
  scale_colour_distiller(name = "AQI",   palette="Spectral",  limits=c(0,300)) +
  geom_polygon(data=ca,aes(x=long,y=lat, group = group),
               colour='black', fill=NA, lwd=.0005) +
  labs(title = "PM2.5 Fitted Values") +
  facet_wrap(~year, nrow = 3) 


library(rdist)
CApmdf <- as.data.frame(CApm)
P <- list() # init list
years <- 2010:2020
for(i in seq_along(years)) { # for each day
  aqi_year <- filter(CApmold,
                     year == years[i]) # filter by year
  # obs_dists <- aqi_year %>% # take the data
  #   dplyr::select(Longitude,Latitude) %>% # extract coords.
  #   dist() %>% # comp. dists.
  #   as.matrix() # conv. to matrix
  obs_dists <- as.matrix(rdist(coordinates(aqi_year)))
  obs_dists[obs_dists == 0] <- .0001
  obs_dists.inv <- 1/obs_dists # weight matrix
  diag(obs_dists.inv) <- 0 # 0 on diag
  P[[i]] <- Moran.I(aqi_year$residuals, # run Moran's I
                    obs_dists.inv) %>%
    do.call("cbind", .) # conv. to df
}

MoransIsummary <- do.call("rbind",P) %>% summary(digits = 2)
#moran's I null hypothesis of no spatial correllation is rejected.

save(glmpreds, glmres, P, MoransIsummary, file = "glmplots.rdata")
#save(glmpreds, glmres, P, MoransIsummary, file = "glmplotsoz.rdata")