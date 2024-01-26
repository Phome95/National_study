knitr::opts_chunk$set(echo = TRUE)

library(tidyverse)  # Modern data science workflow
library(sf)
library(sp)
library(rgeos)
library(tmap)
library(tmaptools)
library(spgwr)
library(grid)
library(gridExtra)
library("dplyr")
library("ggpubr")
library(stats)
# Change the presentation of decimal numbers to 4 and avoid scientific notation
options(prompt="R> ", digits=4, scipen=999)
#read in mean biodiversity data #missing biodiversity data needs investigating
bio.data = read.csv("C:/Users/ezxph2/OneDrive - The University of Nottingham/PhD//DATA/Objective2/bio_mean.csv")
bio.data
ET_data =  st_read(file.path("C:/Users/ezxph2/OneDrive - The University of Nottingham/PhD/DATA/Objective2/gdf.shp"))
st_geometry_type(ET_data)
st_crs(ET_data)
Output.Areas2 <- read_sf("C:/Users/ezxph2/OneDrive - The University of Nottingham/PhD/DATA/Objective2/gdf.shp")
Output.Areas2 <- st_sf(Output.Areas2)
glimpse(Output.Areas2)
colnames(Output.Areas2)

Output.Areas2 = st_as_sf(Output.Areas2, crs=st_crs(Output.Areas2))
Output.Areas2
#map of Health index
qtm(Output.Areas2, fill = "HI_2021",title = "2021 Health Index")

#missing data, replace with median value
summary(Output.Areas2)
# install and/or load package
pacman::p_load(naniar)
# percent of ALL data frame values that are missing
pct_miss(Output.Areas2)
gg_miss_var(Output.Areas2, show_pct = TRUE)
#fill in missing vale with mean of its contiguous neighbours
index <- st_touches(Output.Areas2, Output.Areas2)

output <- output %>% 
  mutate(greenspace = ifelse(is.na(greenspace),
                             apply(index, 1, function(i){mean(.$greenspace[i])}),
                             greenspace))

output <- output %>% 
  mutate(protected_ = ifelse(is.na(protected_),
                             apply(index, 1, function(i){mean(.$protected_[i])}),
                             protected_))


summary(output)
Output.Areas2 =   output
summary(Output.Areas2)
# merge two data frames by LAD code
Output.Areas2<- merge(Output.Areas2,bio.data,by="LAD_code")
#rename total popultion column
names(Output.Areas2)[names(Output.Areas2) == "All person"] <- "population"
colnames(Output.Areas2)
#correlations
ET_data <- subset(Output.Areas2, select = -c(geometry))
ET_data$geometry <- NULL
#only get numeric data
ET_data = (ET_data[sapply(ET_data,is.numeric)])
library(corrplot)
testRes = cor.mtest(ET_data, conf.level = 0.95)
## leave blank on non-significant coefficient
## add significant correlation coefficients
corrplot(cor(ET_data,use="complete.obs"), p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)
corrplot(cor(ET_data[sapply(ET_data,is.numeric)],use="complete.obs"), method = 'square', order = 'FPC', type = 'lower', diag = FALSE)
#missing data
sum(is.na(Output.Areas2))
colSums(is.na(Output.Areas2))

colnames(Output.Areas2)


model = lm(HI_2021 ~ bmean_recordCount+ bmean_speciesCount+ pmean_recordCount + pmean_speciesCount +mmean_recordCount + Mmean_speciesCount+ areakm + GDP_millio + Pop_densit + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c + L1Observat + L2Observat + L3Observat
            + L4Observat + app_Observ +gdp_capita+ population,data = Output.Areas2)

print(shapiro.test(Output.Areas2$HI_2021))
par(mfrow=c(2,2))
plot(model)
summary(model)
#load car package
library(car)
library(performance)
library(see)

# Check model assumptions
pkgs <- c(
  "flextable", "performance", "see", "lmtest", "ggplot2",
  "qqplotr", "ggrepel", "patchwork", "boot","rempsyc"
)
install_if_not_installed(pkgs)
check_model(model)

#import
invisible(lapply(pkgs, library, character.only = TRUE))

nice_table(nice_assumptions(model), col.format.p = 2:4)
View(nice_assumptions(model))
#check VIF values
print(check_collinearity(model2))
#remove independent variables that produce a variance inflation factor above 5, but keep some to maintain information


model = lm(HI_2021 ~ bmean_recordCount+ bmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + Pop_densit + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c 
           + L4Observat + app_Observ +gdp_capita,data = Output.Areas2)


model2 = lm(HI_2021 ~ bmean_recordCount+ bmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + Pop_densit + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c 
            + L4Observat + app_Observ +gdp_capita,data = LAD)

model3 = lm(HI_2021 ~ mmean_recordCount+ Mmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + Pop_densit + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c 
            + L4Observat + app_Observ +gdp_capita,data = LAD)

plot(model3)
check_model(model2)
View(nice_assumptions(model3))
nice_table(nice_assumptions(model3), col.format.p = 2:4)
summary(model3)
#check VIF values
print(check_collinearity(model3))
library(AICcmodavg)

#define list of models
models <- list(model2, model3)

#specify model names
mod.names <- c('model2', 'model3')

#calculate AIC of each model
aictab(cand.set = models, modnames = mod.names)

#try gwr model
resids <- residuals(model2)
map.resids <- cbind(Output.Areas2, resids)
head(map.resids)
names(map.resids)[19] <- "resids"
qtm(map.resids, fill = "resids")
#gwr
OA.Census2 <- st_as_sf(Output.Areas2)
coords = st_coordinates(OA.Census2)
GWRbandwidth <- gwr.sel(HI_2021 ~ bmean_recordCount+ bmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c ,data = Output.Areas2 ,coords=coords, adapt = T)
#run gwr model
gwr.model = gwr(HI_2021 ~ bmean_recordCount+ bmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c ,
                data = OA.Census2,
                coords=coords,
                adapt=GWRbandwidth,
                hatmatrix=TRUE,
                se.fit=TRUE)
gwr.model
results <-as.data.frame(gwr.model$SDF)
names(results)
results
LMZ.F3GWR.test(gwr.model)
#map results
gwr.map <- cbind(OA.Census2, as.matrix(results))
gwr.map2 <- st_as_sf(gwr.map)
qtm(gwr.map, fill = "pred" , title = "Health index")

#tests for gwr
BFC99.gwr.test(gwr.model)
BFC02.gwr.test(gwr.model, approx=FALSE)
BFC02.gwr.test(gwr.model, approx=FALSE)
LMZ.F3GWR.test(gwr.model)

gwr.morantest(gwr.model, lw, zero.policy = FALSE)


#spatial auto correlation
library(spdep)
neighbours <- poly2nb(Output.Areas2)
neighbours

#one area with no neighbours (isle of white, number 44, so remove to undertake spatial weights)
LAD <- Output.Areas2[-c(42), ]
print(Output.Areas2["LAD22NM"], n=44)

#define neighbours

nb <- poly2nb(LAD, queen=TRUE)
nb[[1]]
#assign weights to nearest polygons, equal weight = w, 
lw <- nb2listw(nb, style="W", zero.policy=TRUE)
print(lw,zero.policy=TRUE)
#compute neighbourhood health index value
health.lag <- lag.listw(lw, LAD$HI_2021)
#create regression model
# Create a regression model
M <- lm(health.lag ~ LAD$HI_2021)

# Plot the data
plot(health.lag ~ LAD$HI_2021, pch=20, asp=1, las=1)

moran.test(LAD$HI_2021, lw)
#compute morans I with monte carlo simulation

MC<- moran.mc(LAD$HI_2021, lw, nsim=599)

# View results (including p-value)
MC
#spatial autocorrelation
#and for future use, write the residuals out to a column in your dataframe

LAD$model2_resids <- residuals(model2)
LAD_SP <- as(LAD,"Spatial")
coordsW <- coordinates(LAD_SP)
plot(coordsW)
#Now we need to generate a spatial weights matrix (remember from the lecture a couple of weeks ago). We'll start with a simple binary matrix of queen's case neighbours

#create a neighbours list of queens contiguity (LAD boundaries need fixing)
LADnb <- poly2nb(LAD_SP, queen=T)

#or nearest neighbours
knn_LAD <- knearneigh(coordsW, k=4)
LAD_knn <- knn2nb(knn_LAD)

#plot them
plot(LADnb, coordinates(coordsW), col="red")
plot(LAD_knn, coordinates(coordsW), col="blue")
#add a map underneath
centroids <- st_centroid(st_geometry(LAD))
plot(st_geometry(LAD), border = "grey60", reset = FALSE)
plot(LADnb, coords = centroids, add=T, col = "red")
#create a spatial weights matrix object from these weights
LAD.queens_weight <- nb2listw(LADnb, style="C") #lad boundaries need fixing
LAD.knn_4_weight <- nb2listw(LAD_knn, style="C")



#try to run linear models divided up by spatial regime
gwr.data<-gwr.model$SDF
gwr.data
names(gwr.data)
dat.dist<-dist(gwr.data@data[, 3:17])
#The hclust() function does basic hierarchical clustering
#according to several different methods, we'll use Ward's method
clust.dat<-hclust(dat.dist, method="ward.D")

#And we'll plot the dendrogram, or tree plot
par(mfrow=c(1,1))
plot(clust.dat)
#I only want a few groups, so I'll cut the tree so I get 4 clusters
gwr.data$clus<-cutree(clust.dat, k=4)


#i'll use table() to get the frequencies of each cluster
table(gwr.data$clus)
#to get it to plot right, we have to convert the cluster number
#to a factor variable
gwr.data$b.cf<-as.factor(gwr.data$clus)
spplot(LAD_SP)

spplot(gwr.data,"b.cf", col.regions=brewer.pal(4, "Accent"), main="GWR Spatial Regimes", col="transparent",sp.layout = LAD_SP)
#We can get the means of each coefficient for each of the regime areas
b.means<-aggregate(gwr.data@data[, 3:17], by=list(gwr.data$clus), mean)
b.means

#The main idea behind examining spatial regimes is that we are looking for areas where the model may be working differently, we've done this, now we will
#fit separate ols models to each regime
Output.Areas2$cluster<-gwr.data$clus
lm1<-lm(HI_2021 ~ Median_dis + bird_recor + bird_speci + mammal_rec + mammal_spe + plant_reco + plant_spec + areakm + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c ,data = Output.Areas2, subset=Output.Areas2$cluster==1)
lm2<-lm(HI_2021 ~ Median_dis + bird_recor + bird_speci + mammal_rec + mammal_spe + plant_reco + plant_spec + areakm + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c ,data = Output.Areas2, subset=Output.Areas2$cluster==2)
lm3<-lm(HI_2021 ~ Median_dis + bird_recor + bird_speci + mammal_rec + mammal_spe + plant_reco + plant_spec + areakm + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c ,data = Output.Areas2, subset=Output.Areas2$cluster==3)
lm4<-lm(HI_2021 ~ Median_dis + bird_recor + bird_speci + mammal_rec + mammal_spe + plant_reco + plant_spec + areakm + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c ,data = Output.Areas2, subset=Output.Areas2$cluster==4)

#examine each of the spatial regimes
summary(lm1)
summary(lm2)
summary(lm3)
summary(lm4)

#fit a spatial lag model as alternative
library(spatialreg)
#Maximum Likelihood (ML) estimation of the spatial lag model is carried out with the lagsarlm() function. 
#The required arguments are a regression “formula,” 
#a data set and a listw spatial weights object. The default method uses Ord’s eigenvalue decomposition of the spatial weights matrix
#spatial weights
library(spdep)
nb <- LADnb
#assign weights to nearest polygons, equal weight = w, 
lw <- nb2listw(nb, style="W", zero.policy=TRUE)
#fit spatial lag model
fit_2_lag <- lagsarlm(HI_2021 ~ mmean_recordCount+ Mmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + Pop_densit + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c 
                         + L4Observat + app_Observ +gdp_capita,data = LAD,lw)
summary(fit_2_lag)#compare AIC values to ols model, very similar (1932), rho value is negative and small and not significant, means on spatial model may be better
summary(model2)
AIC(model2)
AIC(fit_2_lag)
####Here, we use Monte Carlo simulation to obtain simulated distributions of the various impacts
W <- as(lw, "CsparseMatrix")
trMC <- trW(W, type="MC")
#get impacts from spatial lag model, r = number of simulations
im<-impacts(fit_2_lag, tr = trMC, R=100)

sums<-summary(im,  zstats=T)

data.frame(sums$res)
#To print the p values
data.frame(sums$pzmat)

#use spatial error model
#We coerce the sf object into a new sp object
LAD_SP <- as(LAD,"Spatial")
#Then we create a list of neighbours using the Queen criteria
LADnb <- poly2nb(LAD_SP, queen=T)
#create matrix
wm_s <- nb2mat(LADnb, style='B')
#create list of weights = equal weights
rwm_s <- mat2listw(wm_s, style='W')

fit_3_err <- errorsarlm(HI_2021 ~ mmean_recordCount+ Mmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + Pop_densit + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c 
                        + L4Observat + app_Observ +gdp_capita,data = LAD, rwm_s)
summary(fit_3_err)
summary(fit_2_lag)
summary(model3)
AIC(fit_3_err)
#Save AIC values
AICs<-c(AIC(model3),AIC(fit_2_lag), AIC(fit_3_err))
labels<-c("OLS", "SLM","SEM" )

flextable(data.frame(Models=labels, AIC=round(AICs, 2)))
#present results
fits <- list()
fits$ols <- model3
fits$slag <- fit_2_lag
fits$serr <- fit_3_err

summary(fit_3_err)
fit_3_err
#try spatial rf
names(LAD)
library(spatialRF)
library(randomForest)
library(randomForestExplainer)
# Extract predictor variables and response variable
# Modify this based on your data structure"
predictor <- LAD[, c("mmean_recordCount","Mmean_speciesCount", "pmean_recordCount", "pmean_speciesCount" ,"areakm", "Pop_densit" , "unemployme","greenspace","protected_","NDVI_mean","bluespace_","elevation_","temp_c" 
                      ,"L4Observat","app_Observ","gdp_capita")]
response <- LAD$HI_2021

#names of the response variable and the predictors
dependent.variable.name <- "HI_2021"
predictor.variable.names <- colnames(predictor)
# Calculate distance matrix
distance_matrix <- st_distance(LAD$geometry)
#coordinates of the cases
xy <- LAD[, c("X.x", "Y")]

#distance matrix
distance.matrix <- distance_matrix

#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0, 100, 200, 400, 800)
###reduce multicolinarity
preference.order <- c("mmean_recordCount","Mmean_speciesCount", "pmean_recordCount", "pmean_speciesCount" ,"areakm", "Pop_densit" , "unemployme","greenspace","protected_","NDVI_mean","bluespace_","elevation_","temp_c" 
                      ,"L4Observat","app_Observ","gdp_capita")

predictor.variable.names <- spatialRF::auto_cor(
  x = LAD[, predictor.variable.names],
  cor.threshold = 0.6,
  preference.order = preference.order
) %>% 
  spatialRF::auto_vif(
    vif.threshold = 10,
    preference.order = preference.order
  )
#random seed for reproducibility
random.seed <- 1

# Split the data into training and testing sets
set.seed(123)
train_indices <- sample(1:nrow(LAD), 0.8 * nrow(LAD))
train_data <- LAD[train_indices, ]
test_data <- LAD[-train_indices, ]
# Train a non spatial random forest model
str(LAD)
model.non.spatial <- spatialRF::rf(
  data = LAD,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  verbose = FALSE,
  n.cores = 4
)
# extracting positions of NA values 
print ("Row and Col positions of NA values") 
which(is.na(LAD), arr.ind=TRUE)
#impute na values with 0
columns_to_impute <- c("A", "C")
# Adjust parameters as needed
#########
fit_3_err <- errorsarlm(HI_2021 ~ mmean_recordCount+ Mmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + Pop_densit + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c 
                        + L4Observat + app_Observ +gdp_capita,data = LAD, rwm_s)
fit_2_lag <- lagsarlm(HI_2021 ~ mmean_recordCount+ Mmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + Pop_densit + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c 
                      + L4Observat + app_Observ +gdp_capita,data = LAD,lw)
model3 = model3
library(stargazer)

view(stargazer(model2,fit_2_lag,fit_3_err, title="Regression Results", align=TRUE,type="text") )

modelsummary(fits, stars = T, fmt=5)
plot_summs(fits,inner_ci_level = .9)
library(lme4)
library(jtools)
library(huxtable)
export_summs(fits, scale = FALSE)

health2021 = read.csv("C:/Users/ezxph2/OneDrive - The University of Nottingham/PhD//DATA/Objective2/health2021.csv")
colnames(Output.Areas2)
health2021 <- subset(health2021, select = -c(geometry))
colnames(health2021)
Output.Areas3<- merge(LAD,health2021,by="LAD_code")
#spatial autocorrelation for physical acitivity
names(Output.Areas3)[names(Output.Areas3) == 'Physical.activity..L1.'] <- 'activity'
#compute neighbourhood health index value


#try with physical activity

physical.lag <- lag.listw(lw, Output.Areas3$activity)
#create regression model
# Create a regression model
M <- lm(health.lag ~ Output.Areas3$activity)

# Plot the data
plot(physical.lag ~ Output.Areas3$activity, pch=20, asp=1, las=1)

moran.test(Output.Areas3$activity, lw)
#compute morans I with monte carlo simulation

MC<- moran.mc(Output.Areas3$activity, lw, nsim=599)

# View results (including p-value)
MC
qtm(Output.Areas3, fill = "activity",title = "Physical Activity % of Survey (Weighted)")
##
#fit models
m
model4= lm(activity ~ mmean_recordCount+ Mmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + Pop_densit + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c 
                    + L4Observat + app_Observ +gdp_capita,data = Output.Areas3)

fit_3_err1 <- errorsarlm(activity ~ mmean_recordCount+ Mmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + Pop_densit + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c 
                        + L4Observat + app_Observ +gdp_capita,data = Output.Areas3, rwm_s)
fit_2_lag1 <- lagsarlm(activity ~ mmean_recordCount+ Mmean_speciesCount+ pmean_recordCount + pmean_speciesCount + areakm + Pop_densit + unemployme + greenspace + protected_+NDVI_mean +bluespace_ + elevation_ + temp_c 
                      + L4Observat + app_Observ +gdp_capita,data = Output.Areas3,lw)
export_summs(model4,fit_2_lag1,fit_3_err1, scale = FALSE)
AICs<-c(AIC(model4),AIC(fit_2_lag1), AIC(fit_3_err1))
labels<-c("OLS", "SLM","SEM" )

flextable(data.frame(Models=labels, AIC=round(AICs, 2)))
hist(Output.Areas3$activity)
