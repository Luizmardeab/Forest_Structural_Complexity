#PCA and Random forest analysis of complexity metrics in detecting
#field measured old-growth structural attributes
# install.packages('pacman')
rm(list = ls()); gc()
library(pacman) 


p_load(corrplot,ggplot2,cluster,
       devtools,patchwork,
       tidyverse,MASS,dplyr,
       exactextractr,psych,reshape2,
       randomForest, caret,Rdimtools,
       caretEnsemble,iml,fastshap,pdp,
       lmerTest,emmeans,Matrix,sjPlot,lme4,
       multcomp,multcompView, sf,data.table,dplyr
)
pckgs = c("iml","corrplot","car","cluster",
  "randomForest", "caret",
  "caretEnsemble",
  "lmerTest","emmeans"
)

for (i in pckgs){
  print(citation(i))
}

# Set work directory
setwd(setwd(paste0(getwd())))

# Get Full dataset with Field and LiDAR metrics
shapepath<-'C2_plot_FSC_Full.shp'
fplots<-st_read(shapepath)

names(fplots)
nrow(fplots)
unique(fplots$samp_id)
sort(unique(fplots$meas_yr))


### Separate the data into Study are and outside
AOI = st_read("G:/Conservation Solution Lab/People/Luizmar/PhD_Luizmar/Van_Is.shp")
st_crs(AOI)==st_crs(fplots)
AOI = st_transform(AOI,st_crs(fplots))
st_crs(AOI)==st_crs(fplots)
AOI$AOI = "AOI"

fplots = st_join(fplots,AOI[c("AOI")])
fplots$AOI[is.na(fplots$AOI)] = "Out"
nrow(fplots)

##
plot_inf = read.csv("G:/Conservation Solution Lab/People/Luizmar/PhD_Luizmar/Chapter_2/C2_plot_FSC.csv")
nrow(plot_inf)
nrow(fplots)
names(fplots)
names(plot_inf)
unique(fplots$samp_id)
fplots = fplots%>%left_join(plot_inf[c(1,3,16)], by= "samp_id")
###

LDR = st_read("G:/Conservation Solution Lab/People/Luizmar/PhD_Luizmar/Chapter_2/BC_LiDAR_UA.shp")
st_crs(LDR)==st_crs(fplots)
LDR = st_transform(LDR,st_crs(fplots))
st_crs(LDR)==st_crs(fplots)
names(LDR)

fplots = st_join(fplots,LDR[c("year")])
nrow(fplots)

## Remove plots with recorded Fire or Harvest post measurment (Hermozila et al. 2019)
fplots$harvest[is.na(fplots$harvest)] = 0
fplots$fire[is.na(fplots$fire)] = 0
fplots$harvest <- as.numeric(fplots$harvest)
fplots$fire <- as.numeric(fplots$fire)
fplots$meas_yr <- as.numeric(fplots$meas_yr)
fplots[c("meas_yr","harvest","fire")]
fplots = fplots%>%filter(!(meas_yr <= harvest | meas_yr <= fire))


### Consolidated cutblocks
CCUT = st_read("G:/Conservation Solution Lab/People/Luizmar/PhD_Luizmar/Chapter_3/VI_CCutBlocks.shp")
st_crs(CCUT)==st_crs(fplots)
CCUT = st_transform(CCUT,st_crs(fplots))
st_crs(CCUT)==st_crs(fplots)
names(CCUT)
CCUT$HARVEST__1
CCUT$PERCENT_CL

fplots = st_join(fplots,CCUT[c("HARVEST__1","PERCENT_CL")])
nrow(fplots)
fplots$HARVEST__1 <- as.numeric(fplots$HARVEST__1)
fplots$HARVEST__1[is.na(fplots$HARVEST__1)]=0
fplots = fplots%>%filter((meas_yr > HARVEST__1))
nrow(fplots)

## REmove duplicates
fplots$duplcts=duplicated(fplots$geometry)
duplic = unique(fplots$samp_id[fplots$duplcts==TRUE])
table(fplots$duplcts)

fplots = fplots %>%
  group_by(geometry) %>%
  mutate(max_yr = max(year)) %>%
  ungroup()%>%data.frame() %>%
  mutate(del = case_when(
    samp_id %in% duplic & year < max_yr ~ TRUE,
    !samp_id %in% duplic ~ FALSE,
    TRUE ~ FALSE  # catch-all for anything else
  ))

names(fplots)

for (i in duplic){
  #keep largest "year)
  print(fplots[fplots$samp_id==i,c(1:3,23,24,52,55)])
}

fplots = fplots%>%filter(del==FALSE, !is.na(CExyz))

fplots$duplcts=duplicated(fplots$samp_id)
table(fplots$duplcts)
duplic = unique(fplots$samp_id[fplots$duplcts==TRUE])
for (i in duplic){
  #keep largest "year)
  print(fplots[fplots$samp_id==i,c(1:3,23,24,52,55)])
}

fplots = fplots%>%filter(duplcts==FALSE)
nrow(fplots)
## Remove plots harvested after field data collection and before LiDar
fplots = fplots%>%filter(year>HARVEST__1)

# Remove BEC variants not present on the Island
data = fplots%>%filter(!BEC%in%c("CWHdm", "CWHds1", "CWHms1"))
# View result
head(data)
nrow(data)
sort(unique(data$meas_yr))

# data<-data%>%filter(meas_yr>1990)
hist(2019-data$meas_yr+data$age)
nrow(data)
data$age<-2025-data$meas_yr+data$age
nrow(data)


# nrow(data)
names(data)[c(3:15,27,28,30,31,32,38,39,40)]<- c("Age","LateSuc_sp","Species_Div","MaxHgt","Special_Att","MaxDBH","DBH_CV" ,
                                                 "Hght_CV","DeadFallen_bio","DeadStand_bio","Live_bio","LargeTreesDen",'TreeDen',
                                                 "CExyz","CVHeight","FD","CRugos","MADHeight","FHD","CANCOV",
                                                 "Hght_Max")
field<-names(data)[c(3:15)]#[c(37:40,42:49)]
lidar<-names(data)[c(27,28,30,31,32,38,39,40)]

data<-data[!is.na(data$Hght_Max),]
data<-data[!data$si==0,]

data<-data[!is.na(data),]
nrow(data)

# Calculate the difference between 'hgt' and 'zmax'
diff_hgt_zmax <- data$MaxHgt - data$Hght_Max
# Calculate IQR-based outliers on the difference
q1 <- quantile(diff_hgt_zmax, 0.25, na.rm = TRUE)
q3 <- quantile(diff_hgt_zmax, 0.75, na.rm = TRUE)
iqr <- q3 - q1

# Define outliers as points below Q1 - 1.5*IQR or above Q3 + 1.5*IQR
lower_bound <- q1 - 1.5 * iqr
upper_bound <- q3 + 1.5 * iqr
outliers_diff <- diff_hgt_zmax < lower_bound | diff_hgt_zmax > upper_bound

data$outlier<-outliers_diff
table(data$outlier)

plot(data$MaxHgt, data$Hght_Max, 
     main = "Difference (hgt - zmax) with Outliers Highlighted", 
     xlab = "Field (hgt)", 
     ylab = "Lidar (zmax)", 
     col = ifelse(data$outlier, "red", "black"),  # Color points based on 'outlier' column
     pch = 19)


data<-data[,c("samp_id", "meas_yr",field, lidar,"outlier",'si',"BEC","AOI","year")]

################# Plotting data histogram ##############################
# Assuming 'data' is your data frame and 'field' is a vector of column names
data[data$Live_bio==0,]
data<-data[!data$Hght_Max==0,]
data = data[!is.na(data$Hght_Max),]
data = data[!is.na(data$CVHeight),]
data = data[!is.na(data$CRugos),]
data = data[!is.na(data$CExyz),]

# Calculate the Maturity Index
nrow(data)
unique(data$year)
data <- data %>%filter(outlier=='FALSE')%>%
  mutate(
    IMAT = rowSums(
      across(all_of(field[-c(1,13)]), ~ . / quantile(., 0.9, na.rm = TRUE))/length(field[-c(1,13)]), 
      na.rm = TRUE
    )
  )

names(data)
field<-c(field[1:13],"si","IMAT")

# Obtain mean and standard deviation for scaling variables
data = data%>%st_drop_geometry()

data = data[data$Age!="-Inf",]
scmean<-scsd<-min<-max<-n <-c()
for (i in c(field, lidar)){
  scmean<-c(scmean,round(mean(data[,i], na.rm=TRUE),3))
  scsd<-c(scsd,round(sd(data[,i],na.rm=TRUE),3))
  min<-c(min,round(min(data[,i],na.rm=TRUE),3))
  max<-c(max,round(max(data[,i],na.rm=TRUE),3))
  n <- c(n, sum(!is.na(data[[i]])))  # Count non-NA values
}

scstat<-data.frame(Var=c(field, lidar), Mean=scmean, SD=scsd, MIN=min,MAX=max, N=n)
nrow(data)
data %>% filter(si!=0)%>%summary()
data %>% filter(si!=0)%>%nrow()
scstat$SE<-scstat$SD/sqrt(scstat$N)


scstat%>%mutate(col1 = paste(round(Mean,2)," \u00b1 ",round(SE,2), sep=""),
                col2 = paste("(",round(MIN,2),",",round(MAX,2),")", sep=""))%>%
  dplyr::select(Var, col1, col2) 

# Initialize an empty list to store histogram plots
namelist<-c(field,lidar)#
hst_plt <- list()
z <- 0

# Loop through each column name in namelist
for (i in 1:length(namelist)) {
  print(namelist[i])
  
  # Generate histogram plot for each variable
  pt1 <- data %>% filter(si!=0)%>%
    ggplot(aes(x = .data[[namelist[i]]])) +  # Access column by name
    geom_histogram(bins = 20, fill = "#999999", color = "black") +
    xlab(namelist[i]) + 
    ylab("Count")
  
  # Add the plot to the list
  z <- z + 1
  hst_plt[[z]] <- pt1
}

# Split the list of plots
first_14_plots <- hst_plt[1:15]
last_7_plots <- hst_plt[16:21]

# Arrange the first 14 plots in a 3x4 grid
first_plot_image <- wrap_plots(first_14_plots, ncol = 3, nrow = 5) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")&
  theme(legend.position="bottom",plot.margin = unit(c(0,30,0,0), "pt"),
        strip.text.x = element_text( size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold.italic" ),
        legend.direction='vertical', legend.title=element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", linewidth=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, linewidth=0.5),
        legend.text=element_text(size=12),
        panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
        panel.background=element_rect(fill="white", colour="grey", linewidth=1, linetype="solid"),
        panel.grid.major=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0,face='bold'), 
        axis.title.x=element_text(size=12),#axis.ticks.x=element_blank(),
        axis.text.x=element_text(size=12, angle=20, hjust = 1), 
        axis.title.y=element_text(size=12), 
        axis.text.y=element_text(size=12), legend.key=element_rect(fill="transparent", colour="transparent")) 

#ggsave("Paper_Images/Field_histo_15.png", first_plot_image, width = 12, height = 16, dpi = 300)


# Arrange the first 14 plots in a 3x4 grid
second_plot_image <- wrap_plots(last_7_plots, ncol = 2, nrow = 3) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")&
  theme(legend.position="bottom",plot.margin = unit(c(0,30,0,0), "pt"),
        strip.text.x = element_text( size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold.italic" ),
        legend.direction='vertical', legend.title=element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", linewidth=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, linewidth=0.5),
        legend.text=element_text(size=12),
        panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
        panel.background=element_rect(fill="white", colour="grey", linewidth=0.5, linetype="solid"),
        panel.grid.major=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0,face='bold'), 
        axis.title.x=element_text(size=12),#axis.ticks.x=element_blank(),
        axis.text.x=element_text(size=10, angle=20, hjust = 1), 
        axis.title.y=element_text(size=12), 
        axis.text.y=element_text(size=10), legend.key=element_rect(fill="transparent", colour="transparent")) 


#ggsave("Paper_Images/Lidar_histo_15.png", second_plot_image, width = 9, height = 12, dpi = 300)


##Create age classes before normalizing
# Define age classes
data<- data %>% filter(si!=0)
data <- data %>%
  mutate(
    Stand_Age = case_when(
      Age <= 70 ~ "<70yrs", 
      Age > 70 & Age <= 140 ~ "70 - 140yrs", 
      Age > 140 & Age < 250 ~ "140 - 250yrs",
      Age >= 250 ~ ">250yrs"
    ),
    Stand_Age = fct_relevel(
      Stand_Age, "<70yrs", "70 - 140yrs", "140 - 250yrs", ">250yrs"
    )
  )

# Define productivity classes based on 'si' value
data <- data %>%
  mutate(Productivity = if_else(si >= 21, "High-Prod", "Low-Prod"))%>%
  mutate(
    Productivity = fct_relevel(
      Productivity,
      "Low-Prod", "High-Prod"
    )
  )
table(data$Productivity)
# Combine Stand_Age and Productivity to create new age+productivity classes
data <- data %>%
  mutate(Age_Productivity = paste(Productivity,Stand_Age,  sep = "_"))

# Order Age_Productivity from young low productivity to old-growth high productivity
data <- data %>%
  mutate(
    Age_Productivity = fct_relevel(
      Age_Productivity,
      "Low-Prod_<70yrs", "High-Prod_<70yrs",
      "Low-Prod_70 - 140yrs", "High-Prod_70 - 140yrs",
      "Low-Prod_140 - 250yrs", "High-Prod_140 - 250yrs",
      "Low-Prod_>250yrs", "High-Prod_>250yrs"
    )
  )

table(data$Age_Productivity)
data<-data%>%
  mutate(
    Old_prod = if_else(Age>=250,Age_Productivity,Stand_Age)
  )%>%mutate(Old_prod=fct_relevel(
    Old_prod, "<70yrs", "70 - 140yrs", "140 - 250yrs",
    "Low-Prod_>250yrs", "High-Prod_>250yrs"
  ))

table(data$Age_Productivity, data$BEC)
# Perform k-means clustering on the ogi index to create 4 clusters
set.seed(19)  # For reproducibility
kmeans_result <- kmeans(data$IMAT, centers = 4, nstart = 20)

# Add the cluster labels as a new column
data <- data %>%
  mutate(ogi_cluster = factor(kmeans_result$cluster, 
                              labels = c(  "Moderate","Low","Very-high","High")))%>%#c("Moderate", "High", "Very-high", "Low")
  mutate(ogi_cluster=fct_relevel(
    ogi_cluster,"Low", "Moderate", "High","Very-high"
  ))

data%>%group_by(ogi_cluster)%>%summarise(test=mean(CExyz))


 ## BEC groups
data <- data %>%
  mutate(BECsubzone = gsub("[0-9]", "", BEC))

data <- data %>%
  mutate(BECzone = str_sub(str_remove_all(BEC, "[0-9]"), 1, 3))



##Create PCA of field measurments:
# PCA Scores (first 3 PCs)
data = data%>%filter(AOI=='AOI') #Keep only data inside the study area

pca_result <- prcomp(data[field[-c(14,15)]])

# Variance explained by each PCA
explained_variance <- summary(pca_result)$importance[2, ]
print("Variance Explained by each PCA:")
print(explained_variance)

# Get variable loadings (contributions of variables to PCs)
loadings <- pca_result$rotation
print("Loadings (contributions of each variable to PCs):")
print(loadings[, 1:3]) # Contributions to the first 3 PCs


data[,c("PC1",'PC2','PC3')] <- pca_result$x[, 1:3]
# Apply normalization to the `field` column(s) and combine with `lidar` column(s)

names(data)[c(17,18,19,20,21,22,25,29)] <- c("CVHeight", "FD","CRugos","MADHeight", "FHD", "CANCOV","SI", "IMAT")
lidar<- c("CExyz" ,"CVHeight", "FD","CRugos","MADHeight", "FHD", "CANCOV","Hght_Max" )
field[14]<-"IMAT"
field[15]<-"SI"
data = data[!is.na(data$MADHeight),]
M<- data[,c(field,lidar,"PC1",'PC2','PC3')] %>%#select(all_of(c(field, lidar,"PC1",'PC2','PC3')))%>%
  data.frame() %>% drop_na() %>% cor(method = "spearman") %>% round(2)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- data[,c(field, lidar,"PC1",'PC2','PC3')] %>%#select(all_of(c(field, lidar,"ogi")))%>%
  data.frame()%>%cor.mtest()
head(p.mat[, 1:5])

### Correlation plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# Plot the correlation matrix subset
corrplot(M[c(1:15),16:21], method="color", col=col(200),  
         type="full", # Type of plot
         addCoef.col="black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, # Text label color and rotation
         p.mat=p.mat[c(1:15),16:21], sig.level=c(0.05), insig="blank", 
         diag=TRUE, # Hide correlation coefficient on the principal diagonal,
         cl.ratio=0.4, # Adjust legend size
         tl.cex=0.9,
         number.cex=0.8)

# Plot the correlation matrix Lidar with PCA
corrplot(M[16:21,24:26], method="color", col=col(200),  
         type="full", # Type of plot
         addCoef.col="black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, # Text label color and rotation
         p.mat=p.mat[16:21,24:26], sig.level=c(0.05), insig="blank", 
         diag=TRUE, # Hide correlation coefficient on the principal diagonal,
         cl.ratio=0.8) # Adjust legend size

# Plot the correlation matrix field with PCA
corrplot(M[1:15,24:26], method="color", col=col(200),  
         type="full", # Type of plot
         addCoef.col="black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, # Text label color and rotation
         p.mat=p.mat[1:15,24:26], sig.level=c(0.05), insig="blank", 
         diag=TRUE, # Hide correlation coefficient on the principal diagonal,
         cl.ratio=0.8) # Adjust legend size
# Adjust the plot parameters to increase legend width
par(oma=c(0, 0, 0, 0)) # Add space on the right for a wider legend

corrplot(M[1:15, 24:26], 
         method="color", 
         col=col(200),  
         type="full", 
         addCoef.col="black", 
         tl.col="black", 
         tl.srt=45, 
         p.mat=p.mat[1:15, 24:26], 
         sig.level=0.05, 
         insig="blank", 
         diag=TRUE,
         cl.ratio=0.8) # Adjust legend size

# Plot the correlation matrix subset
corrplot(M[1:15,1:15], method="color", col=col(200),  
         type = 'lower',
         addCoef.col="black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, # Text label color and rotation
         p.mat=p.mat[1:15,1:15], sig.level=c(0.05), insig="blank", 
         diag=FALSE # Hide correlation coefficient on the principal diagonal,
)

# Plot the correlation matrix subset
corrplot(M[c(16,19,20,21,17,18),c(16,19,20,21,17,18)], method="color", col=col(200),  
         type="lower", # Type of plot
         addCoef.col="black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, # Text label color and rotation
         p.mat=p.mat[c(16,19,20,21,17,18),c(16,19,20,21,17,18)], sig.level=c(0.05), insig="blank", 
         diag=FALSE, # Hide correlation coefficient on the principal diagonal,
         cl.ratio=0.4, # Adjust legend size
         tl.cex=0.9,
         number.cex=0.8)
unique(data$year)



dat<-data%>%filter(AOI=="AOI", meas_yr>= 1990, outlier==FALSE)#, year<2020,AOI=="AOI", !BEC%in%c("CWHdm", "CWHds1", "CWHms1"))#

write.csv(dat, "data_Full.csv")
########################### End Correlation ####################################
# Random forest analysis
# Load pre-processed data
# data<-read.csv("data_Full.csv")
nrFolds <- 5 #Create fold to separe the data into calibration and validation
folds <- rep_len(1:nrFolds, nrow(dat))
folds <- sample(folds, nrow(dat))
################ Random Forest with caret
set.seed(7)
metric <- "MAE"
#Define custom RF model based on mtry and ntree tuning
customRF<-list(type="Regression", library="randomForest", loop=NULL)
customRF$parameters<-data.frame(parameter=c("mtry", "ntree"), class=rep("numeric", 2), label=c("mtry", "ntree"))
customRF$grid<-function(x, y, len=NULL, search="grid") {}
customRF$fit<-function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry=param$mtry, ntree=param$ntree, ...)
}
#Predict label
customRF$predict<-function(modelFit, newdata, preProc=NULL, submodels=NULL)
  predict(modelFit, newdata)
#Predict prob
customRF$prob<-function(modelFit, newdata, preProc=NULL, submodels=NULL)
  predict(modelFit, newdata, type="prob")
customRF$sort<-function(x) x[order(x[,1]),]
customRF$levels<-function(x) x$classes

control<-trainControl(method="repeatedcv", number=2, repeats=2)
tunegrid<-expand.grid(mtry=c(12,6,1), ntree=c(250, 500, 750))

# Define prediction function for the Random Forest model
pred_wrapper <- function(model, newdata) {
  predict(model, newdata)
}


# Check the result
head(dat)

# Initialize performance data frame
cm_per <- data.frame(
  targvar = character(0),
  rsq = numeric(0),
  sd_rsq = numeric(0),
  rmse = numeric(0),
  nrmse = numeric(0),
  sd_rmse = numeric(0),
  mae = numeric(0),
  stringsAsFactors = FALSE
)

# Empty list of models
cm_models <- list()

# Initialize a data frame to store SHAP values
shap_results <- data.frame()
datimp <- data.frame()
summary(dat)

for (i in 1:length(lidar)) {
  # Create dataset
  dataset <- cbind(dat[lidar[i]], dat[c(field[-c(14,15)])])  # Adjust columns as needed
  
  # Tune parameters
  set.seed(7)
  fit.rf1 <- train(as.formula(paste(lidar[i], "~ .")), data = dataset, method = customRF, 
                   metric = metric,
                   preProc = c("center", "scale", "BoxCox"),
                   tuneGrid = tunegrid, trControl = control)
  
  # Get the best model based on performance
  best_model <- fit.rf1$bestTune
  print(best_model)
  
  # Fit the global model
  glob_fit <- randomForest(as.formula(paste(lidar[i], "~ .")), data = dataset, 
                           mtry = best_model$mtry, ntree = best_model$ntree)
  
  # Calculate metrics
  cm_per <- rbind(cm_per, data.frame(
    targvar = lidar[i],
    rsq = mean(glob_fit$rsq),
    sd_rsq = sd(glob_fit$rsq),
    rmse = sqrt(mean(glob_fit$mse)),
    nrmse = sqrt(mean(glob_fit$mse)) / (max(glob_fit$y)-min(glob_fit$y))*100,
    sd_rmse = sd(sqrt(glob_fit$mse)),
    mae = MAE(glob_fit$y, glob_fit$predicted),
    stringsAsFactors = FALSE
  ))
  
  cm_models[[i]] <- glob_fit
  gc()
  
  for (k in 1:nrFolds) {
    # Actual split of the data
    shap_res <- data.frame()
    fold <- which(folds == k)
    data.train <- dataset[-fold, ]
    data.test <- dataset[fold, ]
    
    # Use best parameter to fit a global model on training data
    b_fit <- randomForest(as.formula(paste(lidar[i], "~ .")), data = data.train, 
                          mtry = best_model$mtry, ntree = best_model$ntree)
    
    # Calculate variable importance using IncNodePurity
    varim <- importance(b_fit) %>%
      as.data.frame() %>%
      rownames_to_column(var = "Variable") %>%
      mutate(importance = round(IncNodePurity / max(IncNodePurity) * 100, 2))%>%mutate(targvar=lidar[i], fold=k)
    
    predictor <- Predictor$new(b_fit, data = data.train[-c(1)], y = data.train[,1])
    
    shapley_values <- lapply(1:nrow(dataset), function(i) {
      Shapley$new(predictor, x.interest = dataset[i, -1])$results
    })
    
    # Combine all Shapley values into a single data frame
    shapley_df <- do.call(rbind, shapley_values)
    
    
    # For overall variable importance, calculate mean absolute SHAP values for each feature
    shapley_df <- shapley_df%>%separate(feature.value, into = c("feature", "value"), sep = "=", remove = FALSE) %>%
      mutate(value = as.numeric(value)) %>%
      dplyr::select(-feature.value)%>%mutate(targvar=lidar[i], fold=k)
    
    if (i==1&k==1){
      shap_results<-shapley_df
      datimp<-data.frame(varim)
      cm_per_f<- data.frame(
        targvar = lidar[i],
        rsq = mean(b_fit$rsq),
        sd_rsq = sd(b_fit$rsq),
        rmse = sqrt(mean(b_fit$mse)),
        nrmse = sqrt(mean(b_fit$mse)) / (max(b_fit$y)-min(b_fit$y))*100,
        nrmse_sd = sqrt(mean(b_fit$mse)) / sd(b_fit$y)*100,
        sd_rmse = sd(sqrt(b_fit$mse)),
        mae = MAE(b_fit$y, b_fit$predicted))
    }else{
      shap_results<-rbind(shap_results,shapley_df)
      datimp<-rbind(datimp,varim)
      cm_per_f <- rbind(cm_per_f, data.frame(
        targvar = lidar[i],
        rsq = mean(b_fit$rsq),
        sd_rsq = sd(b_fit$rsq),
        rmse = sqrt(mean(b_fit$mse)),
        nrmse = sqrt(mean(b_fit$mse)) / (max(b_fit$y)-min(b_fit$y))*100,
        nrmse_sd = sqrt(mean(b_fit$mse)) / sd(b_fit$y)*100,
        sd_rmse = sd(sqrt(b_fit$mse)),
        mae = MAE(b_fit$y, b_fit$predicted),
        stringsAsFactors = FALSE
      ))
    }
    gc()
    print (paste('var:', i, ' and fold:', k, sep=""))
  }
}
# Optionally, save the SHAP values to a CSV file for later use
# write.csv(shap_results, "shap_Full.csv", row.names = FALSE)

############Model Performance summary############
# Assuming SD is a column that provides the standard deviation for each metric
# If not, use RSQ_sd, RMSE_sd, and MAE_sd in place of SD
rnd <-2
model_perf<-cm_per_f %>% 
  group_by(targvar) %>% 
  summarise(
    RSQ_mean = mean(rsq), RSQ_sd = sd(rsq), 
    RMSE_mean = mean(nrmse), RMSE_sd = sd(nrmse),
    MAE_mean = mean(mae), MAE_sd = sd(mae), 
    n = n()
  ) %>%
  mutate(
    RSQ_se = RSQ_sd / sqrt(n),
    RMSE_se = RMSE_sd / sqrt(n),
    MAE_se = MAE_sd / sqrt(n),
    RSQ_CI = paste0(round(RSQ_mean,rnd), " (? ", round(RSQ_se, rnd),")", sep=""),
    RMSE_CI = paste0(round(RMSE_mean,rnd), " (? ", round(RMSE_se, rnd),")", sep=""),
    MAE_CI = paste0(round(MAE_mean,rnd), " (? ", round(MAE_se, rnd),")", sep=""),
    
    RSQ_lower.ci = RSQ_mean - qt(1 - (0.05 / 2), n - 1) * RSQ_se,
    RSQ_upper.ci = RSQ_mean + qt(1 - (0.05 / 2), n - 1) * RSQ_se,
    
    RMSE_lower.ci = RMSE_mean - qt(1 - (0.05 / 2), n - 1) * RMSE_se,
    RMSE_upper.ci = RMSE_mean + qt(1 - (0.05 / 2), n - 1) * RMSE_se,
    
    MAE_lower.ci = MAE_mean - qt(1 - (0.05 / 2), n - 1) * MAE_se,
    MAE_upper.ci = MAE_mean + qt(1 - (0.05 / 2), n - 1) * MAE_se
  ) %>%
  as.data.frame()

model_perf[,c("targvar","RSQ_CI", "RMSE_CI","MAE_CI")]
# write.csv(model_perf, "VI_LiDAR/Paper_Images/models_performance_R15.csv", row.names = FALSE)
write.csv(model_perf, "models_performance_Full.csv", row.names = FALSE)

##################### PCA Analysis ############################
# PCA Scores (first 3 PCs)
summary(data)
data = data[!is.na(data$CVHeight),]
for(j in c(1:6)){
  pca_result <- prcomp(data[lidar[-c(j, 7,8)]])
  print(paste("Target Variable:", lidar[j]))
  # Variance explained by each PCA
  explained_variance <- summary(pca_result)$importance[2, ]
  print("Variance Explained by each PCA:")
  print(explained_variance)
  
  
  pca_scores <- pca_result$x[, 1:3]
  for (i in 1:length(lidar)){
    dataset_pca <- cbind(target = data[lidar[i]], pca_scores)
  
    # Linear regression using the first 3 PCs
    lm_fit <- lm(as.formula(paste(lidar[i], "~ PC1 + PC2 + PC3 ")), data = dataset_pca)
  
    # Regression summary (variance explained by each PCA)
    summary(lm_fit)
  
    # Get R-squared and coefficients
    rsq <- round(summary(lm_fit)$r.squared,2)
    coefficients <- summary(lm_fit)$coefficients
    print(paste(lidar[i], " R-squared:", rsq))
    # print("Coefficients of each PCA in predicting target variable:")
    # print(coefficients)
  }
}


############Variable importance############
# Display the updated SHAP results
shap_results<-read.csv("shap_Full.csv")
head(shap_results)
i<-1
for(i in 1:length(field)){
  shap_results$unnorm_value[shap_results$feature==field[i]]<-round(shap_results$value[shap_results$feature==field[i]]*scstat$SD[scstat$Var==field[i]]+scstat$Mean[scstat$Var==field[i]],2)
  
}

# Summary statistics
summary_stats <- shap_results%>%
  group_by(feature,value,targvar, fold) %>%
  summarise(mean_shap = mean(abs(phi)), Direction=if_else(mean(phi)>0,"Positive","Negative"))%>%data.frame()
print(summary_stats)
summary_stats$normalized <- summary_stats$mean_shap / sum(summary_stats$mean_shap)

# # Melt the SHAP values data frame to long format
# VAriable importance using Random Forest output
cmplot<-c()

dat<-data.frame(
  targ_var = datimp$targvar,
  predictors = datimp$Variable,
  importance = datimp$IncNodePurity,
  percimp = datimp$importance
)

for (i in 1:length(lidar)){  
  # Plot the importance values
  VI<-dat%>%filter(targ_var==lidar[i])%>%
    group_by(predictors)%>%
    summarise( mean=mean(percimp),SD=sd(percimp), n=n())%>%
    mutate(se= SD / sqrt(n),
           lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
           upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)%>%as.data.frame()
  
  p1<-VI%>% mutate(Var = fct_reorder(predictors, mean)) %>%
    ggplot(
      aes(x=Var,y=mean)) + 
    geom_col() +
    geom_errorbar(aes(ymin=lower.ci , ymax=upper.ci), width=.3, size=0.5) +
    scale_color_manual(values =c("#006666","#009933","#339999"))+
    scale_fill_manual(values =c("#006666","#009933","#339999"))+
    #geom_point(position=pd,size=3, shape=3) +
    coord_flip() +
    ggtitle(paste(lidar[i]))+
    ylab("Relative Influence (%)") +
    xlab("") + #labs (title ="c)")
    theme_bw() +
    theme(legend.position="", legend.direction = "vertical",legend.title=element_blank(),
          legend.background = element_rect(fill="white", linetype="solid", size=0.25, colour="black"),
          legend.spacing.y=unit(0, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
          legend.text=element_text(size=10),
          panel.background=element_rect(fill="white", colour="grey", size=0.5, linetype="solid"),
          panel.grid.major=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
          panel.grid.minor=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
          plot.title=element_text(color="black", size=14, hjust=0), axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10), axis.text.x=element_text(size=10, hjust = 1),
          axis.text.y=element_text(size=10), legend.key=element_rect(fill="transparent", colour="transparent"))
  p1
  cmplot[[i]]<-p1
}


# VAriable importance using SHAP values
shapplot<-c()

for (i in 1:length(lidar)){  
  # Plot the importance values
  VI<-summary_stats%>%#filter(feature!="MaxHgt")%>%
    mutate(predictors=feature, value=normalized/max(normalized)*100)%>%#
    filter(targvar==lidar[i])%>%
    group_by(predictors)%>%
    summarise(mean=mean(abs(value)), SD=sd(abs(value)), n=n())%>%data.frame()%>%
    mutate(se= SD / sqrt(n),
           lower.ci = mean - qt(1 - (0.1 / 2), n - 1) * se,
           upper.ci = mean + qt(1 - (0.1 / 2), n - 1) * se)%>%as.data.frame()
  levels(VI$Var)
  p1<-VI%>% mutate(Var = fct_reorder(predictors, mean)) %>%
    ggplot(
      aes(x=Var,y=mean)) + 
    geom_col(fill="#009933", color="#009933") +
    geom_errorbar(aes(ymin=lower.ci , ymax=upper.ci), width=.3, size=0.5, color="Black") +
    # scale_color_manual(values =c("#006666","#009933","#339999"))+
    # scale_fill_manual(values =c("#006666","#009933","#339999"))+
    # #geom_point(position=pd,size=3, shape=3) +
    coord_flip() +
    ggtitle(paste(lidar[i]))+
    ylab("mean(|SHAP value|)") +
    xlab("") + #labs (title ="c)")
    theme_bw() +
    theme(legend.position=c(0.8,0.2), legend.direction = "vertical",legend.title=element_blank(),
          legend.background = element_rect(fill="white", linetype="solid", size=0.25, colour="black"),
          legend.spacing.y=unit(0, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
          legend.text=element_text(size=12),
          panel.background=element_rect(fill="white", colour="grey", size=0.5, linetype="solid"),
          panel.grid.major=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
          panel.grid.minor=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
          plot.title=element_text(color="black", size=14, hjust=0), axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12), axis.text.x=element_text(size=12, hjust = 1),
          axis.text.y=element_text(size=12), legend.key=element_rect(fill="transparent", colour="transparent"))
  p1
  shapplot[[i]]<-p1
}

# Assuming shapplot is your list of six ggplot objects
plot_list <-  shapplot

# Arrange in 2 columns and 4 rows, with labels from (a) to (g)
combined_plot <-wrap_plots(
  plot_list[[1]], plot_list[[2]],
  plot_list[[3]], plot_list[[4]],
  plot_list[[5]], plot_list[[6]],
  #plot_list[[7]], plot_list[[8]],
  # plot_list[[7]], plot_spacer(),
  ncol = 3
) + plot_annotation(tag_levels = "a", tag_suffix = ")")+
  plot_layout(guides = "collect")& theme(legend.position="bottom", legend.direction = "horizontal")

# Display the combined plot
print(combined_plot)

# Save the combined plot to a file
ggsave("Paper_Images/FSCs_Full.png", plot = combined_plot, width = 12, height = 9, dpi = 300)


 # for geom_quasirandom (nicer than jitter)
summary(abs(shap_results$value[shap_results$targvar=="CExyz"]))
feature_order <- rev(c("MaxHgt", "Live_bio", "MaxDBH",
                       "Age", "DeadStand_bio", "TreeDen",
                       "Species_Div", "DeadFallen_bio", "DBH_CV",
                       "Special_Att", "Hght_CV", "LargeTreesDen",
                       "LateSuc_sp"))
p2 = shap_results %>%
  filter(targvar == "CExyz") %>%#, feature!="MaxHgt"
  mutate(value = if_else(abs(value) > 1, 1, abs(value)))%>%
  mutate(feature = factor(feature, levels = feature_order)) %>%
  ggplot(aes(x = feature, y = phi*10)) +
  geom_violin(aes(group = feature), fill = "transparent", color = "transparent",
              trim = FALSE, scale = "width") +
  geom_quasirandom(aes(color = abs(value)), alpha = 0.6, size = 1.5, width = 0.2) +
  scale_color_gradient2(
    low = "yellow", mid = "white", high = "#009933",
    midpoint = 0.5, limits = c(0, 1), oob = scales::squish,
    breaks = c(0, 1),
    labels = c("Low", "High")
  ) +
  coord_flip() +
  labs(
    y = "SHAP value",
    x = "",
    color = "Prediction value"
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.2, 0.15),
    legend.direction = "vertical",
    legend.title = element_text(size=12),
    legend.background = element_rect(fill="white", linetype="solid", size=0.25, colour="black"),
    legend.spacing.y = unit(0, "mm"),
    panel.border = element_rect(colour="black", fill=NA, size=0.5),
    legend.text = element_text(size=11),
    panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"),
    panel.grid.major = element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
    panel.grid.minor = element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
    plot.title = element_text(color="black", size=14, hjust=0),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12, hjust = 0.5),
    axis.text.x = element_text(size=12, hjust = 0.5),
    axis.text.y = element_text(size=12),
    legend.key = element_rect(fill="transparent", colour="transparent")
  )

p2

plot_list[[1]]
ggsave("Paper_Images/CExyz_Full.png", plot = plot_list[[1]], width = 6, height = 6, dpi = 300)

p2 <- p2 +
  theme(
    axis.text.y = element_blank(),   # remove left labels
    axis.title.y = element_blank(),
    axis.ticks.y.left = element_blank(),  # remove left ticks
    axis.ticks.y.right = element_line(),  # show right ticks
    axis.text.y.right  = element_text(size=12)  # optional: style right labels
  )

p1p2 = wrap_plots(plot_list[[1]], p2, ncol = 2) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") +
  plot_layout(widths = c(1,1))
ggsave("Paper_Images/CExyz_SHAP.png", plot = p1p2 , width = 10, height = 7, dpi = 300)


############Partial Dependence plots ############
#Using shap values
unique(shap_results$feature)
head(shap_results)
feature_order <- c("MaxHgt", "Live_bio", "MaxDBH",
                       "Age", "DeadStand_bio", "TreeDen",
                       "Species_Div", "DeadFallen_bio", "DBH_CV",
                       "Special_Att", "Hght_CV", "LargeTreesDen",
                       "LateSuc_sp")

x_lim = c(45,900, 250, 210,100,1000)

pdp_full<-c()
z=0
for (i in feature_order[1:6]){
  z=z+1
  
  p1 = shap_results %>%
    filter(targvar == "CExyz", feature==i)%>%
  ggplot(aes(x = unnorm_value, y = phi*10)) +
    geom_point( fill = "#009933", color = "#009933", alpha=0.5) +
    # stat_summary(fun.y = mean, geom = "line", col = "red", size = 1, linetype='dashed')+
    geom_vline(xintercept =x_lim[z] , linetype = "dashed", color = "grey", size = 1, alpha=0.85)+
    geom_smooth(method="loess", size=0.9,se = TRUE,fill="red",color="red") +
    labs(
      y = "SHAP value",
      x = i#,
      # color = "Prediction value"
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.2, 0.15),
      legend.direction = "vertical",
      legend.title = element_text(size=12),
      legend.background = element_rect(fill="white", linetype="solid", size=0.25, colour="black"),
      legend.spacing.y = unit(0, "mm"),
      panel.border = element_rect(colour="black", fill=NA, size=0.5),
      legend.text = element_text(size=11),
      panel.background = element_rect(fill="white", colour="grey", size=0.5, linetype="solid"),
      panel.grid.major = element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
      panel.grid.minor = element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
      plot.title = element_text(color="black", size=14, hjust=0),
      axis.title.x = element_text(size=12),
      axis.title.y = element_text(size=12, hjust = 0.5),
      axis.text.x = element_text(size=12, hjust = 0.5),
      axis.text.y = element_text(size=12),
      legend.key = element_rect(fill="transparent", colour="transparent")
    )
  p1

  pdp_full[[z]]<-p1
}

pdpts = wrap_plots(pdp_full[1:6])+
  plot_annotation(tag_levels = "a", tag_suffix = ")")

ggsave("Paper_Images/PDP_SHAP.png", plot = pdpts , width = 12, height = 9, dpi = 300)

pdpts = wrap_plots(pdp_full[7:13])+
  plot_annotation(tag_levels = "a", tag_suffix = ")")

ggsave("Paper_Images/PDP_SHAP_lower.png", plot = pdpts , width = 12, height = 12, dpi = 300)



##############################################################################

###### SHAP analysis for each group ###############
lidar
cm_models[[1]]
dataset <- cbind(data[lidar[1]], data[c(field[-c(14,15)])])  # Adjust columns as needed

predictor <- Predictor$new(cm_models[[1]], data = dataset[-c(1)], y = dataset[,1])

shapley_values <- lapply(1:nrow(dataset), function(i) {
  Shapley$new(predictor, x.interest = dataset[i, -1])$results
})

# Combine all Shapley values into a single data frame
shapley_df <- do.call(rbind, shapley_values)

# Sum absolute Shapley values for each feature
feature_importance <- shapley_df %>%
  group_by(feature) %>%
  summarise(importance = mean(abs(phi))) %>%
  arrange(desc(importance))%>%data.frame()



table(data$ogi_cluster,data$Stand_Age)

print(feature_importance)
stage<- unique(data$ogi_cluster)
plotlistII<-c()
z=0
for(i in stage) {
  set.seed(7)
  shapley <- Shapley$new(predictor, x.interest = dataset[data$ogi_cluster==i, -1])  # Customize rows as needed
  shapley$results  # View SHAP values
  # Plot SHAP values for individual instances
  shapley$plot()
  
  # For overall variable importance, calculate mean absolute SHAP values for each feature
  shapley_importance <- shapley$results %>%
    # group_by(feature) %>%
    mutate(mean = mean(abs(phi)),Dir= ifelse(phi > 0, "Positive", "Negative")) %>%
    arrange(desc(mean))%>%data.frame()%>%mutate(targvar=lidar[1])
  
  
  p2<-shapley_importance%>% mutate(Var = fct_reorder(feature, mean),
      Dir = factor(Dir, levels = c("Positive", "Negative")))  %>%# Explicitly set both levels)
    ggplot(
      aes(x=Var,y=mean),colour = "black") + 
    geom_col(aes(fill =Dir)) +
    # geom_errorbar(aes(ymin=lower.ci , ymax=upper.ci), width=.3, size=0.5) +
    #scale_color_manual(values =c("Positive"="#006666","Negative"="#339999"))+
    scale_fill_manual(values =c("Positive"="#006666","Negative"="#339999"))+ #"#009933"
    #geom_point(position=pd,size=3, shape=3) +
    coord_flip() +
    ggtitle(paste(i, sep=""))+
    ylab("SHAP values (%)") +
    xlab("") + #labs (title ="c)")
    theme_bw() +
    theme(legend.position="",legend.title=element_blank(),
          legend.background = element_rect(fill="white", linetype="solid", linewidth=0.25, colour="black"),
          legend.spacing.y=unit(0, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
          legend.text=element_text(size=12),
          panel.background=element_rect(fill="white", colour="grey", linewidth=0.5, linetype="solid"),
          panel.grid.major=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
          panel.grid.minor=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
          plot.title=element_text(color="black", size=14, hjust=0), axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12), axis.text.x=element_text(size=12, hjust = 1),
          axis.text.y=element_text(size=12), legend.key=element_rect(fill="transparent", colour="transparent"))
  p2
  z=(z+1)
  plotlistII[[z]]<-p2
}


# Remove legends from all plots except plot 3
plotlistII[[1]] <- plotlistII[[1]] + theme(legend.position = "none")
plotlistII[[2]] <- plotlistII[[2]] + theme(legend.position = "none")
plotlistII[[4]] <- plotlistII[[4]] + theme(legend.position = "none")

# Ensure plot 3 has the legend at the bottom
plotlistII[[1]] <- plotlistII[[1]] + theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.direction = "horizontal"  # Arrange legend items horizontally
  )

# Combine the plots and collect only plot 3's legend
final_plot <- patchwork::wrap_plots(
  (plotlistII[[2]] | plotlistII[[1]] |
     plotlistII[[3]] | plotlistII[[4]]), 
  guides = "collect"
) & theme( axis.text.y =element_text(size=12),
          axis.text.x =element_text(size=12),legend.text=element_text(size=14),
          plot.title=element_text(color="black", size=14, hjust=0)) # Force legend collection at the bottom

# Print final plot
print(final_plot)

ggsave("Paper_Images/CIMAT_Shap.jpeg", height=6, width=18, units="in", dpi = 300)

#####################################################################################
######################### Plot Stand Age and IMAT vs FSC ###########################
# install.packages("ggpmisc")
library(ggpmisc)
data = read.csv("data_Full.csv")
data$Pred = predict(cm_models[[1]],data)


dat <- data.frame(matrix(NA, nrow = nrow(data), ncol = 23))
data$si = data$SI
field2=field
field2[15]="si"
z<-0
for(i in c(lidar,field2)){
  z=z+1
  names(dat)[z]<-i
  dat[,i]<-data[,i]*scstat$SD[scstat$Var==i]+scstat$Mean[scstat$Var==i]
}

dat$Pred = data$Pred


M<- dat[,c("DeadFallen_bio","MaxHgt","TreeDen","Age","IMAT","Pred","Live_bio","Species_Div", "si")] %>%#select(all_of(c(field, lidar,"PC1",'PC2','PC3')))%>%
  data.frame() %>% drop_na() %>% cor(method = "spearman") %>% round(2)
p.mat <- dat[,c("DeadFallen_bio","MaxHgt","TreeDen","Age","IMAT","Pred","Live_bio","Species_Div", "si")] %>%#select(all_of(c(field, lidar,"ogi")))%>%
  data.frame()%>%cor.mtest()


size = 12
p1 <- dat %>%
  ggplot(aes(scale(CExyz), scale(Pred))) +
  geom_point(shape = 21, color = "black", fill = "#56B4E9", size = 1) +
  geom_smooth(method = "lm", color = "#56B4E9", fill = "#56B4E9", size = 0.5) +
  labs(
    title = "",
    x = "Observed - FSC",
    y = "Predicted - FSC"
  ) +
  scale_x_continuous(breaks = sort(c(pretty(scale(dat$CExyz))))) +
  coord_cartesian(ylim = c(-3.5, 2)) +
  theme(
    axis.title.y = element_text(size = size, colour = "black"),
    axis.title.x = element_text(size = size, colour = "black"),
    axis.text.x  = element_text(size = size, colour = "black"),
    axis.text.y  = element_text(size = size, colour = "black")
  ) +
  # equation + R? label
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,        # simple linear regression
    parse = TRUE,
    label.x.npc = "right",  # place in top-right corner
    label.y.npc = 0.95
  )

p1


p2 <- dat %>%
  ggplot(aes(IMAT, scale(Pred))) +
  geom_point(shape = 21,  color = "black",fill = "#56B4E9", size = 1) +
  geom_smooth(method = "loess",  color = "#56B4E9",fill = "#56B4E9", size = 0.5) +
  geom_text(x = 0.25, y =1.75, label = paste("Spearman cor: ", M[42], ", pval: ", sprintf("%.3f", p.mat[42]), sep = ""),
            hjust = -0.1, vjust = -0.5, color = "black", size = 5) +
  labs(title = "", x = "IMAT - Maturity Index", y = paste("FSC - ", i, sep=""), size = "Standard Error") +
  scale_x_continuous(breaks = sort(c(pretty(dat$IMAT))))

p2 = p2+
  coord_cartesian(ylim = c(-3.5, 2)) +  
  theme(axis.title.y=element_text(size=size, colour="black"), 
        axis.title.x=element_text(size=size, colour="black"), 
        axis.text.x=element_text(size=size, colour="black"), 
        axis.text.y=element_text(size=size, colour="black"))

p2



pal=c('#A1DE93',"#679A7D")
p3 <- dat%>%mutate(Productivity=if_else(si>20,"High","Low"))%>%
  ggplot( aes(log(Age),scale(Pred)))+
  geom_rect(aes(xmin = 5.15, xmax = 5.55, ymin = -5, ymax = 2),
            fill = "lightblue", alpha = 0.75, inherit.aes = FALSE) +
  geom_point(shape = 21, aes( fill=Productivity),size=1,color = "black") +
  geom_smooth(method = "loess", aes(color=Productivity, fill=Productivity),size=0.5)+
  scale_color_manual(values = pal)  +
  scale_fill_manual(values = pal)  +
  geom_text(x =3.25, y = 1.75, label = paste("Spearman cor: ",M[33],", pval: ",sprintf("%.3f", p.mat[33]),sep=""),
            hjust   = -0.1, vjust=-0.5, color="black", size = 5)+
  labs(color="Productivity",fill="Productivity", shape="")+
  labs(title="" , x = "Log(Stand Age)", y = paste("FSC - ", i, sep=""), size="Standart Error")+
  scale_x_continuous(breaks = sort(c(pretty(log(dat$Age)))))
p3 = p3+
  theme(legend.position = c(0.75,0.1),legend.direction = "horizontal",
        legend.title=element_text(size=size),legend.text = element_text(size=size),
        legend.background = element_rect(fill="white", linetype="solid", size=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"),
        axis.title.y=element_text(size=size, colour="black"), 
        axis.title.x=element_text(size=size, colour="black"), 
        axis.text.x=element_text(size=size, colour="black"), 
        axis.text.y=element_text(size=size, colour="black"))

p3



p2 <- p2 + coord_cartesian(ylim = c(-4, 2))
p3 <- p3 + coord_cartesian(ylim = c(-4, 2))
wrap_plots((p1|p3))+
  plot_annotation(tag_levels = list(c("a)",'b)', "c)", "d)", "e)", "f)")))&
  theme(#legend.position="bottom",
    strip.text.x = element_text( size = size, color = "black", face = "bold"),
    strip.text.y = element_text(size = size, color = "black", face = "bold.italic" ),
    legend.direction='vertical', legend.title=element_text(size=size),legend.text = element_text(size=size),
    legend.background = element_rect(fill="white", linetype="solid", size=0.5, colour="black"),
    legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
    panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
    panel.background=element_rect(fill="white", colour="grey", size=1, linetype="solid"),
    panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
    panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
    plot.title=element_text(color="black", size=14, hjust=0,face='bold'), #axis.title.x=element_blank(),axis.ticks.x=element_blank(),
    #axis.text.x=element_blank(),#element_text(size=10, angle=45, hjust = 1), 
    #axis.title.y=element_text(size=12), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    legend.key=element_rect(fill="transparent", colour="transparent")) 
ggsave(paste("Paper_Images/Age_IMAT_",i,"_full.jpeg",sep=""), height=7, width=12, units="in", dpi=300)

p1&
  theme(#legend.position="bottom",
    strip.text.x = element_text( size = size, color = "black", face = "bold"),
    strip.text.y = element_text(size = size, color = "black", face = "bold.italic" ),
    legend.direction='vertical', legend.title=element_text(size=size),legend.text = element_text(size=size),
    legend.background = element_rect(fill="white", linetype="solid", size=0.5, colour="black"),
    legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
    panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
    panel.background=element_rect(fill="white", colour="grey", size=1, linetype="solid"),
    panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
    panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
    plot.title=element_text(color="black", size=14, hjust=0,face='bold'), #axis.title.x=element_blank(),axis.ticks.x=element_blank(),
    #axis.text.x=element_blank(),#element_text(size=10, angle=45, hjust = 1), 
    #axis.title.y=element_text(size=12), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    legend.key=element_rect(fill="transparent", colour="transparent")) 

ggsave("Paper_Images/PredvsObs_CExyz.jpeg", height=5, width=5, units="in", dpi=300)