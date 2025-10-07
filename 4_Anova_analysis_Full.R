# Load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  performance,
  ggeffects,
  partR2,
  cluster,        # clustering tools
  dendextend,     # dendrogram visualization
  ggeffects,      # marginal effects
  ggplot2,        # plotting
  emmeans,        # estimated marginal means
  multcompView    # compact letter display for multiple comparisons
)

# List of your original datasets and their names
files <- c("data_Full.csv", "data_75.csv", "data_50.csv")

# Create an empty list to store all six datasets
all_datasets <- list()

for(f in files){
  
  # Load the data
  df <- read.csv(f)
  
  # Extract LiDAR density name from file (between "data_" and ".csv")
  lidar_den <- gsub("data_|\\.csv", "", f)
  
  # Create dataset for meas_yr >= 1990
  df_1990 <- df[df$meas_yr >= 1990, ]
  df_1990$LiDAR_Den <- lidar_den
  df_1990$Filter_Year <- 1990
  
  # Create dataset for meas_yr >= 2000
  df_2000 <- df[df$meas_yr >= 2000, ]
  df_2000$LiDAR_Den <- lidar_den
  df_2000$Filter_Year <- 2000
  
  # Store in list
  all_datasets[[paste0(lidar_den,"_1990")]] <- df_1990
  all_datasets[[paste0(lidar_den,"_2000")]] <- df_2000
}

# Check
names(all_datasets)
# 1. Make sure all datasets have the same columns
all_cols <- unique(unlist(lapply(all_datasets, colnames)))
# 2. Add missing columns to each dataset and fill with NA
all_datasets_fixed <- lapply(all_datasets, function(df) {
  missing_cols <- setdiff(all_cols, colnames(df))
  if(length(missing_cols) > 0){
    df[missing_cols] <- NA
  }
  # Ensure the same column order
  df <- df[all_cols]
  return(df)
})

# 3. Combine datasets
data_all <- do.call(rbind, all_datasets_fixed)

# Optional: reset row names
rownames(data_all) <- NULL

# Check
str(data_all)

data_all <- transform(data_all,
                    StandAge_s = scale(Age,  center=TRUE, scale=TRUE),
                    Productivity_s = scale(SI, center=TRUE, scale=TRUE),
                    MaturityIndex_s = scale(IMAT, center=TRUE, scale=TRUE)
)

data_all = data_all%>%
  mutate(Productivity = case_when(
    Productivity=="Low-Prod"~"Low",
    Productivity=="High-Prod"~"High"))%>%
  mutate(
  Stand_Age = fct_relevel(Stand_Age, "<80yrs", "80 - 140yrs", "140 - 250yrs", ">250yrs"),
  Productivity = fct_relevel(Productivity,"Low", "High"),
  ogi_cluster=fct_relevel(ogi_cluster,"Low", "Moderate", "High","Very-high")
)

#Check if the data subset has signficant effect
test = lmer(scale(CExyz) ~ ogi_cluster*Filter_Year*LiDAR_Den +
       (1 | BEC),
     data = data_all, REML = TRUE)

anova(test)
rm(dfA)


### Mixed effects models
dfA <- data.frame()  # initialize empty results df
z=0
m_int <- list()      # initialize model list
for (i in unique(data_all$Filter_Year)) {
  for (j in unique(data_all$LiDAR_Den)) {
    mydata <- data_all %>%
      filter(AOI=="AOI", Filter_Year==i, LiDAR_Den==j)
    
    if (nrow(mydata) == 0) next  # skip empty combos
    
    ## ----------- Models -------------
    m1 <- lmer(scale(CExyz) ~ Stand_Age * Productivity + (1 | BEC), data=mydata, REML=TRUE)
    m1_red <- lm(scale(CExyz) ~ Stand_Age * Productivity, data=mydata)
    anv_df1 <- anova(m1, m1_red)
    anv1 <- anova(m1, ddf="Kenward-Roger")
    r2_1 <- r2_nakagawa(m1)
    icc_1 <- performance::icc(m1)
    
    m2 <- lmer(scale(CExyz) ~ StandAge_s * Productivity_s + (1 | BEC), data=mydata, REML=TRUE)
    m2_red <- lm(scale(CExyz) ~ StandAge_s * Productivity_s, data=mydata)
    anv_df2 <- anova(m2, m2_red)
    anv2 <- anova(m2, ddf="Kenward-Roger")
    r2_2 <- r2_nakagawa(m2)
    icc_2 <- performance::icc(m2)
    
    m3 <- lmer(scale(CExyz) ~ ogi_cluster * Productivity + (1 | BEC), data=mydata, REML=TRUE)
    m3_red <- lm(scale(CExyz) ~ ogi_cluster * Productivity, data=mydata)
    anv_df3 <- anova(m3, m3_red)
    anv3 <- anova(m3, ddf="Kenward-Roger")
    r2_3 <- r2_nakagawa(m3)
    icc_3 <- performance::icc(m3)
    
    m4 <- lmer(scale(CExyz) ~ MaturityIndex_s * Productivity_s + (1 | BEC), data=mydata, REML=TRUE)
    m4_red <- lm(scale(CExyz) ~ MaturityIndex_s * Productivity_s, data=mydata)
    anv_df4 <- anova(m4, m4_red)
    anv4 <- anova(m4, ddf="Kenward-Roger")
    r2_4 <- r2_nakagawa(m4)
    icc_4 <- performance::icc(m4)
    
    ## ----------- Results table -------------
    dfB <- data.frame(
      Model_type = c("Discrete", "Continuous", "Discrete", "Continuous"),
      Target = c("Stand age", "Stand age","IMAT","IMAT"),
      Target_Significance = c(anv1$`Pr(>F)`[1], anv2$`Pr(>F)`[1],
                              anv3$`Pr(>F)`[1], anv4$`Pr(>F)`[1]),
      Productivity = c(anv1$`Pr(>F)`[2], anv2$`Pr(>F)`[2],
                       anv3$`Pr(>F)`[2], anv4$`Pr(>F)`[2]),
      Interaction = c(anv1$`Pr(>F)`[3], anv2$`Pr(>F)`[3],
                      anv3$`Pr(>F)`[3], anv4$`Pr(>F)`[3]),
      BEC_variant = c(anv_df1$`Pr(>Chisq)`[2], anv_df2$`Pr(>Chisq)`[2],
                      anv_df3$`Pr(>Chisq)`[2], anv_df4$`Pr(>Chisq)`[2]),
      ICC = c(icc_1$ICC_adjusted, icc_2$ICC_adjusted, icc_3$ICC_adjusted, icc_4$ICC_adjusted),
      R2c = c(r2_1$R2_conditional, r2_2$R2_conditional, r2_3$R2_conditional, r2_4$R2_conditional),
      R2m = c(r2_1$R2_marginal, r2_2$R2_marginal, r2_3$R2_marginal, r2_4$R2_marginal),
      Min_Year = i,
      LiDAR_Den = j
    )
    
    ## format p-values
    format_p <- function(x) ifelse(x < 0.001, "<0.001", round(x, 3))
    dfB[, c("Target_Significance","Productivity","Interaction","BEC_variant")] <- 
      lapply(dfB[, c("Target_Significance","Productivity","Interaction","BEC_variant")], format_p)
    
    ## round ICC + R2
    dfB[, c("ICC","R2c","R2m")] <- round(dfB[, c("ICC","R2c","R2m")], 3)
    
    ## append
    dfA <- rbind(dfA, dfB)
    
    ## store models without overwriting
    m_int <- c(m_int, list(m1, m2, m3, m4))
    
  }
}


## Quick glance at results
anova(m_int[[1]], ddf = "Kenward-Roger")
anova(m_int[[2]], ddf = "Kenward-Roger")

anova(m_int[[3]], ddf = "Kenward-Roger")
anova(m_int[[4]], ddf = "Kenward-Roger")

plot( ggeffect(m_int[[1]], terms =  "Stand_Age") )|
plot( ggeffect(m_int[[1]], terms =  "Productivity") )|
plot( ggeffect(m_int[[1]], terms =  c("Stand_Age","Productivity")) )

plot( ggeffect(m_int[[3]], terms =  "ogi_cluster") )|
  plot( ggeffect(m_int[[3]], terms =  "Productivity") )|
  plot( ggeffect(m_int[[3]], terms =  c("ogi_cluster","Productivity")) )


###############################################################################
############ Summarizing all comparisons ######################################
res_all <- data.frame()

for (i in 1:length(m_int)) {
  
  if (i %in% c(1,5,9,13,17,21)) {
    # Stand_Age*Productivity
    em_stand <- emmeans(m_int[[i]], ~ Stand_Age*Productivity)
    cld_stand <- cld(em_stand, Letters = letters) %>%
      mutate(
        Stand_Age = factor(Stand_Age,
                           levels = c("<70yrs", "70 - 140yrs", "140 - 250yrs", ">250yrs")),
        Productivity = factor(Productivity,
                              levels = c("Low", "High"))
      ) %>% 
      arrange(Productivity)#Productivity, Stand_Age
    
    pred_stand <- ggeffect(m_int[[i]], terms = c("Stand_Age","Productivity")) %>%
      as.data.frame()
    pred_stand$signif <- cld_stand$.group
    pred_stand$Model_type <- "Productivity_AGE"
    pred_stand$Min_Year=dfA$Min_Year[i]
    pred_stand$LiDAR_Den=dfA$LiDAR_Den[i]
    
  } else if (i %in% c(3,7,11,15,19,23)) {
    # IMAT*Productivity
    em_stand <- emmeans(m_int[[i]], ~ ogi_cluster*Productivity)
    cld_stand <- cld(em_stand, Letters = letters) %>%
      mutate(
        ogi_cluster = factor(ogi_cluster,
                             levels = c("Low", "Moderate", "High", "Very-high")),
        Productivity = factor(Productivity,
                              levels = c("Low", "High"))
      ) %>% 
      arrange(Productivity)#Productivity,ogi_cluster 
    
    pred_stand <- ggeffect(m_int[[i]], terms = c("ogi_cluster","Productivity")) %>%
      as.data.frame()
    pred_stand$signif <- cld_stand$.group
    pred_stand$Model_type <- "Productivity_IMAT"
    pred_stand$Min_Year=dfA$Min_Year[i]
    pred_stand$LiDAR_Den=dfA$LiDAR_Den[i]
    
  } else {
    next
  }
  
  res_all <- rbind(res_all, pred_stand)
}

## I created one data set at a time: one for productivity, one for stand age, 
# for for the interaction between productive and stand age,
# and lastly, one for the interaction between productivity and IMAT clusters
# Prod 
# AgIm
# AgImvsProd

Mod_res = rbind(AgImvsProd,AgIm)
Mod_res = rbind(Mod_res,Prod)

## Organize Results
t<-Mod_res
colnames(t)
t[c(2:5)]<-round(t[c(2:5)],2)
t$mstd<-paste(t$predicted," (?",t$std.error,")",sep="")
t$let_sig<-trimws(t$signif)
t$Treatment = paste(t$x,"_",t$group,sep="")

t[, "Result"] <-paste( t[ ,"mstd"], t[ ,"let_sig"],sep="")

unique(t$Model_type)
names(t)
results<- spread(t[t$Model_type=="IMAT*Productivity",c("Min_Year","LiDAR_Den","Treatment","Result3")], Treatment, Result, fill=0) 

results<-results [order(results $Carb_pool),]

write.csv(results,"results.csv")

###############################################################################
# Plotting the models' results

# Stand_Age*Productivity
# 1. Compute estimated marginal means (EMMs) for Stand_Age*Productivity
em_stand <- emmeans(m_int[[1]], ~ Stand_Age*Productivity)
cld_stand <- cld(em_stand, Letters = letters)  # Add letters for significance
cld_stand <- cld_stand %>%
  mutate(
    Stand_Age = factor(Stand_Age, 
                       levels = c("<80yrs", "80 - 140yrs", "140 - 250yrs", ">250yrs")),
    Productivity = factor(Productivity, 
                          levels = c("Low", "High"))
  ) %>% 
  arrange(Productivity,Stand_Age)
# 2. Get ggeffect predictions
pred_stand <- ggeffect(m_int[[1]], terms = c("Stand_Age","Productivity"))

# 3. Merge letters with ggeffect data
pred_stand$signif <- cld_stand$.group
pred_stand = as.data.frame(pred_stand)
pal=c('#A1DE93',"#679A7D")
pt1 = pred_stand%>%
ggplot( aes(x = x, y = predicted, group=group, color=group)) +
  geom_point(size = 4,  position = position_dodge(width = dodge_width)) +#aes(shape=Productivity),
  geom_errorbar(aes(ymin = predicted -std.error, ymax = predicted+ std.error), size = 1, width = 0.25, position = position_dodge(width = dodge_width)) +
  geom_text(aes(label = signif, y = predicted + 0.05), color = "black",
            position=position_dodge(width=1), vjust=-0.5, hjust=0,size=4.5) +
  scale_color_manual(values = pal)  +
  labs(x = "Stand Age", y = "Predicted Response", color='Productivity') +
  theme_minimal(base_size = 14)+
  theme_bw()+
  theme(legend.position=c(0.65,0.10),legend.direction='horizontal', legend.title= element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", linewidth=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        legend.text=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey", linewidth=0.5, linetype="solid"),
        panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0), axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.x=element_text(size=12),#element_blank(),#element_text(size=10, angle=45, hjust = 1), 
        axis.title.y=element_blank(), 
        axis.text.y=element_text(size=12), legend.key=element_rect(fill="transparent", colour="transparent")) 
pt1
# Stand_Age
# 1. Compute estimated marginal means (EMMs) for Stand_Age
em_stand <- emmeans(m_int[[1]], ~ Stand_Age)
cld_stand <- cld(em_stand, Letters = letters)  # Add letters for significance
cld_stand <- cld_stand %>%
  mutate(
    Stand_Age = factor(Stand_Age, 
                       levels = c("<80yrs", "80 - 140yrs", "140 - 250yrs", ">250yrs"))
  ) %>% 
  arrange(Stand_Age)
# 2. Get ggeffect predictions
pred_stand <- ggeffect(m_int[[1]], terms = c("Stand_Age"))

# 3. Merge letters with ggeffect data
pred_stand$signif <- cld_stand$.group
pred_stand = as.data.frame(pred_stand)
pt2 = pred_stand%>%
  # 4. Plot with letters
  ggplot( aes(x = x, y = predicted)) +
  # geom_line(color = "steelblue", size = 1.2) +
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_point(size = 4,  position = position_dodge(width = dodge_width)) +#aes(shape=Productivity),
  geom_errorbar(aes(ymin = predicted -std.error, ymax = predicted+ std.error), size = 1, width = 0.25, position = position_dodge(width = dodge_width)) +
  geom_text(aes(label = signif, y = predicted + 0.05), color = "black",
            position=position_dodge(width=1), vjust=-0.5, hjust=0,size=4.5) +
  labs(x = "Stand Age", y = "Forest Strcutural Complexity (CExyz)", color='Productivity') +
  theme_minimal(base_size = 14)+
  theme_bw()+
  theme(legend.position="bottom",legend.direction='horizontal', legend.title= element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", linewidth=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        legend.text=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey", linewidth=0.5, linetype="solid"),
        panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0), axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.x=element_text(size=12),#element_blank(),#element_text(size=10, angle=45, hjust = 1), 
        axis.title.y=element_text(size=12), 
        axis.text.y=element_text(size=12), legend.key=element_rect(fill="transparent", colour="transparent")) 

pt2
# Productivity
# 1. Compute estimated marginal means (EMMs) for Productivity
em_stand <- emmeans(m_int[[1]], ~ Productivity)
cld_stand <- cld(em_stand, Letters = letters)  # Add letters for significance
cld_stand <- cld_stand %>%
  mutate(
    Productivity = factor(Productivity, 
                          levels = c("Low", "High"))
  ) %>% 
  arrange(Productivity)
# 2. Get ggeffect predictions
pred_stand <- ggeffect(m_int[[1]], terms = c("Productivity"))

# 3. Merge letters with ggeffect data
pred_stand$signif <- cld_stand$.group
pred_stand = as.data.frame(pred_stand)
pt3 = pred_stand%>%
  # 4. Plot with letters
  ggplot( aes(x = x, y = predicted)) +
  # geom_line(color = "steelblue", size = 1.2) +
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_point(size = 4,  position = position_dodge(width = dodge_width)) +#aes(shape=Productivity),
  geom_errorbar(aes(ymin = predicted -std.error, ymax = predicted+ std.error), size = 1, width = 0.25, position = position_dodge(width = dodge_width)) +
  geom_text(aes(label = signif, y = predicted + 0.05), color = "black",
            position=position_dodge(width=1), vjust=-0.5, hjust=0,size=4.5) +
  labs(x = "Stand Age", y = "", color='Productivity') +
  theme_minimal(base_size = 14)+
  theme_bw()+
  theme(legend.position="bottom",legend.direction='horizontal', legend.title= element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", linewidth=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        legend.text=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey", linewidth=0.5, linetype="solid"),
        panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0), axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.x=element_text(size=12),#element_blank(),#element_text(size=10, angle=45, hjust = 1), 
        axis.title.y=element_blank(), 
        axis.text.y=element_text(size=12), legend.key=element_rect(fill="transparent", colour="transparent")) 




wrap_plots(pt2|pt3|pt1) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")

# IMAT*Productivity
# 1. Compute estimated marginal means (EMMs) for IMAT*Productivity
em_stand <- emmeans(m_int[[3]], ~ ogi_cluster*Productivity)
cld_stand <- cld(em_stand, Letters = letters)  # Add letters for significance
cld_stand <- cld_stand %>%
  mutate(
    ogi_cluster = factor(ogi_cluster, 
                       levels = c("Low", "Moderate", "High", "Very-high")),
    Productivity = factor(Productivity, 
                          levels = c("Low", "High"))
  ) %>% 
  arrange(Productivity,ogi_cluster)
# 2. Get ggeffect predictions
pred_stand <- ggeffect(m_int[[3]], terms = c("ogi_cluster","Productivity"))

# 3. Merge letters with ggeffect data
pred_stand$signif <- cld_stand$.group
pred_stand = as.data.frame(pred_stand)
pal=c('#A1DE93',"#679A7D")
pt4 = pred_stand%>%
  # 4. Plot with letters
  ggplot( aes(x = x, y = predicted, group=group, color=group)) +
  geom_point(size = 4,  position = position_dodge(width = dodge_width)) +#aes(shape=Productivity),
  geom_errorbar(aes(ymin = predicted -std.error, ymax = predicted+ std.error), size = 1, width = 0.25, position = position_dodge(width = dodge_width)) +
  geom_text(aes(label = signif, y = predicted + 0.05), color = "black",
            position=position_dodge(width=1), vjust=-0.5, hjust=0,size=4.5) +
  scale_color_manual(values = pal)  +
  labs(x = "Stand Age", y = "Predicted Response", color='Productivity') +
  theme_minimal(base_size = 14)+
  theme_bw()+
  theme(legend.position=c(0.65,0.10),legend.direction='horizontal', legend.title= element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", linewidth=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        legend.text=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey", linewidth=0.5, linetype="solid"),
        panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0), axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.x=element_text(size=12),#element_blank(),#element_text(size=10, angle=45, hjust = 1), 
        axis.title.y=element_blank(), 
        axis.text.y=element_text(size=12), legend.key=element_rect(fill="transparent", colour="transparent")) 
pt4
# IMAT(ogi_cluster)
# 1. Compute estimated marginal means (EMMs) for IMAT(ogi_cluster)
em_stand <- emmeans(m_int[[3]], ~ ogi_cluster)
cld_stand <- cld(em_stand, Letters = letters)  # Add letters for significance
cld_stand <- cld_stand %>%
  mutate(
    ogi_cluster = factor(ogi_cluster, 
                         levels = c("Low", "Moderate", "High", "Very-high"))
  ) %>% 
  arrange(ogi_cluster)
# 2. Get ggeffect predictions
pred_stand <- ggeffect(m_int[[3]], terms = c("ogi_cluster"))

# 3. Merge letters with ggeffect data
pred_stand$signif <- cld_stand$.group
pred_stand = as.data.frame(pred_stand)
pt5 = pred_stand%>%
  # 4. Plot with letters
  ggplot( aes(x = x, y = predicted)) +
  # geom_line(color = "steelblue", size = 1.2) +
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_point(size = 4,  position = position_dodge(width = dodge_width)) +#aes(shape=Productivity),
  geom_errorbar(aes(ymin = predicted -std.error, ymax = predicted+ std.error), size = 1, width = 0.25, position = position_dodge(width = dodge_width)) +
  geom_text(aes(label = signif, y = predicted + 0.05), color = "black",
            position=position_dodge(width=1), vjust=-0.5, hjust=0,size=4.5) +
  labs(x = "Stand Age", y = "Forest Strcutural Complexity (CExyz)", color='Productivity') +
  theme_minimal(base_size = 14)+
  theme_bw()+
  theme(legend.position="bottom",legend.direction='horizontal', legend.title= element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", linewidth=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        legend.text=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey", linewidth=0.5, linetype="solid"),
        panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0), axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.x=element_text(size=12),#element_blank(),#element_text(size=10, angle=45, hjust = 1), 
        axis.title.y=element_text(size=12), 
        axis.text.y=element_text(size=12), legend.key=element_rect(fill="transparent", colour="transparent")) 

pt5
# Productivity
# 1. Compute estimated marginal means (EMMs) for Productivity
em_stand <- emmeans(m_int[[3]], ~ Productivity)
cld_stand <- cld(em_stand, Letters = letters)  # Add letters for significance
cld_stand <- cld_stand %>%
  mutate(
    Productivity = factor(Productivity, 
                          levels = c("Low-Prod", "High-Prod"))
  ) %>% 
  arrange(Productivity)
# 2. Get ggeffect predictions
pred_stand <- ggeffect(m_int[[3]], terms = c("Productivity"))

# 3. Merge letters with ggeffect data
pred_stand$signif <- cld_stand$.group
pred_stand = as.data.frame(pred_stand)
pt6 = pred_stand%>%
  # 4. Plot with letters
  ggplot( aes(x = x, y = predicted)) +
  # geom_line(color = "steelblue", size = 1.2) +
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_point(size = 4,  position = position_dodge(width = dodge_width)) +#aes(shape=Productivity),
  geom_errorbar(aes(ymin = predicted -std.error, ymax = predicted+ std.error), size = 1, width = 0.25, position = position_dodge(width = dodge_width)) +
  geom_text(aes(label = signif, y = predicted + 0.05), color = "black",
            position=position_dodge(width=1), vjust=-0.5, hjust=0,size=4.5) +
  labs(x = "Stand Age", y = "", color='Productivity') +
  theme_minimal(base_size = 14)+
  theme_bw()+
  theme(legend.position="bottom",legend.direction='horizontal', legend.title= element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", linewidth=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        legend.text=element_text(size=12),
        panel.background=element_rect(fill="white", colour="grey", linewidth=0.5, linetype="solid"),
        panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0), axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.x=element_text(size=12),#element_blank(),#element_text(size=10, angle=45, hjust = 1), 
        axis.title.y=element_blank(), 
        axis.text.y=element_text(size=12), legend.key=element_rect(fill="transparent", colour="transparent")) 

pt6

pt2 + coord_cartesian(ylim = c(-2.5, 1.25))
pt3 + coord_cartesian(ylim = c(-2.5, 1.25))
pt1 + coord_cartesian(ylim = c(-2.5, 1.25))
pt5 + coord_cartesian(ylim = c(-2.5, 1.25))
pt6 + coord_cartesian(ylim = c(-2.5, 1.25))
pt4 + coord_cartesian(ylim = c(-2.5, 1.25))

combined_pt = wrap_plots((pt2|pt3|pt1)/(pt5|pt6|pt4)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")+
  coord_cartesian(ylim = c(-2.5, 1.25))

ggsave("Paper_Images/Anova_Full.jpeg",combined_pt, height=8.5, width=13.5, units="in", dpi=300)


###### Correlation plots ############
# names(data)[c(17,18,19,20,21,22,25,29)] <- c("CVHeight", "FD","CRugos","MADHeight", "FHD", "CANCOV","SI", "IMAT")
lidar<- c("CExyz" ,"CVHeight", "FD","CRugos","MADHeight", "FHD", "CANCOV","Hght_Max" )
field<-c("Age","LateSuc_sp","Species_Div",  "MaxHgt"  , "Special_Att" ,  "MaxDBH" ,       
        "DBH_CV" , "Hght_CV" ,"DeadFallen_bio", "DeadStand_bio" , "Live_bio" , "LargeTreesDen", 
        "TreeDen","IMAT", "SI")

# Define color palette
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# Set layout for 2 rows x 3 columns **once**
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))  # adjust margins

for (i in unique(data_all$Filter_Year)) {
  for (j in unique(data_all$LiDAR_Den)) {
    
    mydata <- data_all %>%
      filter(AOI == "AOI", Filter_Year == i, LiDAR_Den == j)
    
    M <- mydata[, c(field, lidar, "PC1",'PC2','PC3')] %>%
      data.frame() %>%
      drop_na() %>%
      cor(method = "spearman") %>%
      round(2)
    
    cor.mtest <- function(mat, ...) {
      mat <- as.matrix(mat)
      n <- ncol(mat)
      p.mat <- matrix(NA, n, n)
      diag(p.mat) <- 0
      for (k in 1:(n - 1)) {
        for (l in (k + 1):n) {
          tmp <- cor.test(mat[, k], mat[, l], ...)
          p.mat[k, l] <- p.mat[l, k] <- tmp$p.value
        }
      }
      colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
      p.mat
    }
    
    p.mat <- mydata[, c(field, lidar, "PC1",'PC2','PC3')] %>% 
      data.frame() %>% 
      cor.mtest()
    
    # Plot correlation matrix
    corrplot(
      M[c(16,19,20,21,17,18), c(16,19,20,21,17,18)],
      method = "color",
      col = col(200),
      type = "lower",
      addCoef.col = "black",
      tl.col = "black",
      tl.srt = 45,
      p.mat = p.mat[c(16,19,20,21,17,18), c(16,19,20,21,17,18)],
      sig.level = 0.05,
      insig = "blank",
      diag = FALSE,
      cl.ratio = 0.4,
      tl.cex = 1.1,
      number.cex = 1.1
    )
    # Add title manually closer to the plot
    title(main = paste0("Year: ", i, " | LiDAR: ", j), line = 1, cex.main = 1.2)
    
  }
}

# Reset layout to 1 plot
par(mfrow = c(1,1))

