# install.packages('pacman')
rm(list = ls()); gc()
library(pacman) 

p_load(curl,progress,
       foreign,gdal,
       future,future.apply,
       geometry,
       lidR,
       raster,MASS,
       terra,readxl,
       devtools,
       tidyverse,
       exactextractr,
       sf,MASS,sp,Rdimtools,VIM
)

setwd(paste0(getwd()))
source_gist("524eade46135f6348140")

## Get all excel file combined
# datafolder<-"TSA/"
# files <- list.files(datafolder, pattern=".xlsx", full.names = TRUE)
# file_names <- list.files(datafolder, pattern=".xlsx", full.names = FALSE)
# 
# 
# 
# 
# ncol(dat)
# 
# for (i in 1:length(files)){
#   if (i==1){
#   dat<-read_excel(files[i], sheet = 2)
#   
#   dat$dat_ori <- stringr::str_remove(string = file_names[i], pattern = ".xlsx")
#   
#   }else{
#     dat1<-read_excel(files[i], sheet = 2)
#     dat1$dat_ori <- stringr::str_remove(string = file_names[i], pattern = ".xlsx")
#     dat<-rbind(dat,dat1)
#    
#   }
#   print(paste(file_names[i],length(drop)))
#   }
# 
# warnings()


#Get x and y coordinates from field plots and sample ID to select only usable plots
shapepath<-'G:/Conservation Solution Lab/People/Luizmar/PhD_Luizmar/Chapter_2/BC_data_UA.shp'
# shapepath<-"BC_data.shp"
fplots<-st_read(shapepath)
# # Plot the catalog with chunks displayed
unique(fplots$SAMP_ID)
plot_a=fplots
############### Prepare tree data ######################################
dat<-read.csv("dat_tree.csv")
names(dat)
dat$samp_id[1:10]

names(plot_a)[1]<-"samp_id"
plot_a$samp_id[1:10]
ldr<-unique(dat$samp_id)
fld<-unique(plot_a$samp_id)
ldr[ldr%in%fld]


dat<-dat[dat$samp_id %in% plot_a$samp_id & dat$meas_yr%in%plot_a$meas_yr,]
unique(dat$samp_id)

dat<-dat[(dat$meas_first=='Y'&dat$meas_last=='Y')|(dat$meas_first=='N'&dat$meas_last=='Y'),]

#Plot area
dat<-dat%>%left_join(plot_a, by = c("samp_id","meas_yr"))
nrow(dat)


dat$plt_area_ha<-1/dat$phf_tree
dat$count<-1


dat$area[is.na(dat$area)]<-dat$plt_area_ha[is.na(dat$area)]

###Check the distribution of a few tree features
#Age
ggplot(dat)+geom_histogram(aes(age_tot),bins = 60)+
  theme(legend.position = "",legend.direction = "horizontal",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12, colour="black"),
        axis.title=element_text(size=12,face="bold"))

#Decay BEC-based loss factor: % decay estimate based on tree species, risk group
ggplot(dat)+geom_histogram(aes(pct_dcy),bins = 60)+
  theme(legend.position = "",legend.direction = "horizontal",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12, colour="black"),
        axis.title=element_text(size=12,face="bold"))

#Species
ggplot(dat[dat$age_tot>250,]) +
  geom_col(aes(x = species, y = count)) +
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12, face = "bold")
  )

unique(dat$species)

print("Tree data prepared")

#################################################################################
### Empirical data processing for estimation of old-growth attributes from empirical data
### Tree height, Basal area, Number of dead standing trees, structural complexity
### biomass, understood density, canopy gap, and age

###############################################################################
# Initials:
##1 Density of trees
Tden<-dat%>%filter(dbh>=util)%>%group_by(samp_id, meas_yr) %>%#summarise(Den = sum(phf_tree)/length(unique(plot_no)),Den1 = mean(stemsha, na.rm = TRUE))%>%
  summarise(area=mean(area),den1=n()/mean(area),TDen = mean(stemsha, na.rm = TRUE)) %>%
  # mutate( TDen = Den1)%>%
  data.frame()%>%right_join(all_plt, by = c("samp_id",'meas_yr'))
nrow(Tden)
unique(duplicated(Tden$samp_id))
which(duplicated(Tden$samp_id))
Tden[c(1159, 1160, 1161),]
hist(Tden$den1)
hist(Tden$TDen)
Tden[Tden$samp_id=="07038CG000006",]
Tden = Tden%>%filter(!is.na(TDen))
Tden[Tden$samp_id=="07038CG000006",]
plot(Tden$den1,Tden$TDen)

##2 Density of large tree (dbh>100)
unique(is.na(dat$phf_tree))
dent <- dat%>%filter(dbh>=util)%>%group_by(samp_id, meas_yr)%>%summarise(count=n())%>%data.frame()
LTs<-dat%>%filter(dbh>=100)%>%group_by(samp_id, meas_yr)%>%
  summarise(LTcount=ifelse(is.na(n()),0,n()))%>%data.frame()%>%
  right_join(dent, by = c("samp_id",'meas_yr'))%>%mutate(LTcount=ifelse(is.na(LTcount),0,LTcount))%>%
  right_join(Tden, by = c("samp_id",'meas_yr'))%>%  mutate(LT_100=LTcount/area, LT_100_den1 =den1*(LTcount/count), LT_100_den =TDen*(LTcount/count))%>%
  mutate(LT_100_den = ifelse(is.na(LT_100_den), 0, LT_100_den))
nrow(LTs)
unique(duplicated(LTs$samp_id))

hist(LTs$LT_100)
hist(LTs$LT_100_den)
hist(LTs$LT_100_den1)
si<-LTs$samp_id[LTs$LT_100_den>200]
unique(dat$age_tot[dat$samp_id%in%si])
##3 Total Biomass: High stand volume or biomass
Vol<-dat%>%filter(dbh>=util)%>%group_by(samp_id, meas_yr)%>%summarise(vol=sum(wsvha,na.rm = TRUE))%>%
  data.frame()%>%right_join(all_plt, by = c("samp_id",'meas_yr'))%>%
  mutate(vol = ifelse(is.na(vol), 0, vol))
nrow(Vol)
unique(duplicated(Vol$samp_id))
which(duplicated(Vol$samp_id))
Vol[c(1159, 1160, 1161),]
# Vol = Vol[-c(1159, 1160, 1161),]
Vol[Vol$samp_id=="01015BG000014",]
hist(Vol$vol)
si<-Vol$samp_id[Vol$vol>6000]
sort(unique(dat$age_tot[dat$samp_id%in%si]))
##4 Live Biomass: High stand volume or biomass
Lbio<-dat%>%filter(ld=='L'&dbh>util)%>%group_by(samp_id, meas_yr)%>%summarise(Lbio=sum(wsvha,na.rm = TRUE))%>%
  data.frame()%>%right_join(all_plt, by = c("samp_id",'meas_yr'))%>%
  mutate(Lbio = ifelse(is.na(Lbio), 0, Lbio))
nrow(Lbio)
unique(duplicated(Lbio$samp_id))
which(duplicated(Lbio$samp_id))
Lbio[c(1158,1159,1160),]
Lbio[Lbio$samp_id=="07038CG000006",]
# Lbio = Lbio[-c(1158,1159,1160),]
hist(Lbio$Lbio)
si<-Lbio$samp_id[Lbio$Lbio>6000]
sort(unique(dat$age_tot[dat$samp_id%in%si]))
##5 Dead standing trees: Vol of dead/dying standing trees
DSbio<-dat%>%filter(ld=='D', st_fl=='S'&dbh>=util)%>%group_by(samp_id, meas_yr)%>%summarise(DSbio=sum(wsvha,na.rm = TRUE))%>%
  data.frame()%>%right_join(all_plt, by = c("samp_id",'meas_yr'))%>%mutate(DSbio = ifelse(is.na(DSbio), 0, DSbio))
nrow(DSbio)
unique(duplicated(DSbio$samp_id))
which(duplicated(DSbio$samp_id))
# DSbio = DSbio[-c(869, 874, 957),]
hist(DSbio$DSbio)
DSbio[DSbio$DSbio>600,]
##6 CWD	Large amount/mass of downed CWD: dead fallen biomass
DFbio<-dat%>%filter(ld=='D', st_fl=='F'&dbh>=util)%>%group_by(samp_id, meas_yr)%>%summarise(DFbio=sum(wsvha,na.rm = TRUE))%>%
  data.frame()%>%right_join(all_plt, by = c("samp_id",'meas_yr'))%>%mutate(DFbio = ifelse(is.na(DFbio), 0, DFbio))
unique(duplicated(DFbio$samp_id))
which(duplicated(DFbio$samp_id))
# DFbio = DFbio[-c(467,473, 821),]
nrow(DFbio)

hist(DFbio$DFbio)
DFbio[DFbio$DFbio>400,]
##7 Decay class	Wide decay class distribution of logs and /or snags
names(dat)
dat[dat$samp_id%in%smp,]
DC<-dat%>%filter(dbh>=util)%>%
  group_by(samp_id, meas_yr)%>%summarise(dc_m=mean(pct_dcy, na.rm = TRUE),dc_sd=sd(pct_dcy, na.rm = TRUE))%>%
  data.frame()%>%mutate(DC_CV=dc_sd/dc_m*100)%>%right_join(all_plt, by = c("samp_id",'meas_yr'))

unique(duplicated(DC$samp_id))
which(duplicated(DC$samp_id))
# DC = DC[-c(1159,1160, 1161),]
nrow(DC)

smp<-unique(DC$samp_id[is.na(DC$DC_CV)])
##8 Vertical Complexity: Coeficient of variation of height
Ht_CV<-dat%>%filter(dbh>=util)%>%group_by(samp_id, meas_yr)%>%summarise(ht_m=mean(height, na.rm = TRUE),ht_sd=sd(height, na.rm = TRUE))%>%
  data.frame()%>%mutate(ht_CV=ht_sd/ht_m*100)%>%right_join(all_plt, by = c("samp_id",'meas_yr'))%>%
  mutate(ht_CV = ifelse(is.na(ht_CV), 0, ht_CV))

unique(duplicated(Ht_CV$samp_id))
which(duplicated(Ht_CV$samp_id))
# Ht_CV = Ht_CV[-c(1159,1160, 1161),]
nrow(Ht_CV)
hist(Ht_CV$ht_CV)

##9 Horizontal complexity:	High variation in tree sizes, presence of several cohorts
dbh_CV<-dat%>%filter(dbh>=util)%>%group_by(samp_id, meas_yr)%>%summarise(dbh_m=mean(dbh, na.rm = TRUE),dbh_sd=sd(dbh, na.rm = TRUE))%>%
  data.frame()%>%mutate(dbh_CV=dbh_sd/dbh_m*100)%>%right_join(all_plt, by = c("samp_id",'meas_yr'))%>%
  mutate(dbh_CV = ifelse(is.na(dbh_CV), 0, dbh_CV))

unique(duplicated(dbh_CV$samp_id))
which(duplicated(dbh_CV$samp_id))
# dbh_CV = dbh_CV[-c(1159,1160, 1161),]
nrow(dbh_CV)
hist(dbh_CV$dbh_CV)

##10 Maximum DBH:	Presence of large trees
dbh_max<-dat%>%filter(dbh>=util)%>%group_by(samp_id, meas_yr)%>%summarise(dbh_m=max(dbh, na.rm = TRUE))%>%
  data.frame()%>%right_join(all_plt, by = c("samp_id",'meas_yr'))

unique(duplicated(dbh_max$samp_id))
which(duplicated(dbh_max$samp_id))
# dbh_max = dbh_max[-c(1159,1160, 1161),]
nrow(dbh_max)
hist(dbh_max$dbh_m)

##11 Special attributes	Special attributes (pit and mound relief, presence of epiphytes presence of cavity-trees, tree hollows)	

Sp_at <- dat%>%filter(dbh>=util) %>%mutate(across(
      c(dam_frk, dam_scr, dam_dtop, dam_mis, dam_frs, dam_btop, dam_cnk, dam_rot, dam_bcnk),
      ~ if_else(. >= 1, 1, 0)))%>%group_by(samp_id, meas_yr) %>% summarise(
    att_m = sum(across(c(dam_frk, dam_scr, dam_dtop, dam_mis, dam_frs, dam_btop, dam_cnk, dam_rot, dam_bcnk)), na.rm = TRUE),
    n = n()
  ) %>%data.frame() %>%mutate(SAt_ab = att_m/n*100) %>%right_join(all_plt, by = c("samp_id", "meas_yr"))

unique(duplicated(Sp_at$samp_id))
which(duplicated(Sp_at$samp_id))
# Sp_at = Sp_at[-c(1159,1160, 1161),]
nrow(Sp_at)

hist(Sp_at$SAt_ab)
##12 Stand age	
dat$age_tot[dat$samp_id == "01015BG000014"&dat$meas_yr == 2004]<-58
dat$age_tot[dat$samp_id == "VIL1_1851_CMI"&dat$meas_yr == 2016]<-0
dat$age_tot[dat$samp_id == "05019BR000501"&dat$meas_yr == 1994]<-459


s_age<-dat%>%filter(dbh>=util)%>%group_by(samp_id, meas_yr)%>%summarise(age=max(age_tot, na.rm = TRUE))%>%
  data.frame()%>%right_join(all_plt, by = c("samp_id",'meas_yr'))

nrow(s_age)
unique(duplicated(s_age$samp_id))
which(duplicated(s_age$samp_id))

# s_age = s_age[-c(1159,1160, 1161),]
hist(as.numeric(dat$age_tot))
hist(log(s_age$age))
hist(s_age$age)

##13 Top height	
t_hgt<-dat%>%filter(dbh>=util)%>%group_by(samp_id, meas_yr)%>%summarise(hgt=max(height, na.rm = TRUE))%>%
  data.frame()%>%right_join(all_plt, by = c("samp_id",'meas_yr'))
nrow(t_hgt)
unique(duplicated(t_hgt$samp_id))
which(duplicated(t_hgt$samp_id))
# t_hgt = t_hgt[-c(1159,1160, 1161),]
hist(t_hgt$hgt)

##14 Tree Species Diversity: shannon diversity function
# Biomass per/species for shannon diversity index
sp_bio<- dat%>%filter(ld=='L'&dbh>util)%>%group_by(samp_id, meas_yr,species)%>%summarise(bio=sum(wsvha,na.rm = TRUE))%>%
  data.frame%>%right_join(Lbio, by=c('samp_id', 'meas_yr'))%>%mutate(sp_pi=round(bio/Lbio,3))%>%
  group_by(samp_id, meas_yr)%>%summarise(Sdiv=sum(-(sp_pi*log(sp_pi)), na.rm = TRUE))%>%data.frame()

nrow(sp_bio)
unique(duplicated(sp_bio$samp_id))
hist(sp_bio$Sdiv)

##15 Late succession species for ESSF and SBS BEC zones ('BG','BA',"YC", "HW", "FD", "CW", "HM", 'SS','MB')
#check which species most commonly grow to 250 years or more
LSS<-dat%>%filter(ld=='L'&dbh>100)%>%filter(samp_id%in%unique(s_age$samp_id[s_age$age>500]))%>%group_by(samp_id, meas_yr, species)%>%
  summarise(n=n(), dbh=max(dbh,na.rm=TRUE), hgt=max(height,na.rm=TRUE))%>%
  data.frame()
unique(LSS$species)

LSS<-dat%>%filter(ld=='L'&dbh>util)%>%group_by(samp_id, meas_yr, species)%>%summarise(sp_bio=sum(wsvha,na.rm = TRUE))%>%
  filter(species%in%c("CW", "FD", "HW","BA", "SS","HM","YC"))%>%group_by(samp_id, meas_yr) %>%
  summarise(bio = sum(sp_bio, na.rm = TRUE), .groups = 'drop') %>%   # Sum the biomass across selected species
  data.frame()%>%right_join(Lbio,by=c('samp_id', 'meas_yr'))%>%mutate(bio = ifelse(is.na(bio), 0, bio))%>%
  mutate(LS_sp=(bio/Lbio*100))%>%mutate(LS_sp = ifelse(is.na(LS_sp), 0, LS_sp))
nrow(LSS)
unique(duplicated(LSS$samp_id))

hist(LSS$LS_sp)
###############################################################################
# List of metric data frames
metric_dfs <- list(s_age = s_age, LSS = LSS[c(1,2,5)], sp_bio = sp_bio, t_hgt = t_hgt, 
                   Sp_at = Sp_at[c(1,2,5)],dbh_max = dbh_max, dbh_CV = dbh_CV[c(1,2,5)], Ht_CV = Ht_CV[c(1,2,5)], #DC = DC[c(1,2,5)],
                   VolDF = DFbio, VolDS = DSbio, Vol = Lbio, LTs = LTs[c(1,2,10)], Tden=Tden[c(1,2,5)])

# Initialize data
data <- all_plt[1:2]

# Perform left joins
for (metric_name in names(metric_dfs)) {
  data <- left_join(data, metric_dfs[[metric_name]], by = c('samp_id', 'meas_yr'))
}


names(filtplots)[c(1,3)]<-names(data)[1:2]

data<-left_join(data.frame(data),filtplots,by=c('samp_id','meas_yr'))
nrow(data)
names(data)


# data_sf <- st_as_sf(data, sf_column_name = "geometry", crs = st_crs(fplots))
# nrow(data_sf)

nrow(data)
summary(data)
write_sf(data,"data_complete.shp" )
