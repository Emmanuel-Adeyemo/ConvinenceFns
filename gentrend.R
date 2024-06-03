library(ExSight.R)
library(asreml)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(janitor)
library(LaunchR)

library("qpcR") 

options(launchR.crop="sorghum")

sorg <- read.table("r3_list_from_bill.txt", header = TRUE) %>% clean_names %>% filter(expt_year == 2022)
r3.years <- read.csv("r3_list_from_bill_part2.txt", sep = "\t", header = T) %>% select(-name)

traits_n <- c("YIELD", "MST", "PLTHT", "TSTWT", "LDGSEV", "DAYFLW", "BORLDG", "SCASC", "RLD SC", "COLTOL")
traits_m <- c("DAYFLW")
sorg$expt_year <- as.factor(sorg$expt_year)

#keeps <- c("YIELD~PE", "MST~PE", "PLTHT~PE", "TSTWT~PE", "LDGSEV~PE", "DAYFLW~PE", "BORLDG~PE")
keeps <- c( "DAYFLW~PE")

year_id <- c(2014,2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022 )
year_id_two <- c(2014)

sorg.split <- sorg %>% filter(expt_year %in% year_id)
sorg.split.two <- sorg %>% filter(expt_year %in% year_id_two)

get.blups <- function(x){
  
  sorg.year <- subset(sorg, expt_year == x)  # change data in
  sorg.data <- ExSight.R::data.loader(experiment.ids = sorg.year$expt_id, traits = traits_m, crop = "sorghum", of.interest = "expt")
  blues1 <- ExSight.R::ExSight.R.blues(data.in = sorg.data, traits = traits_m, full.Exsight = TRUE, parallelize = 28, min.aoi.size = 20)
  blues2 <- ExSight.R::post.ExSight.R.blues(blues1)
  
  hyb.blups <- ExSight.R::ExSight.R.blups(blues2, crop = "sorghum", markers = sorghum_markers, model = "zygote", Gmat = TRUE)
  
  blups.year <- process.blups(hyb.blups, x)
  
  blups.res <- id.blups(blups.year, x)
  
}

process.blups <- function(blups.out, x){
  
  blups.data <- blups.out$blups
  blups.names <- intersect(names(blups.data), keeps)
  blups.match <- blups.data %>% select(name, all_of(blups.names))
  blups.match <- blups.match %>% mutate(year = x)
  blups.full <- blups.match %>% rename(dayflw = "DAYFLW~PE")
  #blups.full <- blups.full %>% mutate(yield_kg = yield*67.27273)
  
}


id.blups <- function(processed.blups, x){
  
  x.year <- subset(r3.years, year == x)
  blups.ided <- processed.blups %>% mutate(geNameToID(name)) %>% inner_join(x.year, by = "id")
}


# group one
nahp.grp.one <- lapply(unique(sorg$expt_year), get.blups)
out.names.one <- paste0("sorg.", unique(sorg.split$expt_year))
names(nahp.grp.one) <- out.names.one
out.group.one <- as.data.frame(do.call(rbind, nahp.grp.one))
write.table(out.group.one, "NAHP_out_group_one.txt", sep = "\t")


# group two
nahp.grp.two <- lapply(unique(sorg.split.two$expt_year), get.blups)
out.names.two <- paste0("sorg.", unique(sorg.split.two$expt_year))
names(nahp.grp.two) <- out.names.two
out.group.two <- as.data.frame(do.call(rbind, nahp.grp.two)) %>% rbind(out.group.one) %>% rename(year = "year.x")
write.table(out.group.two, "NAHP_out_group_two.txt", sep = "\t")





r3.years <- searchModel("NAHP R3", crop = "sorghum") %>% mutate(geIDToName(id))
r3.years %>% filter(str_detect(word, "2018"))


write.table(r3.years, "r3_years_nahp.txt", sep = "\t")




out.names <- paste0("sorg.", year_id)

names(test.result) <- out.names

test.result$sorg.2018

yield_bind <- yield_bind %>% rename(yield = "YIELD~PE", mst = "MST~PE", pltht = "PLTHT~PE", ldgsev = "LDGSEV~PE", dayflw = "DAYFLW~PE", borldg = "BORLDG~PE", twt = "TSTWT~PE", scasc = "SCASC~PE", rldsc = "RLDSC~PE", cldtol = "CLDTOL~PE")


s.2018 <- subset(r3.years, year == 2018)
test.result$sorg.2018 %>% mutate(geNameToID(name)) %>% inner_join(s.2018, by = "id")



exptIDToName(3099913 )




out.group.two$year <- as.factor(out.group.two$year)


nahp.dta <- read.csv("NAHP_BLUPS_USE_THIS.txt", header = TRUE, sep = "\t") %>% clean_names
nahp.dta <- nahp.dta %>% mutate(yield_kg = yield*67.27273)

nahp.dta %>% group_by(year) %>% summarise(mean_yield = mean(yield_kg, na.rm = T), sd_yield = sd(yield_kg, na.rm = T), n = sum(!is.na(yield_kg))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean YIELD")+ xlab("Year") + labs(title = "YIELD") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(6000, 8500))

nahp.dta %>% group_by(year) %>% summarise(mean_yield = mean(mst, na.rm = T), sd_yield = sd(mst, na.rm = T), n = sum(!is.na(mst))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean MOISTURE")+ xlab("Year") + labs(title = "MOISTURE") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(12, 16))


nahp.dta %>% group_by(year) %>% summarise(mean_yield = mean(borldg, na.rm = T), sd_yield = sd(borldg, na.rm = T), n = sum(!is.na(borldg))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean BORLDG")+ xlab("Year") + labs(title = "BORLDG") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(-20, 80))


nahp.dta %>% group_by(year) %>% summarise(mean_yield = mean(pltht, na.rm = T), sd_yield = sd(pltht, na.rm = T), n = sum(!is.na(pltht))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean PLTHT")+ xlab("Year") + labs(title = "PLTHT") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(40, 55))


nahp.dta %>% group_by(year) %>% summarise(mean_yield = mean(ldgsev, na.rm = T), sd_yield = sd(ldgsev, na.rm = T), n = sum(!is.na(ldgsev))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean LDGSEV")+ xlab("Year") + labs(title = "LDGSEV") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(3, 9))


nahp.dta %>% group_by(year) %>% summarise(mean_yield = mean(tstwt, na.rm = T), sd_yield = sd(tstwt, na.rm = T), n = sum(!is.na(tstwt))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean TSTWT")+ xlab("Year") + labs(title = "TSTWT") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(52, 62))



nahp.dta %>% group_by(year) %>% summarise(mean_yield = mean(dayflw, na.rm = T), sd_yield = sd(dayflw, na.rm = T), n = sum(!is.na(dayflw))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean DAYFLW")+ xlab("Year") + labs(title = "DAYFLW") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(50, 75))


cols <- c("yield_kg", "mst", "twt", "pltht", "ldgsev", "dayflw", "borldg")
cols <- as.factor(cols)

plot.fig <- function(x){
  
  nahp.dta %>% group_by(year) %>% summarise(mean_yield = mean(x, na.rm = T), sd_yield = sd(x, na.rm = T), n = n()) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
    ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
    geom_line(color = "steelblue3")+
    geom_point()  + ylab(paste("Mean ", as.character(x)))+ xlab("Year") + labs(title = as.character(x)) +
    geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                  position=position_dodge(0.05)) +
    theme_pubr() + grids(linetype = "dashed")+ 
    theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) #+
    #scale_y_continuous(limits = c(6000, 8500))
}

lapply(cols, plot.fig)

lapply(unique(h2.long$loc_group), plot.fig)

nahp.dta %>% group_by(year) %>% summarise(mean_yield = mean(cols[[1]], na.rm = T), sd_yield = sd(cols[1], na.rm = T), n = n()) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " "))


test.sorg <- sorg %>% filter(expt_year == 2022)


#sorg.year <- subset(sorg.split.two, expt_year == x)  # change data in
sorg.data <- ExSight.R::data.loader(experiment.ids = test.sorg$expt_id, traits = "DAYFLW", crop = "sorghum", of.interest = "db")
blues1 <- ExSight.R::ExSight.R.blues(data.in = sorg.data, traits = "DAYFLW", full.Exsight = TRUE, parallelize = 4)
blues2 <- post.ExSight.R.blues(blues1)

hyb.blups <- ExSight.R::ExSight.R.blups(blues2, crop = "sorghum", markers = sorghum_markers, model = "zygote", Gmat = TRUE)
hyb.blups$blups


hhh <- read.csv("NAHP_out_group_two.txt", header = TRUE, sep = "\t") %>% clean_names
s.2020 <- subset(r3.years, year == 2017)
hyb.blups$blups %>% mutate(geNameToID(name)) %>% inner_join(s.2020, by = "id") %>% rename(yld = "YIELD~PE") %>%summarize(meanY = mean(yld))
hhh %>% group_by(year) %>% summarize(meanY = mean(yield, na.rm = T))
hhh <- hhh %>% select(-yield_kg) %>% mutate(yield_kg = yield * 67.27273)



aa <- blues1$DAYFLW$blues.list
ab <- as.data.frame(do.call(rbind, aa))

blup_res = btsat(predicted.value~AOI+GE_ID(r,REML),ab,OtherArgs='WEIGHT=wts')
over_mean <- mean(blup_res$BLUEs$PREDICTED.VALUE)
blup_dt <- blup_res$BLUPs
blup_dt_mean <- blup_dt %>% mutate(colNm = over_mean + PREDICTED.VALUE) %>% select(-PREDICTED.VALUE) %>% rename(id = "Level")

dfw <- r3.years %>% filter(year == 2022) %>% inner_join(blup_dt_mean, by = "id") 

write.table(dfw, "2022_dfw.txt", sep = "\t")


data.table::setnames(blup_dt_mean,'colNm','rand')



r3.model <- searchModel("NAHP R3", crop = "sorghum") %>% mutate(geIDToName(id))

r3_df <- DataFromJob("NAHP R3", crop = "sorghum", searchJobType = "RSTG", searchYear = 2018:2022, traits = "YIELD", What = "BLUE")

DataFromJob("NAHP R2", crop = "sorghum", traits = "yield", What = "BLUE")






a <- setdiff(unique(r3$expt_id), unique(r3.traits$expt_id))
sorg %>% filter(expt_id %in% a)
a <- sorg %>% filter(expt_year == 2022)

z = blueLoader(a$expt_id)


aoi_2022 <- read.csv("aoi_is_for_2022_r3.txt", sep = "\t", header = T)
blues.22 <- deptBlueLoader(aoi_2022$aoi_id)$ge_stats %>% mutate(traitIDToName(trait_id)) %>% rename(trait_nm = "name")
b22.traits <- blues.22 %>% select(-id, -eu_count, -arv, -ge_adjustment)#%>% filter(trait_nm %in% traits_n) 
b22.dta <- b22.traits %>% mutate(expt_year = 2022)


aoi_2021 <- read.csv("aoi_is_for_2022_r3.txt", sep = "\t", header = T) %>% filter(year == 2021)
blues.21 <- deptBlueLoader(aoi_2021$aoi_id)$ge_stats %>% mutate(traitIDToName(trait_id)) %>% rename(trait_nm = "name")
b21.traits <- blues.21 %>% filter(trait_nm %in% traits_n) %>% select(-id, -eu_count, -arv, -ge_adjustment)
b21.dta <- b21.traits %>% mutate(expt_year = 2021)

# Start

sorg <- read.table("r3_list_from_bill.txt", header = TRUE) %>% clean_names
r3.years <- read.csv("r3_list_from_bill_part2.txt", sep = "\t", header = T) %>% select(-name)

traits_n <- c( "DAYFLW")

r3 <- blueLoader(sorg$expt_id) %>% mutate(traitIDToName(trait_id)) %>% rename(trait_nm = "name")
r3.traits <- r3 %>% filter(trait_nm %in% traits_n) %>% select(-id)
r3.dta <- sorg %>% select(expt_id, expt_year) %>% inner_join(r3.traits, by = "expt_id")

blups.btsat <- function(x, WhatYear, dta.year){
  r3.use <- subset(dta.year, trait_nm == x)
  
  blup_res = btsat(trait_value~AOI_ID+GE_ID(r,REML),r3.use,OtherArgs='WEIGHT=weight')
  over_mean <- mean(blup_res$BLUEs$TRAIT_VALUE)
  blup_dt <- blup_res$BLUPs
  blup_dt_mean <- blup_dt %>% mutate(colNm = over_mean + TRAIT_VALUE) %>% select(-TRAIT_VALUE) %>% rename(id = "Level") %>% select(id, colNm)
  data.table::setnames(blup_dt_mean,'colNm',x)
  
  blup.r3.only <- collect.r3.data.only(blup_dt_mean, WhatYear)
  
}

collect.r3.data.only <- function(blups.result, WhatYear){
  
  yoi <- subset(r3.years, year == WhatYear)
  bp.join <- blups.result %>% inner_join(yoi, by = "id")
}


run.blups <- function(WhatYear){
  
  r3.filter <- b22.dta %>% filter(expt_year == WhatYear)  ################## ###################### change data source
  
  out.blups <- lapply(unique(r3.filter$trait_nm), blups.btsat, WhatYear, r3.filter)
  tnames <- unique(r3.filter$trait_nm)
  names(out.blups) <- tnames
  ynames <- c(tnames, "id")
  
  out.now <- out.blups %>% reduce(full_join, by = "id") %>% select(all_of(ynames)) %>% mutate(year = WhatYear)
  return(out.now)
}



blps.12 <- run.blups(2012)
blps.13 <- run.blups(2013)
blps.14 <- run.blups(2014)
blps.15 <- run.blups(2015)
blps.16 <- run.blups(2016)
blps.17 <- run.blups(2017)
blps.18 <- run.blups(2018)
blps.19 <- run.blups(2019)
blps.20 <- run.blups(2020)
blps.21b <- run.blups(2021)
blps.22 <- run.blups(2022)

write.table(a, "tstte_geid.txt", sep = "\t")

blps.21a <- bp.in %>% select(id, colNm) %>% inner_join(blps.21b, by = "id") %>% select(-YIELD) %>% rename(YIELD = "colNm")

blps.21a <- blps.21a %>% select(-PLTHT) %>% rename(PLTHT = "PLTHT_A")
blps.22 <- blps.22 %>% rename(PLTHT = "PLTHT_A")



nahp.blups <- bind_rows(blps.12, blps.13, blps.14, blps.15, blps.16, blps.17, blps.18, blps.19, blps.20, blps.21a, blps.22)
write.table(nahp.blups, "NAHP_BLUPS_USE_THIS.txt", sep = "\t")

####################################################################################################################

blps.21.join <- blps.21 %>% reduce(full_join, by = "id")

r3.long <- sorg %>% select(expt_id, expt_year) %>% inner_join(r3, by = "expt_id")
r3.17 <- r3.dta %>% filter(expt_year == 2020, trait_nm == "PLTHT")#, 2101984))#, 2103180, 2102754, 2101996) ) #%>% select(-name)

blup_res = btsat(trait_value~AOI_ID+GE_ID(r,REML),r3.17,OtherArgs='WEIGHT=weight')
over_mean <- mean(blup_res$BLUEs$TRAIT_VALUE)
blup_dt <- blup_res$BLUPs
blup_dt_mean <- blup_dt %>% mutate(colNm = over_mean + TRAIT_VALUE) %>% select(-TRAIT_VALUE) %>% rename(id = "Level")
data.table::setnames(blup_dt_mean,'colNm','rand')


blps.20$PLTHT <- blup_dt_mean[1:21, "colNm"]

bp.in <- blup_dt_mean %>% inner_join(s, by = "id")
mean(bp.in$colNm, na.rm=T)


ggplot(r3.17, aes(trait_value)) + 
  geom_histogram(aes(y = (..count..)), binwidth = 5) +
  facet_wrap(~aoi_id, ncol = 4) 















##################################################### NACP ##################################

sorg <- read.table("nacp_list_from_bill.txt", header = TRUE) %>% clean_names
r3.years <- read.csv("nacp_list_from_bill_part2.txt", sep = "\t", header = T) %>% select(-name)

traits_n <- c("YIELD", "MST", "PLTHT", "TSTWT", "LDGSEV", "DAYFLW", "SCASC",  "PLTHT_A", "BRLSTL", "BRLCTL", "BRLNLL" )

r3 <- blueLoader(sorg$expt_id) %>% mutate(traitIDToName(trait_id)) %>% rename(trait_nm = "name")
r3.traits <- r3 %>% filter(trait_nm %in% traits_n) %>% select(-id)
r3.dta <- sorg %>% select(expt_id, expt_year) %>% inner_join(r3.traits, by = "expt_id")

blups.btsat <- function(x, WhatYear, dta.year){
  r3.use <- subset(dta.year, trait_nm == x)
  
  blup_res = btsat(trait_value~AOI_ID+GE_ID(r,REML),r3.use,OtherArgs='WEIGHT=weight')
  over_mean <- mean(blup_res$BLUEs$TRAIT_VALUE)
  blup_dt <- blup_res$BLUPs
  blup_dt_mean <- blup_dt %>% mutate(colNm = over_mean + TRAIT_VALUE) %>% select(-TRAIT_VALUE) %>% rename(id = "Level") %>% select(id, colNm)
  data.table::setnames(blup_dt_mean,'colNm',x)
  
  blup.r3.only <- collect.r3.data.only(blup_dt_mean, WhatYear)
  
}

collect.r3.data.only <- function(blups.result, WhatYear){
  
  yoi <- subset(r3.years, year == WhatYear)
  bp.join <- blups.result %>% inner_join(yoi, by = "id")
}


run.blups <- function(WhatYear){
  
  r3.filter <- r3.dta %>% filter(expt_year == WhatYear)  ################## ###################### change data source
  
  out.blups <- lapply(unique(r3.filter$trait_nm), blups.btsat, WhatYear, r3.filter)
  tnames <- unique(r3.filter$trait_nm)
  names(out.blups) <- tnames
  ynames <- c(tnames, "id")
  
  out.now <- out.blups %>% reduce(full_join, by = "id") %>% select(all_of(ynames)) %>% mutate(year = WhatYear)
  return(out.now)
}



blps.12 <- run.blups(2012)
blps.13 <- run.blups(2013)
blps.14 <- run.blups(2014)
blps.15 <- run.blups(2015)
blps.16 <- run.blups(2016)
blps.17 <- run.blups(2017)
blps.18 <- run.blups(2018)
blps.19 <- run.blups(2019)
blps.20 <- run.blups(2020)
blps.21 <- run.blups(2021)
blps.22 <- run.blups(2022)

blps.21a <- bp.in %>% select(id, colNm) %>% inner_join(blps.21b, by = "id") %>% select(-YIELD) %>% rename(YIELD = "colNm")

blps.21a <- blps.21a %>% select(-PLTHT) %>% rename(PLTHT = "PLTHT_A")
blps.22 <- blps.22 %>% rename(PLTHT = "PLTHT_A")



nacp.blups <- bind_rows(blps.12, blps.13, blps.14, blps.15, blps.16, blps.17, blps.18, blps.19, blps.20, blps.22)
write.table(nacp.blups, "NACP_BLUPS_USE_THIS.txt", sep = "\t")


nacp.dta <- nacp.blups %>% clean_names()
nacp.dta <- read.csv("NACP_BLUPS_USE_THIS.txt", header = TRUE, sep = "\t") %>% clean_names
nacp.dta <- nacp.dta %>% mutate(yield_kg = yield*67.27273)

nacp.dta %>% group_by(year) %>% summarise(mean_yield = mean(yield_kg, na.rm = T), sd_yield = sd(yield_kg, na.rm = T), n = sum(!is.na(yield_kg))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean YIELD")+ xlab("Year") + labs(title = "YIELD") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(3500, 7500))

nacp.dta %>% group_by(year) %>% summarise(mean_yield = mean(mst, na.rm = T), sd_yield = sd(mst, na.rm = T), n = sum(!is.na(mst))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean MOISTURE")+ xlab("Year") + labs(title = "MOISTURE") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(12, 16))


nacp.dta %>% group_by(year) %>% summarise(mean_yield = mean(borldg, na.rm = T), sd_yield = sd(borldg, na.rm = T), n = sum(!is.na(borldg))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean BORLDG")+ xlab("Year") + labs(title = "BORLDG") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(-20, 80))


nacp.dta %>% group_by(year) %>% summarise(mean_yield = mean(pltht, na.rm = T), sd_yield = sd(pltht, na.rm = T), n = sum(!is.na(pltht))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean PLTHT")+ xlab("Year") + labs(title = "PLTHT") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(40, 60))


nacp.dta %>% group_by(year) %>% summarise(mean_yield = mean(ldgsev, na.rm = T), sd_yield = sd(ldgsev, na.rm = T), n = sum(!is.na(ldgsev))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean LDGSEV")+ xlab("Year") + labs(title = "LDGSEV") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(3, 9))


nacp.dta %>% group_by(year) %>% summarise(mean_yield = mean(tstwt, na.rm = T), sd_yield = sd(tstwt, na.rm = T), n = sum(!is.na(tstwt))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean TSTWT")+ xlab("Year") + labs(title = "TSTWT") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(52, 62))



nacp.dta %>% group_by(year) %>% summarise(mean_yield = mean(dayflw, na.rm = T), sd_yield = sd(dayflw, na.rm = T), n = sum(!is.na(dayflw))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean DAYFLW")+ xlab("Year") + labs(title = "DAYFLW") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(64, 72))




####################################################################################################################




##################################################### ARG ##################################

sorg <- read.table("arg_list_from_bill.txt", header = TRUE) %>% clean_names
r3.years <- read.csv("arg_list_from_bill_part2.txt", sep = "\t", header = T) %>% select(-name)

traits_n <- c("YIELD", "MST", "PLTHT", "TSTWT", "LDGSEV", "DAYFLW", "SCASC",  "PLTHT_A", "BRLSTL", "BRLCTL", "BRLNLL" )

r3 <- blueLoader(sorg$expt_id) %>% mutate(traitIDToName(trait_id)) %>% rename(trait_nm = "name")
r3.traits <- r3 %>% filter(trait_nm %in% traits_n) %>% select(-id)
r3.dta <- sorg %>% select(expt_id, expt_year) %>% inner_join(r3.traits, by = "expt_id")

blups.btsat <- function(x, WhatYear, dta.year){
  r3.use <- subset(dta.year, trait_nm == x)
  
  blup_res = btsat(trait_value~AOI_ID+GE_ID(r,REML),r3.use,OtherArgs='WEIGHT=weight')
  over_mean <- mean(blup_res$BLUEs$TRAIT_VALUE)
  blup_dt <- blup_res$BLUPs
  blup_dt_mean <- blup_dt %>% mutate(colNm = over_mean + TRAIT_VALUE) %>% select(-TRAIT_VALUE) %>% rename(id = "Level") %>% select(id, colNm)
  data.table::setnames(blup_dt_mean,'colNm',x)
  
  blup.r3.only <- collect.r3.data.only(blup_dt_mean, WhatYear)
  
}

collect.r3.data.only <- function(blups.result, WhatYear){
  
  yoi <- subset(r3.years, year == WhatYear)
  bp.join <- blups.result %>% inner_join(yoi, by = "id")
}


run.blups <- function(WhatYear){
  
  r3.filter <- r3.dta %>% filter(expt_year == WhatYear)  ################## ###################### change data source
  
  out.blups <- lapply(unique(r3.filter$trait_nm), blups.btsat, WhatYear, r3.filter)
  tnames <- unique(r3.filter$trait_nm)
  names(out.blups) <- tnames
  ynames <- c(tnames, "id")
  
  out.now <- out.blups %>% reduce(full_join, by = "id") %>% select(all_of(ynames)) %>% mutate(year = WhatYear)
  return(out.now)
}



blps.12 <- run.blups(2012)
blps.13 <- run.blups(2013)
blps.14 <- run.blups(2014)
blps.15 <- run.blups(2015)
blps.16 <- run.blups(2016)
blps.17 <- run.blups(2017)
blps.18 <- run.blups(2018)
blps.19 <- run.blups(2019)
blps.20 <- run.blups(2020)
blps.21 <- run.blups(2021)
blps.22 <- run.blups(2022)

blps.21a <- bp.in %>% select(id, colNm) %>% inner_join(blps.21b, by = "id") %>% select(-YIELD) %>% rename(YIELD = "colNm")

blps.21a <- blps.21a %>% select(-PLTHT) %>% rename(PLTHT = "PLTHT_A")
blps.22 <- blps.22 %>% rename(PLTHT = "PLTHT_A")



nacp.blups <- bind_rows(blps.12, blps.13, blps.14, blps.15, blps.16, blps.17, blps.18, blps.19, blps.20, blps.22)
write.table(argt, "ARGT_BLUPS_USE_THIS.txt", sep = "\t")


nacp.dta <- nacp.blups %>% clean_names()
argt <- read.csv("ARGT_BLUPS_USE_THIS.txt", header = TRUE, sep = "\t") %>% clean_names
argt <- argt %>% mutate(yield_kg = yield*67.27273)

argt %>% group_by(year) %>% summarise(mean_yield = mean(yield_kg, na.rm = T), sd_yield = sd(yield_kg, na.rm = T), n = sum(!is.na(yield_kg))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean YIELD")+ xlab("Year") + labs(title = "YIELD") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(6000, 7500))

argt %>% group_by(year) %>% summarise(mean_yield = mean(mst, na.rm = T), sd_yield = sd(mst, na.rm = T), n = sum(!is.na(mst))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean MOISTURE")+ xlab("Year") + labs(title = "MOISTURE") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(12, 20))


nacp.dta %>% group_by(year) %>% summarise(mean_yield = mean(borldg, na.rm = T), sd_yield = sd(borldg, na.rm = T), n = sum(!is.na(borldg))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean BORLDG")+ xlab("Year") + labs(title = "BORLDG") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(-20, 80))


argt %>% group_by(year) %>% summarise(mean_yield = mean(pltht, na.rm = T), sd_yield = sd(pltht, na.rm = T), n = sum(!is.na(pltht))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean PLTHT")+ xlab("Year") + labs(title = "PLTHT") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(50, 70))


nacp.dta %>% group_by(year) %>% summarise(mean_yield = mean(ldgsev, na.rm = T), sd_yield = sd(ldgsev, na.rm = T), n = sum(!is.na(ldgsev))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean LDGSEV")+ xlab("Year") + labs(title = "LDGSEV") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(3, 9))


argt %>% group_by(year) %>% summarise(mean_yield = mean(tstwt, na.rm = T), sd_yield = sd(tstwt, na.rm = T), n = sum(!is.na(tstwt))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean TSTWT")+ xlab("Year") + labs(title = "TSTWT") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(52, 62))



argt %>% group_by(year) %>% summarise(mean_yield = mean(dayflw, na.rm = T), sd_yield = sd(dayflw, na.rm = T), n = sum(!is.na(dayflw))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean DAYFLW")+ xlab("Year") + labs(title = "DAYFLW") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(65, 85))



####################################
sorg <- sorg %>% filter(expt_year != 2019)
nahp.grp.one <- lapply(unique(sorg$expt_year), get.blups)
out.names.one <- paste0("sorg.", unique(sorg$expt_year))
names(nahp.grp.one) <- out.names.one
arg.20 <- nahp.grp.one$sorg.2020$blups %>% rename(PLTHT = "PLTHT~PE", DAYFLW = "DAYFLW~PE", MST = "MST~PE", YIELD = "YIELD~PE", id = "ge_id")
arg.21 <- nahp.grp.one$sorg.2021$blups %>% rename(PLTHT = "PLTHT~PE", DAYFLW = "DAYFLW~PE", MST = "MST~PE", YIELD = "YIELD~PE", id = "ge_id")

arg.20$id <- as.numeric(as.character(arg.20$id))
arg.21$id <- as.numeric(as.character(arg.21$id))

out.group.one <- as.data.frame(do.call(rbind, nahp.grp.one))
write.table(out.group.one, "NAHP_out_group_one.txt", sep = "\t")

take <- c("id", "PLTHT", "DAYFLW", "MST", "YIELD")

x.year <- subset(r3.years, year == 2020)
argt.20 <- arg.20 %>% select(all_of(take)) %>% inner_join(x.year, by = "id")

x.year <- subset(r3.years, year == 2021)
argt.21 <- arg.21 %>% select(all_of(take)) %>% inner_join(x.year, by = "id")

argt <- bind_rows(blps.19, argt.20, argt.21) %>% clean_names()


ssc




gd <- c(124552979,
124552981,
132724532,
132724546,
132724547,
132724570,
132724588,
132724597,
132724599,
132724600,
132724627,
132724634,
132724639,
132724640,
132724679,
132724692,
132724716,
132724722,
132724769,
132724569)



geIDToName(gd)


#####################################################################################################################################


##################################################### LACM ##################################

sorg <- read.table("lacm_list_from_bill.txt", header = TRUE) %>% clean_names %>% filter(expt_year == 2022)
r3.years <- read.csv("lacm_list_from_bill_part2.txt", sep = "\t", header = T) %>% select(-name)

traits_n <- c("YIELD", "MST", "PLTHT","DAYFLW")

r3 <- blueLoader(sorg$expt_id) %>% mutate(traitIDToName(trait_id)) %>% rename(trait_nm = "name")
r3.traits <- r3 %>% filter(trait_nm %in% traits_n) %>% select(-id)
r3.dta <- sorg %>% select(expt_id, expt_year) %>% inner_join(r3.traits, by = "expt_id")

blups.btsat <- function(x, WhatYear, dta.year){
  r3.use <- subset(dta.year, trait_nm == x)
  
  blup_res = btsat(trait_value~AOI_ID+GE_ID(r,REML),r3.use,OtherArgs='WEIGHT=weight')
  over_mean <- mean(blup_res$BLUEs$TRAIT_VALUE)
  blup_dt <- blup_res$BLUPs
  blup_dt_mean <- blup_dt %>% mutate(colNm = over_mean + TRAIT_VALUE) %>% select(-TRAIT_VALUE) %>% rename(id = "Level") %>% select(id, colNm)
  data.table::setnames(blup_dt_mean,'colNm',x)
  
  blup.r3.only <- collect.r3.data.only(blup_dt_mean, WhatYear)
  
}

collect.r3.data.only <- function(blups.result, WhatYear){
  
  yoi <- subset(r3.years, year == WhatYear)
  bp.join <- blups.result %>% inner_join(yoi, by = "id")
}


run.blups <- function(WhatYear){
  
  r3.filter <- r3.dta %>% filter(expt_year == WhatYear)  ################## ###################### change data source
  
  out.blups <- lapply(unique(r3.filter$trait_nm), blups.btsat, WhatYear, r3.filter)
  tnames <- unique(r3.filter$trait_nm)
  names(out.blups) <- tnames
  ynames <- c(tnames, "id")
  
  out.now <- out.blups %>% reduce(full_join, by = "id") %>% select(all_of(ynames)) %>% mutate(year = WhatYear)
  return(out.now)
}



blps.12 <- run.blups(2012)
blps.13 <- run.blups(2013)
blps.14 <- run.blups(2014)
blps.15 <- run.blups(2015)
blps.16 <- run.blups(2016)
blps.17 <- run.blups(2017)
blps.18 <- run.blups(2018)
blps.19 <- run.blups(2019)
blps.20 <- run.blups(2020)
blps.21 <- run.blups(2021)
blps.22 <- run.blups(2022)

blps.12$YIELD <- blps.12$YIELD - 42
blps.14$YIELD <- blps.14$YIELD - 46
blps.16$YIELD <- blps.16$YIELD -20

blps.21a <- bp.in %>% select(id, colNm) %>% inner_join(blps.21b, by = "id") %>% select(-YIELD) %>% rename(YIELD = "colNm")

blps.21a <- blps.21a %>% select(-PLTHT) %>% rename(PLTHT = "PLTHT_A")
blps.22 <- blps.22 %>% rename(PLTHT = "PLTHT_A")



lacm.blups <- bind_rows(blps.12, blps.13, blps.14, blps.15, blps.16, blps.17, blps.18, blps.19, blps.20)#, blps.22)
write.table(lacm.blups, "LACM_BLUPS_USE_THIS.txt", sep = "\t")
lacm <- lacm.blups %>% clean_names


lamc.dta <- lacm.blups %>% clean_names()
lacm <- read.csv("LACM_BLUPS_USE_THIS.txt", header = TRUE, sep = "\t") %>% clean_names
lacm <- lacm %>% mutate(yield_kg = yield*67.27273)

lacm %>% group_by(year) %>% summarise(mean_yield = mean(yield_kg, na.rm = T), sd_yield = sd(yield_kg, na.rm = T), n = sum(!is.na(yield_kg))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean YIELD")+ xlab("Year") + labs(title = "YIELD") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(9000, 14000))

lacm %>% group_by(year) %>% summarise(mean_yield = mean(mst, na.rm = T), sd_yield = sd(mst, na.rm = T), n = sum(!is.na(mst))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean MOISTURE")+ xlab("Year") + labs(title = "MOISTURE") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(12, 20))


nacp.dta %>% group_by(year) %>% summarise(mean_yield = mean(borldg, na.rm = T), sd_yield = sd(borldg, na.rm = T), n = sum(!is.na(borldg))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean BORLDG")+ xlab("Year") + labs(title = "BORLDG") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(-20, 80))


lacm %>% group_by(year) %>% summarise(mean_yield = mean(pltht, na.rm = T), sd_yield = sd(pltht, na.rm = T), n = sum(!is.na(pltht))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean PLTHT")+ xlab("Year") + labs(title = "PLTHT") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(65, 90))


nacp.dta %>% group_by(year) %>% summarise(mean_yield = mean(ldgsev, na.rm = T), sd_yield = sd(ldgsev, na.rm = T), n = sum(!is.na(ldgsev))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean LDGSEV")+ xlab("Year") + labs(title = "LDGSEV") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(3, 9))


argt %>% group_by(year) %>% summarise(mean_yield = mean(tstwt, na.rm = T), sd_yield = sd(tstwt, na.rm = T), n = sum(!is.na(tstwt))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean TSTWT")+ xlab("Year") + labs(title = "TSTWT") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(52, 62))



lacm %>% group_by(year) %>% summarise(mean_yield = mean(dayflw, na.rm = T), sd_yield = sd(dayflw, na.rm = T), n = sum(!is.na(dayflw))) %>% mutate(sd = sd_yield*(1.96/sqrt(n))) %>% mutate(n = paste0("(", n, ")")) %>% mutate(year_n = paste(year, n, sep = " ")) %>%
  ggplot(aes(x=year_n, y=mean_yield, group = 1)) + 
  geom_line(color = "steelblue3", size = 0.7)+
  geom_point()  + ylab("Mean DAYFLW")+ xlab("Year") + labs(title = "DAYFLW") +
  geom_errorbar(aes(ymin=mean_yield-sd, ymax=mean_yield+sd), color = "black", width=.2,
                position=position_dodge(0.05)) +
  theme_pubr() + grids(linetype = "dashed")+ 
  theme(panel.border = element_rect(color = "blue4",fill = NA,size = 1)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(70, 90))












aoi_id_lacm <- read.csv("lacm_aoi.txt", sep = "\t", header = T) %>% filter(year == 2022)

a <- sorg %>% filter(expt_year == 2021)
s <- deptBlueLoader(aoi_id_lacm$aoi_id)$ge_stats

lacm21 <- r3.years %>% filter(year == 2021) %>% inner_join(lacm21, by = "id")



deptBlueLoader(2173218)


s13 <- blueLoader(3248250) %>% filter(trait_id == 11)

blup_res = btsat(trait_value~AOI_ID+GE_ID(r,REML),s13,TraitArgs='WEIGHT=weight')
over_mean <- mean(blup_res$BLUEs$TRAIT_VALUE)
blup_dt <- blup_res$BLUPs
blup_dt_mean <- blup_dt %>% mutate(colNm = over_mean + TRAIT_VALUE) %>% select(-TRAIT_VALUE) %>% rename(id = "Level")

mst <- blup_dt_mean %>% select(colNm) %>% rename(mst = "colNm")

write.table(blup_dt_mean, "lacm_yield_13.txt", sep = "\t")
























aoi_id_nacp <- read.csv("aoi_is_for_2022_nacp.txt", sep = "\t", header = T)

a <- sorg %>% filter(expt_year == 2022)
s <- deptBlueLoader(aoi_id_nacp$aoi_id)$ge_stats

dta.20 <- s %>% mutate(traitIDToName(trait_id)) %>% rename(trait_nm = "name")
dta.20 <- dta.20 %>% filter(trait_nm %in% traits_n) %>% select(-id) %>% mutate(expt_year = 2022)

sorg.yearnow <- sorg.split %>% filter(expt_year == 2019)
sorg.data <- ExSight.R::data.loader(experiment.ids = 4352745, traits = traits_n, crop = "sorghum", of.interest = "db")
blues1 <- ExSight.R::ExSight.R.blues(data.in = sorg.data, traits = traits_n, parallelize = 4, min.aoi.size = 10)
blues2 <- post.ExSight.R.blues(blues1)
hyb.blups <- ExSight.R::ExSight.R.blups(blues2, crop = "sorghum", markers = sorghum_markers, model = "zygote", Gmat = TRUE)
blups.data <- hyb.blups$blups
blups.names <- intersect(names(blups.data), keeps)
blups.match <- blups.data %>% select(name, all_of(blups.names))
blups.match <- blups.match %>% mutate(year = 2019)
blups.full <- blups.match %>% rename(ldgsev = "LDGSEV~PE", dayflw = "DAYFLW~PE", twt = "TSTWT~PE")
x.year <- subset(r3.years, year == 2019)
blups.ided <- blups.full %>% mutate(geNameToID(name)) %>% inner_join(x.year, by = "id")
y19 <- blups.ided
nacp.others <- bind_rows(y14, y15, y16, y17, y18, y19)
write.table(nacp.others, "nacp.others.txt", sep = "\t")



a <- blues1$DAYFLW$blues.list %>% reduce(rbind)
blup_res = btsat(predicted.value~AOI+GE_ID(r,REML),a,TraitArgs='WEIGHT=wts')
over_mean <- mean(blup_res$BLUEs$PREDICTED.VALUE)
blup_dt <- blup_res$BLUPs
blup_dt_mean <- blup_dt %>% mutate(colNm = over_mean + PREDICTED.VALUE) %>% select(-PREDICTED.VALUE) %>% rename(id = "Level")

dfw <- blup_dt_mean %>% select(colNm, id) %>% rename(dfw = "colNm")


lacm21 <- bind_cols(yld, mst, dfw, plt) %>% mutate(year = 2021)
lacm21 <- lacm21 %>% mutate(year = 2021)
lacm_one <- bind_rows(lacm22, lacm21)
write.table(lacm_one, "lacm_two.txt", sep = "\t")




# group one
nacp.grp.one <- lapply(unique(sorg.split$expt_year), get.blups)

out.names.one <- paste0("sorg.", unique(sorg.split$expt_year))
names(nahp.grp.one) <- out.names.one
out.group.one <- as.data.frame(do.call(rbind, nahp.grp.one))


r3.long <- sorg %>% select(expt_id, expt_year) %>% inner_join(r3, by = "expt_id")
r3.17 <- r3.dta %>% filter(expt_year == 2022)#, trait_nm == "PLTHT")#, 2101984))#, 2103180, 2102754, 2101996) ) #%>% select(-name)

blup_res = btsat(trait_value~AOI_ID+GE_ID(r,REML),r3.17,OtherArgs='WEIGHT=weight')
over_mean <- mean(blup_res$BLUEs$TRAIT_VALUE)
blup_dt <- blup_res$BLUPs
blup_dt_mean <- blup_dt %>% mutate(colNm = over_mean + TRAIT_VALUE) %>% select(-TRAIT_VALUE) %>% rename(id = "Level")


a =DataFromStage(STAGES=c('R3'), GEO='NA',searchYear= 2021, crop='sorghum',traits='DAYFLW',RawObs=FALSE)


samcount <- read.csv("samcount.txt", sep = "\t", header = T)
samcount <- samcount %>% mutate(newC = ifelse(count < 6, count, sam))
write.table(samcount, "samcount_new.txt", sep = "\t")
