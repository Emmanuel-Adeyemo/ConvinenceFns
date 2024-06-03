library(dplyr)
library(reshape2)
library(bWGR)
library(btsatr)
library(LaunchR)
library(stringr)


########################################

get.datr <- function(GEO, searchYear, prdYear){
  
  dta = DataFromStage(STAGES = c("R1", "R2"), GEO=GEO,crop='sorghum',traits="YIELD", searchYear = searchYear)
  #dta <- dta %>% filter(geo == GEO)
  dta = dta %>% filter(!str_detect(maturity, "Forage")) # remove forage materials
  
  keep <- c('ge_id','aoi_id','trait_name','trait_value')
  dta.tn <- dta %>% filter(year != prdYear) %>% select(all_of(keep))
  PSet <- dta %>% filter(year == prdYear) %>% select(ge_id) %>% distinct(ge_id)
  
  out <- list(dta = dta, dta.tn = dta.tn, PSet = PSet)
  return(out)
}

GEO = "NAHP"

dta <- get.datr(GEO, 2022:2023, 2023)
dta = dta$dta 
out = c("AFSA", "MMSD")
`%!in%` = Negate(`%in%`)
dta = dta %>% filter(geo %!in% out )
dta = dta %>% rename(yield = "trait_value")

keep = c("aoi_id", "ge_id", "AoiCode", "yield")

dta = dta %>% select(all_of(keep))

id <- as.numeric(dta$ge_id)

Z <- GetZ(id, crop = "sorghum")

#######################################


r1_join <- readRDS("r1_pheno.RData")
Z <- readRDS("Z_r1.RData")


#r1_join = dta

r1_join %>% head
r1_join <- r1_join %>% filter(aoiname != "MMSB02LMARQN", year %in% c(2020,2021))

r1_join$ Loc = paste(r1_join$aoiname,r1_join$year)

Y = dcast(r1_join, ge_id~Loc, mean, value.var='yield' )

Y = dcast(dta, ge_id~AoiCode, mean, value.var='yield' )

rownames(Y) = Y$ge_id
Y = data.matrix(Y[,-1])
Y = Y[,which(colSums(!is.na(Y))>150)]


ii = intersect(rownames(Y),rownames(Z))
Y = Y[ii,]
Z = Z[ii,]

v = apply(Z,2,var, na.rm = TRUE)
Z = Z[,which(v>0.1)]
Z = Z-1
Z[is.na(Z)] = 0
Z %>% dim



#########
E1 = EigenGAU(Z,phi=1.2) # lowered from 1.5 to 1.2 on 07/05/23
E = EigenEVD(E1)
E$U = E$U[,which(E$D>0.01)]
E$D = E$D[which(E$D>0.01)]
Xgen = E$U %*% diag(sqrt(E$D))
rownames(Xgen) = rownames(Z); 
rm(E)
rm(E1)




Z = IMP(Z)
boxplot(Y)


fit = MRR3(Y,Z,cores=8)


GC = fit$GC
rownames(GC) = colnames(GC) = colnames(Y)
disT = dist(1-GC) %>% as.matrix()
#disb = dist(GC) %>% as.matrix()
test = 1-GC
test_dist = test %>% dist() %>% as.matrix()
hc = hclust(dist(1-GC))
plot(hc)


write.table(disT, "Nacp_nahp_distance_matrix.txt", sep="\t")





####################################

traits = c("YIELD")

dta = r1_join

if(length(traits)>0) dta = subset(dta,(trait_name %in% traits))
cat('Converting tall to wide data\n')
# replace reshape2 on 05/02/23

keep = c("gename", "trait_value", "AoiCode", "trait_name")

dta = dta %>% select(all_of(keep))
Ys = suppressWarnings(by(dta, dta$trait_name, function(dta) reshape(dta, direction = 'wide', v.names = 'trait_value', timevar = 'AoiCode',  idvar='gename')[,-2] )) 
Ys = Ys[which(sapply(Ys,ncol)>=3)]  # debug with Rocio 12/13/22
Ys = lapply(Ys, function(x){colnames(x)=gsub('trait_value\\.','',colnames(x));return(x)} )
tmp_trts = names(Ys) # debug with Rocio 12/13/22
Ys = lapply(tmp_trts, function(whichY){ 
  Y = Ys[[whichY]]
  GEnames = Y$gename; AoiName = (colnames(Y)[-1]);
  rownames(Y) = Y$gename; Y = data.matrix(Y[,-1]);
  #ww = which(colMeans(!is.na(Y))>0.5)
  ww = which(colSums(!is.na(Y))>10)
  if(length(ww)>1){
    Y = Y[,ww]; }else if(length(ww)==1){
      Y = matrix(Y[,ww],ncol = 1,dimnames = list(GEnames,AoiName[ww]))  # debug with Rocio 12/13/22
    }else{ Y = NULL}
  return(Y)})
names(Ys) = tmp_trts # debug with Rocio 12/13/22
Ys = Ys[which(!sapply(Ys,is.null))] # debug with Ron 08/22/22
Ys = Ys[which(sapply(Ys,nrow)>10)]




cat(paste0('Fit models (',length(Ys),' traits)\n'))
fits = suppressMessages(
  lapply(Ys, function(Y){
    if(UseGenomics|GCA){
      fit = MRR3(Y,Xgen[rownames(Y),],df0=5,maxit=99,XFA2=XFA2) }
    else{
        fit = MRR3(Y,diag(nrow(Y)),df0=5,maxit=99,XFA2=XFA2) }}))
    

res = lapply(Ys, function(Y){
fit = MRR3(Y, diag(nrow(Y)),df0=5,maxit=99)})


aoi_names = Ys[[1]]
test = res$YIELD$mu %>% data.frame() %>% round(0)

test = res$YIELD$hat

colnames(test) = colnames(aoi_names)
rownames(test) = rownames(aoi_names)
test = as.data.frame(test)

test = test %>% tibble::rownames_to_column("ge_id")

test_res = test %>% tidyr::pivot_longer(!ge_id, names_to = "AoiCode", values_to = "pred_yield")



pheno = aoi_names %>% as.data.frame %>% tibble::rownames_to_column("ge_id") %>% tidyr::pivot_longer(!ge_id, names_to = "AoiCode", values_to = "yield")

join_pheno = merge(pheno, test_res, by = c("ge_id", "AoiCode"))

join_pheno %>% group_by(AoiCode) %>% summarise(corr = cor(pred_yield, yield, use = "p"))
####################

boxplot(test[-1])





####


library(pacman)
p_load(akima,corrgram,gclus,gge,latticeExtra,luzlogr,pagenum,pls,readr,readxl,reshape2)

library(gex)

gex()


library(readr)

dta = read_csv("NACP GDCS YIELD BLUPS 2018-2023.csv")


dta_new = dta %>% tidyr::pivot_longer(cols = "2070A992-01":XS8523, names_to = "ge_id", values_to = "yield")

write_csv(dta_new, "NACP_Heber_Stuff.csv")


#














































require(ggplot2)

require(ggpubr)
require(ggfortify)
require(factoextra)
require(NbClust)
require(ggrepel)

dist.gc <- as.matrix(dist(1-GC))

#' @aliases get optimum K for clustering - elbow method
get.optimum.k <- function(geno_distance, intercept_line){
  
  optimum_k <- fviz_nbclust(geno_distance, kmeans, method = "wss") +
    geom_vline(xintercept = intercept_line, linetype = 2)+  # include intercept number after optimum number has been determined
    labs(subtitle = "Elbow method")
  
  return(optimum_k)
}

get.optimum.k(dist.gc, 2)

#' @aliases get Kmeans
get.kmeans <- function(geno_distance, number_of_clusters){
  
  kmeans_out <- kmeans(geno_distance, number_of_clusters, nstart=2, iter.max=1000, algorithm="Hartigan-Wong") # options:algorith="Hartigan-Wong", "Lloyd", "Forgy"
  
  return(kmeans_out)
}

ans.cl <- get.kmeans(dist.gc, 2)


library(FactoMineR)
pca.res <- PCA(GC)
pca <-pca.res$ind$coord
pca <- as.data.frame(pca)
plot(pca$Dim.1, pca$Dim.2)

cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)

#' @aliases plot kmeans
plot.kmeans <- function(cluster_with_pca, cluster_number){
  
  cluster_with_pca[,cluster_number] <- as.factor(cluster_with_pca[,cluster_number])
  
  kmeans_plot <- ggplot(cluster_with_pca, aes(x=Dim.1, y=Dim.2, shape=cluster_with_pca[,cluster_number], color=cluster_with_pca[,cluster_number])) +
    geom_point()
  
  kmeans_bw <- kmeans_plot + theme_bw()
  
  return(kmeans_bw)
}



plot.kmeans(cl.pca, "cl.ass")


fviz_cluster(ans.cl, data = GC,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#E7B900"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

eigv <- pca.res$eig[,1:2]
var.pct <- round(eigv[,2], 1)

require(ggpubr)


cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  #shape = "Species", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  stat_mean(aes(color = cl.ass), size = 4)

rr <- gdcs %>% group_by(loc) %>% filter(n() > 21)


mu <- gdcs %>% group_by(loc) %>% filter(n() > 21) %>% summarize(MeanY = mean(yield))

summary(mu$MeanY)

mu <- mu %>% mutate(class = ifelse(MeanY >= 130, "High", 
                                   ifelse(MeanY < 130 & MeanY >= 115, "Mid-High", 
                                          ifelse(MeanY < 115 & MeanY >= 100, "Mid-Low",
                                                 ifelse(MeanY < 100, "Low", "Unk")))))

mu <- mu %>% mutate(class = ifelse(MeanY >= 115, "High", "Low"))

cl.pca$class <- mu$class

options(ggrepel.max.overlaps = Inf)

ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "class", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,#subset(cl.pca, class == "High"),
    #aes(label = rownames(subset(cl.pca, class == "Low"))),
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )



ver_all <- cl.pca %>% group_by(class) %>% summarize(n_total = n()) %>% 
  mutate(freq = formattable::percent(n_total / sum(n_total))) %>% 
  mutate(ypos = cumsum(n_total)- 0.5*n_total ) %>%
  as.data.frame()
# Make the plot
ggplot(ver_all, aes(x=2, y=n_total, fill=class)) +
  geom_bar(position = 'fill', stat = 'identity')  +
  xlim(0.5, 2.5) +
  coord_polar(theta = 'y') + 
  labs(x=NULL, y=NULL)+
  theme_void()

ver_all$class <- factor(ver_all$class, levels = c("High", "Mid-High", "Mid-Low", "Low"))
ver_all %>%
  ggplot( aes(x=class, y=n_total)) +
  geom_segment( aes(xend=class, yend=0)) +
  geom_point( size=4, color="orange") +
  coord_flip() +
  theme_bw() +
  xlab("")



# three clusters
ans.cl <- get.kmeans(dist.gc, 3)
cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)
mu <- mu %>% mutate(class = ifelse(MeanY >= 130, "High", 
                                   ifelse(MeanY < 130 & MeanY >= 100, "Mid", 
                                       ifelse(MeanY < 100, "Low", "Unk"))))
cl.pca$class <- mu$class
cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
cl.pca$class <- factor(cl.pca$class, levels = c("Low", "Mid", "High"))

ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "class", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
)# +
  

ver_all <- cl.pca %>% group_by(class) %>% summarize(n_total = n()) %>% 
  mutate(freq = formattable::percent(n_total / sum(n_total))) %>% 
  mutate(ypos = cumsum(n_total)- 0.5*n_total ) %>%
  as.data.frame()
ver_all$class <- factor(ver_all$class, levels = c("High", "Mid", "Low"))
ver_all %>%
  ggplot( aes(x=class, y=n_total)) +
  geom_segment( aes(xend=class, yend=0)) +
  geom_point( size=4, color="orange") +
  coord_flip() +
  theme_bw() +
  xlab("")




# two clusters
ans.cl <- get.kmeans(dist.gc, 2)
cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)
mu <- r1_join %>% group_by(Loc) %>% filter(n() > 100) %>% summarize(MeanY = mean(yield))
mu <- mu %>% mutate(class = ifelse(MeanY >= 100, "High", "Low"))

cl.pca$class <- mu$class
cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
#cl.pca$class <- factor(cl.pca$class, levels = c("Low", "High"))

ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "class", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )



ver_all <- cl.pca %>% group_by(class) %>% summarize(n_total = n()) %>% 
  mutate(freq = formattable::percent(n_total / sum(n_total))) %>% 
  mutate(ypos = cumsum(n_total)- 0.5*n_total ) %>%
  as.data.frame()
ver_all$class <- factor(ver_all$class, levels = c("High", "Low"))
ver_all %>%
  ggplot( aes(x=class, y=n_total)) +
  geom_segment( aes(xend=class, yend=0)) +
  geom_point( size=4, color="orange") +
  coord_flip() +
  theme_bw() +
  xlab("")


#############################################################################

r1_join <- readRDS("r4_nacp.RData")
Z <- readRDS("Z_r4_nacp.RData")

Y = dcast(r1_join, ge_id~aoi_code, mean, value.var='yield' )
rownames(Y) = Y$ge_id
Y = data.matrix(Y[,-1])
Y = Y[,which(colSums(!is.na(Y))>30)]


ii = intersect(rownames(Y),rownames(Z))
Y = Y[ii,]
Z = Z[ii,]

v = apply(Z,2,var, na.rm = TRUE)
Z = Z[,which(v>0.1)]
Z = Z-1
Z[is.na(Z)] = 0
Z %>% dim



fit = MRR3(Y,Z,cores=4)


GC = fit$GC
rownames(GC) = colnames(GC) = colnames(Y)
hc = hclust(as.dist(1-GC))
plot(hc)


dist.gc <- as.matrix(dist(1-GC))

get.optimum.k(dist.gc, 4)

ans.cl <- get.kmeans(dist.gc, 3)


pca.res <- PCA(GC)
pca <-pca.res$ind$coord
pca <- as.data.frame(pca)

cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)

eigv <- pca.res$eig[,1:2]
var.pct <- round(eigv[,2], 1)

cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "cl.ass", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

mu <- r1_join %>% group_by(aoi_code) %>% filter(n() > 30) %>% summarize(MeanY = mean(yield))

summary(mu$MeanY)

mu <- mu %>% mutate(class = ifelse(MeanY >= 110, "High",
                                   ifelse(MeanY >= 80 & MeanY < 110, "Mid", "Low")))

rem <- setdiff(mu$aoi_code, rownames(cl.pca))
`%nin%` = Negate(`%in%`)
mu <- mu %>% filter(aoi_code %nin% rem)
all.equal(mu$aoi_code, rownames(cl.pca))


cl.pca$class <- mu$class
cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
#cl.pca$class <- factor(cl.pca$class, levels = c("Low", "High"))

ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "class", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )









#################################### GDCS #########################################
 
r1_join <- readRDS("cp_gdcs.RData")
Z <- readRDS("Z_cp_gdcs.RData")

Y = dcast(r1_join, ge_id~aoi_code, mean, value.var='yield' )
rownames(Y) = Y$ge_id
Y = data.matrix(Y[,-1])
Y = Y[,which(colSums(!is.na(Y))>30)]


ii = intersect(rownames(Y),rownames(Z))
Y = Y[ii,]
Z = Z[ii,]

v = apply(Z,2,var, na.rm = TRUE)
Z = Z[,which(v>0.1)]
Z = Z-1
Z[is.na(Z)] = 0
Z %>% dim



fit = MRR3(Y,Z,cores=4)


GC = fit$GC
rownames(GC) = colnames(GC) = colnames(Y)
hc = hclust(as.dist(1-GC))
plot(hc)


dist.gc <- as.matrix(dist(1-GC))

get.optimum.k(dist.gc, 4)

ans.cl <- get.kmeans(dist.gc, 3)


pca.res <- PCA(GC)
pca <-pca.res$ind$coord
pca <- as.data.frame(pca)

cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)

eigv <- pca.res$eig[,1:2]
var.pct <- round(eigv[,2], 1)

cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "cl.ass", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

mu <- r1_join %>% group_by(aoi_code)  %>% summarize(MeanY = mean(yield))

summary(mu$MeanY)

mu <- mu %>% mutate(class = ifelse(MeanY >= 95, "High",
                                   ifelse(MeanY >= 80 & MeanY < 95, "Mid", "Low")))

rem <- setdiff(mu$aoi_code, rownames(cl.pca))
`%nin%` = Negate(`%in%`)
mu <- mu %>% filter(aoi_code %nin% rem)
all.equal(mu$aoi_code, rownames(cl.pca))


cl.pca$class <- mu$class
cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
#cl.pca$class <- factor(cl.pca$class, levels = c("Low", "High"))

ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "class", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )




# two clusters
ans.cl <- get.kmeans(dist.gc, 2)


pca.res <- PCA(GC)
pca <-pca.res$ind$coord
pca <- as.data.frame(pca)

cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)

eigv <- pca.res$eig[,1:2]
var.pct <- round(eigv[,2], 1)

cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "cl.ass", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

mu <- r1_join %>% group_by(aoi_code)  %>% summarize(MeanY = mean(yield))

summary(mu$MeanY)

mu <- mu %>% mutate(class = ifelse(MeanY >= 95, "High",
                                   ifelse(MeanY >= 80 & MeanY < 95, "Mid", "Low")))

rem <- setdiff(mu$aoi_code, rownames(cl.pca))
`%nin%` = Negate(`%in%`)
mu <- mu %>% filter(aoi_code %nin% rem)
all.equal(mu$aoi_code, rownames(cl.pca))


cl.pca$class <- mu$class
cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
#cl.pca$class <- factor(cl.pca$class, levels = c("Low", "High"))

ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "class", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )






#################################### NAHP GDCS #########################################

r1_join <- readRDS("hp_gdcs.RData")
Z <- readRDS("Z_hp_gdcs.RData")

remove <- c('MMSY04NCONCN_2020','MMSY04NSOLON_2020','PPSY04FBO01N_2020','PPSY04FSL01N_2020', "MMSB02LMARQN_2020", "MMSB02LMARQN_2021")
`%nin%` = Negate(`%in%`)
r1_join <- r1_join %>% filter(aoi_code %nin% remove)

Y = dcast(r1_join, ge_id~aoi_code, mean, value.var='yield' )
rownames(Y) = Y$ge_id
Y = data.matrix(Y[,-1])
Y = Y[,which(colSums(!is.na(Y))>15)]


ii = intersect(rownames(Y),rownames(Z))
Y = Y[ii,]
Z = Z[ii,]

v = apply(Z,2,var, na.rm = TRUE)
Z = Z[,which(v>0.1)]
Z = Z-1
Z[is.na(Z)] = 0
Z %>% dim



fit = MRR3(Y,Z,cores=4)


GC = fit$GC
rownames(GC) = colnames(GC) = colnames(Y)
hc = hclust(as.dist(1-GC))
plot(hc)


dist.gc <- as.matrix(dist(1-GC))

get.optimum.k(dist.gc, 3)

ans.cl <- get.kmeans(dist.gc, 3)


pca.res <- PCA(GC)
pca <-pca.res$ind$coord
pca <- as.data.frame(pca)

cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)

eigv <- pca.res$eig[,1:2]
var.pct <- round(eigv[,2], 1)

cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "cl.ass", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

mu <- r1_join %>% group_by(aoi_code)  %>% summarize(MeanY = mean(yield))

summary(mu$MeanY)

mu <- mu %>% mutate(class = ifelse(MeanY >= 110, "High", "Low"))

mu <- mu %>% mutate(class = ifelse(MeanY >= 120, "High",
                                   ifelse(MeanY >= 90 & MeanY < 125, "Mid", "Low")))

rem <- setdiff(mu$aoi_code, rownames(cl.pca))
`%nin%` = Negate(`%in%`)
mu <- mu %>% filter(aoi_code %nin% rem)
all.equal(mu$aoi_code, rownames(cl.pca))


cl.pca$class <- mu$class
cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
#cl.pca$class <- factor(cl.pca$class, levels = c("Low", "High"))

ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "class", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )




# two clusters
ans.cl <- get.kmeans(dist.gc, 2)


pca.res <- PCA(GC)
pca <-pca.res$ind$coord
pca <- as.data.frame(pca)

cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)

eigv <- pca.res$eig[,1:2]
var.pct <- round(eigv[,2], 1)

cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "cl.ass", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

mu <- r1_join %>% group_by(aoi_code)  %>% summarize(MeanY = mean(yield))

summary(mu$MeanY)

mu <- mu %>% mutate(class = ifelse(MeanY >= 110, "High", "Low"))

mu <- mu %>% mutate(class = ifelse(MeanY >= 95, "High",
                                   ifelse(MeanY >= 80 & MeanY < 95, "Mid", "Low")))

rem <- setdiff(mu$aoi_code, rownames(cl.pca))
`%nin%` = Negate(`%in%`)
mu <- mu %>% filter(aoi_code %nin% rem)
all.equal(mu$aoi_code, rownames(cl.pca))


cl.pca$class <- mu$class
cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
#cl.pca$class <- factor(cl.pca$class, levels = c("Low", "High"))

ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "class", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )



#################################### LACM GDCS #########################################

r1_join <- readRDS("la_r12.RData")
Z <- readRDS("Z_la_r12.RData")

Y = dcast(r1_join, ge_id~aoi_code, mean, value.var='yield' )
rownames(Y) = Y$ge_id
Y = data.matrix(Y[,-1])
#Y = Y[,which(colSums(!is.na(Y))>12)]


ii = intersect(rownames(Y),rownames(Z))
Y = Y[ii,]
Z = Z[ii,]

v = apply(Z,2,var, na.rm = TRUE)
Z = Z[,which(v>0.1)]
Z = Z-1
Z[is.na(Z)] = 0
Z %>% dim



fit = MRR3(Y,Z,cores=4)


GC = fit$GC
rownames(GC) = colnames(GC) = colnames(Y)
hc = hclust(as.dist(1-GC))
plot(hc)


dist.gc <- as.matrix(dist(1-GC))

get.optimum.k(dist.gc, 3)

ans.cl <- get.kmeans(dist.gc, 3)


pca.res <- PCA(GC)
pca <-pca.res$ind$coord
pca <- as.data.frame(pca)

cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)

eigv <- pca.res$eig[,1:2]
var.pct <- round(eigv[,2], 1)

cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "cl.ass", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

mu <- r1_join %>% group_by(aoi_code)  %>% summarize(MeanY = mean(yield))

summary(mu$MeanY)

mu <- mu %>% mutate(class = ifelse(MeanY >= 150, "High", "Low"))

mu <- mu %>% mutate(class = ifelse(MeanY >= 120, "High",
                                   ifelse(MeanY >= 90 & MeanY < 125, "Mid", "Low")))

rem <- setdiff(mu$aoi_code, rownames(cl.pca))
`%nin%` = Negate(`%in%`)
mu <- mu %>% filter(aoi_code %nin% rem)
all.equal(mu$aoi_code, rownames(cl.pca))


cl.pca$class <- mu$class
cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
#cl.pca$class <- factor(cl.pca$class, levels = c("Low", "High"))

ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "class", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )




# two clusters
ans.cl <- get.kmeans(dist.gc, 2)


pca.res <- PCA(GC)
pca <-pca.res$ind$coord
pca <- as.data.frame(pca)

cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)

eigv <- pca.res$eig[,1:2]
var.pct <- round(eigv[,2], 1)

cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "cl.ass", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

mu <- r1_join %>% group_by(aoi_code)  %>% summarize(MeanY = mean(yield))

summary(mu$MeanY)

mu <- mu %>% mutate(class = ifelse(MeanY >= 150, "High", "Low"))

mu <- mu %>% mutate(class = ifelse(MeanY >= 95, "High",
                                   ifelse(MeanY >= 80 & MeanY < 95, "Mid", "Low")))

rem <- setdiff(mu$aoi_code, rownames(cl.pca))
`%nin%` = Negate(`%in%`)
mu <- mu %>% filter(aoi_code %nin% rem)
all.equal(mu$aoi_code, rownames(cl.pca))


cl.pca$class <- mu$class
cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
#cl.pca$class <- factor(cl.pca$class, levels = c("Low", "High"))

ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "class", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )






##########

cp <- readRDS("cp_gdcs.RData")
dt <- rbind(cp, r1_join)

id <- as.numeric(dt$ge_id)
Z <- GetZ(id, crop = "sorghum")

Y = dcast(dt, ge_id~aoi_code, mean, value.var='yield' )
rownames(Y) = Y$ge_id
Y = data.matrix(Y[,-1])
Y = Y[,which(colSums(!is.na(Y))>15)]


ii = intersect(rownames(Y),rownames(Z))
Y = Y[ii,]
Z = Z[ii,]

v = apply(Z,2,var, na.rm = TRUE)
Z = Z[,which(v>0.1)]
Z = Z-1
Z[is.na(Z)] = 0
Z %>% dim



fit = MRR3(Y,Z,cores=4)


GC = fit$GC
rownames(GC) = colnames(GC) = colnames(Y)
hc = hclust(as.dist(1-GC))
plot(hc)


dist.gc <- as.matrix(dist(1-GC))

get.optimum.k(dist.gc, 3)

ans.cl <- get.kmeans(dist.gc, 2)


pca.res <- PCA(GC)
pca <-pca.res$ind$coord
pca <- as.data.frame(pca)

cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)

eigv <- pca.res$eig[,1:2]
var.pct <- round(eigv[,2], 1)

cl.pca$cl.ass <- as.factor(cl.pca$cl.ass)
options(ggrepel.max.overlaps = Inf)
ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "cl.ass", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )








url <- "https://dsigitlab.research.auto.pioneer.com/alencar.xavier/btsatr.git"
devtools::install_git(url = url)

devtools::install_git(url = "https://dsigitlab.research.auto.pioneer.com/alencar.xavier/btsatr", subdir = "btsatr")

remotes::install_gitlab("alencar.xavier/btsatr")

devtools::install_git(
  "https://dsigitlab.research.auto.pioneer.com/alencar.xavier/btsatr", subdir = "btsatr",
  credentials = git2r::cred_user_pass("emmanuel.adeyemo", getPass::getPass())
)


###
r1_join <- r1_join %>% filter(year == 2019)
r1_join$ Loc = paste(r1_join$aoiname,r1_join$year)
Y = dcast(r1_join, ge_id~Loc, mean, value.var='yield' )
rownames(Y) = Y$ge_id
Y = data.matrix(Y[,-1])
Y = Y[,which(colSums(!is.na(Y))>100)]


ii = intersect(rownames(Y),rownames(Z))
Y = Y[ii,]
Z = Z[ii,]

v = apply(Z,2,var, na.rm = TRUE)
Z = Z[,which(v>0.1)]
Z = Z-1
Z[is.na(Z)] = 0
Z %>% dim


hist(v)


fit = MRR3(Y,Z,cores=4)


GC = fit$GC
rownames(GC) = colnames(GC) = colnames(Y)
hc = hclust(as.dist(1-GC))
plot(hc)

r1_join %>% group_by(aoiname) %>% summarize(MeanY = mean(yield))
fit$h2



sort(fit$h2)




r2_df <- DataFromJob("NAHP R2", crop = "sorghum", traits = "yield", What = "BLUE")
a <- SearchPrism("NAHP R2", crop = "sorghum") %>% filter(Type == "RSTG")# & AOI.Group == NA)
b <- searchModel('NAHP R2','sorghum')

library(btsatr)
dt <- DataFromJob("NAHP", crop = "sorghum", searchYear = 2019, What = "Phenotype", trait = c("PLTHT", "PLTHT_RS"))
yr <- c(2019, 2020, 2021)

dt_for_seyi <- dt %>% dplyr::filter(YEAR %in% yr)
write.table(dt_for_seyi, "Pltht_data_for_Seyi.txt", sep = "\t")

#===========================================================================================================
library(dplyr)
library(janitor)
  
f1 <- read.csv("f1_crosses.txt", sep = "\t", header = T) %>% clean_names



f1_high <- f1 %>% arrange(pedigree, -weight) %>% filter(duplicated(pedigree) == FALSE)
write.table(f1_high, "f1_high.txt", sep = "\t")







test <- read.csv("Daily_GHG_data_ea.csv", header=T, na.strings="") %>% clean_names
test <- test[1:10,]


#collective factor setting




arv <- function(dat){
  
  cols <- colnames(dat)
  for(i in 1:ncol(dat)){
    
    x <- cols[i]
    if(is.numeric(dat[,i])){
      dat[sapply(dat, is.numeric)]
      x_name <- paste0(x, "_arv")
      dat[,x_name] <-  round((100+(sapply(dat[,x], mean) - mean(dat[,x]))/sd(dat[,x])*10),0)
    }
  }
  
 
  return(dat)
}

arv(test)


f3 <- read.csv("F3_predictions.txt", sep = "\t", header = T) %>% clean_names

f3_arv <- arv(f3)
write.table(f3_arv, "F3_Predictions_ARV.txt", sep = "\t", row.names = F)



#############################  add SI to RSGT data ##############################

dta <- read.csv("R2.txt", sep = "\t", header = T) %>% clean_names()

dta$yield <- as.numeric(dta$yield)



dta <- dta %>% mutate(hy_ly = yield_highyld/yield_lowyld) %>% mutate(uearly_SI = (0.45 * yield_arv) + (0.1 * yield_lowyld_arv) 
                                                                     + (0.15 * ldgsev_arv) + (-0.05 * pltht_arv) + (-0.05 * pltht_rs_arv)
                                                                     + (0.2 * borldg_arv))

dta <- dta %>% mutate(yld_ly_ldg = yield_lowyld_arv/ldgsev_arv) %>% mutate(earlyStab_SI = (0.4 * yield_arv) + (0 * yield_lowyld_arv) 
                                                                     + (0.1 * ldgsev_arv) + (-0.15 * pltht_arv) + (-0.15 * pltht_rs_arv)
                                                                     + (0.15 * borldg_arv) + (0.05 * yld_ly_ldg))

dta <- dta %>% mutate(earlyTopE_SI = (0.45 * yield_arv) + (0 * yield_lowyld_arv) 
                                                                           + (0.05 * ldgsev_arv) + (-0.15 * pltht_arv) + (-0.15 * pltht_rs_arv)
                                                                           + (0.15 * borldg_arv) + (0.05 * yld_ly_ldg))

dta <- dta %>% mutate(midStab_SI = (0.4 * yield_arv) + (0 * yield_lowyld_arv) 
                      + (0.1 * ldgsev_arv) + (-0.1 * pltht_arv) + (-0.1 * pltht_rs_arv)
                      + (0.2 * borldg_arv) + (0.1 * yld_ly_ldg))

dta <- dta %>% mutate(midLateStab_SI = (0.55 * yield_arv) + (0 * yield_lowyld_arv) 
                      + (0.1 * ldgsev_arv) + (-0.1 * pltht_arv) + (-0.1 * pltht_rs_arv)
                      + (0.05 * borldg_arv) + (0.05 * yld_ly_ldg) + (0.05 * hy_ly))

dta <- dta %>% mutate(yld_ldg = yield_arv/ldgsev_arv) %>% mutate(midLateTopE_SI = (0.65 * yield_arv) + (0 * yield_lowyld_arv) 
                      + (0.05 * ldgsev_arv) + (-0 * pltht_arv) + (-0 * pltht_rs_arv)
                      + (0.15 * borldg_arv) + (0.1 * yld_ldg) + (0.05 * hy_ly))

write.table(dta, "R2_withSI.txt", sep = "\t")



################################################## Feb 2023
library(dplyr)
library(reshape2)
library(bWGR)
library(btsatr)
library(LaunchR)
library(corrplot)
library(stringr)
library(cluster)

setwd("C:/Users/gzp242/OneDrive - Corteva/Documents/Projects/R2_GDCS")
r2 <- read.csv("sorg_r2_cecs.txt", sep = "\t", header = T)

r2 %>% filter(str_detect(loc, "201999")) %>% tail

Z <- GetZ(unique(r2$ge_id), crop = "sorghum")

Y = dcast(r2, ge_id~loc, mean, value.var='yield' )
rownames(Y) = Y$ge_id
Y = data.matrix(Y[,-1])
Y = Y[,which(colSums(!is.na(Y))>99)] # removes environments with less than 100 entries

# >99 = only five locs were dropped

ii = intersect(rownames(Y),rownames(Z))
Y = Y[ii,]
Z = Z[ii,]

v = apply(Z,2,var, na.rm = TRUE)
Z = Z[,which(v>0.1)]
Z = Z-1
Z[is.na(Z)] = 0
Z %>% dim

fit = MRR3(Y,Z,cores=4)


GC = fit$GC
rownames(GC) = colnames(GC) = colnames(Y)
hc = hclust(as.dist(1-GC))
plot(hc)

h2 <- as.data.frame(fit$h2)
h2 <- cbind(loc = colnames(Y), h2) %>% rename(h2 = "fit$h2")
h2_0.2 <- h2 %>% filter(h2 >= 0.2)

id <-  intersect(rownames(GC),h2_0.2$loc)
GC <- GC[id, id]

GC_loc <- cbind(loc = rownames(GC), GC) %>% data.frame()
r2_for_join <- r2 %>% select(aoi_id, loc) %>% distinct(aoi_id, .keep_all = T)

GC.aoi <- GC_loc %>% inner_join(r2_for_join, by = "loc") 
rownames(GC.aoi) <- GC.aoi$loc 
GC_aoi <- GC.aoi %>% select(-loc)

id_new <- intersect(rownames(GC_aoi), colnames(GC_aoi))
id_new <- c(id_new, "aoi_id")
GC_aoi <- GC_aoi[id_new, id_new]

write.table(GC_aoi, "r2_genetic_correlation.txt", sep = "\t")

r2.data <- read.csv("r2_genetic_correlation.txt", sep = "\t", header = T, row.names = 1) %>% select(-aoi_id) %>% as.matrix()

require(ggplot2)
require(ggpubr)
require(ggfortify)
require(factoextra)
require(NbClust)
require(ggrepel)

#r2.data <- r2.data %>% mutate_if(is.factor, as.numeric) %>% textshape::column_to_rownames("locname")

gc <- 1 - r2.data
dist.gc <- as.matrix(dist(gc))

#' @aliases get optimum K for clustering - elbow method
get.optimum.k <- function(geno_distance, intercept_line){
  
  optimum_k <- fviz_nbclust(geno_distance, kmeans, method = "wss") +
    geom_vline(xintercept = intercept_line, linetype = 2)+  # include intercept number after optimum number has been determined
    labs(subtitle = "Elbow method")
  
  return(optimum_k)
}

get.optimum.k(dist.gc, 2)

#' @aliases get Kmeans
get.kmeans <- function(geno_distance, number_of_clusters){
  
  kmeans_out <- kmeans(geno_distance, number_of_clusters, nstart=2, iter.max=1000, algorithm="Hartigan-Wong") # options:algorith="Hartigan-Wong", "Lloyd", "Forgy"
  
  return(kmeans_out)
}

ans.cl <- get.kmeans(dist.gc, 2)

ans.cl <- get.kmeans(disT, 2)


library(FactoMineR)
pca.res <- PCA(r2.data)
pca <-pca.res$ind$coord
pca <- as.data.frame(pca)
plot(pca$Dim.1, pca$Dim.2)

cl.pca <- cbind(pca, cl.ass = ans.cl$cluster)

#' @aliases plot kmeans
plot.kmeans <- function(cluster_with_pca, cluster_number){
  
  cluster_with_pca[,cluster_number] <- as.factor(cluster_with_pca[,cluster_number])
  
  kmeans_plot <- ggplot(cluster_with_pca, aes(x=Dim.1, y=Dim.2, shape=cluster_with_pca[,cluster_number], color=cluster_with_pca[,cluster_number])) +
    geom_point()
  
  kmeans_bw <- kmeans_plot + theme_bw()
  
  return(kmeans_bw)
}



plot.kmeans(cl.pca, "cl.ass")


fviz_cluster(ans.cl, data = gc,
             palette = c("#4477AA", "#EE9988", "#E7B900"), 
             geom = c("point", "text"),
             ellipse.type = "convex", 
             labelsize = 8,
             ggtheme = theme_bw()
)


fviz_cluster(ans.cl, data = test,
             palette = c("#4477AA", "#EE9988", "#E7B900"), 
             geom = c("point", "text"),
             ellipse.type = "convex", 
             labelsize = 8,
             ggtheme = theme_bw()
)

pm <- eclust(gc,FUNcluster="pam", k=2,hc_metric = "euclidean")
pm.cl <- pam(dist.gc, 2)

fviz_cluster(pm.cl, data = gc,
             palette = c("#4477AA", "#EE9988", "#E7B900"), 
             geom = c("point", "text"),
             ellipse.type = "convex", 
             labelsize = 6,
             #repel = TRUE,
             ggtheme = theme_bw()
)

eigv <- pca.res$eig[,1:2]
var.pct <- round(eigv[,2], 1)

require(ggpubr)


r2_select <- r2 %>% filter(loc %in% rownames(cl.pca))

mu <- r2_select %>% group_by(loc)  %>% summarize(MeanY = mean(yield))

summary(mu$MeanY)


mu <- mu %>% mutate(class = ifelse(MeanY >= 115, "High", "Low"))

cl.pca$class <- mu$class

options(ggrepel.max.overlaps = Inf)

ggscatter(
  cl.pca, x = "Dim.1", y = "Dim.2", 
  color = "cl.ass", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "class", 
  size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", var.pct[1], "% )" ),
  ylab = paste0("Dim 2 (", var.pct[2], "% )" )
) +
  geom_text_repel(
    data = cl.pca,#subset(cl.pca, class == "High"),
    #aes(label = rownames(subset(cl.pca, class == "Low"))),
    aes(label = rownames(cl.pca)),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )







col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(mtx_corr, method="color", col=col(200),  
         type="upper", #order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         #p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)



corr_simple <- function(data=df,sig=0.5){
  #run a correlation and drop the insignificant ones
  #corr <- cor(df_cor)
  corr <- df
  #prepare to drop duplicates and correlations of 1     
  corr[lower.tri(corr,diag=TRUE)] <- NA 
  #drop perfect correlations
  corr[corr == 1] <- NA 
  #turn into a 3-column table
  corr <- as.data.frame(as.table(corr))
  #remove the NA values from above 
  corr <- na.omit(corr) 
  #select significant values  
  corr <- subset(corr, abs(Freq) > sig) 
  #sort by highest correlation
  corr <- corr[order(-abs(corr$Freq)),] 
  #print table
  print(corr)
  #turn corr back into matrix in order to plot with corrplot
  mtx_corr <- reshape2::acast(corr, Var1~Var2, value.var="Freq")
  
  #plot correlations visually
  corrplot(mtx_corr, is.corr=FALSE, tl.col="black", na.label=" ")
}
corr_simple(GC)


corr <- GC
corr[lower.tri(corr,diag=TRUE)] <- NA 
corr[corr == 1] <- NA 
corr <- as.data.frame(as.table(corr))

corr <- na.omit(corr) 
corr <- subset(corr, abs(Freq) > 0.8) 
corr <- corr[order(-abs(corr$Freq)),] 
mtx_corr <- reshape2::acast(corr, Var1~Var2, value.var="Freq")
corrplot(mtx_corr, is.corr=FALSE, tl.col="black", na.label=" ")
