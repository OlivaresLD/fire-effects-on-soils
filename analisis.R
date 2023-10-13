####################################################################
####################################################################
## Script para análisis estadísticos y visualización de gráficos
## de Ks y repelencia para artículo de la Sp.J.S.Sc.
## Desarrollado por: Luis Daniel Olivares Martínez
## 13/04/2023
## Universidad Miguel Hernández
####################################################################
####################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # directorio a la ubicación del archivo

art <- read.csv("_Ks_repe_SITIOS.csv", encoding = "latin1")
art$tipoR_0 <- factor(art$tipoR, levels = c("C(0)","Sup(1)","Sup(2)","Sub(1)","Sub(2)"))
art$tipoR <- factor(art$tipoR, levels = c("C(0)","S(1)","S(2)","G(1)","G(2)"))
art$trat <- factor(art$trat)
art$rec <- factor(art$rec)
repe <- subset(art, !is.na(a.max))
ks <- subset(art, !is.na(Ks.ms))
repe$IRDI[is.na(repe$IRDI)] <- 0

#### incluir todas las repelencias ####
# re <- read.csv("Hidrofob_ensayos.csv", encoding = "latin1")
# ref <- subset(re, re$Sitio != 'Cru' & !re$ensayo %in% art$id)
# write.csv(ref,"repe sin ks asociada.csv", row.names = F, fileEncoding = "latin1") #ya están incluidas en "_Ks_repe_SITIOS.csv"

#### estadística descriptiva ####
library(tidyverse)

boxplot(ks$Ks.mmhr~ks$tipoR)
boxplot(log(ks$Ks.mmhr)~ks$tipoR)
boxplot(ks$a.m~ks$tipoR)
boxplot(ks$PF_Macroporos~ks$tipoR)
boxplot(ks$N.m2_Macroporos~ks$tipoR)
boxplot(log(ks$N.m2_Macroporos)~ks$tipoR)

fi <- which(!grepl('0.5r',ks$id))# filtro control sólo en claros
table(ks$tipoR[fi])
table(ks$tipoR)
boxplot(ks$Ks.mmhr[fi]~ks$tipoR[fi])
boxplot(ks$N.m2_Macroporos[fi]~ks$tipoR[fi])
# ks <- subset(ks, !grepl('0.5r',ks$id)) # aplicando el filtro para tener solo claros en los análisis (a 0.5m de radio dista de la metodología para los otros sitios y baja mucho la Ks)

{resks <- ks %>% group_by(tipoR) %>% dplyr::summarise(Ks.mmh = mean(Ks.mmhr) %>% round(2), Ks_de = sd(Ks.mmhr) %>% round(2), a = mean(a.m) %>% round(2), a_de = sd(a.m) %>% round(2), PF = mean(PF_Macroporos) %>% round(2), PF_de = sd(PF_Macroporos) %>% round(2), N = mean(N.m2_Macroporos) %>% round(2), N_de = sd(N.m2_Macroporos) %>% round(2))
  # tapply(tabla.ks$Ks.mmhr, tabla.ks$Piso, mean)
  # tapply(tabla.ks$Ks.mmhr, tabla.ks$Piso, sd)
  # tapply(tabla.ks$Ks.mmhr, tabla.ks$Piso, max)
  # paste0(round(mean(tabla.ks$Ks.mmhr),2),"±",round(sd(tabla.ks$Ks.mmhr),2))
  resks <- cbind(resks, table(ks$tipoR))
  resks$Var1 <- NULL
  resks$Ks.mmh_e <- (resks$Ks.mmh / resks$Freq ^ 0.5) %>% round(2)
  resks$a_e <- (resks$a / resks$Freq ^ 0.5) %>% round(2)
  resks$PF_e <- (resks$PF / resks$Freq ^ 0.5) %>% round(2)
  resks$N_e <- (resks$N / resks$Freq ^ 0.5) %>% round(1)
}
write.csv(resks, "ks_5sitios.csv", row.names = F, fileEncoding = "latin1")

boxplot(repe$IRDI~repe$tipoR) # IRDI más alto suelos más repelentes considerando todos sus estados de humedad
boxplot(repe$w.min~repe$tipoR)
boxplot(repe$w.max~repe$tipoR)
boxplot(repe$a.max~repe$tipoR)
boxplot(repe$a.105~repe$tipoR)
boxplot(repe$S~repe$tipoR)

{resrepe <- repe %>% group_by(tipoR) %>% dplyr::summarise(irdi = mean(IRDI) %>% round(2), irdi_de = sd(IRDI) %>% round(2), wmin = mean(w.min) %>% round(3), wmin_de = sd(w.min) %>% round(3), wmax = mean(w.max) %>% round(3), wmax_de = sd(w.max) %>% round(3), amax = mean(a.max) %>% round(2), amax_de = sd(a.max) %>% round(2), a105 = mean(a.105) %>% round(2), a105_de = sd(a.105) %>% round(2), S = mean(S, na.rm = T) %>% round(2), S_de = sd(S, na.rm = T) %>% round(2))
  resrepe$S_de <- tapply(repe$S, repe$tipoR, sd) %>% round(2)
  resrepe$Type <- factor(c('Control','Surface','Surface','Ground','Ground'), levels = c('Control','Surface','Ground'))
  resrepe$Recurrence <- factor(c(0,1,2,1,2))
  names(resrepe)[seq(2,13,2)] <- c('IRDI','w.min','w.max','a.max','a.105','S')
  names(resrepe)[seq(3,13,2)] <- paste0(c('IRDI','w.min','w.max',
                                          'a.max','a.105','S'),'_de')
  resrepe <- cbind(resrepe, table(repe$tipoR))
  resrepe$Var1 <- NULL
  resrepe$w.min_e <- resrepe$w.min_de / resrepe$Freq ^ 0.5
  resrepe$a.max_e <- resrepe$a.max_de / resrepe$Freq ^ 0.5
  resrepe$S_e <- resrepe$S_de / resrepe$Freq ^ 0.5
  resrepe$a.105_e <- resrepe$a.105_de / resrepe$Freq ^ 0.5
  resrepe$w.max_e <- resrepe$w.max_de / resrepe$Freq ^ 0.5
  resrepe$IRDI_e <- resrepe$IRDI_de / resrepe$Freq ^ 0.5
  }
write.csv(resrepe, "repe_5sitios.csv", row.names = F, fileEncoding = "latin1")

#### estadísticos ####
library(agricolae)
library(conover.test)
library(fitdistrplus)

xxx <- aggregate(repe$a.105, by= list(P=repe$Piso, S=repe$Sitio), FUN=mean); friedman.test(x ~ S | P, data = xxx[c(-1,-11),]) #considerando los pisos, no significativa
rm(xxx)

#### extracción de los modelos ANOVA Ks y repe ####
vars <- c('Ks.mmhr', 
          'a.m', 
          'PF_Macroporos', 
          'N.m2_Macroporos', 
          'w.min','a.max','S','a.105','w.max','IRDI')
vars2 <- c(rep('ks',4),rep('repe',6))
log_si <- c(rep(1,4),0,1,0,1,0,0)
lom <- length(vars)
mods <- array(NA,c(lom*4,7))
mods <- data.frame(mods, stringsAsFactors = F)
names(mods) <- c('Type', 'Recurrence', 'Type * Rec', 'Formula', 'F_Type','F_Rec','F_T_R')
tryCatch( expr = for(i in 1:lom){ # condicional para anova con log y para anova sin log según log_si en i, comparando con los modelos de trat, rec y trat*rec
 
  if (log_si[i] == 0) {
    ano <- anova(aov(eval(parse(text = vars[i])) ~ trat * rec, data = eval(parse(text =vars2[i]))))
    ano2 <- anova(aov(eval(parse(text = vars[i])) ~ trat, data = eval(parse(text =vars2[i]))))
    ano3 <- anova(aov(eval(parse(text = vars[i])) ~ rec, data = eval(parse(text =vars2[i]))))
    ano4 <- anova(aov(eval(parse(text = vars[i])) ~ trat + rec, data = eval(parse(text =vars2[i]))))
  }else{
    ano <- anova(aov(log(eval(parse(text = vars[i]))) ~ trat * rec, data = eval(parse(text =vars2[i]))))
    ano2 <- anova(aov(log(eval(parse(text = vars[i]))) ~ trat, data = eval(parse(text =vars2[i]))))
    ano3 <- anova(aov(log(eval(parse(text = vars[i]))) ~ rec, data = eval(parse(text =vars2[i]))))
    ano4 <- anova(aov(log(eval(parse(text = vars[i]))) ~ trat + rec, data = eval(parse(text =vars2[i]))))
  }
  
  mods[i,1:3] <- ano[1:3,5]
  mods[i,5:7] <- ano[1:3,4]
  mods[i,4] <- paste(vars[i],'~ Type * Recurrence')
  
  mods[i+lom,1] <- ano2[1,5]
  mods[i+lom,5] <- ano2[1,4]
  mods[i+lom,4] <- paste(vars[i],'~ Type')
  
  mods[i+lom*2,2] <- ano3[1,5]
  mods[i+lom*2,6] <- ano3[1,4]
  mods[i+lom*2,4] <- paste(vars[i],'~ Recurrence')
  
  mods[i+lom*3,1:3] <- ano4[1:3,5]
  mods[i+lom*3,5:7] <- ano4[1:3,4]
  mods[i+lom*3,4] <- paste(vars[i],'~ Type + Recurrence')
  
}, warning=NULL,error=function(e)paste("Ensayo",vars[i],"con problmeas en fila ",i))

mini_mods <- subset(mods, Recurrence < 0.05 | Type < 0.05 | `Type * Rec` < 0.05)
setwd("H:/Mi unidad/oliolivares@ciencias.unam.mx 2022-11-11 17 22/avances papers/resiliencia incendios AN/resultados art/outputs")
write.csv(mini_mods, "Modelos significativos.csv")

##Relación Ks ~ SOM + n_m2_Macroporos
tapply(art$Ks.mmhr, art$tipoR, function(x)mean(x,na.rm=T) %>% round(2))
tapply(art$IRDI, art$tipoR, function(x)mean(x,na.rm=T) %>% round(2))
tapply(art$N.m2_Macroporos, art$tipoR, function(x)mean(x,na.rm=T) %>% round(0))
descdist(ks$N.m2_Macroporos, discrete = FALSE) # entre uniforme y normal

soc_ori <- read.csv("SOC.csv")
soc <- subset(soc_ori, layer == "mineral")
plot(soc$Ks~soc$OM_.+soc$n_m2_M)
plot(soc$S~soc$OM_.)
plot(soc$IRDI~soc$OM_.)
plot(soc$w.min~soc$OM_.)
plot(Ks~OM_.+n_m2_M, data=soc[which(soc$layer == "mineral"),])
anova(aov(log(soc_ori$Ks) ~ soc_ori$type * soc_ori$n_m2_M)) # Significativo
anova(aov(log(soc_ori$Ks) ~ soc_ori$type + soc_ori$OM_. + soc_ori$n_m2_M)) # Significativo
anova(aov(log(soc$Ks) ~ soc$type + soc$OM_. + soc$n_m2_M)) #faltan datos para este modelo
anova(aov(log(soc$Ks) ~ soc$OM_.*soc$n_m2_M)) # No significativo
anova(aov(log(soc$Ks) ~ log(soc$OM_.) + log(soc$n_m2_M))) # No significativo
glm(soc$Ks~soc$OM_.*soc$n_m2_M, family = Gamma) %>% summary()

##hipótesis infiltración ~ recurrencia + tipo de incendio
# ANOVA's y Post Hoc
anova(aov(ks$Ks.mmhr~ks$trat*ks$rec)) # Significativo solo trat, no rec
anova(aov(log(ks$Ks.mmhr) ~ ks$trat*ks$rec)) # Significativos AMBOS sin interacción
anova(aov(ks$a.m~ks$trat*ks$rec)) # Significativos trat e interacción de ambos
anova(aov(log(ks$a.m) ~ ks$trat*ks$rec)) # Significativos trat e interacción de ambos
anova(aov(ks$PF_Macroporos~ks$trat*ks$rec)) # Significativos trat e interacción de ambos
anova(aov(ks$N.m2_Macroporos~ks$trat*ks$rec)) # Significativos rec e interacción de ambos
anova(aov(log(ks$N.m2_Macroporos) ~ ks$trat*ks$rec)) # Significativos AMBOS e interacción

pairwise.t.test(ks$Ks.mmhr, ks$tipoR, p.adj = "holm") # c(a,b,b,b,b)
pairwise.t.test(log(ks$Ks.mmhr), ks$tipoR, p.adj = "holm") # c(a,b,bc,bd,bc) #
# pairwise.t.test(log(ks$Ks.mmhr), ks$tipoR, p.adj = "bon") # c(a,b,b,b,b)
pairwise.t.test(ks$a.m, ks$tipoR, p.adj = "holm") # c(ab,a,ab,a,ac)
pairwise.t.test(log(ks$a.m), ks$tipoR, p.adj = "holm") # c(a,b,bc,bc,bd) #
pairwise.t.test(ks$PF_Macroporos, ks$tipoR, p.adj = "holm") # c(a,b,cb,b,db)
pairwise.t.test(log(ks$PF_Macroporos), ks$tipoR, p.adj = "holm") # c(a,b,bc,bc,bd) #
pairwise.t.test(ks$N.m2_Macroporos, ks$tipoR, p.adj = "holm") # c(a,a,a,a,a)
pairwise.t.test(log(ks$N.m2_Macroporos), ks$tipoR, p.adj = "holm") # c(ab,a,a,ac,ab) #

pairwise.t.test(ks$Ks.mmhr, ks$trat, p.adj = "holm") # c(a,a,a,a,a)
TukeyHSD(aov(ks$a.m ~ ks$trat * ks$rec))
pairwise.t.test(ks$N.m2_Macroporos, ks$rec, p.adj = "holm") # c(a,a,a,a,a)
TukeyHSD(aov(ks$N.m2_Macroporos ~ ks$trat * ks$rec))

c1 <- HSD.test(aov(log(Ks.mmhr) ~ tipoR, data=ks),"tipoR", group=T)$groups
c2 <- HSD.test(aov(log(a.m) ~ tipoR, data=ks),"tipoR", group=T)$groups # 
c3 <- HSD.test(aov(log(PF_Macroporos) ~ tipoR, data=ks),"tipoR", group=T)$groups
c4 <- HSD.test(aov(log(N.m2_Macroporos) ~ tipoR, data=ks),"tipoR", group=T)$groups # 
ks_g_l <- list(c1,c2,c3,c4)
ks_g <- ks_g_l %>% 
  # map(~select(.x)) %>%
  map(rownames_to_column) %>%
  purrr::reduce(dplyr::full_join, by='rowname')
names(ks_g)[seq(2,9,by=2)] <- paste0(c('Ks.mmhr','a.m','PF_Macroporos','N.m2_Macroporos'),'_mean')
names(ks_g)[seq(3,9,by=2)] <- c('Ks.mmhr','a.m','PF_Macroporos','N.m2_Macroporos')
ks_g$rowname <- factor(ks_g$rowname, levels = c("C(0)","Sup(1)","Sup(2)","Sub(1)","Sub(2)"))
ks_g <- with(ks_g, ks_g[order(ks_g$rowname),])
write.csv(ks_g, "ks_groups_corrected.csv", row.names = F, fileEncoding = "latin1")

##hipótesis repelencia ~ recurrencia + tipo de incendio
# Distribución de datos
descdist(repe$w.min, discrete = FALSE) # entre uniforme y normal
descdist(repe$a.max, discrete = FALSE) # uniforme
descdist(repe$S, discrete = FALSE) # entre uniforme y normal
descdist(repe$a.105, discrete = FALSE) # entre uniforme, normal y beta
 # beta dist mod <- glm(y ~ x1 + x2, data = foo, family = binomial(link = "logit"))
descdist(repe$w.max, discrete = FALSE) # entre uniforme, normal y beta
descdist(repe$IRDI, discrete = FALSE) # entre "logística" y lognormal

fitdist(repe$IRDI, "logis")$aic
fitdist(repe$IRDI+1, "lnorm")$aic
fitdist(repe$IRDI, "norm")$aic

fitdist(repe$w.max, "unif")$aic
fitdist(repe$w.max, "norm")$aic
fitdist(repe$w.max+1, "weibull")$aic

fitdist(repe$a.105, "unif")$aic
fitdist(repe$a.105/120, "beta")$aic
fitdist(repe$a.105, "norm")$aic

# GLM
glm(repe$w.min ~ repe$trat*repe$rec, family = gaussian) %>% summary()
glm(repe$a.max ~ repe$trat*repe$rec, family = gaussian) %>% summary()
glm(repe$S ~ repe$trat*repe$rec, family = gaussian) %>% summary()
glm(repe$a.105 ~ repe$trat*repe$rec, family = quasipoisson) %>% summary()
glm(repe$w.max ~ repe$trat*repe$rec, family = quasipoisson) %>% summary()
glm(repe$IRDI ~ repe$trat*repe$rec, family = quasipoisson) %>% summary()

# wilcox.test(x=repe$a.105, repe$PISO.EQ | repe$Sitio) # ver como copartar con varias vías
# wilcox.test(repe$a.105 ~ repe$PISO.EQ+repe$sitio)
# Prueba en bloques no paramétrica (Test de Friedman)
kruskal.test(repe$a.105 ~ repe$tipoR)
kruskal.test(repe$a.105 ~ repe$rec)
conover.test(repe$a.105, repe$rec, method="Bonferroni")
conover.test(repe$a.105, repe$tipoR, method="Bonferroni")

xxx <- aggregate(repe$a.max, by= list(P=repe$trat, S=repe$rec), FUN=mean); friedman.test(x ~ S | P, data=xxx[-1,]) #considerando recurrencia y tipo de incendio. No significativa
rm(xxx)
xxx <- aggregate(repe$w.min, by= list(P=repe$trat, S=repe$rec), FUN=mean); friedman.test(x ~ S | P, data=xxx[-1,]) #considerando recurrencia y tipo de incendio. No significativa
rm(xxx)
xxx <- aggregate(repe$IRDI, by= list(P=repe$trat, S=repe$rec), FUN=function(x)mean(x,na.rm=T)); friedman.test(x ~ S | P, data=xxx[-1,]) #considerando recurrencia y tipo de incendio. No significativa
rm(xxx)
xxx <- aggregate(repe$S, by= list(P=repe$trat, S=repe$rec), FUN=mean); friedman.test(x ~ S | P, data=xxx[-1,]) #considerando recurrencia y tipo de incendio. No significativa
rm(xxx)
xxx <- aggregate(repe$w.max, by= list(P=repe$trat, S=repe$rec), FUN=mean); friedman.test(x ~ S | P, data=xxx[-1,]) #considerando recurrencia y tipo de incendio. No significativa
rm(xxx)
xxx <- aggregate(repe$a.105, by= list(P=repe$trat, S=repe$rec), FUN=mean); friedman.test(x ~ S | P, data=xxx[-1,]) #considerando recurrencia y tipo de incendio. No significativa
rm(xxx)
xxx <- aggregate(repe$Y90.max, by= list(P=repe$trat, S=repe$rec), FUN=mean); friedman.test(x ~ S | P, data=xxx[-1,]) #considerando recurrencia y tipo de incendio. No significativa
rm(xxx)

# ANOVAs
anova(aov(repe$w.min~repe$trat)) # Significativa
anova(aov(repe$w.min~repe$trat*repe$rec)) # Significativo solo trat, no rec
anova(aov(log(repe$w.min+1)~repe$trat*repe$rec)) # Significativo solo trat, no rec

anova(aov(repe$a.max~repe$trat*repe$rec)) # Significativos AMBOS, no interacción
anova(aov(log(repe$a.max)~repe$trat*repe$rec)) # Significativos AMBOS, no interacción

anova(aov(repe$a.105~repe$trat*repe$rec)) # Significativo solo trat, no rec
anova(aov(log(repe$a.105)~repe$trat*repe$rec)) # Significativos AMBOS, no interacción

anova(aov(repe$S~repe$trat*repe$rec)) # Significativo solo trat, no rec
anova(aov(repe$w.max~repe$trat*repe$rec)) # Significativo solo trat, no rec
anova(aov(repe$IRDI~repe$trat*repe$rec)) # Significativos trat e interacción

anova(aov(repe$w.min ~ repe$tipoR)) # Significativo
anova(aov(repe$a.max ~ repe$tipoR)) # Significativo
anova(aov(repe$S ~ repe$tipoR)) # Significativo
anova(aov(repe$a.105 ~ repe$tipoR)) # Significativo
anova(aov(repe$w.max ~ repe$tipoR)) # Significativo 
anova(aov(repe$IRDI ~ repe$tipoR)) # Significativo

# Post Hoc ANOVAs
TukeyHSD(aov(repe$IRDI ~ repe$tipoR))
TukeyHSD(aov(repe$IRDI ~ repe$trat * repe$rec))
# LSD.test(aov(repe$IRDI ~ repe$trat * repe$rec), "tipoR", p.adj= "bon")
# HSD.test(aov(repe$IRDI ~ repe$tipoR), "tipoR")
# LSD.test(aov(repe$IRDI ~ repe$tipoR), "tipoR", p.adj= "bon")

TukeyHSD(aov(repe$a.max ~ repe$tipoR))
TukeyHSD(aov(repe$a.max ~ repe$trat + repe$rec))
LSD.test(aov(repe$a.max ~ repe$trat * repe$rec), "tipoR", p.adj= "bon")
HSD.test(aov(repe$a.max ~ repe$tipoR), "tipoR", p.adj= "bon")

pairwise.t.test(repe$w.min, repe$tipoR, p.adj = "bonf")
pairwise.t.test(log(repe$a.max), repe$tipoR, p.adj = "holm")
pairwise.t.test(repe$S, repe$tipoR, p.adj = "holm")
pairwise.t.test(log(repe$a.105), repe$tipoR, p.adj = "holm")
pairwise.t.test(repe$w.max, repe$tipoR, p.adj = "holm")
pairwise.t.test(repe$IRDI, repe$tipoR, p.adj = "holm")

c1 <- HSD.test(aov(w.min ~ tipoR, data=repe),"tipoR", group=T)$groups # , console=T
c2 <- HSD.test(aov(log(a.max) ~ tipoR, data=repe),"tipoR", group=T)$groups # 
c3 <- HSD.test(aov(S ~ tipoR, data=repe),"tipoR", group=T)$groups # , console=T
c4 <- HSD.test(aov(log(a.105) ~ tipoR, data=repe),"tipoR", group=T)$groups # 
c5 <- HSD.test(aov(w.max ~ tipoR, data=repe),"tipoR", group=T)$groups # , console=T
c6 <- HSD.test(aov(IRDI ~ tipoR, data=repe),"tipoR", group=T)$groups # , console=T

repe_g_l <- list(c1,c2,c3,c4,c5,c6)
repe_g <- repe_g_l %>% 
  # map(~select(.x)) %>%
  map(rownames_to_column) %>%
  purrr::reduce(dplyr::full_join, by='rowname')
# repe_g2 <- cbind(c1,c2,c3,c4,c5,c6)
names(repe_g)[seq(2,13,by=2)] <- paste0(c('w.min','a.max','S','a.105','w.max','IRDI'),'_mean')
names(repe_g)[seq(3,13,by=2)] <- c('w.min','a.max','S','a.105','w.max','IRDI')
repe_g$rowname <- factor(repe_g$rowname, levels = c("C(0)","Sup(1)","Sup(2)","Sub(1)","Sub(2)"))
repe_g <- with(repe_g, repe_g[order(repe_g$rowname),])
rm(list = c('c1','c2','c3','c4','c5','c6'))

#### Grafiquirris ggplot2 ####
library(ggpubr)
# ggplot(resrepe, aes(x = tipoR))+
#   geom_errorbar(aes(ymin = IRDI - 5, ymax = IRDI + IRDI_de), 
#                   linewidth = 0.55, width = 0.2, color = "gray38")+
#   geom_col(aes(y = IRDI, fill = Type))+
#   scale_fill_manual(values = c("aquamarine4","salmon4","cornsilk3"))+
#   scale_y_continuous(limits = c(0,110))+
#   geom_text(label = repe_g$IRDI, y = 110, colour = "gray11", size = 4.5)+
#   xlab(NULL)+
#   theme_light()+
#   theme(legend.position = "none", axis.text.x=element_blank())
#   # theme(legend.position = "bottom")

pA <- ggplot(resrepe, aes(x = tipoR))+
  scale_fill_manual(values = c("aquamarine4","salmon4","cornsilk3"))+
  xlab(NULL)+
  theme_light() +
  theme(axis.text.x=element_blank() ) #, 
        # legend.position = "none")

pB <- ggplot(resrepe, aes(x = tipoR))+
  scale_fill_manual(values = c("aquamarine4","salmon4","cornsilk3"))+
  xlab(NULL)+
  theme_light()+
  theme(legend.position = "bottom", axis.text.x = element_text(size=rel(1.3)))

# Versión con Desviaciones Estándar
w.min <- pA+
  geom_errorbar(aes(ymin = w.min - .1, ymax = w.min + w.min_de), 
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = w.min, fill = Type))+
  coord_cartesian(ylim = c(0,1.50))+
  geom_text(label = repe_g$w.min, y = 1.45, colour = "gray11", size = 4.5)
a.max <- pA+
  geom_errorbar(aes(ymin = a.max - 5, ymax = a.max + a.max_de), 
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = a.max, fill = Type))+
  coord_cartesian(ylim = c(60,125))+
  geom_text(label = repe_g$a.max, y = 123, colour = "gray11", size = 4.5)
S <- pA+
  geom_errorbar(aes(ymin = S - 5, ymax = S + S_de),
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = S, fill = Type))+
  coord_cartesian(ylim = c(0,90))+
  geom_text(label = repe_g$S, y = 85, colour = "gray11", size = 4.5)
a.105 <- pA+
  geom_errorbar(aes(ymin = a.105 - 5, ymax = a.105 + a.105_de), 
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = a.105, fill = Type))+
  coord_cartesian(ylim = c(60,125))+
  geom_text(label = repe_g$a.105, y = 117, colour = "gray11", size = 4.5)
w.max <- pA+
  geom_errorbar(aes(ymin = w.max - 0.05, ymax = w.max + w.max_de), 
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = w.max, fill = Type))+
  coord_cartesian(ylim = c(0,1.0))+
  geom_text(label = repe_g$w.max, y = 0.80, colour = "gray11", size = 4.5)
IRDI <- pA+
  geom_errorbar(aes(ymin = IRDI - 5, ymax = IRDI + IRDI_de), 
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = IRDI, fill = Type))+
  coord_cartesian(ylim = c(0,110))+
  geom_text(label = repe_g$IRDI, y = 110, colour = "gray11", size = 4.5)

# Versión con Errores Estándar
w.min <- pA+
  geom_errorbar(aes(ymin = w.min - .1, ymax = w.min + w.min_e), 
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = w.min, fill = Type))+
  coord_cartesian(ylim = c(0,1.2))+
  geom_text(label = repe_g$w.min, y = 1.16, colour = "gray11", size = 4.5)
a.max <- pA+
  geom_errorbar(aes(ymin = a.max - 5, ymax = a.max + a.max_e), 
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = a.max, fill = Type))+
  coord_cartesian(ylim = c(60,120))+
  geom_text(label = repe_g$a.max, y = 118, colour = "gray11", size = 4.5)
S <- pA+
  geom_errorbar(aes(ymin = S - 5, ymax = S + S_e),
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = S, fill = Type))+
  coord_cartesian(ylim = c(0,90))+
  geom_text(label = repe_g$S, y = 85, colour = "gray11", size = 4.5)
a.105 <- pB+
  geom_errorbar(aes(ymin = a.105 - 2, ymax = a.105 + a.105_e),
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = a.105, fill = Type))+
  coord_cartesian(ylim = c(80,110))+
  geom_text(label = repe_g$a.105, y = 108, colour = "gray11", size = 4.5)
w.max <- pA+
  geom_errorbar(aes(ymin = w.max - 0.05, ymax = w.max + w.max_e), 
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = w.max, fill = Type))+
  coord_cartesian(ylim = c(0,.5))+
  geom_text(label = repe_g$w.max, y = 0.40, colour = "gray11", size = 4.5)
IRDI <- pA+
  geom_errorbar(aes(ymin = IRDI - 5, ymax = IRDI + IRDI_e), 
                linewidth = 0.55, width = 0.2, color = "gray38")+
  geom_col(aes(y = IRDI, fill = Type))+
  coord_cartesian(ylim = c(0,100))+
  geom_text(label = repe_g$IRDI, y = 97, colour = "gray11", size = 4.5)

## plot imagen y guardado en sistema
setwd("H:/Mi unidad/oliolivares@ciencias.unam.mx 2022-11-11 17 22/avances papers/resiliencia incendios AN/resultados art/outputs")
png("Repe_5sites.png", units="cm", width=10, height=17, res=200)
ggarrange(IRDI, w.min,S, a.max, a.105, common.legend = T, legend = "bottom", ncol = 1)
dev.off()
tiff("Repe_5sites.tiff", units="cm", width=10, height=17, res=200)
ggarrange(IRDI, w.min,S, a.max, a.105, common.legend = T, legend = "bottom", ncol = 1)
dev.off()

# plot(log(ks$Ks.mmhr)~ks$IRDI, col=c("green","blue","red","gray55")[as.numeric(ks$tipoR)]); legend("bottomright",legend=levels(ks$tipoR), col = c("green","blue","red","gray55"),cex= 0.8,pch=21)
# 
# plot(log(ks$Ks.mmhr)~ks$a.max, col=c("green","blue","red","gray55")[as.numeric(ks$tipoR)]); legend("bottomright",legend=levels(ks$tipoR), col = c("green","blue","red","gray55"),cex= 0.8,pch=21)