# Analyses axe 1
#version finale

#loader les données
setwd("~/Documents/Maîtrise/Données Axe 1")

lonicera <- read.table("Lonicera.txt", header=TRUE)
lonicera$site <- as.factor(lonicera$site)
lonicera$inoculum <- as.factor(lonicera$inoculum)

#communautés fongiques :
# ~~~ À changer ~~~
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/fungal metacom.RData")
metacom <- metacom[-41,]
colnames(metacom) <- paste("champi",1:ncol(metacom),sep="")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# One-way ANOVAs de la performance des plants ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(datasets)
library(ggplot2)
library(multcompView)
library(dplyr)
library(carData)
library(car)

#Cette analyse se fait avec le jeu de données de performance originel, et non
#avec les moyennes qui seront utilisées plus tard
#loader les données ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
croissance <- read.table("croissance.txt", header=TRUE)
croissance$site <- as.factor(croissance$site)

#ANOVAs initiales sur données brutes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d.racines.aov <- aov(croissance$diam.racine ~ croissance$site)
summary(d.racines.aov)
#p = 0.272
racines.aov <- aov(croissance$l.racines ~ croissance$site)
summary(racines.aov)
#p = 0.0177 *
tiges.aov <- aov(croissance$masse.tige ~ croissance$site)
summary(tiges.aov)
#p = 0.0108 *

# Vérifier les conditions d'application ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#normalité des résidus

#diamètre des racines ~~~~~~~~~~~~~~~
shapiro.test(resid(d.racines.aov))
#p-value = 6.266e-14
#résidus ne suivent pas une distribution normale
hist(croissance$diam.racine.moy)
qqPlot(croissance$diam.racine.moy)

#longueur des racines ~~~~~~~~~~~~~~~
shapiro.test(resid(racines.aov))
#p-value = 0.002748
#résidus ne suivent pas une distribution normale
hist(croissance$l.racines)
qqPlot(croissance$l.racines)

#transformation racine carrée
croissance$l.racines.V <- sqrt((croissance$l.racines))
hist(croissance$l.racines.V)
qqPlot(croissance$l.racines.V)
#semble mieux


#masse des tiges ~~~~~~~~~~~~~~~
shapiro.test(resid(tiges.aov))
#p-value = 1.617e-07
#résidus ne suivent pas une distribution normale
hist(croissance$masse.tige)
qqPlot(croissance$masse.tige)

#transformation racine carrée
croissance$masse.tige.V <- sqrt((croissance$masse.tige))
hist(croissance$masse.tige.V)
qqPlot(croissance$masse.tige.V)
#semble mieux

#ANOVAs sur données tranformées racine carré ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
racinesV.aov <- aov(croissance$l.racines.V ~ croissance$site)
summary(racinesV.aov)
#p = 0.127
tigesV.aov <- aov(croissance$masse.tige.V ~ croissance$site)
summary(tigesV.aov)
#p = 0.0567

#test de normalité sur données transformées ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shapiro.test(resid(racinesV.aov))
#toujours pas
shapiro.test(resid(tigesV.aov))
#toujours pas

#égalité des variances ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bartlett.test(croissance$l.racines, croissance$site)
#variances ne sont pas homogènes
bartlett.test(croissance$masse.tige, croissance$site)
#variances ne sont pas hommogènes

#Tests post-hoc de Tukey ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Tuckey
d.racines.T <- TukeyHSD(d.racines.aov)
racines.T <- TukeyHSD(racines.aov)
tiges.T <- TukeyHSD(tiges.aov)

#attribuer les lettres
d.racines.cld <- multcompLetters4(d.racines.aov, d.racines.T)
print(d.racines.cld)

l.racines.cld <- multcompLetters4(racines.aov, racines.T)
print(l.racines.cld)

tiges.cld <- multcompLetters4(tiges.aov, tiges.T)
print(tiges.cld)


# Graphique longueur des racines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#tableau résumé avec lettres et 3e quartile
Tk.L <- group_by(croissance, site) %>%
  summarise(mean=mean(l.racines), quant = quantile(l.racines, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
l.racines.cld <- as.data.frame.list(l.racines.cld$`croissance$site`)
Tk.L$cld <- l.racines.cld$Letters
print(Tk.L)

#grahphique longueur
ggplot(croissance, aes(x=site, y=l.racines)) + 
  geom_boxplot(aes(fill = factor(..middle..)), show.legend = FALSE) +
  scale_fill_brewer(palette = "Purples") +
  xlab("Site") +
  ylab("Longueur totale des racines (m)") +
  geom_label(data = Tk.L, aes(x = site, y = quant, label = cld), 
             size = 3, nudge_x=0.15, nudge_y =0.07)



# Graphique masse des tiges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#tableau résumé avec lettres et 3e quartile
Tk.T <- group_by(croissance, site) %>%
  summarise(mean=mean(masse.tige), quant = quantile(masse.tige, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
m.tiges.cld <- as.data.frame.list(tiges.cld$`croissance$site`)
Tk.T$cld <- m.tiges.cld$Letters
print(Tk.T)

#grahphique longueur
ggplot(croissance, aes(x=site, y=masse.tige)) + 
  geom_boxplot(aes(fill = factor(..middle..)), show.legend = FALSE) +
  scale_fill_brewer(palette = "Blues") +
  xlab("Site") +
  ylab("Masse des tiges (g)") +
  geom_label(data = Tk.T, aes(x = site, y = quant, label = cld), 
             size = 3, nudge_x=0.18, nudge_y=0.008)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analyse factorielle multiple (AFM) des données de physicochimie et de colonisation racinaire ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(vegan)
library(FactoMineR)
#loader la fonction 'screestick'

#séparer les données en matrice réponse et matrices explicatives
lonicera.perfo <- lonicera[, 3:5] #s
lonicera.colo <- lonicera[, 6:10] #c
lonicera.pc <- lonicera[, 11:16] #s


#MFA
# MFA on 2 groups of variables: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Regroup the 2 tables
tab2 <- data.frame(lonicera.colo, lonicera.pc)
dim(tab2)
# Number of variables in each group
(grn <- c(ncol(lonicera.colo), ncol(lonicera.pc)))

# Compute the MFA without multiple plots
t2.mfa <- MFA(tab2,
              group = grn,
              type = c("c", "s"),
              ncp = 2,
              name.group = c("Colonisation", "Physicochimie"), graph = FALSE)

t2.mfa

# Plot the results
plot(t2.mfa,
     choix = "axes",
     habillage = "group",
     shadowtext = TRUE) 

plot(
  t2.mfa,
  choix = "ind",
  partial = "all",
  habillage = "group")

plot(t2.mfa,
     choix = "var",
     habillage = "group",
     shadowtext = TRUE) 

plot(t2.mfa, choix = "group")

# Eigenvalues, scree plot and broken stick model
ev <- t2.mfa$eig[, 1]
names(ev) <- paste("MFA", 1 : length(ev)) 
screestick(ev, las = 2)

# RV coefficients with tests (p-values above the diagonal of the matrix)
#lonicera.perfo, lonicera.colo, lonicera.pc, com.hel
rvp <- t2.mfa$group$RV
rvp[1, 2] <- coeffRV(lonicera.colo, scale(lonicera.pc))$p.value 
round(rvp[-4, -4], 4)


