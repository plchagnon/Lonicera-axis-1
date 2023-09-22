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
# GLM des données de performance des plants ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Cette analyse se fait avec le jeu de données de performance originel, et non
#avec les moyennes qui seront utilisées plus tard
#loader les données ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
croissance <- read.table("croissance.txt", header=TRUE)
croissance$site <- as.factor(croissance$site)
croissance$pot <- as.factor(croissance$pot)
croissance$inoculum <- as.factor(croissance$inoculum)

#réf : https://www.datacamp.com/tutorial/generalized-linear-models


library(energy)
library(multcomp)
library(multcompView)
library(datasets)
library(ggplot2)
library(dplyr)


#Modèle loi normale

x=rep(1:5,each=10)
y=3.4*.8*x+rnorm(50,0,.4)
mod=glm(y~x,family="gaussian")
r=residuals(mod)
hist(r)
plot(r~x)
abline(h=0,lty=2)

#Masse des tiges ~ sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#régression normale
tige.glm <- glm(masse.tige ~ site, data = croissance, family = "gaussian")

#vérifier les conditions
hist(residuals(tige.glm))
plot(residuals(tige.glm) ~ croissance$site)
abline(h=0,lty=2)
qqPlot(residuals(tige.glm))

# tests post-hoc
summary(tige.glm)
tige.T <- glht(tige.glm, mcp(site="Tukey"))
plot(tige.T)

#attribuer les lettres
tige.cld <- cld(tige.T, level=0.05)
print(tige.cld)


# Graphique masse des tiges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#tableau résumé avec lettres et 3e quartile
Tk.T <- group_by(croissance, site) %>%
  summarise(mean=mean(masse.tige), quant = quantile(masse.tige, probs = 0.75))

# extracting the compact letter display and adding to the Tk table
Tk.T$cld <- tige.cld$mcletters$Letters
print(Tk.T)

#grahphique longueur
ggplot(croissance, aes(x=site, y=masse.tige)) + 
  geom_boxplot(aes(fill = factor(..middle..)), show.legend = FALSE) +
  scale_fill_brewer(palette = "Greens") +
  xlab("Site") +
  ylab("Masse sèche des tiges (g)") +
  geom_label(data = Tk.T, aes(x = site, y = quant, label = cld), 
             size = 3, nudge_x=0.15, nudge_y =0.007)


#Longueur des racines ~ sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#régression normale
racine.glm <- glm(l.racines ~ site, data = croissance, family = "gaussian")

#vérifier les conditions
hist(residuals(racine.glm))

plot(residuals(racine.glm) ~ croissance$site)
abline(h=0,lty=2)

qqPlot(residuals(racine.glm))

# tests post-hoc
summary(racine.glm)
racine.T <- glht(racine.glm, mcp(site="Tukey"))
plot(racine.T)

#attribuer les lettres
racine.cld <- cld(racine.T, level=0.05)
print(racine.cld)

# Graphique longueur des racines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#tableau résumé avec lettres et 3e quartile
Tk.R <- group_by(croissance, site) %>%
  summarise(mean=mean(l.racines), quant = quantile(l.racines, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
racine.cld <- as.data.frame.list(racine.cld$`croissance$site`)
Tk.R$cld <- racine.cld$Letters
print(Tk.R)

#grahphique longueur
ggplot(croissance, aes(x=site, y=l.racines)) + 
  geom_boxplot(aes(fill = factor(..middle..)), show.legend = FALSE) +
  scale_fill_brewer(palette = "Purples") +
  xlab("Site") +
  ylab("Longueur totale des racines (m)") +
  geom_label(data = Tk.R, aes(x = site, y = quant, label = cld), 
             size = 3, nudge_x=0.15, nudge_y =0.07)

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


