# Analyses axe 1
#version finale

#loader les données
setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")
#comm : métacommunauté (rangées = sites, colonnes = champignons) pour 31 échantillons de sol
#taxa : taxonomie (autant de rangées que de champignons, et rangs taxonomiques en colonnes)
#ASV : abondance de 7 espèces de champignons jugées pertinentes par régression elastic net, pour 32 échantillons de sol
#lonicera : variables de performance, de colonisation racinaire et de physicochimie, pour 49 échantillons de sol
#croissance : variables de performance, pour 314 pots de Lonicera
#perfo.comm : variables de performance, pour les 32 échantillons ayant des données fongiques
#lonicera.comm : variables de performance, de colonisation racinaire, de physicochimie, et abondance des ASV de champignons pertinentes, pour les 31 échantillons ayant des données fongiques


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GLM des données de performance des plants ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#loader les données ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#réf : https://www.datacamp.com/tutorial/generalized-linear-models

library(energy)
library(multcomp)
library(multcompView)
library(datasets)
library(ggplot2)
library(dplyr)
library(car)


#Créer un modèle de loi normale
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
# a l'air normal

# tests post-hoc
summary(tige.glm)
#valeur p du site 1: ??
#valeur p du site 2: 0.05894
#valeur p du site 3: 0.00536 ** 
#valeur p du site 4: 0.00634 **
#valeur p du site 5: 0.16636

tige.T <- glht(tige.glm, mcp(site="Tukey"))
plot(tige.T)

#attribuer les lettres
tige.cld <- cld(tige.T, level=0.05)
print(tige.cld)


# Graphique masse des tiges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#tableau résumé avec lettres et 3e quartile
Tk.T <- group_by(croissance, site) %>%
  summarise(mean=mean(masse.tige), quant = quantile(masse.tige, probs = 0.75)) %>%
  arrange(desc(mean))

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
Tk.R$cld <- racine.cld$mcletters$Letters
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
# GLM des données de physicochimie ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#loader les données ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")


library(energy)
library(multcomp)
library(multcompView)
library(datasets)
library(ggplot2)
library(dplyr)
library(car)


#Créer un modèle de loi normale
x=rep(1:5,each=10)
y=3.4*.8*x+rnorm(50,0,.4)
mod=glm(y~x,family="gaussian")
r=residuals(mod)
hist(r)
plot(r~x)
abline(h=0,lty=2)

#Phosphore ~ sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#régression normale
phos.glm <- glm(phosphore ~ site, data = lonicera, family = "gaussian")

#vérifier les conditions
hist(residuals(phos.glm))
plot(residuals(phos.glm) ~ lonicera$site)
abline(h=0,lty=2)
qqPlot(residuals(phos.glm))
# pas tout-à-fait normal

#transformation log
lonicera$phosphore.log <- log(lonicera$phosphore)
#régression racine carrée
phos.glm.l <- glm(phosphore.log ~ site, data = lonicera, family = "gaussian")

#vérifier les conditions
hist(residuals(phos.glm.r))
plot(residuals(phos.glm.r) ~ lonicera$site)
abline(h=0,lty=2)
qqPlot(residuals(phos.glm.r))
# pas tout-à-fait normal
# tests post-hoc
summary(phos.glm.r)
#valeur p du site 1: ??
#valeur p du site 2: 4.87e-06 ***
#valeur p du site 3: 2.30e-06 ***
#valeur p du site 4: 7.24e-08 ***
#valeur p du site 5: 6.22e-07 ***

tige.T <- glht(tige.glm, mcp(site="Tukey"))
plot(tige.T)

#attribuer les lettres
tige.cld <- cld(tige.T, level=0.05)
print(tige.cld)


# Graphique masse des tiges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#tableau résumé avec lettres et 3e quartile
Tk.T <- group_by(croissance, site) %>%
  summarise(mean=mean(masse.tige), quant = quantile(masse.tige, probs = 0.75)) %>%
  arrange(desc(mean))

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analyse factorielle multiple (AFM) des données de physicochimie et de colonisation racinaire ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#loader les données
setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")

library(vegan)
library(FactoMineR)
#loader la fonction 'screestick'

#séparer les données en matrice réponse et matrices explicatives
lonicera.perfo <- lonicera[, 3:5] #s
lonicera.colo <- lonicera[, 6:10] #c
lonicera.pc <- lonicera[, 11:16] #s
lonicera.seq <- lonicera.comm[17:23]

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#perMANOVA ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ref: https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/

library(vegan)

#importer les données 
setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")

#données de site
Sites <- rep(c(1,2,3,4,5),times=c(10,10,10,9,10))

lonicera.perfo <- lonicera[, 3:5] # unités diff
lonicera.colo <- lonicera[, 6:10] # mêmes unités ; abondances
lonicera.pc <- lonicera[, 11:16]  # unités diff


#perMANOVA des données de physico-chimie ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Verify multivariate homogeneity of within-group covariance
pc.d1 <- dist(lonicera.pc)
# Test of homogeneity of within-cell dispersions
(pc.cell.MHV <- betadisper(pc.d1, Sites)) 
permutest(pc.cell.MHV) #valeur p est de 0.006 < 0.05

#permanova avec adonis2
pc.adonis.scale <- adonis2(scale(lonicera.pc) ~ Sites, 
                           method = "euc",
                           by = "term")
# p = 0.001 ***

#graphique - cadrage de type 1
pc.rda <- rda(lonicera.pc ~ Sites)
plot(pc.rda,
     scaling = 1,
     display = "wa",
     main = "MANOVA permutationnelle, cadrage 1 - scores 'wa'",
     xlim=c(-2.4, 1.5),
     ylim=c(-1.2, 2.2))

ordispider(pc.rda, Sites, scaling = 1,
           label = TRUE,
           col = "blue"
)

pc.sc1 <-
  scores(pc.rda, scaling = 1,
         display = "species")
arrows(0, 0, pc.sc1[, 1] * 0.3, pc.sc1[, 2] * 0.3, length = 0.1,
       angle = 10,
       col = "red"
) 
text(
  pc.sc1[, 1] * 0.3, pc.sc1[, 2] * 0.3,
  labels = rownames(pc.sc1), pos = 4,
  cex = 0.8,
  col = "red"
)

#graphique - cadrage de type 2
plot(pc.rda,
     scaling = 2,
     display = "wa",
     main = "MANOVA permutationnelle, cadrage 2 - scores 'wa'",
     ylim=c(-2, 3))

ordispider(pc.rda, Sites, scaling = 2,
           label = TRUE,
           col = "blue"
)

pc.sc2 <-
  scores(pc.rda, scaling = 2,
         display = "species")
arrows(0, 0, pc.sc2[, 1] * 0.3, pc.sc2[, 2] * 0.3, length = 0.1,
       angle = 10,
       col = "red"
) 
text(
  pc.sc2[, 1] * 0.3, pc.sc2[, 2] * 0.3,
  labels = rownames(pc.sc2), pos = 4,
  cex = 0.8,
  col = "red"
)


#perMANOVA des données de colonisation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Verify multivariate homogeneity of within-group covariance
colo.d1 <- dist(lonicera.colo)
# Test of homogeneity of within-cell dispersions
(colo.cell.MHV <- betadisper(colo.d1, Sites)) 
permutest(colo.cell.MHV) #valeur p est de 0.396 > 0.05

colo.adonis.scale <- adonis2(scale(lonicera.colo) ~ Sites, 
                             method = "euc",
                             by = "term")
# p = 0.001 ***

#graphique - cadrage de type 1

colo.rda <- rda(lonicera.colo ~ Sites)

plot(colo.rda,
     scaling = 1,
     display = "wa",
     main = "MANOVA permutationnelle, cadrage 1 - scores 'wa'",
     xlim = c(-3.5, 2.5),  
     ylim = c(-2.2, 4.2))

ordispider(colo.rda, Sites, scaling = 1,
           label = TRUE,
           col = "blue"
)

colo.sc1 <-
  scores(colo.rda, scaling = 1,
         display = "species")
arrows(0, 0, colo.sc1[, 1] * 0.3, colo.sc1[, 2] * 0.3, length = 0.1,
       angle = 10,
       col = "red"
) 
text(
  colo.sc1[, 1] * 0.3, colo.sc1[, 2] * 0.3,
  labels = rownames(colo.sc1), pos = 4,
  cex = 0.8,
  col = "red"
)

#graphique - cadrage de type 2
plot(colo.rda,
     scaling = 2,
     display = "wa",
     main = "MANOVA permutationnelle, cadrage 2 - scores 'wa'")

ordispider(colo.rda, Sites, scaling = 2,
           label = TRUE,
           col = "blue"
)

colo.sc2 <-
  scores(colo.rda, scaling = 2,
         display = "species")
arrows(0, 0, colo.sc2[, 1] * 0.3, colo.sc2[, 2] * 0.3, length = 0.1,
       angle = 10,
       col = "red"
) 
text(
  colo.sc2[, 1] * 0.3, colo.sc2[, 2] * 0.3,
  labels = rownames(colo.sc2), pos = 4,
  cex = 0.8,
  col = "red"
)



#perMANOVA des données de performance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Verify multivariate homogeneity of within-group covariance
perf.d1 <- dist(lonicera.perfo)
# Test of homogeneity of within-cell dispersions
(perf.cell.MHV <- betadisper(perf.d1, Sites)) 
permutest(perf.cell.MHV) #valeur p est de 0.085 > 0.05

#permanova avec adonis2
perf.adonis.scale <- adonis2(scale(lonicera.perfo) ~ Sites, 
                             method = "euc",
                             by = "term")
# p = 0.114

#pas significatif
#pas de graphiques (ne marchent pas)

#perMANOVA des données de séquençage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sites.comm <- rep(c(1,2,3,4,5),times=c(6,7,8,4,7))

# Verify multivariate homogeneity of within-group covariance
seq.d1 <- dist(lonicera.seq)
# Test of homogeneity of within-cell dispersions
(seq.cell.MHV <- betadisper(seq.d1, Sites.comm)) 
permutest(seq.cell.MHV) #valeur p est de 0.551 > 0.05

#permanova avec adonis 2
seq.adonis.scale <- adonis2(scale(lonicera.seq) ~ Sites.comm, 
                            method = "euc",
                            by = "term")
# p = 0.119
#PAS significatif
#graphiques : il faudrait enlever les sites 49, 23 et possiblement 50

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Varpart ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(vegan)
library(eulerr)
library(remotes)
library(colorRamps)
library(RColorBrewer)

#importer les données 
setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")


lonicera.perfo <- lonicera[, 3:4] # unités diff
lonicera.colo <- lonicera[, 6:10] # mêmes unités ; abondances
lonicera.pc <- lonicera[, 11:16]  # unités diff

loni.perfo <- lonicera.comm[, 3:4]
loni.colo <- lonicera.comm[, 6:10]
loni.pc <- lonicera.comm[, 11:16]
loni.seq <- lonicera.comm[17:23]

showvarparts(3, bg = c("red", "blue", "yellow"))

# Biomasse en fonction de la colonisation et physicochimie ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(lonicera.part.2 <- varpart(lonicera.perfo$l.racines, lonicera.colo, lonicera.pc)) 

plot(lonicera.part.2, 
     digits = 2, 
     Xnames = c("Colonisation", "Physicochimie"),
     bg = c("green", "yellow"))

#meilleur diagramme

venn.cp <- euler(c("Colonisation" = 0,
                "Physicochimie" = lonicera.part.2$part$fract$Adj.R.squared[2],
                "Residuals" = lonicera.part.2$part$indfract$Adj.R.squared[4],
                "Colonisation&Physicochimie" = lonicera.part.2$part$fract$Adj.R.squared[3]))
                
plot(venn.cp)


# performance en fonction de la physicochimie et des champignons du sol
(lonicera.part.seq <- varpart(loni.perfo, loni.seq, loni.pc)) 
plot(lonicera.part.seq, 
     digits = 2, 
     Xnames = c("Physicochimie", "Champignons"),
     bg = c("yellow", "pink"))

#meilleur diagramme

venn.chp <- euler(c("Physicochimie" = lonicera.part.seq$part$fract$Adj.R.squared[1],
                   "Champignons" = 0,
                   "Residuals" = lonicera.part.seq$part$indfract$Adj.R.squared[4],
                   "Physicochimie&Champignons" = 0))

plot(venn.chp)

# biomasse en fonction de la p-c, colonisation et les espèces de champignon sélectionnées ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(loni.part.all <- varpart(loni.perfo$masse.tige, loni.pc, loni.colo, loni.seq)) 

plot(loni.part.all, 
     digits = 2, 
     Xnames = c("Physicochimie", "Colonisation", "Champignons"),
     bg = c("yellow", "green", "pink"))

#meilleur diagramme

fit3 <- c("Physicochimie" = 0,
          "Colonisation" = 0,
          "Champignons" = loni.part.all$part$indfract$Adj.R.square[3],
          "Residuals" = loni.part.all$part$indfract$Adj.R.square[8],
          "Physicochimie&Colonisation" = loni.part.all$part$indfract$Adj.R.square[4],
          "Physicochimie&Champignons" = 0,
          "Colonisation&Champignons" = loni.part.all$part$indfract$Adj.R.square[5],
          "Colonisation&Champignons&Physicochimie" = 0)

venn3 <- euler(fit3, shape = "ellipse")
plot(venn, shape = "ellipse")

plot(venn3,
     fills = c("dodgerblue4", "darkgoldenrod1", "cornsilk4", "coral1"),
     edges = FALSE,
     fontsize = 8,
     quantities = list(fontsize = 8))

plot(venn3,
     fills = brewer.pal(4, "Spectral"),
     edges = FALSE,
     fontsize = 8,
     quantities = list(fontsize = 8))


plot(venn,
     fills = blue2red(4),
     edges = FALSE,
     fontsize = 8,
     quantities = list(fontsize = 8))

#Longueur des racines en fonction des 3 autres ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(loni.part.all.r <- varpart(loni.perfo$l.racines, loni.pc, loni.colo, loni.seq)) 

plot(loni.part.all.r, 
     digits = 2, 
     Xnames = c("Physicochimie", "Colonisation", "Champignons"),
     bg = c("yellow", "green", "pink"))

#meilleur diagramme

venn <- euler(c("Physicochimie" = 0,
                "Colonisation" = 0,
                "Champignons" = loni.part.all.r$part$fract$Adj.R.square[3],
                "Residuals" = loni.part.all.r$part$indfract$Adj.R.square[8],
                "Physicochimie&Colonisation" = 0,
                "Physicochimie&Champignons" = 0,
                "Colonisation&Champignons" = 0,
                "Colonisation&Champignons&Physicochimie" = loni.part.all.r$part$indfract$Adj.R.square[7]
))


venn <- euler(c("X1" = model$part$fract$Adj.R.square[2],
                "X2" = model$part$fract$Adj.R.square[2],
                "X3" = model$part$fract$Adj.R.square[3],
                "Residuals" = model$part$indfract$Adj.R.square[8],
                "X1&X2" = model$part$indfract$Adj.R.square[4],
                "X1&X3" = model$part$indfract$Adj.R.square[6],
                "X2&X3" = model$part$indfract$Adj.R.square[5],
                "X2&X3&X1" = model$part$indfract$Adj.R.square[7]
))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLS reg ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Ref: https://www.statology.org/partial-least-squares-in-r/
#Ref2 : https://support.minitab.com/en-us/minitab/21/help-and-how-to/statistical-modeling/regression/how-to/partial-least-squares/interpret-the-results/key-results/
#Voir exemple dans "Exemple_test.R"

library("pls")

#importer les données 
setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")



# Variable réponse : masse des tiges avec séquençage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit PCR model (Fit Partial Least Squares Model) ~~~~~~~~~~~~~~~~~~~~~~~
#scale : centrer-réduire
#validation=”CV”: This tells R to use k-fold cross-validation to evaluate the performance of the model. 
#Note that this uses k=10 folds by default. Also note that you can specify “LOOCV” instead to perform leave-one-out cross-validation.
lonicera.comm.reduit <- lonicera.comm[,-c(1:2, 4:5)]
pls.tige <- plsr(masse.tige~., data=lonicera.comm.reduit, scale=TRUE, validation="CV")

#Choose the Number of PLS Components ~~~~~~~~~~~~~~~~~~~~~~~
#by looking at the test root mean squared error (test RMSE) calculated by the k-fold cross-validation
#view summary of model fitting
summary(pls.tige)


#1) table : VALIDATION: RMSEP
#This table tells us the test RMSE calculated by the k-fold cross validation. We can see the following:

#If we only use the intercept term in the model, the test RMSE is 0.01647
#If we add in the first PLS component, the test RMSE goes up to 0.01815
#If we add in the second PLS component, the test RMSE goes up to 0.01878
#We can see that adding additional PLS components also increases the test RMSE
# Je ne sais pas combien de composants sont optimaux

#2) Table : TRAINING: % variance explained
#percentage of the variance in the response variable explained by the PLS components

#Using just the first PLS component :  explain X% of the variation in the response variable
#adding in the second PLS component, we can explain X% of the variation in the response variable.
##adding in the second PLS component, we can explain X% of the variation in the response variable.
#Note that we’ll always be able to explain more variance by using more PLS components, but we can see that adding in more than two PLS components doesn’t actually increase the percentage of explained variance by much.


#visualize cross-validation plots (the test RMSE)
validationplot(pls.tige)
validationplot(pls.tige, val.type="MSEP")
validationplot(pls.tige, val.type="R2")


#Use the Final Model to Make Predictions ~~~~~~~~~~~~~~~~~~~~~~~

#use the final model with two PLS components to make predictions on new observations.
#split the original dataset into a training and testing set and use the final model with two PLS components to make predictions on the testing set

#define training and testing sets
train.pls.tige <- lonicera.comm.reduit[1:20, c("masse.tige", "arb", "coils", "ves", "nonmyc", "dse", "pH", "ammonium", "nitrate", "phosphore", "wsa", "mo", "Calonectria_sp.", "Rhizophaguss_irregularis", "Rhizophaguss_intraradices", "Claroideoglomu_sp.", "Rhizophaguss_intraradices", "Podilas_humilis", "Minutisphaeras_aspera")]
y_test.pls.tige <- lonicera.comm.reduit[20:nrow(lonicera.comm.reduit), c("masse.tige")]
test.pls.tige <- lonicera.comm.reduit[20:nrow(lonicera.comm.reduit), c("arb", "coils", "ves", "nonmyc", "dse", "pH", "ammonium", "nitrate", "phosphore", "wsa", "mo", "Calonectria_sp.", "Rhizophaguss_irregularis", "Rhizophaguss_intraradices", "Claroideoglomu_sp.", "Rhizophaguss_intraradices", "Podilas_humilis", "Minutisphaeras_aspera")]

#use model to make predictions on a test set
pls.tige.2 <- plsr(masse.tige~., data=train.pls.tige, scale=TRUE, validation="CV")
pcr_pred.pls <- predict(pls.tige.2, test.pls.tige, ncomp=2)

#calculate RMSE
sqrt(mean((pcr_pred.pls - y_test.pls.tige)^2))

#[1] 0.02515554



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(vegan)
library(FactoMineR)
library(ggplot2)

#importer les données 
setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")

lonicera.perfo <- lonicera[, 3:4] # unités diff
lonicera.colo <- lonicera[, 6:10] # mêmes unités ; abondances
lonicera.pc <- lonicera[, 11:16]  # unités diff

loni.perfo <- lonicera.comm[, 3:4]
loni.colo <- lonicera.comm[, 6:10]
loni.pc <- lonicera.comm[, 11:16]
loni.seq <- lonicera.comm[17:23]

#Échantillons par sites
Sites <- rep(c(1,2,3,4,5),times=c(10,10,10,9,10))
Sites.ellipse <- c(1, 2, 3, 4, 5)


# Colonisation des racines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

colo.pca <- rda(lonicera.colo, scale = TRUE)
summary(colo.pca, axes = 0)
biplot(colo.pca, scaling = 1)
biplot(colo.pca, scaling = 2)

#plus facile à visionner
PCA(lonicera.colo, scale=TRUE)

#en couleur
myc.site <- colo.pca$CA$u[, -c(3, 4, 5)]
myc.site <- cbind(myc.site, Sites)

plot(colo.pca, type="none", xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4))
points(myc.site, col=myc.site[,3])
ordiellipse(myc.site, groups = myc.site[,3], col=Sites.ellipse, label=TRUE)

#Graphique avec ggplot2
myc.site <- as.data.frame(myc.site)
myc.site$Sites <- as.factor(myc.site$Sites)

ggplot(data=myc.site, aes(x=PC1, y=PC2, colour = Sites))+
  geom_point()+
  xlab("Premier axe en compostantes principales")+
  ylab("Deuxième axe en compostantes principales")+
  stat_ellipse(type = "t", lwd = 0.4, level = 0.75)

# Paramètres du sol ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sol.pca <- rda(lonicera.pc, scale=TRUE)
summary(sol.pca, axes=0)
biplot(sol.pca, scaling = 1)
biplot(sol.pca, scaling = 2)

#plus facile à visualiser
PCA(lonicera.pc, scale=TRUE)

#en couleur
sol.site <- sol.pca$CA$u[, -c(3, 4, 5, 6)]
sol.site <- cbind(sol.site, Sites)

plot(sol.pca, type="none", xlim = c(-0.4, 0.4), ylim = c(-0.45, 0.45))
points(sol.site, col=sol.site[,3])
ordiellipse(sol.site, groups = sol.site[,3], col=Sites.ellipse, label=TRUE)

#Graphique avec ggplot2
sol.site <- as.data.frame(sol.site)
sol.site$Sites <- as.factor(sol.site$Sites)

ggplot(data=sol.site, aes(x=PC1, y=PC2, colour = Sites))+
  geom_point()+
  xlab("Premier axe en compostantes principales")+
  ylab("Deuxième axe en compostantes principales")+
  stat_ellipse(type = "t", lwd = 0.4, level = 0.75)

#Performance des plants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PAS utilisé ecm

plante.pca <- rda(lonicera.perfo, scale = TRUE)
summary(plante.pca, axes=0)
biplot(plante.pca, scaling = 1)
biplot(plante.pca, scaling = 2)

#plus facile à visionner
PCA(lonicera.pl, scale=TRUE)

#en couleur
plante.site <- plante.pca$CA$u[, -c(3)]
plante.site <- cbind(plante.site, Sites)

plot(plante.pca, type="none", xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4))
points(plante.site, col=plante.site[,3])
ordiellipse(plante.site, groups = plante.site[,3], col=Sites.ellipse, label=TRUE)

#Graphique avec ggplot2
plante.site <- as.data.frame(plante.site)
plante.site$Sites <- as.factor(plante.site$Sites)


ggplot(data=myc.site, aes(x=PC1, y=PC2, colour = Sites))+
  geom_point()+
  xlab("Premier axe en compostantes principales")+
  ylab("Deuxième axe en compostantes principales")+
  stat_ellipse(type = "t", lwd = 0.4, level = 0.75)


#Champignons du sol ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sites.comm <- rep(c(1,2,3,4,5),times=c(6,7,8,4,7))

#transfo 1
loni.seq.hel <- decostand(loni.seq, "hellinger")

#transfo 2
tr4=function(x){(x/sum(x))^(1/4)}
loni.seq.hel <- t(apply(loni.seq,1,tr4))

seq.pca <- rda(loni.seq.hel, scale = TRUE)
summary(seq.pca, axes=0)
biplot(seq.pca, scaling = 1)
biplot(seq.pca, scaling = 2)

#plus facile à visionner
PCA(loni.seq.hel, scale=TRUE)

#en couleur
seq.site <- seq.pca$CA$u[, -c(3:7)]
seq.site <- cbind(seq.site, Sites.comm)

plot(seq.pca, type="none", xlim = c(-0.6, 0.4), ylim = c(-0.5, 0.3))
points(seq.site, col=seq.site[,3])
ordiellipse(seq.site, groups = seq.site[,3], col=Sites.ellipse, label=TRUE)

#Graphique avec ggplot2
seq.site <- as.data.frame(seq.site)
seq.site$Sites.comm <- as.factor(seq.site$Sites.comm)


ggplot(data=seq.site, aes(x=PC1, y=PC2, colour = Sites.comm))+
  geom_point()+
  xlab("Premier axe en compostantes principales")+
  ylab("Deuxième axe en compostantes principales")+
  stat_ellipse(type = "t", lwd = 0.4, level = 0.75)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Régession linéaire mixte ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(lme4)
library(base)

#importer les données 
setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")


#Scale
lonicera.comm.s <- scale(lonicera.comm[, 3:23], center = TRUE, scale = TRUE)
lonicera.comm.s <- as.data.frame(lonicera.comm.s)
lonicera.comm.s <- cbind(lonicera.comm$site, lonicera.comm.s)
names(lonicera.comm.s)[names(lonicera.comm.s) == 'lonicera.comm$site'] <- 'site'

#Avec glmer de lme4
#version sans effet de site
model.lmer <- lmer(masse.tige ~ arb + coils + ves + nonmyc + dse + pH + ammonium + nitrate + phosphore + wsa + mo + Calonectria_sp. + Rhizophaguss_irregularis + Rhizophaguss_intraradices + Claroideoglomu_sp. + Rhizophaguss_intraradices + Podilas_humilis + Minutisphaeras_aspera + (1|site), data = lonicera.comm.s)

#summary
summary(model.lmer)
print(model.lmer, correlation=TRUE)
plot(model.lmer)

#version avec effet de site
model.lmer <- lmer(masse.tige ~ arb + coils + ves + nonmyc + dse + pH + ammonium + nitrate + phosphore + wsa + mo + Calonectria_sp. + Rhizophaguss_irregularis + Rhizophaguss_intraradices + Claroideoglomu_sp. + Rhizophaguss_intraradices + Podilas_humilis + Minutisphaeras_aspera, data = lonicera.comm.s)

#summary
summary(model.lmer)
print(model.lmer, correlation=TRUE)
plot(model.lmer)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Régession globale ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(stats)

#importer les données 
setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")

#régression globale en fonction de toutes les variables
test.glm <- glm(masse.tige ~ arb + coils + ves + nonmyc + dse + pH + ammonium + nitrate + phosphore + wsa + mo + Calonectria_sp. + Rhizophaguss_irregularis + Rhizophaguss_intraradices + Claroideoglomu_sp. + Rhizophaguss_intraradices + Podilas_humilis + Minutisphaeras_aspera, data = lonicera.comm, family = "gaussian")
summary(test.glm)

newdata = data.frame(disp=200, hp= 100)
predict(test.glm, type="response")

#régression globale sans les variables fongiques
test.glm2 <- glm(masse.tige ~ arb + coils + ves + nonmyc + dse + pH + ammonium + nitrate + phosphore + wsa + mo, data = lonicera, family = "gaussian")
summary(test.glm2)

#régression avec juste une variable
test.glm3 <- glm(masse.tige ~ phosphore, data = lonicera, family = "gaussian")
summary(test.glm3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Histogramme ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(RColorBrewer)

#importer les données 
setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")

#créer les barres d'erreur
ecart <- NULL
for(t in 1:50)
{ecart <- append(ecart, sd(croissance$masse.tige[croissance$inoculum == t]))}
ecart <- ecart[-40]
cbind(lonicera, ecart)
                     
#ggplot : avec les moyennes
ggplot(lonicera, aes(x = site, y = masse.tige, color=site))+
  geom_col(position=position_dodge2(width = 0.9), show.legend = FALSE)+
  geom_errorbar(aes(ymin=masse.tige-ecart, ymax=masse.tige+ecart), width=0.4, position =  position_dodge2(width = 0.9), show.legend = FALSE)+
  xlab("Sites")+
  ylab("Masse des tiges (g)")

# à faire : barres d'erreur et changer les couleurs

ggplot(croissance, aes(x = inoculum, y = masse.tige, color=site))+
  geom_col(position="dodge", show.legend = FALSE)

ggplot(lonicera, aes(x = site, y = masse.tige, fill=inoculum))+
  geom_bar(position="dodge", stat = "identity", show.legend = FALSE)+
  geom_errorbar(aes(ymin=masse.tige-ecart, ymax=masse.tige+ecart), width=0.4, position =  position_dodge2(width = 0.9), show.legend = FALSE)+
  xlab("Sites")+
  ylab("Masse des tiges (g)")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section tests ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#VarPart sur tout
a = -0.00557
b = -0.05649
c = 0.33373
d = 0.06491
e =  0.32944
f = -0.02993
g = -0.16380

a+d+f+g
b+d+e+g
c+e+f+g


#Varpart sur les données d'abondance relative
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Autre.Rdata")

# biomasse en fonction de la p-c, colonisation et les espèces de champignon sélectionnées ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(loni.part.rel <- varpart(loni.perfo$masse.tige, loni.pc, loni.colo, ASV.6)) 

plot(loni.part.rel, 
     digits = 2, 
     Xnames = c("Physicochimie", "Colonisation", "Champignons"),
     bg = c("yellow", "green", "pink"))

#meilleur diagramme
venn <- euler(c("Physicochimie" = loni.part.rel$part$indfract$Adj.R.square[1],
                "Colonisation" = 0,
                "Champignons" = loni.part.rel$part$indfract$Adj.R.square[3],
                "Residuals" = loni.part.rel$part$indfract$Adj.R.square[8],
                "Physicochimie&Colonisation" = loni.part.rel$part$indfract$Adj.R.square[4],
                "Physicochimie&Champignons" = 0,
                "Colonisation&Champignons" = loni.part.rel$part$indfract$Adj.R.square[5],
                "Colonisation&Champignons&Physicochimie" = 0
))

plot(venn, shape = "ellipse")

plot(venn,
     fills = c("dodgerblue4", "darkgoldenrod1", "cornsilk4", "coral1"),
     edges = FALSE,
     fontsize = 8,
     quantities = list(fontsize = 8))

plot(venn,
     fills = brewer.pal(4, "Spectral"),
     edges = FALSE,
     fontsize = 8,
     quantities = list(fontsize = 8))

