## Analyses Lonicera

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Préparation des données ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#importer les données 
setwd("~/Documents/Maîtrise/Données Axe 1")
lonicera <- read.table("Lonicera.txt", header=TRUE)
lonicera$site <- as.factor(lonicera$site)
lonicera$inoculum <- as.factor(lonicera$inoculum)
#enlever le dernière ligne (stérile)
lonicera <- lonicera[-50,]

#séparer les données en matrice réponse et matrices explicatives
lonicera.perfo <- lonicera[, 3:5] # unités diff
lonicera.colo <- lonicera[, 6:10] # mêmes unités ; abondances
lonicera.pc <- lonicera[, 11:16]  # unités diff

#données de site
Sites <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#perMANOVA ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ref: https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/

library(vegan)

#importer les données 
setwd("~/Documents/Maîtrise/Données Axe 1")
lonicera <- read.table("Lonicera.txt", header=TRUE)
lonicera$site <- as.factor(lonicera$site)
lonicera$inoculum <- as.factor(lonicera$inoculum)
#enlever le dernière ligne (stérile)
lonicera <- lonicera[-50,]

#données de site
Sites <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5)

lonicera.perfo <- lonicera[, 3:5] # unités diff
lonicera.colo <- lonicera[, 6:10] # mêmes unités ; abondances
lonicera.pc <- lonicera[, 11:16]  # unités diff
rownames(lonicera.pc) <- c(1:39, 41:50)
colnames(lonicera.pc) <- c(1:31)

# Verify multivariate homogeneity of within-group covariance
pc.d1 <- dist(lonicera.pc)
# Test of homogeneity of within-cell dispersions
(pc.cell.MHV <- betadisper(pc.d1, Sites)) 
permutest(pc.cell.MHV) #valeur p est de 0.006 < 0.05


#perMANOVA
pc.adonis.scale <- adonis2(scale(lonicera.pc) ~ Sites, 
                     method = "euc",
                     by = "term")
# p = 0.001 ***

#graphique - cadrage de type 1
pc.rda <- rda(lonicera.pc ~ Sites)
plot(pc.rda,
     scaling = 1,
     display = "wa",
     main = "MANOVA permutationnelle, cadrage 1 - scores 'wa'")

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
     main = "MANOVA permutationnelle, cadrage 2 - scores 'wa'")

ordispider(pc.rda, Sites, scaling = 2,
           label = TRUE,
           col = "blue"
)

pc.sc1 <-
  scores(pc.rda, scaling = 2,
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


#Exemple écologie numérique ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#données
load("/Users/coralie/Documents/Maîtrise/Analyse quantitative des données/En R/NEwR-2ed_code_data/NEwR2-Data/Doubs.RData")
#espèces de poissons
head(spe)

#matrice de distance de Hellinger
spe.hel <- decostand(spe, "hellinger")

#création des facteurs
# Creation of a factor 'elevation' (3 levels, 9 sites each)
ele.fac <- gl(3, 9, labels = c("high", "mid", "low")) 
# Creation of a factor mimicking 'pH'
pH.fac <-
  as.factor(c(1, 2, 3, 2, 3, 1, 3, 2, 1, 2, 1, 3, 3, 2, 1, 1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 3))
# Is the two-way factorial design balanced?
table(ele.fac, pH.fac)

#perMANOVA

perma.ex <- adonis2(spe.hel[1:27, ] ~ ele.fac * pH.fac, 
        method = "euc",
        by = "term"
)
#l'élévation est significatif, mais pas le pH ou l'interaction

#plot?
plot(perma.ex)
#dear god


#Avec mes données : longueur des racines et sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~
#(ne fonctionne pas)

setwd("~/Documents/Maîtrise/Données Axe 1")

lonicera <- read.table("Lonicera.txt", header=TRUE)
lonicera$site <- as.factor(lonicera$site)
lonicera$inoculum <- as.factor(lonicera$inoculum)
#enlever le dernière ligne (stérile)
lonicera <- lonicera[-50,]

#Préparer les données
racines.hel <- decostand(lonicera$l.racines, "hellinger")

#perMANOVA

perma.racines <- adonis2(racines.hel ~ lonicera$pH * lonicera$phosphore, 
                    method = "euc",
                    by = "term"
)
perma.racines

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Régression Elastic Net ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#réf: https://hackernoon.com/an-introduction-to-ridge-lasso-and-elastic-net-regression-cca60b4b934f

# also allows us to tune the alpha parameter where 𝞪 = 0 corresponds to ridge and 𝞪 = 1 to lasso. 
#Simply put, if you plug in 0 for alpha, the penalty function reduces to the L1 (ridge) term and if we set alpha to 1 we get the L2 (lasso) term. Therefore we can choose an alpha value between 0 and 1 to optimize the elastic net. 
#Effectively this will shrink some coefficients and set some to 0 for sparse selection.

#Exemple de hackenoon ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Praparing the data ~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(caret)
library(glmnet)

data("Boston", package = "MASS")

#set a seed so you can reproduce the results
set.seed(1212)

#split the data into training and test data
sample_size <- floor(0.75 * nrow(Boston))
training_index <- sample(seq_len(nrow(Boston)), size = sample_size)
train <- Boston[training_index, ]
test <- Boston[-training_index, ]

# We also should create two objects to store predictor (x) and response variables (y, median value)
Predictorx <- model.matrix(medv~., train)[,-1]
Responsey <- train$medv

# Performing Elastic Net regression ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#tune parameters to identify the best alpha and lambda values
#We will tune the model by iterating over a number of alpha and lambda pairs and we can see which pair has the lowest associated error

model.net <- train(medv ~., data = train, method = "glmnet",trControl = trainControl("cv", number = 10),tuneLength = 10)
model.net$bestTune
# alpha     lambda
#5   0.1 0.08468587

coef(model.net$finalModel, model.net$bestTune$lambda)
x.test.net <- model.matrix(medv ~., test)[,-1]
predictions.net <- model.net %>% predict(x.test.net)

data.frame(RMSE.net = RMSE(predictions.net, test$medv),Rsquare.net = R2(predictions.net, test$medv))
# RMSE = 5.240535 and R² = 0.768797



# Avec mes données ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Importer les données ~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(caret)
library(glmnet)
library(prospectr)

setwd("~/Documents/Maîtrise/Données Axe 1")

lonicera <- read.table("Lonicera.txt", header=TRUE)
lonicera$site <- as.factor(lonicera$site)
lonicera$inoculum <- as.factor(lonicera$inoculum)
#enlever le dernière ligne (stérile)
lonicera <- lonicera[-50,]

# Préparer les données ~~~~~~~~~~~~~~~~~~~~~~~

#split the data into training and test data
sample_size <- floor(0.75 * nrow(lonicera))
lonicera.num <- data.matrix(lonicera)
training_index <- kenStone(lonicera.num, k=sample_size, metric="euclid")
train <- lonicera[training_index$model, ]
test <- lonicera[training_index$test, ]

## Première régression elastic net : en fonction de la masse des tiges ~~~~~~~~~~~~~~~~~~~~~~~

# Create two objects to store predictor (x) and response variables (y, median value)
Predictor.x <- model.matrix(masse.tige~ arb + coils + ves + nonmyc + dse + pH + ammonium + nitrate + phosphore, train)[,-1]
Response.y <- train$masse.tige


#tune parameters to identify the best alpha and lambda values
#We will tune the model by iterating over a number of alpha and lambda pairs and we can see which pair has the lowest associated error
model.net <- train(masse.tige~ arb + coils + ves + nonmyc + dse + pH + ammonium + nitrate + phosphore, data = train, method = "glmnet",trControl = trainControl("cv", number = 10),tuneLength = 10)
model.net$bestTune
# alpha     lambda
# 63   0.7  0.005558789
#C'est différent chaque fois
#12JN2023: 81   0.9 0.005994475

coef(model.net$finalModel, model.net$bestTune$lambda)
#changent d'une fois à l'autre
#12JN2023: le phosphore (0.0002481043) et l'ordonnée à l'origine (0.0670081665) sont les seuls paramètres significatifs

x.test.net <- model.matrix(masse.tige~ arb + coils + ves + nonmyc + dse + pH + ammonium + nitrate + phosphore, test)[,-1]
predictions.net <- model.net %>% predict(x.test.net)

data.frame(RMSE.net = RMSE(predictions.net, test$masse.tige),Rsquare.net = R2(predictions.net, test$masse.tige))
# RMSE = 0.01330959 and R² = 8.385322e-07
# RMSE = 0.0211265 and R² = 0.08918796
#12JN2023: RMSE = 0.01876104 and R² = 0.1371562


# Deuxième régression elastic net : en fonction de la longueur des racines ~~~~~~~~~~~~~~~~~~~~~~~

# Create two objects to store predictor (x) and response variables (y, median value)
Predictor.x2 <- model.matrix(l.racines~ arb + coils + ves + nonmyc + dse + pH + ammonium + nitrate + phosphore, train)[,-1]
Response.y2 <- train$l.racines

#tune parameters to identify the best alpha and lambda values
#We will tune the model by iterating over a number of alpha and lambda pairs and we can see which pair has the lowest associated error
model.net2 <- train(l.racines~ arb + coils + ves + nonmyc + dse + pH + ammonium + nitrate + phosphore, data = train, method = "glmnet",trControl = trainControl("cv", number = 10),tuneLength = 10)
model.net2$bestTune
# alpha     lambda
#8   0.1  0.1333897
#C'est différent chaque fois

coef(model.net2$finalModel, model.net2$bestTune$lambda)
#12JN2023 : tout sauf ves?
#10 x 1 sparse Matrix of class "dgCMatrix"
#s1
#(Intercept)  4.055643143
#arb         -0.014486167
#coils       -0.007330374
#ves          .          
#nonmyc      -0.026253299
#dse         -0.011520079
#pH           0.059424659
#ammonium    -0.080433777
#nitrate      0.018581527
#phosphore       0.059714695

x.test.net2 <- model.matrix(l.racines~ arb + coils + ves + nonmyc + dse + pH + ammonium + nitrate + phosphore, test)[,-1]
predictions.net2 <- model.net2 %>% predict(x.test.net2)

data.frame(RMSE.net = RMSE(predictions.net2, test$l.racines),Rsquare.net = R2(predictions.net2, test$l.racines))
# RMSE = 0.8274268 and R² = 0.03673776


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ordinations : PCA #1 ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(vegan)
library(FactoMineR)
library(ggplot2)

#fichiers de données
setwd("~/Documents/Maîtrise/Données Axe 1")


lonicera <- read.table("Lonicera.txt", header=TRUE)
lonicera$site <- as.factor(lonicera$site)
lonicera$inoculum <- as.factor(lonicera$inoculum)
#enlever le dernière ligne (stérile)
lonicera <- lonicera[-50,]

#Couleurs par sites
Sites <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5)
Sites.ellipse <- c(1, 2, 3, 4, 5)


# Colonisation des racines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lonicera.myc <- lonicera[, -c(1:5, 11:16)]

myc.pca <- rda(lonicera.myc, scale = TRUE)
summary(myc.pca)
biplot(myc.pca, scaling = 1)
biplot(myc.pca, scaling = 2)

#plus facile à visionner
PCA(lonicera.myc, scale=TRUE)

#en couleur
myc.site <- myc.pca$CA$u[, -c(3, 4, 5)]
myc.site <- cbind(myc.site, Sites)

plot(myc.pca, type="none", xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4))
points(myc.site, col=myc.site[,3])
ordiellipse(myc.site, groups = myc.site[,3], col=Sites.ellipse, label=TRUE)

#Graphique avec ggplot2
myc.site <- as.data.frame(myc.site)
myc.site$Sites <- as.factor(myc.site$Sites)

ggplot(data=myc.site, aes(x=PC1, y=PC2, colour = Sites))+
  geom_point()+
  stat_ellipse()

# Paramètres du sol ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lonicera.sol <- lonicera[, -c(1:10)]
sol.pca <- rda(lonicera.sol, scale=TRUE)
summary(sol.pca)
biplot(sol.pca, scaling = 1)
biplot(sol.pca, scaling = 2)

#plus facile à visualiser
PCA(lonicera.sol, scale=TRUE)

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
  stat_ellipse()

#Performance des plants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lonicera.pl <- lonicera[, -c(1:2, 6:16)]
plante.pca <- rda(lonicera.pl, scale = TRUE)
summary(plante.pca)
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

ggplot(data=plante.site, aes(x=PC1, y=PC2, colour = Sites))+
  geom_point()+
  stat_ellipse()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLS reg ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Ref: https://www.statology.org/partial-least-squares-in-r/
#Ref2 : https://support.minitab.com/en-us/minitab/21/help-and-how-to/statistical-modeling/regression/how-to/partial-least-squares/interpret-the-results/key-results/
#Voir exemple dans "Exemple_test.R"

library("pls")

#fichiers de données
setwd("~/Documents/Maîtrise/Données Axe 1")


lonicera <- read.table("Lonicera.txt", header=TRUE)
lonicera$site <- as.factor(lonicera$site)
lonicera$inoculum <- as.factor(lonicera$inoculum)
#enlever le dernière ligne (stérile)
lonicera <- lonicera[-50,]


# Variable réponse : masse des tiges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit PCR model (Fit Partial Least Squares Model) ~~~~~~~~~~~~~~~~~~~~~~~
#scale : centrer-réduire
#validation=”CV”: This tells R to use k-fold cross-validation to evaluate the performance of the model. 
#Note that this uses k=10 folds by default. Also note that you can specify “LOOCV” instead to perform leave-one-out cross-validation.
model.pls <- plsr(masse.tige~arb+coils+ves+nonmyc+dse+pH+ammonium+nitrate+phosphore+wsa+mo, data=lonicera, scale=TRUE, validation="CV")


#Choose the Number of PLS Components ~~~~~~~~~~~~~~~~~~~~~~~
#by looking at the test root mean squared error (test RMSE) calculated by the k-fold cross-validation
#view summary of model fitting
summary(model.pls)
#Je dirais 3 composants

#1) table : VALIDATION: RMSEP
#This table tells us the test RMSE calculated by the k-fold cross validation. We can see the following:

#If we only use the intercept term in the model, the test RMSE is 0.01647
#If we add in the first PLS component, the test RMSE goes up to 0.01815
#If we add in the second PLS component, the test RMSE goes up to 0.01878
#We can see that adding additional PLS components also increases the test RMSE
# Je ne sais pas combien de composants sont optimaux

#2) Table : TRAINING: % variance explained
#percentage of the variance in the response variable explained by the PLS components

#Using just the first PLS component :  explain 14.66% of the variation in the response variable
#adding in the second PLS component, we can explain 17.91% of the variation in the response variable.
##adding in the second PLS component, we can explain 19.19% of the variation in the response variable.
#Note that we’ll always be able to explain more variance by using more PLS components, but we can see that adding in more than two PLS components doesn’t actually increase the percentage of explained variance by much.


#visualize cross-validation plots (the test RMSE)
validationplot(model.pls) #2 ou3
validationplot(model.pls, val.type="MSEP")
validationplot(model.pls, val.type="R2")


#Use the Final Model to Make Predictions ~~~~~~~~~~~~~~~~~~~~~~~

#use the final model with two PLS components to make predictions on new observations.
#split the original dataset into a training and testing set and use the final model with two PLS components to make predictions on the testing set

#define training and testing sets
train.pls <- lonicera[1:25, c("masse.tige", "arb", "coils", "ves", "nonmyc", "dse", "pH", "ammonium", "nitrate", "phosphore", "wsa", "mo")]
y_test.pls <- lonicera[26:nrow(lonicera), c("masse.tige")]
test.pls <- lonicera[26:nrow(lonicera), c("arb", "coils", "ves", "nonmyc", "dse", "pH", "ammonium", "nitrate", "phosphore", "wsa", "mo")]

#use model to make predictions on a test set
model.pls.2 <- plsr(masse.tige~arb+coils+ves+nonmyc+dse+pH+ammonium+nitrate+phosphore+wsa+mo, data=train.pls, scale=TRUE, validation="CV")
pcr_pred.pls <- predict(model.pls.2, test.pls, ncomp=2)

#calculate RMSE
sqrt(mean((pcr_pred.pls - y_test.pls)^2))

#[1] 0.02515554


# Variable réponse : longeur de racines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit PCR model (Fit Partial Least Squares Model) ~~~~~~~~~~~~~~~~~~~~~~~
#scale : centrer-réduire
#validation=”CV”: This tells R to use k-fold cross-validation to evaluate the performance of the model. 
#Note that this uses k=10 folds by default. Also note that you can specify “LOOCV” instead to perform leave-one-out cross-validation.
model.pls.r <- plsr(l.racines~arb+coils+ves+nonmyc+dse+pH+ammonium+nitrate+phosphore+wsa+mo, data=lonicera, scale=TRUE, validation="CV")


#Choose the Number of PLS Components ~~~~~~~~~~~~~~~~~~~~~~~
#by looking at the test root mean squared error (test RMSE) calculated by the k-fold cross-validation
#view summary of model fitting
summary(model.pls.r)
#1) table : VALIDATION: RMSEP
#This table tells us the test RMSE calculated by the k-fold cross validation. We can see the following:

#If we only use the intercept term in the model, the test RMSE is 0.8225
#If we add in the first PLS component, the test RMSE drops to 0.8063
#If we add in the second PLS component, the test RMSE goes back up to 0.8121
#We can see that adding additional PLS components also increases the test RMSE (0.8324)
# 1 composant est optimal ?

#2) Table : TRAINING: % variance explained
#percentage of the variance in the response variable explained by the PLS components

#Using just the first PLS component :  explain 26.33% of the variation in the response variable
#adding in the second PLS component, we can explain 34.64% of the variation in the response variable.
##adding in the second PLS component, we can explain 38.71% of the variation in the response variable.
#Puis 39.54, 39.65, etc.
#Note that we’ll always be able to explain more variance by using more PLS components, but we can see that adding in more than two PLS components doesn’t actually increase the percentage of explained variance by much.


#visualize cross-validation plots (the test RMSE)
validationplot(model.pls.r) #il y a un pic
validationplot(model.pls.r, val.type="MSEP")
validationplot(model.pls.r, val.type="R2")


#Use the Final Model to Make Predictions ~~~~~~~~~~~~~~~~~~~~~~~

#use the final model with two PLS components to make predictions on new observations.
#split the original dataset into a training and testing set and use the final model with two PLS components to make predictions on the testing set

#define training and testing sets
train.pls.r <- lonicera[1:25, c("l.racines", "arb", "coils", "ves", "nonmyc", "dse", "pH", "ammonium", "nitrate", "phosphore", "wsa", "mo")]
y_test.pls.r <- lonicera[26:nrow(lonicera), c("l.racines")]
test.pls.r <- lonicera[26:nrow(lonicera), c("arb", "coils", "ves", "nonmyc", "dse", "pH", "ammonium", "nitrate", "phosphore", "wsa", "mo")]

#use model to make predictions on a test set
model.pls.r.2 <- plsr(l.racines~arb+coils+ves+nonmyc+dse+pH+ammonium+nitrate+phosphore+wsa+mo, data=train.pls.r, scale=TRUE, validation="CV")
pcr_pred.pls.r <- predict(model.pls.r.2, test.pls.r, ncomp=2)

#calculate RMSE
sqrt(mean((pcr_pred.pls.r - y_test.pls.r)^2))

#[1] 1.213277



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANOVAs one-way ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ATTENTION ## les noms des fichiers de données ne sont pas à jour !!

setwd("~/Documents/Maîtrise/Données Axe 1")

lonicera <- read.table("Lonicera.txt", header=TRUE)
lonicera$site <- as.factor(lonicera$site)
lonicera$inoculum <- as.factor(lonicera$inoculum)
#enlever le dernière ligne (stérile)
lonicera <- lonicera[-50,]

library(datasets)
library(ggplot2)
library(multcompView)
library(dplyr)

d.racines.aov <- aov(lonicera$diam.racine ~ lonicera$site)
#p = 0.0357*
racines.aov <- aov(lonicera$l.racines ~ lonicera$site)
#p = 0.0023 **
tiges.aov <- aov(lonicera$masse.tige ~ lonicera$site)
#p = 0.000111 ***

#Analyse des lettres ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Tuckey
d.racines.T <- TukeyHSD(d.racines.aov)
l.racines.T <- TukeyHSD(l.racines.aov)
tiges.T <- TukeyHSD(m.tiges.aov)

#attitrer les lettres
d.racines.cld <- multcompLetters4(d.racines.aov, d.racines.T)
print(d.racines.cld)

l.racines.cld <- multcompLetters4(l.racines.aov, l.racines.T)
print(l.racines.cld)

tiges.cld <- multcompLetters4(m.tiges.aov, tiges.T)
print(tiges.cld)


#longueur des racines
#tableau résumé avec lettres et 3e quartile
Tk.L <- group_by(plantes.graph, site) %>%
  summarise(mean=mean(l.racines), quant = quantile(l.racines, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
l.racines.cld.2 <- as.data.frame.list(l.racines.cld$`plantes$site`)
Tk.L$cld <- l.racines.cld.2$Letters
print(Tk.L)

#grahphique longueur
ggplot(plantes.graph, aes(x=site, y=l.racines)) + 
  geom_boxplot(aes(fill = factor(..middle..)), show.legend = FALSE) +
  scale_fill_brewer(palette = "Purples") +
  xlab("Site") +
  ylab("Longueur totale des racines (m)") +
  geom_label(data = Tk.L, aes(x = site, y = quant, label = cld), 
             size = 3, nudge_x=0.15, nudge_y =0.07)



#tiges
#tableau résumé avec lettres et 3e quartile
Tk.T <- group_by(plantes.graph, site) %>%
  summarise(mean=mean(masse.tige), quant = quantile(masse.tige, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
m.tiges.cld.2 <- as.data.frame.list(tiges.cld$`plantes$site`)
Tk.T$cld <- m.tiges.cld.2$Letters
print(Tk.T)

#grahphique longueur
ggplot(plantes.graph, aes(x=site, y=masse.tige)) + 
  geom_boxplot(aes(fill = factor(..middle..)), show.legend = FALSE) +
  scale_fill_brewer(palette = "Blues") +
  xlab("Site") +
  ylab("Masse des tiges (g)") +
  geom_label(data = Tk.T, aes(x = site, y = quant, label = cld), 
             size = 3, nudge_x=0.18, nudge_y=0.008)


#diamètre des racines
#tableau résumé avec lettres et 3e quartile
Tk.D <- group_by(plantes.graph, site) %>%
  summarise(mean=mean(diam.racine.moy), quant = quantile(diam.racine.moy, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
d.racines.cld2 <- as.data.frame.list(d.racines.cld$`plantes$site`)
Tk.D$cld <- d.racines.cld2$Letters
print(Tk.D)

#grahphique longueur
ggplot(plantes.graph, aes(x=site, y=diam.racine.moy)) + 
  geom_boxplot(aes(fill = factor(..middle..)), show.legend = FALSE) +
  scale_fill_brewer(palette = "Greens") +
  xlab("Site") +
  ylab("Diamètre moyen des racines (mm)") +
  geom_label(data = Tk.D, aes(x = site, y = quant, label = cld), 
             size = 3, nudge_x=0.18, nudge_y=0.008)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Varpart ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(vegan)


#fichiers de données
setwd("~/Documents/Maîtrise/Données Axe 1")


lonicera <- read.table("Lonicera.txt", header=TRUE)
lonicera$site <- as.factor(lonicera$site)
lonicera$inoculum <- as.factor(lonicera$inoculum)
#enlever le dernière ligne (stérile)
lonicera <- lonicera[-50,]

# longueur des racines en fonction des nutriments
(lonicera.part.all <- varpart(lonicera.perfo, lonicera.colo, lonicera.pc)) 
lonicera.part.all
plot(lonicera.part.all, 
     digits = 2, 
     Xnames = c("Colonisation", "Physicochimie"),
     bg = c("yellow", "pink"))
showvarparts(2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Genralized linear model ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#réf : https://www.datacamp.com/tutorial/generalized-linear-models

library(carData)
library(car)
library(energy)
library(multcomp)

#fichiers de données
setwd("~/Documents/Maîtrise/Données Axe 1")

lonicera <- read.table("Lonicera.txt", header=TRUE)
lonicera$site <- as.factor(lonicera$site)
lonicera$inoculum <- as.factor(lonicera$inoculum)
#enlever le dernière ligne (stérile)
lonicera <- lonicera[-50,]


#Masse des tiges ~ sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#régression normale
tige.glm <- glm(masse.tige ~ site, data = lonicera, family = "gaussian")

#vérifier les conditions
plot(tige.glm)

#normalité
shapiro.test(resid(tige.glm))
#non normale (p-value = 0.01687)

hist(lonicera$masse.tige)
qqPlot(lonicera$masse.tige)

#suit une distribution poisson ?
#régression poisson
tige.glm.pois <- glm(masse.tige ~ site, data = lonicera, family = poisson(link = "log"))
#seulement pour des décomptes (?) - ne fonctionne pas

#binomiale
tige.glm.bi <- glm(masse.tige ~ site, data = lonicera, family = binomial(link = "logit"))

#revenons à la régression normale pour l'instant
tige.glm <- glm(masse.tige ~ site, data = lonicera, family = "gaussian")

#test post-hoc
tige.tukey <- glht(tige.glm, mcp(site="Tukey"))
plot(tige.tukey)


#Longueur de racines ~ sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#régression normale
racine.glm <- glm(l.racines ~ site, data = lonicera, family = "gaussian")

racine.tukey <- glht(racine.glm, mcp(site="Tukey"))
plot(racine.tukey)


# performance ~ propriétés de sol ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

racine.sol.glm <- glm(l.racines ~ pH+ammonium+phosphore+nitrate+wsa+mo, data = lonicera, family = "gaussian")
plot(racine.sol.glm)
summary(racine.sol.glm)
#phospohre est significatif



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generalized linear mixed model ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Ref: https://search.r-project.org/CRAN/refmans/lme4/html/glmer.html
#ref2: https://rdrr.io/cran/lme4/man/glmer.html

library(lme4)
#site = facteur aléatoire
#per ~ propriété sol 1, 2, 3 et contrôler effet du site
model.test <- glmer(l.racines  ~ site + (1|phosphore), data = lonicera, family = gaussian)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analyse factorielle multiple (AFM) ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(vegan)
library(FactoMineR)
#load fonction screestick

#fichiers de données
setwd("~/Documents/Maîtrise/Données Axe 1")

lonicera <- read.table("Lonicera.txt", header=TRUE)
lonicera$site <- as.factor(lonicera$site)
lonicera$inoculum <- as.factor(lonicera$inoculum)
#enlever le dernière ligne (stérile)
lonicera <- lonicera[-50,]

#séparer les données en matrice réponse et matrices explicatives
lonicera.perfo <- lonicera[, 3:5] #s
lonicera.colo <- lonicera[, 6:10] #c
lonicera.pc <- lonicera[, 11:16] #s
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/fungal metacom.RData")
com.hel <- decostand(metacom, "hellinger") #c
com.hel <- com.hel[-c(41),]

#MFA
# MFA on 4 groups of variables: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Regroup the 3 tables (Hellinger-transformed species physiographic variables, chemical variables)
tab4 <- data.frame(lonicera.perfo, lonicera.colo, lonicera.pc, com.hel)
dim(tab4)
# Number of variables in each group
(grn <- c(ncol(lonicera.perfo), ncol(lonicera.colo), ncol(lonicera.pc), ncol(com.hel)))

# Compute the MFA without multiple plots
t4.mfa <- MFA(tab4,
              group = grn,
              type = c("s", "c", "s", "c"),
              ncp = 2,
              name.group = c("Performance", "Colonisation", "Physicochimie", "Communauté fongique"), graph = FALSE)

t4.mfa

# Plot the results
plot(t4.mfa,
     choix = "axes",
     habillage = "group",
     shadowtext = TRUE) 

plot(
       t4.mfa,
       choix = "ind",
       partial = "all",
       habillage = "group")

plot(t4.mfa,
     choix = "var",
     habillage = "group",
     shadowtext = TRUE) 

plot(t4.mfa, choix = "group")

# Eigenvalues, scree plot and broken stick model
ev <- t4.mfa$eig[, 1]
names(ev) <- paste("MFA", 1 : length(ev)) 
screestick(ev, las = 2)

# RV coefficients with tests (p-values above the diagonal of the matrix)
#lonicera.perfo, lonicera.colo, lonicera.pc, com.hel
rvp <- t4.mfa$group$RV
rvp[1, 2] <- coeffRV(scale(lonicera.perfo), lonicera.colo)$p.value
rvp[1, 3] <- coeffRV(scale(lonicera.perfo), scale(lonicera.pc))$p.value
rvp[1, 4] <- coeffRV(scale(lonicera.perfo), com.hel)$p.value 
rvp[2, 3] <- coeffRV(lonicera.colo, scale(lonicera.pc))$p.value 
rvp[2, 4] <- coeffRV(lonicera.colo, com.hel)$p.value 
rvp[3, 4] <- coeffRV(scale(lonicera.pc), com.hel)$p.value 
round(rvp[-6, -6], 6)


# MFA sans les metacommunautés: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Regroup the 3 tables (Hellinger-transformed species physiographic variables, chemical variables)
tab3 <- data.frame(lonicera.perfo, lonicera.colo, lonicera.pc)
dim(tab3)
# Number of variables in each group
(grn3 <- c(ncol(lonicera.perfo), ncol(lonicera.colo), ncol(lonicera.pc)))

# Compute the MFA without multiple plots
t3.mfa <- MFA(tab3,
              group = grn3,
              type = c("s", "c", "s"),
              ncp = 2,
              name.group = c("Performance", "Colonisation", "Physicochimie"), graph = FALSE)

t3.mfa

# Plot the results
plot(t3.mfa,
     choix = "axes",
     habillage = "group",
     shadowtext = TRUE) 

plot(
  t3.mfa,
  choix = "ind",
  partial = "all",
  habillage = "group")

plot(t3.mfa,
     choix = "var",
     habillage = "group",
     shadowtext = TRUE) 

plot(t3.mfa, choix = "group")
# Eigenvalues, scree plot and broken stick model
ev3 <- t3.mfa$eig[, 1]
names(ev3) <- paste("MFA", 1 : length(ev3)) 
screestick(ev3, las = 2)

