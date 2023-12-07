#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Analyses d'ADN =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Roulé sur Calcul Canada

#Préparation des données

#appeler les packages
library("Rcpp")
library("dada2")
packageVersion("dada2")

#appeler les fichiers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("~/Documents/Maîtrise/Analyse_données génomiques/Données exemple")

## path
path <- "~/Documents/Maîtrise/Analyse_données génomiques/Données exemple" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Préparation des données génomiques ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##string manipulation to get matched lists of the forward and reverse fastq files

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Remove primers
prim <- c(22, 20)

## visualizing the quality profiles :
#of the forward reads
plotQualityProfile(fnFs[1:2])

#gray-scale is a heat map of the frequency of each quality score at each base position. 
#The mean quality score at each position is shown by the green line, 
#and the quartiles of the quality score distribution by the orange lines. 
#The red line shows the scaled proportion of reads that extend to at least that position 
#(this is more useful for other sequencing technologies, as Illumina reads are typically all the same length, hence the flat red line).
#The forward reads are good quality. We generally advise trimming the last few nucleotides to avoid less well-controlled errors that can arise there. 
#trimming the last 10 nucleotides
truncF <- 280

#of the reverse reads
plotQualityProfile(fnRs[1:2])

#The reverse reads are of significantly worse quality, especially at the end, which is common in Illumina sequencing
#truncate the reverse reads at position  where the quality distribution crashes.
truncR <- 140

# est-ce qu'il y a un overlap ?


## Filter and trim ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

#standard filtering parameters
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncF,truncR), trimleft=prim,
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

## Learn the Error Rates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#visualize the estimated error rates:
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


## Sample Inference ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# apply the core sample inference algorithm
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#inspect the return objects :
dadaFs[[1]]
dadaRs[[1]]

#inférer 104 et 108 séquences à partir des 7135 et 2753 du début


## Merge paired reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
head(mergers[[2]])
head(mergers[[3]])


## Construct sequence table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Sequences that are much longer or shorter than expected may be the result of non-specific priming. You can remove non-target-length sequences from your sequence table


## Remove chimeras ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 94 bimeras out of 161 input sequences.
dim(seqtab.nochim)
# [1]  3 67
# 3 et 63 avec MaxEE à 1

sum(seqtab.nochim)/sum(seqtab)
#0.8159767
#when we account for the abundances of those variants 
#we see they account for about 19% of the merged sequence reads.


## Track reads through the pipeline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: 
#e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

#Outside of filtering, there should no step in which a majority of reads are lost.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Assignation taxonomique ====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Avec le document dev.fasta
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/Maîtrise/Analyse_données génomiques/Données exemple/sh_general_release_dynamic_s_all_29.11.2022.fasta", multithread=TRUE)
#R encountered a fatal error

#Avec le document .fasta
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/Maîtrise/Analyse_données génomiques/Données exemple/sh_general_release_dynamic_s_all_29.11.2022.fasta", multithread=TRUE)

library(bioseq)
read.fasta("~/Documents/Maîtrise/Analyse_données génomiques/Données exemple/sh_general_release_dynamic_s_all_29.11.2022.fasta")

## NE FONCTIONNE PAS
# Pas assez de mémoire vive ; utiliser calcul canada

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Cluster sequences by similarity =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library("bioseq")
help(seq_cluster)

#test 1
seq.test <- seq_cluster(seqtab.nochim, threshold = 0.05, method = "complete")

## NE FONCTIONNE PAS : à faire sur calcul canada

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Traiter le jeu de données =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
load("/Users/coralie/Documents/Maîtrise/Analyse_données génomiques/Données finales/Coralie ITS.RData")

#Modifier les noms de colonnes
rownames(taxa)=paste0("ASV",1:ncol(comm))
colnames(comm) = rownames(taxa)

#Modifier les numéros de rangées
rangee <- c(1:39, 41:50)
rownames(lonicera) = rangee
#rangee.moins <- c(1, 2, 3, 5, 8, 9, 11, 12, 13, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 29, 32, 33, 34, 37, 41, 42, 43, 45, 47, 49, 50)


#données de performance sans les échantillons manquants des métacommunautés
perfo.comm <- lonicera [,1:5]
perfo.comm <- perfo.comm[rownames(comm),]



#toutes les variables sans les échantillons manquants des métacommunautés
lonicera.comm <- lonicera[rownames(comm),]

save.image("Lonicera.Coralie.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Isoler les espèces importantes par Elastic Net =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#réf: https://hackernoon.com/an-introduction-to-ridge-lasso-and-elastic-net-regression-cca60b4b934f

#Importer les données ~~~~~~~~~~~~~~~~~~~~~~~

setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")
#comm : métacommunauté (rangées = sites, colonnes = champignons) pour 31 échantillons de sol
#taxa : taxonomie (autant de rangées que de champignons, et rangs taxonomiques en colonnes)
#lonicera : variables de performance, de colonisation racinaire et de physicochimie, pour 49 échantillons de sol
#croissance : variables de performance, pour 314 pots de Lonicera
#perfo.comm : variables de performance, pour les 31 échantillons ayant des données fongiques
#lonicera.comm : variables de performance, de colonisation racinaire et de physicochimie, pour les 31 échantillons ayant des données fongiques

library(tidyverse)
library(caret)
library(glmnet)
library(prospectr)
library(vegan)
library(cgwtools)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Elastic net #1 : masse des tiges =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Préparer le fichier contenant un vecteur de performance et les communautés fongiques ~~~~~~~~~~~~~~~~~~~~~~~

#avec la masse des tiges
comm.tige <- as.data.frame(comm)
comm.tige <- cbind(perfo.comm$masse.tige, comm.tige)
names(comm.tige)[names(comm.tige) == 'perfo.comm$masse.tige'] <- 'masse.tige'

#régression élastic net

#split the data into training and test data
sample_size <- floor(0.9 * nrow(comm.tige))
comm.tige.num <- data.matrix(comm.tige)
training_index <- kenStone(comm.tige.num, k=sample_size, metric="euclid")
train <- comm.tige[training_index$model, ]
test <- comm.tige[training_index$test, ]

# Create two objects to store predictor (x) and response variables (y, median value)
Predictor.x <- model.matrix(masse.tige~., train)[,-1]
Response.y <- train[,1]

#tune parameters to identify the best alpha and lambda values
#We will tune the model by iterating over a number of alpha and lambda pairs and we can see which pair has the lowest associated error
model.net <- train(masse.tige~., train, method = "glmnet",trControl = trainControl("cv", number = 10),tuneLength = 10)
#attention, plusieurs messages d'erreur avec le fichier de pratique
model.net$bestTune
#avec 75% des données dans train
#alpha      lambda
#71   0.8 0.007908121
#avec 90% des données dans train
#   alpha    lambda
#45   0.5 0.0111314

#varImp de caret
varImp(model.net)
varImp(model.net, scale=FALSE)

#trouver les espèces à contribution significative
all.coef <- coef(model.net$finalModel, model.net$bestTune$lambda)
champignons <- all.coef[!all.coef[,1]==0,]
#7 espèces
# ASV 107, 185, 188, 192, 270, 593 et 820

x.test.net <- model.matrix(masse.tige~., test)[,-1]
predictions.net <- model.net %>% predict(x.test.net)

data.frame(RMSE.net = RMSE(predictions.net, test$masse.tige),Rsquare.net = R2(predictions.net, test$masse.tige))
#avec 90% des données dans train
#RMSE.net Rsquare.net
#1 0.008256518   0.1619261


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Elastic net #2 : longueur des racines =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Probably best to remove all items from the environment first

setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")
#comm : métacommunauté (rangées = sites, colonnes = champignons) pour 31 échantillons de sol
#taxa : taxonomie (autant de rangées que de champignons, et rangs taxonomiques en colonnes)
#lonicera : variables de performance, de colonisation racinaire et de physicochimie, pour 49 échantillons de sol
#croissance : variables de performance, pour 314 pots de Lonicera
#perfo.comm : variables de performance, pour les 31 échantillons ayant des données fongiques
#lonicera.comm : variables de performance, de colonisation racinaire et de physicochimie, pour les 31 échantillons ayant des données fongiques

library(tidyverse)
library(caret)
library(glmnet)
library(prospectr)

#Préparer le fichier contenant un vecteur de performance et les communautés fongiques ~~~~~~~~~~~~~~~~~~~~~~~

#avec la longueur des racines
comm.racine <- as.data.frame(comm)
comm.racine <- cbind(perfo.comm$l.racines, comm.racine)
names(comm.racine)[names(comm.racine) == 'perfo.comm$l.racines'] <- 'l.racines'


#régression élastic net

#split the data into training and test data
sample_size <- floor(0.9 * nrow(comm.racine))
comm.racine.num <- data.matrix(comm.racine)
training_index <- kenStone(comm.racine.num, k=sample_size, metric="euclid")
train <- comm.racine[training_index$model, ]
test <- comm.racine[training_index$test, ]


# Create two objects to store predictor (x) and response variables (y, median value)
Predictor.x <- model.matrix(l.racines~., train)[,-1]
Response.y <- train[,1]


#tune parameters to identify the best alpha and lambda values
#We will tune the model by iterating over a number of alpha and lambda pairs and we can see which pair has the lowest associated error
model.net <- train(l.racines~., train, method = "glmnet",trControl = trainControl("cv", number = 10),tuneLength = 10)
#attention, plusieurs messages d'erreur avec le fichier de pratique
model.net$bestTune
#avec 75% des données dans train
#   alpha    lambda
#81   0.9 0.6092961

#avec 90% des données dans train
#alpha    lambda
#81   0.9 0.5562389

#trouver les espèces à contribution significative
all.coef <- coef(model.net$finalModel, model.net$bestTune$lambda)
champignons.racine <- all.coef[!all.coef[,1]==0,]

#combiner les deux fichiers
champignons <- cbind(champignons, all.coef[!all.coef[,1]==0,])

#aucune espèce

x.test.net <- model.matrix(l.racines~., test)[,-1]
predictions.net <- model.net %>% predict(x.test.net)

data.frame(RMSE.net = RMSE(predictions.net, test$l.racines),Rsquare.net = R2(predictions.net, test$l.racines))
#   RMSE.net Rsquare.net
#1 0.1753813          NA


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Elastic net #3 : diamètre des racines =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Probably best to remove all items from the environment first

setwd("~/Documents/Maîtrise/Données Axe 1")
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")
#comm : métacommunauté (rangées = sites, colonnes = champignons) pour 31 échantillons de sol
#taxa : taxonomie (autant de rangées que de champignons, et rangs taxonomiques en colonnes)
#lonicera : variables de performance, de colonisation racinaire et de physicochimie, pour 49 échantillons de sol
#croissance : variables de performance, pour 314 pots de Lonicera
#perfo.comm : variables de performance, pour les 31 échantillons ayant des données fongiques
#lonicera.comm : variables de performance, de colonisation racinaire et de physicochimie, pour les 31 échantillons ayant des données fongiques
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Autre.Rdata")

library(tidyverse)
library(caret)
library(glmnet)
library(prospectr)
library(vegan)

#Préparer le fichier contenant un vecteur de performance et les communautés fongiques ~~~~~~~~~~~~~~~~~~~~~~~

#avec la longueur des racines
comm.diam <- as.data.frame(comm)
comm.diam <- cbind(perfo.comm$diam.racine, comm.diam)
names(comm.diam)[names(comm.diam) == 'perfo.comm$diam.racine'] <- 'd.racines'


#régression élastic net

#split the data into training and test data
sample_size <- floor(0.9 * nrow(comm.diam))
comm.diam.num <- data.matrix(comm.diam)
training_index <- kenStone(comm.diam.num, k=sample_size, metric="euclid")
train <- comm.diam[training_index$model, ]
test <- comm.diam[training_index$test, ]


# Create two objects to store predictor (x) and response variables (y, median value)
Predictor.x <- model.matrix(d.racines~., train)[,-1]
Response.y <- train[,1]


#tune parameters to identify the best alpha and lambda values
#We will tune the model by iterating over a number of alpha and lambda pairs and we can see which pair has the lowest associated error
model.net <- train(d.racines~., train, method = "glmnet",trControl = trainControl("cv", number = 10),tuneLength = 10)
#attention, plusieurs messages d'erreur avec le fichier de pratique
model.net$bestTune
#avec 90% des données dans train
#alpha     lambda
#90     1 0.01143887

#trouver les espèces à contribution significative
all.coef <- coef(model.net$finalModel, model.net$bestTune$lambda)
champignons.diam <- all.coef[!all.coef[,1]==0,]
#aucune espèce

#combiner les deux fichiers
champignons <- cbind(champignons, all.coef[!all.coef[,1]==0,])

x.test.net <- model.matrix(d.racines~., test)[,-1]
predictions.net <- model.net %>% predict(x.test.net)

data.frame(RMSE.net = RMSE(predictions.net, test$l.racines),Rsquare.net = R2(predictions.net, test$l.racines))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Ajout des espèces de champignons sélectionnées =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

champignons
#7 espèces
# ASV 107, 185, 188, 192, 270, 593 et 820
ASV <- cbind(comm$ASV107, comm$ASV185, comm$ASV188, comm$ASV192, comm$ASV270, comm$ASV593, comm$ASV820)
rownames(ASV) = rownames(comm)
ASV.names <- c(taxa[107,6], paste0(taxa[185,6], taxa[185,7]), paste0(taxa[188,6], taxa[188,7]), taxa[192,6], paste0(taxa[270,6], taxa[270,7]), paste0(taxa[593,6], taxa[593,7]), paste0(taxa[820,6], taxa[820,7]))
ASV.names <- c("Calonectria_sp.", "Rhizophaguss_irregularis", "Rhizophaguss_intraradices", "Claroideoglomu_sp.", "Rhizophaguss_intraradices", "Podilas_humilis", "Minutisphaeras_aspera")
colnames(ASV) = ASV.names
ASV = as.data.frame(ASV)

#préparer un nouveau fichier RData
setwd("~/Documents/Maîtrise/Données Axe 1")

## not run to avoid creating detritus
save(champignons,file='Autre.Rdata')
resave(ASV,file='Autre.Rdata')
resave(ASV.names,file='Autre.Rdata')
#check your work
lsdata('Autre.RData')

#ajouter les ASV à lonicera.comm
lonicera.comm <- cbind(lonicera.comm, ASV)

#comm : métacommunauté (rangées = sites, colonnes = champignons) pour 31 échantillons de sol
#taxa : taxonomie (autant de rangées que de champignons, et rangs taxonomiques en colonnes)
#ASV : abondance de 7 espèces de champignons jugées pertinentes par régression elastic net, pour 31 échantillons de sol
#lonicera : variables de performance, de colonisation racinaire et de physicochimie, pour 49 échantillons de sol
#croissance : variables de performance, pour 314 pots de Lonicera
#perfo.comm : variables de performance, pour les 31 échantillons ayant des données fongiques
#lonicera.comm : variables de performance, de colonisation racinaire, de physicochimie, et abondance des ASV de champignons pertinentes, pour les 31 échantillons ayant des données fongiques


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Forest plot =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(forestplot)
library(dplyr)

load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Autre.RData")


champignons <- as.data.frame(champignons)
new.champignons <- cbind(rownames(champignons), champignons$champignons)
new.champignons <- new.champignons[-1,]
colnames(new.champignons) = (c("ASV", "Coefficient"))
new.champignons <- as.data.frame(new.champignons)
new.champignons$ASV <- as.factor(new.champignons$ASV)
new.champignons$Coefficient <- as.numeric(new.champignons$Coefficient)

resave(new.champignons,file='Autre.Rdata')
resave(champignons,file='Autre.Rdata')


#Standard deviation error
sd(new.champignons$Coefficient)

ggplot(new.champignons, aes(x = ASV, y=Coefficient)) + 
  geom_point() +
  geom_hline(yintercept = 0, col = "gray", size = 0.4, linetype="dashed")

#arguments non-utilisés
#geom_errorbar(aes(ymin=Coefficient-sd, ymax=Coefficient+sd), width=.2) +
#  geom_errorbar(mapping = NULL, data = NULL, stat = “identity”, position = “identity”)
#geom_errorbarh( height )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Elastic net #4 : masse des tiges sur données transformées Hellinger =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")

#Préparer le fichier contenant un vecteur de performance et les communautés fongiques ~~~~~~~~~~~~~~~~~~~~~~~

comm.tige <- as.data.frame(comm)
comm.tige.hell <- decostand(comm.tige, "hellinger")
comm.tige.hell <- cbind(perfo.comm$masse.tige, comm.tige.hell)
names(comm.tige.hell)[names(comm.tige.hell) == 'perfo.comm$masse.tige'] <- 'masse.tige'

#régression élastic net

#split the data into training and test data
sample_size <- floor(0.9 * nrow(comm.tige.hell))
comm.tigeH.num <- data.matrix(comm.tige.hell)
training_index <- kenStone(comm.tigeH.num, k=sample_size, metric="euclid")
train <- comm.tige.hell[training_index$model, ]
test <- comm.tige.hell[training_index$test, ]

# Create two objects to store predictor (x) and response variables (y, median value)
Predictor.x <- model.matrix(masse.tige~., train)[,-1]
Response.y <- train[,1]

#tune parameters to identify the best alpha and lambda values
#We will tune the model by iterating over a number of alpha and lambda pairs and we can see which pair has the lowest associated error
model.net <- train(masse.tige~., train, method = "glmnet",trControl = trainControl("cv", number = 10),tuneLength = 10)
#attention, plusieurs messages d'erreur avec le fichier de pratique
model.net$bestTune
#avec 90% des données dans train
#   alpha    lambda
#88     1 0.004851214
#VARIE D'UNE FOIS À L'AUTRE
#36   0.4 0.01120693

#varImp de caret
varImp(model.net)
varImp(model.net, scale=FALSE)

#trouver les espèces à contribution significative
all.coef <- coef(model.net$finalModel, model.net$bestTune$lambda)
champignons4 <- all.coef[!all.coef[,1]==0,]
#9 espèces
# ASV 61, 76, 107, 117, 203, 535, 593, 614, 820
#3 en commun avec elastic net sur données brutes : 820, 593 et 107

#deuxième fois : 16 espèces (!)
#ASV 820, 203, 614, 535, 61, 1122, 593, 107, 76, 87, 117, 408, 235, 28, 135, 348

x.test.net <- model.matrix(masse.tige~., test)[,-1]
predictions.net <- model.net %>% predict(x.test.net)

data.frame(RMSE.net = RMSE(predictions.net, test$masse.tige),Rsquare.net = R2(predictions.net, test$masse.tige))
#avec 90% des données dans train
#     RMSE.net Rsquare.net
#1 0.01214996  0.05472083


#Ajouter les espèces sélectionnées à un fichier 
champignons4
ASV.4 <- cbind(comm$ASV28, comm$ASV61, comm$ASV76, comm$ASV87, comm$ASV107, comm$ASV117, comm$ASV135, comm$ASV203, comm$ASV235, comm$ASV348, comm$ASV408, comm$ASV535, comm$ASV593, comm$ASV614, comm$ASV820, comm$ASV1122)
rownames(ASV.4) = rownames(comm)
ASV.names4 <- c(paste0(taxa[28,6], taxa[28,7]),
                paste0(taxa[61,6], taxa[61,7]),
                paste0(taxa[76,6], taxa[76,7]),
                paste0(taxa[87,6], taxa[87,7]),
                paste0(taxa[107,6], taxa[107,7]),
                paste0(taxa[117,6], taxa[117,7]),
                paste0(taxa[135,6], taxa[135,7]),
                paste0(taxa[203,6], taxa[203,7]),
                paste0(taxa[235,6], taxa[235,7]), 
                paste0(taxa[348,6], taxa[348,7]), 
                paste0(taxa[408,6], taxa[408,7]), 
                paste0(taxa[535,6], taxa[535,7]), 
                paste0(taxa[593,6], taxa[593,7]), 
                paste0(taxa[614,6], taxa[614,7]), 
                paste0(taxa[820,6], taxa[820,7]),
                paste0(taxa[1122,6], taxa[1122,7]))
colnames(ASV.4) = ASV.names4
ASV.4 = as.data.frame(ASV.4)
ASV.hell <- ASV.4

#ajouter au fichier RData
setwd("~/Documents/Maîtrise/Données Axe 1")
resave(ASV.hell,file='Autre.Rdata')
resave(ASV.names4,file='Autre.Rdata')


#graphique
champignons4 <- as.data.frame(champignons4)

new.champignons4 <- cbind(rownames(champignons4), champignons4$champignons4)
new.champignons4 <- as.data.frame(new.champignons4)
new.champignons4 <- new.champignons4[-1,]
colnames(new.champignons4) = (c("ASV", "Coefficient"))

new.champignons4$ASV <- as.factor(new.champignons4$ASV)
new.champignons4$Coefficient <- as.numeric(new.champignons4$Coefficient)

resave(new.champignons4,file='Autre.Rdata')

#graphique
ggplot(new.champignons4, aes(x = ASV, y=Coefficient)) + 
  geom_point() +
  geom_hline(yintercept = 0, col = "gray", size = 0.4, linetype="dashed")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Elastic net #5 : masse des tiges sur données transformées chi-carré =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Lonicera.Coralie.RData")

#Préparer le fichier contenant un vecteur de performance et les communautés fongiques ~~~~~~~~~~~~~~~~~~~~~~~

comm.tige <- as.data.frame(comm)
comm.tige.chi <- decostand(comm.tige, "chi.square")
comm.tige.chi <- cbind(perfo.comm$masse.tige, comm.tige.chi)
names(comm.tige.chi)[names(comm.tige.chi) == 'perfo.comm$masse.tige'] <- 'masse.tige'

#régression élastic net

#split the data into training and test data
sample_size <- floor(0.9 * nrow(comm.tige.chi))
comm.tigeC.num <- data.matrix(comm.tige.chi)
training_index <- kenStone(comm.tigeC.num, k=sample_size, metric="euclid")
train <- comm.tige.chi[training_index$model, ]
test <- comm.tige.chi[training_index$test, ]

# Create two objects to store predictor (x) and response variables (y, median value)
Predictor.x <- model.matrix(masse.tige~., train)[,-1]
Response.y <- train[,1]

#tune parameters to identify the best alpha and lambda values
#We will tune the model by iterating over a number of alpha and lambda pairs and we can see which pair has the lowest associated error
model.net <- train(masse.tige~., train, method = "glmnet",trControl = trainControl("cv", number = 10),tuneLength = 10)
#attention, plusieurs messages d'erreur avec le fichier de pratique
model.net$bestTune
#avec 90% des données dans train
#   alpha    lambda
#45     0.5 0.01095148

#varImp de caret
varImp(model.net)
varImp(model.net, scale=FALSE)

#trouver les espèces à contribution significative
all.coef <- coef(model.net$finalModel, model.net$bestTune$lambda)
champignons5 <- all.coef[!all.coef[,1]==0,]
#9 espèces
# ASV 593, 203, 820, 270, 185, 348, 188, 192, 535
#6 en commun avec elastic net sur données brutes : 593, 203, 820, 270, 185, 188 et 192

x.test.net <- model.matrix(masse.tige~., test)[,-1]
predictions.net <- model.net %>% predict(x.test.net)

data.frame(RMSE.net = RMSE(predictions.net, test$masse.tige),Rsquare.net = R2(predictions.net, test$masse.tige))
#avec 90% des données dans train
#     RMSE.net Rsquare.net
#1 0.009633471   0.1219671

#Ajouter les espèces sélectionnées à un fichier 
champignons5
ASV.5 <- cbind(comm$ASV185, comm$ASV188, comm$ASV192, comm$ASV203, comm$ASV270, comm$ASV348, comm$ASV535, comm$ASV593, comm$ASV820)
rownames(ASV.5) = rownames(comm)
ASV.names5 <- c(paste0(taxa[185,6], taxa[185,7]), 
               paste0(taxa[188,6], taxa[188,7]), 
               paste0(taxa[192,6]), 
               paste0(taxa[203,6], taxa[203,7]), 
               paste0(taxa[270,6], taxa[270,7]), 
               paste0(taxa[348,6], taxa[348,7]), 
               paste0(taxa[535,6], taxa[535,7]), 
               paste0(taxa[593,6], taxa[593,7]), 
               paste0(taxa[820,6], taxa[820,7]))
#ASV.names <- c("Calonectria_sp.", "Rhizophaguss_irregularis", "Rhizophaguss_intraradices", "Claroideoglomu_sp.", "Rhizophaguss_intraradices", "Podilas_humilis", "Minutisphaeras_aspera")
colnames(ASV.5) = ASV.names5
ASV.5 = as.data.frame(ASV.5)
ASV.chi <- ASV.5

#ajouter au fichier RData
setwd("~/Documents/Maîtrise/Données Axe 1")

resave(ASV.chi,file='Autre.Rdata')
resave(ASV.names5,file='Autre.Rdata')
resave(champignons5,file='Autre.Rdata')

lsdata('Autre.RData')

#graphique
champignons5 <- as.data.frame(champignons5)

new.champignons5 <- cbind(rownames(champignons5), champignons5$champignons5)
new.champignons5 <- as.data.frame(new.champignons5)
new.champignons5 <- new.champignons5[-1,]
colnames(new.champignons5) = (c("ASV", "Coefficient"))

new.champignons5$ASV <- as.factor(new.champignons5$ASV)
new.champignons5$Coefficient <- as.numeric(new.champignons5$Coefficient)

resave(new.champignons5,file='Autre.Rdata')

#graphique
ggplot(new.champignons5, aes(x = ASV, y=Coefficient)) + 
  geom_point() +
  geom_hline(yintercept = 0, col = "gray", size = 0.4, linetype="dashed")

saveRDS(new.champignons5, file = "Autre.RData")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Elastic net #6 : masse des tiges sur abondances relatives =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load("/Users/coralie/Documents/Maîtrise/Données Axe 1/Autre.Rdata")


rel=function(x){x/sum(x)}
comm.rel=t(apply(comm,1,rel))

#corrélation
comm.rel <- as.matrix(comm.rel)
comm <- as.matrix(comm)

cor.test(c(comm), c(comm.rel))


#Préparer le fichier contenant un vecteur de performance et les communautés fongiques ~~~~~~~~~~~~~~~~~~~~~~~

comm.rel <- as.data.frame(comm.rel)
comm.rel <- cbind(perfo.comm$masse.tige, comm.rel)
names(comm.rel)[names(comm.rel) == 'perfo.comm$masse.tige'] <- 'masse.tige'

#régression élastic net

#split the data into training and test data
sample_size <- floor(0.9 * nrow(comm.rel))
comm.rel.num <- data.matrix(comm.rel)
training_index <- kenStone(comm.rel.num, k=sample_size, metric="euclid")
train <- comm.rel[training_index$model, ]
test <- comm.rel[training_index$test, ]

# Create two objects to store predictor (x) and response variables (y, median value)
Predictor.x <- model.matrix(masse.tige~., train)[,-1]
Response.y <- train[,1]

#tune parameters to identify the best alpha and lambda values
#We will tune the model by iterating over a number of alpha and lambda pairs and we can see which pair has the lowest associated error
model.net <- train(masse.tige~., train, method = "glmnet",trControl = trainControl("cv", number = 10),tuneLength = 10)
#attention, plusieurs messages d'erreur avec le fichier de pratique
model.net$bestTune
#avec 90% des données dans train
#   alpha    lambda
#62     0.7 0.007345175

#varImp de caret
varImp(model.net)
varImp(model.net, scale=FALSE)

#trouver les espèces à contribution significative
all.coef <- coef(model.net$finalModel, model.net$bestTune$lambda)
champignons6 <- all.coef[!all.coef[,1]==0,]
#6 espèces seulement
# ASV 593, 820, 860, 270, 185, 107
#5 en commun avec elastic net sur données brutes : 593, 820, 270, 185, 107
#188 et 192 ont disparut
# et 860 est nouveaux

x.test.net <- model.matrix(masse.tige~., test)[,-1]
predictions.net <- model.net %>% predict(x.test.net)

data.frame(RMSE.net = RMSE(predictions.net, test$masse.tige),Rsquare.net = R2(predictions.net, test$masse.tige))
#avec 90% des données dans train
#     RMSE.net Rsquare.net
#1 0.008566437   0.2076059

#Ajouter les espèces sélectionnées à un fichier 
champignons6
ASV.6 <- cbind(comm.rel$ASV107, comm.rel$ASV185, comm.rel$ASV270, comm.rel$ASV593, comm.rel$ASV820, comm.rel$ASV860)
rownames(ASV.6) = rownames(comm.rel)
ASV.names6 <- c(paste0(taxa[107,6], taxa[107,7]),
                paste0(taxa[185,6], taxa[185,7]), 
                paste0(taxa[270,6], taxa[270,7]), 
                paste0(taxa[593,6], taxa[593,7]), 
                paste0(taxa[820,6], taxa[820,7]),
                paste0(taxa[860,6], taxa[860,7]))
colnames(ASV.6) = ASV.names6
ASV.6 = as.data.frame(ASV.6)


#ajouter au fichier RData
setwd("~/Documents/Maîtrise/Données Axe 1")

resave(ASV.6,file='Autre.Rdata')
resave(ASV.names6,file='Autre.Rdata')
resave(champignons6,file='Autre.Rdata')

lsdata('Autre.RData')

#graphique
champignons6 <- as.data.frame(champignons6)

new.champignons6 <- cbind(rownames(champignons6), champignons6$champignons6)
new.champignons6 <- as.data.frame(new.champignons6)
new.champignons6 <- new.champignons6[-1,]
new.champignons6 <- cbind(ASV.names6, new.champignons6)
colnames(new.champignons6) = (c("Taxo", "ASV", "Coefficient"))

new.champignons6$ASV <- as.factor(new.champignons6$ASV)
new.champignons6$Coefficient <- as.numeric(new.champignons6$Coefficient)

resave(new.champignons6,file='Autre.Rdata')

#graphique
ggplot(new.champignons6, aes(x = ASV, y=Coefficient)) + 
  geom_point() +
  geom_hline(yintercept = 0, col = "gray", size = 0.4, linetype="dashed")

saveRDS(new.champignons5, file = "Autre.RData")
