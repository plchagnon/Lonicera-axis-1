#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Analyses d'ADN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

## NE FONCTIONNE PAS

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Isoler les espèces importantes par Elastic Net =====
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Données d'exemple ; à remplacer

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

## Première régression elastic net : en fonction de la masse des tiges ~~~~~~~~~~~~~~~~~~~~~~~

#Préparer le fichier contenant un vecteur de performance et les communautés fongiques ~~~~~~~~~~~~~~~~~~~~~~~
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/fungal metacom.RData")
metacom <- metacom[-41,]
colnames(metacom) <- paste("champi",1:ncol(metacom),sep="")
#avec la masse des tiges
metacom <- as.data.frame(metacom)
metacom <- cbind(lonicera$masse.tige, metacom)
names(metacom)[names(metacom) == 'lonicera$masse.tige'] <- 'masse.tige'

#régression élastic net

#split the data into training and test data
sample_size <- floor(0.75 * nrow(metacom))
metacom.num <- data.matrix(metacom)
training_index <- kenStone(metacom.num, k=sample_size, metric="euclid")
train <- metacom[training_index$model, ]
test <- metacom[training_index$test, ]


# Create two objects to store predictor (x) and response variables (y, median value)
Predictor.x <- model.matrix(masse.tige~., train)[,-1]
Response.y <- train[,1]


#tune parameters to identify the best alpha and lambda values
#We will tune the model by iterating over a number of alpha and lambda pairs and we can see which pair has the lowest associated error
model.net <- train(masse.tige~., train, method = "glmnet",trControl = trainControl("cv", number = 10),tuneLength = 10)
#attention, plusieurs messages d'erreur avec le fichier de pratique
model.net$bestTune
coef(model.net$finalModel, model.net$bestTune$lambda)
#changent d'une fois à l'autre
#champignons 12, 22 et 31

x.test.net <- model.matrix(masse.tige~., test)[,-1]
predictions.net <- model.net %>% predict(x.test.net)

data.frame(RMSE.net = RMSE(predictions.net, test$masse.tige),Rsquare.net = R2(predictions.net, test$masse.tige))


## Deuxième régression elastic net : en fonction de la longueur des racines ~~~~~~~~~~~~~~~~~~~~~~~
#Probably best to remove all items from the environment first

#Préparer le fichier contenant un vecteur de performance et les communautés fongiques ~~~~~~~~~~~~~~~~~~~~~~~
load("/Users/coralie/Documents/Maîtrise/Données Axe 1/fungal metacom.RData")
metacom <- metacom[-41,]
colnames(metacom) <- paste("champi",1:ncol(metacom),sep="")
#avec la masse des tiges
metacom <- as.data.frame(metacom)
metacom <- cbind(lonicera$l.racines, metacom)
names(metacom)[names(metacom) == 'lonicera$l.racines'] <- 'l.racines'

#régression élastic net

#split the data into training and test data
sample_size <- floor(0.75 * nrow(metacom))
metacom.num <- data.matrix(metacom)
training_index <- kenStone(metacom.num, k=sample_size, metric="euclid")
train <- metacom[training_index$model, ]
test <- metacom[training_index$test, ]


# Create two objects to store predictor (x) and response variables (y, median value)
Predictor.x <- model.matrix(l.racines~., train)[,-1]
Response.y <- train[,1]


#tune parameters to identify the best alpha and lambda values
#We will tune the model by iterating over a number of alpha and lambda pairs and we can see which pair has the lowest associated error
model.net <- train(l.racines~., train, method = "glmnet",trControl = trainControl("cv", number = 10),tuneLength = 10)
#attention, plusieurs messages d'erreur avec le fichier de pratique
model.net$bestTune
coef(model.net$finalModel, model.net$bestTune$lambda)
#presque tous les champignons
#sauf 2, 6, 7, 19 et 24

x.test.net <- model.matrix(l.racines~., test)[,-1]
predictions.net <- model.net %>% predict(x.test.net)

data.frame(RMSE.net = RMSE(predictions.net, test$l.racines),Rsquare.net = R2(predictions.net, test$l.racines))


