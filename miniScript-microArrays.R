
# numéro à changer pour chaque essai
setwd("~/m1meg-ue1-tp5-r-4")

# Chargment des packages dans l'environnement présent
library(Matrix)
library(lattice)
library(fdrtool)
library(rpart)
library(Biobase)
library(Biostrings)
library(mouse4302cdf)
library(ggplot2)
library(limma)
library(affy)
library(affyPLM)
library(simpleaffy)
library(genefilter)

sessionInfo()

# 1 - Récupération des données
celpath = "/srv/data/meg-m1-ue1/DataTP5/"

# import des données d'intensité de fluorescence brutes des sondes depuis les fichiers CEL , et ajout dans un objet R AffyBatch
data = ReadAffy(celfile.path=celpath)

# Cet objet AffyBatch est accessible via des méthodes contenues dans le package affy 
# Récupération de l'annotation des échantillons
ph = data@phenoData
ph@data

# 2 - Contrôle qualité des données de puces

# Attribution de noms informatifs aux échantillons
ph@data[ ,1] = c("Testis_E13_1","Testis_E13_2","Testis_E13_3","Testis_P0_1","Testis_P0_2","Testis_P0_3","Ovary_E13_1","Ovary_E13_2","Ovary_E13_3", "Ovary_P0_1", "Ovary_P0_2", "Ovary_P0_3")
ph@data

# Histogramme de la distribution des données (avec ggplot2)

# Récupérations des intensités PM (sondes Perfect Match dans les probe sets)
pmexp = pm(data)

# initialisation de 3 vecteurs
sampleNames = vector()
nologs = vector()

# remplissage des vecteurs avec les données
for (i in 1:12)
{
  sampleNames = c(sampleNames,rep(ph@data[i,1],dim(pmexp)[1]))
  nologs = c(nologs, pmexp[,i])
}
nologs[1:10]

# Après avoir rempli les vecteurs, nous combinons sampleNames et nologs ou logs dans 2 dataframes
nologData = data.frame(nologInt=nologs,sampleName=sampleNames)

nologData[1:5,]

# Maintenant nous pouvons créer l'histogramme: 
dataHist1 = ggplot(nologData, aes(nologInt, colour = sampleName)) 
dataHist1 + geom_density()

# 3 - Normalisation des données
library("gcrma")
data.gcrma = gcrma(data)
# création d'une data farame avec seulement les valeurs d'expression normalisées
normexpr.gcrma = data.frame(exprs(data.gcrma))
normexpr.gcrma[1:5,]


# Box plots des valeurs normalisées
# initialisation des vecteurs
sampleNames = vector()
normlogs = vector()

# remplissage des vecteurs
for (i in 1:12)
{
  sampleNames = c(sampleNames,rep(ph@data[i,1],dim(data.gcrma)[1]))
  normlogs = c(normlogs,normexpr.gcrma[,i])
}

normlogs[1:10]

# création de la dataframe avec les 2 vecteurs
normData = data.frame(norm_logInt=normlogs,sampleName=sampleNames)
normData[1:5,]

# création des graphs
dataBox3 = ggplot(normData, aes(sampleName,norm_logInt))
dataBox3 + geom_boxplot() + ylim(2,16) + ggtitle("after normalization")




