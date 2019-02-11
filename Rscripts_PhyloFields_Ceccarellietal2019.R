############################################################################
## Functions to generate range map polygons using the alphahull R package ##
############################################################################

# Set your working directory (dir)
setwd(dir)

# Load packages
library(sp)
library(alphahull)
library(maptools)

# Load the ConvexHull function
source ("ConvexHull.R")

# Import points shape files of each species and find point pairs with equal spatial coordinates
data.shape <- readShapePoints("CompilacionTotal_63sp_Acronimo_Tsor.shp")
data.shape2 <- data.shape[-zerodist(data.shape )[,1],]
data.shape <- data.shape2

# Test multiple alpha values that best fit with your occurence data (i.e. a= 8, 10, 12, 14, 15, 18)
a=c(8, 10, 12, 14, 15, 18)
par(mfrow=c(2,3))
for (alpha in a) {
  ch <- ConvexHull(data.shape, alpha=alpha)
  plot(ch, axes= TRUE)
  box()
  points(data.shape, pch=19) 
}

#alpha value
b <- ConvexHull(data.shape, alpha=14)

#write Spatial data Shape
writeSpatialShape(b, "Tsor")

#########################################################################################################################################
## Functions to import phylogenetic and spatial data, create presence-absence matrices and calculate diversity and phylogenetic fields ##
#########################################################################################################################################

# Load packages
library(maptools)
library(letsR)
library(picante)
library(foreign)
library(rgdal)
library(ape)

#Load the phylogentic data (in parenthetic format)
arbol <- read.tree("arbol_check")

# Load 'divfield_pams', 'randomrange_phylofields' and 'compare_phylofields' functions
source ("div_phylo_fields.R")
source ("randomrange_pf2.R")
source ("compare_PFs.R")

#Load your spatial data
####Points####
#get the data (as a table, not as a shapefile)
puntos.table <- read.dbf("Puntos_63spFINAL.dbf")
#create the PAM (presenceAbsence object)
puntos.table.pam1 <- lets.presab.points(puntos.table[,1:2],puntos.table[,3],resol = 1)
puntos.table.pam05 <- lets.presab.points(puntos.table[,1:2],puntos.table[,3],resol = 0.5)
#save only the presence-absence matrix
puntos.pam.clean <- puntos.table.pam1$P[,-c(1,2)]
puntos.pam.clean05 <- puntos.table.pam05$P[,-c(1,2)]

#calculate the diversity fields and matrices
#change the WD to save the matrices
#1x1°
setwd("SpxSp_puntos")
puntos.divfields <- divfield_pams(pam_obs = puntos.pam.clean)
#get the sp_assemblages
sp_assemb(spp.divfields = puntos.divfields)
#calculate the phylogenetic fields
puntos.phylofields <- phylofields_psv(spp.divfields = puntos.divfields,tree = arbol)
#0.5°#
setwd("puntos_05")
puntos.divfields05 <- divfield_pams(pam_obs = puntos.pam.clean05)
#get the sp_assemblages
sp_assemb(spp.divfields = puntos.divfields05)
#calculate the phylogenetic fields
puntos.phylofields05 <- phylofields_psv(spp.divfields = puntos.divfields05,tree = arbol)


####Polygons####
poligonos <- readOGR("Poligonos_57spFINAL.shp")
poligonos.pam1 <- lets.presab(poligonos,resol = 1)
poligonos.pam05 <- lets.presab(poligonos,resol = 0.5)

#save only the pres-ab matrix
poligonos.pam.clean <- poligonos.pam1$P[,-c(1,2)]
poligonos.pam.clean05 <- poligonos.pam05$P[,-c(1,2)]

#calculate the diversity fields and matrices
#change the WD to save the matrices
#1x1°
setwd("SpxSp")
poligonos.divfields <- divfield_pams(pam_obs = poligonos.pam.clean)
#get the sp_assemblages
sp_assemb(spp.divfields = poligonos.divfields)
#calculate the phylogenetic fields
poligonos.phylofields <- phylofields_psv(spp.divfields = poligonos.divfields,tree = arbol)
#0.5°
setwd("poligonos_pam05")
poligonos.divfields05 <- divfield_pams(pam_obs = poligonos.pam.clean05)
#get the sp_assemblages
sp_assemb(spp.divfields = poligonos.divfields05)
#calculate the phylogenetic fields
poligonos.phylofields05 <- phylofields_psv(spp.divfields = poligonos.divfields05,tree = arbol)

#### randomizations####
#Points#
#0.5°
randomrange_phylofields(PAM = puntos.pam.clean05, spp.divfields = puntos.divfields05, path = "SpxSp_puntos/puntos_05/", tree = arbol, sims = 1000)
compare_phylofields(observed = "SpxSp_puntos/puntos_05/puntos05_phylo_fields.txt",nulls = "SpxSp_puntos/puntos_05/puntos05_random_phylofields_constrained.txt")
#1°
randomrange_phylofields(PAM = puntos.pam.clean, spp.divfields = puntos.divfields, path = "SpxSp_puntos/", tree = arbol, sims = 1000)
setwd("SpxSp_puntos")
compare_phylofields(observed = "puntos1_phylo_fields.txt",nulls = "puntos1_random_phylofields_constrained.txt")

#Polygons#
#0.5°
randomrange_phylofields(PAM = poligonos.pam.clean05, spp.divfields = poligonos.divfields05, path = "SpxSp/poligonos_pam05/", tree = arbol, sims = 1000)
setwd("SpxSp/poligonos_pam05")
compare_phylofields(observed = "poligonos05_phylo_fields.txt",nulls = "poligonos05_random_phylofields_constrained.txt")
#1°
randomrange_phylofields(PAM = poligonos.pam.clean, spp.divfields = poligonos.divfields, path = "SpxSp/", tree = arbol, sims = 1000)
setwd("SpxSp")
compare_phylofields(observed = "poligonos1_phylo_fields.txt",nulls = "poligonos1_random_phylofields_constrained.txt")

##########################################################
## Phylogenetic Generalized Least Squares (PGLS) method ##
##########################################################

#Load packages
require(ape)
require(caper)
require(geiger)

##Load the phylogentic data (in parenthetic format)
triato.tree <- read.tree("arbol_check")

#read the data with species name (SCINAME),	Phylogenetic Fields,	climatic variables and	PSV information as separate columns
triato.data.all <- read.table("PSVvsBio4_7_allsp.txt", header=TRUE, sep="")
triato.data.sig05 <- read.table("PSVvsBio4_7_sigsp_05.txt", header=TRUE, sep="")
triato.data.sig1 <- read.table("PSVvsBio4_7_sigsp_1.txt", header=TRUE, sep="")

#Create a "comparative data" object (caper package)
triatoall.compdata <- comparative.data(triato.tree,triato.data.all,SCINAME,vcv=T)
triatosig05.compdata <- comparative.data(triato.tree,triato.data.sig05,SCINAME,vcv=T)
triatosig1.compdata <- comparative.data(triato.tree,triato.data.sig1,SCINAME,vcv=T)

#PGLS method
triatoall1.pgls.Bio4 <- pgls(PSV1~Bio4mean,triatoall.compdata,lambda='ML')
triatoall1.pgls.Bio7 <- pgls(PSV1~Bio7mean,triatoall.compdata,lambda='ML')
triatoall05.pgls.Bio4 <- pgls(PSV05~Bio4mean,triatoall.compdata,lambda='ML')
triatoall05.pgls.Bio7 <- pgls(PSV05~Bio7mean,triatoall.compdata,lambda='ML')

triatosig1.pgls.Bio4 <- pgls(PSV1~Bio4mean,triatosig1.compdata,lambda='ML')
triatosig1.pgls.Bio7 <- pgls(PSV1~Bio7mean,triatosig1.compdata,lambda='ML')
triatosig05.pgls.Bio4 <- pgls(PSV05~Bio4mean,triatosig05.compdata,lambda='ML')
triatosig05.pgls.Bio7 <- pgls(PSV05~Bio7mean,triatosig05.compdata,lambda='ML')