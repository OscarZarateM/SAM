#S0_Data_preparation############################
#######################################################
library(sf)
library(tidyverse)
library(readxl)
library(sp)
library(terra)
library(iNEXT)
library(vegan)
library(DESeq2)
library(factoextra)
library(reshape2)
#######################################################
###Read no hits
OTU_t = read_excel("wanda23_ITS_c98.unite.SAM_samples.FunGuild.xlsx")
names(OTU_t)
##Taxonomy
Fungi_C = data.frame(OTU_t$OTU, OTU_t$taxonomy)
write(Fungi_C$OTU_t.taxonomy, "taxo.txt")
taxo = read.csv("taxo.txt", sep = ";", header = FALSE)
taxo

delete_before_double_underscore <- function(x) {
  sub(".*__", "", x)}
taxo = apply(taxo, 2, delete_before_double_underscore)
Fungi_C = as.data.frame(cbind(Fungi_C, taxo))

names(Fungi_C) = c("OTU", "Taxo","Kingdom", "Phylum", "Class", "Order", "Genus", "Species", "GS")
##Traits
Fungi_T = data.frame(OTU_t$OTU, OTU_t$trophicMode, OTU_t$guild, OTU_t$confidenceRanking)


OTU_t = select(OTU_t, -c("taxonomy", "taxonomicLevel", "trophicMode", "guild", "confidenceRanking", "notes","trait", "growthForm", "citationSource","SAM-KALJA-ITS"))

##Fungi traits


#### VT
##Remove Kalja
VT = read_excel("OTU_table.xlsx")
dim(VT)
reads_VT = sum(colSums(VT[,-1]))
VT = filter(VT, Sample !="SAM-KALJA-Wanda")
dim(VT)
reads_VT = sum(colSums(VT[,-1]))
VT = VT[!VT$Sample == "SAM-KALJA-Wanda",]
AMF_Y = t(VT[,-1])
colnames(AMF_Y) = VT$Sample
dim(AMF_Y)
##VT taxo
AMF_C = read_excel("wanda23_oscar_data_pioneer.v3.i97.a95.xlsx", sheet = "Sheet1", col_names = c("VT", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
AMF_C
## Vt guilds
VT_guilds = read_excel("guilds.xlsx")
AMF_T =left_join(AMF_C[,c("VT", "Family")], VT_guilds)
## OTU_t no hits table
Fungi_Y =OTU_t
## taxo no hits taxonomy
Fungi_C
## comm AMF community matrix
AMF_Y




#### X data
#### Soil
soil_t = read_excel("soil.xlsx", skip = 4)
### Spatial
points = st_read("SAM_points.gpkg")
lines = st_read("SAM_lines.gpkg")
### DEM
DEM = rast("DEM.tif")
#######################################################
#######################################################
#### II.Check #########################################
str(OTU_t)
str(soil_t)
str(points)
###Soil pHKCL
soil_t$pHKCl = as.numeric(soil_t$pHKCl)
### Soil No. = OTU_T$sample
soil_t$Sample = paste("SAM",soil_t$No.,"Wanda", sep = "-")
soil_t$No. = NULL
### points$layer
points$layer = gsub("SAM_", "", points$layer)
points$site = as.numeric(gsub("_points", "", points$layer))
points$layer =NULL
#######################################################
#######################################################
### III. Extract Elevation from DEMs ##################
elv = extract(DEM, points, ID =TRUE)
points = cbind(points, elv)
### IV. Join spatial-soil ############################
X = merge(soil_t, points, by.x = "Sample", by.y = "ID" )
str(X)
X$path = NULL
X$ID.1 = NULL
str(X)
### Distance st_distance###############################
d = data.frame()
for (point in X$Sample) {
  f = X[ which( X$Sample == point), "site"]
  p = X[ which( X$Sample == point), "geom"]
  #lines
  l =  subset(lines, layer == f )
  m_distance = as.numeric(st_distance(p, l))
  d =  rbind(d ,c(point, m_distance))
}
colnames(d) = c("Sample", "m_distance" )
###join with X
X = merge(X, d)
str(X)
X$m_distance = as.numeric(X$m_distance)
X
save(AMF_Y, AMF_T, AMF_C, X, file = "SAM_data.rda")
save(Fungi_Y, Fungi_C ,X, file = "Fungi_data.rda")
