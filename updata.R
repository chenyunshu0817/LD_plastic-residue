#library(BiodiversityR)
library(ggplot2)
library(vegan)
library(reshape2)
library(Hmisc)
library(plotrix)
library(phyloseq)
library(MASS)
library(bioDist) 
library(igraph)
library(car)
library(coin)
library(edgeR)
library(formatR)
library(gridExtra)
library(gplots)
library(indicspecies)
library(sciplot)
library(ape)
library(grid)
library(RVAideMemoire)
library(gridBase)
library(TukeyC)
library(corrplot)
library(userfriendlyscience)
library(caret)
library(multcompView)
library(ggpmisc)
library(agricolae)
library(ggrepel)
library(ggcor)
library(dplyr)
library(plyr)
library(stringr)
library(nlme)
library(lmerTest)
library(lme4)
library(Hmisc)
library(minpack.lm)
library(stats4)
library(spaa)
library(tidyr)
library("pulsar")
library(tidyverse)
library(ggClusterNet)


source("vennDia.R")
source("CorrDF.R")
source("plotOTU.R")
source("cor.mtest.R")
source("variance_functions.R")
source("maPalette.R")
source("triangle_shape.R")
source("star_shape.R")
source("zi_pi.r")
debuggingState(on=FALSE)
options(scipen=10)






###################################
##### Import and prepare data #####
###################################

#######################
##### Bacteria #####
#######################

##### Import Data #####
otu_B <- read.table("B_otutab.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
otu_B <- as.matrix(otu_B)
rownames(otu_B) <- paste("b", rownames(otu_B), sep="")

##### Import design file #####
design_B <- read.table("metadata.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
design_B <- design_B[c(which(design_B$Treatment1 == "CK"),which(design_B$Treatment1 == "PE"),which(design_B$Treatment1 == "DE")),]
design_B$Group <-factor(design_B$Treatment1,c("CK","PE","DE"))
design_B$Treatment1<-factor(design_B$Treatment1,c("CK","PE","DE"))
design_B$Season<-factor(design_B$Season,c("one","two"))
design_B$Site<-factor(design_B$Site,c("soil","root","end","film"))
design_B$Concentration <- factor(design_B$Concentration,c("0","800","4000"))
design_B$Block <- factor(design_B$Block,c("F1.1","F1.2","F1.3","F1.4"))

str(design_B)

##### Import Taxonomy #####
tax_B <- read.table("B_taxonomy.txt", row.names=1, sep="\t", header=T ,stringsAsFactors=F,quote="")
rownames(tax_B) <- paste("b",rownames(tax_B),sep="")


# create separate taxonomy label specifying classes of Proteobacteria
tax_B$labels <- tax_B$Phylum
tax_B[ rownames(tax_B)[tax_B$Class=="Alphaproteobacteria" ], ]$labels <- "Alphaproteobacteria"
tax_B[ rownames(tax_B)[tax_B$Class=="Betaproteobacteria" ], ]$labels <- "Betaproteobacteria"
tax_B[ rownames(tax_B)[tax_B$Class=="Gammaproteobacteria" ], ]$labels <- "Gammaproteobacteria"
tax_B[ rownames(tax_B)[tax_B$Class=="Deltaproteobacteria" ], ]$labels <- "Deltaproteobacteria"
table(tax_B$labels)

### Defining bOTU colors by phylum (using the taxonomy file)
tax_B$cols <- tax_B$labels
table(tax_B$cols)

##### Store Cyanobacteria and mitrochrondia sequeces #####
unique(tax_B$Kingdom)
table(tax_B$Kingdom)
r1 <- rownames(tax_B[tax_B$Kingdom=="Archaea",])

table(tax_B$Phylum)
r2 <- rownames(tax_B[tax_B$Phylum=="Cyanobacteria",])

otus_remove_B <- c(r1,r2)

## Remove these from otu table, tax table   
otu_B <- otu_B[-which(rownames(otu_B) %in% otus_remove_B),]
otu_B <- otu_B[,rownames(design_B)]
otu_B <- otu_B[rowSums(otu_B) > 0,]
nrow(otu_B)

tax_B <- tax_B[rownames(otu_B),]

design_B <- droplevels(design_B[rownames(design_B) %in% colnames(otu_B),])
design_B <- design_B[colnames(otu_B),]

dim(otu_B)
dim(tax_B)
dim(design_B)



###### B sequence and OTU counts ######
sum(colSums(otu_B))
sort(colSums(otu_B))
median(colSums(otu_B))

nrow(tax_B)
table(tax_B$Kingdom)

## Order taxonmy file by OTU
otu_order_B <- match(rownames(otu_B), rownames(tax_B))
tax_B <- tax_B[otu_order_B,]


##### Define sample types #####
onesamples <- rownames(design_B)[which(design_B$Season == "one")]
twosamples <- rownames(design_B)[which(design_B$Season == "two")]

soilsamples <- rownames(design_B)[which(design_B$Site == "soil")]
rootsamples <- rownames(design_B)[which(design_B$Site == "root")]
particlesamples <- rownames(design_B)[which(design_B$Site == "film")]
endsamples <- rownames(design_B)[which(design_B$Site == "end")]

##### Total number of soil / root OTUs #####
soil_B_otu <- otu_B[,soilsamples]
soil_B_otu <- soil_B_otu[rowSums(soil_B_otu) > 0,]
nrow(soil_B_otu)

root_B_otu <- otu_B[,rootsamples]
root_B_otu <- root_B_otu[rowSums(root_B_otu) > 0,]
nrow(root_B_otu)

end_B_otu <- otu_B[,endsamples]
end_B_otu <- end_B_otu[rowSums(end_B_otu) > 0,]
nrow(end_B_otu)

particle_B_otu <- otu_B[,particlesamples]
particle_B_otu <- particle_B_otu[rowSums(particle_B_otu) > 0,]
nrow(particle_B_otu)



#######################
##### Fungi #####
#######################

#### Import OTU table #####
otu_F <- read.table("F_otutab.txt",row.names=1,sep="\t",header=T,blank.lines.skip=F,check.names=F)
otu_F <- as.matrix(otu_F)
rownames(otu_F) <- paste("f",rownames(otu_F),sep="")


##### Import design file #####
design_F <- read.table("metadata.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
design_F <- design_F[c(which(design_F$Treatment1 == "CK"),which(design_F$Treatment1 == "PE"),which(design_F$Treatment1 == "DE")),]
design_F$Group <-factor(design_F$Treatment1,c("CK","PE","DE"))
design_F$Treatment1<-factor(design_F$Treatment1,c("CK","PE","DE"))
design_F$Season<-factor(design_F$Season,c("one","two"))
design_F$Site<-factor(design_F$Site,c("soil","root","end","film"))
design_F$Concentration <- factor(design_F$Concentration,c("0","800","4000"))
design_F$Block <- factor(design_F$Block,c("F1.1","F1.2","F1.3","F1.4"))

str(design_F)


##### Import Taxonomy #####
tax_F <- read.table("F_taxonomy.txt", row.names=1, sep="\t", header=T ,stringsAsFactors=F,quote="")
rownames(tax_F) <- paste("f",rownames(tax_F),sep="")
tax_F$labels <- tax_F$Phylum

tax_F[ rownames(tax_F)[tax_F$Class=="Sordariomycetes" ], ]$labels <- "Sordariomycetes"
tax_F[ rownames(tax_F)[tax_F$Class=="Eurotiomycetes" ], ]$labels <- "Eurotiomycetes"
tax_F[ rownames(tax_F)[tax_F$Class=="Dothideomycetes" ], ]$labels <- "Dothideomycetes"
tax_F[ rownames(tax_F)[tax_F$Class=="Orbiliomycetes" ], ]$labels <- "Orbiliomycetes"
tax_F[ rownames(tax_F)[tax_F$Class=="Agaricomycetes" ], ]$labels <- "Agaricomycetes"
tax_F[ rownames(tax_F)[tax_F$Class=="Tremellomycetes" ], ]$labels <- "Tremellomycetes"


table(tax_F$labels)


### Defining fOTU colors by phylum (using the taxonomy file)

tax_F$cols <- tax_F$labels
table(tax_F$cols)



##### Store Plantae, Protists and sequences of unknown origin and remove #####
#OTUs classified in kingdom Plantae or Protista are removed, as well as those whose kingdom is unassigned

unique(tax_F$Kingdom)
table(tax_F$Kingdom)

r3 <- rownames(tax_F[tax_F$Phylum == "Unassigned",])

otus_remove_F <- c(r3)

##### Calculate % sequences removed #####
otu_F <- otu_F[-which(rownames(otu_F) %in% otus_remove_F),]
otu_F <- otu_F[,rownames(design_F)]
otu_F <- otu_F[rowSums(otu_F) > 0,]
nrow(otu_F)

tax_F <- tax_F[rownames(otu_F),]

design_F <- droplevels(design_F[rownames(design_F) %in% colnames(otu_F),])
design_F <- design_F[colnames(otu_F),]

dim(otu_F)
dim(tax_F)
dim(design_F)


###### F sequence and OTU counts ######
sum(colSums(otu_F))
sort(colSums(otu_F))
median(colSums(otu_F))

nrow(tax_F)
table(tax_F$Kingdom)

## Order taxonmy file by OTU
otu_order_F <- match(rownames(otu_F), rownames(tax_F))
tax_F <- tax_F[otu_order_F, ]

##### Define sample types #####
##### Define sample types #####
onesamples <- rownames(design_F)[which(design_F$Season == "one")]
twosamples <- rownames(design_F)[which(design_F$Season == "two")]

soilsamples <- rownames(design_F)[which(design_F$Site == "soil")]
rootsamples <- rownames(design_F)[which(design_F$Site == "root")]
particlesamples <- rownames(design_F)[which(design_F$Site == "film")]
endsamples <- rownames(design_F)[which(design_F$Site == "end")]

##### Total number of soil / root OTUs #####
soil_F_otu <- otu_F[,soilsamples]
soil_F_otu <- soil_F_otu[rowSums(soil_F_otu) > 0,]
nrow(soil_F_otu)

root_F_otu <- otu_F[,rootsamples]
root_F_otu <- root_F_otu[rowSums(root_F_otu) > 0,]
nrow(root_F_otu)

end_F_otu <- otu_F[,endsamples]
end_F_otu <- end_F_otu[rowSums(end_F_otu) > 0,]
nrow(end_F_otu)

particle_F_otu <- otu_F[,particlesamples]
particle_F_otu <- particle_F_otu[rowSums(particle_F_otu) > 0,]
nrow(particle_F_otu)




#######################
##### Protist #####
#######################

#### Import OTU table #####
otu_P <- read.table("P_otutab.txt",row.names=1,sep="\t",header=T,blank.lines.skip=F,check.names=F)
otu_P <- as.matrix(otu_P)
rownames(otu_P) <- paste("p",rownames(otu_P),sep="")


##### Import design file #####
design_P <- read.table("P_metadata.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
design_P <- design_P[c(which(design_P$Treatment1 == "CK"),which(design_P$Treatment1 == "PE"),which(design_P$Treatment1 == "DE")),]
design_P$Group <-factor(design_P$Treatment1,c("CK","PE","DE"))
design_P$Treatment1<-factor(design_P$Treatment1,c("CK","PE","DE"))
design_P$Season<-factor(design_P$Season,c("one","two"))
design_P$Site<-factor(design_P$Site,c("soil","root","end","film"))
design_P$Concentration <- factor(design_P$Concentration,c("0","800","4000"))
design_P$Block <- factor(design_P$Block,c("F1.1","F1.2","F1.3","F1.4"))

str(design_P)


##### Import Taxonomy #####
tax_P <- read.table("P_taxonomy.txt", row.names=1, sep="\t", header=T ,stringsAsFactors=F,quote="")
rownames(tax_P) <- paste("p",rownames(tax_P),sep="")
tax_P$labels <- tax_P$Phylum

tax_P[ rownames(tax_P)[tax_P$Class=="Metazoa" ], ]$labels <- "Metazoa"
tax_P[ rownames(tax_P)[tax_P$Class=="Choanoflagellida" ], ]$labels <- "Choanoflagellida"
tax_P[ rownames(tax_P)[tax_P$Class=="Mesomycetozoa" ], ]$labels <- "Mesomycetozoa"
table(tax_P$labels)

table(tax_P$labels)


### Defining OTU colors by phylum (using the taxonomy file)
tax_P$cols <- tax_P$labels
table(tax_P$cols)

unique(tax_P$Kingdom)
table(tax_P$Kingdom)

# unassignphy_otu <- otu_P[rownames(otu_P) %in% r6,]; dim(unassignphy_otu); length(r6)
# sum(colSums(unassignphy_otu))/sum(colSums(otu_P))*100
otu_P <- otu_P[,rownames(design_P)]
otu_P <- otu_P[rowSums(otu_P) > 0,]
nrow(otu_P)

tax_P <- tax_P[rownames(otu_P),]

design_P <- droplevels(design_P[rownames(design_P) %in% colnames(otu_P),])
design_P <- design_P[colnames(otu_P),]

dim(otu_P)
dim(tax_P)
dim(design_P)

###### F sequence and OTU counts ######
sum(colSums(otu_P))
sort(colSums(otu_P))
median(colSums(otu_P))

nrow(tax_P)
table(tax_P$Kingdom)

## Order taxonmy file by OTU
otu_order_P <- match(rownames(otu_P), rownames(tax_P))
tax_P <- tax_P[otu_order_P, ]

##### Define sample types #####
ponesamples <- rownames(design_P)[which(design_P$Season == "one")]
ptwosamples <- rownames(design_P)[which(design_P$Season == "two")]

psoilsamples <- rownames(design_P)[which(design_P$Site == "soil")]
prootsamples <- rownames(design_P)[which(design_P$Site == "root")]
pparticlesamples <- rownames(design_P)[which(design_P$Site == "film")]
pendsamples <- rownames(design_P)[which(design_P$Site == "end")]

##### Total number of soil / root OTUs #####
soil_P_otu <- otu_P[,psoilsamples]
soil_P_otu <- soil_P_otu[rowSums(soil_P_otu) > 0,]
nrow(soil_P_otu)

root_P_otu <- otu_P[,prootsamples]
root_P_otu <- root_P_otu[rowSums(root_P_otu) > 0,]
nrow(root_P_otu)

end_P_otu <- otu_P[,pendsamples]
end_P_otu <- end_P_otu[rowSums(end_P_otu) > 0,]
nrow(end_P_otu)

particle_P_otu <- otu_P[,pparticlesamples]
particle_P_otu <- particle_P_otu[rowSums(particle_P_otu) > 0,]
nrow(particle_P_otu)


##### Alpha #####
set.seed(619)
cp_otu_B <- as.data.frame(t(rrarefy(t(otu_B), min(colSums(otu_B)))))
cp_B_otu_shannon <- vegan::diversity(t(otu_B), index = "shannon")
design_B$Shannon <- cp_B_otu_shannon
cp_B_otu_simpson <- vegan::diversity(t(otu_B), index = "simpson")
design_B$Simpson <- cp_B_otu_simpson
chao1 <- t(estimateR(t(cp_otu_B)))
chao1 <- as.data.frame(chao1)
design_B$chao1 <- chao1$S.chao1
design_B$ACE <- chao1$S.ACE


# line_alpha_bacteria <- lmer(Shannon~Site + Type + Concentration + Treatment1 + Treatment2 + Treatment3+ (1|Block), REML=TRUE,data=design_B)
# Anova(line_alpha_bacteria,type=2,ddf="lme4",test="F")

set.seed(619)
cp_otu_F <- as.data.frame(t(rrarefy(t(otu_F), min(colSums(otu_F)))))
cp_F_otu_shannon <- vegan::diversity(t(otu_F), index = "shannon")
design_F$Shannon <- cp_F_otu_shannon
cp_F_otu_simpson <- vegan::diversity(t(otu_F), index = "simpson")
design_F$Simpson <- cp_F_otu_simpson
chao1 <- t(estimateR(t(cp_otu_F)))
chao1 <- as.data.frame(chao1)
design_F$chao1 <- chao1$S.chao1
design_F$ACE <- chao1$S.ACE


# line_alpha_fungi <- lmer(Shannon~Site + Type + Concentration + Treatment1 + Treatment2 + Treatment3+ (1|Block), REML=TRUE,data=design_F)
# Anova(line_alpha_fungi,type=2,ddf="lme4",test="F")


set.seed(619)
cp_otu_P <- as.data.frame(t(rrarefy(t(otu_P), min(colSums(otu_P)))))
cp_P_otu_shannon <- vegan::diversity(t(otu_P), index = "shannon")
design_P$Shannon <- cp_P_otu_shannon
cp_P_otu_simpson <- vegan::diversity(t(otu_P), index = "simpson")
design_P$Simpson <- cp_P_otu_simpson
chao1 <- t(estimateR(t(cp_otu_P)))
chao1 <- as.data.frame(chao1)
design_P$chao1 <- chao1$S.chao1
design_P$ACE <- chao1$S.ACE



# line_alpha_protist <- lmer(Shannon~Site + Type + Concentration + Treatment1 + Treatment2 + Treatment3+ (1|Block), REML=TRUE,data=design_P)
# Anova(line_alpha_protist,type=2,ddf="lme4",test="F")



#####The beta diversity analysis of Bacteria #####
## Apply sequence count threshold to bacteria soil community and TMM normalize counts

## Bacteria one year 
otu_B_one <- otu_B[, onesamples ]
otu_B_one <- otu_B_one[rowSums(otu_B_one) > 0,]
dim(otu_B_one)

tax_B_one <- tax_B[rownames(otu_B_one),]
design_B_one <- droplevels(design_B[onesamples,])

edgeR_B_one <- DGEList(counts=otu_B_one, 
                       group=design_B_one$Site,
                       genes=tax_B_one)

edgeR_B_one <- calcNormFactors(edgeR_B_one)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_one_B <- cpm(edgeR_B_one, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk soil bacteria community into phyloseq objects
## for further analysis
phy_B_one <- otu_table(otu_norm_one_B,taxa_are_rows=T)
phy_tax_one_B <-tax_table(as.matrix(tax_B_one))
phy_design_one_B <- sample_data(design_B_one)
physeq_one_norm_B <- phyloseq(phy_B_one,phy_tax_one_B,phy_design_one_B)
sample_data(physeq_one_norm_B)$Site<- factor(sample_data(physeq_one_norm_B)$Site,levels=c("soil","root","end","film"))


## Create bray-curtis dissimiliartiy matrix
all_dis_B_one <- vegdist(t(otu_table(physeq_one_norm_B)),method="bray")

paov_all_B_one <- adonis2(all_dis_B_one ~Block+Site+Concentration+Type+Treatment1, data=design_B_one, permutations=9999)

paov_all_B_one  

pcoa_norm_B_one <- ordinate(physeq_one_norm_B,"PCoA","bray")
pcoa_B_one <- plot_ordination(physeq_one_norm_B, pcoa_norm_B_one, type="sites", color="Site", shape="Treatment1")
pcoa_B_one <- pcoa_B_one+
  geom_point(size=4)+
  scale_shape_manual(values = c(19,17,18))+ 
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("Bacteria one community")


## Alpha
CK  <- design_B[design_B$Treatment1=="CK", 9 ]
PE  <- design_B[design_B$Treatment1=="PE", 9 ]
DE  <- design_B[design_B$Treatment1=="DE", 9 ]


my_comparisons = list( c('CK', 'PE'), c("CK", "DE"), c("PE", "DE"))

Shannon_B <- ggplot(design_B, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=Shannon, fill =Treatment1))+geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+labs(x="Type", y="Bacteria Shannon index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")

CK  <- design_B[design_B$Treatment1=="CK", 11 ]
PE  <- design_B[design_B$Treatment1=="PE", 11 ]
DE  <- design_B[design_B$Treatment1=="DE", 11 ]

Chao1_B <- ggplot(design_B, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=chao1, fill =Treatment1))+
  geom_violin(trim=FALSE,color="white") +geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Bacteria Chao1 index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")

grid.newpage()
grid.arrange(Shannon_B,Chao1_B,ncol = 2)



Shannon_B <- ggplot(design_B, aes(x=factor(Site,levels = c("soil","root","end","film")), y=Shannon, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Bacteria Shannon index")

Chao1_B <- ggplot(design_B, aes(x=factor(Site,levels = c("soil","root","end","film")), y=chao1, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Bacteria Chao1 index")

grid.newpage()
grid.arrange(Shannon_B,Chao1_B,ncol = 2)




Shannon_B_one <- ggplot(design_B_one, aes(x=factor(Site,levels = c("soil","root","end","film")), y=Shannon, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Bacteria one soil Shannon index")

Chao1_B_one <- ggplot(design_B_one, aes(x=factor(Site,levels = c("soil","root","end","film")), y=chao1, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Bacteria one soil Chao1 index")

grid.newpage()
grid.arrange(Shannon_B_one,Chao1_B_one,ncol = 2)


## Bacteria two year 
otu_B_two <- otu_B[, twosamples ]
otu_B_two <- otu_B_two[rowSums(otu_B_two) > 0,]
dim(otu_B_two)

tax_B_two <- tax_B[rownames(otu_B_two),]
design_B_two <- droplevels(design_B[twosamples,])

edgeR_B_two <- DGEList(counts=otu_B_two, 
                       group=design_B_two$Treatment1,
                       genes=tax_B_two)

edgeR_B_two <- calcNormFactors(edgeR_B_two)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_two_B <- cpm(edgeR_B_two, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk soil bacteria community into phyloseq objects
## for further analysis
phy_B_two <- otu_table(otu_norm_two_B,taxa_are_rows=T)
phy_tax_two_B <-tax_table(as.matrix(tax_B_two))
phy_design_two_B <- sample_data(design_B_two)
physeq_two_norm_B <- phyloseq(phy_B_two,phy_tax_two_B,phy_design_two_B)
sample_data(physeq_two_norm_B)$Type<- factor(sample_data(physeq_two_norm_B)$Type,levels=c("PBAT","PE"))


## Create bray-curtis dissimiliartiy matrix
all_dis_B_two <- vegdist(t(otu_table(physeq_two_norm_B)),method="bray")

paov_all_B_two <- adonis2(all_dis_B_two ~Block+Site+Concentration+Type+Treatment1, data=design_B_two, permutations=9999)

paov_all_B_two  

pcoa_norm_B_two <- ordinate(physeq_two_norm_B,"PCoA","bray")
pcoa_B_two <- plot_ordination(physeq_two_norm_B, pcoa_norm_B_two, type="sites", color="Site", shape="Treatment1")
pcoa_B_two <- pcoa_B_two+
  geom_point(size=4)+
  scale_shape_manual(values = c(19,17,18))+ 
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("Bacteria two community")


## Alpha
Shannon_B_two <- ggplot(design_B_two, aes(x=factor(Site,levels = c("soil","root","end","film")), y=Shannon, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Bacteria two soil Shannon index")

Chao1_B_two <- ggplot(design_B_two, aes(x=factor(Site,levels = c("soil","root","end","film")), y=chao1, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Bacteria two soil Chao1 index")

grid.newpage()
grid.arrange(Shannon_B_two,Chao1_B_two,ncol = 2)



## Fungi one year 
otu_F_one <- otu_F[, onesamples ]
otu_F_one <- otu_F_one[rowSums(otu_F_one) > 0,]
dim(otu_F_one)

tax_F_one <- tax_F[rownames(otu_F_one),]
design_F_one <- droplevels(design_F[onesamples,])

edgeR_F_one <- DGEList(counts=otu_F_one, 
                       group=design_F_one$Type,
                       genes=tax_F_one)

edgeR_F_one <- calcNormFactors(edgeR_F_one)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_one_F <- cpm(edgeR_F_one, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk soil bacteria community into phyloseq objects
## for further analysis
phy_F_one <- otu_table(otu_norm_one_F,taxa_are_rows=T)
phy_tax_one_F <-tax_table(as.matrix(tax_F_one))
phy_design_one_F <- sample_data(design_F_one)
physeq_one_norm_F <- phyloseq(phy_F_one,phy_tax_one_F,phy_design_one_F)
sample_data(physeq_one_norm_F)$Type<- factor(sample_data(physeq_one_norm_F)$Type,levels=c("PBAT","PE"))


## Create bray-curtis dissimiliartiy matrix
all_dis_F_one <- vegdist(t(otu_table(physeq_one_norm_F)),method="bray")

paov_all_F_one <- adonis2(all_dis_F_one ~Block+Site+Concentration+Type+Treatment1, data=design_F_one, permutations=9999)

paov_all_F_one  

pcoa_norm_F_one <- ordinate(physeq_one_norm_F,"PCoA","bray")
pcoa_F_one <- plot_ordination(physeq_one_norm_F, pcoa_norm_F_one, type="sites", color="Site", shape="Treatment1")
pcoa_F_one <- pcoa_F_one+
  geom_point(size=4)+
  scale_shape_manual(values = c(19,17,18))+ 
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("Fungi one community")


## Alpha
CK  <- design_F[design_F$Treatment1=="CK", 9 ]
PE  <- design_F[design_F$Treatment1=="PE", 9 ]
DE  <- design_F[design_F$Treatment1=="DE", 9 ]


my_comparisons = list( c('CK', 'PE'), c("CK", "DE"), c("PE", "DE"))

Shannon_F <- ggplot(design_F, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=Shannon, fill =Treatment1))+geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+labs(x="Type", y="Fungi Shannon index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")

CK  <- design_F[design_F$Treatment1=="CK", 11 ]
PE  <- design_F[design_F$Treatment1=="PE", 11 ]
DE  <- design_F[design_F$Treatment1=="DE", 11 ]

Chao1_F <- ggplot(design_F, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=chao1, fill =Treatment1))+geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Fungi Chao1 index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")

grid.newpage()
grid.arrange(Shannon_F,Chao1_F,ncol = 2)




Shannon_F <- ggplot(design_F, aes(x=factor(Site,levels = c("soil","root","end","film")), y=Shannon, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Fungi one soil Shannon index")

Chao1_F <- ggplot(design_F, aes(x=factor(Site,levels = c("soil","root","end","film")), y=chao1, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Fungi one soil Chao1 index")

grid.newpage()
grid.arrange(Shannon_F,Chao1_F,ncol = 2)


Shannon_F_one <- ggplot(design_F_one, aes(x=factor(Site,levels = c("soil","root","end","film")), y=Shannon, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Fungi one soil Shannon index")

Chao1_F_one <- ggplot(design_F_one, aes(x=factor(Site,levels = c("soil","root","end","film")), y=chao1, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Fungi one soil Chao1 index")

grid.newpage()
grid.arrange(Shannon_F_one,Chao1_F_one,ncol = 2)


## Fungi two year 
otu_F_two <- otu_F[, twosamples ]
otu_F_two <- otu_F_two[rowSums(otu_F_two) > 0,]
dim(otu_F_two)

tax_F_two <- tax_F[rownames(otu_F_two),]
design_F_two <- droplevels(design_F[twosamples,])

edgeR_F_two <- DGEList(counts=otu_F_two, 
                       group=design_F_two$Treatment1,
                       genes=tax_F_two)

edgeR_F_two <- calcNormFactors(edgeR_F_two)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_two_F <- cpm(edgeR_F_two, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk soil bacteria community into phyloseq objects
## for further analysis
phy_F_two <- otu_table(otu_norm_two_F,taxa_are_rows=T)
phy_tax_two_F <-tax_table(as.matrix(tax_F_two))
phy_design_two_F <- sample_data(design_F_two)
physeq_two_norm_F <- phyloseq(phy_F_two,phy_tax_two_F,phy_design_two_F)
sample_data(physeq_two_norm_F)$Site<- factor(sample_data(physeq_two_norm_F)$Site,levels=c("soil","root","end","film"))


## Create bray-curtis dissimiliartiy matrix
all_dis_F_two <- vegdist(t(otu_table(physeq_two_norm_F)),method="bray")

paov_all_F_two <- adonis2(all_dis_F_two ~Block+Site+Concentration+Type+Treatment1, data=design_F_two, permutations=9999)

paov_all_F_two  

pcoa_norm_F_two <- ordinate(physeq_two_norm_F,"PCoA","bray")
pcoa_F_two <- plot_ordination(physeq_two_norm_F, pcoa_norm_F_two, type="sites",color="Site", shape="Treatment1")
pcoa_F_two <- pcoa_F_two+
  geom_point(size=4)+
  scale_shape_manual(values = c(19,17,18))+ 
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("Fungi two community")


## Alpha
Shannon_F_two <- ggplot(design_F_two, aes(x=factor(Site,levels = c("soil","root","end","film")), y=Shannon, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Fungi two soil Shannon index")

Chao1_F_two <- ggplot(design_F_two, aes(x=factor(Site,levels = c("soil","root","end","film")), y=chao1, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Fungi two soil Chao1 index")

grid.newpage()
grid.arrange(Shannon_F_two,Chao1_F_two,ncol = 2)



## Protist one year 
otu_P_one <- otu_P[, ponesamples ]
otu_P_one <- otu_P_one[rowSums(otu_P_one) > 0,]
dim(otu_P_one)

tax_P_one <- tax_P[rownames(otu_P_one),]
design_P_one <- droplevels(design_P[ponesamples,])

edgeR_P_one <- DGEList(counts=otu_P_one, 
                       group=design_P_one$Site,
                       genes=tax_P_one)

edgeR_P_one <- calcNormFactors(edgeR_P_one)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_one_P <- cpm(edgeR_P_one, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk soil bacteria community into phyloseq objects
## for further analysis
phy_P_one <- otu_table(otu_norm_one_P,taxa_are_rows=T)
phy_tax_one_P <-tax_table(as.matrix(tax_P_one))
phy_design_one_P <- sample_data(design_P_one)
physeq_one_norm_P <- phyloseq(phy_P_one,phy_tax_one_P,phy_design_one_P)
sample_data(physeq_one_norm_P)$Site<- factor(sample_data(physeq_one_norm_P)$Site,levels=c("soil","root","film"))


## Create bray-curtis dissimiliartiy matrix
all_dis_P_one <- vegdist(t(otu_table(physeq_one_norm_P)),method="bray")

paov_all_P_one <- adonis2(all_dis_P_one ~Block+Site+Concentration+Type+Treatment1, data=design_P_one, permutations=9999)

paov_all_P_one  

pcoa_norm_P_one <- ordinate(physeq_one_norm_P,"PCoA","bray")
pcoa_P_one <- plot_ordination(physeq_one_norm_P, pcoa_norm_P_one, type="sites", color="Site", shape="Treatment1")
pcoa_P_one <- pcoa_P_one+
  geom_point(size=4)+
  scale_shape_manual(values = c(19,17,18))+ 
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("Protist one community")


## Alpha
Shannon_P_one <- ggplot(design_P_one, aes(x=factor(Site,levels = c("soil","root","film")), y=Shannon, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Fungi one soil Shannon index")

Chao1_P_one <- ggplot(design_P_one, aes(x=factor(Site,levels = c("soil","root","film")), y=chao1, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Fungi one soil Chao1 index")

grid.newpage()
grid.arrange(Shannon_P_one,Chao1_P_one,ncol = 2)


## Fungi two year 
otu_P_two <- otu_P[, ptwosamples ]
otu_P_two <- otu_P_two[rowSums(otu_P_two) > 0,]
dim(otu_P_two)

tax_P_two <- tax_P[rownames(otu_P_two),]
design_P_two <- droplevels(design_P[ptwosamples,])

edgeR_P_two <- DGEList(counts=otu_P_two, 
                       group=design_P_two$Site,
                       genes=tax_P_two)

edgeR_P_two <- calcNormFactors(edgeR_P_two)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_two_P <- cpm(edgeR_P_two, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk soil bacteria community into phyloseq objects
## for further analysis
phy_P_two <- otu_table(otu_norm_two_P,taxa_are_rows=T)
phy_tax_two_P <-tax_table(as.matrix(tax_P_two))
phy_design_two_P <- sample_data(design_P_two)
physeq_two_norm_P <- phyloseq(phy_P_two,phy_tax_two_P,phy_design_two_P)
sample_data(physeq_two_norm_P)$Site<- factor(sample_data(physeq_two_norm_P)$Site,levels=c("soil","root","end","film"))


## Create bray-curtis dissimiliartiy matrix
all_dis_P_two <- vegdist(t(otu_table(physeq_two_norm_P)),method="bray")

paov_all_P_two <- adonis2(all_dis_P_two ~Block+Site+Treatment1+Concentration+Type, data=design_P_two, permutations=9999)

paov_all_P_two  

pcoa_norm_P_two <- ordinate(physeq_two_norm_P,"PCoA","bray")
pcoa_P_two <- plot_ordination(physeq_two_norm_P, pcoa_norm_P_two, type="sites", color="Site", shape="Treatment1")
pcoa_P_two <- pcoa_P_two+
  geom_point(size=4)+
  scale_shape_manual(values = c(19,17,18))+ 
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("Protist two community")


## Alpha
CK  <- design_P[design_P$Treatment1=="CK", 9 ]
PE  <- design_P[design_P$Treatment1=="PE", 9 ]
DE  <- design_P[design_P$Treatment1=="DE", 9 ]


my_comparisons = list( c('CK', 'PE'), c("CK", "DE"), c("PE", "DE"))

Shannon_P <- ggplot(design_P, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=Shannon, fill =Treatment1))+geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+labs(x="Type", y="Protist Shannon index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")

CK  <- design_P[design_P$Treatment1=="CK", 11 ]
PE  <- design_P[design_P$Treatment1=="PE", 11 ]
DE  <- design_P[design_P$Treatment1=="DE", 11 ]

Chao1_P <- ggplot(design_P, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=chao1, fill =Treatment1))+geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Protist Chao1 index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")

grid.newpage()
grid.arrange(Shannon_P,Chao1_P,ncol = 2)




Shannon_P_two <- ggplot(design_P_two, aes(x=factor(Site,levels = c("soil","root","end","film")), y=Shannon, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Fungi two soil Shannon index")

Chao1_P_two <- ggplot(design_P_two, aes(x=factor(Site,levels = c("soil","root","end","film")), y=chao1, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Type", y="Fungi two soil Chao1 index")

grid.newpage()
grid.arrange(Shannon_P_two,Chao1_P_two,ncol = 2)

grid.newpage()
grid.arrange(Shannon_B_one,Chao1_B_one,Shannon_B_two,Chao1_B_two,Shannon_F_one,Chao1_F_one,Shannon_F_two,Chao1_F_two,Shannon_P_one,Chao1_P_one,Shannon_P_two,Chao1_P_two,ncol = 4)

grid.newpage()
grid.arrange(pcoa_B_one,pcoa_F_one,pcoa_P_one,pcoa_B_two,pcoa_F_two,pcoa_P_two,ncol = 3)

grid.newpage()
grid.arrange(Shannon_B,Shannon_F,Shannon_P,Chao1_B,Chao1_F,Chao1_P,ncol = 3)


##### Supplementary Figure S3: phyla relative abundance plots #####

##### Bacteria

## Express B OTU counts as relative abunance percent
otu_B_RA <- t(t(otu_B)/colSums(otu_B)) * 100
colSums(otu_B_RA)
nrow(otu_B_RA)

## Get names of bacteria phyla present (use 'labels' as this specifies class within Proteobacteria)
PHYLAnames_B <- names(sort(table(tax_B[,"labels"]), decr=T))
length(PHYLAnames_B)
sort(table(tax_B[,"labels"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(otu_B_RA)
for (i in PHYLAnames_B){
  x <- array(colSums(otu_B_RA[rownames(tax_B)[which(tax_B$labels == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(PHYLAnames_B)
colnames(y) <- paste(colnames(otu_B_RA))
PHYLUM_mat_B <- y
PHYLUM_mat_B[,1:5]
colSums(PHYLUM_mat_B)
PHYLUM_mat_B_mean <- sort(apply(PHYLUM_mat_B,1,mean),decr=T)
PHYLUM_mat_B <- PHYLUM_mat_B[names(PHYLUM_mat_B_mean),]

## Express B OTU counts as relative abunance percent
design_B <- unite(design_B, Treatment2, c(Site, Treatment1))
design_B <- unite(design_B, Treatment2, c(Season,Treatment2))

one_soil_CK<- rownames(design_B[design_B$Treatment2=="one_soil_CK",])
one_root_CK<- rownames(design_B[design_B$Treatment2=="one_root_CK",])
one_end_CK<- rownames(design_B[design_B$Treatment2=="one_end_CK",])
two_soil_CK<- rownames(design_B[design_B$Treatment2=="two_soil_CK",])
two_root_CK<- rownames(design_B[design_B$Treatment2=="two_root_CK",])
two_end_CK<- rownames(design_B[design_B$Treatment2=="two_end_CK",])
one_soil_PE<- rownames(design_B[design_B$Treatment2=="one_soil_PE",])
one_root_PE<- rownames(design_B[design_B$Treatment2=="one_root_PE",])
one_end_PE<- rownames(design_B[design_B$Treatment2=="one_end_PE",])
one_film_PE <- rownames(design_B[design_B$Treatment2=="one_film_PE",])
two_soil_PE <- rownames(design_B[design_B$Treatment2=="two_soil_PE",])
two_root_PE <- rownames(design_B[design_B$Treatment2=="two_root_PE",])
two_end_PE <- rownames(design_B[design_B$Treatment2=="two_end_PE",])
two_film_PE <- rownames(design_B[design_B$Treatment2=="two_film_PE",])
one_soil_DE <- rownames(design_B[design_B$Treatment2=="one_soil_DE",])
one_root_DE <- rownames(design_B[design_B$Treatment2=="one_root_DE",])
one_end_DE <- rownames(design_B[design_B$Treatment2=="one_end_DE",])
one_film_DE <- rownames(design_B[design_B$Treatment2=="one_film_DE",])
two_soil_DE <- rownames(design_B[design_B$Treatment2=="two_soil_DE",])
two_root_DE <- rownames(design_B[design_B$Treatment2=="two_root_DE",])
two_end_DE <- rownames(design_B[design_B$Treatment2=="two_end_DE",])
two_film_DE <- rownames(design_B[design_B$Treatment2=="two_film_DE",])




## Express B OTU counts as relative abunance percent
design_B <- unite(design_B, Treatment2, c(Site, Season))

one_soil <- rownames(design_B[design_B$Treatment2=="soil_one",])
one_root <- rownames(design_B[design_B$Treatment2=="root_one",])
one_end <- rownames(design_B[design_B$Treatment2=="end_one",])
two_soil <- rownames(design_B[design_B$Treatment2=="soil_two",])
two_root <- rownames(design_B[design_B$Treatment2=="root_two",])
two_end <- rownames(design_B[design_B$Treatment2=="end_two",])

one_film <- rownames(design_B[design_B$Treatment2=="film_one",])
two_film <- rownames(design_B[design_B$Treatment2=="film_two",])



## Soil phyla abundances
PHYLUM_mat_B_film <- PHYLUM_mat_B[,two_film]
colSums(PHYLUM_mat_B_film)
PHYLUM_mat_B_film_mean <- sort(apply(PHYLUM_mat_B_film,1,mean),decr=T)
PHYLUM_mat_B_film <- PHYLUM_mat_B_film[names(PHYLUM_mat_B_film_mean),]
PHYLUM_mat_B_film_se <- apply(PHYLUM_mat_B_film,1,se)[names(PHYLUM_mat_B_film_mean)]

length(PHYLUM_mat_B_film_mean[PHYLUM_mat_B_film_mean > 0])

design_B_one_film <- design_B[two_film,]


PHYLUM_mat_B_film <- t(PHYLUM_mat_B_film)
PHYLUM_mat_B_film <- as.data.frame(PHYLUM_mat_B_film)
PHYLUM_mat_B_film$Treatment1 <- design_B_one_film$Treatment1

for(i in PHYLUM_mat_B_film[,1:12]) {
  fit <- aov(i~Treatment1,data = PHYLUM_mat_B_film)
  print(summary(fit))
  out <- HSD.test(fit,"Treatment1")
  print(out$groups)
}




## Use hierarchical clustering to order samples
norm_clu_B <- hclust(all_dis_B, "average")
norm_clu_B$height
Beta_labs_B <- norm_clu_B$labels
Beta_order_B <- norm_clu_B$order
Beta_labs_B <- Beta_labs_B[Beta_order_B]
PHYLUM_mat_B <- PHYLUM_mat_B[ ,Beta_labs_B]


## Soil phyla abundances
PHYLUM_mat_B_soil <- PHYLUM_mat_B[,onesoilsamples]
colSums(PHYLUM_mat_B_soil)
PHYLUM_mat_B_soil_mean <- sort(apply(PHYLUM_mat_B_soil,1,mean),decr=T)
PHYLUM_mat_B_soil <- PHYLUM_mat_B_soil[names(PHYLUM_mat_B_soil_mean),]
PHYLUM_mat_B_soil_se <- apply(PHYLUM_mat_B_soil,1,se)[names(PHYLUM_mat_B_soil_mean)]

length(PHYLUM_mat_B_soil_mean[PHYLUM_mat_B_soil_mean > 0])

## Root phyla abundances
PHYLUM_mat_B_root <- PHYLUM_mat_B[,rootsamples]
colSums(PHYLUM_mat_B_root)
PHYLUM_mat_B_root_mean <- apply(PHYLUM_mat_B_root,1,mean)[names(PHYLUM_mat_B_soil_mean)]
PHYLUM_mat_B_root_se <- apply(PHYLUM_mat_B_root,1,se)[names(PHYLUM_mat_B_soil_mean)]

length(PHYLUM_mat_B_root_mean[PHYLUM_mat_B_root_mean > 0])

## endophytes phyla abundances
PHYLUM_mat_B_endophytes <- PHYLUM_mat_B[,endsamples]
colSums(PHYLUM_mat_B_endophytes)
PHYLUM_mat_B_endophytes_mean <- sort(apply(PHYLUM_mat_B_endophytes,1,mean),decr=T)
PHYLUM_mat_B_endophytes <- PHYLUM_mat_B_endophytes[names(PHYLUM_mat_B_endophytes_mean),]
PHYLUM_mat_B_endophytes_se <- apply(PHYLUM_mat_B_endophytes,1,se)[names(PHYLUM_mat_B_endophytes_mean)]

length(PHYLUM_mat_B_endophytes_mean[PHYLUM_mat_B_endophytes_mean > 0])

## endophytes phyla abundances
PHYLUM_mat_B_endophytes <- PHYLUM_mat_B[,endsamples]
colSums(PHYLUM_mat_B_endophytes)
PHYLUM_mat_B_endophytes_mean <- sort(apply(PHYLUM_mat_B_endophytes,1,mean),decr=T)
PHYLUM_mat_B_endophytes <- PHYLUM_mat_B_endophytes[names(PHYLUM_mat_B_endophytes_mean),]
PHYLUM_mat_B_endophytes_se <- apply(PHYLUM_mat_B_endophytes,1,se)[names(PHYLUM_mat_B_endophytes_mean)]

length(PHYLUM_mat_B_endophytes_mean[PHYLUM_mat_B_endophytes_mean > 0])

## particle phyla abundances
PHYLUM_mat_B_particle <- PHYLUM_mat_B[,particlesamples]
colSums(PHYLUM_mat_B_particle)
PHYLUM_mat_B_particle_mean <- sort(apply(PHYLUM_mat_B_particle,1,mean),decr=T)
PHYLUM_mat_B_particle <- PHYLUM_mat_B_particle[names(PHYLUM_mat_B_particle_mean),]
PHYLUM_mat_B_particle_se <- apply(PHYLUM_mat_B_particle,1,se)[names(PHYLUM_mat_B_particle_mean)]

length(PHYLUM_mat_B_particle_mean[PHYLUM_mat_B_particle_mean > 0])


## Make matrix of phyla abundances by cropping system
PHYLUM_mat_B_mean_type <-as.matrix(cbind(`one_soil_CK`=apply(PHYLUM_mat_B[,one_soil_CK],1,mean),
                                         `one_root_CK`=apply(PHYLUM_mat_B[,one_root_CK],1,mean),
                                         `one_end_CK`=apply(PHYLUM_mat_B[,one_end_CK],1,mean),
                                         `two_soil_CK`=apply(PHYLUM_mat_B[,two_soil_CK],1,mean),
                                         `two_root_CK`=apply(PHYLUM_mat_B[,two_root_CK],1,mean),
                                         `two_end_CK`=apply(PHYLUM_mat_B[,two_end_CK],1,mean),
                                         `one_soil_PE`=apply(PHYLUM_mat_B[,one_soil_PE],1,mean),
                                         `one_root_PE`=apply(PHYLUM_mat_B[,one_root_PE],1,mean),
                                         `one_end_PE`=apply(PHYLUM_mat_B[,one_end_PE],1,mean),
                                         `one_film_PE`=apply(PHYLUM_mat_B[,one_film_PE],1,mean),
                                         `two_soil_PE`=apply(PHYLUM_mat_B[,two_soil_PE],1,mean),
                                         `two_root_PE`=apply(PHYLUM_mat_B[,two_root_PE],1,mean),
                                         `two_end_PE`=apply(PHYLUM_mat_B[,two_end_PE],1,mean),
                                         `two_film_PE`=apply(PHYLUM_mat_B[,two_film_PE],1,mean),
                                         `one_soil_DE`=apply(PHYLUM_mat_B[,one_soil_DE],1,mean),
                                         `one_root_DE`=apply(PHYLUM_mat_B[,one_root_DE],1,mean),
                                         `one_end_DE`=apply(PHYLUM_mat_B[,one_end_DE],1,mean),
                                         `one_film_DE`=apply(PHYLUM_mat_B[,one_film_DE],1,mean),
                                         `two_soil_DE`=apply(PHYLUM_mat_B[,two_soil_DE],1,mean),
                                         `two_root_DE`=apply(PHYLUM_mat_B[,two_root_DE],1,mean),
                                         `two_end_DE`=apply(PHYLUM_mat_B[,two_end_DE],1,mean),
                                         `two_film_DE`=apply(PHYLUM_mat_B[,two_film_DE],1,mean)))



new_B_cols<-c("one_soil_CK","one_soil_PE","one_soil_DE","one_root_CK","one_root_PE","one_root_DE","one_end_CK","one_end_PE","one_end_DE","one_film_PE","one_film_DE","two_soil_CK","two_soil_PE","two_soil_DE","two_root_CK","two_root_PE","two_root_DE","two_end_CK","two_end_PE","two_end_DE","two_film_PE","two_film_DE")

PHYLUM_mat_B_mean_type <- PHYLUM_mat_B_mean_type[, c(new_B_cols)]



### Defining bOTU colors by phylum (using the taxonomy file)
tax_B$cols <- tax_B$labels
table(tax_B$cols)

# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_B, 1, mean) < 1)
low_count_phyla_B <- rownames(PHYLUM_mat_B_mean_type)[sort(apply(PHYLUM_mat_B, 1, mean), decr=T) < 1]
# attribute grey color
for(i in low_count_phyla_B){
  tax_B[ rownames(tax_B)[tax_B$Phylum==paste(i) ], ]$cols <- "lightgrey"
}
table(tax_B$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
abundant_phyla_B <- rownames(PHYLUM_mat_B)[sort(apply(PHYLUM_mat_B, 1, mean), decr=T) > 1]
abundant_phyla_B
tax_B[ rownames(tax_B)[tax_B$labels=="Betaproteobacteria" ], ]$cols <- "lightcoral"
tax_B[ rownames(tax_B)[tax_B$labels=="Gammaproteobacteria" ], ]$cols <- "palegreen3"
tax_B[ rownames(tax_B)[tax_B$labels=="Alphaproteobacteria" ], ]$cols <- "palegreen1"
tax_B[ rownames(tax_B)[tax_B$labels=="Acidobacteria" ], ]$cols <- "palegreen4"
tax_B[ rownames(tax_B)[tax_B$labels=="Actinobacteria" ], ]$cols <- "dodgerblue"
tax_B[ rownames(tax_B)[tax_B$labels=="Bacteroidetes" ], ]$cols <- "sandybrown"
tax_B[ rownames(tax_B)[tax_B$labels=="Chloroflexi" ], ]$cols <- "plum1"
tax_B[ rownames(tax_B)[tax_B$labels=="Deltaproteobacteria" ], ]$cols <- "palegreen2"
tax_B[ rownames(tax_B)[tax_B$labels=="Firmicutes" ], ]$cols <- "steelblue1"
tax_B[ rownames(tax_B)[tax_B$labels=="Gemmatimonadetes" ], ]$cols <- "steelblue4"
tax_B[ rownames(tax_B)[tax_B$labels=="Saccharibacteria" ], ]$cols <- "palevioletred3"

## collaps OTU colors to prepare Phylum level colors
label_cols_B <- tax_B[, c("labels", "cols") ]
library(plyr)
PHYLA_label_cols_B <- ddply(label_cols_B, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_B) <- PHYLA_label_cols_B[,1]
PHYLA_label_cols_B <- PHYLA_label_cols_B[c(abundant_phyla_B, low_count_phyla_B),]
PHYLA_label_cols_B

## Legend for Phylum colors
PHYLA_label_cols_B_legend <- PHYLA_label_cols_B[1:12,]
PHYLA_label_cols_B_legend[12,1] <- "other"
rownames(PHYLA_label_cols_B_legend)[12] <- "other"
PHYLA_label_cols_B_legend


##### Plot Supplementary Figure S3
par(oma=c(0,0,0,0), mar=c(6,4,1,5), xpd=NA)
phylum_bar_B <- barplot(as.matrix(PHYLUM_mat_B_mean_type), col=PHYLA_label_cols_B[rownames(PHYLUM_mat_B_mean_type),]$cols,
                        ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_B, labels=colnames(PHYLUM_mat_B_mean_type), col.axis="black", las=2, cex.axis=1)
title(ylab="Relative abundance (%)")
title(main="Bacteria Community")
legend(26, 30, bty="n", cex=0.7, x.intersp= 0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_B_legend$labels), 
       fill=rev(PHYLA_label_cols_B_legend$cols), 
       border=rev(PHYLA_label_cols_B_legend$cols) )



##### Fungi

## Express B OTU counts as relative abunance percent
otu_F_RA <- t(t(otu_F)/colSums(otu_F)) * 100
colSums(otu_F_RA)
nrow(otu_F_RA)

## Get names of fungi phyla present (use 'labels' as this specifies class within Proteobacteria)
PHYLAnames_F <- names(sort(table(tax_F[,"labels"]), decr=T))
length(PHYLAnames_F)
sort(table(tax_F[,"labels"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(otu_F_RA)
for (i in PHYLAnames_F){
  x <- array(colSums(otu_F_RA[rownames(tax_F)[which(tax_F$labels == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(PHYLAnames_F)
colnames(y) <- paste(colnames(otu_F_RA))
PHYLUM_mat_F <- y
PHYLUM_mat_F[,1:5]
colSums(PHYLUM_mat_F)
PHYLUM_mat_F_mean <- sort(apply(PHYLUM_mat_F,1,mean),decr=T)
PHYLUM_mat_F <- PHYLUM_mat_F[names(PHYLUM_mat_F_mean),]

## Use hierarchical clustering to order samples
norm_clu_F <- hclust(all_dis_F, "average")
norm_clu_F$height
Beta_labs_F <- norm_clu_F$labels
Beta_order_F <- norm_clu_F$order
Beta_labs_F <- Beta_labs_F[Beta_order_F]
PHYLUM_mat_F <- PHYLUM_mat_F[ ,Beta_labs_F]



## Express B OTU counts as relative abunance percent
design_F <- unite(design_F, Treatment2, c(Site, Season))

one_soil <- rownames(design_F[design_F$Treatment2=="soil_one",])
one_root <- rownames(design_F[design_F$Treatment2=="root_one",])
one_end <- rownames(design_F[design_F$Treatment2=="end_one",])
two_soil <- rownames(design_F[design_F$Treatment2=="soil_two",])
two_root <- rownames(design_F[design_F$Treatment2=="root_two",])
two_end <- rownames(design_F[design_F$Treatment2=="end_two",])

one_film <- rownames(design_F[design_F$Treatment2=="film_one",])
two_film <- rownames(design_F[design_F$Treatment2=="film_two",])



## Soil phyla abundances
PHYLUM_mat_F_end <- PHYLUM_mat_F[,two_end]
colSums(PHYLUM_mat_F_end)
PHYLUM_mat_F_end_mean <- sort(apply(PHYLUM_mat_F_end,1,mean),decr=T)
PHYLUM_mat_F_end <- PHYLUM_mat_F_end[names(PHYLUM_mat_F_end_mean),]
PHYLUM_mat_F_end_se <- apply(PHYLUM_mat_F_end,1,se)[names(PHYLUM_mat_F_end_mean)]

length(PHYLUM_mat_F_end_mean[PHYLUM_mat_F_end_mean > 0])

design_F_one_end <- design_F[two_end,]


PHYLUM_mat_F_end <- t(PHYLUM_mat_F_end)
PHYLUM_mat_F_end <- as.data.frame(PHYLUM_mat_F_end)
PHYLUM_mat_F_end$Treatment1 <- design_F_one_end$Treatment1

for(i in PHYLUM_mat_F_end[,1:11]) {
  fit <- aov(i~Treatment1,data = PHYLUM_mat_F_end)
  print(summary(fit))
  out <- HSD.test(fit,"Treatment1")
  print(out$groups)
}





## Soil phyla abundances
PHYLUM_mat_F_soil <- PHYLUM_mat_F[,soilsamples]
colSums(PHYLUM_mat_F_soil)
PHYLUM_mat_F_soil_mean <- sort(apply(PHYLUM_mat_F_soil,1,mean),decr=T)
PHYLUM_mat_F_soil <- PHYLUM_mat_F_soil[names(PHYLUM_mat_F_soil_mean),]
PHYLUM_mat_F_soil_se <- apply(PHYLUM_mat_F_soil,1,se)[names(PHYLUM_mat_F_soil_mean)]

length(PHYLUM_mat_F_soil_mean[PHYLUM_mat_F_soil_mean > 0])

## Root phyla abundances
PHYLUM_mat_F_root <- PHYLUM_mat_F[,rootsamples]
colSums(PHYLUM_mat_F_root)
PHYLUM_mat_F_root_mean <- apply(PHYLUM_mat_F_root,1,mean)[names(PHYLUM_mat_F_soil_mean)]
PHYLUM_mat_F_root_se <- apply(PHYLUM_mat_F_root,1,se)[names(PHYLUM_mat_F_soil_mean)]

length(PHYLUM_mat_F_root_mean[PHYLUM_mat_F_root_mean > 0])

## endophytes phyla abundances
PHYLUM_mat_F_endophytes <- PHYLUM_mat_F[,endophytessamples]
colSums(PHYLUM_mat_F_endophytes)
PHYLUM_mat_F_endophytes_mean <- sort(apply(PHYLUM_mat_F_endophytes,1,mean),decr=T)
PHYLUM_mat_F_endophytes <- PHYLUM_mat_F_endophytes[names(PHYLUM_mat_F_endophytes_mean),]
PHYLUM_mat_F_endophytes_se <- apply(PHYLUM_mat_F_endophytes,1,se)[names(PHYLUM_mat_F_endophytes_mean)]

length(PHYLUM_mat_F_endophytes_mean[PHYLUM_mat_F_endophytes_mean > 0])

## endophytes phyla abundances
PHYLUM_mat_F_endophytes <- PHYLUM_mat_F[,endophytessamples]
colSums(PHYLUM_mat_F_endophytes)
PHYLUM_mat_F_endophytes_mean <- sort(apply(PHYLUM_mat_F_endophytes,1,mean),decr=T)
PHYLUM_mat_F_endophytes <- PHYLUM_mat_F_endophytes[names(PHYLUM_mat_F_endophytes_mean),]
PHYLUM_mat_F_endophytes_se <- apply(PHYLUM_mat_F_endophytes,1,se)[names(PHYLUM_mat_F_endophytes_mean)]

length(PHYLUM_mat_F_endophytes_mean[PHYLUM_mat_F_endophytes_mean > 0])

## particle phyla abundances
PHYLUM_mat_F_particle <- PHYLUM_mat_F[,particlesamples]
colSums(PHYLUM_mat_F_particle)
PHYLUM_mat_F_particle_mean <- sort(apply(PHYLUM_mat_F_particle,1,mean),decr=T)
PHYLUM_mat_F_particle <- PHYLUM_mat_F_particle[names(PHYLUM_mat_F_particle_mean),]
PHYLUM_mat_F_particle_se <- apply(PHYLUM_mat_F_particle,1,se)[names(PHYLUM_mat_F_particle_mean)]

length(PHYLUM_mat_F_particle_mean[PHYLUM_mat_F_particle_mean > 0])


## Make matrix of phyla abundances by cropping system
PHYLUM_mat_F_mean_type <-as.matrix(cbind(`one_soil_CK`=apply(PHYLUM_mat_F[,one_soil_CK],1,mean),
                                         `one_root_CK`=apply(PHYLUM_mat_F[,one_root_CK],1,mean),
                                         `one_end_CK`=apply(PHYLUM_mat_F[,one_end_CK],1,mean),
                                         `two_soil_CK`=apply(PHYLUM_mat_F[,two_soil_CK],1,mean),
                                         `two_root_CK`=apply(PHYLUM_mat_F[,two_root_CK],1,mean),
                                         `two_end_CK`=apply(PHYLUM_mat_F[,two_end_CK],1,mean),
                                         `one_soil_PE`=apply(PHYLUM_mat_F[,one_soil_PE],1,mean),
                                         `one_root_PE`=apply(PHYLUM_mat_F[,one_root_PE],1,mean),
                                         `one_end_PE`=apply(PHYLUM_mat_F[,one_end_PE],1,mean),
                                         `one_film_PE`=apply(PHYLUM_mat_F[,one_film_PE],1,mean),
                                         `two_soil_PE`=apply(PHYLUM_mat_F[,two_soil_PE],1,mean),
                                         `two_root_PE`=apply(PHYLUM_mat_F[,two_root_PE],1,mean),
                                         `two_end_PE`=apply(PHYLUM_mat_F[,two_end_PE],1,mean),
                                         `two_film_PE`=apply(PHYLUM_mat_F[,two_film_PE],1,mean),
                                         `one_soil_DE`=apply(PHYLUM_mat_F[,one_soil_DE],1,mean),
                                         `one_root_DE`=apply(PHYLUM_mat_F[,one_root_DE],1,mean),
                                         `one_end_DE`=apply(PHYLUM_mat_F[,one_end_DE],1,mean),
                                         `one_film_DE`=apply(PHYLUM_mat_F[,one_film_DE],1,mean),
                                         `two_soil_DE`=apply(PHYLUM_mat_F[,two_soil_DE],1,mean),
                                         `two_root_DE`=apply(PHYLUM_mat_F[,two_root_DE],1,mean),
                                         `two_end_DE`=apply(PHYLUM_mat_F[,two_end_DE],1,mean),
                                         `two_film_DE`=apply(PHYLUM_mat_F[,two_film_DE],1,mean)))



new_F_cols<-c("one_soil_CK","one_soil_PE","one_soil_DE","one_root_CK","one_root_PE","one_root_DE","one_end_CK","one_end_PE","one_end_DE","one_film_PE","one_film_DE","two_soil_CK","two_soil_PE","two_soil_DE","two_root_CK","two_root_PE","two_root_DE","two_end_CK","two_end_PE","two_end_DE","two_film_PE","two_film_DE")

PHYLUM_mat_F_mean_type <- PHYLUM_mat_F_mean_type[, c(new_F_cols)]


# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_F, 1, mean) < 1)
low_count_phyla_F <- rownames(PHYLUM_mat_F)[apply(PHYLUM_mat_F, 1, mean) < 1]
# attribute grey color
for(i in low_count_phyla_F){
  tax_F[ rownames(tax_F)[tax_F$Phylum==paste(i) ], ]$cols <- "lightgrey"
}
table(tax_F$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
phyla_F <- names(sort(apply(PHYLUM_mat_F, 1, mean), decr=T) )
phyla_F
tax_F[ rownames(tax_F)[tax_F$Phylum=="Ascomycota" ], ]$cols <- "dodgerblue2"
tax_F[ rownames(tax_F)[tax_F$Class=="Sordariomycetes" ], ]$cols <- "deepskyblue"
tax_F[ rownames(tax_F)[tax_F$Class=="Dothideomycetes" ], ]$cols <- "lightskyblue"
tax_F[ rownames(tax_F)[tax_F$Class=="Eurotiomycetes" ], ]$cols <- "cornflowerblue"
tax_F[ rownames(tax_F)[tax_F$Class=="Orbiliomycetes" ], ]$cols <- "skyblue1"
tax_F[ rownames(tax_F)[tax_F$Phylum=="Basidiomycota" ], ]$cols <- "firebrick1"
tax_F[ rownames(tax_F)[tax_F$Class=="Agaricomycetes" ], ]$cols <- "tomato"
tax_F[ rownames(tax_F)[tax_F$Class=="Tremellomycetes" ], ]$cols <- "indianred1"
tax_F[ rownames(tax_F)[tax_F$Phylum=="Chytridiomycota" ], ]$cols <- "goldenrod2"
tax_F[ rownames(tax_F)[tax_F$Phylum=="Mortierellomycota" ], ]$cols <- "seagreen3"
tax_F[ rownames(tax_F)[tax_F$Phylum=="Rozellomycota" ], ]$cols <- "mediumorchid1"


## collaps OTU colors to prepare Phylum colors
PHYLA_label_cols_F <- tax_F[,c("labels", "cols")]
library(plyr)
PHYLA_label_cols_F <- ddply(PHYLA_label_cols_F, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_F) <- PHYLA_label_cols_F[,1]
PHYLA_label_cols_F <- PHYLA_label_cols_F[phyla_F, ]
PHYLA_label_cols_F

## Legend for Phylum colors
PHYLA_label_cols_F_legend <- PHYLA_label_cols_F[1:12,]
PHYLA_label_cols_F_legend[12, 1] <- "other"
rownames(PHYLA_label_cols_F_legend)[12] <- "other"
PHYLA_label_cols_F_legend



##### Plot Supplementary Figure S3
par(oma=c(0,0,0,0), mar=c(6,4,1,5), xpd=NA)
phylum_bar_F <- barplot(as.matrix(PHYLUM_mat_F_mean_type), col=PHYLA_label_cols_F[rownames(PHYLUM_mat_F_mean_type),]$cols,ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_F, labels=colnames(PHYLUM_mat_F_mean_type), col.axis="black", las=2, cex.axis=0.9)
title(ylab="Relative abundance (%)")
title(main="FUNGI Community")
legend(26, 100, bty="n", cex=0.7, x.intersp=0.1, y.intersp=1, 
       legend=rev(PHYLA_label_cols_F_legend$labels),
       fill=rev(PHYLA_label_cols_F_legend$cols),
       border=rev(PHYLA_label_cols_F_legend$cols) )



##### Protist

## Express B OTU counts as relative abunance percent
otu_P_RA <- t(t(otu_P)/colSums(otu_P)) * 100
colSums(otu_P_RA)
nrow(otu_P_RA)

## Get names of fungi phyla present (use 'labels' as this specifies class within Proteobacteria)
PHYLAnames_P <- names(sort(table(tax_P[,"labels"]), decr=T))
length(PHYLAnames_P)
sort(table(tax_P[,"labels"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(otu_P_RA)
for (i in PHYLAnames_P){
  x <- array(colSums(otu_P_RA[rownames(tax_P)[which(tax_P$labels == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(PHYLAnames_P)
colnames(y) <- paste(colnames(otu_P_RA))
PHYLUM_mat_P <- y
PHYLUM_mat_P[,1:5]
colSums(PHYLUM_mat_P)
PHYLUM_mat_P_mean <- sort(apply(PHYLUM_mat_P,1,mean),decr=T)
PHYLUM_mat_P <- PHYLUM_mat_P[names(PHYLUM_mat_P_mean),]

## Use hierarchical clustering to order samples
norm_clu_P <- hclust(all_dis_P, "average")
norm_clu_P$height
Beta_labs_P <- norm_clu_P$labels
Beta_order_P <- norm_clu_P$order
Beta_labs_P <- Beta_labs_P[Beta_order_P]
PHYLUM_mat_P <- PHYLUM_mat_P[ ,Beta_labs_P]






## Express B OTU counts as relative abunance percent
design_P <- unite(design_P, Treatment2, c(Site, Season))

one_soil <- rownames(design_P[design_P$Treatment2=="soil_one",])
one_root <- rownames(design_P[design_P$Treatment2=="root_one",])
one_end <- rownames(design_P[design_P$Treatment2=="end_one",])
two_soil <- rownames(design_P[design_P$Treatment2=="soil_two",])
two_root <- rownames(design_P[design_P$Treatment2=="root_two",])
two_end <- rownames(design_P[design_P$Treatment2=="end_two",])

one_film <- rownames(design_P[design_P$Treatment2=="film_one",])
two_film <- rownames(design_P[design_P$Treatment2=="film_two",])



## Soil phyla abundances
PHYLUM_mat_P_film <- PHYLUM_mat_P[,one_root]
colSums(PHYLUM_mat_P_film)
PHYLUM_mat_P_film_mean <- sort(apply(PHYLUM_mat_P_film,1,mean),decr=T)
PHYLUM_mat_P_film <- PHYLUM_mat_P_film[names(PHYLUM_mat_P_film_mean),]
PHYLUM_mat_P_film_se <- apply(PHYLUM_mat_P_film,1,se)[names(PHYLUM_mat_P_film_mean)]

length(PHYLUM_mat_P_film_mean[PHYLUM_mat_P_film_mean > 0])

design_P_one_film <- design_P[one_root,]


PHYLUM_mat_P_film <- t(PHYLUM_mat_P_film)
PHYLUM_mat_P_film <- as.data.frame(PHYLUM_mat_P_film)
PHYLUM_mat_P_film$Treatment1 <- design_P_one_film$Treatment1

for(i in PHYLUM_mat_P_film[,1:7]) {
  fit <- aov(i~Treatment1,data = PHYLUM_mat_P_film)
  print(summary(fit))
  out <- HSD.test(fit,"Treatment1")
  print(out$groups)
}



## Soil phyla abundances
PHYLUM_mat_P_soil <- PHYLUM_mat_P[,soilsamples]
colSums(PHYLUM_mat_P_soil)
PHYLUM_mat_P_soil_mean <- sort(apply(PHYLUM_mat_P_soil,1,mean),decr=T)
PHYLUM_mat_P_soil <- PHYLUM_mat_P_soil[names(PHYLUM_mat_P_soil_mean),]
PHYLUM_mat_P_soil_se <- apply(PHYLUM_mat_P_soil,1,se)[names(PHYLUM_mat_P_soil_mean)]

length(PHYLUM_mat_P_soil_mean[PHYLUM_mat_P_soil_mean > 0])

## Root phyla abundances
PHYLUM_mat_P_root <- PHYLUM_mat_P[,rootsamples]
colSums(PHYLUM_mat_P_root)
PHYLUM_mat_P_root_mean <- apply(PHYLUM_mat_P_root,1,mean)[names(PHYLUM_mat_P_soil_mean)]
PHYLUM_mat_P_root_se <- apply(PHYLUM_mat_P_root,1,se)[names(PHYLUM_mat_P_soil_mean)]

length(PHYLUM_mat_P_root_mean[PHYLUM_mat_P_root_mean > 0])

## endophytes phyla abundances
PHYLUM_mat_P_endophytes <- PHYLUM_mat_P[,endophytessamples]
colSums(PHYLUM_mat_P_endophytes)
PHYLUM_mat_P_endophytes_mean <- sort(apply(PHYLUM_mat_P_endophytes,1,mean),decr=T)
PHYLUM_mat_P_endophytes <- PHYLUM_mat_P_endophytes[names(PHYLUM_mat_P_endophytes_mean),]
PHYLUM_mat_P_endophytes_se <- apply(PHYLUM_mat_P_endophytes,1,se)[names(PHYLUM_mat_P_endophytes_mean)]

length(PHYLUM_mat_P_endophytes_mean[PHYLUM_mat_P_endophytes_mean > 0])

## endophytes phyla abundances
PHYLUM_mat_P_endophytes <- PHYLUM_mat_P[,endophytessamples]
colSums(PHYLUM_mat_P_endophytes)
PHYLUM_mat_P_endophytes_mean <- sort(apply(PHYLUM_mat_P_endophytes,1,mean),decr=T)
PHYLUM_mat_P_endophytes <- PHYLUM_mat_P_endophytes[names(PHYLUM_mat_P_endophytes_mean),]
PHYLUM_mat_P_endophytes_se <- apply(PHYLUM_mat_P_endophytes,1,se)[names(PHYLUM_mat_P_endophytes_mean)]

length(PHYLUM_mat_P_endophytes_mean[PHYLUM_mat_P_endophytes_mean > 0])

## particle phyla abundances
PHYLUM_mat_P_particle <- PHYLUM_mat_P[,particlesamples]
colSums(PHYLUM_mat_P_particle)
PHYLUM_mat_P_particle_mean <- sort(apply(PHYLUM_mat_P_particle,1,mean),decr=T)
PHYLUM_mat_P_particle <- PHYLUM_mat_P_particle[names(PHYLUM_mat_P_particle_mean),]
PHYLUM_mat_P_particle_se <- apply(PHYLUM_mat_P_particle,1,se)[names(PHYLUM_mat_P_particle_mean)]

length(PHYLUM_mat_P_particle_mean[PHYLUM_mat_P_particle_mean > 0])


## Make matrix of phyla abundances by cropping system
PHYLUM_mat_P_mean_type <-as.matrix(cbind(`one_soil_CK`=apply(PHYLUM_mat_P[,one_soil_CK],1,mean),
                                         `one_root_CK`=apply(PHYLUM_mat_P[,one_root_CK],1,mean),
                                         `two_soil_CK`=apply(PHYLUM_mat_P[,two_soil_CK],1,mean),
                                         `two_root_CK`=apply(PHYLUM_mat_P[,two_root_CK],1,mean),
                                         `two_end_CK`=apply(PHYLUM_mat_P[,two_end_CK],1,mean),
                                         `one_soil_PE`=apply(PHYLUM_mat_P[,one_soil_PE],1,mean),
                                         `one_root_PE`=apply(PHYLUM_mat_P[,one_root_PE],1,mean),
                                         `one_film_PE`=apply(PHYLUM_mat_P[,one_film_PE],1,mean),
                                         `two_soil_PE`=apply(PHYLUM_mat_P[,two_soil_PE],1,mean),
                                         `two_root_PE`=apply(PHYLUM_mat_P[,two_root_PE],1,mean),
                                         `two_end_PE`=apply(PHYLUM_mat_P[,two_end_PE],1,mean),
                                         `two_film_PE`=apply(PHYLUM_mat_P[,two_film_PE],1,mean),
                                         `one_soil_DE`=apply(PHYLUM_mat_P[,one_soil_DE],1,mean),
                                         `one_root_DE`=apply(PHYLUM_mat_P[,one_root_DE],1,mean),
                                         `one_film_DE`=apply(PHYLUM_mat_P[,one_film_DE],1,mean),
                                         `two_soil_DE`=apply(PHYLUM_mat_P[,two_soil_DE],1,mean),
                                         `two_root_DE`=apply(PHYLUM_mat_P[,two_root_DE],1,mean),
                                         `two_end_DE`=apply(PHYLUM_mat_P[,two_end_DE],1,mean),
                                         `two_film_DE`=apply(PHYLUM_mat_P[,two_film_DE],1,mean)))


new_P_cols<-c("one_soil_CK","one_soil_PE","one_soil_DE","one_root_CK","one_root_PE","one_root_DE","one_film_PE","one_film_DE","two_soil_CK","two_soil_PE","two_soil_DE","two_root_CK","two_root_PE","two_root_DE","two_end_CK","two_end_PE","two_end_DE","two_film_PE","two_film_DE")

PHYLUM_mat_P_mean_type <- PHYLUM_mat_P_mean_type[, c(new_P_cols)]



# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_P, 1, mean) < 1)
low_count_phyla_P <- rownames(PHYLUM_mat_P)[apply(PHYLUM_mat_P, 1, mean) < 1]
# attribute grey color
for(i in low_count_phyla_P){
  tax_P[ rownames(tax_P)[tax_P$labels==paste(i) ], ]$cols <- "lightgrey"
}
table(tax_P$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
phyla_P <- names(sort(apply(PHYLUM_mat_P, 1, mean), decr=T) )
phyla_P
tax_P[ rownames(tax_P)[tax_P$labels=="Opisthokonta" ], ]$cols <- "#78fee0"
tax_P[ rownames(tax_P)[tax_P$labels=="Metazoa" ], ]$cols <- "#F38181"
tax_P[ rownames(tax_P)[tax_P$labels=="Alveolata" ], ]$cols <- "#ffad60"
tax_P[ rownames(tax_P)[tax_P$labels=="Amoebozoa" ], ]$cols <- "#ceefe4"
tax_P[ rownames(tax_P)[tax_P$labels=="Archaeplastida" ], ]$cols <- "#9dd3a8"
tax_P[ rownames(tax_P)[tax_P$labels=="Rhizaria" ], ]$cols <- "#96ceb4"
tax_P[ rownames(tax_P)[tax_P$labels=="Stramenopiles" ], ]$cols <- "#feda77"


## collaps OTU colors to prepare Phylum colors
PHYLA_label_cols_P <- tax_P[,c("labels", "cols")]
library(plyr)
PHYLA_label_cols_P <- ddply(PHYLA_label_cols_P, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_P) <- PHYLA_label_cols_P[,1]
PHYLA_label_cols_P <- PHYLA_label_cols_P[phyla_P, ]
PHYLA_label_cols_P

## Legend for Phylum colors
PHYLA_label_cols_P_legend <- PHYLA_label_cols_P[1:8,]
PHYLA_label_cols_P_legend[8, 1] <- "other"
rownames(PHYLA_label_cols_P_legend)[8] <- "other"
PHYLA_label_cols_P_legend



##### Plot Supplementary Figure S3
par(oma=c(0,0,0,0), mar=c(6,4,1,5), xpd=NA)
phylum_bar_P <- barplot(as.matrix(PHYLUM_mat_P_mean_type), col=PHYLA_label_cols_P[rownames(PHYLUM_mat_P_mean_type),]$cols,ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_P, labels=colnames(PHYLUM_mat_P_mean_type), col.axis="black", las=2, cex.axis=0.9)
title(ylab="Relative abundance (%)")
title(main="Protist Community")
legend(30, 30, bty="n", cex=0.7, x.intersp=0.1, y.intersp=1, 
       legend=rev(PHYLA_label_cols_P_legend$labels),
       fill=rev(PHYLA_label_cols_P_legend$cols),
       border=rev(PHYLA_label_cols_P_legend$cols))




## Kingdom network type ##
## PE
otu_B_PE <- otu_B[, c(soilPEsamples,rootPEsamples,PEsamples) ]
keep_otu_B_PE <- which(rowSums(otu_B_PE >= 8) >= 12)

otu_B_PE <- otu_B_PE[keep_otu_B_PE,]

dim(otu_B_PE)

tax_PE_B <- tax_B[rownames(otu_B_PE),]
design_B_PE <- droplevels(design_B[c(soilPEsamples,rootPEsamples,PEsamples),])

edgeR_B_PE <- DGEList(counts=otu_B_PE, 
                      group=design_B_PE$Type,
                      genes=tax_PE_B)

edgeR_B_PE <- calcNormFactors(edgeR_B_PE)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_PE_B <- cpm(edgeR_B_PE, normalized.lib.sizes=T, log=F)



otu_B_DE <- otu_B[, c(soilDEsamples,rootDEsamples,DEsamples)]
keep_otu_B_DE <- which(rowSums(otu_B_DE >= 8) >= 12)

otu_B_DE <- otu_B_DE[keep_otu_B_DE,]

dim(otu_B_DE)

tax_DE_B <- tax_B[rownames(otu_B_DE),]
design_B_DE <- droplevels(design_B[c(soilDEsamples,rootDEsamples,DEsamples),])

edgeR_B_DE <- DGEList(counts=otu_B_DE, 
                      group=design_B_DE$Type,
                      genes=tax_DE_B)

edgeR_B_DE <- calcNormFactors(edgeR_B_DE)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_DE_B <- cpm(edgeR_B_DE, normalized.lib.sizes=T, log=F)



otu_B_CK <- otu_B[, c(soilCKsamples,rootCKsamples) ]
keep_otu_B_CK <- which(rowSums(otu_B_CK >= 8) >= 12)

otu_B_CK <- otu_B_CK[keep_otu_B_CK,]

dim(otu_B_CK)

tax_CK_B <- tax_B[rownames(otu_B_CK),]
design_B_CK <- droplevels(design_B[c(soilCKsamples,rootCKsamples),])

edgeR_B_CK <- DGEList(counts=otu_B_CK, 
                      group=design_B_CK$Type,
                      genes=tax_CK_B)

edgeR_B_CK <- calcNormFactors(edgeR_B_CK)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_CK_B <- cpm(edgeR_B_CK, normalized.lib.sizes=T, log=F)




otu_F_PE <- otu_F[, c(soilPEsamples,rootPEsamples,PEsamples) ]
keep_otu_F_PE <- which(rowSums(otu_F_PE >= 6) >= 8)

otu_F_PE <- otu_F_PE[keep_otu_F_PE,]

dim(otu_F_PE)

tax_PE_F <- tax_F[rownames(otu_F_PE),]
design_F_PE <- droplevels(design_F[c(soilPEsamples,rootPEsamples,PEsamples),])

edgeR_F_PE <- DGEList(counts=otu_F_PE, 
                      group=design_F_PE$Type,
                      genes=tax_PE_F)

edgeR_F_PE <- calcNormFactors(edgeR_F_PE)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_PE_F <- cpm(edgeR_F_PE, normalized.lib.sizes=T, log=F)



otu_F_DE <- otu_F[, c(soilDEsamples,rootDEsamples,DEsamples) ]
keep_otu_F_DE <- which(rowSums(otu_F_DE >= 6) >= 8)

otu_F_DE <- otu_F_DE[keep_otu_F_DE,]

dim(otu_F_DE)

tax_DE_F <- tax_F[rownames(otu_F_DE),]
design_F_DE <- droplevels(design_F[c(soilDEsamples,rootDEsamples,DEsamples),])

edgeR_F_DE <- DGEList(counts=otu_F_DE, 
                      group=design_F_DE$Type,
                      genes=tax_DE_F)

edgeR_F_DE <- calcNormFactors(edgeR_F_DE)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_DE_F <- cpm(edgeR_F_DE, normalized.lib.sizes=T, log=F)



otu_F_CK <- otu_F[, c(soilCKsamples,rootCKsamples) ]
keep_otu_F_CK <- which(rowSums(otu_F_CK >= 6) >= 8)

otu_F_CK <- otu_F_CK[keep_otu_F_CK,]

dim(otu_F_CK)

tax_CK_F <- tax_F[rownames(otu_F_CK),]
design_F_CK <- droplevels(design_F[c(soilCKsamples,rootCKsamples),])

edgeR_F_CK <- DGEList(counts=otu_F_CK, 
                      group=design_F_CK$Type,
                      genes=tax_CK_F)

edgeR_F_CK <- calcNormFactors(edgeR_F_CK)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_CK_F <- cpm(edgeR_F_CK, normalized.lib.sizes=T, log=F)



otu_P_PE <- otu_P[, c(psoilPEsamples,prootPEsamples,pPEsamples) ]
keep_otu_P_PE <- which(rowSums(otu_P_PE >= 8) >= 12)

otu_P_PE <- otu_P_PE[keep_otu_P_PE,]

dim(otu_P_PE)

tax_PE_P <- tax_P[rownames(otu_P_PE),]
design_P_PE <- droplevels(design_P[c(psoilPEsamples,prootPEsamples,pPEsamples),])

edgeR_P_PE <- DGEList(counts=otu_P_PE, 
                      group=design_P_PE$Type,
                      genes=tax_PE_P)

edgeR_P_PE <- calcNormFactors(edgeR_P_PE)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_PE_P <- cpm(edgeR_P_PE, normalized.lib.sizes=T, log=F)



otu_P_DE <- otu_P[, c(psoilDEsamples,prootDEsamples,pDEsamples) ]
keep_otu_P_DE <- which(rowSums(otu_P_DE >= 8) >= 12)

otu_P_DE <- otu_P_DE[keep_otu_P_DE,]

dim(otu_P_DE)

tax_DE_P <- tax_P[rownames(otu_P_DE),]
design_P_DE <- droplevels(design_P[c(psoilDEsamples,prootDEsamples,pDEsamples),])

edgeR_P_DE <- DGEList(counts=otu_P_DE, 
                      group=design_P_DE$Type,
                      genes=tax_DE_P)

edgeR_P_DE <- calcNormFactors(edgeR_P_DE)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_DE_P <- cpm(edgeR_P_DE, normalized.lib.sizes=T, log=F)



otu_P_CK <- otu_P[, c(psoilCKsamples,prootCKsamples) ]
keep_otu_P_CK <- which(rowSums(otu_P_CK >= 8) >= 12)

otu_P_CK <- otu_P_CK[keep_otu_P_CK,]

dim(otu_P_CK)

tax_CK_P <- tax_P[rownames(otu_P_CK),]
design_P_CK <- droplevels(design_P[c(psoilCKsamples,prootCKsamples),])

edgeR_P_CK <- DGEList(counts=otu_P_CK, 
                      group=design_P_CK$Type,
                      genes=tax_CK_P)

edgeR_P_CK <- calcNormFactors(edgeR_P_CK)

## Get TMM normalized counts expressed as relative abundance counts per million 
otu_norm_CK_P <- cpm(edgeR_P_CK, normalized.lib.sizes=T, log=F)







##### Soil meta co-occurrence network creation and analysis and defining keystone OTUs #####

## Combine OTU counts of both kingdoms together
otu_norm_PE_combine <- rbind(otu_norm_PE_B, otu_norm_PE_F, otu_norm_PE_P)

tax_PE_B <- tax_PE_B[,-c(5:9)]
tax_PE_F <- tax_PE_F[,-c(5:9)]
tax_PE_P <- tax_PE_P[,-c(5:10)]

otu_cp_PE_tax <- rbind(tax_PE_B, tax_PE_F, tax_PE_P)

## Perform Spearman correlation of all OTU pairs
all_PE_cor <- rcorr(t(otu_norm_PE_combine), type=c("spearman"))

## Create data frame of co-occurring OTUs
all_cor_PE_df <- CorrDF(all_PE_cor$r, all_PE_cor$P)
all_cor_PE_df$padj <- p.adjust(all_cor_PE_df$p, method="none")

all_cor_PE_df_padj <- all_cor_PE_df[which(all_cor_PE_df$cor > 0.8 | all_cor_PE_df$cor < -0.8),]
all_cor_PE_df_padj <- all_cor_PE_df_padj[which(all_cor_PE_df_padj$padj < 0.001),]

## Make node attribute table
# indic_PE_combine <- c(indic_B_PE, indic_F_PE)

# PE_r_values_combine <- rbind(B_PE_r_values, F_PE_r_values)

nodeattrib_PE_combine <- data.frame(node=union(all_cor_PE_df_padj$from,all_cor_PE_df_padj$to))
nodeattrib_PE_combine$Kingdom <- 0

for (i in as.character(nodeattrib_PE_combine$node))
{
  if (i %in% rownames(otu_cp_PE_tax) == TRUE)
  {nodeattrib_PE_combine[nodeattrib_PE_combine$node==i,"Kingdom"] <- paste(otu_cp_PE_tax[i,1:1])}
  else
  { nodeattrib_PE_combine[nodeattrib_PE_combine$node==i,"Kingdom"]<- "NA"}
}

rownames(nodeattrib_PE_combine) <- as.character(nodeattrib_PE_combine$node)

all_PE_net <- graph_from_data_frame(all_cor_PE_df_padj,direct=F,vertices=nodeattrib_PE_combine)



## Number of nodes
length(V(all_PE_net))

## Number of bacteria and fungi nodes
length(grep("bASV*",names(V(all_PE_net))))
length(grep("fASV*",names(V(all_PE_net))))
length(grep("pASV*",names(V(all_PE_net))))

## Number of edges in network
length(E(all_PE_net))

## Connections 
bb_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_soil <- bb_occur_soil[bb_occur_soil$cor>0,]
nrow(bb_occur_soil)

ff_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("fASV*",from) & grepl("fASV*",to)),])
ff_occur_soil <- ff_occur_soil[ff_occur_soil$cor>0,]
nrow(ff_occur_soil)

pp_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("pASV*",from) & grepl("pASV*",to)),])
pp_occur_soil <- pp_occur_soil[pp_occur_soil$cor>0,]
nrow(pp_occur_soil)

bf_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("fASV*",from) & grepl("bASV*",to)),])
bf_occur_soil <- bf_occur_soil[bf_occur_soil$cor>0,]
nrow(bf_occur_soil)

bp_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("pASV*",from) & grepl("bASV*",to)),])
bp_occur_soil <- bp_occur_soil[bp_occur_soil$cor>0,]
nrow(bp_occur_soil)

pf_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("pASV*",from) & grepl("fASV*",to)),])
pf_occur_soil <- pf_occur_soil[pf_occur_soil$cor>0,]
nrow(pf_occur_soil)

bb_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_soil <- bb_occur_soil[bb_occur_soil$cor<0,]
nrow(bb_occur_soil)

ff_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("fASV*",from) & grepl("fASV*",to)),])
ff_occur_soil <- ff_occur_soil[ff_occur_soil$cor<0,]
nrow(ff_occur_soil)

pp_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("pASV*",from) & grepl("pASV*",to)),])
pp_occur_soil <- pp_occur_soil[pp_occur_soil$cor<0,]
nrow(pp_occur_soil)

bf_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("fASV*",from) & grepl("bASV*",to)),])
bf_occur_soil <- bf_occur_soil[bf_occur_soil$cor<0,]
nrow(bf_occur_soil)

bp_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("pASV*",from) & grepl("bASV*",to)),])
bp_occur_soil <- bp_occur_soil[bp_occur_soil$cor<0,]
nrow(bp_occur_soil)

pf_occur_soil <- droplevels(all_cor_PE_df_padj[with(all_cor_PE_df_padj, grepl("pASV*",from) & grepl("fASV*",to)),])
pf_occur_soil <- pf_occur_soil[pf_occur_soil$cor<0,]
nrow(pf_occur_soil)



mean(E(all_PE_net))

graph.density(all_PE_net)

all_PE_net_cfg <- cluster_fast_greedy(as.undirected(all_PE_net))
all_PE_net_modules <- sort(table(membership(all_PE_net_cfg)),decr=T)

all_PE_modules_20 <- all_PE_net_modules[1:20]
sum(all_PE_modules_20)/sum(all_PE_net_modules)
rm20_plot <- all_PE_modules_20
names(rm20_plot) <- as.factor(1:20)

transitivity(all_PE_net)

all_PE_net_modularity <- modularity(cluster_fast_greedy(all_PE_net))

# ec <- edge.betweenness.community(all_PE_net)
# print(modularity(ec))

all_PE_net_betweenness <- betweenness(all_PE_net,normalized = T)
all_PE_net_degree <- sort(igraph::degree(all_PE_net),decreasing = T)

# all_PE_topological <- cbind(all_PE_net_betweenness,all_PE_net_degree)
# 
# all_PE_topological_data <- as.data.frame(all_PE_topological)
# 
# all_PE_topological <- ggplot(data = all_PE_topological_data,aes(x=all_PE_net_betweenness,y=all_PE_net_degree))+geom_point(colour = "purple", size = 3)+geom_text(aes(label =rownames(all_PE_topological_data)))


plot(all_PE_net_betweenness,all_PE_net_degree)


## Calculate network node degrees/max/min
all_PE_deg <- sort(igraph::degree(all_PE_net,mode="all"),decr=T)
max(all_PE_deg)
mean(all_PE_deg)

## Calculate network node degrees/max/min
all_PE_deg <- sort(igraph::degree(all_PE_net,mode="all"),decr=T)
max(all_PE_deg)
mean(all_PE_deg)




## Combine OTU counts of both kingdoms together
otu_norm_DE_combine <- rbind(otu_norm_DE_B, otu_norm_DE_F, otu_norm_DE_P)

tax_DE_B <- tax_DE_B[,-c(5:9)]
tax_DE_F <- tax_DE_F[,-c(5:9)]
tax_DE_P <- tax_DE_P[,-c(5:10)]

otu_cp_DE_tax <- rbind(tax_DE_B, tax_DE_F, tax_DE_P)

## Perform Spearman correlation of all OTU pairs
all_DE_cor <- rcorr(t(otu_norm_DE_combine), type=c("spearman"))

## Create data frame of co-occurring OTUs
all_cor_DE_df <- CorrDF(all_DE_cor$r, all_DE_cor$P)
all_cor_DE_df$padj <- p.adjust(all_cor_DE_df$p, method="none")

all_cor_DE_df_padj <- all_cor_DE_df[which(all_cor_DE_df$cor > 0.8 | all_cor_DE_df$cor < -0.8),]
all_cor_DE_df_padj <- all_cor_DE_df_padj[which(all_cor_DE_df_padj$padj < 0.001),]

## Make node attribute table
# indic_DE_combine <- c(indic_B_DE, indic_F_DE)

# DE_r_values_combine <- rbind(B_DE_r_values, F_DE_r_values)

nodeattrib_DE_combine <- data.frame(node=union(all_cor_DE_df_padj$from,all_cor_DE_df_padj$to))
nodeattrib_DE_combine$Kingdom <- 0

for (i in as.character(nodeattrib_DE_combine$node))
{
  if (i %in% rownames(otu_cp_DE_tax) == TRUE)
  {nodeattrib_DE_combine[nodeattrib_DE_combine$node==i,"Kingdom"] <- paste(otu_cp_DE_tax[i,1:1])}
  else
  { nodeattrib_DE_combine[nodeattrib_DE_combine$node==i,"Kingdom"]<- "NA"}
}

rownames(nodeattrib_DE_combine) <- as.character(nodeattrib_DE_combine$node)

all_DE_net <- graph_from_data_frame(all_cor_DE_df_padj,direct=F,vertices=nodeattrib_DE_combine)


## Number of nodes
length(V(all_DE_net))

## Number of bacteria and fungi nodes
length(grep("bASV*",names(V(all_DE_net))))
length(grep("fASV*",names(V(all_DE_net))))
length(grep("pASV*",names(V(all_DE_net))))

## Number of edges in network
length(E(all_DE_net))

## Connections 
bb_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_soil <- bb_occur_soil[bb_occur_soil$cor>0,]
nrow(bb_occur_soil)

ff_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("fASV*",from) & grepl("fASV*",to)),])
ff_occur_soil <- ff_occur_soil[ff_occur_soil$cor>0,]
nrow(ff_occur_soil)

pp_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("pASV*",from) & grepl("pASV*",to)),])
pp_occur_soil <- pp_occur_soil[pp_occur_soil$cor>0,]
nrow(pp_occur_soil)

bf_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("fASV*",from) & grepl("bASV*",to)),])
bf_occur_soil <- bf_occur_soil[bf_occur_soil$cor>0,]
nrow(bf_occur_soil)

bp_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("pASV*",from) & grepl("bASV*",to)),])
bp_occur_soil <- bp_occur_soil[bp_occur_soil$cor>0,]
nrow(bp_occur_soil)

pf_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("pASV*",from) & grepl("fASV*",to)),])
pf_occur_soil <- pf_occur_soil[pf_occur_soil$cor>0,]
nrow(pf_occur_soil)

bb_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_soil <- bb_occur_soil[bb_occur_soil$cor<0,]
nrow(bb_occur_soil)

ff_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("fASV*",from) & grepl("fASV*",to)),])
ff_occur_soil <- ff_occur_soil[ff_occur_soil$cor<0,]
nrow(ff_occur_soil)

pp_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("pASV*",from) & grepl("pASV*",to)),])
pp_occur_soil <- pp_occur_soil[pp_occur_soil$cor<0,]
nrow(pp_occur_soil)

bf_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("fASV*",from) & grepl("bASV*",to)),])
bf_occur_soil <- bf_occur_soil[bf_occur_soil$cor<0,]
nrow(bf_occur_soil)

bp_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("pASV*",from) & grepl("bASV*",to)),])
bp_occur_soil <- bp_occur_soil[bp_occur_soil$cor<0,]
nrow(bp_occur_soil)

pf_occur_soil <- droplevels(all_cor_DE_df_padj[with(all_cor_DE_df_padj, grepl("pASV*",from) & grepl("fASV*",to)),])
pf_occur_soil <- pf_occur_soil[pf_occur_soil$cor<0,]
nrow(pf_occur_soil)



mean(E(all_DE_net))

graph.density(all_DE_net)

all_DE_net_cfg <- cluster_fast_greedy(as.undirected(all_DE_net))
all_DE_net_modules <- sort(table(membership(all_DE_net_cfg)),decr=T)

all_DE_modules_20 <- all_DE_net_modules[1:20]
sum(all_DE_modules_20)/sum(all_DE_net_modules)
rm20_plot <- all_DE_modules_20
names(rm20_plot) <- as.factor(1:20)

transitivity(all_DE_net)

all_DE_net_modularity <- modularity(cluster_fast_greedy(all_DE_net))

# ec <- edge.betweenness.community(all_DE_net)
# print(modularity(ec))

all_DE_net_betweenness <- betweenness(all_DE_net,normalized = T)
all_DE_net_degree <- sort(igraph::degree(all_DE_net),decreasing = T)

# all_DE_topological <- cbind(all_DE_net_betweenness,all_DE_net_degree)
# 
# all_DE_topological_data <- as.data.frame(all_DE_topological)
# 
# all_DE_topological <- ggplot(data = all_DE_topological_data,aes(x=all_DE_net_betweenness,y=all_DE_net_degree))+geom_point(colour = "purple", size = 3)+geom_text(aes(label =rownames(all_DE_topological_data)))


plot(all_DE_net_betweenness,all_DE_net_degree)


## Calculate network node degrees/max/min
all_DE_deg <- sort(igraph::degree(all_DE_net,mode="all"),decr=T)
max(all_DE_deg)
mean(all_DE_deg)

## Calculate network node degrees/max/min
all_DE_deg <- sort(igraph::degree(all_DE_net,mode="all"),decr=T)
max(all_DE_deg)
mean(all_DE_deg)




## Combine OTU counts of both kingdoms together
otu_norm_CK_combine <- rbind(otu_norm_CK_B, otu_norm_CK_F, otu_norm_CK_P)

tax_CK_B <- tax_CK_B[,-c(5:9)]
tax_CK_F <- tax_CK_F[,-c(5:9)]
tax_CK_P <- tax_CK_P[,-c(5:10)]

otu_cp_CK_tax <- rbind(tax_CK_B, tax_CK_F, tax_CK_P)

## Perform Spearman correlation of all OTU pairs
all_CK_cor <- rcorr(t(otu_norm_CK_combine), type=c("spearman"))

## Create data frame of co-occurring OTUs
all_cor_CK_df <- CorrDF(all_CK_cor$r, all_CK_cor$P)
all_cor_CK_df$padj <- p.adjust(all_cor_CK_df$p, method="none")

all_cor_CK_df_padj <- all_cor_CK_df[which(all_cor_CK_df$cor > 0.8 | all_cor_CK_df$cor < -0.8),]
all_cor_CK_df_padj <- all_cor_CK_df_padj[which(all_cor_CK_df_padj$padj < 0.001),]

## Make node attribute table
# indic_CK_combine <- c(indic_B_CK, indic_F_CK)

# CK_r_values_combine <- rbind(B_CK_r_values, F_CK_r_values)

nodeattrib_CK_combine <- data.frame(node=union(all_cor_CK_df_padj$from,all_cor_CK_df_padj$to))
nodeattrib_CK_combine$Kingdom <- 0

for (i in as.character(nodeattrib_CK_combine$node))
{
  if (i %in% rownames(otu_cp_CK_tax) == TRUE)
  {nodeattrib_CK_combine[nodeattrib_CK_combine$node==i,"Kingdom"] <- paste(otu_cp_CK_tax[i,1:1])}
  else
  { nodeattrib_CK_combine[nodeattrib_CK_combine$node==i,"Kingdom"]<- "NA"}
}

rownames(nodeattrib_CK_combine) <- as.character(nodeattrib_CK_combine$node)

all_CK_net <- graph_from_data_frame(all_cor_CK_df_padj,direct=F,vertices=nodeattrib_CK_combine)


# sub_graph <- list()
# for (i in names(all_indic_CK)) {
#   sample_i <- all_indic_CK[i]
#   select_node <- rownames(sample_i)[which(sample_i != 0.00000)]
#   sub_graph[[i]] <- induced_subgraph(all_CK_net, select_node)
# }
# sub_graph
# plot_network(sub_graph$CK)
# plot_network(sub_graph$PBAT)
# plot_network(sub_graph$PLA)


## Number of nodes
length(V(all_CK_net))

## Number of bacteria and fungi nodes
length(grep("bASV*",names(V(all_CK_net))))
length(grep("fASV*",names(V(all_CK_net))))
length(grep("pASV*",names(V(all_CK_net))))

## Number of edges in network
length(E(all_CK_net))

## Connections 
bb_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_soil <- bb_occur_soil[bb_occur_soil$cor>0,]
nrow(bb_occur_soil)

ff_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("fASV*",from) & grepl("fASV*",to)),])
ff_occur_soil <- ff_occur_soil[ff_occur_soil$cor>0,]
nrow(ff_occur_soil)

pp_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("pASV*",from) & grepl("pASV*",to)),])
pp_occur_soil <- pp_occur_soil[pp_occur_soil$cor>0,]
nrow(pp_occur_soil)

bf_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("fASV*",from) & grepl("bASV*",to)),])
bf_occur_soil <- bf_occur_soil[bf_occur_soil$cor>0,]
nrow(bf_occur_soil)

bp_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("pASV*",from) & grepl("bASV*",to)),])
bp_occur_soil <- bp_occur_soil[bp_occur_soil$cor>0,]
nrow(bp_occur_soil)

pf_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("pASV*",from) & grepl("fASV*",to)),])
pf_occur_soil <- pf_occur_soil[pf_occur_soil$cor>0,]
nrow(pf_occur_soil)

bb_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_soil <- bb_occur_soil[bb_occur_soil$cor<0,]
nrow(bb_occur_soil)

ff_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("fASV*",from) & grepl("fASV*",to)),])
ff_occur_soil <- ff_occur_soil[ff_occur_soil$cor<0,]
nrow(ff_occur_soil)

pp_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("pASV*",from) & grepl("pASV*",to)),])
pp_occur_soil <- pp_occur_soil[pp_occur_soil$cor<0,]
nrow(pp_occur_soil)

bf_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("fASV*",from) & grepl("bASV*",to)),])
bf_occur_soil <- bf_occur_soil[bf_occur_soil$cor<0,]
nrow(bf_occur_soil)

bp_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("pASV*",from) & grepl("bASV*",to)),])
bp_occur_soil <- bp_occur_soil[bp_occur_soil$cor<0,]
nrow(bp_occur_soil)

pf_occur_soil <- droplevels(all_cor_CK_df_padj[with(all_cor_CK_df_padj, grepl("pASV*",from) & grepl("fASV*",to)),])
pf_occur_soil <- pf_occur_soil[pf_occur_soil$cor<0,]
nrow(pf_occur_soil)



mean(E(all_CK_net))

graph.density(all_CK_net)

all_CK_net_cfg <- cluster_fast_greedy(as.undirected(all_CK_net))
all_CK_net_modules <- sort(table(membership(all_CK_net_cfg)),decr=T)

all_CK_modules_20 <- all_CK_net_modules[1:20]
sum(all_CK_modules_20)/sum(all_CK_net_modules)
rm20_plot <- all_CK_modules_20
names(rm20_plot) <- as.factor(1:20)

transitivity(all_CK_net)

all_CK_net_modularity <- modularity(cluster_fast_greedy(all_CK_net))

# ec <- edge.betweenness.community(all_CK_net)
# print(modularity(ec))

all_CK_net_betweenness <- betweenness(all_CK_net,normalized = T)
all_CK_net_degree <- sort(igraph::degree(all_CK_net),decreasing = T)

# all_CK_topological <- cbind(all_CK_net_betweenness,all_CK_net_degree)
# 
# all_CK_topological_data <- as.data.frame(all_CK_topological)
# 
# all_CK_topological <- ggplot(data = all_CK_topological_data,aes(x=all_CK_net_betweenness,y=all_CK_net_degree))+geom_point(colour = "purple", size = 3)+geom_text(aes(label =rownames(all_CK_topological_data)))


plot(all_CK_net_betweenness,all_CK_net_degree)


## Calculate network node degrees/max/min
all_CK_deg <- sort(igraph::degree(all_CK_net,mode="all"),decr=T)
max(all_CK_deg)
mean(all_CK_deg)

## Calculate network node degrees/max/min
all_CK_deg <- sort(igraph::degree(all_CK_net,mode="all"),decr=T)
max(all_CK_deg)
mean(all_CK_deg)


## Network properties##
igraph_character <- data.frame(c(length(grep("bASV*",names(V(all_CK_net)))), length(grep("fASV*",names(V(all_CK_net)))), length(grep("fASV*",names(V(all_CK_net)))),mean(E(all_CK_net)), graph.density(all_CK_net), transitivity(all_CK_net),length(all_CK_net_modules),all_CK_net_modularity,centralization.degree(all_CK_net)$centralization),
                               c(length(grep("bASV*",names(V(all_PE_net)))), length(grep("fASV*",names(V(all_PE_net)))), length(grep("pASV*",names(V(all_PE_net)))),mean(E(all_PE_net)),graph.density(all_PE_net), transitivity(all_PE_net),length(all_PE_net_modules),all_PE_net_modularity,centralization.degree(all_PE_net)$centralization),
                               c(length(grep("bASV*",names(V(all_DE_net)))), length(grep("fASV*",names(V(all_DE_net)))), length(grep("pASV*",names(V(all_DE_net)))), mean(E(all_DE_net)),graph.density(all_DE_net), transitivity(all_DE_net),length(all_DE_net_modules),all_DE_net_modularity,centralization.degree(all_DE_net)$centralization))



names(igraph_character) <- c("CK","PE","DE")
row.names(igraph_character) <- c("B_node","F_node","P_node","mean","density","transitivity","modules","modularity","centralization")

igraph_character <- as.data.frame(t(igraph_character))
igraph_character$type <- c("CK","PE","DE")

# igraph_character$site <- c(rep(c("soil","root","end","film"),3))
# igraph_character$site <- factor(igraph_character$site,c("soil","root","end","film"))

mean.edge <- ggplot(data = igraph_character, mapping = aes(x = factor(type,levels = c("CK","PE","DE")), y = mean, fill = type)) + geom_bar(stat = 'identity', position = 'dodge')

density <- ggplot(data = igraph_character, mapping = aes(x = factor(type,levels = c("CK","PE","DE")), y =density , fill = type)) + geom_bar(stat = 'identity', position = 'dodge')

transitivity <- ggplot(data = igraph_character, mapping = aes(x = factor(type,levels = c("CK","PE","DE")), y =transitivity  , fill = type)) + geom_bar(stat = 'identity', position = 'dodge')

centralization <- ggplot(data = igraph_character, mapping = aes(x = factor(type,levels = c("CK","PE","DE")), y =centralization  , fill = type)) + geom_bar(stat = 'identity', position = 'dodge')

grid.newpage()
grid.arrange(mean.edge,density,transitivity,centralization,ncol = 2)


write_graph(all_CK_net,"all_CK_net.graphml","graphml")
write_graph(all_PE_net,"all_PE_net.graphml","graphml")
write_graph(all_DE_net,"all_DE_net.graphml","graphml")


all_CK_net_betweenness <- evcent(all_CK_net)$vector
all_PE_net_betweenness <- evcent(all_PE_net)$vector
all_DE_net_betweenness <- evcent(all_DE_net)$vector




all_soil_topological <- cbind(all_CK_net_degree,all_PE_net_degree,all_DE_net_degree)

all_soil_topological <- as.data.frame(all_soil_topological)
colnames(all_soil_topological) <- c("CK","PE","DE")

all_soil_topological1 <- melt(all_soil_topological)

CK  <- all_soil_topological1[all_soil_topological1$variable=="CK", 2 ]
PE  <- all_soil_topological1[all_soil_topological1$variable=="PE", 2 ]
DE  <- all_soil_topological1[all_soil_topological1$variable=="DE", 2 ]
my_comparisons = list( c('CK', 'PE'), c("CK", "DE"), c("PE", "DE"))

all_degree <- ggplot(data = all_soil_topological1,aes(x=factor(variable,levels = c("CK","PE","DE")),y=value))+ geom_boxplot(outlier.shape = NA,aes(fill=variable))+geom_jitter(width =0.2,shape = 21,size=2.5)+  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test")


all_soil_betweenness <- cbind(all_CK_net_betweenness,all_PE_net_betweenness,all_DE_net_betweenness)

all_soil_betweenness <- as.data.frame(all_soil_betweenness)
colnames(all_soil_betweenness) <- c("CK","PE","DE")

all_soil_betweenness1 <- melt(all_soil_betweenness)

CK  <- all_soil_betweenness1[all_soil_betweenness1$variable=="CK", 2 ]
PE  <- all_soil_betweenness1[all_soil_betweenness1$variable=="PE", 2 ]
DE  <- all_soil_betweenness1[all_soil_betweenness1$variable=="DE", 2 ]
my_comparisons = list( c('CK', 'PE'), c("CK", "DE"), c("PE", "DE"))

all_betweenness <- ggplot(data = all_soil_betweenness1,aes(x=factor(variable,levels = c("CK","PE","DE")),y=value))+ geom_boxplot(outlier.shape = NA,aes(fill=variable))+geom_jitter(width =0.2,shape = 21,size=2.5)+  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test")



##### env #####
## Bacteria
###  Env
env <- read.table("env1.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")

env_design <- design_B[rownames(env),]
env$Treatment1 <- env_design$Treatment1

env_one <- env[c(1:24),]
env_two <- env[c(25:48),]

for(i in env[,1:7]) {
  fit <- aov(i~Treatment,data = env)
  print(summary(fit))
  out <- HSD.test(fit,"Treatment")
  print(out$groups)
}

env <- env[,-8]


for(i in env_one[,1:12]) {
  fit <- aov(i~Treatment1,data = env_one)
  print(summary(fit))
  out <- HSD.test(fit,"Treatment1")
  print(out$groups)
}


for(i in env_two[,1:10]) {
  fit <- aov(i~Treatment1,data = env_two)
  print(summary(fit))
  out <- HSD.test(fit,"Treatment1")
  print(out$groups)
}




### RDA
env_log <- log10(env)

### Bacteria RDA
B_soil_otu_rda <- decostand(cp_otu_B_soil,method = "hellinger")
B_soil_otu_rda <- t(B_soil_otu_rda)
B_soil_otu_rda <- as.data.frame(B_soil_otu_rda, scale = FALSE)

B_soil_rda <- rda(B_soil_otu_rda~.,env_log)
B_soil_rda_ill <- summary(B_soil_rda)

vif.cca(B_soil_rda)
anova(B_soil_rda, by = "term", permutations=999)

RsquareAdj(B_soil_rda)
RsquareAdj(rda(B_soil_otu_rda~PH+C+N+P+MC+MN+MP,data=env_log))$adj.r.squared

st=as.data.frame(B_soil_rda_ill$sites)[1:2]
yz=as.data.frame(B_soil_rda_ill$biplot)[1:2]*3

B_soil_rda_plot <- ggplot(data = st,aes(RDA1,RDA2)) +
  geom_point(aes(color = env_design$Treatment1,size=4))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, linewidth=0.6,colour = "blue")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x="RDA1 7.27%",y="RDA2 4.54%")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank())




### RDA one
env_one <- env_one[,-11]

env_one <- env

env_one_log <- log10(env_one)
env_one_design <- design_B[rownames(env_one),]

### Bacteria RDA


B_soil_one_otu_rda <- decostand(cp_otu_B_soil[,rownames(env_one)],method = "hellinger")
B_soil_one_otu_rda <- t(B_soil_one_otu_rda)
B_soil_one_otu_rda <- as.data.frame(B_soil_one_otu_rda, scale = FALSE)

B_soil_one_rda <- rda(B_soil_one_otu_rda~.,env_one_log)
B_soil_one_rda_ill <- summary(B_soil_one_rda)

vif.cca(B_soil_one_rda)
anova(B_soil_one_rda, by = "term", permutations=999)

RsquareAdj(B_soil_one_rda)
RsquareAdj(rda(B_soil_one_otu_rda~PH+C+N+P+MC+MN+MP+BG+CB+GAP,data=env_one_log))$adj.r.squared

st=as.data.frame(B_soil_one_rda_ill$sites)[1:2]
yz=as.data.frame(B_soil_one_rda_ill$biplot)[1:2]*3

B_soil_one_rda_plot <- ggplot(data = st,aes(RDA1,RDA2)) +
  geom_point(aes(color = env_one_design$Treatment1,size=4))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, linewidth=0.6,colour = "blue")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x="RDA1 6.32%",y="RDA2 4.14%")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank())



### RDA two
env_two <- env_two[,-11]
env_two_log <- log10(env_two)
env_two_design <- design_B[rownames(env_two),]

### Bacteria RDA
B_soil_two_otu_rda <- decostand(cp_otu_B_soil[,rownames(env_two)],method = "hellinger")
B_soil_two_otu_rda <- t(B_soil_two_otu_rda)
B_soil_two_otu_rda <- as.data.frame(B_soil_two_otu_rda, scale = FALSE)

B_soil_two_rda <- cca(B_soil_two_otu_rda~.,env_two_log)
B_soil_two_rda_ill <- summary(B_soil_two_rda)

vif.cca(B_soil_two_rda)
anova(B_soil_two_rda, by = "term", permutations=999)

RsquareAdj(B_soil_two_rda)

st=as.data.frame(B_soil_two_rda_ill$sites)[1:2]
yz=as.data.frame(B_soil_two_rda_ill$biplot)[1:2]*3

B_soil_two_rda_plot <- ggplot(data = st,aes(RDA1,RDA2)) +
  geom_point(aes(color = env_two_design$Treatment1,size=4))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, linewidth=0.6,colour = "blue")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x="RDA1 7.81%",y="RDA2 5.75%")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank())



### Fungi RDA
F_soil_otu_rda <- decostand(cp_otu_F_soil,method = "hellinger")
F_soil_otu_rda <- t(F_soil_otu_rda)
F_soil_otu_rda <- as.data.frame(F_soil_otu_rda, scale = FALSE)

F_soil_rda <- rda(F_soil_otu_rda~.,env_log)
F_soil_rda_ill <- summary(F_soil_rda)

vif.cca(F_soil_rda)
anova(F_soil_rda, by = "term", permutations=999)

RsquareAdj(F_soil_rda)
RsquareAdj(rda(F_soil_otu_rda~PH+C+N+P+MC+MN+MP,data=env_log))$adj.r.squared

st=as.data.frame(F_soil_rda_ill$sites)[1:2]
yz=as.data.frame(F_soil_rda_ill$biplot)[1:2]*3

F_soil_rda_plot <- ggplot(data = st,aes(RDA1,RDA2)) +
  geom_point(aes(color = env_design$Treatment1,size=4))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "blue")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x="RDA1 9.56%",y="RDA2 4.65%")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank())



### RDA one
env_one <- env_one[,-11]
env_one_log <- log10(env_one)
env_one_design <- design_B[rownames(env_one),]

### Bacteria RDA
F_soil_one_otu_rda <- decostand(cp_otu_F_soil[,rownames(env_one)],method = "hellinger")
F_soil_one_otu_rda <- t(F_soil_one_otu_rda)
F_soil_one_otu_rda <- as.data.frame(F_soil_one_otu_rda, scale = FALSE)

F_soil_one_rda <- cca(F_soil_one_otu_rda~.,env_one_log)
F_soil_one_rda_ill <- summary(F_soil_one_rda)

vif.cca(F_soil_one_rda)
anova(F_soil_one_rda, by = "term", permutations=999)

RsquareAdj(F_soil_one_rda)

st=as.data.frame(F_soil_one_rda_ill$sites)[1:2]
yz=as.data.frame(F_soil_one_rda_ill$biplot)[1:2]*3

F_soil_one_rda_plot <- ggplot(data = st,aes(RDA1,RDA2)) +
  geom_point(aes(color = env_one_design$Treatment1,size=4))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, linewidth=0.6,colour = "blue")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x="RDA1 2.11%",y="RDA2 1.43%")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank())



### RDA two
env_two <- env_two[,-11]
env_two_log <- log10(env_two)
env_two_design <- design_F[rownames(env_two),]

### FUNGI RDA
F_soil_two_otu_rda <- decostand(cp_otu_F_soil[,rownames(env_two)],method = "hellinger")
F_soil_two_otu_rda <- t(F_soil_two_otu_rda)
F_soil_two_otu_rda <- as.data.frame(F_soil_two_otu_rda, scale = FALSE)

F_soil_two_rda <- rda(F_soil_two_otu_rda~.,env_two_log)
F_soil_two_rda_ill <- summary(F_soil_two_rda)

vif.cca(F_soil_two_rda)
anova(F_soil_two_rda, by = "term", permutations=999)

RsquareAdj(F_soil_two_rda)
RsquareAdj(rda(F_soil_two_otu_rda~PH+C+N+P+MC+MN+MP,data=env_two_log))$adj.r.squared

st=as.data.frame(F_soil_two_rda_ill$sites)[1:2]
yz=as.data.frame(F_soil_two_rda_ill$biplot)[1:2]*3

F_soil_two_rda_plot <- ggplot(data = st,aes(RDA1,RDA2)) +
  geom_point(aes(color = env_two_design$Treatment1,size=4))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, linewidth=0.6,colour = "blue")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x="RDA1 2.38%",y="RDA2 1.62%")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank())




### RDA one
env_one <- env_one[,-11]
env_one_log <- log10(env_one)
env_one_design <- design_B[rownames(env_one),]

### Bacteria RDA
P_soil_one_otu_rda <- decostand(cp_otu_P_soil[,rownames(env_one)],method = "hellinger")
P_soil_one_otu_rda <- t(P_soil_one_otu_rda)
P_soil_one_otu_rda <- as.data.frame(P_soil_one_otu_rda, scale = FALSE)

P_soil_one_rda <- rda(P_soil_one_otu_rda~.,env_one_log)
P_soil_one_rda_ill <- summary(P_soil_one_rda)

vif.cca(P_soil_one_rda)
anova(P_soil_one_rda, by = "term", permutations=999)

RsquareAdj(P_soil_one_rda)
RsquareAdj(rda(P_soil_one_otu_rda~PH+C+N+P+MC+MN+MP,data=env_one_log))$adj.r.squared

st=as.data.frame(P_soil_one_rda_ill$sites)[1:2]
yz=as.data.frame(P_soil_one_rda_ill$biplot)[1:2]*3

P_soil_one_rda_plot <- ggplot(data = st,aes(RDA1,RDA2)) +
  geom_point(aes(color = env_one_design$Treatment1,size=4))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, linewidth=0.6,colour = "blue")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x="RDA1 16.19%",y="RDA2 8.17%")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank())



### RDA two
env_two <- env_two[,-11]
env_two_log <- log10(env_two)
env_two_design <- design_P[rownames(env_two),]

### Protist RDA
P_soil_two_otu_rda <- decostand(cp_otu_P_soil[,rownames(env_two)],method = "hellinger")
P_soil_two_otu_rda <- t(P_soil_two_otu_rda)
P_soil_two_otu_rda <- as.data.frame(P_soil_two_otu_rda, scale = FALSE)

P_soil_two_rda <- rda(P_soil_two_otu_rda~.,env_two_log)
P_soil_two_rda_ill <- summary(P_soil_two_rda)

vif.cca(P_soil_two_rda)
anova(P_soil_two_rda, by = "term", permutations=999)

RsquareAdj(P_soil_two_rda)
RsquareAdj(rda(P_soil_two_otu_rda~PH+C+N+P+MC+MN+MP,data=env_two_log))$adj.r.squared

st=as.data.frame(P_soil_two_rda_ill$sites)[1:2]
yz=as.data.frame(P_soil_two_rda_ill$biplot)[1:2]*3

P_soil_two_rda_plot <- ggplot(data = st,aes(RDA1,RDA2)) +
  geom_point(aes(color = env_two_design$Treatment1,size=4))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, linewidth=0.6,colour = "blue")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x="RDA1 12.30%",y="RDA2 10.65%")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank())


grid.newpage()
grid.arrange(B_soil_one_rda_plot,B_soil_two_rda_plot,F_soil_one_rda_plot,F_soil_two_rda_plot,P_soil_one_rda_plot,P_soil_two_rda_plot,nrow = 2)


#envfit  
B_one_ef=envfit(B_soil_one_rda,env_one_log,permu=999)
B_one_ef$vectors$r#
B_one_ef$vectors$pvals


B_one_rda_envfit <- data.frame(tax=colnames(env_one_log),B.r=B_one_ef$vectors$r,B.p=B_one_ef$vectors$pvals)
B_one_com=arrange(B_one_rda_envfit,B.r)
head(B_one_com,n=3)
B_one_com[c(3)]=B_one_com[c(3)]>0.05

rownames(B_one_com) <- B_one_com$tax
head(B_one_com,n=3)

B_one_com$tax = factor(B_one_com$tax,order = T,levels = row.names(B_one_com))

B_one_com_p2 <- ggplot(B_one_com, aes(x =tax, y = B.r),size=2) +
  geom_bar(stat = 'identity', width = 0.8,color="black",fill="red") +
  scale_fill_manual(guide = FALSE)+
  geom_text(aes(y = B.r+0.005, label = ifelse(B.p==TRUE,"","*")),
            size = 5, fontface = "bold") +
  xlab("Environmental factor")+
  ylab(expression(r^"2"))+
  scale_y_continuous(expand = c(0,0))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))
B_one_com_p2



B_two_ef=envfit(B_soil_two_rda,env_two_log,permu=999)
B_two_ef$vectors$r
B_two_ef$vectors$pvals


B_two_rda_envfit <- data.frame(tax=colnames(env_two_log),B.r=B_two_ef$vectors$r,B.p=B_two_ef$vectors$pvals)
B_two_com=arrange(B_two_rda_envfit,B.r)
head(B_two_com,n=3)
B_two_com[c(3)]=B_two_com[c(3)]>0.05

rownames(B_two_com) <- B_two_com$tax
head(B_two_com,n=3)

B_two_com$tax = factor(B_two_com$tax,order = T,levels = row.names(B_two_com))

B_two_com_p2 <- ggplot(B_two_com, aes(x =tax, y = B.r),size=2) +
  geom_bar(stat = 'identity', width = 0.8,color="black",fill="red") +
  scale_fill_manual(guide = FALSE)+
  geom_text(aes(y = B.r+0.005, label = ifelse(B.p==TRUE,"","*")),
            size = 5, fontface = "bold") +
  xlab("Environmental factor")+
  ylab(expression(r^"2"))+
  scale_y_continuous(expand = c(0,0))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))
B_two_com_p2


F_one_ef=envfit(F_soil_one_rda,env_one_log,permu=999)
F_one_ef$vectors$r
F_one_ef$vectors$pvals


F_one_rda_envfit <- data.frame(tax=colnames(env_one_log),B.r=F_one_ef$vectors$r,B.p=F_one_ef$vectors$pvals)
F_one_com=arrange(F_one_rda_envfit,B.r)
head(F_one_com,n=3)
F_one_com[c(3)]=F_one_com[c(3)]>0.05

rownames(F_one_com) <- F_one_com$tax
head(F_one_com,n=3)

F_one_com$tax = factor(F_one_com$tax,order = T,levels = row.names(F_one_com))

F_one_com_p2 <- ggplot(F_one_com, aes(x =tax, y = B.r),size=2) +
  geom_bar(stat = 'identity', width = 0.8,color="black",fill="red") +
  scale_fill_manual(guide = FALSE)+
  geom_text(aes(y = B.r+0.005, label = ifelse(B.p==TRUE,"","*")),
            size = 5, fontface = "bold") +
  xlab("Environmental factor")+
  ylab(expression(r^"2"))+
  scale_y_continuous(expand = c(0,0))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))
F_one_com_p2



F_two_ef=envfit(F_soil_two_rda,env_two_log,permu=999)
F_two_ef$vectors$r
F_two_ef$vectors$pvals


F_two_rda_envfit <- data.frame(tax=colnames(env_two_log),B.r=F_two_ef$vectors$r,B.p=F_two_ef$vectors$pvals)
F_two_com=arrange(F_two_rda_envfit,B.r)
head(F_two_com,n=3)
F_two_com[c(3)]=F_two_com[c(3)]>0.05

rownames(F_two_com) <- F_two_com$tax
head(F_two_com,n=3)

F_two_com$tax = factor(F_two_com$tax,order = T,levels = row.names(F_two_com))

F_two_com_p2 <- ggplot(F_two_com, aes(x =tax, y = B.r),size=2) +
  geom_bar(stat = 'identity', width = 0.8,color="black",fill="red") +
  scale_fill_manual(guide = FALSE)+
  geom_text(aes(y = B.r+0.005, label = ifelse(B.p==TRUE,"","*")),
            size = 5, fontface = "bold") +
  xlab("Environmental factor")+
  ylab(expression(r^"2"))+
  scale_y_continuous(expand = c(0,0))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))
F_two_com_p2

P_one_ef=envfit(P_soil_one_rda,env_one_log,permu=999)
P_one_ef$vectors$r
P_one_ef$vectors$pvals


P_one_rda_envfit <- data.frame(tax=colnames(env_one_log),B.r=P_one_ef$vectors$r,B.p=P_one_ef$vectors$pvals)
P_one_com=arrange(P_one_rda_envfit,B.r)
head(P_one_com,n=3)
P_one_com[c(3)]=P_one_com[c(3)]>0.05

rownames(P_one_com) <- P_one_com$tax
head(P_one_com,n=3)

P_one_com$tax = factor(P_one_com$tax,order = T,levels = row.names(P_one_com))

P_one_com_p2 <- ggplot(P_one_com, aes(x =tax, y = B.r),size=2) +
  geom_bar(stat = 'identity', width = 0.8,color="black",fill="red") +
  scale_fill_manual(guide = FALSE)+
  geom_text(aes(y = B.r+0.005, label = ifelse(B.p==TRUE,"","*")),
            size = 5, fontface = "bold") +
  xlab("Environmental factor")+
  ylab(expression(r^"2"))+
  scale_y_continuous(expand = c(0,0))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))
P_one_com_p2



P_two_ef=envfit(P_soil_two_rda,env_two_log,permu=999)
P_two_ef$vectors$r
P_two_ef$vectors$pvals


P_two_rda_envfit <- data.frame(tax=colnames(env_two_log),B.r=P_two_ef$vectors$r,B.p=P_two_ef$vectors$pvals)
P_two_com=arrange(P_two_rda_envfit,B.r)
head(P_two_com,n=3)
P_two_com[c(3)]=P_two_com[c(3)]>0.05

rownames(P_two_com) <- P_two_com$tax
head(P_two_com,n=3)

P_two_com$tax = factor(P_two_com$tax,order = T,levels = row.names(P_two_com))

P_two_com_p2 <- ggplot(P_two_com, aes(x =tax, y = B.r),size=2) +
  geom_bar(stat = 'identity', width = 0.8,color="black",fill="red") +
  scale_fill_manual(guide = FALSE)+
  geom_text(aes(y = B.r+0.005, label = ifelse(B.p==TRUE,"","*")),
            size = 5, fontface = "bold") +
  xlab("Environmental factor")+
  ylab(expression(r^"2"))+
  scale_y_continuous(expand = c(0,0))+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 45))
P_two_com_p2

grid.newpage()
grid.arrange(B_one_com_p2,F_one_com_p2,P_one_com_p2,B_two_com_p2,F_two_com_p2,P_two_com_p2,nrow = 2)



##### Metagenom ####
metadata1 <- read.table("metadata1.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
## KEGG
Kegg_ko <- read.table("Kegg_tab.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")
Kegg_ko[is.na(Kegg_ko)] <- 0
Kegg_tax <- read.table("Kegg_tax.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")

edgeR_kegg <- DGEList(counts=Kegg_ko, 
                      group=metadata1$Treatment1)

edgeR_kegg <- calcNormFactors(edgeR_kegg)

otu_norm_kegg <- cpm(edgeR_kegg, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk root bacteria community into phyloseq objects
## for further analysis
phy_kegg <- otu_table(Kegg_ko,taxa_are_rows=T)
phy_tax_kegg <-tax_table(as.matrix(Kegg_tax))
phy_design_kegg <- sample_data(metadata1)
physeq_kegg <- phyloseq(phy_kegg,phy_design_kegg)


## Create bray-curtis dissimiliartiy matrix
all_dis_kegg <- vegdist(t(otu_table(physeq_kegg)),method="bray")


## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_kegg <- adonis2(all_dis_kegg ~Treatment1, data=metadata1, permutations=9999)

paov_all_kegg  



## NR
NR <- read.table("NR_tab.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
NR_tax <- read.table("NR_tax.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")

NR_tax_B <- NR_tax[which(NR_tax$Kingdom=="Bacteria"),]
NR_B <- NR[rownames(NR_tax_B),]

# create separate taxonomy label specifying classes of Proteobacteria
NR_tax_B$labels <- NR_tax_B$Phylum
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$Class=="Alphaproteobacteria" ], ]$labels <- "Alphaproteobacteria"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$Class=="Betaproteobacteria" ], ]$labels <- "Betaproteobacteria"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$Class=="Gammaproteobacteria" ], ]$labels <- "Gammaproteobacteria"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$Class=="Deltaproteobacteria" ], ]$labels <- "Deltaproteobacteria"
table(NR_tax_B$labels)

### Defining bOTU colors by phylum (using the taxonomy file)
NR_tax_B$cols <- NR_tax_B$labels
table(NR_tax_B$cols)

edgeR_NR <- DGEList(counts=NR, 
                    group=metadata1$Treatment1,
                    genes=NR_tax)

edgeR_NR <- calcNormFactors(edgeR_NR)

otu_norm_NR <- cpm(edgeR_NR, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk root bacteria community into phyloseq objects
## for further analysis
phy_NR <- otu_table(NR,taxa_are_rows=T)
phy_tax_NR <-tax_table(as.matrix(NR_tax))
phy_design_NR <- sample_data(metadata1)
physeq_NR <- phyloseq(phy_NR,phy_design_kegg,phy_tax_NR)


## Create bray-curtis dissimiliartiy matrix
all_dis_NR <- vegdist(t(otu_table(physeq_NR)),method="bray")

## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_NR <- adonis2(all_dis_NR ~Treatment1, data=metadata1, permutations=9999)

paov_all_NR  



## CAZy
CAZy <- read.table("CAZy_tab.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
CAZy_tax <- read.table("CAZy_tax.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)

edgeR_CAZy <- DGEList(counts=CAZy, 
                      group=metadata1$type,
                      genes=CAZy_tax)

edgeR_CAZy <- calcNormFactors(edgeR_CAZy)

otu_norm_CAZy <- cpm(edgeR_CAZy, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk root bacteria community into phyloseq objects
## for further analysis
phy_CAZy <- otu_table(CAZy,taxa_are_rows=T)
phy_tax_CAZy <-tax_table(as.matrix(CAZy_tax))
phy_design_CAZy <- sample_data(metadata1)
physeq_CAZy <- phyloseq(phy_CAZy,phy_design_CAZy,phy_tax_CAZy)


## Create bray-curtis dissimiliartiy matrix
all_dis_CAZy <- vegdist(t(otu_table(physeq_CAZy)),method="bray")

## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_CAZy <- adonis2(all_dis_CAZy ~Treatment1, data=metadata1, permutations=9999)

paov_all_CAZy  



## COG
COG <- read.table("COG_tab.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
COG_tax <- read.table("COG_tax.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)

edgeR_COG <- DGEList(counts=COG, 
                     group=metadata1$type)

edgeR_COG <- calcNormFactors(edgeR_COG)

otu_norm_COG <- cpm(edgeR_COG, normalized.lib.sizes=T, log=F)

## Input TMM normalized counts, taxonomy, and design of bulk root bacteria community into phyloseq objects
## for further analysis
phy_COG <- otu_table(COG,taxa_are_rows=T)
phy_design_COG <- sample_data(metadata1)
physeq_COG <- phyloseq(phy_COG,phy_design_COG)


## Create bray-curtis dissimiliartiy matrix
all_dis_COG <- vegdist(t(otu_table(physeq_COG)),method="bray")

## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_COG <- adonis2(all_dis_COG ~Treatment1, data=metadata1, permutations=9999)

paov_all_COG




##kegg_diversity
pcoa_kegg <- ordinate(physeq_kegg,"PCoA","bray")
pcoa_all_kegg <- plot_ordination(physeq_kegg, pcoa_kegg, type="sites", color ="Treatment1")
pcoa_all_kegg1 <- pcoa_all_kegg+
  geom_point(size=5)+
  xlab(paste("PCo 1", paste("(",round(pcoa_all_kegg$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_all_kegg$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("KEGG")



set.seed(619)
cp_kegg <- as.data.frame(t(rrarefy(t(Kegg_ko), min(colSums(Kegg_ko)))))
Kegg_shannon <- vegan::diversity(t(cp_kegg), index = "shannon")

chao1 <- t(estimateR(t(cp_kegg)))
chao1 <- as.data.frame(chao1)

metadata1$shannon <- Kegg_shannon
metadata1$chao1 <- chao1$S.chao1

CK  <- metadata1[metadata1$shannon=="CK", 8 ]
PE  <- metadata1[metadata1$shannon=="PE", 8 ]
DE  <- metadata1[metadata1$shannon=="DE", 8 ]

my_comparisons = list( c('CK', 'PE'), c("CK", "DE"), c("PE", "DE"))

Kegg_shannon <- ggplot(metadata1, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=shannon, fill =Treatment1))+geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#3CB371","#DAA520","gray30"))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+  labs(x="Types", y="Kegg Shannon index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")


## Kegg_de
CK_samples <- rownames(metadata1)[which(metadata1$Treatment1 == "CK")]
PE_samples <- rownames(metadata1)[which(metadata1$Treatment1 == "PE")]
DE_samples <- rownames(metadata1)[which(metadata1$Treatment1 == "DE")]

Kegg_de_tax <- read.table("kegg_de_tax.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
Kegg_ko_de <- cp_kegg[rownames(Kegg_de_tax),]

Kegg_ko_de_mat <- as.data.frame(t(Kegg_ko_de))
par(mfrow=c(1,2), mar=c(0.5,3.5,2,0))

CS_cols <- c("#2c7fbe","#f7d368","#8bc224")

metadata1$Treatment1<-factor(metadata1$Treatment1,c("CK","PE","DE"))

bargraph.CI(metadata1$Treatment1, Kegg_ko_de_mat$K16320,
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)

bargraph.CI(metadata1$Treatment1,Kegg_ko_de_mat$K23359,
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)



bargraph.CI(metadata1$Treatment1,Kegg_ko_de_mat$K00448,
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)


bargraph.CI(metadata1$Treatment1,Kegg_ko_de_mat$K10676,
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)




## Get names of bacteria phyla present (use 'labels' as this specifies class within Proteobacteria)
Pathway_Kegg_de <- names(sort(table(Kegg_de_tax[,"Pathway"]), decr=T))
length(Pathway_Kegg_de)
sort(table(Kegg_de_tax[,"Pathway"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_ko_de)
for (i in Pathway_Kegg_de){
  x <- array(colSums(Kegg_ko_de[rownames(Kegg_de_tax)[which(Kegg_de_tax$Pathway == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(Pathway_Kegg_de)
colnames(y) <- paste(colnames(Kegg_ko_de))
Pathway_mat_Kegg_de <- y
Pathway_mat_Kegg_de[,1:5]
colSums(Pathway_mat_Kegg_de)
PHYLUM_mat_B_mean <- sort(apply(Pathway_mat_Kegg_de,1,mean),decr=F)
Pathway_mat_Kegg_de <- Pathway_mat_Kegg_de[names(PHYLUM_mat_B_mean),]


Pathway_mat_Kegg_de_ra <- t(t(Pathway_mat_Kegg_de)/colSums(Pathway_mat_Kegg_de)) *100

PHYLUM_mat_flowering_mean <- sort(apply(Pathway_mat_Kegg_de_ra,1,mean),decr=T)
Pathway_mat_Kegg_de_ra <- Pathway_mat_Kegg_de_ra[names(PHYLUM_mat_flowering_mean),]



Pathway_mat_Kegg_de_ra <- as.data.frame(t(Pathway_mat_Kegg_de))


bargraph.CI(metadata1$Treatment1, rowSums(Pathway_mat_Kegg_de_ra), 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025, ylim=c(39000,43000), border=F)




Pathway_treat_de_ra <-as.matrix(cbind(`CK`=apply(Pathway_mat_Kegg_de_ra[,CK_samples],1,mean),
                                      `PE`=apply(Pathway_mat_Kegg_de_ra[,PE_samples],1,mean),
                                      `DE`=apply(Pathway_mat_Kegg_de_ra[,DE_samples],1,mean)))


## Define soil samples by management system
Pathway_treat <-as.matrix(cbind(`CK`=apply(Pathway_mat_Kegg_de[,CK_samples],1,mean),
                                `PE`=apply(Pathway_mat_Kegg_de[,PE_samples],1,mean),
                                `DE`=apply(Pathway_mat_Kegg_de[,DE_samples],1,mean)))

Pathway_treat_mean <- sort(apply(Pathway_treat,1,mean),decr=T)
Pathway_treat <- Pathway_treat[names(Pathway_treat_mean),]

Pathway_treat <- as.matrix(Pathway_treat)

Pathway_color <- c("lightcoral","darkolivegreen2","paleturquoise2","lightgoldenrod1","darkorange3","orchid1","darkslateblue","royalblue1","yellow","limegreen","purple","tomato","white","palegreen3","steelblue3","cornflowerblue","skyblue3","sienna4","mediumslateblue","orange","orchid")

##### Plot Supplementary Figure S3
par(oma=c(0,0,0,0), mar=c(6,4,1,5), xpd=NA)
phylum_bar_B <- barplot(as.matrix(Pathway_treat), col=Pathway_color, xaxt="n", border=NA, las=2)
legend(20,4000, bty="n", cex=0.7, x.intersp= 0.1, y.intersp=1,
       legend=rev(rownames(Pathway_treat)), 
       fill=rev(Pathway_color), 
       border=rev(Pathway_color))

library(gg.gap)
library(ggplot2)
data <-
  data.frame(x = c("CK","PE","DE"),
             y = c(40416.00,41434.67,42636.50))

p1 = ggplot(data, aes(x = x, y = y, fill = x)) +
  geom_bar(stat = 'identity', position = position_dodge(),show.legend = FALSE) +
  theme_bw() +
  labs(x = NULL, y = NULL)+
  theme_minimal()+theme(legend.position = "none")


p1 = p1 + theme_classic()+coord_flip()
p2 = gg.gap(plot = p1,
            segments = c(1500,39000),
            tick_width = 1000,
            rel_heights = c(0.2,0,0.2,0,1),
            ylim = c(0, 43000)
)

p2



##NR_diversity
pcoa_NR <- ordinate(physeq_NR,"PCoA","bray")
pcoa_all_NR <- plot_ordination(physeq_NR, pcoa_NR, type="sites", color ="Treatment1")
pcoa_all_NR1 <- pcoa_all_NR+
  geom_point(size=5)+
  xlab(paste("PCo 1", paste("(",round(pcoa_all_NR$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_all_NR$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("NR")


set.seed(619)
cp_NR <- as.data.frame(t(rrarefy(t(NR), min(colSums(NR)))))
NR_shannon <- vegan::diversity(t(cp_NR), index = "shannon")

chao1 <- t(estimateR(t(cp_NR)))
chao1 <- as.data.frame(chao1)

metadata1$NR_shannon <- NR_shannon
metadata1$NR_chao1 <- chao1$S.chao1


NR_diversity <- ggplot(metadata1, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=NR_shannon, fill =Treatment1))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#3CB371","#DAA520","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Types", y="NR shannon index")


NR_shannon <- ggplot(metadata1, aes(x=factor(Treatment2,levels = c("soil_NO","soil_ADD","film_ADD","root_NO","root_ADD")), y=NR_shannon, fill =Treatment2))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("white","#3CB371","#DAA520","gray30","blue"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Types", y="NR Shannon index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")





## Soil bacteria: Calculate percentage of sequences classified into each phylum
NR_B_RA <- t(t(NR_B)/colSums(NR_B)) * 100
colSums(NR_B_RA)
nrow(NR_B_RA)

## Get names of bacteria phyla present (use 'labels' as this specifies class within Proteobacteria)
PHYLAnames_B <- names(sort(table(NR_tax_B[,"labels"]), decr=T))
length(PHYLAnames_B)
sort(table(NR_tax_B[,"labels"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(NR_B_RA)
for (i in PHYLAnames_B){
  x <- array(colSums(NR_B_RA[rownames(NR_tax_B)[which(NR_tax_B$labels == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}


## Create matrix
rownames(y) <- paste(PHYLAnames_B)
colnames(y) <- paste(colnames(NR_B_RA))
PHYLUM_mat_B <- y
PHYLUM_mat_B[,1:5]
colSums(PHYLUM_mat_B)
PHYLUM_mat_B_mean <- sort(apply(PHYLUM_mat_B,1,mean),decr=T)
PHYLUM_mat_B <- PHYLUM_mat_B[names(PHYLUM_mat_B_mean),]


new_B_cols<-c("CK1","CK2","CK3","PE1","PE2","PE3","PBAT1","PBAT2","PBAT3","PLA1","PLA2","PLA3")

PHYLUM_mat_B_mean_type <- PHYLUM_mat_B[, c(new_B_cols)]


# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_B_mean_type, 1, mean) < 1)
low_count_phyla_B <- rownames(PHYLUM_mat_B_mean_type)[sort(apply(PHYLUM_mat_B_mean_type, 1, mean), decr=T) < 1]
# attribute grey color
for(i in low_count_phyla_B){
  NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels==paste(i) ], ]$cols <- "lightgrey"
}
table(NR_tax_B$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
abundant_phyla_B <- rownames(PHYLUM_mat_B_mean_type)[sort(apply(PHYLUM_mat_B_mean_type, 1, mean), decr=T) > 1]
abundant_phyla_B
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Betaproteobacteria" ], ]$cols <- "lightcoral"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Gammaproteobacteria" ], ]$cols <- "palegreen3"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Alphaproteobacteria" ], ]$cols <- "palegreen1"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Acidobacteria" ], ]$cols <- "palegreen4"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Actinobacteria" ], ]$cols <- "dodgerblue"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Bacteroidetes" ], ]$cols <- "sandybrown"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Chloroflexi" ], ]$cols <- "plum1"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Deltaproteobacteria" ], ]$cols <- "palegreen2"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Planctomycetes" ], ]$cols <- "orchid"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Gemmatimonadetes" ], ]$cols <- "steelblue4"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Verrucomicrobia" ], ]$cols <- "steelblue1"
NR_tax_B[ rownames(NR_tax_B)[NR_tax_B$labels=="Proteobacteria" ], ]$cols <- "orange"


## collaps OTU colors to prepare Phylum level colors
label_cols_B <- NR_tax_B[, c("labels", "cols") ]
library(plyr)
PHYLA_label_cols_B <- ddply(label_cols_B, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_B) <- PHYLA_label_cols_B[,1]
PHYLA_label_cols_B <- PHYLA_label_cols_B[c(abundant_phyla_B, low_count_phyla_B),]
PHYLA_label_cols_B

## Legend for Phylum colors
PHYLA_label_cols_B_legend <- PHYLA_label_cols_B[1:12,]
PHYLA_label_cols_B_legend[11,1] <- "other"
rownames(PHYLA_label_cols_B_legend)[12] <- "other"
PHYLA_label_cols_B_legend


##### Plot Supplementary Figure S3
par(oma=c(0,0,0,0), mar=c(6,4,1,5), xpd=NA)
phylum_bar_B <- barplot(as.matrix(PHYLUM_mat_B_mean_type), col=PHYLA_label_cols_B[rownames(PHYLUM_mat_B_mean_type),]$cols, ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_B, labels=colnames(PHYLUM_mat_B_mean_type), col.axis="black", las=2, cex.axis=1)
title(ylab="Relative abundance (%)")
title(main="Bacteria Community")
legend(14.5, 100, bty="n", cex=0.7, x.intersp= 0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_B_legend$labels), 
       fill=rev(PHYLA_label_cols_B_legend$cols), 
       border=rev(PHYLA_label_cols_B_legend$cols) )







##CAZy_diversity
pcoa_CAZy <- ordinate(physeq_CAZy,"PCoA","bray")
pcoa_all_CAZy <- plot_ordination(physeq_CAZy, pcoa_CAZy, type="sites", color ="Treatment1")
pcoa_all_CAZy1 <- pcoa_all_CAZy+
  geom_point(size=5)+
  xlab(paste("PCo 1", paste("(",round(pcoa_all_CAZy$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_all_CAZy$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("CAZy")


set.seed(619)
cp_CAZy <- as.data.frame(t(rrarefy(t(CAZy), min(colSums(CAZy)))))
CAZy_shannon <- vegan::diversity(t(cp_CAZy), index = "shannon")

chao1 <- t(estimateR(t(cp_CAZy)))
chao1 <- as.data.frame(chao1)

metadata1$CAZy_shannon <- CAZy_shannon
metadata1$CAZy_chao1 <- chao1$S.chao1

CK  <- metadata1[metadata1$shannon=="CK", 10 ]
PE  <- metadata1[metadata1$shannon=="PE", 10 ]
DE  <- metadata1[metadata1$shannon=="DE", 10 ]

my_comparisons = list( c('CK', 'PE'), c("CK", "DE"), c("PE", "DE"))

CAZy_shannon <- ggplot(metadata1, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=CAZy_shannon, fill =Treatment1))+geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#3CB371","#DAA520","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+labs(x="Types", y="CAZy shannon index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



##COG_diversity
pcoa_COG <- ordinate(physeq_COG,"PCoA","bray")
pcoa_all_COG <- plot_ordination(physeq_COG, pcoa_COG, type="sites", color ="Treatment1")
pcoa_all_COG1 <- pcoa_all_COG+
  geom_point(size=5)+
  xlab(paste("PCo 1", paste("(",round(pcoa_all_COG$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_all_COG$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("COG")

set.seed(619)
cp_COG <- as.data.frame(t(rrarefy(t(COG), min(colSums(COG)))))
COG_shannon <- vegan::diversity(t(cp_COG), index = "shannon")

chao1 <- t(estimateR(t(cp_COG)))
chao1 <- as.data.frame(chao1)

metadata1$COG_shannon <- COG_shannon
metadata1$COG_chao1 <- chao1$S.chao1

CK  <- metadata1[metadata1$shannon=="CK", 12 ]
PE  <- metadata1[metadata1$shannon=="PE", 12 ]
DE  <- metadata1[metadata1$shannon=="DE", 12 ]

my_comparisons = list( c('CK', 'PE'), c("CK", "DE"), c("PE", "DE"))

COG_diversity <- ggplot(metadata1, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=COG_shannon, fill =Treatment1))+geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#3CB371","#DAA520","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  labs(x="Types", y="COG shannon index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



grid.newpage()
grid.arrange(Kegg_shannon,CAZy_shannon,COG_diversity,ncol = 3)

grid.newpage()
grid.arrange(pcoa_all_kegg1,pcoa_all_CAZy1,pcoa_all_COG1,ncol = 3)


for(i in metadata1[,4:9]) {
  fit <- aov(i~type,data = metadata1)
  print(summary(fit))
  out <- HSD.test(fit,"type")
  print(out$groups)
}



##### sfd #####
sfd <- read.csv("metadata2.csv",header = T,row.names = 1)

B_KO_SFD <- ggplot(sfd, aes(x = KO_shannon, y = B_Shannon)) +
  geom_point(aes(color = Treatment1),shape=16,size=4) +
  geom_smooth(method = "lm", color = '#08141E') +  scale_color_manual(values = c('CK' = '#2c7fbe', 'PE' = '#f7d368','DE' = '#8bc224'))+stat_poly_eq(aes(label = paste(eq.label, ..adj.rr.label.., sep = '~~~~')), parse = T,family = "SH")

B_COG_SFD <- ggplot(sfd, aes(x = COG_shannon, y = B_Shannon)) +
  geom_point(aes(color = Treatment1),shape=16,size=4) +
  geom_smooth(method = "lm", color = '#08141E') +  scale_color_manual(values = c('CK' = '#2c7fbe', 'PE' = '#f7d368','DE' = '#8bc224'))+stat_poly_eq(aes(label = paste(eq.label, ..adj.rr.label.., sep = '~~~~')), parse = T,family = "SH")

B_CAZy_SFD <- ggplot(sfd, aes(x = CAZy_shannon, y = B_Shannon)) +
  geom_point(aes(color = Treatment1),shape=16,size=4) +
  geom_smooth(method = "lm", se = T, color = '#08141E') +  scale_color_manual(values = c('CK' = '#2c7fbe', 'PE' = '#f7d368','DE' = '#8bc224'))+stat_poly_eq(aes(label = paste(eq.label, ..adj.rr.label.., sep = '~~~~')), parse = T,family = "SH")


F_KO_SFD <- ggplot(sfd, aes(x = KO_shannon, y = F_Shannon)) +
  geom_point(aes(color = Treatment1),shape=16,size=4) +
  geom_smooth(method = "lm", se = T, color = '#08141E') +  scale_color_manual(values = c('CK' = '#2c7fbe', 'PE' = '#f7d368','DE' = '#8bc224'))+stat_poly_eq(aes(label = paste(eq.label, ..adj.rr.label.., sep = '~~~~')), parse = T,family = "SH")

F_COG_SFD <- ggplot(sfd, aes(x = COG_shannon, y = F_Shannon)) +
  geom_point(aes(color = Treatment1),shape=16,size=4) +
  geom_smooth(method = "lm",se = T, color = '#08141E') +  scale_color_manual(values = c('CK' = '#2c7fbe', 'PE' = '#f7d368','DE' = '#8bc224'))+stat_poly_eq(aes(label = paste(eq.label, ..adj.rr.label.., sep = '~~~~')), parse = T,family = "SH")

F_CAZy_SFD <- ggplot(sfd, aes(x = CAZy_shannon, y = F_Shannon)) +
  geom_point(aes(color = Treatment1),shape=16,size=4) +
  geom_smooth(method = "lm", se = T, color = '#08141E') +  scale_color_manual(values = c('CK' = '#2c7fbe', 'PE' = '#f7d368','DE' = '#8bc224'))+stat_poly_eq(aes(label = paste(eq.label, ..adj.rr.label.., sep = '~~~~')), parse = T,family = "SH")



P_KO_SFD <- ggplot(sfd, aes(x = KO_shannon, y = P_Shannon)) +
  geom_point(aes(color = Treatment1),shape=16,size=4) +
  geom_smooth(method = "lm", se = T, color = '#08141E') +  scale_color_manual(values = c('CK' = '#2c7fbe', 'PE' = '#f7d368','DE' = '#8bc224'))+stat_poly_eq(aes(label = paste(eq.label, ..adj.rr.label.., sep = '~~~~')), parse = T,family = "SH")

P_COG_SFD <- ggplot(sfd, aes(x = COG_shannon, y = P_Shannon)) +
  geom_point(aes(color = Treatment1),shape=16,size=4) +
  geom_smooth(method = "lm",se = T, color = '#08141E') +  scale_color_manual(values = c('CK' = '#2c7fbe', 'PE' = '#f7d368','DE' = '#8bc224'))+stat_poly_eq(aes(label = paste(eq.label, ..adj.rr.label.., sep = '~~~~')), parse = T,family = "SH")

P_CAZy_SFD <- ggplot(sfd, aes(x = CAZy_shannon, y = P_Shannon)) +
  geom_point(aes(color = Treatment1),shape=16,size=4) +
  geom_smooth(method = "lm", se = T, color = '#08141E') +  scale_color_manual(values = c('CK' = '#2c7fbe', 'PE' = '#f7d368','DE' = '#8bc224'))+stat_poly_eq(aes(label = paste(eq.label, ..adj.rr.label.., sep = '~~~~')), parse = T,family = "SH")


grid.newpage()
grid.arrange(B_KO_SFD ,B_COG_SFD ,B_CAZy_SFD,F_KO_SFD ,F_COG_SFD ,F_CAZy_SFD,P_KO_SFD ,P_COG_SFD ,P_CAZy_SFD,ncol = 3)


library(patchwork)
B_KO_SFD+B_COG_SFD +B_CAZy_SFD+F_KO_SFD+F_COG_SFD+F_CAZy_SFD+P_KO_SFD+P_COG_SFD +P_CAZy_SFD+plot_layout(guides = "collect")





## C/N/S/A metaboilte
## C metabilte
C_tax <- read.table("C_tax.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
Kegg_Ctab <- cp_kegg[rownames(C_tax),]

Kegg_Ctab$KEGG_Name <- C_tax$KEGG_Name

Kegg_Ctab <- melt (Kegg_Ctab)

func <- function(x)(c(n = length(x),mean=mean(x,na.rm = T)))
Kegg_Ctab <- aggregate(value ~ KEGG_Name + variable, data=Kegg_Ctab,func)

ggplot(Kegg_Ctab, aes(x = variable, y = KEGG_Name, size = value[,"mean"], color=variable)) + geom_point()


## N metabilte
N_tax <- read.table("N_tax.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
Kegg_Ntab <- cp_kegg[rownames(N_tax),]

Kegg_Ntab$KEGG_Name <- N_tax$KEGG_Name

Kegg_Ntab <- melt (Kegg_Ntab)

func <- function(x)(c(n = length(x),mean=mean(x,na.rm = T)))
Kegg_Ntab <- aggregate(value ~ KEGG_Name + variable, data=Kegg_Ntab,func)

ggplot(Kegg_Ntab, aes(x = variable, y = KEGG_Name, size = value[,"mean"], color=variable)) + geom_point()




##### lefse #####
library(microeco)
library(dplyr)
library(tidyverse)
library(magrittr)


## Bacteria all
keep_cp_otu_B <- which(rowSums(cp_otu_B) > 0)
cp_otu_B <- cp_otu_B[keep_cp_otu_B,]

nrow(cp_otu_B)

tax_B <- tax_B[rownames(cp_otu_B),]
tax_B$otu <- rownames(tax_B)

B_dataset <- microtable$new(sample_table = design_B,
                            otu_table = cp_otu_B, 
                            tax_table = tax_B)


B_lefse <- trans_diff$new(dataset = B_dataset, 
                          method = "lefse", 
                          group = "Treatment1", 
                          alpha = 0.05, taxa_level = "otu",
                          lefse_subgroup = NULL,
                          p_adjust_method = "none")

B_lefse_plot <- B_lefse$plot_diff_bar(use_number = 1:30, 
                                      width = 0.8, heatmap_y = "Taxa",
                                      group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() + 
  ggsci::scale_fill_npg()


## Fungi all
keep_cp_otu_F <- which(rowSums(cp_otu_F) > 0)
cp_otu_F <- cp_otu_F[keep_cp_otu_F,]

nrow(cp_otu_F)

tax_F <- tax_F[rownames(cp_otu_F),]
tax_F$otu <- rownames(tax_F)

F_dataset <- microtable$new(sample_table = design_F,
                            otu_table = cp_otu_F, 
                            tax_table = tax_F)


F_lefse <- trans_diff$new(dataset = F_dataset, 
                          method = "lefse", 
                          group = "Treatment1", 
                          alpha = 0.05, taxa_level = "otu",
                          lefse_subgroup = NULL,
                          p_adjust_method = "none")

F_lefse_plot <- F_lefse$plot_diff_bar(use_number = 1:30, 
                                      width = 0.8, heatmap_y = "Taxa",
                                      group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() + 
  ggsci::scale_fill_npg()


## Protist all
keep_cp_otu_P <- which(rowSums(cp_otu_P) > 0)
cp_otu_P <- cp_otu_P[keep_cp_otu_P,]

nrow(cp_otu_P)

tax_P <- tax_P[rownames(cp_otu_P),]
tax_P$otu <- rownames(tax_P)

design_P <- droplevels(design_P[soilsamples,])

P_dataset <- microtable$new(sample_table = design_P,
                            otu_table = cp_otu_P, 
                            tax_table = tax_P)


P_lefse <- trans_diff$new(dataset = P_dataset, 
                          method = "lefse", 
                          group = "Treatment1", 
                          alpha = 0.05, taxa_level = "otu",
                          lefse_subgroup = NULL,
                          p_adjust_method = "none")

P_lefse_plot <- P_lefse$plot_diff_bar(use_number = 1:30, 
                                      width = 0.8, heatmap_y = "Taxa",
                                      group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() + 
  ggsci::scale_fill_npg()


grid.newpage()
grid.arrange(B_lefse_plot,F_lefse_plot,P_lefse_plot,nrow = 1)



## cladogram
# bacteria
B_cladogram_ASV <- c(B_lefse[["plot_diff_bar_taxa"]])
cp_otu_B_cladogram <- cp_otu_B[B_cladogram_ASV,]


PE <- rownames(design_B[design_B$Treatment1=="PE",])
DE <- rownames(design_B[design_B$Treatment1=="DE",])
CK <- rownames(design_B[design_B$Treatment1=="CK",])


B_cladogram_ASV_tab <-as.matrix(cbind(`CK`=apply(cp_otu_B_cladogram[,CK],1,sum),
                                      `PE`=apply(cp_otu_B_cladogram[,PE],1,sum),
                                      `DE`=apply(cp_otu_B_cladogram[,DE],1,sum)))

B_cladogram_ASV_tax <- tax_B[B_cladogram_ASV,]

## fungi
F_cladogram_ASV <- c(F_lefse[["plot_diff_bar_taxa"]])
cp_otu_F_cladogram <- cp_otu_F[F_cladogram_ASV,]


PE <- rownames(design_F[design_F$Treatment1=="PE",])
DE <- rownames(design_F[design_F$Treatment1=="DE",])
CK <- rownames(design_F[design_F$Treatment1=="CK",])


F_cladogram_ASV_tab <-as.matrix(cbind(`CK`=apply(cp_otu_F_cladogram[,CK],1,sum),
                                      `PE`=apply(cp_otu_F_cladogram[,PE],1,sum),
                                      `DE`=apply(cp_otu_F_cladogram[,DE],1,sum)))

F_cladogram_ASV_tax <- tax_F[F_cladogram_ASV,]


## protist
P_cladogram_ASV <- c(P_lefse[["plot_diff_bar_taxa"]])
cp_otu_P_cladogram <- cp_otu_P[P_cladogram_ASV,]


PE <- rownames(design_P[design_P$Treatment1=="PE",])
DE <- rownames(design_P[design_P$Treatment1=="DE",])
CK <- rownames(design_P[design_P$Treatment1=="CK",])


P_cladogram_ASV_tab <-as.matrix(cbind(`CK`=apply(cp_otu_P_cladogram[,CK],1,sum),
                                      `PE`=apply(cp_otu_P_cladogram[,PE],1,sum),
                                      `DE`=apply(cp_otu_P_cladogram[,DE],1,sum)))

P_cladogram_ASV_tax <- tax_P[P_cladogram_ASV,]



## Bacteria soil
cp_otu_B_soil <- cp_otu_B[,soilsamples]

keep_cp_otu_B_soil <- which(rowSums(cp_otu_B_soil) > 0)
cp_otu_B_soil <- cp_otu_B_soil[keep_cp_otu_B_soil,]

nrow(cp_otu_B_soil)

tax_B_soil <- tax_B[rownames(cp_otu_B_soil),]
tax_B_soil$otu <- rownames(tax_B_soil)

design_B_soil <- droplevels(design_B[soilsamples,])

B_soil_dataset <- microtable$new(sample_table = design_B_soil,
                                 otu_table = cp_otu_B_soil, 
                                 tax_table = tax_B_soil)


B_soil_lefse <- trans_diff$new(dataset = B_soil_dataset, 
                               method = "lefse", 
                               group = "Treatment1", 
                               alpha = 0.05, taxa_level = "otu",
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

B_soil_lefse_plot <- B_soil_lefse$plot_diff_bar(use_number = 1:30, 
                                                width = 0.8, heatmap_y = "Taxa",
                                                group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() + 
  ggsci::scale_fill_npg()



## Bacteria root
cp_otu_B_root <- cp_otu_B[,rootsamples]

keep_cp_otu_B_root <- which(rowSums(cp_otu_B_root) > 0)
cp_otu_B_root <- cp_otu_B_root[keep_cp_otu_B_root,]

nrow(cp_otu_B_root)

tax_B_root <- tax_B[rownames(cp_otu_B_root),]
tax_B_root$otu <- rownames(tax_B_root)

design_B_root <- droplevels(design_B[rootsamples,])

B_root_dataset <- microtable$new(sample_table = design_B_root,
                                 otu_table = cp_otu_B_root, 
                                 tax_table = tax_B_root)


B_root_lefse <- trans_diff$new(dataset = B_root_dataset, 
                               method = "lefse", 
                               group = "Treatment1", taxa_level = "otu",
                               alpha = 0.01, 
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

B_root_lefse_plot <- B_root_lefse$plot_diff_bar(use_number = 1:30, 
                                                width = 0.8, 
                                                group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Bacteria end
cp_otu_B_end <- cp_otu_B[,endsamples]

keep_cp_otu_B_end <- which(rowSums(cp_otu_B_end) > 0)
cp_otu_B_end <- cp_otu_B_end[keep_cp_otu_B_end,]

nrow(cp_otu_B_end)

tax_B_end <- tax_B[rownames(cp_otu_B_end),]
tax_B_end$otu <- rownames(tax_B_end)

design_B_end <- droplevels(design_B[endsamples,])

B_end_dataset <- microtable$new(sample_table = design_B_end,
                                otu_table = cp_otu_B_end, 
                                tax_table = tax_B_end)


B_end_lefse <- trans_diff$new(dataset = B_end_dataset, 
                              method = "lefse", 
                              group = "Treatment1", taxa_level = "otu",
                              alpha = 0.05, 
                              lefse_subgroup = NULL,
                              p_adjust_method = "none")

B_end_lefse_plot <- B_end_lefse$plot_diff_bar(use_number = 1:30, 
                                              width = 0.8, 
                                              group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Bacteria film
cp_otu_B_film <- cp_otu_B[,particlesamples]

keep_cp_otu_B_film <- which(rowSums(cp_otu_B_film) > 0)
cp_otu_B_film <- cp_otu_B_film[keep_cp_otu_B_film,]

nrow(cp_otu_B_film)

tax_B_film <- tax_B[rownames(cp_otu_B_film),]
tax_B_film$otu <- rownames(tax_B_film)

design_B_film <- droplevels(design_B[particlesamples,])

B_film_dataset <- microtable$new(sample_table = design_B_film,
                                 otu_table = cp_otu_B_film, 
                                 tax_table = tax_B_film)


B_film_lefse <- trans_diff$new(dataset = B_film_dataset, 
                               method = "lefse", 
                               group = "Treatment1", 
                               alpha = 0.05, taxa_level = "otu",
                               lefse_subgroup = NULL)

B_film_lefse_plot <- B_film_lefse$plot_diff_bar(use_number = 1:30, 
                                                width = 0.8, 
                                                group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()


## Fungi soil
cp_otu_F_soil <- cp_otu_F[,soilsamples]

keep_cp_otu_F_soil <- which(rowSums(cp_otu_F_soil) > 0)
cp_otu_F_soil <- cp_otu_F_soil[keep_cp_otu_F_soil,]

nrow(cp_otu_F_soil)

tax_F_soil <- tax_F[rownames(cp_otu_F_soil),]
tax_F_soil$otu <- rownames(tax_F_soil)

design_F_soil <- droplevels(design_F[soilsamples,])

F_soil_dataset <- microtable$new(sample_table = design_F_soil,
                                 otu_table = cp_otu_F_soil, 
                                 tax_table = tax_F_soil)


F_soil_lefse <- trans_diff$new(dataset = F_soil_dataset, 
                               method = "lefse", 
                               group = "Treatment1", taxa_level = "otu",
                               alpha = 0.05, 
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

F_soil_lefse_plot <- F_soil_lefse$plot_diff_bar(use_number = 1:30, 
                                                width = 0.8, 
                                                group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Bacteria root
cp_otu_F_root <- cp_otu_F[,rootsamples]

keep_cp_otu_F_root <- which(rowSums(cp_otu_F_root) > 0)
cp_otu_F_root <- cp_otu_F_root[keep_cp_otu_F_root,]

nrow(cp_otu_F_root)

tax_F_root <- tax_F[rownames(cp_otu_F_root),]
tax_F_root$otu <- rownames(tax_F_root)

design_F_root <- droplevels(design_F[rootsamples,])

F_root_dataset <- microtable$new(sample_table = design_F_root,
                                 otu_table = cp_otu_F_root, 
                                 tax_table = tax_F_root)


F_root_lefse <- trans_diff$new(dataset = F_root_dataset, 
                               method = "lefse", 
                               group = "Treatment1", 
                               alpha = 0.01, taxa_level = "otu",
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

F_root_lefse_plot <- F_root_lefse$plot_diff_bar(use_number = 1:30, 
                                                width = 0.8, 
                                                group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Bacteria end
cp_otu_F_end <- cp_otu_F[,endsamples]

keep_cp_otu_F_end <- which(rowSums(cp_otu_F_end) > 0)
cp_otu_F_end <- cp_otu_F_end[keep_cp_otu_F_end,]

nrow(cp_otu_F_end)

tax_F_end <- tax_F[rownames(cp_otu_F_end),]
tax_F_end$otu <- rownames(tax_F_end)

design_F_end <- droplevels(design_F[endsamples,])

F_end_dataset <- microtable$new(sample_table = design_F_end,
                                otu_table = cp_otu_F_end, 
                                tax_table = tax_F_end)


F_end_lefse <- trans_diff$new(dataset = F_end_dataset, 
                              method = "lefse", 
                              group = "Treatment1", taxa_level = "otu",
                              alpha = 0.05, 
                              lefse_subgroup = NULL,
                              p_adjust_method = "none")

F_end_lefse_plot <- F_end_lefse$plot_diff_bar(use_number = 1:30, 
                                              width = 0.8, 
                                              group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Bacteria film
cp_otu_F_film <- cp_otu_F[,particlesamples]

keep_cp_otu_F_film <- which(rowSums(cp_otu_F_film) > 0)
cp_otu_F_film <- cp_otu_F_film[keep_cp_otu_F_film,]

nrow(cp_otu_F_film)

tax_F_film <- tax_F[rownames(cp_otu_F_film),]
tax_F_film$otu <- rownames(tax_F_film)

design_F_film <- droplevels(design_F[particlesamples,])

F_film_dataset <- microtable$new(sample_table = design_F_film,
                                 otu_table = cp_otu_F_film, 
                                 tax_table = tax_F_film)


F_film_lefse <- trans_diff$new(dataset = F_film_dataset, 
                               method = "lefse", 
                               group = "Treatment1", taxa_level = "otu",
                               alpha = 0.05, 
                               lefse_subgroup = NULL)

F_film_lefse_plot <- F_film_lefse$plot_diff_bar(use_number = 1:30, 
                                                width = 0.8, 
                                                group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Protist soil
cp_otu_P_soil <- cp_otu_P[,soilsamples]

keep_cp_otu_P_soil <- which(rowSums(cp_otu_P_soil) > 0)
cp_otu_P_soil <- cp_otu_P_soil[keep_cp_otu_P_soil,]

nrow(cp_otu_P_soil)

tax_P_soil <- tax_P[rownames(cp_otu_P_soil),]
tax_P_soil$otu <- rownames(tax_P_soil)

design_P_soil <- droplevels(design_P[soilsamples,])

P_soil_dataset <- microtable$new(sample_table = design_P_soil,
                                 otu_table = cp_otu_P_soil, 
                                 tax_table = tax_P_soil)


P_soil_lefse <- trans_diff$new(dataset = P_soil_dataset, 
                               method = "lefse", 
                               group = "Treatment1", 
                               alpha = 0.05, taxa_level = "otu",
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

P_soil_lefse_plot <- P_soil_lefse$plot_diff_bar(use_number = 1:30, 
                                                width = 0.8, 
                                                group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Protist root
cp_otu_P_root <- cp_otu_P[,rootsamples]

keep_cp_otu_P_root <- which(rowSums(cp_otu_P_root) > 0)
cp_otu_P_root <- cp_otu_P_root[keep_cp_otu_P_root,]

nrow(cp_otu_P_root)

tax_P_root <- tax_P[rownames(cp_otu_P_root),]
tax_P_root$otu <- rownames(tax_P_root)

design_P_root <- droplevels(design_P[rootsamples,])

P_root_dataset <- microtable$new(sample_table = design_P_root,
                                 otu_table = cp_otu_P_root, 
                                 tax_table = tax_P_root)


P_root_lefse <- trans_diff$new(dataset = P_root_dataset, 
                               method = "lefse", 
                               group = "Treatment1", taxa_level = "otu",
                               alpha = 0.01, 
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

P_root_lefse_plot <- P_root_lefse$plot_diff_bar(use_number = 1:30, 
                                                width = 0.8, 
                                                group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Protist end
cp_otu_P_end <- cp_otu_P[,pendsamples]

keep_cp_otu_P_end <- which(rowSums(cp_otu_P_end) > 0)
cp_otu_P_end <- cp_otu_P_end[keep_cp_otu_P_end,]

nrow(cp_otu_P_end)

tax_P_end <- tax_P[rownames(cp_otu_P_end),]
tax_P_end$otu <- rownames(tax_P_end)

design_P_end <- droplevels(design_P[pendsamples,])

P_end_dataset <- microtable$new(sample_table = design_P_end,
                                otu_table = cp_otu_P_end, 
                                tax_table = tax_P_end)


P_end_lefse <- trans_diff$new(dataset = P_end_dataset, 
                              method = "lefse", 
                              group = "Treatment1", taxa_level = "otu",
                              alpha = 0.05, 
                              lefse_subgroup = NULL,
                              p_adjust_method = "none")

P_end_lefse_plot <- P_end_lefse$plot_diff_bar(use_number = 1:30, 
                                              width = 0.8, 
                                              group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Protist film
cp_otu_P_film <- cp_otu_P[,particlesamples]

keep_cp_otu_P_film <- which(rowSums(cp_otu_P_film) > 0)
cp_otu_P_film <- cp_otu_P_film[keep_cp_otu_P_film,]

nrow(cp_otu_P_film)

tax_P_film <- tax_P[rownames(cp_otu_P_film),]
tax_P_film$otu <- rownames(tax_P_film)

design_P_film <- droplevels(design_P[particlesamples,])

P_film_dataset <- microtable$new(sample_table = design_P_film,
                                 otu_table = cp_otu_P_film, 
                                 tax_table = tax_P_film)


P_film_lefse <- trans_diff$new(dataset = P_film_dataset, 
                               method = "lefse", taxa_level = "otu",
                               group = "Treatment1", 
                               alpha = 0.05, 
                               lefse_subgroup = NULL)

P_film_lefse_plot <- P_film_lefse$plot_diff_bar(use_number = 1:30, 
                                                width = 0.8, 
                                                group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()


grid.newpage()
grid.arrange(B_soil_lefse_plot ,B_root_lefse_plot,B_end_lefse_plot,B_film_lefse_plot,ncol = 4)

grid.newpage()
grid.arrange(F_soil_lefse_plot ,F_root_lefse_plot,F_end_lefse_plot,F_film_lefse_plot,ncol = 4)

grid.newpage()
grid.arrange(P_soil_lefse_plot ,P_root_lefse_plot,P_end_lefse_plot,P_film_lefse_plot,ncol = 4)



## cladogram
# bacteria
B_cladogram_ASV <- c(B_soil_lefse[["plot_diff_bar_taxa"]],B_root_lefse[["plot_diff_bar_taxa"]],B_end_lefse[["plot_diff_bar_taxa"]],B_film_lefse[["plot_diff_bar_taxa"]])
B_cladogram_ASV <- unique(B_cladogram_ASV)

cp_otu_B_cladogram <- cp_otu_B[B_cladogram_ASV,]


PE <- rownames(design_B[design_B$Treatment1=="PE",])
DE <- rownames(design_B[design_B$Treatment1=="DE",])
CK <- rownames(design_B[design_B$Treatment1=="CK",])


B_cladogram_ASV_tab <-as.matrix(cbind(`CK`=apply(cp_otu_B_cladogram[,CK],1,sum),
                                      `PE`=apply(cp_otu_B_cladogram[,PE],1,sum),
                                      `DE`=apply(cp_otu_B_cladogram[,DE],1,sum)))

B_cladogram_ASV_tax <- tax_B[B_cladogram_ASV,]

write.table(B_cladogram_ASV_tax,paste0("B_cladogram_ASV_tax.txt"),sep="\t",quote=F)
write.table(B_cladogram_ASV_tab,paste0("B_cladogram_ASV.txt"),sep="\t",quote=F)


# fungi
F_cladogram_ASV <- c(F_soil_lefse[["plot_diff_bar_taxa"]],F_root_lefse[["plot_diff_bar_taxa"]],F_end_lefse[["plot_diff_bar_taxa"]],F_film_lefse[["plot_diff_bar_taxa"]])
F_cladogram_ASV <- unique(F_cladogram_ASV)

cp_otu_F_cladogram <- cp_otu_F[F_cladogram_ASV,]


PE <- rownames(design_F[design_F$Treatment1=="PE",])
DE <- rownames(design_F[design_F$Treatment1=="DE",])
CK <- rownames(design_F[design_F$Treatment1=="CK",])


F_cladogram_ASV_tab <-as.matrix(cbind(`CK`=apply(cp_otu_F_cladogram[,CK],1,sum),
                                      `PE`=apply(cp_otu_F_cladogram[,PE],1,sum),
                                      `DE`=apply(cp_otu_F_cladogram[,DE],1,sum)))

F_cladogram_ASV_tax <- tax_F[F_cladogram_ASV,]

write.table(F_cladogram_ASV_tax,paste0("F_cladogram_ASV_tax.txt"),sep="\t",quote=F)
write.table(F_cladogram_ASV_tab,paste0("F_cladogram_ASV.txt"),sep="\t",quote=F)


# protist
P_cladogram_ASV <- c(P_soil_lefse[["plot_diff_bar_taxa"]],P_root_lefse[["plot_diff_bar_taxa"]],P_end_lefse[["plot_diff_bar_taxa"]],P_film_lefse[["plot_diff_bar_taxa"]])
P_cladogram_ASV <- unique(P_cladogram_ASV)

cp_otu_P_cladogram <- cp_otu_P[P_cladogram_ASV,]


PE <- rownames(design_P[design_P$Treatment1=="PE",])
DE <- rownames(design_P[design_P$Treatment1=="DE",])
CK <- rownames(design_P[design_P$Treatment1=="CK",])


P_cladogram_ASV_tab <-as.matrix(cbind(`CK`=apply(cp_otu_P_cladogram[,CK],1,sum),
                                      `PE`=apply(cp_otu_P_cladogram[,PE],1,sum),
                                      `DE`=apply(cp_otu_P_cladogram[,DE],1,sum)))

P_cladogram_ASV_tax <- tax_P[P_cladogram_ASV,]

write.table(P_cladogram_ASV_tax,paste0("P_cladogram_ASV_tax.txt"),sep="\t",quote=F)
write.table(P_cladogram_ASV_tab,paste0("P_cladogram_ASV.txt"),sep="\t",quote=F)




##### Metagenom C/N/S/P #####
## C gene
Kegg_C <- read.table("KO_C.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")

# Kegg_C <- Kegg_C*1000000

Kegg_C_tax <- read.table("C_tax_gene.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")


## C fixcation
KO_C <- names(sort(table(Kegg_C_tax[,"KO"]), decr=T))
length(KO_C)
sort(table(Kegg_C_tax[,"KO"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_C)
for (i in KO_C){
  x <- array(colSums(Kegg_C[rownames(Kegg_C_tax)[which(Kegg_C_tax$KO == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(KO_C)
colnames(y) <- paste(colnames(Kegg_C))
KO_C_sum <- y

a <- c(	"K01602","K03841","K01807","K09709","K00242","K01007","K01595","K01678","K01676","K01681","K01682")

Kegg_C_fix <- KO_C_sum[a,]

b <- c(	"K00193","K00156","K00134","K00873","K00927","K01689","K01834","K15633","K15634","K15635","K00150","K00844","K00845","K00850","K00886","K01623","K01624","K01803","K01810","K06859","K11645","K13810","K15916","K16306","K16370")
Kegg_C_de <- KO_C_sum[b,]


##barplot 
metadata1$Treatment1 <- factor(metadata1$Treatment1,c("CK","PE","DE"))

CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])

PHYLUM_mat_C_fix <-as.matrix(cbind(`CK`=apply(Kegg_C_fix[,CK],1,sum),
                                   `PE`=apply(Kegg_C_fix[,PE],1,sum),
                                   `DE`=apply(Kegg_C_fix[,DE],1,sum)))

PHYLUM_mat_C_fix <- as.data.frame(t(Kegg_C_fix))

bargraph.CI(metadata1$Treatment1, rowSums(PHYLUM_mat_C_fix)*100, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)

phylum_bar_C <- barplot(as.matrix(PHYLUM_mat_C_fix),  xaxt="n", border=NA, las=2)


PHYLUM_mat_C_bar <- as.matrix(t(PHYLUM_mat_C_fix))


barplot(PHYLUM_mat_C_bar, col=colors()[c("#2c7fbe","#8bc224","#f7d368")] ,  border="white", font.axis=2, beside=T, font.lab=2)

par(mfrow=c(2,4), mar=c(0.5,3.5,2,0))

CS_cols <- c("#2c7fbe","#8bc224","#f7d368")

bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_fix$K03841)*100, 
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_fix$K01681)*100, 
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_fix$K01807)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_fix$K09709)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_fix$K00242)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_fix$K01007)*100, 
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_fix$K01595)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_fix$K01676)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)



PHYLUM_mat_C_de <-as.matrix(cbind(`CK`=apply(Kegg_C_de[,CK],1,sum),
                                  `PE`=apply(Kegg_C_de[,PE],1,sum),
                                  `DE`=apply(Kegg_C_de[,DE],1,sum)))

PHYLUM_mat_C_de <- as.data.frame(t(Kegg_C_de ))

bargraph.CI(metadata1$Treatment1, rowSums(PHYLUM_mat_C_de)*100, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)

phylum_bar_C <- barplot(as.matrix(PHYLUM_mat_C_de),  xaxt="n", border=NA, las=2)


PHYLUM_mat_C_bar <- as.matrix(t(Kegg_C_de))


CS_cols <- c("#2c7fbe","#8bc224","#f7d368")

bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_de$K00845)*100, 
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_de$K00156)*100, 
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)



## collaps OTU colors to prepare Phylum level colors
label_cols_C <- Kegg_C_tax[, c("Kingdom", "cols") ]
library(plyr)
PHYLA_label_cols_C <- ddply(label_cols_C, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_C) <- PHYLA_label_cols_C[,1]
PHYLA_label_cols_C

##### Plot Supplementary Figure S3
par(oma=c(0,0,0,0), mar=c(10,4,5,5), xpd=NA)
phylum_bar_C <- barplot(as.matrix(PHYLUM_mat_C_fix),  xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_C, labels=colnames(PHYLUM_mat_C_mean_type), col.axis="black", las=2, cex.axis=2)
title(main="C function gene")


PHYLUM_mat_C_fix_ra <- PHYLUM_mat_C_fix*10000



edgeR_kegg_C <- DGEList(counts=Kegg_C, 
                        group=metadata1$Treatment1)

## Input TMM normalized counts, taxonomy, and design of bulk root bacteria community into phyloseq objects
## for further analysis
phy_kegg_C <- otu_table(Kegg_C,taxa_are_rows=T)
phy_tax_kegg_C <-tax_table(as.matrix(Kegg_C_tax))
phy_design_kegg_C <- sample_data(metadata1)
physeq_kegg_C <- phyloseq(phy_kegg_C,phy_design_kegg_C)


## Create bray-curtis dissimiliartiy matrix
all_dis_kegg_C <- vegdist(t(otu_table(physeq_kegg_C)),method="bray")


## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_kegg_C <- adonis2(all_dis_kegg_C ~Treatment1, data=metadata1, permutations=9999)

paov_all_kegg_C  


##kegg_diversity
pcoa_kegg_C <- ordinate(physeq_kegg_C,"PCoA","bray")
pcoa_all_kegg_C <- plot_ordination(physeq_kegg_C, pcoa_kegg_C, type="sites", color ="Treatment1")
pcoa_all_kegg1_C <- pcoa_all_kegg_C+
  geom_point(size=5)+
  xlab(paste("PCo 1", paste("(",round(pcoa_all_kegg_C$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_all_kegg_C$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("KEGG_C")



microeco_kegg_C<- microtable$new(sample_table = metadata1, otu_table = Kegg_C, tax_table = Kegg_C_tax)

range(microeco_kegg_C$sample_sums())
microeco_kegg_C$rarefy_samples(sample.size =min(colSums(microeco_kegg_C$otu_table)))
microeco_kegg_C$sample_sums() %>% range
print(microeco_kegg_C)
microeco_kegg_C$cal_abund()
microeco_kegg_C$cal_alphadiv()
microeco_kegg_C$cal_betadiv()
microeco_kegg_C


lefse_microeco_kegg_C<- trans_diff$new(dataset = microeco_kegg_C, method = "lefse", group = "Treatment1", alpha = 0.01, taxa_level = "Function",group_order = c("CK","PE","DE"),p_adjust_method = "none")

###write.csv(lefse_microeco_B_LD$res_abund,"diff_abund_B_otu_lefse.csv")
##write.csv(lefse_microeco_B_LD$res_diff,"diff_diff_B_otu_lefse.csv")
# we show 40 taxa with the highest LDA (log10)
g1 <- lefse_microeco_kegg_C$plot_diff_bar(use_number = 1:100,color_values = c("#D95F02","#1B9E77" ,"lightcoral"), group_order = c("CK","PE","DE"))
g2 <- lefse_microeco_kegg_C$plot_diff_abund(color_values = c("#D95F02","#1B9E77" ,"lightcoral"),group_order = c("CK","PE","DE"),select_taxa = lefse_microeco_kegg_C$plot_diff_bar_taxa,add_sig = TRUE)
g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
gridExtra::grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2, 1.7))

diff_otu_lefse_kegg_C <- lefse_microeco_kegg_C$res_diff


###class
lefse_microeco_kegg_C_Class<- trans_diff$new(dataset = microeco_kegg_C, method = "lefse", group = "Treatment1", alpha = 0.01, taxa_level = "Function",group_order = c("CK","PE","DE"),p_adjust_method = "none")

lefse_kegg_C_class <- lefse_microeco_kegg_C_Class$plot_diff_bar(use_number = 1:10,group_order = c("CK","PE","DE"),color_values = c("#D95F02","#1B9E77" ,"lightcoral"))






##barplot 
CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])


Kegg_C <- t(t(Kegg_C)/colSums(Kegg_C)) * 100


PHYLAnames_C <- names(sort(table(Kegg_C_tax[,"KO"]), decr=T))
length(Kegg_C)
sort(table(Kegg_C_tax[,"KO"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_C)
for (i in PHYLAnames_C){
  x <- array(colSums(Kegg_C[rownames(Kegg_C_tax)[which(Kegg_C_tax$KO == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}


## Create matrix
rownames(y) <- paste(PHYLAnames_C)
colnames(y) <- paste(colnames(Kegg_C))
PHYLUM_mat_C <- y
PHYLUM_mat_C[,1:5]
colSums(PHYLUM_mat_C)
PHYLUM_mat_C_mean <- sort(apply(PHYLUM_mat_C,1,mean),decr=T)
PHYLUM_mat_C <- PHYLUM_mat_C[names(PHYLUM_mat_C_mean),]


PHYLUM_mat_C_mean_type <-as.matrix(cbind(`CK`=apply(PHYLUM_mat_C[,CK],1,mean),
                                         `PE`=apply(PHYLUM_mat_C[,PE],1,mean),
                                         `DE`=apply(PHYLUM_mat_C[,DE],1,mean)))


a <- c("K00029", "K00123","K00164","K00242","K01595","K01690")

PHYLUM_mat_C_bar <- PHYLUM_mat_C[a,]
PHYLUM_mat_C_bar <- as.data.frame(t(PHYLUM_mat_C_bar))


par(mfrow=c(1,6), mar=c(0.5,3.5,2,0))

CS_cols <- c("#2c7fbe","#8bc224","#f7d368")

bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_bar$K00029),
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_bar$K00123), 
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_bar$K00164), 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_bar$K00242), 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_bar$K01595), 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_C_bar$K01690), 
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)




barplot(PHYLUM_mat_C_bar, col=colors()[c(23,89,12)] ,  border="white", font.axis=2, beside=T, font.lab=2)


b <- c("K01676", "K00126","K04073","K00850","K09709","K01682")

PHYLUM_mat_C_bar <- PHYLUM_mat_C_mean_type[b,]
PHYLUM_mat_C_bar <- as.matrix(t(PHYLUM_mat_C_bar))


barplot(PHYLUM_mat_C_bar, col=colors()[c(23,89,12)] ,  border="white", font.axis=2, beside=T, font.lab=2)



# Phyla with MEAN abundances higher than 1% relative abundances
Kegg_C_tax$cols <- Kegg_C_tax$Kingdom
Kegg_C_tax[ rownames(Kegg_C_tax)[Kegg_C_tax$Kingdom=="Bacteria" ], ]$cols <- "peachpuff3"
Kegg_C_tax[ rownames(Kegg_C_tax)[Kegg_C_tax$Kingdom=="unclassified" ], ]$cols <- "orchid3"
Kegg_C_tax[ rownames(Kegg_C_tax)[Kegg_C_tax$Kingdom=="Archaea" ], ]$cols <- "gold1"
Kegg_C_tax[ rownames(Kegg_C_tax)[Kegg_C_tax$Kingdom=="Eukaryota" ], ]$cols <- "lightsalmon4"
Kegg_C_tax[ rownames(Kegg_C_tax)[Kegg_C_tax$Kingdom=="Unclassified" ], ]$cols <- "tan1"
Kegg_C_tax[ rownames(Kegg_C_tax)[Kegg_C_tax$Kingdom=="unknown" ], ]$cols <- "palegreen4"


## collaps OTU colors to prepare Phylum level colors
label_cols_C <- Kegg_C_tax[, c("Kingdom", "cols") ]
library(plyr)
PHYLA_label_cols_C <- ddply(label_cols_C, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_C) <- PHYLA_label_cols_C[,1]
PHYLA_label_cols_C

##### Plot Supplementary Figure S3
par(oma=c(0,0,0,0), mar=c(10,4,5,5), xpd=NA)
phylum_bar_C <- barplot(as.matrix(PHYLUM_mat_C_mean_type),  xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_C, labels=colnames(PHYLUM_mat_C_mean_type), col.axis="black", las=2, cex.axis=2)
title(main="C function gene")




## N gene
Kegg_N <- read.table("KO_N.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")
Kegg_N_tax <- read.table("N_tax_gene.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")

edgeR_kegg_N <- DGEList(counts=Kegg_N, 
                        group=metadata1$Treatment1)


## N cycle
KO_N <- names(sort(table(Kegg_N_tax[,"Function"]), decr=T))
length(KO_N)
sort(table(Kegg_N_tax[,"Function"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_N)
for (i in KO_N){
  x <- array(colSums(Kegg_N[rownames(Kegg_N_tax)[which(Kegg_N_tax$Function == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(KO_N)
colnames(y) <- paste(colnames(Kegg_N))
KO_N_sum <- y

a <- c(	"K17229","K21307","K21308","K21309","K00380","K00381","K17224","K18277","K21307","K17226")

Kegg_N_fix <- KO_N_sum[a,]


##barplot 
metadata1$Treatment1 <- factor(metadata1$Treatment1,c("CK","PE","DE"))

CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])

PHYLUM_mat_N_fix <-as.matrix(cbind(`CK`=apply(Kegg_N[,CK],1,sum),
                                   `PE`=apply(Kegg_N[,PE],1,sum),
                                   `DE`=apply(Kegg_N[,DE],1,sum)))

PHYLUM_mat_N_fix <- as.data.frame(t(Kegg_N_fix))

bargraph.CI(metadata1$Treatment1, rowSums(PHYLUM_mat_N_fix)*100, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)

phylum_bar_N <- barplot(as.matrix(PHYLUM_mat_N_fix),  xaxt="n", border=NA, las=2)


PHYLUM_mat_N_bar <- as.matrix(t(PHYLUM_mat_N_fix))


barplot(PHYLUM_mat_N_bar, col=colors()[c("#2c7fbe","#8bc224","#f7d368")] ,  border="white", font.axis=2, beside=T, font.lab=2)

par(mfrow=c(2,4), mar=c(0.5,3.5,2,0))

CS_cols <- c("#2c7fbe","#f7d368","#8bc224")

bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_N_fix$K17224)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
# bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K18277)*100, 
#             las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
#             err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K21307)*100, 
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K17226)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)

bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K17229)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K21308)*100,
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K21309)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K00380)*100, 
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K00381)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)





## Input TMM normalized counts, taxonomy, and design of bulk root bacteria community into phyloseq objects
## for further analysis
phy_kegg_N <- otu_table(Kegg_N,taxa_are_rows=T)
phy_tax_kegg_N <-tax_table(as.matrix(Kegg_N_tax))
phy_design_kegg_N <- sample_data(metadata1)
physeq_kegg_N <- phyloseq(phy_kegg_N,phy_design_kegg_N)


## Create bray-curtis dissimiliartiy matrix
all_dis_kegg_N <- vegdist(t(otu_table(physeq_kegg_N)),method="bray")


## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_kegg_N <- adonis2(all_dis_kegg_N ~Treatment1, data=metadata1, permutations=9999)

paov_all_kegg_N  


##kegg_diversity
pcoa_kegg_N <- ordinate(physeq_kegg_N,"PCoA","bray")
pcoa_all_kegg_N <- plot_ordination(physeq_kegg_N, pcoa_kegg_N, type="sites", color ="Treatment1")
pcoa_all_kegg1_N <- pcoa_all_kegg_N+
  geom_point(size=5)+
  xlab(paste("PCo 1", paste("(",round(pcoa_all_kegg_N$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_all_kegg_N$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("KEGG_N")


##barplot 
PHYLAnames_N <- names(sort(table(Kegg_N_tax[,"Function"]), decr=T))
length(Kegg_N)
sort(table(Kegg_N_tax[,"Function"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_N)
for (i in PHYLAnames_N){
  x <- array(colSums(Kegg_N[rownames(Kegg_N_tax)[which(Kegg_N_tax$Function == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}


## Create matrix
rownames(y) <- paste(PHYLAnames_N)
colnames(y) <- paste(colnames(Kegg_N))
PHYLUM_mat_N <- y
PHYLUM_mat_N[,1:5]
colSums(PHYLUM_mat_N)
PHYLUM_mat_N_mean <- sort(apply(PHYLUM_mat_N,1,mean),decr=T)
PHYLUM_mat_N <- PHYLUM_mat_N[names(PHYLUM_mat_N_mean),]

PHYLUM_mat_N <- t(PHYLUM_mat_N)
PHYLUM_mat_N <- as.data.frame(PHYLUM_mat_N)
PHYLUM_mat_N$Treatment <- metadata1$Treatment1

for(i in PHYLUM_mat_N[,1:27]) {
  fit <- aov(i~Treatment,data = PHYLUM_mat_N)
  print(summary(fit))
  out <- HSD.test(fit,"Treatment")
  print(out$groups)
}







PHYLUM_mat_N_mean_type <-as.matrix(cbind(`CK`=apply(PHYLUM_mat_N[,CK],1,mean),
                                         `PE`=apply(PHYLUM_mat_N[,PE],1,mean),
                                         `DE`=apply(PHYLUM_mat_N[,DE],1,mean)))


# Phyla with MEAN abundances higher than 1% relative abundances
Kegg_N_tax$cols <- Kegg_N_tax$Function
Kegg_N_tax[ rownames(Kegg_N_tax)[Kegg_N_tax$Kingdom=="Bacteria" ], ]$cols <- "peachpuff3"
Kegg_N_tax[ rownames(Kegg_N_tax)[Kegg_N_tax$Kingdom=="unclassified" ], ]$cols <- "orchid3"
Kegg_N_tax[ rownames(Kegg_N_tax)[Kegg_N_tax$Kingdom=="Archaea" ], ]$cols <- "gold1"
Kegg_N_tax[ rownames(Kegg_N_tax)[Kegg_N_tax$Kingdom=="Eukaryota" ], ]$cols <- "lightsalmon4"
Kegg_N_tax[ rownames(Kegg_N_tax)[Kegg_N_tax$Kingdom=="Unclassified" ], ]$cols <- "tan1"


## collaps OTU colors to prepare Phylum level colors
label_cols_N <- Kegg_N_tax[, c("Kingdom", "cols") ]
library(plyr)
PHYLA_label_cols_N <- ddply(label_cols_N, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_N) <- PHYLA_label_cols_N[,1]
PHYLA_label_cols_N

##### Plot Supplementary Figure S3
par(oma=c(0,0,0,0), mar=c(10,4,5,5), xpd=NA)
phylum_bar_N <- barplot(as.matrix(PHYLUM_mat_N_mean_type), col=PHYLA_label_cols_N[rownames(PHYLUM_mat_N),]$cols, xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_N, labels=colnames(PHYLUM_mat_C_mean_type), col.axis="black", las=2, cex.axis=2)
title(main="N function gene")



## S gene
Kegg_S <- read.table("KO_S.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")
Kegg_S_tax <- read.table("S_tax_gene.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")

edgeR_kegg_S <- DGEList(counts=Kegg_S, 
                        group=metadata1$Treatment1)


## S cycle
KO_S <- names(sort(table(Kegg_S_tax[,"KO"]), decr=T))
length(KO_S)
sort(table(Kegg_S_tax[,"KO"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_S)
for (i in KO_S){
  x <- array(colSums(Kegg_S[rownames(Kegg_S_tax)[which(Kegg_S_tax$KO == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(KO_S)
colnames(y) <- paste(colnames(Kegg_S))
KO_S_sum <- y

a <- c(	"K17229","K21307","K21308","K21309","K00380","K00381","K17224","K18277","K21307","K17226")

Kegg_S_fix <- KO_S_sum[a,]


##barplot 
metadata1$Treatment1 <- factor(metadata1$Treatment1,c("CK","PE","DE"))

CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])

PHYLUM_mat_S_fix <-as.matrix(cbind(`CK`=apply(Kegg_S_fix[,CK],1,sum),
                                   `PE`=apply(Kegg_S_fix[,PE],1,sum),
                                   `DE`=apply(Kegg_S_fix[,DE],1,sum)))

PHYLUM_mat_S_fix <- as.data.frame(t(Kegg_S_fix))

bargraph.CI(metadata1$Treatment1, rowSums(PHYLUM_mat_S_fix)*100, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)

phylum_bar_S <- barplot(as.matrix(PHYLUM_mat_S_fix),  xaxt="n", border=NA, las=2)


PHYLUM_mat_S_bar <- as.matrix(t(PHYLUM_mat_S_fix))


barplot(PHYLUM_mat_S_bar, col=colors()[c("#2c7fbe","#8bc224","#f7d368")] ,  border="white", font.axis=2, beside=T, font.lab=2)

par(mfrow=c(2,4), mar=c(0.5,3.5,2,0))

CS_cols <- c("#2c7fbe","#f7d368","#8bc224")

bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K17224)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
# bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K18277)*100, 
#             las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
#             err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K21307)*100, 
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K17226)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)

bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K17229)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K21308)*100,
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K21309)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K00380)*100, 
            las=2, ylab="", cex.lab=.5, col=CS_cols,cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_S_fix$K00381)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)



## Input TMM normalized counts, taxonomy, and design of bulk root bacteria community into phyloseq objects
## for further analysis
phy_kegg_S <- otu_table(Kegg_S,taxa_are_rows=T)
phy_tax_kegg_S <-tax_table(as.matrix(Kegg_S_tax))
phy_design_kegg_S <- sample_data(metadata1)
physeq_kegg_S <- phyloseq(phy_kegg_S,phy_design_kegg_S)


## Create bray-curtis dissimiliartiy matrix
all_dis_kegg_S <- vegdist(t(otu_table(physeq_kegg_S)),method="bray")


## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_kegg_S <- adonis2(all_dis_kegg_S ~Treatment1, data=metadata1, permutations=9999)

paov_all_kegg_S  


##kegg_diversity
pcoa_kegg_S <- ordinate(physeq_kegg_S,"PCoA","bray")
pcoa_all_kegg_S <- plot_ordination(physeq_kegg_S, pcoa_kegg_S, type="sites", color ="Treatment1")
pcoa_all_kegg1_S <- pcoa_all_kegg_S+
  geom_point(size=5)+
  xlab(paste("PCo 1", paste("(",round(pcoa_all_kegg_S$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_all_kegg_S$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("KEGG_S")


##barplot 
PHYLAnames_S <- names(sort(table(Kegg_S_tax[,"Function"]), decr=T))
length(Kegg_S)
sort(table(Kegg_S_tax[,"Function"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_S)
for (i in PHYLAnames_S){
  x <- array(colSums(Kegg_S[rownames(Kegg_S_tax)[which(Kegg_S_tax$Function == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}


## Create matrix
rownames(y) <- paste(PHYLAnames_S)
colnames(y) <- paste(colnames(Kegg_S))
PHYLUM_mat_S <- y
PHYLUM_mat_S[,1:5]
colSums(PHYLUM_mat_S)
PHYLUM_mat_S_mean <- sort(apply(PHYLUM_mat_S,1,mean),decr=T)
PHYLUM_mat_S <- PHYLUM_mat_S[names(PHYLUM_mat_S_mean),]


PHYLUM_mat_S <- t(PHYLUM_mat_S)
PHYLUM_mat_S <- as.data.frame(PHYLUM_mat_S)
PHYLUM_mat_S$Treatment <- metadata1$Treatment1

for(i in PHYLUM_mat_S[,1:37]) {
  fit <- aov(i~Treatment,data = PHYLUM_mat_S)
  a <- summary(fit)
  print(a[[1]][["Pr(>F)"]])
}



newdata=c()
for (i in 1:dim(newid)[1]){
  family <- newid[i,]
  if ((family[1] < family[2]) | (family[1] < family[3])){
    newdata=rbind(newdata,family)
  }
}
show(newdata)




PHYLUM_mat_S_mean_type <-as.matrix(cbind(`CK`=apply(PHYLUM_mat_S[,CK],1,mean),
                                         `PE`=apply(PHYLUM_mat_S[,PE],1,mean),
                                         `DE`=apply(PHYLUM_mat_S[,DE],1,mean)))



# Phyla with MEAN abundances higher than 1% relative abundances
Kegg_S_tax$cols <- Kegg_S_tax$Kingdom
Kegg_S_tax[ rownames(Kegg_S_tax)[Kegg_S_tax$Kingdom=="Bacteria" ], ]$cols <- "peachpuff3"
Kegg_S_tax[ rownames(Kegg_S_tax)[Kegg_S_tax$Kingdom=="unclassified" ], ]$cols <- "orchid3"
Kegg_S_tax[ rownames(Kegg_S_tax)[Kegg_S_tax$Kingdom=="Archaea" ], ]$cols <- "gold1"
Kegg_S_tax[ rownames(Kegg_S_tax)[Kegg_S_tax$Kingdom=="Eukaryota" ], ]$cols <- "lightsalmon4"
Kegg_S_tax[ rownames(Kegg_S_tax)[Kegg_S_tax$Kingdom=="unknown" ], ]$cols <- "tan1"


## collaps OTU colors to prepare Phylum level colors
label_cols_S <- Kegg_S_tax[, c("Kingdom", "cols") ]
library(plyr)
PHYLA_label_cols_S <- ddply(label_cols_S, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_S) <- PHYLA_label_cols_S[,1]
PHYLA_label_cols_S

##### Plot Supplementary Figure S3
par(oma=c(0,0,0,0), mar=c(10,4,5,5), xpd=NA)
phylum_bar_S <- barplot(as.matrix(PHYLUM_mat_S_mean_type), col=PHYLA_label_cols_S[rownames(PHYLUM_mat_S),]$cols, xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_S, labels=colnames(PHYLUM_mat_S_mean_type), col.axis="black", las=2, cex.axis=2)
title(main="S function gene")



## P gene
Kegg_P <- read.table("KO_P.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")
Kegg_P_tax <- read.table("P_tax_gene.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")

edgeR_kegg_P <- DGEList(counts=Kegg_P, 
                        group=metadata1$Treatment1)



## P cycle
KO_P <- names(sort(table(Kegg_P_tax[,"KO"]), decr=T))
length(KO_P)
sort(table(Kegg_P_tax[,"KO"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_P)
for (i in KO_P){
  x <- array(colSums(Kegg_P[rownames(Kegg_P_tax)[which(Kegg_P_tax$KO == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(KO_P)
colnames(y) <- paste(colnames(Kegg_P))
KO_P_sum <- y

a <- c(	"K11912","K02457")

Kegg_P_fix <- KO_P_sum[a,]

##barplot 
metadata1$Treatment1 <- factor(metadata1$Treatment1,c("CK","PE","DE"))

CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])

PHYLUM_mat_P_fix <-as.matrix(cbind(`CK`=apply(Kegg_P_fix[,CK],1,sum),
                                   `PE`=apply(Kegg_P_fix[,PE],1,sum),
                                   `DE`=apply(Kegg_P_fix[,DE],1,sum)))

PHYLUM_mat_P_fix <- as.data.frame(t(Kegg_P_fix))

bargraph.CI(metadata1$Treatment1, rowSums(PHYLUM_mat_P_fix)*100, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)

phylum_bar_P <- barplot(as.matrix(PHYLUM_mat_P_fix),  xaxt="n", border=NA, las=2)


PHYLUM_mat_P_bar <- as.matrix(t(PHYLUM_mat_P_fix))


barplot(PHYLUM_mat_P_bar, col=colors()[c("#2c7fbe","#8bc224","#f7d368")] ,  border="white", font.axis=2, beside=T, font.lab=2)

par(mfrow=c(1,2), mar=c(0.5,3.5,2,0))

CS_cols <- c("#2c7fbe","#8bc224","#f7d368")

bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_P_fix$K11912)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)
bargraph.CI(metadata1$Treatment1, (PHYLUM_mat_P_fix$K02457)*100, 
            las=2, ylab="", cex.lab=.5,col=CS_cols, cex.axis=.7, cex.names=.7,
            err.width=.025,  border=F)





## Input TMM normalized counts, taxonomy, and design of bulk root bacteria community into phyloseq objects
## for further analysis
phy_kegg_P <- otu_table(Kegg_P,taxa_are_rows=T)
phy_tax_kegg_P <-tax_table(as.matrix(Kegg_P_tax))
phy_design_kegg_P <- sample_data(metadata1)
physeq_kegg_P <- phyloseq(phy_kegg_P,phy_design_kegg_P)


## Create bray-curtis dissimiliartiy matrix
all_dis_kegg_P <- vegdist(t(otu_table(physeq_kegg_P)),method="bray")


## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_kegg_P <- adonis2(all_dis_kegg_P ~Treatment1, data=metadata1, permutations=9999)

paov_all_kegg_P  


##kegg_diversity
pcoa_kegg_P <- ordinate(physeq_kegg_P,"PCoA","bray")
pcoa_all_kegg_P <- plot_ordination(physeq_kegg_P, pcoa_kegg_P, type="sites", color ="Treatment1")
pcoa_all_kegg1_P <- pcoa_all_kegg_P+
  geom_point(size=5)+
  xlab(paste("PCo 1", paste("(",round(pcoa_all_kegg_C$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_all_kegg_C$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("KEGG_P")


##barplot 
PHYLAnames_P <- names(sort(table(Kegg_P_tax[,"Kingdom"]), decr=T))
length(Kegg_P)
sort(table(Kegg_P_tax[,"Kingdom"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_P)
for (i in PHYLAnames_P){
  x <- array(colSums(Kegg_P[rownames(Kegg_P_tax)[which(Kegg_P_tax$Kingdom == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}


## Create matrix
rownames(y) <- paste(PHYLAnames_P)
colnames(y) <- paste(colnames(Kegg_P))
PHYLUM_mat_P <- y
PHYLUM_mat_P[,1:5]
colSums(PHYLUM_mat_P)
PHYLUM_mat_P_mean <- sort(apply(PHYLUM_mat_P,1,mean),decr=T)
PHYLUM_mat_P <- PHYLUM_mat_P[names(PHYLUM_mat_P_mean),]

PHYLUM_mat_P_mean_type <-as.matrix(cbind(`CK`=apply(PHYLUM_mat_P[,CK],1,mean),
                                         `PE`=apply(PHYLUM_mat_P[,PE],1,mean),
                                         `DE`=apply(PHYLUM_mat_P[,DE],1,mean)))



# Phyla with MEAN abundances higher than 1% relative abundances
Kegg_P_tax$cols <- Kegg_P_tax$Kingdom
Kegg_P_tax[ rownames(Kegg_P_tax)[Kegg_P_tax$Kingdom=="Bacteria" ], ]$cols <- "peachpuff3"
Kegg_P_tax[ rownames(Kegg_P_tax)[Kegg_P_tax$Kingdom=="unclassified" ], ]$cols <- "orchid3"
Kegg_P_tax[ rownames(Kegg_P_tax)[Kegg_P_tax$Kingdom=="Archaea" ], ]$cols <- "gold1"
Kegg_P_tax[ rownames(Kegg_P_tax)[Kegg_P_tax$Kingdom=="Eukaryota" ], ]$cols <- "lightsalmon4"
Kegg_P_tax[ rownames(Kegg_P_tax)[Kegg_P_tax$Kingdom=="Unclassified" ], ]$cols <- "tan1"


## collaps OTU colors to prepare Phylum level colors
label_cols_P <- Kegg_P_tax[, c("Kingdom", "cols") ]
library(plyr)
PHYLA_label_cols_P <- ddply(label_cols_P, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_P) <- PHYLA_label_cols_P[,1]
PHYLA_label_cols_P

##### Plot Supplementary Figure S3
par(oma=c(0,0,0,0), mar=c(10,4,5,5), xpd=NA)
phylum_bar_P <- barplot(as.matrix(PHYLUM_mat_P_mean_type), col=PHYLA_label_cols_P[rownames(PHYLUM_mat_P),]$cols, xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_P, labels=colnames(PHYLUM_mat_P_mean_type), col.axis="black", las=2, cex.axis=2)
title(main="P function gene")


par(mfrow=c(2,2))




##### neutral model #####
library(Hmisc)
library(minpack.lm)
library(stats4)
library(spaa)


## cp_otu_B_soil_CK 
cp_otu_B_soil_CK <- cp_otu_B[,soilCKsamples]

keep_cp_otu_B_soil_CK <- which(rowSums(cp_otu_B_soil_CK ) >= 0)
cp_otu_B_soil_CK <- cp_otu_B_soil_CK[keep_cp_otu_B_soil_CK,]

nrow(cp_otu_B_soil_CK)

tax_B_soil_CK <- tax_B[rownames(cp_otu_B_soil_CK),]
design_B_soil_CK <- droplevels(design_B[soilCKsamples,])


spp<-t(cp_otu_B_soil_CK)

B_soil_CK_niche_width <- niche.width(cp_otu_B_soil_CK, method = 'levins')
B_soil_CK_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  



bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'



library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_B_soil_PE 
cp_otu_B_soil_PE <- cp_otu_B[,soilPEsamples]

keep_cp_otu_B_soil_PE <- which(rowSums(cp_otu_B_soil_PE ) >= 0)
cp_otu_B_soil_PE <- cp_otu_B_soil_PE[keep_cp_otu_B_soil_PE,]

nrow(cp_otu_B_soil_PE)

tax_B_soil_PE <- tax_B[rownames(cp_otu_B_soil_PE),]
design_B_soil_PE <- droplevels(design_B[soilPEsamples,])


spp<-t(cp_otu_B_soil_PE)

B_soil_PE_niche_width <- niche.width(cp_otu_B_soil_PE, method = 'levins')
B_soil_PE_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_B_soil_DE 
cp_otu_B_soil_DE <- cp_otu_B[,soilDEsamples]

keep_cp_otu_B_soil_DE <- which(rowSums(cp_otu_B_soil_DE ) >= 0)
cp_otu_B_soil_DE <- cp_otu_B_soil_DE[keep_cp_otu_B_soil_DE,]

nrow(cp_otu_B_soil_DE)

tax_B_soil_DE <- tax_B[rownames(cp_otu_B_soil_DE),]
design_B_soil_DE <- droplevels(design_B[soilDEsamples,])


spp<-t(cp_otu_B_soil_DE)

B_soil_DE_niche_width <- niche.width(cp_otu_B_soil_DE, method = 'levins')
B_soil_DE_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)





## cp_otu_B_root_CK 
cp_otu_B_root_CK <- cp_otu_B[,rootCKsamples]

keep_cp_otu_B_root_CK <- which(rowSums(cp_otu_B_root_CK ) >= 0)
cp_otu_B_root_CK <- cp_otu_B_root_CK[keep_cp_otu_B_root_CK,]

nrow(cp_otu_B_root_CK)

tax_B_root_CK <- tax_B[rownames(cp_otu_B_root_CK),]
design_B_root_CK <- droplevels(design_B[rootCKsamples,])


spp<-t(cp_otu_B_root_CK)

B_root_CK_niche_width <- niche.width(cp_otu_B_root_CK, method = 'levins')
B_root_CK_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_B_root_PE 
cp_otu_B_root_PE <- cp_otu_B[,rootPEsamples]

keep_cp_otu_B_root_PE <- which(rowSums(cp_otu_B_root_PE ) >= 0)
cp_otu_B_root_PE <- cp_otu_B_root_PE[keep_cp_otu_B_root_PE,]

nrow(cp_otu_B_root_PE)

tax_B_root_PE <- tax_B[rownames(cp_otu_B_root_PE),]
design_B_root_PE <- droplevels(design_B[rootPEsamples,])


spp<-t(cp_otu_B_root_PE)

B_root_PE_niche_width <- niche.width(cp_otu_B_root_PE, method = 'levins')
B_root_PE_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_B_root_DE 
cp_otu_B_root_DE <- cp_otu_B[,rootDEsamples]

keep_cp_otu_B_root_DE <- which(rowSums(cp_otu_B_root_DE ) >= 0)
cp_otu_B_root_DE <- cp_otu_B_root_DE[keep_cp_otu_B_root_DE,]

nrow(cp_otu_B_root_DE)

tax_B_root_DE <- tax_B[rownames(cp_otu_B_root_DE),]
design_B_root_DE <- droplevels(design_B[rootDEsamples,])


spp<-t(cp_otu_B_root_DE)

B_root_DE_niche_width <- niche.width(cp_otu_B_root_DE, method = 'levins')
B_root_DE_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_B_end_CK 
cp_otu_B_end_CK <- cp_otu_B[,endCKsamples]

keep_cp_otu_B_end_CK <- which(rowSums(cp_otu_B_end_CK ) >= 0)
cp_otu_B_end_CK <- cp_otu_B_end_CK[keep_cp_otu_B_end_CK,]

nrow(cp_otu_B_end_CK)

tax_B_end_CK <- tax_B[rownames(cp_otu_B_end_CK),]
design_B_end_CK <- droplevels(design_B[endCKsamples,])


spp<-t(cp_otu_B_end_CK)

B_end_CK_niche_width <- niche.width(cp_otu_B_end_CK, method = 'levins')
B_end_CK_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_B_end_PE 
cp_otu_B_end_PE <- cp_otu_B[,endPEsamples]

keep_cp_otu_B_end_PE <- which(rowSums(cp_otu_B_end_PE ) >= 0)
cp_otu_B_end_PE <- cp_otu_B_end_PE[keep_cp_otu_B_end_PE,]

nrow(cp_otu_B_end_PE)

tax_B_end_PE <- tax_B[rownames(cp_otu_B_end_PE),]
design_B_end_PE <- droplevels(design_B[endPEsamples,])


spp<-t(cp_otu_B_end_PE)

B_end_PE_niche_width <- niche.width(cp_otu_B_end_PE, method = 'levins')
B_end_PE_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_B_end_DE 
cp_otu_B_end_DE <- cp_otu_B[,endDEsamples]

keep_cp_otu_B_end_DE <- which(rowSums(cp_otu_B_end_DE ) >= 0)
cp_otu_B_end_DE <- cp_otu_B_end_DE[keep_cp_otu_B_end_DE,]

nrow(cp_otu_B_end_DE)

tax_B_end_DE <- tax_B[rownames(cp_otu_B_end_DE),]
design_B_end_DE <- droplevels(design_B[endDEsamples,])


spp<-t(cp_otu_B_end_DE)

B_end_DE_niche_width <- niche.width(cp_otu_B_end_DE, method = 'levins')
B_end_DE_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_B_film_PE 
cp_otu_B_film_PE <- cp_otu_B[,PEsamples]

keep_cp_otu_B_film_PE <- which(rowSums(cp_otu_B_film_PE ) >= 0)
cp_otu_B_film_PE <- cp_otu_B_film_PE[keep_cp_otu_B_film_PE,]

nrow(cp_otu_B_film_PE)

tax_B_film_PE <- tax_B[rownames(cp_otu_B_film_PE),]
design_B_film_PE <- droplevels(design_B[PEsamples,])


spp<-t(cp_otu_B_film_PE)

B_film_PE_niche_width <- niche.width(cp_otu_B_film_PE, method = 'levins')
B_film_PE_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_B_film_DE 
cp_otu_B_film_DE <- cp_otu_B[,DEsamples]

keep_cp_otu_B_film_DE <- which(rowSums(cp_otu_B_film_DE ) >= 0)
cp_otu_B_film_DE <- cp_otu_B_film_DE[keep_cp_otu_B_film_DE,]

nrow(cp_otu_B_film_DE)

tax_B_film_DE <- tax_B[rownames(cp_otu_B_film_DE),]
design_B_film_DE <- droplevels(design_B[DEsamples,])


spp<-t(cp_otu_B_film_DE)

B_film_DE_niche_width <- niche.width(cp_otu_B_film_DE, method = 'levins')
B_film_DE_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)



B_soil_niche_width <- data.frame(CK=t(B_soil_CK_niche_width),PE=t(B_soil_PE_niche_width),DE=t(B_soil_DE_niche_width))

CK <- B_soil_CK_niche_width
PE <- B_soil_CK_niche_width
DE <- B_soil_DE_niche_width

my_comparisons = list( c("CK", "PE"), c("CK", "DE"), c("PE", "DE") )

B_soil_niche_width <- melt(B_soil_niche_width)

B_soil_width_box <- ggplot(B_soil_niche_width, aes(x=factor(variable,levels = c("CK","PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#d2a128","#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



B_root_niche_width <- data.frame(CK=t(B_root_CK_niche_width),PE=t(B_root_PE_niche_width),DE=t(B_root_DE_niche_width))

CK <- B_soil_CK_niche_width
PE <- B_soil_CK_niche_width
DE <- B_soil_DE_niche_width

my_comparisons = list( c("CK", "PE"), c("CK", "DE"), c("PE", "DE") )

B_root_niche_width <- melt(B_root_niche_width)

B_root_width_box <- ggplot(B_root_niche_width, aes(x=factor(variable,levels = c("CK","PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#d2a128","#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



B_end_niche_width <- data.frame(CK=t(B_end_CK_niche_width),PE=t(B_end_PE_niche_width),DE=t(B_end_DE_niche_width))

CK <- B_soil_CK_niche_width
PE <- B_soil_CK_niche_width
DE <- B_soil_DE_niche_width

my_comparisons = list( c("CK", "PE"), c("CK", "DE"), c("PE", "DE") )

B_end_niche_width <- melt(B_end_niche_width)

B_end_width_box <- ggplot(B_end_niche_width, aes(x=factor(variable,levels = c("CK","PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#d2a128","#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



B_film_niche_width <- data.frame(PE=t(B_film_PE_niche_width),DE=t(B_film_DE_niche_width))

PE <- B_soil_CK_niche_width
DE <- B_soil_DE_niche_width

my_comparisons = list( c("PE", "DE") )

B_film_niche_width <- melt(B_film_niche_width)

B_film_width_box <- ggplot(B_film_niche_width, aes(x=factor(variable,levels = c("PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



grid.newpage()
grid.arrange(B_soil_width_box,B_root_width_box,B_end_width_box,B_film_width_box,ncol = 4)



## cp_otu_F_soil_CK 
cp_otu_F_soil_CK <- cp_otu_F[,soilCKsamples]

keep_cp_otu_F_soil_CK <- which(rowSums(cp_otu_F_soil_CK ) >= 0)
cp_otu_F_soil_CK <- cp_otu_F_soil_CK[keep_cp_otu_F_soil_CK,]

nrow(cp_otu_F_soil_CK)

tax_F_soil_CK <- tax_F[rownames(cp_otu_F_soil_CK),]
design_F_soil_CK <- droplevels(design_F[soilCKsamples,])


spp<-t(cp_otu_F_soil_CK)

F_soil_CK_niche_width <- niche.width(cp_otu_F_soil_CK, method = 'levins')
F_soil_CK_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_F_soil_PE 
cp_otu_F_soil_PE <- cp_otu_F[,soilPEsamples]

keep_cp_otu_F_soil_PE <- which(rowSums(cp_otu_F_soil_PE ) >= 0)
cp_otu_F_soil_PE <- cp_otu_F_soil_PE[keep_cp_otu_F_soil_PE,]

nrow(cp_otu_F_soil_PE)

tax_F_soil_PE <- tax_F[rownames(cp_otu_F_soil_PE),]
design_F_soil_PE <- droplevels(design_F[soilPEsamples,])


spp<-t(cp_otu_F_soil_PE)

F_soil_PE_niche_width <- niche.width(cp_otu_F_soil_PE, method = 'levins')
F_soil_PE_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_F_soil_DE 
cp_otu_F_soil_DE <- cp_otu_F[,soilDEsamples]

keep_cp_otu_F_soil_DE <- which(rowSums(cp_otu_F_soil_DE ) >= 0)
cp_otu_F_soil_DE <- cp_otu_F_soil_DE[keep_cp_otu_F_soil_DE,]

nrow(cp_otu_F_soil_DE)

tax_F_soil_DE <- tax_F[rownames(cp_otu_F_soil_DE),]
design_F_soil_DE <- droplevels(design_F[soilDEsamples,])


spp<-t(cp_otu_F_soil_DE)

F_soil_DE_niche_width <- niche.width(cp_otu_F_soil_DE, method = 'levins')
F_soil_DE_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)





## cp_otu_F_root_CK 
cp_otu_F_root_CK <- cp_otu_F[,rootCKsamples]

keep_cp_otu_F_root_CK <- which(rowSums(cp_otu_F_root_CK ) >= 0)
cp_otu_F_root_CK <- cp_otu_F_root_CK[keep_cp_otu_F_root_CK,]

nrow(cp_otu_F_root_CK)

tax_F_root_CK <- tax_F[rownames(cp_otu_F_root_CK),]
design_F_root_CK <- droplevels(design_F[rootCKsamples,])


spp<-t(cp_otu_F_root_CK)

F_root_CK_niche_width <- niche.width(cp_otu_F_root_CK, method = 'levins')
F_root_CK_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_F_root_PE 
cp_otu_F_root_PE <- cp_otu_F[,rootPEsamples]

keep_cp_otu_F_root_PE <- which(rowSums(cp_otu_F_root_PE ) >= 0)
cp_otu_F_root_PE <- cp_otu_F_root_PE[keep_cp_otu_F_root_PE,]

nrow(cp_otu_F_root_PE)

tax_F_root_PE <- tax_F[rownames(cp_otu_F_root_PE),]
design_F_root_PE <- droplevels(design_F[rootPEsamples,])


spp<-t(cp_otu_F_root_PE)

F_root_PE_niche_width <- niche.width(cp_otu_F_root_PE, method = 'levins')
F_root_PE_niche_width

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 


draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_F_root_DE 
cp_otu_F_root_DE <- cp_otu_F[,rootDEsamples]

keep_cp_otu_F_root_DE <- which(rowSums(cp_otu_F_root_DE ) >= 0)
cp_otu_F_root_DE <- cp_otu_F_root_DE[keep_cp_otu_F_root_DE,]

nrow(cp_otu_F_root_DE)

tax_F_root_DE <- tax_F[rownames(cp_otu_F_root_DE),]
design_F_root_DE <- droplevels(design_F[rootDEsamples,])


spp<-t(cp_otu_F_root_DE)

F_root_DE_niche_width <- niche.width(cp_otu_F_root_DE, method = 'levins')
F_root_DE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_F_end_CK 
cp_otu_F_end_CK <- cp_otu_F[,endCKsamples]

keep_cp_otu_F_end_CK <- which(rowSums(cp_otu_F_end_CK ) >= 0)
cp_otu_F_end_CK <- cp_otu_F_end_CK[keep_cp_otu_F_end_CK,]

nrow(cp_otu_F_end_CK)

tax_F_end_CK <- tax_F[rownames(cp_otu_F_end_CK),]
design_F_end_CK <- droplevels(design_F[endCKsamples,])


spp<-t(cp_otu_F_end_CK)

F_end_CK_niche_width <- niche.width(cp_otu_F_end_CK, method = 'levins')
F_end_CK_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_F_end_PE 
cp_otu_F_end_PE <- cp_otu_F[,endPEsamples]

keep_cp_otu_F_end_PE <- which(rowSums(cp_otu_F_end_PE ) >= 0)
cp_otu_F_end_PE <- cp_otu_F_end_PE[keep_cp_otu_F_end_PE,]

nrow(cp_otu_F_end_PE)

tax_F_end_PE <- tax_F[rownames(cp_otu_F_end_PE),]
design_F_end_PE <- droplevels(design_F[endPEsamples,])


spp<-t(cp_otu_F_end_PE)

F_end_PE_niche_width <- niche.width(cp_otu_F_end_PE, method = 'levins')
F_end_PE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_F_end_DE 
cp_otu_F_end_DE <- cp_otu_F[,endDEsamples]

keep_cp_otu_F_end_DE <- which(rowSums(cp_otu_F_end_DE ) >= 0)
cp_otu_F_end_DE <- cp_otu_F_end_DE[keep_cp_otu_F_end_DE,]

nrow(cp_otu_F_end_DE)

tax_F_end_DE <- tax_F[rownames(cp_otu_F_end_DE),]
design_F_end_DE <- droplevels(design_F[endDEsamples,])


spp<-t(cp_otu_F_end_DE)

F_end_DE_niche_width <- niche.width(cp_otu_F_end_DE, method = 'levins')
F_end_DE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_F_film_PE 
cp_otu_F_film_PE <- cp_otu_F[,PEsamples]

keep_cp_otu_F_film_PE <- which(rowSums(cp_otu_F_film_PE ) >= 0)
cp_otu_F_film_PE <- cp_otu_F_film_PE[keep_cp_otu_F_film_PE,]

nrow(cp_otu_F_film_PE)

tax_F_film_PE <- tax_F[rownames(cp_otu_F_film_PE),]
design_F_film_PE <- droplevels(design_F[PEsamples,])


spp<-t(cp_otu_F_film_PE)

F_film_PE_niche_width <- niche.width(cp_otu_F_film_PE, method = 'levins')
F_film_PE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_F_film_DE 
cp_otu_F_film_DE <- cp_otu_F[,DEsamples]

keep_cp_otu_F_film_DE <- which(rowSums(cp_otu_F_film_DE ) >= 0)
cp_otu_F_film_DE <- cp_otu_F_film_DE[keep_cp_otu_F_film_DE,]

nrow(cp_otu_F_film_DE)

tax_F_film_DE <- tax_F[rownames(cp_otu_F_film_DE),]
design_F_film_DE <- droplevels(design_F[DEsamples,])


spp<-t(cp_otu_F_film_DE)

F_film_DE_niche_width <- niche.width(cp_otu_F_film_DE, method = 'levins')
F_film_DE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)


F_soil_niche_width <- data.frame(CK=t(F_soil_CK_niche_width),PE=t(F_soil_PE_niche_width),DE=t(F_soil_DE_niche_width))

CK <- F_soil_CK_niche_width
PE <- F_soil_CK_niche_width
DE <- F_soil_DE_niche_width

my_comparisons = list( c("CK", "PE"), c("CK", "DE"), c("PE", "DE") )

F_soil_niche_width <- melt(F_soil_niche_width)

F_soil_width_box <- ggplot(F_soil_niche_width, aes(x=factor(variable,levels = c("CK","PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#d2a128","#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



F_root_niche_width <- data.frame(CK=t(F_root_CK_niche_width),PE=t(F_root_PE_niche_width),DE=t(F_root_DE_niche_width))

CK <- F_soil_CK_niche_width
PE <- F_soil_CK_niche_width
DE <- F_soil_DE_niche_width

my_comparisons = list( c("CK", "PE"), c("CK", "DE"), c("PE", "DE") )

F_root_niche_width <- melt(F_root_niche_width)

F_root_width_box <- ggplot(F_root_niche_width, aes(x=factor(variable,levels = c("CK","PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#d2a128","#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



F_end_niche_width <- data.frame(CK=t(F_end_CK_niche_width),PE=t(F_end_PE_niche_width),DE=t(F_end_DE_niche_width))

CK <- F_soil_CK_niche_width
PE <- F_soil_CK_niche_width
DE <- F_soil_DE_niche_width

my_comparisons = list( c("CK", "PE"), c("CK", "DE"), c("PE", "DE") )

F_end_niche_width <- melt(F_end_niche_width)

F_end_width_box <- ggplot(F_end_niche_width, aes(x=factor(variable,levels = c("CK","PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#d2a128","#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



F_film_niche_width <- data.frame(PE=t(F_film_PE_niche_width),DE=t(F_film_DE_niche_width))

PE <- F_soil_CK_niche_width
DE <- F_soil_DE_niche_width

my_comparisons = list( c("PE", "DE") )

F_film_niche_width <- melt(F_film_niche_width)

F_film_width_box <- ggplot(F_film_niche_width, aes(x=factor(variable,levels = c("PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



grid.newpage()
grid.arrange(F_soil_width_box,F_root_width_box,F_end_width_box,F_film_width_box,ncol = 4)




## cp_otu_P_soil_CK 
cp_otu_P_soil_CK <- cp_otu_P[,soilCKsamples]

keep_cp_otu_P_soil_CK <- which(rowSums(cp_otu_P_soil_CK ) >= 0)
cp_otu_P_soil_CK <- cp_otu_P_soil_CK[keep_cp_otu_P_soil_CK,]

nrow(cp_otu_P_soil_CK)

tax_P_soil_CK <- tax_P[rownames(cp_otu_P_soil_CK),]
design_P_soil_CK <- droplevels(design_P[soilCKsamples,])


spp<-t(cp_otu_P_soil_CK)

P_soil_CK_niche_width <- niche.width(cp_otu_P_soil_CK, method = 'levins')
P_soil_CK_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_P_soil_PE 
cp_otu_P_soil_PE <- cp_otu_P[,soilPEsamples]

keep_cp_otu_P_soil_PE <- which(rowSums(cp_otu_P_soil_PE ) >= 0)
cp_otu_P_soil_PE <- cp_otu_P_soil_PE[keep_cp_otu_P_soil_PE,]

nrow(cp_otu_P_soil_PE)

tax_P_soil_PE <- tax_P[rownames(cp_otu_P_soil_PE),]
design_P_soil_PE <- droplevels(design_P[soilPEsamples,])


spp<-t(cp_otu_P_soil_PE)

P_soil_PE_niche_width <- niche.width(cp_otu_P_soil_PE, method = 'levins')
P_soil_PE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_P_soil_DE 
cp_otu_P_soil_DE <- cp_otu_P[,soilDEsamples]

keep_cp_otu_P_soil_DE <- which(rowSums(cp_otu_P_soil_DE ) >= 0)
cp_otu_P_soil_DE <- cp_otu_P_soil_DE[keep_cp_otu_P_soil_DE,]

nrow(cp_otu_P_soil_DE)

tax_P_soil_DE <- tax_P[rownames(cp_otu_P_soil_DE),]
design_P_soil_DE <- droplevels(design_P[soilDEsamples,])


spp<-t(cp_otu_P_soil_DE)

P_soil_DE_niche_width <- niche.width(cp_otu_P_soil_DE, method = 'levins')
P_soil_DE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)





## cp_otu_P_root_CK 
cp_otu_P_root_CK <- cp_otu_P[,rootCKsamples]

keep_cp_otu_P_root_CK <- which(rowSums(cp_otu_P_root_CK ) >= 0)
cp_otu_P_root_CK <- cp_otu_P_root_CK[keep_cp_otu_P_root_CK,]

nrow(cp_otu_P_root_CK)

tax_P_root_CK <- tax_P[rownames(cp_otu_P_root_CK),]
design_P_root_CK <- droplevels(design_P[rootCKsamples,])


spp<-t(cp_otu_P_root_CK)

P_root_CK_niche_width <- niche.width(cp_otu_P_root_CK, method = 'levins')
P_root_CK_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_P_root_PE 
cp_otu_P_root_PE <- cp_otu_P[,rootPEsamples]

keep_cp_otu_P_root_PE <- which(rowSums(cp_otu_P_root_PE ) >= 0)
cp_otu_P_root_PE <- cp_otu_P_root_PE[keep_cp_otu_P_root_PE,]

nrow(cp_otu_P_root_PE)

tax_P_root_PE <- tax_P[rownames(cp_otu_P_root_PE),]
design_P_root_PE <- droplevels(design_P[rootPEsamples,])


spp<-t(cp_otu_P_root_PE)

P_root_PE_niche_width <- niche.width(cp_otu_P_root_PE, method = 'levins')
P_root_PE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_P_root_DE 
cp_otu_P_root_DE <- cp_otu_P[,rootDEsamples]

keep_cp_otu_P_root_DE <- which(rowSums(cp_otu_P_root_DE ) >= 0)
cp_otu_P_root_DE <- cp_otu_P_root_DE[keep_cp_otu_P_root_DE,]

nrow(cp_otu_P_root_DE)

tax_P_root_DE <- tax_P[rownames(cp_otu_P_root_DE),]
design_P_root_DE <- droplevels(design_P[rootDEsamples,])


spp<-t(cp_otu_P_root_DE)

P_root_DE_niche_width <- niche.width(cp_otu_P_root_DE, method = 'levins')
P_root_DE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_P_end_CK 
cp_otu_P_end_CK <- cp_otu_P[,pendCKsamples]

keep_cp_otu_P_end_CK <- which(rowSums(cp_otu_P_end_CK ) >= 0)
cp_otu_P_end_CK <- cp_otu_P_end_CK[keep_cp_otu_P_end_CK,]

nrow(cp_otu_P_end_CK)

tax_P_end_CK <- tax_P[rownames(cp_otu_P_end_CK),]
design_P_end_CK <- droplevels(design_P[pendCKsamples,])


spp<-t(cp_otu_P_end_CK)

P_end_CK_niche_width <- niche.width(cp_otu_P_end_CK, method = 'levins')
P_end_CK_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_P_end_PE 
cp_otu_P_end_PE <- cp_otu_P[,pendPEsamples]

keep_cp_otu_P_end_PE <- which(rowSums(cp_otu_P_end_PE ) >= 0)
cp_otu_P_end_PE <- cp_otu_P_end_PE[keep_cp_otu_P_end_PE,]

nrow(cp_otu_P_end_PE)

tax_P_end_PE <- tax_P[rownames(cp_otu_P_end_PE),]
design_P_end_PE <- droplevels(design_P[pendPEsamples,])


spp<-t(cp_otu_P_end_PE)

P_end_PE_niche_width <- niche.width(cp_otu_P_end_PE, method = 'levins')
P_end_PE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_P_end_DE 
cp_otu_P_end_DE <- cp_otu_P[,pendDEsamples]

keep_cp_otu_P_end_DE <- which(rowSums(cp_otu_P_end_DE ) >= 0)
cp_otu_P_end_DE <- cp_otu_P_end_DE[keep_cp_otu_P_end_DE,]

nrow(cp_otu_P_end_DE)

tax_P_end_DE <- tax_P[rownames(cp_otu_P_end_DE),]
design_P_end_DE <- droplevels(design_P[pendDEsamples,])


spp<-t(cp_otu_P_end_DE)

P_end_DE_niche_width <- niche.width(cp_otu_P_end_DE, method = 'levins')
P_end_DE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)




## cp_otu_P_film_PE 
cp_otu_P_film_PE <- cp_otu_P[,PEsamples]

keep_cp_otu_P_film_PE <- which(rowSums(cp_otu_P_film_PE ) >= 0)
cp_otu_P_film_PE <- cp_otu_P_film_PE[keep_cp_otu_P_film_PE,]

nrow(cp_otu_P_film_PE)

tax_P_film_PE <- tax_P[rownames(cp_otu_P_film_PE),]
design_P_film_PE <- droplevels(design_P[PEsamples,])


spp<-t(cp_otu_P_film_PE)

P_film_PE_niche_width <- niche.width(cp_otu_P_film_PE, method = 'levins')
P_film_PE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'

library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)


P_soil_niche_width <- data.frame(CK=t(P_soil_CK_niche_width),PE=t(P_soil_PE_niche_width),DE=t(P_soil_DE_niche_width))

CK <- P_soil_CK_niche_width
PE <- P_soil_CK_niche_width
DE <- P_soil_DE_niche_width

my_comparisons = list( c("CK", "PE"), c("CK", "DE"), c("PE", "DE") )

P_soil_niche_width <- melt(P_soil_niche_width)

P_soil_width_box <- ggplot(P_soil_niche_width, aes(x=factor(variable,levels = c("CK","PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#d2a128","#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



P_root_niche_width <- data.frame(CK=t(P_root_CK_niche_width),PE=t(P_root_PE_niche_width),DE=t(P_root_DE_niche_width))

CK <- P_soil_CK_niche_width
PE <- P_soil_CK_niche_width
DE <- P_soil_DE_niche_width

my_comparisons = list( c("CK", "PE"), c("CK", "DE"), c("PE", "DE") )

P_root_niche_width <- melt(P_root_niche_width)

P_root_width_box <- ggplot(P_root_niche_width, aes(x=factor(variable,levels = c("CK","PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#d2a128","#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



P_end_niche_width <- data.frame(CK=t(P_end_CK_niche_width),PE=t(P_end_PE_niche_width),DE=t(P_end_DE_niche_width))

CK <- P_soil_CK_niche_width
PE <- P_soil_CK_niche_width
DE <- P_soil_DE_niche_width

my_comparisons = list( c("CK", "PE"), c("CK", "DE"), c("PE", "DE") )

P_end_niche_width <- melt(P_end_niche_width)

P_end_width_box <- ggplot(P_end_niche_width, aes(x=factor(variable,levels = c("CK","PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#d2a128","#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



P_film_niche_width <- data.frame(PE=t(P_film_PE_niche_width),DE=t(P_film_DE_niche_width))

PE <- P_soil_CK_niche_width
DE <- P_soil_DE_niche_width

my_comparisons = list( c("PE", "DE") )

P_film_niche_width <- melt(P_film_niche_width)

P_film_width_box <- ggplot(P_film_niche_width, aes(x=factor(variable,levels = c("PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#41ab6d","#297ec2"))+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



grid.newpage()
grid.arrange(P_soil_width_box,P_root_width_box,P_end_width_box,P_film_width_box,ncol = 4)


grid.newpage()
grid.arrange(B_soil_width_box,B_root_width_box,B_end_width_box,B_film_width_box,F_soil_width_box,F_root_width_box,F_end_width_box,F_film_width_box,P_soil_width_box,P_root_width_box,P_end_width_box,P_film_width_box,ncol = 4)




#####  type network #####
CKsamples <- rownames(design_B)[which(design_B$Group == "CK")]
PEsamples <- rownames(design_B)[which(design_B$Group == "PE")]
DEsamples <- rownames(design_B)[which(design_B$Group == "DE")]

pCKsamples <- rownames(design_P)[which(design_P$Group == "CK")]
pPEsamples <- rownames(design_P)[which(design_P$Group == "PE")]
pDEsamples <- rownames(design_P)[which(design_P$Group == "DE")]


## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
cp_otu_B_CK <- cp_otu_B[,CKsamples]

keep_cp_otu_B_CK <- which(rowSums(cp_otu_B_CK >= 1) >= 4)
cp_otu_B_CK <- cp_otu_B_CK[keep_cp_otu_B_CK,]

nrow(cp_otu_B_CK)

tax_B_CK <- tax_B[rownames(cp_otu_B_CK),]
design_B_CK <- droplevels(design_B[CKsamples,])

B_CK_otu_cor <- rcorr(t(cp_otu_B_CK), type=c("spearman"))

## Create data frame of co-occurring OTUs
B_CK_cor_df <- CorrDF(B_CK_otu_cor$r, B_CK_otu_cor$P)
B_CK_cor_df$padj <- p.adjust(B_CK_cor_df$p, method="none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
B_CK_cor_df_padj <- B_CK_cor_df[which(B_CK_cor_df$cor > 0.7| B_CK_cor_df$cor < -0.7),]
B_CK_cor_df_padj <- B_CK_cor_df_padj[which(B_CK_cor_df_padj$padj < 0.001),]

## Make node attribute table
B_nodeattrib_CK <- data.frame(node = union(B_CK_cor_df_padj$from,B_CK_cor_df_padj$to))
B_nodeattrib_CK$Phylum <- 0

for (i in as.character(B_nodeattrib_CK$node))
{
  if (i %in% rownames(tax_B) == TRUE)
  {B_nodeattrib_CK[B_nodeattrib_CK$node==i,"Phylum"] <- paste(tax_B[i,9:9])}
  else
  { B_nodeattrib_CK[B_nodeattrib_CK$node==i,"Phylum"]<- "NA"}
}


## Create co-occurrence network with igraph
B_CK_net <- graph_from_data_frame(B_CK_cor_df_padj, direct=F, vertices=B_nodeattrib_CK)

## Calculate relative abudnance of OTU nodes
B_CK_ra <- apply(cp_otu_B_CK,1,mean)
B_CK_ra <- B_CK_ra[V(B_CK_net)$name]

## Network properties ##

## Number of nodes in network
length(V(B_CK_net))

## Number of edges in network
length(E(B_CK_net))

## Connections 
bb_occur_PE <- droplevels(B_CK_cor_df_padj[with(B_CK_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor>0,]
nrow(bb_occur_PE)

bb_occur_PE <- droplevels(B_CK_cor_df_padj[with(B_CK_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor<0,]
nrow(bb_occur_PE)

mean(E(B_CK_net))

graph.density(B_CK_net)

set.seed(619)
V(B_CK_net)$modularity <- membership(cluster_fast_greedy(B_CK_net))

B_CK_modularity <- V(B_CK_net)$modularity

B_CK_net_cfg <- cluster_fast_greedy(as.undirected(B_CK_net))
B_CK_net_modules <- sort(table(membership(B_CK_net_cfg)),decr=T)

transitivity(B_CK_net)

B_CK_net_betweenness <- betweenness(B_CK_net,normalized = T)
B_CK_net_degree <- igraph::degree(B_CK_net)

B_CK_net_modularity <- modularity(cluster_fast_greedy(B_CK_net))

plot(B_CK_net_betweenness,B_CK_net_degree)

B_CK_topological <- cbind(B_CK_net_betweenness,B_CK_net_degree,B_CK_modularity)

B_CK_topological_data <- as.data.frame(B_CK_topological)
colnames(B_CK_topological_data) <- c("betweenness","degree","modularity")

adjacency_B_CK <- get.adjacency(B_CK_net,attr = "cor",sparse = FALSE)

B_CK_zi_pi <- zi.pi(B_CK_topological_data, adjacency_B_CK, degree = 'degree', modularity_class = 'modularity')
head(B_CK_zi_pi)

B_CK_zi_pi <- na.omit(B_CK_zi_pi) 
B_CK_zi_pi[which(B_CK_zi_pi$within_module_connectivities < 2.5 & B_CK_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
B_CK_zi_pi[which(B_CK_zi_pi$within_module_connectivities < 2.5 & B_CK_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
B_CK_zi_pi[which(B_CK_zi_pi$within_module_connectivities > 2.5 & B_CK_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
B_CK_zi_pi[which(B_CK_zi_pi$within_module_connectivities > 2.5 & B_CK_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

B_CK_zi_pi_plot <- ggplot(B_CK_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)


## Calculate network node degrees/max/min
B_CK_deg <- sort(igraph::degree(B_CK_net,mode="all"),decr=T)
max(B_CK_deg)
mean(B_CK_deg)

## Set node colors based upon sensitivity to management system
unique(V(B_CK_net)$Phylum)
V(B_CK_net)$color <- V(B_CK_net)$Phylum

V(B_CK_net)$frame.color <- V(B_CK_net)$color

## Set node shape
V(B_CK_net)$shape <- "circle"

## Set node size
V(B_CK_net)$size <- 3



## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
cp_otu_B_PE <- cp_otu_B[,PEsamples]

keep_cp_otu_B_PE <- which(rowSums(cp_otu_B_PE >= 1) >= 4)
cp_otu_B_PE <- cp_otu_B_PE[keep_cp_otu_B_PE,]

nrow(cp_otu_B_PE)

tax_B_PE <- tax_B[rownames(cp_otu_B_PE),]
design_B_PE <- droplevels(design_B[PEsamples,])

B_PE_otu_cor <- rcorr(t(cp_otu_B_PE), type=c("spearman"))

## Create data frame of co-occurring OTUs
B_PE_cor_df <- CorrDF(B_PE_otu_cor$r, B_PE_otu_cor$P)
B_PE_cor_df$padj <- p.adjust(B_PE_cor_df$p, method="none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
B_PE_cor_df_padj <- B_PE_cor_df[which(B_PE_cor_df$cor > 0.7| B_PE_cor_df$cor < -0.7),]
B_PE_cor_df_padj <- B_PE_cor_df_padj[which(B_PE_cor_df_padj$padj < 0.001),]

## Make node attribute table
B_nodeattrib_PE <- data.frame(node = union(B_PE_cor_df_padj$from,B_PE_cor_df_padj$to))
B_nodeattrib_PE$Phylum <- 0

for (i in as.character(B_nodeattrib_PE$node))
{
  if (i %in% rownames(tax_B) == TRUE)
  {B_nodeattrib_PE[B_nodeattrib_PE$node==i,"Phylum"] <- paste(tax_B[i,9:9])}
  else
  { B_nodeattrib_PE[B_nodeattrib_PE$node==i,"Phylum"]<- "NA"}
}


## Create co-occurrence network with igraph
B_PE_net <- graph_from_data_frame(B_PE_cor_df_padj, direct=F, vertices=B_nodeattrib_PE)

## Calculate relative abudnance of OTU nodes
B_PE_ra <- apply(cp_otu_B_PE,1,mean)
B_PE_ra <- B_PE_ra[V(B_PE_net)$name]

## Network properties ##

## Number of nodes in network
length(V(B_PE_net))

## Number of edges in network
length(E(B_PE_net))

## Connections 
bb_occur_PE <- droplevels(B_PE_cor_df_padj[with(B_PE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor>0,]
nrow(bb_occur_PE)

bb_occur_PE <- droplevels(B_PE_cor_df_padj[with(B_PE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor<0,]
nrow(bb_occur_PE)

mean(E(B_PE_net))

graph.density(B_PE_net)

set.seed(619)
V(B_PE_net)$modularity <- membership(cluster_fast_greedy(B_PE_net))

B_PE_modularity <- V(B_PE_net)$modularity

B_PE_net_cfg <- cluster_fast_greedy(as.undirected(B_PE_net))
B_PE_net_modules <- sort(table(membership(B_PE_net_cfg)),decr=T)

transitivity(B_PE_net)

B_PE_net_betweenness <- betweenness(B_PE_net,normalized = T)
B_PE_net_degree <- igraph::degree(B_PE_net)

B_PE_net_modularity <- modularity(cluster_fast_greedy(B_PE_net))

plot(B_PE_net_betweenness,B_PE_net_degree)

B_PE_topological <- cbind(B_PE_net_betweenness,B_PE_net_degree,B_PE_modularity)

B_PE_topological_data <- as.data.frame(B_PE_topological)
colnames(B_PE_topological_data) <- c("betweenness","degree","modularity")

adjacency_B_PE <- get.adjacency(B_PE_net,attr = "cor",sparse = FALSE)

B_PE_zi_pi <- zi.pi(B_PE_topological_data, adjacency_B_PE, degree = 'degree', modularity_class = 'modularity')
head(B_PE_zi_pi)

B_PE_zi_pi <- na.omit(B_PE_zi_pi) 
B_PE_zi_pi[which(B_PE_zi_pi$within_module_connectivities < 2.5 & B_PE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
B_PE_zi_pi[which(B_PE_zi_pi$within_module_connectivities < 2.5 & B_PE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
B_PE_zi_pi[which(B_PE_zi_pi$within_module_connectivities > 2.5 & B_PE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
B_PE_zi_pi[which(B_PE_zi_pi$within_module_connectivities > 2.5 & B_PE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

B_PE_zi_pi_plot <- ggplot(B_PE_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)


## Calculate network node degrees/max/min
B_PE_deg <- sort(igraph::degree(B_PE_net,mode="all"),decr=T)
max(B_PE_deg)
mean(B_PE_deg)

## Set node colors based upon sensitivity to management system
unique(V(B_PE_net)$Phylum)
V(B_PE_net)$color <- V(B_PE_net)$Phylum

V(B_PE_net)$frame.color <- V(B_PE_net)$color

## Set node shape
V(B_PE_net)$shape <- "circle"

## Set node size
V(B_PE_net)$size <- 3




## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
cp_otu_B_DE <- cp_otu_B[,DEsamples]

keep_cp_otu_B_DE <- which(rowSums(cp_otu_B_DE >= 1) >= 4)
cp_otu_B_DE <- cp_otu_B_DE[keep_cp_otu_B_DE,]

nrow(cp_otu_B_DE)

tax_B_DE <- tax_B[rownames(cp_otu_B_DE),]
design_B_DE <- droplevels(design_B[DEsamples,])

B_DE_otu_cor <- rcorr(t(cp_otu_B_DE), type=c("spearman"))

## Create data frame of co-occurring OTUs
B_DE_cor_df <- CorrDF(B_DE_otu_cor$r, B_DE_otu_cor$P)
B_DE_cor_df$padj <- p.adjust(B_DE_cor_df$p, method="none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
B_DE_cor_df_padj <- B_DE_cor_df[which(B_DE_cor_df$cor > 0.7| B_DE_cor_df$cor < -0.7),]
B_DE_cor_df_padj <- B_DE_cor_df_padj[which(B_DE_cor_df_padj$padj < 0.001),]

## Make node attribute table
B_nodeattrib_DE <- data.frame(node = union(B_DE_cor_df_padj$from,B_DE_cor_df_padj$to))
B_nodeattrib_DE$Phylum <- 0

for (i in as.character(B_nodeattrib_DE$node))
{
  if (i %in% rownames(tax_B) == TRUE)
  {B_nodeattrib_DE[B_nodeattrib_DE$node==i,"Phylum"] <- paste(tax_B[i,9:9])}
  else
  { B_nodeattrib_DE[B_nodeattrib_DE$node==i,"Phylum"]<- "NA"}
}


## Create co-occurrence network with igraph
B_DE_net <- graph_from_data_frame(B_DE_cor_df_padj, direct=F, vertices=B_nodeattrib_DE)

## Calculate relative abudnance of OTU nodes
B_DE_ra <- apply(cp_otu_B_DE,1,mean)
B_DE_ra <- B_DE_ra[V(B_DE_net)$name]

## Network properties ##

## Number of nodes in network
length(V(B_DE_net))

## Number of edges in network
length(E(B_DE_net))

## Connections 
bb_occur_PE <- droplevels(B_DE_cor_df_padj[with(B_DE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor>0,]
nrow(bb_occur_PE)

bb_occur_PE <- droplevels(B_DE_cor_df_padj[with(B_DE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor<0,]
nrow(bb_occur_PE)

mean(E(B_DE_net))

graph.density(B_DE_net)

set.seed(619)
V(B_DE_net)$modularity <- membership(cluster_fast_greedy(B_DE_net))

B_DE_modularity <- V(B_DE_net)$modularity

B_DE_net_cfg <- cluster_fast_greedy(as.undirected(B_DE_net))
B_DE_net_modules <- sort(table(membership(B_DE_net_cfg)),decr=T)

transitivity(B_DE_net)

B_DE_net_betweenness <- betweenness(B_DE_net,normalized = T)
B_DE_net_degree <- igraph::degree(B_DE_net)

B_DE_net_modularity <- modularity(cluster_fast_greedy(B_DE_net))

plot(B_DE_net_betweenness,B_DE_net_degree)

B_DE_topological <- cbind(B_DE_net_betweenness,B_DE_net_degree,B_DE_modularity)

B_DE_topological_data <- as.data.frame(B_DE_topological)
colnames(B_DE_topological_data) <- c("betweenness","degree","modularity")

adjacency_B_DE <- get.adjacency(B_DE_net,attr = "cor",sparse = FALSE)

B_DE_zi_pi <- zi.pi(B_DE_topological_data, adjacency_B_DE, degree = 'degree', modularity_class = 'modularity')
head(B_DE_zi_pi)

B_DE_zi_pi <- na.omit(B_DE_zi_pi) 
B_DE_zi_pi[which(B_DE_zi_pi$within_module_connectivities < 2.5 & B_DE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
B_DE_zi_pi[which(B_DE_zi_pi$within_module_connectivities < 2.5 & B_DE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
B_DE_zi_pi[which(B_DE_zi_pi$within_module_connectivities > 2.5 & B_DE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
B_DE_zi_pi[which(B_DE_zi_pi$within_module_connectivities > 2.5 & B_DE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

B_DE_zi_pi_plot <- ggplot(B_DE_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)


## Calculate network node degrees/max/min
B_DE_deg <- sort(igraph::degree(B_DE_net,mode="all"),decr=T)
max(B_DE_deg)
mean(B_DE_deg)

## Set node colors based upon sensitivity to management system
unique(V(B_DE_net)$Phylum)
V(B_DE_net)$color <- V(B_DE_net)$Phylum

V(B_DE_net)$frame.color <- V(B_DE_net)$color

## Set node shape
V(B_DE_net)$shape <- "circle"

## Set node size
V(B_DE_net)$size <- 3



## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
cp_otu_F_CK <- cp_otu_F[,CKsamples]

keep_cp_otu_F_CK <- which(rowSums(cp_otu_F_CK >= 1) >= 4)
cp_otu_F_CK <- cp_otu_F_CK[keep_cp_otu_F_CK,]

nrow(cp_otu_F_CK)

tax_F_CK <- tax_F[rownames(cp_otu_F_CK),]
design_F_CK <- droplevels(design_F[CKsamples,])

F_CK_otu_cor <- rcorr(t(cp_otu_F_CK), type=c("spearman"))

## Create data frame of co-occurring OTUs
F_CK_cor_df <- CorrDF(F_CK_otu_cor$r, F_CK_otu_cor$P)
F_CK_cor_df$padj <- p.adjust(F_CK_cor_df$p, method="none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
F_CK_cor_df_padj <- F_CK_cor_df[which(F_CK_cor_df$cor > 0.7| F_CK_cor_df$cor < -0.7),]
F_CK_cor_df_padj <- F_CK_cor_df_padj[which(F_CK_cor_df_padj$padj < 0.001),]

## Make node attribute table
F_nodeattrib_CK <- data.frame(node = union(F_CK_cor_df_padj$from,F_CK_cor_df_padj$to))
F_nodeattrib_CK$Phylum <- 0

for (i in as.character(F_nodeattrib_CK$node))
{
  if (i %in% rownames(tax_F) == TRUE)
  {F_nodeattrib_CK[F_nodeattrib_CK$node==i,"Phylum"] <- paste(tax_F[i,9:9])}
  else
  { F_nodeattrib_CK[F_nodeattrib_CK$node==i,"Phylum"]<- "NA"}
}


## Create co-occurrence network with igraph
F_CK_net <- graph_from_data_frame(F_CK_cor_df_padj, direct=F, vertices=F_nodeattrib_CK)

## Calculate relative abudnance of OTU nodes
F_CK_ra <- apply(cp_otu_F_CK,1,mean)
F_CK_ra <- F_CK_ra[V(F_CK_net)$name]

## Network properties ##

## Number of nodes in network
length(V(F_CK_net))

## Number of edges in network
length(E(F_CK_net))

## Connections 
bb_occur_PE <- droplevels(F_CK_cor_df_padj[with(F_CK_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor>0,]
nrow(bb_occur_PE)

bb_occur_PE <- droplevels(F_CK_cor_df_padj[with(F_CK_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor<0,]
nrow(bb_occur_PE)

mean(E(F_CK_net))

graph.density(F_CK_net)

set.seed(619)
V(F_CK_net)$modularity <- membership(cluster_fast_greedy(F_CK_net))

F_CK_modularity <- V(F_CK_net)$modularity

F_CK_net_cfg <- cluster_fast_greedy(as.undirected(F_CK_net))
F_CK_net_modules <- sort(table(membership(F_CK_net_cfg)),decr=T)

transitivity(F_CK_net)

F_CK_net_betweenness <- betweenness(F_CK_net,normalized = T)
F_CK_net_degree <- igraph::degree(F_CK_net)

F_CK_net_modularity <- modularity(cluster_fast_greedy(F_CK_net))

plot(F_CK_net_betweenness,F_CK_net_degree)

F_CK_topological <- cbind(F_CK_net_betweenness,F_CK_net_degree,F_CK_modularity)

F_CK_topological_data <- as.data.frame(F_CK_topological)
colnames(F_CK_topological_data) <- c("betweenness","degree","modularity")

adjacency_F_CK <- get.adjacency(F_CK_net,attr = "cor",sparse = FALSE)

F_CK_zi_pi <- zi.pi(F_CK_topological_data, adjacency_F_CK, degree = 'degree', modularity_class = 'modularity')
head(F_CK_zi_pi)

F_CK_zi_pi <- na.omit(F_CK_zi_pi) 
F_CK_zi_pi[which(F_CK_zi_pi$within_module_connectivities < 2.5 & F_CK_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
F_CK_zi_pi[which(F_CK_zi_pi$within_module_connectivities < 2.5 & F_CK_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
F_CK_zi_pi[which(F_CK_zi_pi$within_module_connectivities > 2.5 & F_CK_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
F_CK_zi_pi[which(F_CK_zi_pi$within_module_connectivities > 2.5 & F_CK_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

F_CK_zi_pi_plot <- ggplot(F_CK_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)


## Calculate network node degrees/max/min
F_CK_deg <- sort(igraph::degree(F_CK_net,mode="all"),decr=T)
max(F_CK_deg)
mean(F_CK_deg)

## Set node colors based upon sensitivity to management system
unique(V(F_CK_net)$Phylum)
V(F_CK_net)$color <- V(F_CK_net)$Phylum

V(F_CK_net)$frame.color <- V(F_CK_net)$color

## Set node shape
V(F_CK_net)$shape <- "circle"

## Set node size
V(F_CK_net)$size <- 3



## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
cp_otu_F_PE <- cp_otu_F[,PEsamples]

keep_cp_otu_F_PE <- which(rowSums(cp_otu_F_PE >= 1) >= 4)
cp_otu_F_PE <- cp_otu_F_PE[keep_cp_otu_F_PE,]

nrow(cp_otu_F_PE)

tax_F_PE <- tax_F[rownames(cp_otu_F_PE),]
design_F_PE <- droplevels(design_F[PEsamples,])

F_PE_otu_cor <- rcorr(t(cp_otu_F_PE), type=c("spearman"))

## Create data frame of co-occurring OTUs
F_PE_cor_df <- CorrDF(F_PE_otu_cor$r, F_PE_otu_cor$P)
F_PE_cor_df$padj <- p.adjust(F_PE_cor_df$p, method="none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
F_PE_cor_df_padj <- F_PE_cor_df[which(F_PE_cor_df$cor > 0.7| F_PE_cor_df$cor < -0.7),]
F_PE_cor_df_padj <- F_PE_cor_df_padj[which(F_PE_cor_df_padj$padj < 0.001),]

## Make node attribute table
F_nodeattrib_PE <- data.frame(node = union(F_PE_cor_df_padj$from,F_PE_cor_df_padj$to))
F_nodeattrib_PE$Phylum <- 0

for (i in as.character(F_nodeattrib_PE$node))
{
  if (i %in% rownames(tax_F) == TRUE)
  {F_nodeattrib_PE[F_nodeattrib_PE$node==i,"Phylum"] <- paste(tax_F[i,9:9])}
  else
  { F_nodeattrib_PE[F_nodeattrib_PE$node==i,"Phylum"]<- "NA"}
}


## Create co-occurrence network with igraph
F_PE_net <- graph_from_data_frame(F_PE_cor_df_padj, direct=F, vertices=F_nodeattrib_PE)

## Calculate relative abudnance of OTU nodes
F_PE_ra <- apply(cp_otu_F_PE,1,mean)
F_PE_ra <- F_PE_ra[V(F_PE_net)$name]

## Network properties ##

## Number of nodes in network
length(V(F_PE_net))

## Number of edges in network
length(E(F_PE_net))

## Connections 
bb_occur_PE <- droplevels(F_PE_cor_df_padj[with(F_PE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor>0,]
nrow(bb_occur_PE)

bb_occur_PE <- droplevels(F_PE_cor_df_padj[with(F_PE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor<0,]
nrow(bb_occur_PE)

mean(E(F_PE_net))

graph.density(F_PE_net)

set.seed(619)
V(F_PE_net)$modularity <- membership(cluster_fast_greedy(F_PE_net))

F_PE_modularity <- V(F_PE_net)$modularity

F_PE_net_cfg <- cluster_fast_greedy(as.undirected(F_PE_net))
F_PE_net_modules <- sort(table(membership(F_PE_net_cfg)),decr=T)

transitivity(F_PE_net)

F_PE_net_betweenness <- betweenness(F_PE_net,normalized = T)
F_PE_net_degree <- igraph::degree(F_PE_net)

F_PE_net_modularity <- modularity(cluster_fast_greedy(F_PE_net))

plot(F_PE_net_betweenness,F_PE_net_degree)

F_PE_topological <- cbind(F_PE_net_betweenness,F_PE_net_degree,F_PE_modularity)

F_PE_topological_data <- as.data.frame(F_PE_topological)
colnames(F_PE_topological_data) <- c("betweenness","degree","modularity")

adjacency_F_PE <- get.adjacency(F_PE_net,attr = "cor",sparse = FALSE)

F_PE_zi_pi <- zi.pi(F_PE_topological_data, adjacency_F_PE, degree = 'degree', modularity_class = 'modularity')
head(F_PE_zi_pi)

F_PE_zi_pi <- na.omit(F_PE_zi_pi) 
F_PE_zi_pi[which(F_PE_zi_pi$within_module_connectivities < 2.5 & F_PE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
F_PE_zi_pi[which(F_PE_zi_pi$within_module_connectivities < 2.5 & F_PE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
F_PE_zi_pi[which(F_PE_zi_pi$within_module_connectivities > 2.5 & F_PE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
F_PE_zi_pi[which(F_PE_zi_pi$within_module_connectivities > 2.5 & F_PE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

F_PE_zi_pi_plot <- ggplot(F_PE_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)


## Calculate network node degrees/max/min
F_PE_deg <- sort(igraph::degree(F_PE_net,mode="all"),decr=T)
max(F_PE_deg)
mean(F_PE_deg)

## Set node colors based upon sensitivity to management system
unique(V(F_PE_net)$Phylum)
V(F_PE_net)$color <- V(F_PE_net)$Phylum

V(F_PE_net)$frame.color <- V(F_PE_net)$color

## Set node shape
V(F_PE_net)$shape <- "circle"

## Set node size
V(F_PE_net)$size <- 3




## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
cp_otu_F_DE <- cp_otu_F[,DEsamples]

keep_cp_otu_F_DE <- which(rowSums(cp_otu_F_DE >= 1) >= 4)
cp_otu_F_DE <- cp_otu_F_DE[keep_cp_otu_F_DE,]

nrow(cp_otu_F_DE)

tax_F_DE <- tax_F[rownames(cp_otu_F_DE),]
design_F_DE <- droplevels(design_F[DEsamples,])

F_DE_otu_cor <- rcorr(t(cp_otu_F_DE), type=c("spearman"))

## Create data frame of co-occurring OTUs
F_DE_cor_df <- CorrDF(F_DE_otu_cor$r, F_DE_otu_cor$P)
F_DE_cor_df$padj <- p.adjust(F_DE_cor_df$p, method="none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
F_DE_cor_df_padj <- F_DE_cor_df[which(F_DE_cor_df$cor > 0.7| F_DE_cor_df$cor < -0.7),]
F_DE_cor_df_padj <- F_DE_cor_df_padj[which(F_DE_cor_df_padj$padj < 0.001),]

## Make node attribute table
F_nodeattrib_DE <- data.frame(node = union(F_DE_cor_df_padj$from,F_DE_cor_df_padj$to))
F_nodeattrib_DE$Phylum <- 0

for (i in as.character(F_nodeattrib_DE$node))
{
  if (i %in% rownames(tax_F) == TRUE)
  {F_nodeattrib_DE[F_nodeattrib_DE$node==i,"Phylum"] <- paste(tax_F[i,9:9])}
  else
  { F_nodeattrib_DE[F_nodeattrib_DE$node==i,"Phylum"]<- "NA"}
}


## Create co-occurrence network with igraph
F_DE_net <- graph_from_data_frame(F_DE_cor_df_padj, direct=F, vertices=F_nodeattrib_DE)

## Calculate relative abudnance of OTU nodes
F_DE_ra <- apply(cp_otu_F_DE,1,mean)
F_DE_ra <- F_DE_ra[V(F_DE_net)$name]

## Network properties ##

## Number of nodes in network
length(V(F_DE_net))

## Number of edges in network
length(E(F_DE_net))

## Connections 
bb_occur_PE <- droplevels(F_DE_cor_df_padj[with(F_DE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor>0,]
nrow(bb_occur_PE)

bb_occur_PE <- droplevels(F_DE_cor_df_padj[with(F_DE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor<0,]
nrow(bb_occur_PE)

mean(E(F_DE_net))

graph.density(F_DE_net)

set.seed(619)
V(F_DE_net)$modularity <- membership(cluster_fast_greedy(F_DE_net))

F_DE_modularity <- V(F_DE_net)$modularity

F_DE_net_cfg <- cluster_fast_greedy(as.undirected(F_DE_net))
F_DE_net_modules <- sort(table(membership(F_DE_net_cfg)),decr=T)

transitivity(F_DE_net)

F_DE_net_betweenness <- betweenness(F_DE_net,normalized = T)
F_DE_net_degree <- igraph::degree(F_DE_net)

F_DE_net_modularity <- modularity(cluster_fast_greedy(F_DE_net))

plot(F_DE_net_betweenness,F_DE_net_degree)

F_DE_topological <- cbind(F_DE_net_betweenness,F_DE_net_degree,F_DE_modularity)

F_DE_topological_data <- as.data.frame(F_DE_topological)
colnames(F_DE_topological_data) <- c("betweenness","degree","modularity")

adjacency_F_DE <- get.adjacency(F_DE_net,attr = "cor",sparse = FALSE)

F_DE_zi_pi <- zi.pi(F_DE_topological_data, adjacency_F_DE, degree = 'degree', modularity_class = 'modularity')
head(F_DE_zi_pi)

F_DE_zi_pi <- na.omit(F_DE_zi_pi) 
F_DE_zi_pi[which(F_DE_zi_pi$within_module_connectivities < 2.5 & F_DE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
F_DE_zi_pi[which(F_DE_zi_pi$within_module_connectivities < 2.5 & F_DE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
F_DE_zi_pi[which(F_DE_zi_pi$within_module_connectivities > 2.5 & F_DE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
F_DE_zi_pi[which(F_DE_zi_pi$within_module_connectivities > 2.5 & F_DE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

F_DE_zi_pi_plot <- ggplot(F_DE_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)


## Calculate network node degrees/max/min
F_DE_deg <- sort(igraph::degree(F_DE_net,mode="all"),decr=T)
max(F_DE_deg)
mean(F_DE_deg)

## Set node colors based upon sensitivity to management system
unique(V(F_DE_net)$Phylum)
V(F_DE_net)$color <- V(F_DE_net)$Phylum

V(F_DE_net)$frame.color <- V(F_DE_net)$color

## Set node shape
V(F_DE_net)$shape <- "circle"

## Set node size
V(F_DE_net)$size <- 3



## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
cp_otu_P_CK <- cp_otu_P[,pCKsamples]

keep_cp_otu_P_CK <- which(rowSums(cp_otu_P_CK >= 1) >= 4)
cp_otu_P_CK <- cp_otu_P_CK[keep_cp_otu_P_CK,]

nrow(cp_otu_P_CK)

tax_P_CK <- tax_P[rownames(cp_otu_P_CK),]
design_P_CK <- droplevels(design_P[CKsamples,])

P_CK_otu_cor <- rcorr(t(cp_otu_P_CK), type=c("spearman"))

## Create data frame of co-occurring OTUs
P_CK_cor_df <- CorrDF(P_CK_otu_cor$r, P_CK_otu_cor$P)
P_CK_cor_df$padj <- p.adjust(P_CK_cor_df$p, method="none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
P_CK_cor_df_padj <- P_CK_cor_df[which(P_CK_cor_df$cor > 0.7| P_CK_cor_df$cor < -0.7),]
P_CK_cor_df_padj <- P_CK_cor_df_padj[which(P_CK_cor_df_padj$padj < 0.001),]

## Make node attribute table
P_nodeattrib_CK <- data.frame(node = union(P_CK_cor_df_padj$from,P_CK_cor_df_padj$to))
P_nodeattrib_CK$Phylum <- 0

for (i in as.character(P_nodeattrib_CK$node))
{
  if (i %in% rownames(tax_P) == TRUE)
  {P_nodeattrib_CK[P_nodeattrib_CK$node==i,"Phylum"] <- paste(tax_P[i,9:9])}
  else
  { P_nodeattrib_CK[P_nodeattrib_CK$node==i,"Phylum"]<- "NA"}
}


## Create co-occurrence network with igraph
P_CK_net <- graph_from_data_frame(P_CK_cor_df_padj, direct=F, vertices=P_nodeattrib_CK)

## Calculate relative abudnance of OTU nodes
P_CK_ra <- apply(cp_otu_P_CK,1,mean)
P_CK_ra <- P_CK_ra[V(P_CK_net)$name]

## Network properties ##

## Number of nodes in network
length(V(P_CK_net))

## Number of edges in network
length(E(P_CK_net))

## Connections 
bb_occur_PE <- droplevels(P_CK_cor_df_padj[with(P_CK_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor>0,]
nrow(bb_occur_PE)

bb_occur_PE <- droplevels(P_CK_cor_df_padj[with(P_CK_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor<0,]
nrow(bb_occur_PE)

mean(E(P_CK_net))

graph.density(P_CK_net)

set.seed(619)
V(P_CK_net)$modularity <- membership(cluster_fast_greedy(P_CK_net))

P_CK_modularity <- V(P_CK_net)$modularity

P_CK_net_cfg <- cluster_fast_greedy(as.undirected(P_CK_net))
P_CK_net_modules <- sort(table(membership(P_CK_net_cfg)),decr=T)

transitivity(P_CK_net)

P_CK_net_betweenness <- betweenness(P_CK_net,normalized = T)
P_CK_net_degree <- igraph::degree(P_CK_net)

P_CK_net_modularity <- modularity(cluster_fast_greedy(P_CK_net))

plot(P_CK_net_betweenness,P_CK_net_degree)

P_CK_topological <- cbind(P_CK_net_betweenness,P_CK_net_degree,P_CK_modularity)

P_CK_topological_data <- as.data.frame(P_CK_topological)
colnames(P_CK_topological_data) <- c("betweenness","degree","modularity")

adjacency_P_CK <- get.adjacency(P_CK_net,attr = "cor",sparse = FALSE)

P_CK_zi_pi <- zi.pi(P_CK_topological_data, adjacency_P_CK, degree = 'degree', modularity_class = 'modularity')
head(P_CK_zi_pi)

P_CK_zi_pi <- na.omit(P_CK_zi_pi) 
P_CK_zi_pi[which(P_CK_zi_pi$within_module_connectivities < 2.5 & P_CK_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
P_CK_zi_pi[which(P_CK_zi_pi$within_module_connectivities < 2.5 & P_CK_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
P_CK_zi_pi[which(P_CK_zi_pi$within_module_connectivities > 2.5 & P_CK_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
P_CK_zi_pi[which(P_CK_zi_pi$within_module_connectivities > 2.5 & P_CK_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

P_CK_zi_pi_plot <- ggplot(P_CK_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)


## Calculate network node degrees/max/min
P_CK_deg <- sort(igraph::degree(P_CK_net,mode="all"),decr=T)
max(P_CK_deg)
mean(P_CK_deg)

## Set node colors based upon sensitivity to management system
unique(V(P_CK_net)$Phylum)
V(P_CK_net)$color <- V(P_CK_net)$Phylum

V(P_CK_net)$frame.color <- V(P_CK_net)$color

## Set node shape
V(P_CK_net)$shape <- "circle"

## Set node size
V(P_CK_net)$size <- 3



## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
cp_otu_P_PE <- cp_otu_P[,pPEsamples]

keep_cp_otu_P_PE <- which(rowSums(cp_otu_P_PE >= 1) >= 4)
cp_otu_P_PE <- cp_otu_P_PE[keep_cp_otu_P_PE,]

nrow(cp_otu_P_PE)

tax_P_PE <- tax_P[rownames(cp_otu_P_PE),]
design_P_PE <- droplevels(design_P[PEsamples,])

P_PE_otu_cor <- rcorr(t(cp_otu_P_PE), type=c("spearman"))

## Create data frame of co-occurring OTUs
P_PE_cor_df <- CorrDF(P_PE_otu_cor$r, P_PE_otu_cor$P)
P_PE_cor_df$padj <- p.adjust(P_PE_cor_df$p, method="none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
P_PE_cor_df_padj <- P_PE_cor_df[which(P_PE_cor_df$cor > 0.7| P_PE_cor_df$cor < -0.7),]
P_PE_cor_df_padj <- P_PE_cor_df_padj[which(P_PE_cor_df_padj$padj < 0.001),]

## Make node attribute table
P_nodeattrib_PE <- data.frame(node = union(P_PE_cor_df_padj$from,P_PE_cor_df_padj$to))
P_nodeattrib_PE$Phylum <- 0

for (i in as.character(P_nodeattrib_PE$node))
{
  if (i %in% rownames(tax_P) == TRUE)
  {P_nodeattrib_PE[P_nodeattrib_PE$node==i,"Phylum"] <- paste(tax_P[i,9:9])}
  else
  { P_nodeattrib_PE[P_nodeattrib_PE$node==i,"Phylum"]<- "NA"}
}


## Create co-occurrence network with igraph
P_PE_net <- graph_from_data_frame(P_PE_cor_df_padj, direct=F, vertices=P_nodeattrib_PE)

## Calculate relative abudnance of OTU nodes
P_PE_ra <- apply(cp_otu_P_PE,1,mean)
P_PE_ra <- P_PE_ra[V(P_PE_net)$name]

## Network properties ##

## Number of nodes in network
length(V(P_PE_net))

## Number of edges in network
length(E(P_PE_net))

## Connections 
bb_occur_PE <- droplevels(P_PE_cor_df_padj[with(P_PE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor>0,]
nrow(bb_occur_PE)

bb_occur_PE <- droplevels(P_PE_cor_df_padj[with(P_PE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor<0,]
nrow(bb_occur_PE)

mean(E(P_PE_net))

graph.density(P_PE_net)

set.seed(619)
V(P_PE_net)$modularity <- membership(cluster_fast_greedy(P_PE_net))

P_PE_modularity <- V(P_PE_net)$modularity

P_PE_net_cfg <- cluster_fast_greedy(as.undirected(P_PE_net))
P_PE_net_modules <- sort(table(membership(P_PE_net_cfg)),decr=T)

transitivity(P_PE_net)

P_PE_net_betweenness <- betweenness(P_PE_net,normalized = T)
P_PE_net_degree <- igraph::degree(P_PE_net)

P_PE_net_modularity <- modularity(cluster_fast_greedy(P_PE_net))

plot(P_PE_net_betweenness,P_PE_net_degree)

P_PE_topological <- cbind(P_PE_net_betweenness,P_PE_net_degree,P_PE_modularity)

P_PE_topological_data <- as.data.frame(P_PE_topological)
colnames(P_PE_topological_data) <- c("betweenness","degree","modularity")

adjacency_P_PE <- get.adjacency(P_PE_net,attr = "cor",sparse = FALSE)

P_PE_zi_pi <- zi.pi(P_PE_topological_data, adjacency_P_PE, degree = 'degree', modularity_class = 'modularity')
head(P_PE_zi_pi)

P_PE_zi_pi <- na.omit(P_PE_zi_pi) 
P_PE_zi_pi[which(P_PE_zi_pi$within_module_connectivities < 2.5 & P_PE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
P_PE_zi_pi[which(P_PE_zi_pi$within_module_connectivities < 2.5 & P_PE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
P_PE_zi_pi[which(P_PE_zi_pi$within_module_connectivities > 2.5 & P_PE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
P_PE_zi_pi[which(P_PE_zi_pi$within_module_connectivities > 2.5 & P_PE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

P_PE_zi_pi_plot <- ggplot(P_PE_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)


## Calculate network node degrees/max/min
P_PE_deg <- sort(igraph::degree(P_PE_net,mode="all"),decr=T)
max(P_PE_deg)
mean(P_PE_deg)

## Set node colors based upon sensitivity to management system
unique(V(P_PE_net)$Phylum)
V(P_PE_net)$color <- V(P_PE_net)$Phylum

V(P_PE_net)$frame.color <- V(P_PE_net)$color

## Set node shape
V(P_PE_net)$shape <- "circle"

## Set node size
V(P_PE_net)$size <- 3




## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
cp_otu_P_DE <- cp_otu_P[,pDEsamples]

keep_cp_otu_P_DE <- which(rowSums(cp_otu_P_DE >= 1) >= 4)
cp_otu_P_DE <- cp_otu_P_DE[keep_cp_otu_P_DE,]

nrow(cp_otu_P_DE)

tax_P_DE <- tax_P[rownames(cp_otu_P_DE),]
design_P_DE <- droplevels(design_P[DEsamples,])

P_DE_otu_cor <- rcorr(t(cp_otu_P_DE), type=c("spearman"))

## Create data frame of co-occurring OTUs
P_DE_cor_df <- CorrDF(P_DE_otu_cor$r, P_DE_otu_cor$P)
P_DE_cor_df$padj <- p.adjust(P_DE_cor_df$p, method="none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
P_DE_cor_df_padj <- P_DE_cor_df[which(P_DE_cor_df$cor > 0.7| P_DE_cor_df$cor < -0.7),]
P_DE_cor_df_padj <- P_DE_cor_df_padj[which(P_DE_cor_df_padj$padj < 0.001),]

## Make node attribute table
P_nodeattrib_DE <- data.frame(node = union(P_DE_cor_df_padj$from,P_DE_cor_df_padj$to))
P_nodeattrib_DE$Phylum <- 0

for (i in as.character(P_nodeattrib_DE$node))
{
  if (i %in% rownames(tax_P) == TRUE)
  {P_nodeattrib_DE[P_nodeattrib_DE$node==i,"Phylum"] <- paste(tax_P[i,9:9])}
  else
  { P_nodeattrib_DE[P_nodeattrib_DE$node==i,"Phylum"]<- "NA"}
}


## Create co-occurrence network with igraph
P_DE_net <- graph_from_data_frame(P_DE_cor_df_padj, direct=F, vertices=P_nodeattrib_DE)

## Calculate relative abudnance of OTU nodes
P_DE_ra <- apply(cp_otu_P_DE,1,mean)
P_DE_ra <- P_DE_ra[V(P_DE_net)$name]

## Network properties ##

## Number of nodes in network
length(V(P_DE_net))

## Number of edges in network
length(E(P_DE_net))

## Connections 
bb_occur_PE <- droplevels(P_DE_cor_df_padj[with(P_DE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor>0,]
nrow(bb_occur_PE)

bb_occur_PE <- droplevels(P_DE_cor_df_padj[with(P_DE_cor_df_padj, grepl("bASV*",from) & grepl("bASV*",to)),])
bb_occur_PE <- bb_occur_PE[bb_occur_PE$cor<0,]
nrow(bb_occur_PE)

mean(E(P_DE_net))

graph.density(P_DE_net)

set.seed(619)
V(P_DE_net)$modularity <- membership(cluster_fast_greedy(P_DE_net))

P_DE_modularity <- V(P_DE_net)$modularity

P_DE_net_cfg <- cluster_fast_greedy(as.undirected(P_DE_net))
P_DE_net_modules <- sort(table(membership(P_DE_net_cfg)),decr=T)

transitivity(P_DE_net)

P_DE_net_betweenness <- betweenness(P_DE_net,normalized = T)
P_DE_net_degree <- igraph::degree(P_DE_net)

P_DE_net_modularity <- modularity(cluster_fast_greedy(P_DE_net))

plot(P_DE_net_betweenness,P_DE_net_degree)

P_DE_topological <- cbind(P_DE_net_betweenness,P_DE_net_degree,P_DE_modularity)

P_DE_topological_data <- as.data.frame(P_DE_topological)
colnames(P_DE_topological_data) <- c("betweenness","degree","modularity")

adjacency_P_DE <- get.adjacency(P_DE_net,attr = "cor",sparse = FALSE)

P_DE_zi_pi <- zi.pi(P_DE_topological_data, adjacency_P_DE, degree = 'degree', modularity_class = 'modularity')
head(P_DE_zi_pi)

P_DE_zi_pi <- na.omit(P_DE_zi_pi) 
P_DE_zi_pi[which(P_DE_zi_pi$within_module_connectivities < 2.5 & P_DE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
P_DE_zi_pi[which(P_DE_zi_pi$within_module_connectivities < 2.5 & P_DE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
P_DE_zi_pi[which(P_DE_zi_pi$within_module_connectivities > 2.5 & P_DE_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
P_DE_zi_pi[which(P_DE_zi_pi$within_module_connectivities > 2.5 & P_DE_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

P_DE_zi_pi_plot <- ggplot(P_DE_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'),
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)


## Calculate network node degrees/max/min
P_DE_deg <- sort(igraph::degree(P_DE_net,mode="all"),decr=T)
max(P_DE_deg)
mean(P_DE_deg)

## Set node colors based upon sensitivity to management system
unique(V(P_DE_net)$Phylum)
V(P_DE_net)$color <- V(P_DE_net)$Phylum

V(P_DE_net)$frame.color <- V(P_DE_net)$color

## Set node shape
V(P_DE_net)$shape <- "circle"

## Set node size
V(P_DE_net)$size <- 3



Type_igraph_character <- data.frame(c(length(V(B_CK_net)), length(E(B_CK_net)), mean(E(B_CK_net)),graph.density(B_CK_net), transitivity(B_CK_net), centralization.degree(B_CK_net)$centralization), 
                                    c(length(V(B_PE_net)), length(E(B_PE_net)), mean(E(B_PE_net)),graph.density(B_PE_net), transitivity(B_PE_net), centralization.degree(B_PE_net)$centralization),
                                    c(length(V(B_DE_net)), length(E(B_DE_net)), mean(E(B_DE_net)),graph.density(B_DE_net), transitivity(B_DE_net), centralization.degree(B_DE_net)$centralization),
                                    c(length(V(F_CK_net)), length(E(F_CK_net)), mean(E(F_CK_net)),graph.density(F_CK_net), transitivity(F_CK_net),centralization.degree(F_CK_net)$centralization), 
                                    c(length(V(F_PE_net)), length(E(F_PE_net)), mean(E(F_PE_net)),graph.density(F_PE_net), transitivity(F_PE_net),centralization.degree(F_PE_net)$centralization),
                                    c(length(V(F_DE_net)), length(E(F_DE_net)), mean(E(F_DE_net)),graph.density(F_DE_net), transitivity(F_DE_net),centralization.degree(F_DE_net)$centralization),
                                    c(length(V(P_CK_net)), length(E(P_CK_net)), mean(E(P_CK_net)),graph.density(P_CK_net), transitivity(P_CK_net),centralization.degree(P_CK_net)$centralization), 
                                    c(length(V(P_PE_net)), length(E(P_PE_net)), mean(E(P_PE_net)),graph.density(P_PE_net), transitivity(P_PE_net),centralization.degree(P_PE_net)$centralization),
                                    c(length(V(P_DE_net)), length(E(P_DE_net)), mean(E(P_DE_net)),graph.density(P_DE_net), transitivity(P_DE_net),centralization.degree(P_DE_net)$centralization))

names(Type_igraph_character) <- c("B_CK","B_PE","B_DE","F_CK","F_PE","F_DE","P_CK","P_PE","P_DE")
row.names(Type_igraph_character) <- c("node","edge","mean.edge","density","transitivity","modularity")

Type_igraph_character <- as.data.frame(t(Type_igraph_character))

Type_igraph_character$Kingdom <- c(rep("Bacteria",3),rep("Fungi",3),rep("Protist",3))
Type_igraph_character$type <- c(rep(c("CK","PE","DE"),3))
Type_igraph_character$type <- factor(Type_igraph_character$type,c("CK","PE","DE"))

B_sub_node <- ggplot(data = Type_igraph_character, mapping = aes(x = factor(Kingdom,levels = c("Bacteria","Fungi","Protist")), y = node, fill = type)) + scale_fill_manual(values=c("#2c7fbe","#8bc224","#f7d368"))+geom_bar(stat = 'identity', position = 'dodge')

B_sub_edge <- ggplot(data = Type_igraph_character, mapping = aes(x = factor(Kingdom,levels = c("Bacteria","Fungi","Protist")), y = edge, fill = type)) + scale_fill_manual(values=c("#2c7fbe","#8bc224","#f7d368"))+geom_bar(stat = 'identity', position = 'dodge')

B_sub_mean.edge <- ggplot(data = Type_igraph_character, mapping = aes(x = factor(Kingdom,levels = c("Bacteria","Fungi","Protist")), y = mean.edge, fill = type)) + scale_fill_manual(values=c("#2c7fbe","#8bc224","#f7d368"))+geom_bar(stat = 'identity', position = 'dodge')

B_sub_density <- ggplot(data = Type_igraph_character, mapping = aes(x = factor(Kingdom,levels = c("Bacteria","Fungi","Protist")), y =density , fill = type)) + scale_fill_manual(values=c("#2c7fbe","#8bc224","#f7d368"))+geom_bar(stat = 'identity', position = 'dodge')

B_sub_transitivity <- ggplot(data = Type_igraph_character, mapping = aes(x = factor(Kingdom,levels = c("Bacteria","Fungi","Protist")), y =transitivity  , fill = type)) + scale_fill_manual(values=c("#2c7fbe","#8bc224","#f7d368"))+geom_bar(stat = 'identity', position = 'dodge')

B_sub_modularity <- ggplot(data = Type_igraph_character, mapping = aes(x = factor(Kingdom,levels = c("Bacteria","Fungi","Protist")), y =modularity  , fill = type)) +scale_fill_manual(values=c("#2c7fbe","#8bc224","#f7d368"))+ geom_bar(stat = 'identity', position = 'dodge')

grid.newpage()
grid.arrange(B_sub_node,B_sub_edge,B_sub_mean.edge,B_sub_density,B_sub_transitivity,B_sub_modularity,ncol = 3)



## Note:the permutations of the network layouts is very time consuming and processor intensive
## Set layouts of the networks with Fruchterman & Reingold algorithim
set.seed(619)
coords_B_CK <- layout_(B_CK_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_B_CK,paste0("coords_B_CK.txt"),sep="\t",row.names=F,col.names=F,quote=F)
#
set.seed(619)
coords_B_PE <- layout_(B_PE_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_B_PE,paste0("coords_B_PE.txt"),sep="\t",row.names=F,col.names=F,quote=F)
#
set.seed(619)
coords_B_DE <- layout_(B_DE_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_B_DE,paste0("coords_B_DE.txt"),sep="\t",row.names=F,col.names=F,quote=F)

set.seed(619)
coords_F_CK <- layout_(F_CK_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_F_CK,paste0("coords_F_CK.txt"),sep="\t",row.names=F,col.names=F,quote=F)
#
set.seed(619)
coords_F_PE <- layout_(F_PE_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_F_PE,paste0("coords_F_PE.txt"),sep="\t",row.names=F,col.names=F,quote=F)
#
set.seed(619)
coords_F_DE <- layout_(F_DE_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_F_DE,paste0("coords_F_DE.txt"),sep="\t",row.names=F,col.names=F,quote=F)


set.seed(619)
coords_P_CK <- layout_(P_CK_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_P_CK,paste0("coords_P_CK.txt"),sep="\t",row.names=F,col.names=F,quote=F) 

set.seed(619)
coords_P_PE <- layout_(P_PE_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_P_PE,paste0("coords_P_PE.txt"),sep="\t",row.names=F,col.names=F,quote=F) 

set.seed(619)
coords_P_DE <- layout_(P_DE_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_P_DE,paste0("coords_P_DE.txt"),sep="\t",row.names=F,col.names=F,quote=F) 


## Import pre-calculated FR layout coordinates to save time 
coords_B_CK <- as.matrix(read.table("coords_B_CK.txt"))
dimnames(coords_B_CK) <-  NULL

coords_B_PE <- as.matrix(read.table("coords_B_PE.txt"))
dimnames(coords_B_PE) <-  NULL

coords_B_DE <- as.matrix(read.table("coords_B_DE.txt"))
dimnames(coords_B_DE) <-  NULL

coords_F_CK <- as.matrix(read.table("coords_F_CK.txt"))
dimnames(coords_F_CK) <-  NULL

coords_F_PE <- as.matrix(read.table("coords_F_PE.txt"))
dimnames(coords_F_PE) <-  NULL

coords_F_DE <- as.matrix(read.table("coords_F_DE.txt"))
dimnames(coords_F_DE) <-  NULL

coords_P_CK <- as.matrix(read.table("coords_P_CK.txt"))
dimnames(coords_P_CK) <-  NULL

coords_P_PE <- as.matrix(read.table("coords_P_PE.txt"))
dimnames(coords_P_PE) <-  NULL

coords_P_DE <- as.matrix(read.table("coords_P_DE.txt"))
dimnames(coords_P_DE) <-  NULL


## Plot Supplementary Figure S13
# pdf(paste0("network_2.pdf"),height=7,width=7)
par(mfrow=c(1,3),mar=c(0,0,0,0))
plot(B_CK_net, vertex.label=NA,layout=coords_B_CK)
plot(B_PE_net, vertex.label=NA,  layout=coords_B_PE)
plot(B_DE_net, vertex.label=NA, layout=coords_B_DE)
plot(F_CK_net, vertex.label=NA,layout=coords_F_CK)
plot(F_PE_net, vertex.label=NA,  layout=coords_F_PE)
plot(F_DE_net, vertex.label=NA, layout=coords_F_DE)

plot(P_CK_net, vertex.label=NA,layout=coords_P_CK)
plot(P_PE_net, vertex.label=NA,  layout=coords_P_PE)
plot(P_DE_net, vertex.label=NA, layout=coords_P_DE)

dev.off()


## Network Robustness
B_ps = phyloseq(sample_data(design_B),
                otu_table(as.matrix(cp_otu_B), taxa_are_rows=TRUE),
                tax_table(as.matrix(tax_B)))

B_res_rd = Robustness.Random.removal(ps = B_ps,
                                     r.threshold= 0.7,
                                     p.threshold=0.001,
                                     method = "spearman")

B_res_rd[[1]]

B_res_ta = Robustness.Targeted.removal(ps = B_ps,
                                       degree = TRUE,
                                       zipi = FALSE,
                                       r.threshold= 0.7,
                                       p.threshold=0.001,
                                       method = "spearman")
B_res_ta[[1]]

res = Vulnerability.micro(ps = B_ps,
                          Top = 500,
                          degree = TRUE,
                          zipi = FALSE,
                          r.threshold= 0.7,
                          p.threshold=0.001,
                          method = "spearman")

B_ps_Vul = res[[1]]
B_ps_Vul

tab.r = network.pip(
  ps = B_ps,
  N = 200,
  # ra = 0.05,
  big = FALSE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.7,
  p.threshold = 0.001,
  maxnode = 2,
  method = "spearman",
  label = FALSE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1
)

dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab

res6 = natural.con.microp (
  ps = B_ps,
  corg = cortab,
  norm = TRUE,
  end = 150,
  start = 0
)
B_con = res6[[1]]
B_con



## F
F_ps = phyloseq(sample_data(design_F),
                otu_table(as.matrix(cp_otu_F), taxa_are_rows=TRUE),
                tax_table(as.matrix(tax_F)))

F_res_rd = Robustness.Random.removal(ps = F_ps,
                                     r.threshold= 0.7,
                                     p.threshold=0.001,
                                     method = "spearman")

F_res_rd[[1]]

F_res_ta = Robustness.Targeted.removal(ps = F_ps,
                                       degree = TRUE,
                                       zipi = FALSE,
                                       r.threshold= 0.7,
                                       p.threshold=0.001,
                                       method = "spearman")
F_res_ta[[1]]

res = Vulnerability.micro(ps = F_ps,
                          Top = 500,
                          degree = TRUE,
                          zipi = FALSE,
                          r.threshold= 0.7,
                          p.threshold=0.001,
                          method = "spearman")

F_ps_Vul = res[[1]]
F_ps_Vul

tab.r = network.pip(
  ps = F_ps,
  N = 200,
  # ra = 0.05,
  big = FALSE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.7,
  p.threshold = 0.001,
  maxnode = 2,
  method = "spearman",
  label = FALSE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1
)

dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab

res6 = natural.con.microp (
  ps = F_ps,
  corg = cortab,
  norm = TRUE,
  end = 150,
  start = 0
)
F_con = res6[[1]]
F_con


## P_soil
P_ps = phyloseq(sample_data(design_P),
                otu_table(as.matrix(cp_otu_P), taxa_are_rows=TRUE),
                tax_table(as.matrix(tax_P)))

P_res_rd = Robustness.Random.removal(ps = P_ps,
                                     r.threshold= 0.7,
                                     p.threshold=0.001,
                                     method = "spearman")

P_res_rd[[1]]

P_res_ta = Robustness.Targeted.removal(ps = P_ps,
                                       degree = TRUE,
                                       zipi = FALSE,
                                       r.threshold= 0.8,
                                       p.threshold=0.05,
                                       method = "spearman")
P_res_ta[[1]]

res = Vulnerability.micro(ps = P_ps,
                          Top = 500,
                          degree = TRUE,
                          zipi = FALSE,
                          r.threshold= 0.7,
                          p.threshold=0.001,
                          method = "spearman")

P_ps_Vul = res[[1]]
P_ps_Vul

tab.r = network.pip(
  ps = P_ps,
  N = 200,
  # ra = 0.05,
  big = FALSE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.7,
  p.threshold = 0.001,
  maxnode = 2,
  method = "spearman",
  label = FALSE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1
)

dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab

res6 = natural.con.microp (
  ps = P_ps,
  corg = cortab,
  norm = TRUE,
  end = 150,
  start = 0
)
P_con = res6[[1]]
P_con


grid.newpage()
grid.arrange(B_res_rd[[1]],B_res_ta[[1]],B_ps_Vul,B_con,F_res_rd[[1]],F_res_ta[[1]],F_ps_Vul,F_con,P_res_rd[[1]],P_res_ta[[1]],P_ps_Vul,P_con,ncol = 4)

grid.newpage()
grid.arrange(B_ps_Vul,F_ps_Vul,P_ps_Vul,ncol = 3)

grid.newpage()
grid.arrange(B_con,F_con,P_con,ncol = 3)





##### netural model #####

## B_CK
cp_otu_B_CK <- cp_otu_B[,CKsamples]

keep_cp_otu_B_CK <- which(rowSums(cp_otu_B_CK ) >= 0)
cp_otu_B_CK <- cp_otu_B_CK[keep_cp_otu_B_CK,]

nrow(cp_otu_B_CK)

tax_B_CK <- tax_B[rownames(cp_otu_B_CK),]
design_B_CK <- droplevels(design_B[CKsamples,])


spp<-t(cp_otu_B_CK)

B_CK_niche_width <- niche.width(cp_otu_B_CK, method = 'levins')
B_CK_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 


draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)


## B_PE
cp_otu_B_PE <- cp_otu_B[,PEsamples]

keep_cp_otu_B_PE <- which(rowSums(cp_otu_B_PE ) >= 0)
cp_otu_B_PE <- cp_otu_B_PE[keep_cp_otu_B_PE,]

nrow(cp_otu_B_PE)

tax_B_PE <- tax_B[rownames(cp_otu_B_PE),]
design_B_PE <- droplevels(design_B[PEsamples,])


spp<-t(cp_otu_B_PE)

B_PE_niche_width <- niche.width(cp_otu_B_PE, method = 'levins')
B_PE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)



## B_DE
cp_otu_B_DE <- cp_otu_B[,DEsamples]

keep_cp_otu_B_DE <- which(rowSums(cp_otu_B_DE ) >= 0)
cp_otu_B_DE <- cp_otu_B_DE[keep_cp_otu_B_DE,]

nrow(cp_otu_B_DE)

tax_B_DE <- tax_B[rownames(cp_otu_B_DE),]
design_B_DE <- droplevels(design_B[DEsamples,])


spp<-t(cp_otu_B_DE)

B_DE_niche_width <- niche.width(cp_otu_B_DE, method = 'levins')
B_DE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)

B_niche_width <- data.frame(CK=t(B_CK_niche_width),PE=t(B_PE_niche_width),DE=t(B_DE_niche_width))

B_niche_width <- melt(B_niche_width)

B_width_box <- ggplot(B_niche_width, aes(x=factor(variable,levels = c("CK","PE","DE")), y=value, fill =variable))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  scale_fill_manual(values = c("#d2a128","#41ab6d","#297ec2"))


## F_CK
cp_otu_F_CK <- cp_otu_F[,CKsamples]

keep_cp_otu_F_CK <- which(rowSums(cp_otu_F_CK ) >= 0)
cp_otu_F_CK <- cp_otu_F_CK[keep_cp_otu_F_CK,]

nrow(cp_otu_F_CK)

tax_F_CK <- tax_F[rownames(cp_otu_F_CK),]
design_F_CK <- droplevels(design_F[CKsamples,])


spp<-t(cp_otu_F_CK)

F_CK_niche_width <- niche.width(cp_otu_F_CK, method = 'levins')
F_CK_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)


## F_PE
cp_otu_F_PE <- cp_otu_F[,PEsamples]

keep_cp_otu_F_PE <- which(rowSums(cp_otu_F_PE ) >= 0)
cp_otu_F_PE <- cp_otu_F_PE[keep_cp_otu_F_PE,]

nrow(cp_otu_F_PE)

tax_F_PE <- tax_F[rownames(cp_otu_F_PE),]
design_F_PE <- droplevels(design_F[PEsamples,])


spp<-t(cp_otu_F_PE)

F_PE_niche_width <- niche.width(cp_otu_F_PE, method = 'levins')
F_PE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)



## F_DE
cp_otu_F_DE <- cp_otu_F[,DEsamples]

keep_cp_otu_F_DE <- which(rowSums(cp_otu_F_DE ) >= 0)
cp_otu_F_DE <- cp_otu_F_DE[keep_cp_otu_F_DE,]

nrow(cp_otu_F_DE)

tax_F_DE <- tax_F[rownames(cp_otu_F_DE),]
design_F_DE <- droplevels(design_F[DEsamples,])


spp<-t(cp_otu_F_DE)

F_DE_niche_width <- niche.width(cp_otu_F_DE, method = 'levins')
F_DE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)


## P_CK
cp_otu_P_CK <- cp_otu_P[,pCKsamples]

keep_cp_otu_P_CK <- which(rowSums(cp_otu_P_CK ) >= 0)
cp_otu_P_CK <- cp_otu_P_CK[keep_cp_otu_P_CK,]

nrow(cp_otu_P_CK)

tax_P_CK <- tax_P[rownames(cp_otu_P_CK),]
design_P_CK <- droplevels(design_P[pCKsamples,])


spp<-t(cp_otu_P_CK)

P_CK_niche_width <- niche.width(cp_otu_P_CK, method = 'levins')
P_CK_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)


## P_PE
cp_otu_P_PE <- cp_otu_P[,pPEsamples]

keep_cp_otu_P_PE <- which(rowSums(cp_otu_P_PE ) >= 0)
cp_otu_P_PE <- cp_otu_P_PE[keep_cp_otu_P_PE,]

nrow(cp_otu_P_PE)

tax_P_PE <- tax_P[rownames(cp_otu_P_PE),]
design_P_PE <- droplevels(design_P[pPEsamples,])


spp<-t(cp_otu_P_PE)

P_PE_niche_width <- niche.width(cp_otu_P_PE, method = 'levins')
P_PE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)



## P_DE
cp_otu_P_DE <- cp_otu_P[,pDEsamples]

keep_cp_otu_P_DE <- which(rowSums(cp_otu_P_DE ) >= 0)
cp_otu_P_DE <- cp_otu_P_DE[keep_cp_otu_P_DE,]

nrow(cp_otu_P_DE)

tax_P_DE <- tax_P[rownames(cp_otu_P_DE),]
design_P_DE <- droplevels(design_P[pDEsamples,])


spp<-t(cp_otu_P_DE)

P_DE_niche_width <- niche.width(cp_otu_P_DE, method = 'levins')
P_DE_niche_width



N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'


library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 





draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)



#####  all network #####

B_lefse_tax <- read.table("B_cladogram_ASV.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)

F_lefse_tax <- read.table("F_cladogram_ASV.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)

P_lefse_tax <- read.table("P_ASV.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)


## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
keep_cp_otu_B <- which(rowSums(cp_otu_B >= 1) >= 4)
cp_otu_B <- cp_otu_B[keep_cp_otu_B,]

nrow(cp_otu_B)

tax_B <- tax_B[rownames(cp_otu_B),]

B_occor <- rcorr(t(cp_otu_B),type=c("spearman"))

## Create data frame of co-occurring OTUs
B_cor_df <- CorrDF(B_occor$r,B_occor$P)
B_cor_df$padj <- p.adjust(B_cor_df$p, method = "none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
B_cor_df_padj <- B_cor_df[which(B_cor_df$cor > 0.7| B_cor_df$cor < -0.7),]
B_cor_df_padj <- B_cor_df_padj[which(B_cor_df_padj$padj < 0.001),]

B_nodeattrib <- data.frame(node = union(B_cor_df_padj$from,B_cor_df_padj$to))

for (i in as.character(B_nodeattrib$node))
{
  if (i %in% rownames(B_lefse_tax) == TRUE)
  {B_nodeattrib[B_nodeattrib$node==i,"Phylum"] <- paste(B_lefse_tax[i,4:4])}
  else
  { B_nodeattrib[B_nodeattrib$node==i,"Phylum"]<- "NA"}
}

## Create co-occurrence network with igraph
B_net <- graph_from_data_frame(B_cor_df_padj,direct=F, vertices = B_nodeattrib)

## Calculate relative abundance of OTU nodes
B_ra <- apply(cp_otu_B,1,mean)
B_ra <- B_ra[V(B_net)$name]

## Network properties ##

## Number of nodes in network
length(V(B_net))

## Number of edges in network
length(E(B_net))

mean(E(B_net))

graph.density(B_net)

B_net_cfg <- cluster_fast_greedy(as.undirected(B_net))
B_net_modules <- sort(table(membership(B_net_cfg)),decr=T)

B_modules_20 <- B_net_modules[1:20]
sum(B_modules_20)/sum(B_net_modules)
rm20_plot <- B_modules_20
names(rm20_plot) <- as.factor(1:20)

transitivity(B_net)

B_net_betweenness <- betweenness(B_net,normalized = T)
B_net_degree <- igraph::degree(B_net)
B_net_closeness <- closeness(B_net)

## Set node colors based upon sensitivity to management system
cs <- c("PE","DE","CK")
unique(V(B_net)$Phylum)
V(B_net)$color <- V(B_net)$Phylum
V(B_net)$color[!V(B_net)$color %in% cs] <- "gray30"
V(B_net)$color[V(B_net)$color == "PE"] <- "#f7d368"
V(B_net)$color[V(B_net)$color == "CK"] <- "#2c7fbe"
V(B_net)$color[V(B_net)$color == "DE"] <- "#8bc224"



## Set vector of nodes responding to certain management systems
rownames(B_nodeattrib) <- B_nodeattrib$node
rootB_nodes <- rownames(B_nodeattrib[B_nodeattrib$Phylum %in% cs,])

## Set node shape
V(B_net)$shape <- "circle"

## Set node size
V(B_net)$size <- V(B_net)$name
V(B_net)$size[!V(B_net)$size %in% rootB_nodes] <- 3
V(B_net)$size[V(B_net)$size %in% rootB_nodes] <- 6
rootB_nodesizes <- as.numeric(V(B_net)$size)


CK_nodes_B <- rownames(B_nodeattrib[B_nodeattrib$Phylum=="CK",])
PE_nodes_B <- rownames(B_nodeattrib[B_nodeattrib$Phylum=="PE",])
DE_nodes_B <- rownames(B_nodeattrib[B_nodeattrib$Phylum=="DE",])

cs_nodes_all <- c(CK_nodes_B,PE_nodes_B,DE_nodes_B)

sum(B_modules_20)/sum(B_net_modules)
sm20_plot <- B_modules_20
names(sm20_plot) <- as.factor(1:20)
B_modules_cs  <- table(factor(membership(B_net_cfg)[cs_nodes_all],levels=names(B_net_modules)))
B_modules_cs_20 <- B_modules_cs[names(B_modules_20)]
rmcs20_plot <- B_modules_cs_20
names(rmcs20_plot) <- as.factor(1:20)

## Make vector of nodes in top 20 modules
B_modules_points <- membership(B_net_cfg)[membership(B_net_cfg) %in% names(B_modules_20)]

B_points <- NULL
for(i in B_modules_points){
  B_rootx <- which(names(B_modules_20)==i)
  B_points <- c(B_points,B_rootx)
}
names(B_points) <- names(B_modules_points)

B_mods_list_cs <- list()
for (i in names(B_modules_cs_20)){
  x1 <- names(membership(B_net_cfg)[membership(B_net_cfg) == i])
  x2 <- x1[x1 %in% cs_nodes_all]
  B_mods_list_cs[[i]] <- as.numeric(V(B_net)[x2])
}


## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
keep_cp_otu_F <- which(rowSums(cp_otu_F >= 1) >= 4)
cp_otu_F <- cp_otu_F[keep_cp_otu_F,]

nrow(cp_otu_F)

tax_F <- tax_F[rownames(cp_otu_F),]

F_occor <- rcorr(t(cp_otu_F),type=c("spearman"))

## Create data frame of co-occurring OTUs
F_cor_df <- CorrDF(F_occor$r,F_occor$P)
F_cor_df$padj <- p.adjust(F_cor_df$p, method = "none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
F_cor_df_padj <- F_cor_df[which(F_cor_df$cor > 0.7| F_cor_df$cor < -0.7),]
F_cor_df_padj <- F_cor_df_padj[which(F_cor_df_padj$padj < 0.001),]

F_nodeattrib <- data.frame(node = union(F_cor_df_padj$from,F_cor_df_padj$to))

for (i in as.character(F_nodeattrib$node))
{
  if (i %in% rownames(F_lefse_tax) == TRUE)
  {F_nodeattrib[F_nodeattrib$node==i,"Phylum"] <- paste(F_lefse_tax[i,4:4])}
  else
  { F_nodeattrib[F_nodeattrib$node==i,"Phylum"]<- "NA"}
}

## Create co-occurrence network with igraph
F_net <- graph_from_data_frame(F_cor_df_padj,direct=F, vertices = F_nodeattrib)

## Calculate relative abundance of OTU nodes
F_ra <- apply(cp_otu_F,1,mean)
F_ra <- F_ra[V(F_net)$name]

## Network properties ##

## Number of nodes in network
length(V(F_net))

## Number of edges in network
length(E(F_net))

mean(E(F_net))

graph.density(F_net)

F_net_cfg <- cluster_fast_greedy(as.undirected(F_net))
F_net_modules <- sort(table(membership(F_net_cfg)),decr=T)

F_modules_20 <- F_net_modules[1:20]
sum(F_modules_20)/sum(F_net_modules)
rm20_plot <- F_modules_20
names(rm20_plot) <- as.factor(1:20)

transitivity(F_net)

F_net_betweenness <- betweenness(F_net,normalized = T)
F_net_degree <- igraph::degree(F_net)
F_net_closeness <- closeness(F_net)

## Set node colors based upon sensitivity to management system
cs <- c("PE","DE","CK")
unique(V(F_net)$Phylum)
V(F_net)$color <- V(F_net)$Phylum
V(F_net)$color[!V(F_net)$color %in% cs] <- "gray30"
V(F_net)$color[V(F_net)$color == "PE"] <- "#f7d368"
V(F_net)$color[V(F_net)$color == "CK"] <- "#2c7fbe"
V(F_net)$color[V(F_net)$color == "DE"] <- "#8bc224"



## Set vector of nodes responding to certain management systems
rownames(F_nodeattrib) <- F_nodeattrib$node
rootF_nodes <- rownames(F_nodeattrib[F_nodeattrib$Phylum %in% cs,])

## Set node shape
V(F_net)$shape <- "circle"

## Set node size
V(F_net)$size <- V(F_net)$name
V(F_net)$size[!V(F_net)$size %in% rootF_nodes] <- 3
V(F_net)$size[V(F_net)$size %in% rootF_nodes] <- 6
rootF_nodesizes <- as.numeric(V(F_net)$size)


CK_nodes_F <- rownames(F_nodeattrib[F_nodeattrib$Phylum=="CK",])
PE_nodes_F <- rownames(F_nodeattrib[F_nodeattrib$Phylum=="PE",])
DE_nodes_F <- rownames(F_nodeattrib[F_nodeattrib$Phylum=="DE",])

cs_nodes_all <- c(CK_nodes_F,PE_nodes_F,DE_nodes_F)

sum(F_modules_20)/sum(F_net_modules)
sm20_plot <- F_modules_20
names(sm20_plot) <- as.factor(1:20)
F_modules_cs  <- table(factor(membership(F_net_cfg)[cs_nodes_all],levels=names(F_net_modules)))
F_modules_cs_20 <- F_modules_cs[names(F_modules_20)]
rmcs20_plot <- F_modules_cs_20
names(rmcs20_plot) <- as.factor(1:20)

## Make vector of nodes in top 20 modules
F_modules_points <- membership(F_net_cfg)[membership(F_net_cfg) %in% names(F_modules_20)]

F_points <- NULL
for(i in F_modules_points){
  F_rootx <- which(names(F_modules_20)==i)
  F_points <- c(F_points,F_rootx)
}
names(F_points) <- names(F_modules_points)

F_mods_list_cs <- list()
for (i in names(F_modules_cs_20)){
  x1 <- names(membership(F_net_cfg)[membership(F_net_cfg) == i])
  x2 <- x1[x1 %in% cs_nodes_all]
  F_mods_list_cs[[i]] <- as.numeric(V(F_net)[x2])
}


## Create bulk soil bacteria co-occurrence network with igraph package
## Perform pairwise Spearman correlations on bulk soil bacteria community on TMM normalized CPM counts
keep_cp_otu_P <- which(rowSums(cp_otu_P >= 1) >= 4)
cp_otu_P <- cp_otu_P[keep_cp_otu_P,]

nrow(cp_otu_P)

tax_P <- tax_P[rownames(cp_otu_P),]

P_occor <- rcorr(t(cp_otu_P),type=c("spearman"))

## Create data frame of co-occurring OTUs
P_cor_df <- CorrDF(P_occor$r,P_occor$P)
P_cor_df$padj <- p.adjust(P_cor_df$p, method = "none")

## Subset data frame for co-occurring OTUs with Spearman's rho > 0.5 and a p-value < 0.001
P_cor_df_padj <- P_cor_df[which(P_cor_df$cor > 0.7| P_cor_df$cor < -0.7),]
P_cor_df_padj <- P_cor_df_padj[which(P_cor_df_padj$padj < 0.001),]

P_nodeattrib <- data.frame(node = union(P_cor_df_padj$from,P_cor_df_padj$to))

for (i in as.character(P_nodeattrib$node))
{
  if (i %in% rownames(P_lefse_tax) == TRUE)
  {P_nodeattrib[P_nodeattrib$node==i,"Phylum"] <- paste(P_lefse_tax[i,4:4])}
  else
  { P_nodeattrib[P_nodeattrib$node==i,"Phylum"]<- "NA"}
}

## Create co-occurrence network with igraph
P_net <- graph_from_data_frame(P_cor_df_padj,direct=F, vertices = P_nodeattrib)

## Calculate relative abundance of OTU nodes
P_ra <- apply(cp_otu_P,1,mean)
P_ra <- P_ra[V(P_net)$name]

## Network properties ##

## Number of nodes in network
length(V(P_net))

## Number of edges in network
length(E(P_net))

mean(E(P_net))

graph.density(P_net)

P_net_cfg <- cluster_fast_greedy(as.undirected(P_net))
P_net_modules <- sort(table(membership(P_net_cfg)),decr=T)

P_modules_20 <- P_net_modules[1:20]
sum(P_modules_20)/sum(P_net_modules)
rm20_plot <- P_modules_20
names(rm20_plot) <- as.factor(1:20)

transitivity(P_net)

P_net_betweenness <- betweenness(P_net,normalized = T)
P_net_degree <- igraph::degree(P_net)
P_net_closeness <- closeness(P_net)

## Set node colors based upon sensitivity to management system
cs <- c("PE","DE","CK")
unique(V(P_net)$Phylum)
V(P_net)$color <- V(P_net)$Phylum
V(P_net)$color[!V(P_net)$color %in% cs] <- "gray30"
V(P_net)$color[V(P_net)$color == "PE"] <- "#f7d368"
V(P_net)$color[V(P_net)$color == "CK"] <- "#2c7fbe"
V(P_net)$color[V(P_net)$color == "DE"] <- "#8bc224"



## Set vector of nodes responding to certain management systems
rownames(P_nodeattrib) <- P_nodeattrib$node
rootP_nodes <- rownames(P_nodeattrib[P_nodeattrib$Phylum %in% cs,])

## Set node shape
V(P_net)$shape <- "circle"

## Set node size
V(P_net)$size <- V(P_net)$name
V(P_net)$size[!V(P_net)$size %in% rootP_nodes] <- 3
V(P_net)$size[V(P_net)$size %in% rootP_nodes] <- 6
rootP_nodesizes <- as.numeric(V(P_net)$size)


CK_nodes_P <- rownames(P_nodeattrib[P_nodeattrib$Phylum=="CK",])
PE_nodes_P <- rownames(P_nodeattrib[P_nodeattrib$Phylum=="PE",])
DE_nodes_P <- rownames(P_nodeattrib[P_nodeattrib$Phylum=="DE",])

cs_nodes_all <- c(CK_nodes_P,PE_nodes_P,DE_nodes_P)

sum(P_modules_20)/sum(P_net_modules)
sm20_plot <- P_modules_20
names(sm20_plot) <- as.factor(1:20)
P_modules_cs  <- table(factor(membership(P_net_cfg)[cs_nodes_all],levels=names(P_net_modules)))
P_modules_cs_20 <- P_modules_cs[names(P_modules_20)]
rmcs20_plot <- P_modules_cs_20
names(rmcs20_plot) <- as.factor(1:20)

## Make vector of nodes in top 20 modules
P_modules_points <- membership(P_net_cfg)[membership(P_net_cfg) %in% names(P_modules_20)]

P_points <- NULL
for(i in P_modules_points){
  P_rootx <- which(names(P_modules_20)==i)
  P_points <- c(P_points,P_rootx)
}
names(P_points) <- names(P_modules_points)

P_mods_list_cs <- list()
for (i in names(P_modules_cs_20)){
  x1 <- names(membership(P_net_cfg)[membership(P_net_cfg) == i])
  x2 <- x1[x1 %in% cs_nodes_all]
  P_mods_list_cs[[i]] <- as.numeric(V(P_net)[x2])
}



#####  Figure S13 : Individual Co-occurence networks #####

## Note:the permutations of the network layouts is very time consuming and processor intensive
## Set layouts of the networks with Fruchterman & Reingold algorithim
set.seed(619)
coords_B <- layout_(B_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_B,paste0("coords_B.txt"),sep="\t",row.names=F,col.names=F,quote=F)

set.seed(619)
coords_F <- layout_(F_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_F,paste0("coords_F.txt"),sep="\t",row.names=F,col.names=F,quote=F)

set.seed(619)
coords_P <- layout_(P_net,with_fr(niter=9999, grid="nogrid"))
write.table(coords_P,paste0("coords_P.txt"),sep="\t",row.names=F,col.names=F,quote=F)

## Import pre-calculated FR layout coordinates to save time 
coords_B <- as.matrix(read.table("coords_B.txt"))
dimnames(coords_B) <-  NULL

coords_F <- as.matrix(read.table("coords_F.txt"))
dimnames(coords_F) <-  NULL

coords_P <- as.matrix(read.table("coords_P.txt"))
dimnames(coords_P) <-  NULL

## Plot Supplementary Figure S13
# pdf(paste0("network_2.pdf"),height=7,width=7)
par(mfrow=c(1,3),mar=c(0,0,0,0))

plot(B_net, vertex.label=NA, vertex.size= rootB_nodesizes, layout=coords_B)
plot(F_net, vertex.label=NA, vertex.size= rootF_nodesizes, layout=coords_F)
plot(P_net, vertex.label=NA, vertex.size= rootP_nodesizes, layout=coords_P)

dev.off()




## N cycle phylum
## N gene
Kegg_N <- read.table("KO_N.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")
Kegg_N_tax <- read.table("N_tax_gene.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")

edgeR_kegg_N <- DGEList(counts=Kegg_N, 
                        group=metadata1$Treatment1)


## N cycle Phylum
KO_N <- names(sort(table(Kegg_N_tax[,"Phylum"]), decr=T))
length## N cycle
KO_N <- names(sort(table(Kegg_N_tax[,"Phylum"]), decr=T))
length(KO_N)
sort(table(Kegg_N_tax[,"Phylum"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_N)
for (i in KO_N){
  x <- array(colSums(Kegg_N[rownames(Kegg_N_tax)[which(Kegg_N_tax$Phylum == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(KO_N)
colnames(y) <- paste(colnames(Kegg_N))
KO_N_sum <- y


##barplot 
metadata1$Treatment1 <- factor(metadata1$Treatment1,c("CK","PE","DE"))

CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])

PHYLUM_mat_N_fix <-as.matrix(cbind(`CK`=apply(KO_N_sum[,CK],1,sum),
                                   `PE`=apply(KO_N_sum[,PE],1,sum),
                                   `DE`=apply(KO_N_sum[,DE],1,sum)))

PHYLUM_mat_16s_mean <- sort(apply(PHYLUM_mat_N_fix,1,mean),decr=T)
PHYLUM_mat_N_fix <- PHYLUM_mat_N_fix[names(PHYLUM_mat_16s_mean),]

PHYLUM_mat_N_fix <- PHYLUM_mat_N_fix[-c(4,9),]
PHYLUM_mat_N_fix <- PHYLUM_mat_N_fix[c(1:9),]


PHYLUM_mat_N_fix <- PHYLUM_mat_N_fix*1000

PHYLUM_mat_N_fix<-melt (PHYLUM_mat_N_fix)

a <- ggplot(PHYLUM_mat_N_fix, aes(x = Var2, y = Var1, size = value, color=Var2)) +  geom_point()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray"),
        panel.border = element_rect(colour="black",fill=NA))


## N gene
Kegg_N <- read.table("KO_N.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")
Kegg_N_tax <- read.table("N_tax_gene.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F,quote = "")

edgeR_kegg_N <- DGEList(counts=Kegg_N, 
                        group=metadata1$Treatment1)

## N cycle Genus
KO_N <- names(sort(table(Kegg_N_tax[,"Genus"]), decr=T))
length## N cycle
KO_N <- names(sort(table(Kegg_N_tax[,"Genus"]), decr=T))
length(KO_N)
sort(table(Kegg_N_tax[,"Genus"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_N)
for (i in KO_N){
  x <- array(colSums(Kegg_N[rownames(Kegg_N_tax)[which(Kegg_N_tax$Genus == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(KO_N)
colnames(y) <- paste(colnames(Kegg_N))
KO_N_sum <- y


##barplot 
metadata1$Treatment1 <- factor(metadata1$Treatment1,c("CK","PE","DE"))

CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])

PHYLUM_mat_N_fix <-as.matrix(cbind(`CK`=apply(KO_N_sum[,CK],1,sum),
                                   `PE`=apply(KO_N_sum[,PE],1,sum),
                                   `DE`=apply(KO_N_sum[,DE],1,sum)))

PHYLUM_mat_16s_mean <- sort(apply(PHYLUM_mat_N_fix,1,mean),decr=T)
PHYLUM_mat_N_fix <- PHYLUM_mat_N_fix[names(PHYLUM_mat_16s_mean),]


PHYLUM_mat_N_fix <- PHYLUM_mat_N_fix[-c(1,2),]
PHYLUM_mat_N_fix <- PHYLUM_mat_N_fix[c(1:20),]


PHYLUM_mat_N_fix <- PHYLUM_mat_N_fix*1000

PHYLUM_mat_N_fix<-melt (PHYLUM_mat_N_fix)

b <- ggplot(PHYLUM_mat_N_fix, aes(x = Var2, y = Var1, size = value, color=Var2)) + geom_point()+scale_fill_manual(values = c("#2c7fbe","#f7d368","#8bc224"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray"),
        panel.border = element_rect(colour="black",fill=NA))

grid.newpage()
grid.arrange(a ,b ,ncol = 2)




## N cycle process 

r1 <- rownames(Kegg_N_tax[Kegg_N_tax$Phylum == "unclassified",])
r2 <- rownames(Kegg_N_tax[Kegg_N_tax$Phylum == "unknown",])

otus_remove <- c(r1,r2)

Kegg_N <- Kegg_N[-which(rownames(Kegg_N) %in% otus_remove),]

Kegg_N_tax <- Kegg_N_tax[rownames(Kegg_N),]


## N fixation
fixation <- c("nifH","nifK","norC","nosZ")
Kegg_N_fix_tax  <-  Kegg_N_tax %>% filter(Function %in% fixation)
Kegg_N_fix_tab <- Kegg_N[rownames(Kegg_N_fix_tax),]

metadata1$Treatment1 <- factor(metadata1$Treatment1,c("CK","PE","DE"))

CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])

Kegg_N_fix_tab <-as.matrix(cbind(`CK`=apply(Kegg_N_fix_tab[,CK],1,sum),
                                 `PE`=apply(Kegg_N_fix_tab[,PE],1,sum),
                                 `DE`=apply(Kegg_N_fix_tab[,DE],1,sum)))


Kegg_N_fix_tab <- t(t(Kegg_N_fix_tab)/colSums(Kegg_N_fix_tab)) * 100

KO_N_fix <- names(sort(table(Kegg_N_fix_tax[,"Phylum"]), decr=T))

length(KO_N_fix)
sort(table(Kegg_N_fix_tax[,"Phylum"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_N_fix_tab)
for (i in KO_N_fix){
  x <- array(colSums(Kegg_N_fix_tab[rownames(Kegg_N_fix_tax)[which(Kegg_N_fix_tax$Phylum == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(KO_N_fix)
colnames(y) <- paste(colnames(Kegg_N_fix_tab))
KO_N_sum_fix <- y


# Phyla with MEAN abundances higher than 1% relative abundances
Kegg_N_fix_tax$labels <- Kegg_N_fix_tax$Phylum
Kegg_N_fix_tax$cols <- Kegg_N_fix_tax$Phylum

Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="Acidobacteria" ], ]$cols <- "palegreen4"
Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="Bacteroidetes" ], ]$cols <- "sandybrown"
Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="Nitrospirae" ], ]$cols <- "gray"
Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="Ignavibacteriae" ], ]$cols <- "gray"
Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="CandidatusKapabacteria" ], ]$cols <- "orchid"
Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="Gemmatimonadetes" ], ]$cols <- "steelblue4"
Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="Verrucomicrobia" ], ]$cols <- "steelblue1"
Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="Proteobacteria" ], ]$cols <- "palegreen3"
Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="Calditrichaeota" ], ]$cols <- "gray"
Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="CandidatusBlackallbacteria" ], ]$cols <- "gray"
Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="CandidatusRokubacteria" ], ]$cols <- "gray"
Kegg_N_fix_tax[ rownames(Kegg_N_fix_tax)[Kegg_N_fix_tax$labels=="Calditrichaeota" ], ]$cols <- "gray"


## collaps OTU colors to prepare Phylum level colors
label_cols_B <- Kegg_N_fix_tax[, c("labels", "cols") ]
library(plyr)
PHYLA_label_cols_B <- ddply(label_cols_B, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_B) <- PHYLA_label_cols_B[,1]
PHYLA_label_cols_B


##### Plot Supplementary Figure S3
par(mfrow=c(1,5), mar=c(6,4,1,5), xpd=NA)
phylum_bar_B <- barplot(as.matrix(KO_N_sum_fix), col=PHYLA_label_cols_B[rownames(KO_N_sum_fix),]$cols, ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_B, labels=colnames(PHYLUM_mat_B_mean_type), col.axis="black", las=2, cex.axis=1)
title(ylab="Relative abundance (%)")
title(main="Bacteria Community")
legend(3, 100, bty="n", cex=0.7, x.intersp= 0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_B$labels), 
       fill=rev(PHYLA_label_cols_B$cols), 
       border=rev(PHYLA_label_cols_B$cols) )



## N Nitrification
Nitrification <- c("pmoA-amoA","pmoB-amoB","pmoC-amoC","hao")
Kegg_N_nit_tax  <-  Kegg_N_tax %>% filter(Function %in% Nitrification)
Kegg_N_nit_tab <- Kegg_N[rownames(Kegg_N_nit_tax),]

metadata1$Treatment1 <- factor(metadata1$Treatment1,c("CK","PE","DE"))

CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])

Kegg_N_nit_tab <-as.matrix(cbind(`CK`=apply(Kegg_N_nit_tab[,CK],1,sum),
                                 `PE`=apply(Kegg_N_nit_tab[,PE],1,sum),
                                 `DE`=apply(Kegg_N_nit_tab[,DE],1,sum)))


Kegg_N_nit_tab <- t(t(Kegg_N_nit_tab)/colSums(Kegg_N_nit_tab)) * 100

KO_N_nit <- names(sort(table(Kegg_N_nit_tax[,"Phylum"]), decr=T))

length(KO_N_nit)
sort(table(Kegg_N_nit_tax[,"Phylum"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_N_nit_tab)
for (i in KO_N_nit){
  x <- array(colSums(Kegg_N_nit_tab[rownames(Kegg_N_nit_tax)[which(Kegg_N_nit_tax$Phylum == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(KO_N_nit)
colnames(y) <- paste(colnames(Kegg_N_nit_tab))
KO_N_sum_nit <- y


# Phyla with MEAN abundances higher than 1% relative abundances
Kegg_N_nit_tax$labels <- Kegg_N_nit_tax$Phylum
Kegg_N_nit_tax$cols <- Kegg_N_nit_tax$Phylum


Kegg_N_nit_tax[ rownames(Kegg_N_nit_tax)[Kegg_N_nit_tax$labels=="Actinobacteria" ], ]$cols <- "dodgerblue"
Kegg_N_nit_tax[ rownames(Kegg_N_nit_tax)[Kegg_N_nit_tax$labels=="Nitrospirae" ], ]$cols <- "gray"
Kegg_N_nit_tax[ rownames(Kegg_N_nit_tax)[Kegg_N_nit_tax$labels=="Proteobacteria" ], ]$cols <- "palegreen3"
Kegg_N_nit_tax[ rownames(Kegg_N_nit_tax)[Kegg_N_nit_tax$labels=="Crenarchaeota" ], ]$cols <- "gray"
Kegg_N_nit_tax[ rownames(Kegg_N_nit_tax)[Kegg_N_nit_tax$labels=="Thaumarchaeota" ], ]$cols <- "gray"


## collaps OTU colors to prepare Phylum level colors
label_cols_B <- Kegg_N_nit_tax[, c("labels", "cols") ]
library(plyr)
PHYLA_label_cols_B <- ddply(label_cols_B, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_B) <- PHYLA_label_cols_B[,1]
PHYLA_label_cols_B

##### Plot Supplementary Figure S3
phylum_bar_B <- barplot(as.matrix(KO_N_sum_nit), col=PHYLA_label_cols_B[rownames(KO_N_sum_nit),]$cols, ylim=c(0,100), xaxt="n", border=NA, las=2)

legend(3, 100, bty="n", cex=0.7, x.intersp= 0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_B$labels), 
       fill=rev(PHYLA_label_cols_B$cols), 
       border=rev(PHYLA_label_cols_B$cols) )



## N Denitrification
Denitrification <- c("narB","nirA","nirB","nirS")
Kegg_N_den_tax  <-  Kegg_N_tax %>% filter(Function %in% Denitrification)
Kegg_N_den_tab <- Kegg_N[rownames(Kegg_N_den_tax),]

metadata1$Treatment1 <- factor(metadata1$Treatment1,c("CK","PE","DE"))

CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])

Kegg_N_den_tab <-as.matrix(cbind(`CK`=apply(Kegg_N_den_tab[,CK],1,sum),
                                 `PE`=apply(Kegg_N_den_tab[,PE],1,sum),
                                 `DE`=apply(Kegg_N_den_tab[,DE],1,sum)))


Kegg_N_den_tab <- t(t(Kegg_N_den_tab)/colSums(Kegg_N_den_tab)) * 100

KO_N_den <- names(sort(table(Kegg_N_den_tax[,"Phylum"]), decr=T))

length(KO_N_den)
sort(table(Kegg_N_den_tax[,"Phylum"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_N_den_tab)
for (i in KO_N_den){
  x <- array(colSums(Kegg_N_den_tab[rownames(Kegg_N_den_tax)[which(Kegg_N_den_tax$Phylum == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(KO_N_den)
colnames(y) <- paste(colnames(Kegg_N_den_tab))
KO_N_sum_den <- y


# Phyla with MEAN abundances higher than 1% relative abundances
Kegg_N_den_tax$labels <- Kegg_N_den_tax$Phylum
Kegg_N_den_tax$cols <- Kegg_N_den_tax$Phylum


Kegg_N_den_tax[ rownames(Kegg_N_den_tax)[Kegg_N_den_tax$labels=="Actinobacteria" ], ]$cols <- "dodgerblue"
Kegg_N_den_tax[ rownames(Kegg_N_den_tax)[Kegg_N_den_tax$labels=="Verrucomicrobia" ], ]$cols <- "gray"
Kegg_N_den_tax[ rownames(Kegg_N_den_tax)[Kegg_N_den_tax$labels=="Proteobacteria" ], ]$cols <- "palegreen3"
Kegg_N_den_tax[ rownames(Kegg_N_den_tax)[Kegg_N_den_tax$labels=="Acidobacteria" ], ]$cols <- "palegreen4"
Kegg_N_den_tax[ rownames(Kegg_N_den_tax)[Kegg_N_den_tax$labels=="Bacteroidetes" ], ]$cols <- "sandybrown"
Kegg_N_den_tax[ rownames(Kegg_N_den_tax)[Kegg_N_den_tax$labels=="Planctomycetes" ], ]$cols <- "gray"
Kegg_N_den_tax[ rownames(Kegg_N_den_tax)[Kegg_N_den_tax$labels=="Bacteroidetes" ], ]$cols <- "sandybrown"
Kegg_N_den_tax[ rownames(Kegg_N_den_tax)[Kegg_N_den_tax$labels=="Cyanobacteria" ], ]$cols <- "gray"
Kegg_N_den_tax[ rownames(Kegg_N_den_tax)[Kegg_N_den_tax$labels=="Chloroflexi" ], ]$cols <- "plum1"



## collaps OTU colors to prepare Phylum level colors
label_cols_B <- Kegg_N_den_tax[, c("labels", "cols") ]
library(plyr)
PHYLA_label_cols_B <- ddply(label_cols_B, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_B) <- PHYLA_label_cols_B[,1]
PHYLA_label_cols_B

##### Plot Supplementary Figure S3
phylum_bar_B <- barplot(as.matrix(KO_N_sum_den), col=PHYLA_label_cols_B[rownames(KO_N_sum_den),]$cols, ylim=c(0,100), xaxt="n", border=NA, las=2)

legend(3, 100, bty="n", cex=0.7, x.intersp= 0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_B$labels), 
       fill=rev(PHYLA_label_cols_B$cols), 
       border=rev(PHYLA_label_cols_B$cols) )




## N DNRA
DNRA <- c("napA","napB","nrfA")
Kegg_N_DNRA_tax  <-  Kegg_N_tax %>% filter(Function %in% DNRA)
Kegg_N_DNRA_tab <- Kegg_N[rownames(Kegg_N_DNRA_tax),]

metadata1$Treatment1 <- factor(metadata1$Treatment1,c("CK","PE","DE"))

CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])

Kegg_N_DNRA_tab <-as.matrix(cbind(`CK`=apply(Kegg_N_DNRA_tab[,CK],1,sum),
                                  `PE`=apply(Kegg_N_DNRA_tab[,PE],1,sum),
                                  `DE`=apply(Kegg_N_DNRA_tab[,DE],1,sum)))


Kegg_N_DNRA_tab <- t(t(Kegg_N_DNRA_tab)/colSums(Kegg_N_DNRA_tab)) * 100

KO_N_DNRA <- names(sort(table(Kegg_N_DNRA_tax[,"Phylum"]), decr=T))

length(KO_N_DNRA)
sort(table(Kegg_N_DNRA_tax[,"Phylum"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_N_DNRA_tab)
for (i in KO_N_DNRA){
  x <- array(colSums(Kegg_N_DNRA_tab[rownames(Kegg_N_DNRA_tax)[which(Kegg_N_DNRA_tax$Phylum == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(KO_N_DNRA)
colnames(y) <- paste(colnames(Kegg_N_DNRA_tab))
KO_N_sum_DNRA <- y


# Phyla with MEAN abundances higher than 1% relative abundances
Kegg_N_DNRA_tax$labels <- Kegg_N_DNRA_tax$Phylum
Kegg_N_DNRA_tax$cols <- Kegg_N_DNRA_tax$Phylum

Kegg_N_DNRA_tax[ rownames(Kegg_N_DNRA_tax)[Kegg_N_DNRA_tax$labels=="Verrucomicrobia" ], ]$cols <- "gray"
Kegg_N_DNRA_tax[ rownames(Kegg_N_DNRA_tax)[Kegg_N_DNRA_tax$labels=="Proteobacteria" ], ]$cols <- "palegreen3"
Kegg_N_DNRA_tax[ rownames(Kegg_N_DNRA_tax)[Kegg_N_DNRA_tax$labels=="Acidobacteria" ], ]$cols <- "palegreen4"
Kegg_N_DNRA_tax[ rownames(Kegg_N_DNRA_tax)[Kegg_N_DNRA_tax$labels=="CandidatusOmnitrophica" ], ]$cols <- "gray"
Kegg_N_DNRA_tax[ rownames(Kegg_N_DNRA_tax)[Kegg_N_DNRA_tax$labels=="CandidatusRokubacteria" ], ]$cols <- "gray"
Kegg_N_DNRA_tax[ rownames(Kegg_N_DNRA_tax)[Kegg_N_DNRA_tax$labels=="Gemmatimonadetes" ], ]$cols <- "steelblue4"
Kegg_N_DNRA_tax[ rownames(Kegg_N_DNRA_tax)[Kegg_N_DNRA_tax$labels=="Chloroflexi" ], ]$cols <- "plum1"



## collaps OTU colors to prepare Phylum level colors
label_cols_B <- Kegg_N_DNRA_tax[, c("labels", "cols") ]
library(plyr)
PHYLA_label_cols_B <- ddply(label_cols_B, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_B) <- PHYLA_label_cols_B[,1]
PHYLA_label_cols_B

##### Plot Supplementary Figure S3
phylum_bar_B <- barplot(as.matrix(KO_N_sum_DNRA), col=PHYLA_label_cols_B[rownames(KO_N_sum_DNRA),]$cols, ylim=c(0,100), xaxt="n", border=NA, las=2)

legend(3, 100, bty="n", cex=0.7, x.intersp= 0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_B$labels), 
       fill=rev(PHYLA_label_cols_B$cols), 
       border=rev(PHYLA_label_cols_B$cols) )



## N ANRA
ANRA <- c("nasA")
Kegg_N_ANRA_tax  <-  Kegg_N_tax %>% filter(Function %in% ANRA)
Kegg_N_ANRA_tab <- Kegg_N[rownames(Kegg_N_ANRA_tax),]

metadata1$Treatment1 <- factor(metadata1$Treatment1,c("CK","PE","DE"))

CK <- rownames(metadata1[metadata1$Treatment1=="CK",])
PE <- rownames(metadata1[metadata1$Treatment1=="PE",])
DE <- rownames(metadata1[metadata1$Treatment1=="DE",])

Kegg_N_ANRA_tab <-as.matrix(cbind(`CK`=apply(Kegg_N_ANRA_tab[,CK],1,sum),
                                  `PE`=apply(Kegg_N_ANRA_tab[,PE],1,sum),
                                  `DE`=apply(Kegg_N_ANRA_tab[,DE],1,sum)))


Kegg_N_ANRA_tab <- t(t(Kegg_N_ANRA_tab)/colSums(Kegg_N_ANRA_tab)) * 100

KO_N_ANRA <- names(sort(table(Kegg_N_ANRA_tax[,"Phylum"]), decr=T))

length(KO_N_ANRA)
sort(table(Kegg_N_ANRA_tax[,"Phylum"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(Kegg_N_ANRA_tab)
for (i in KO_N_ANRA){
  x <- array(colSums(Kegg_N_ANRA_tab[rownames(Kegg_N_ANRA_tax)[which(Kegg_N_ANRA_tax$Phylum == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(KO_N_ANRA)
colnames(y) <- paste(colnames(Kegg_N_ANRA_tab))
KO_N_sum_ANRA <- y


# Phyla with MEAN abundances higher than 1% relative abundances
Kegg_N_ANRA_tax$labels <- Kegg_N_ANRA_tax$Phylum
Kegg_N_ANRA_tax$cols <- Kegg_N_ANRA_tax$Phylum


Kegg_N_ANRA_tax[ rownames(Kegg_N_ANRA_tax)[Kegg_N_ANRA_tax$labels=="Verrucomicrobia" ], ]$cols <- "gray"
Kegg_N_ANRA_tax[ rownames(Kegg_N_ANRA_tax)[Kegg_N_ANRA_tax$labels=="Proteobacteria" ], ]$cols <- "palegreen3"
Kegg_N_ANRA_tax[ rownames(Kegg_N_ANRA_tax)[Kegg_N_ANRA_tax$labels=="Acidobacteria" ], ]$cols <- "palegreen4"
Kegg_N_ANRA_tax[ rownames(Kegg_N_ANRA_tax)[Kegg_N_ANRA_tax$labels=="Gemmatimonadetes" ], ]$cols <- "steelblue4"
Kegg_N_ANRA_tax[ rownames(Kegg_N_ANRA_tax)[Kegg_N_ANRA_tax$labels=="Chloroflexi" ], ]$cols <- "plum1"
Kegg_N_ANRA_tax[ rownames(Kegg_N_ANRA_tax)[Kegg_N_ANRA_tax$labels=="CandidatusOmnitrophica" ], ]$cols <- "gray"
Kegg_N_ANRA_tax[ rownames(Kegg_N_ANRA_tax)[Kegg_N_ANRA_tax$labels=="CandidatusRokubacteria" ], ]$cols <- "gray"



## collaps OTU colors to prepare Phylum level colors
label_cols_B <- Kegg_N_ANRA_tax[, c("labels", "cols") ]
library(plyr)
PHYLA_label_cols_B <- ddply(label_cols_B, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_B) <- PHYLA_label_cols_B[,1]
PHYLA_label_cols_B

##### Plot Supplementary Figure S3
phylum_bar_B <- barplot(as.matrix(KO_N_sum_ANRA), col=PHYLA_label_cols_B[rownames(KO_N_sum_ANRA),]$cols, ylim=c(0,100), xaxt="n", border=NA, las=2)

legend(3, 100, bty="n", cex=0.7, x.intersp= 0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_B$labels), 
       fill=rev(PHYLA_label_cols_B$cols), 
       border=rev(PHYLA_label_cols_B$cols) )




#####  Calculate NTI #####
library(ape)
library(picante)
library(NST)
library(Rcpp)


## Bacteria
B_tree <- read.tree('B_tree.nwk')
B_tree_node <- B_tree[["tip.label"]]
cp_otu_B_tree <- cp_otu_B[c(B_tree_node),]
B_tree <- prune.sample(t(cp_otu_B_tree), B_tree)

## PD
B_comm <- as.matrix(t(cp_otu_B_tree))
B_pd <- pd(B_comm, B_tree,include.root=FALSE)
design_B$pd <- B_pd$PD

my_comparisons = list( c('CK', 'PE'), c("CK", "DE"), c("PE", "DE"))

CK  <- design_B[design_B$Treatment1=="CK", 13 ]
PE  <- design_B[design_B$Treatment1=="PE", 13 ]
DE  <- design_B[design_B$Treatment1=="DE", 13 ]

pd_B <- ggplot(design_B, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=pd, fill =Treatment1))+geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+labs(x="Type", y="Bacteria pd index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")



set.seed(619)
B_pnst <- pNST(comm = t(cp_otu_B_tree), tree = B_tree, group = design_B, phylo.shuffle = TRUE, taxo.null.model = NULL,pd.wd = tempdir(), abundance.weighted = TRUE, rand = 1000, nworker = 4, SES = TRUE, RC = FALSE)

B_betaMNTD <- B_pnst$index.pair
head(B_betaMNTD)
write.csv(B_betaMNTD,"B_betaMNTD.csv")

## Fungi
F_tree <- read.tree('F_otus.nwk')

F_tree_node <- F_tree[["tip.label"]]
a <- setdiff(F_tree_node,rownames(cp_otu_F))
F_tree<- drop.tip(F_tree,a)

F_tree1 <- prune.sample(t(cp_otu_F), F_tree)

## PD
F_comm <- as.matrix(t(cp_otu_F))
F_pd <- pd(F_comm, F_tree,include.root=FALSE)
design_F$pd <- F_pd$PD

my_comparisons = list( c('CK', 'PE'), c("CK", "DE"), c("PE", "DE"))

CK  <- design_F[design_F$Treatment1=="CK", 13 ]
PE  <- design_F[design_F$Treatment1=="PE", 13 ]
DE  <- design_F[design_F$Treatment1=="DE", 13 ]

pd_F <- ggplot(design_F, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=pd, fill =Treatment1))+geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+labs(x="Type", y="Fungi pd index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")

set.seed(619)
F_pnst <- pNST(comm = t(cp_otu_F), tree = F_tree1, group = design_F, phylo.shuffle = TRUE, taxo.null.model = NULL,pd.wd = tempdir(), abundance.weighted = TRUE, rand = 1000, nworker = 4, SES = TRUE, RC = FALSE)

F_betaMNTD <- F_pnst$index.pair
head(F_betaMNTD)
write.csv(F_betaMNTD,"F_betaMNTD.csv")


## Protist
P_tree <- read.tree('P_otus.nwk')
P_tree_node <- P_tree[["tip.label"]]
cp_otu_P_tree <- cp_otu_P[c(P_tree_node),]
P_tree <- prune.sample(t(cp_otu_P_tree), P_tree)

## PD
P_comm <- as.matrix(t(cp_otu_P_tree))
P_pd <- pd(P_comm, P_tree,include.root=FALSE)
design_P$pd <- P_pd$PD

my_comparisons = list( c('CK', 'PE'), c("CK", "DE"), c("PE", "DE"))

CK  <- design_P[design_P$Treatment1=="CK", 13 ]
PE  <- design_P[design_P$Treatment1=="PE", 13 ]
DE  <- design_P[design_P$Treatment1=="DE", 13 ]

pd_P <- ggplot(design_P, aes(x=factor(Treatment1,levels = c("CK","PE","DE")), y=pd, fill =Treatment1))+geom_violin(trim=FALSE,color="white") +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  scale_fill_manual(values=c("#DAA520","#3CB371","gray30"))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+labs(x="Type", y="Protist pd index")+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test")

grid.newpage()
grid.arrange(pd_B,pd_F,pd_P,nrow = 1)


##### lefse #####
library(microeco)
library(dplyr)
library(tidyverse)
library(magrittr)


## Bacteria all
keep_cp_otu_B <- which(rowSums(cp_otu_B) > 0)
cp_otu_B <- cp_otu_B[keep_cp_otu_B,]

nrow(cp_otu_B)

tax_B <- tax_B[rownames(cp_otu_B),]
tax_B$otu <- rownames(tax_B)

B_dataset <- microtable$new(sample_table = design_B,
                                 otu_table = cp_otu_B, 
                                 tax_table = tax_B)


B_lefse <- trans_diff$new(dataset = B_dataset, 
                               method = "lefse", 
                               group = "Treatment1", 
                               alpha = 0.05, taxa_level = "otu",
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

B_lefse_plot <- B_lefse$plot_diff_bar(use_number = 1:30, 
                                                width = 0.8, heatmap_y = "Taxa",
                                                group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() + 
  ggsci::scale_fill_npg()


## Fungi all
keep_cp_otu_F <- which(rowSums(cp_otu_F) > 0)
cp_otu_F <- cp_otu_F[keep_cp_otu_F,]

nrow(cp_otu_F)

tax_F <- tax_F[rownames(cp_otu_F),]
tax_F$otu <- rownames(tax_F)

F_dataset <- microtable$new(sample_table = design_F,
                            otu_table = cp_otu_F, 
                            tax_table = tax_F)


F_lefse <- trans_diff$new(dataset = F_dataset, 
                          method = "lefse", 
                          group = "Treatment1", 
                          alpha = 0.05, taxa_level = "otu",
                          lefse_subgroup = NULL,
                          p_adjust_method = "none")

F_lefse_plot <- F_lefse$plot_diff_bar(use_number = 1:30, 
                                      width = 0.8, heatmap_y = "Taxa",
                                      group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() + 
  ggsci::scale_fill_npg()


## Protist all
keep_cp_otu_P <- which(rowSums(cp_otu_P) > 0)
cp_otu_P <- cp_otu_P[keep_cp_otu_P,]

nrow(cp_otu_P)

tax_P <- tax_P[rownames(cp_otu_P),]
tax_P$otu <- rownames(tax_P)

design_P <- droplevels(design_P[soilsamples,])

P_dataset <- microtable$new(sample_table = design_P,
                            otu_table = cp_otu_P, 
                            tax_table = tax_P)


P_lefse <- trans_diff$new(dataset = P_dataset, 
                          method = "lefse", 
                          group = "Treatment1", 
                          alpha = 0.05, taxa_level = "otu",
                          lefse_subgroup = NULL,
                          p_adjust_method = "none")

P_lefse_plot <- P_lefse$plot_diff_bar(use_number = 1:30, 
                                      width = 0.8, heatmap_y = "Taxa",
                                      group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() + 
  ggsci::scale_fill_npg()


grid.newpage()
grid.arrange(B_lefse_plot,F_lefse_plot,P_lefse_plot,nrow = 1)



## cladogram
# bacteria
B_cladogram_ASV <- c(B_lefse[["plot_diff_bar_taxa"]])
cp_otu_B_cladogram <- cp_otu_B[B_cladogram_ASV,]


PE <- rownames(design_B[design_B$Treatment1=="PE",])
DE <- rownames(design_B[design_B$Treatment1=="DE",])
CK <- rownames(design_B[design_B$Treatment1=="CK",])


B_cladogram_ASV_tab <-as.matrix(cbind(`CK`=apply(cp_otu_B_cladogram[,CK],1,sum),
                                      `PE`=apply(cp_otu_B_cladogram[,PE],1,sum),
                                      `DE`=apply(cp_otu_B_cladogram[,DE],1,sum)))

B_cladogram_ASV_tax <- tax_B[B_cladogram_ASV,]

## fungi
F_cladogram_ASV <- c(F_lefse[["plot_diff_bar_taxa"]])
cp_otu_F_cladogram <- cp_otu_F[F_cladogram_ASV,]


PE <- rownames(design_F[design_F$Treatment1=="PE",])
DE <- rownames(design_F[design_F$Treatment1=="DE",])
CK <- rownames(design_F[design_F$Treatment1=="CK",])


F_cladogram_ASV_tab <-as.matrix(cbind(`CK`=apply(cp_otu_F_cladogram[,CK],1,sum),
                                      `PE`=apply(cp_otu_F_cladogram[,PE],1,sum),
                                      `DE`=apply(cp_otu_F_cladogram[,DE],1,sum)))

F_cladogram_ASV_tax <- tax_F[F_cladogram_ASV,]


## protist
P_cladogram_ASV <- c(P_lefse[["plot_diff_bar_taxa"]])
cp_otu_P_cladogram <- cp_otu_P[P_cladogram_ASV,]


PE <- rownames(design_P[design_P$Treatment1=="PE",])
DE <- rownames(design_P[design_P$Treatment1=="DE",])
CK <- rownames(design_P[design_P$Treatment1=="CK",])


P_cladogram_ASV_tab <-as.matrix(cbind(`CK`=apply(cp_otu_P_cladogram[,CK],1,sum),
                                      `PE`=apply(cp_otu_P_cladogram[,PE],1,sum),
                                      `DE`=apply(cp_otu_P_cladogram[,DE],1,sum)))

P_cladogram_ASV_tax <- tax_P[P_cladogram_ASV,]



## Bacteria soil
cp_otu_B_soil <- cp_otu_B[,soilsamples]

keep_cp_otu_B_soil <- which(rowSums(cp_otu_B_soil) > 0)
cp_otu_B_soil <- cp_otu_B_soil[keep_cp_otu_B_soil,]

nrow(cp_otu_B_soil)

tax_B_soil <- tax_B[rownames(cp_otu_B_soil),]
tax_B_soil$otu <- rownames(tax_B_soil)

design_B_soil <- droplevels(design_B[soilsamples,])

B_soil_dataset <- microtable$new(sample_table = design_B_soil,
                          otu_table = cp_otu_B_soil, 
                          tax_table = tax_B_soil)


B_soil_lefse <- trans_diff$new(dataset = B_soil_dataset, 
                        method = "lefse", 
                        group = "Treatment1", 
                        alpha = 0.05, taxa_level = "otu",
                        lefse_subgroup = NULL,
                        p_adjust_method = "none")

B_soil_lefse_plot <- B_soil_lefse$plot_diff_bar(use_number = 1:30, 
                    width = 0.8, heatmap_y = "Taxa",
                    group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() + 
  ggsci::scale_fill_npg()



## Bacteria root
cp_otu_B_root <- cp_otu_B[,rootsamples]

keep_cp_otu_B_root <- which(rowSums(cp_otu_B_root) > 0)
cp_otu_B_root <- cp_otu_B_root[keep_cp_otu_B_root,]

nrow(cp_otu_B_root)

tax_B_root <- tax_B[rownames(cp_otu_B_root),]
tax_B_root$otu <- rownames(tax_B_root)

design_B_root <- droplevels(design_B[rootsamples,])

B_root_dataset <- microtable$new(sample_table = design_B_root,
                                 otu_table = cp_otu_B_root, 
                                 tax_table = tax_B_root)


B_root_lefse <- trans_diff$new(dataset = B_root_dataset, 
                               method = "lefse", 
                               group = "Treatment1", taxa_level = "otu",
                               alpha = 0.01, 
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

B_root_lefse_plot <- B_root_lefse$plot_diff_bar(use_number = 1:30, 
                           width = 0.8, 
                           group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Bacteria end
cp_otu_B_end <- cp_otu_B[,endsamples]

keep_cp_otu_B_end <- which(rowSums(cp_otu_B_end) > 0)
cp_otu_B_end <- cp_otu_B_end[keep_cp_otu_B_end,]

nrow(cp_otu_B_end)

tax_B_end <- tax_B[rownames(cp_otu_B_end),]
tax_B_end$otu <- rownames(tax_B_end)

design_B_end <- droplevels(design_B[endsamples,])

B_end_dataset <- microtable$new(sample_table = design_B_end,
                                 otu_table = cp_otu_B_end, 
                                 tax_table = tax_B_end)


B_end_lefse <- trans_diff$new(dataset = B_end_dataset, 
                               method = "lefse", 
                               group = "Treatment1", taxa_level = "otu",
                               alpha = 0.05, 
                               lefse_subgroup = NULL,
                              p_adjust_method = "none")

B_end_lefse_plot <- B_end_lefse$plot_diff_bar(use_number = 1:30, 
                           width = 0.8, 
                           group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Bacteria film
cp_otu_B_film <- cp_otu_B[,particlesamples]

keep_cp_otu_B_film <- which(rowSums(cp_otu_B_film) > 0)
cp_otu_B_film <- cp_otu_B_film[keep_cp_otu_B_film,]

nrow(cp_otu_B_film)

tax_B_film <- tax_B[rownames(cp_otu_B_film),]
tax_B_film$otu <- rownames(tax_B_film)

design_B_film <- droplevels(design_B[particlesamples,])

B_film_dataset <- microtable$new(sample_table = design_B_film,
                                 otu_table = cp_otu_B_film, 
                                 tax_table = tax_B_film)


B_film_lefse <- trans_diff$new(dataset = B_film_dataset, 
                               method = "lefse", 
                               group = "Treatment1", 
                               alpha = 0.05, taxa_level = "otu",
                               lefse_subgroup = NULL)

B_film_lefse_plot <- B_film_lefse$plot_diff_bar(use_number = 1:30, 
                           width = 0.8, 
                           group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()


## Fungi soil
cp_otu_F_soil <- cp_otu_F[,soilsamples]

keep_cp_otu_F_soil <- which(rowSums(cp_otu_F_soil) > 0)
cp_otu_F_soil <- cp_otu_F_soil[keep_cp_otu_F_soil,]

nrow(cp_otu_F_soil)

tax_F_soil <- tax_F[rownames(cp_otu_F_soil),]
tax_F_soil$otu <- rownames(tax_F_soil)

design_F_soil <- droplevels(design_F[soilsamples,])

F_soil_dataset <- microtable$new(sample_table = design_F_soil,
                                 otu_table = cp_otu_F_soil, 
                                 tax_table = tax_F_soil)


F_soil_lefse <- trans_diff$new(dataset = F_soil_dataset, 
                               method = "lefse", 
                               group = "Treatment1", taxa_level = "otu",
                               alpha = 0.05, 
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

F_soil_lefse_plot <- F_soil_lefse$plot_diff_bar(use_number = 1:30, 
                           width = 0.8, 
                           group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Bacteria root
cp_otu_F_root <- cp_otu_F[,rootsamples]

keep_cp_otu_F_root <- which(rowSums(cp_otu_F_root) > 0)
cp_otu_F_root <- cp_otu_F_root[keep_cp_otu_F_root,]

nrow(cp_otu_F_root)

tax_F_root <- tax_F[rownames(cp_otu_F_root),]
tax_F_root$otu <- rownames(tax_F_root)

design_F_root <- droplevels(design_F[rootsamples,])

F_root_dataset <- microtable$new(sample_table = design_F_root,
                                 otu_table = cp_otu_F_root, 
                                 tax_table = tax_F_root)


F_root_lefse <- trans_diff$new(dataset = F_root_dataset, 
                               method = "lefse", 
                               group = "Treatment1", 
                               alpha = 0.01, taxa_level = "otu",
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

F_root_lefse_plot <- F_root_lefse$plot_diff_bar(use_number = 1:30, 
                           width = 0.8, 
                           group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Bacteria end
cp_otu_F_end <- cp_otu_F[,endsamples]

keep_cp_otu_F_end <- which(rowSums(cp_otu_F_end) > 0)
cp_otu_F_end <- cp_otu_F_end[keep_cp_otu_F_end,]

nrow(cp_otu_F_end)

tax_F_end <- tax_F[rownames(cp_otu_F_end),]
tax_F_end$otu <- rownames(tax_F_end)

design_F_end <- droplevels(design_F[endsamples,])

F_end_dataset <- microtable$new(sample_table = design_F_end,
                                otu_table = cp_otu_F_end, 
                                tax_table = tax_F_end)


F_end_lefse <- trans_diff$new(dataset = F_end_dataset, 
                              method = "lefse", 
                              group = "Treatment1", taxa_level = "otu",
                              alpha = 0.05, 
                              lefse_subgroup = NULL,
                              p_adjust_method = "none")

F_end_lefse_plot <- F_end_lefse$plot_diff_bar(use_number = 1:30, 
                          width = 0.8, 
                          group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Bacteria film
cp_otu_F_film <- cp_otu_F[,particlesamples]

keep_cp_otu_F_film <- which(rowSums(cp_otu_F_film) > 0)
cp_otu_F_film <- cp_otu_F_film[keep_cp_otu_F_film,]

nrow(cp_otu_F_film)

tax_F_film <- tax_F[rownames(cp_otu_F_film),]
tax_F_film$otu <- rownames(tax_F_film)

design_F_film <- droplevels(design_F[particlesamples,])

F_film_dataset <- microtable$new(sample_table = design_F_film,
                                 otu_table = cp_otu_F_film, 
                                 tax_table = tax_F_film)


F_film_lefse <- trans_diff$new(dataset = F_film_dataset, 
                               method = "lefse", 
                               group = "Treatment1", taxa_level = "otu",
                               alpha = 0.05, 
                               lefse_subgroup = NULL)

F_film_lefse_plot <- F_film_lefse$plot_diff_bar(use_number = 1:30, 
                           width = 0.8, 
                           group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Protist soil
cp_otu_P_soil <- cp_otu_P[,psoilsamples]

keep_cp_otu_P_soil <- which(rowSums(cp_otu_P_soil) > 0)
cp_otu_P_soil <- cp_otu_P_soil[keep_cp_otu_P_soil,]

nrow(cp_otu_P_soil)

tax_P_soil <- tax_P[rownames(cp_otu_P_soil),]
tax_P_soil$otu <- rownames(tax_P_soil)

design_P_soil <- droplevels(design_P[psoilsamples,])

P_soil_dataset <- microtable$new(sample_table = design_P_soil,
                                 otu_table = cp_otu_P_soil, 
                                 tax_table = tax_P_soil)


P_soil_lefse <- trans_diff$new(dataset = P_soil_dataset, 
                               method = "lefse", 
                               group = "Treatment1", 
                               alpha = 0.05, taxa_level = "otu",
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

P_soil_lefse_plot <- P_soil_lefse$plot_diff_bar(use_number = 1:30, 
                           width = 0.8, 
                           group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Protist root
cp_otu_P_root <- cp_otu_P[,rootsamples]

keep_cp_otu_P_root <- which(rowSums(cp_otu_P_root) > 0)
cp_otu_P_root <- cp_otu_P_root[keep_cp_otu_P_root,]

nrow(cp_otu_P_root)

tax_P_root <- tax_P[rownames(cp_otu_P_root),]
tax_P_root$otu <- rownames(tax_P_root)

design_P_root <- droplevels(design_P[rootsamples,])

P_root_dataset <- microtable$new(sample_table = design_P_root,
                                 otu_table = cp_otu_P_root, 
                                 tax_table = tax_P_root)


P_root_lefse <- trans_diff$new(dataset = P_root_dataset, 
                               method = "lefse", 
                               group = "Treatment1", taxa_level = "otu",
                               alpha = 0.01, 
                               lefse_subgroup = NULL,
                               p_adjust_method = "none")

P_root_lefse_plot <- P_root_lefse$plot_diff_bar(use_number = 1:30, 
                           width = 0.8, 
                           group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Protist end
cp_otu_P_end <- cp_otu_P[,pendsamples]

keep_cp_otu_P_end <- which(rowSums(cp_otu_P_end) > 0)
cp_otu_P_end <- cp_otu_P_end[keep_cp_otu_P_end,]

nrow(cp_otu_P_end)

tax_P_end <- tax_P[rownames(cp_otu_P_end),]
tax_P_end$otu <- rownames(tax_P_end)

design_P_end <- droplevels(design_P[pendsamples,])

P_end_dataset <- microtable$new(sample_table = design_P_end,
                                otu_table = cp_otu_P_end, 
                                tax_table = tax_P_end)


P_end_lefse <- trans_diff$new(dataset = P_end_dataset, 
                              method = "lefse", 
                              group = "Treatment1", taxa_level = "otu",
                              alpha = 0.05, 
                              lefse_subgroup = NULL,
                              p_adjust_method = "none")

P_end_lefse_plot <- P_end_lefse$plot_diff_bar(use_number = 1:30, 
                          width = 0.8, 
                          group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()



## Protist film
cp_otu_P_film <- cp_otu_P[,particlesamples]

keep_cp_otu_P_film <- which(rowSums(cp_otu_P_film) > 0)
cp_otu_P_film <- cp_otu_P_film[keep_cp_otu_P_film,]

nrow(cp_otu_P_film)

tax_P_film <- tax_P[rownames(cp_otu_P_film),]
tax_P_film$otu <- rownames(tax_P_film)

design_P_film <- droplevels(design_P[particlesamples,])

P_film_dataset <- microtable$new(sample_table = design_P_film,
                                 otu_table = cp_otu_P_film, 
                                 tax_table = tax_P_film)


P_film_lefse <- trans_diff$new(dataset = P_film_dataset, 
                               method = "lefse", taxa_level = "otu",
                               group = "Treatment1", 
                               alpha = 0.05, 
                               lefse_subgroup = NULL)

P_film_lefse_plot <- P_film_lefse$plot_diff_bar(use_number = 1:30, 
                           width = 0.8, 
                           group_order = c("CK", "PE", "DE")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()


grid.newpage()
grid.arrange(B_soil_lefse_plot ,B_root_lefse_plot,B_end_lefse_plot,B_film_lefse_plot,ncol = 4)

grid.newpage()
grid.arrange(F_soil_lefse_plot ,F_root_lefse_plot,F_end_lefse_plot,F_film_lefse_plot,ncol = 4)

grid.newpage()
grid.arrange(P_soil_lefse_plot ,P_root_lefse_plot,P_end_lefse_plot,P_film_lefse_plot,ncol = 4)
                 
