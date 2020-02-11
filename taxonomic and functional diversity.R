

Appendix A. R codes for calculating taxonomic diversity index
# Set Working Directory
setwd("N:/github/functional_diversity_training")
# Load the package
library(vegan)   
# Read the data file
fish<-read.csv("data/site_matrix_invented.csv",header=T, row.names=1,sep=",")
tax<-read.csv("data/tax_matrix_invented.csv",header=T, row.names=1, sep=",")  
#Transform the dataframe for diversity calculation
#(369 fish * 56 basin transform to 56 basin*369 fish)
t_fish<-t(fish)
# Taxonomic distances from a classification table with variable step lengths
taxdis <- taxa2dist(tax, varstep=TRUE)                      
# Calculate taxonomic diversity and distinctness
mod <- taxondive(t_fish, taxdis)
# List the result
mod

Appendix B. R codes for calculating functional diversity index

# Load the package
library(vegan)  
library(FD) 
# Read the data file
fish<-read.csv("data/site_matrix_invented.csv",header=T, row.names=1,sep=",")
tax<-read.csv("data/tax_matrix_invented.csv",header=T, row.names=1, sep=",")  
trait<-read.csv("data/traits_matrix_invented.csv",header=T, row.names=1, sep=",")
#Prepare trait data for FRic calculation
gower.f<-gowdis(trait)   # Gower’s distance analysis
pcoa.g<-pcoa(gower.f,rn=NULL )   # PCoA analysis
Fun_pcoa<-pcoa.g$vectors[,1:3]    # the first three axes of PCoA were kept as 
# the synthetic functional traits
#Transform the dataframe for diversity calculation
t_fish<-t(fish)                   #(369 fish * 56 basin transform to 56 basin*369 fish)
# Calculate functional diversity indices
FD.alpha<-dbFD(Fun_pcoa, t_fish)
# List the result
FD.alpha$FRic

Appendix C. R codes for calculating β-diversity index, and its component

# Load the packages
library(vegan)  
library(FD) 
library(betapart)
# Read the data file
fish<-read.csv("data/site_matrix_invented.csv",header=T, row.names=1,sep=",")
tax<-read.csv("data/tax_matrix_invented.csv",header=T, row.names=1, sep=",")  
trait<-read.csv("data/traits_matrix_invented.csv",header=T, row.names=1, sep=",")
#Prepare trait data for calculation
gower.f<-gowdis(trait)   # Gower’s distance analysis
pcoa.f<-pcoa(gower.f,rn=NULL )   # PCoA analysis
Fun_pcoa<-pcoa.f$vectors[,1:3]    # keep the first three axes of PCoA as
# the synthetic functional traits
#Prepare taxonomic data for calculation
gower.t<-gowdis(tax,ord = c("classic"))    # Gower’s distance analysis
pcoa.t<-pcoa(gower.t,rn=NULL )         # PCoA analysis
tax_pcoa<-pcoa.t$vectors[,1:2]          # the first 2 axes of PCoA were kept as
# the synthetic taxonomic traits
#Transform the dataframe for diversity calculation
t_fish<-t(fish)                   #(N fish * M basin transform to M basin*N fish)
#Analysis the beta indice between sub-basins
SR_beta<- beta.pair(t_fish, index.family="jaccard")
FD_beta<-functional.beta.pair(t_fish, Fun_pcoa, index.family="jaccard")
TD_beta<-functional.beta.pair(t_fish, tax_pcoa, index.family="jaccard")
#Get the results
SR.jtu<-as.matrix(SR_beta$beta.jtu)
SR.jne<-as.matrix(SR_beta$beta.jne)
SR.jac<-as.matrix(SR_beta$beta.jac)
FD.jtu<-as.matrix(FD_beta$funct.beta.jtu)
FD.jne<-as.matrix(FD_beta$funct.beta.jne)
FD.jac<-as.matrix(FD_beta$funct.beta.jac)
TD.jtu<-as.matrix(TD_beta$funct.beta.jtu)
TD.jne<-as.matrix(TD_beta$funct.beta.jne)
TD.jac<-as.matrix(TD_beta$funct.beta.jac)
#average beta between sub-basins, for M sub-basins, the summary should be divided by (M-1).
Beta_bs<-cbind(colSums(FD.jac)/ (14-1), colSums(FD.jtu)/ (14-1), colSums(FD.jne)/ (14-1), 
               colSums(TD.jac)/ (14-1), colSums(TD.jtu)/ (14-1), colSums(TD.jne)/ (14-1),
               colSums(SR.jac)/ (14-1), colSums(SR.jtu)/ (14-1), colSums(SR.jne)/ (14-1))
colnames(Beta_bs)=c('avg_FD.jac', 'avg_FD.jtu', 'avg_FD.jne',  'avg_TD.jac', 'avg_TD.jtu', 'avg_TD.jne', 
                    'avg_SR.jac',  'avg_SR.jtu',  'avg_SR.jne')
Beta_bs
# Load the package
#install.packages ("pheatmap")
library(pheatmap)
# Heatmap of beta-diversity between sub-basins
pheatmap(FD.jtu,main="FD.jtu")
pheatmap(FD.jne,main="FD.jne")
pheatmap(FD.jac,main="FD.jac")
pheatmap(TD.jtu,main="TD.jtu")
pheatmap(TD.jne,main="TD.jne")
pheatmap(TD.jac,main="TD.jac")
pheatmap(SR.jtu,main="SR.jtu")
pheatmap(SR.jne,main="SR.jne")
pheatmap(SR.jac,main="SR.jac")
#Calculate the proportion of turnover and nestedness components in the overall dissimilarity
#of functional beta diversity and taxonomic beta diversity
B1<- FD.jtu*FD.jac^-1
B2<- FD.jne*FD.jac^-1
B2[lower.tri(B2)]=0
B1[upper.tri(B1)]=0
FD_beta.r<- B1+B2
pheatmap(FD_beta.r, cluster_rows = F,  cluster_cols = F, show_rownames = T, show_colnames = T)
B1<- TD.jtu*TD.jac^-1
B2<- TD.jne*TD.jac^-1
B2[lower.tri(B2)]=0
B1[upper.tri(B1)]=0
TD_beta.r<- B1+B2
pheatmap(TD_beta.r, cluster_rows = F,  cluster_cols = F, show_rownames= T, show_colnames = T)
B1<- SR.jtu*SR.jac^-1
B2<- SR.jne*SR.jac^-1
B2[lower.tri(B2)]=0
B1[upper.tri(B1)]=0
SR_beta.r<- B1+B2
pheatmap(SR_beta.r, cluster_rows = F,  cluster_cols = F, show_rownames= T, show_colnames = T)



##### 4 communities in a 2D functional space (convex hulls are rectangles)
traits.test<-cbind( c(1,1,1,2,2,3,3,4,4,5,5) , c(1,2,4,1,2,3,5,1,4,3,5) )
dimnames(traits.test)<-list(paste("sp",1:11,sep="") , c("Trait 1","Trait 2") ) 

comm.test<-matrix(0,4,11,dimnames=list( c("A","B","C","D") , paste("sp",1:11,sep="") ) )
comm.test["A",c(1,2,4,5)]<-1
comm.test["B",c(1,3,8,9)]<-1
comm.test["C",c(6,7,10,11)]<-1
comm.test["D",c(2,4,7,9)]<-1

plot(5,5,xlim=c(0,6), ylim=c(0,6), type="n", xlab="Trait 1",ylab="Trait 2")
points(traits.test[,1],traits.test[,2], pch=21,cex=1.5,bg="black")
rect(1,1,4,4, col="#458B0050", border="#458B00") ; text(2.5,2.5,"B",col="#458B00",cex=1.5)	
polygon(c(2,1,3,4), c(1,2,5,4), col="#DA70D650", border="#DA70D6") ; 
text(2.5,3,"D",col="#DA70D6",cex=1.5)	
rect(1,1,2,2, col="#FF000050" , border="#FF0000") ; text(1.5,1.5,"A",col="#FF0000",cex=1.5)	
rect(3,3,5,5, col="#1E90FF50", border="#1E90FF") ; text(4,4.2,"C",col="#1E90FF",cex=1.5)	

test.pair<-functional.beta.pair(x=comm.test, traits=traits.test, index.family = "jaccard" )
lapply(test.pair,round,2)

#PCA

traits.pca <- prcomp(trait[,c(1:7)], center = TRUE,scale. = TRUE)
summary (traits.pca)
str (traits.pca)

mtcars.pca <- prcomp(mtcars[,c(1:7,10,11)], center = TRUE,scale. = TRUE)
summary(mtcars.pca)

library(devtools)
install_github("vqv/ggbiplot")

library(ggbiplot)

ggbiplot(mtcars.pca)
ggbiplot (traits.pca)


ggbiplot(mtcars.pca, labels=rownames(mtcars))

ggbiplot(traits.pca, labels=rownames(trait))

mtcars.country <- c(rep("Japan", 3), rep("US",4), rep("Europe", 7),rep("US",3), "Europe", rep("Japan", 3), rep("US",4), rep("Europe", 3), "US", rep("Europe", 3))

ggbiplot(mtcars.pca,ellipse=TRUE,  labels=rownames(mtcars), groups=mtcars.country)

traits.demerspelag <- c("pelagic","reef-associated", "pelagic", "pelagic", "pelagic", "pelagic", "demersal", "pelagic", "demersal", "reef-associated")

ggbiplot(traits.pca,ellipse=TRUE,  labels=rownames(trait), groups=traits.demerspelag)



traits.bodyshape <- c("elongated", "fusiform_normal", "elongated", "fusiform_normal", "fusiform_normal", "short_andor_deep", "fusiform_normal", "fusiform_normal", "other", "short_andor_deep")

ggbiplot(traits.pca,ellipse=TRUE,  labels=rownames(trait), groups=traits.bodyshape)

ggbiplot(traits.pca,ellipse=TRUE,circle=TRUE, labels=rownames(trait), groups=traits.bodyshape)

ggbiplot(traits.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1, labels=rownames(trait), groups=traits.bodyshape)


#To install the package run:
  
install.packages("devtools")
require(devtools)
install_github("ibartomeus/fundiv")
require(fundiv)

#To calculate dendogram indexes:

FD_dendro <- FD_dendro(S = trait, A = t_fish, Cluster.method = "average", ord = "podani",
                       Weigthedby = "abundance")
FD_dendro

FD_all <- FDindexes(S = trait, A = t_fish, Distance.method= "gower", ord= "podani", 
                    Cluster.method= "average", corr= "cailliez", Weigthedby = "abundance")
FD_all

results_ibartomeus <- write.csv(FD_all, "results_ibartomeus.csv")

#We can see that both families of indexes are correlated

plot(FD_all$FDpg ~ FD_all$Frich)

#Note that Functional richness indexes are highly correlated with richness, hence you may also want to know if those indexes are higher or lower than expected by its richness levels (See Rader et al. 2014 D&D for details http://onlinelibrary.wiley.com/doi/10.1111/ddi.12221/abstract)

null <- null.FD(S = trait, A = t_fish, it = 100, w = NA)
null

#No much significant results in this dataset, but is not surprising given:
  
plot(FD_all$FDpg ~ FD_all$n_sp)

#You can calculate standardized FD indexes as in Rader et al. 2014 by

(null$FD - null$null_meanFD) / null$null_sdFD

#Indexes implemented in Clark et al. 2012 (Plos One) are also available (FD_Clark), but in my opinion they perform very similar that the proposed FDw, but are computational time consuming because they need to create a dendogram for each community, instead of using only a general dendogram containing all species. This last approach is better according to the literature...

#FUNCTIONAL GROUPS APPROACH

traits_grouped <- trait %>% 
  mutate (Loo_class= cut(Loo, c(0,25,50,100,1500), labels=c("short", "medium", "long", "very_long"))) %>% 
  mutate (K_class= cut(K, c(0,0.5,1,1.5,2.5), labels=c("low_growth", "medium_growth", "high_growth", "very_high_growth"))) %>%
  mutate (W_class= cut(Winfinity, c(0,100,1000,10000,6000000), labels=c("light_weigth", "medium_weigth", "heavy_weigth", "very_heavy_weigth"))) %>%
  mutate (age_class= cut(tmax, c(0,5,15,25,55), labels=c("short_life", "medium_life", "long_life", "verylong_life"))) %>% 
  mutate (mature_class= cut(tm, c(0,1,5,10,25), labels=c("fast_mature", "medium_mature", "long_mature", "verylong_mature"))) %>%
  mutate (mortality_class= cut(M, c(0,0.5,1,2,3), labels=c("low_mortal", "medium_mortal", "high_mortal", "veryhigh_mortal"))) %>%
  mutate (maturelength_class=cut(Lm, c(0,10,20,100,500), labels=c("short_maturelength", "medium_maturelength", "long_maturelength", "verylong_maturelength"))) %>%
  mutate (deeprange_class=cut(depthrangedeep, c(0,250,500,2500,5000), labels=c("shallow", "medium_deep", "deep", "very_deep")))

colnames (traits_grouped)

functional_groups <- traits_grouped %>% 
  group_by(Loo_class,K_class, W_class, age_class, mature_class, mortality_class, maturelength_class, deeprange_class, demerspelag, FeedingType, bodyshape) %>% 
  summarise(n=n())








