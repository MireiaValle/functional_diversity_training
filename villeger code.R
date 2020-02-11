# FUNCTIONAL SPACE

##############################################################################################
# IMPORTING DATASETS

# importing biomass of species in assemblages from a txt file: species names are in the first column and will be used as row names
weight_fruits_baskets<-read.table("weight_fruits_baskets.txt", header=T, row.names =1)
weight_fruits_baskets<-as.matrix(weight_fruits_baskets)
weight_fruits_baskets

#importing presence-absence fish
fish<-read.csv("data/site_matrix_invented.csv",header=T, row.names=1,sep=",")
t_fish <- t(fish)

##tax<-read.csv("data/tax_matrix_invented.csv",header=T, row.names=1, sep=",")  


############################################

# importing species raw trait values from a txt file: row names are in the column named 'Species'
traits_raw_fruits<-read.table("traits_raw_fruits.txt", header=T, row.names = "Species")
traits_raw_fruits

#traits fish
trait<-read.csv("data/traits_matrix_invented.csv",header=T, row.names=1, sep=",")

# checking that species names are the same in the two matrices
sum( row.names(trait) %in% colnames(t_fish) ) == ncol(t_fish)

# checking that species names are the same in the two matrices
sum( row.names(traits_raw_fruits) %in% colnames(weight_fruits_baskets) ) == ncol(weight_fruits_baskets)


# computing all functional spaces based on dendrogram or PCoA (up to 10 axes)
# NOTE: you need to be connected to internet for a correct functioning of quality_funct_space function

my_path<-"N:/github/Tutorial_FD"
# loading functions
setwd( paste(my_path,"/functions", sep="") ) # folder where R functions have been saved
source("quality_funct_space.R")
source("plot_funct_space.R")
source("multidimFD.R")
source("multidimFbetaD.R")

install.packages ("ape")
install.packages("clue")
install.packages("cluster")
install.packages("geometry")
install.packages("gtools")

setwd ("N:/github/functional_diversity_training/FDresults")   # setting working directory for results

qual_funct_space_clear<-quality_funct_space(trait, traits_weights=NULL, nbdim=10, metric="Gower", dendro=TRUE, plot="quality_funct_space_MarineFish")

qual_funct_space_clear$meanSD # => best space has 5 dimensions + have a look to the "quality_funct_space_MarineFish.jpeg" file in the /results folder
# mean Squared-Deviation of 0.074876 means that average deviation between Euclidean distance and Gower's distance is of (0.074)^0.5=0.27, so it can be seen like an average error of 3%

# species coordinates in the best space
coord_fish_3D<-qual_funct_space_clear$details_funct_space$mat_coord[,1:3]

################################# HASTA AQUÍ CON LA EXTENT NUEVA
##################################
# plot of best functional dendrogram
plot(qual_funct_space_clear$details_funct_space$best_tree, sub="", main="UPGMA", xlab="", h=-1)

##############################################################################################
##############################################################################################

# MULTIDIMENSIONAL FUNCTIONAL DIVERISTY INDICES

# computing Functional diversity indices with plots of FD indices put in a subfolder named plot_FD

my_path<-"N:/github/functional_diversity_training"

FD_assemblages<-multidimFD(coord_fish_3D, t_fish , check_species_pool=TRUE, verb=TRUE,
                           nm_asb_plot=row.names(t_fish), folder_plot= paste(my_path,"FDresults/plot_FD", sep=""),
                           Faxes_plot=colnames(coord_fish_3D)[1:3], Faxes_nm_plot=colnames(coord_fish_3D)[1:3],
                           plot_pool=TRUE, col_bg="grey90", col_sp_pool="grey30",  pch_sp_pool="+", cex_sp_pool=1, 
                           pch_sp=21, col_sp="#1E90FF", transp=50 )

# printing results = rounded FD indices values
results <- round(FD_assemblages,3)
write.csv(results, "results_Villeger.csv")

#plotting results

rich_fd <- subset (FD_assemblages, select=c("Nb_sp", "FRic", "FDiv", "FEve", "FDis", "FSpe", "FOri"))
rich_fd <- as.data.frame(rich_fd )

#FRic=Functional richness = convex hull volume

ggplot(rich_fd, aes(x =Nb_sp, y = FRic )) +
  geom_point()

#FDiv=Functional Divergence

ggplot(rich_fd, aes(x =Nb_sp, y = FDiv)) +
  geom_point()

#FEve=Functional Evenness

ggplot(rich_fd, aes(x =Nb_sp, y = FEve)) +
  geom_point()

#FDis=Functional Dispersion: abundance-weighted mean distance to abundance-weighted centroid scaled by maximum value possible given species pool (i.e. the two most distant species have half of total weight)

ggplot(rich_fd, aes(x =Nb_sp, y = FDis, na.rm = TRUE)) +
  geom_point()

#FSpe=Functional Specialization: abundance-weighted mean distance to centroid of the global pool of species scaled by maximum value possible given species pool (i.e. an assmeblage hosting only the most specialized species)

ggplot(rich_fd, aes(x =Nb_sp, y = FSpe)) +
  geom_point()

#FOri=Functional Originality : abundance-weighted mean distance to nearest neighbour in the global pool of species scaled by maximum value possible given species pool (i.e. an assmeblage hosting only the most original species)
ggplot(rich_fd, aes(x =Nb_sp, y = FOri)) +
  geom_point()

#categories for congruence calculations see Martín-Regalado et al 2019

summary (rich_fd$Nb_sp)

#selecting the category ranges manually: 
#rich_fd$categoryNb_sp[rich_fd$Nb_sp < 500] <- "low"
#rich_fd$categoryNb_sp[rich_fd$Nb_sp > 501 & rich_fd$Nb_sp < 1000] <- "middle"
#rich_fd$categoryNb_sp[rich_fd$Nb_sp > 1001] <- "high"

#selecting the category ranges by quantile:
xs<-quantile(rich_fd$Nb_sp,c(0,1/3,2/3,1))
xs[1]<- xs[1]-.00005
df1 <- rich_fd %>% mutate(category=cut(Nb_sp, breaks=xs, labels=c("low","middle","high")))
boxplot(df1$Nb_sp~df1$category,col=3:5)

#selecting the category ranges by quantile:
xsfd<-quantile(rich_fd$FRic,c(0,1/3,2/3,1),na.rm=TRUE)
xsfd[1]<- xsfd[1]-.00005
df2 <- df1 %>% mutate(category2=cut(FRic, breaks=xsfd, labels=c("low","middle","high")))
boxplot(df2$FRic~df2$category2,col=3:5)

#CALCULATING CONGRUENCE LEVEL

#high if identical categories
high if Nb_sp == high & FRic == high or Nb_sp == middle & FRic == middle or Nb_sp == low & FRic == low

#moderate if contiguous
moderate if Nb_sp == high & FRic == middle or Nb_sp == middle & FRic == high or Nb_sp == moderate & FRic == low or Nb_sp == low & FRic == moderate

#low if extreme
low if Nb_sp == high & FRic == low or Nb_sp == low & FRic == high 


#FUNCIONAL DIvergency
loiczid <- rownames(rich_fd)
rich_fd <- cbind(loiczid=loiczid, rich_fd)

zoom_FD_map <- subs (zoom_basemap, rich_fd, by= "loiczid", which= "FDiv")

plot (zoom_FD_map, zlim=c(0.8,1), col=rainbow(5))

#FUNCIONAL richness MAP

zoom_FDric_map <- subs (zoom_basemap, rich_fd, by= "loiczid", which= "FRic")

plot (zoom_FDric_map, zlim=c(0,1))


zoom_spp2_map <- subs (zoom_basemap, rich_fd, by= "loiczid", which= "Nb_sp")
plot (zoom_spp2_map)

# look to the folder "../results/plot_FD"
##############################################################################################
# MULTIDIMENSIONAL FUNCTIONAL BETA DIVERISTY INDICES

# occurences from species weights
occ_fish_assemblages<-matrix_zoom_assemblage_ok_100
occ_fish_assemblages[which(occ_fish_assemblages>0)]<-1

# indices computation for all pairs of asssemblages, plot only on 3 first axes
install.packages ("betapart")
FbetaD_assemblages<-multidimFbetaD ( coord=coord_fish_5D,  occ=occ_fish_assemblages,  check_species_pool=FALSE, verb=TRUE,
                                     nm_asb_plot=row.names(occ_fish_assemblages) , folder_plot= paste(my_path,"/results/plot_FbetaD", sep="") , 
                                     Faxes_plot=colnames(coord_fish_5D)[1:3] , Faxes_nm_plot=colnames(coord_fish_5D)[1:3] )

# printing results
FbetaD_assemblages

# look to the folder "../results/plot_FbetaD"

##############################################################################################
# saving R objects

setwd( paste(my_path,"/results", sep="") ) # setting working directory for results

save(qual_funct_space_clear, file="qual_funct_space_clear")


save(zoom_traits_ok_NA, file="zoom_traits_ok_NA")
save(species_codes, file="species_codes")

save(coord_fish_5D, file="coord_fish_5D")
save(FD_assemblages, file="FD_assemblages")

save(weight_fish_baskets, file="weight_fish_baskets")
save(FbetaD_baskets, file="FbetaD_baskets")

##############################################################################################
# Going one step further: example of how running a null model
# here we want to test H0: Basket_4 results from a random sorting of species in the regional species pool (i.e. all species) given its species richness (SR). To keep the example simple, will focus only on FRic index so taking into account only species composition, not species weights.


# picking assemblages of SR(basket_4) species at random among the 25 fish
nbrep<-99 # number of replicates
SR_basket_4<-FD_baskets["basket_4","Nb_sp"]
SR_basket_4 # 8 species`

# empty matrix to store simulated species occurences
basket_4_H0<-matrix(0, nbrep, ncol(weight_fish_baskets), dimnames=list(1:nbrep, colnames(weight_fish_baskets) ) )
for (k in 1:nbrep)
{
  basket_4_H0[k, sample(  colnames(weight_fish_baskets), SR_basket_4) ]<-1 # random sorting of species
}# end of k

# computing FD indices on these assemblages, check_species_pool=FALSE since by chance some species could be never picked but this is not an issue
FD_basket_4_H0<-multidimFD(coord_fish_5D, basket_4_H0,  check_species_pool=FALSE )

# comparing observed and expected values under H0 using SES and p-value metrics
SES_FRic_basket_4<- (FD_baskets["basket_4","FRic"]-mean(FD_basket_4_H0[,"FRic"]) ) / sd(FD_basket_4_H0[,"FRic"])
SES_FRic_basket_4 # SES<(-1) means that observed FRic is lower than expected

pvalue_FRic_basket_4<- length(which(FD_baskets["basket_4","FRic"]<=FD_basket_4_H0[,"FRic"]))/ ( length(FD_basket_4_H0[,"FRic"]) +1 )
pvalue_FRic_basket_4 # p-value >0.975 => FRic is significantly lower than expected under H0


##############################################################################################
# END
##############################################################################################
