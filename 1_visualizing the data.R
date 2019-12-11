#READING THE TWO DIFFERENT DATASETS----------------------------------------------------------------------------------#

##1) presence-absence matrix
species_hs_cc_50_random_spread <- read.csv('/home/valle/github/functional_diversity_training/data/presence_absence_matrix_cc_50_random.csv')

#rename the colnames
names(species_hs_cc_50_random_spread) <- gsub(x = names(species_hs_cc_50_random_spread),
                                              pattern = "\\.",
                                              replacement = "-")
#verify is changes have been correctly done
as.tibble (species_hs_cc_50_random_spread)

#we name the rows with the loiczid value
species_hs_cc_50_random_spread <- column_to_rownames(species_hs_cc_50_random_spread, var = "loiczid")

#verify is changes have been correctly done
as.tibble (species_hs_cc_50_random_spread)

#deleting the X column that is created when we read the csv
species_hs_cc_50_random_spread <- species_hs_cc_50_random_spread[-1]
as.tibble (species_hs_cc_50_random_spread)
str (species_hs_cc_50_random_spread)

##2)Traits matrix
#we assign am_sid to row names while reading
traits_sp_cc_50_random<-read.csv("/home/valle/github/functional_diversity_training/data/traits_matrix_cc_50_random.csv", header=T, row.names = "am_sid")

head (traits_sp_cc_50_random)
str (traits_sp_cc_50_random)
#deleting the X column that is created when we read the csv
traits_sp_cc_50_random <- traits_sp_cc_50_random[-1]
as.tibble (traits_sp_cc_50_random)

## checking that species names are the same in the two matrices
sum(row.names(traits_sp_cc_50_random) %in% colnames(species_hs_cc_50_random_spread)) == ncol(species_hs_cc_50_random_spread)

#VISUALIZING THE DATA-------------------------------------------------------------------------------------------------#

#count up by families
counts_families <- traits_sp_cc_50_random %>% group_by(family) %>% count()

counts_families$family <- factor(counts_families$family, levels = counts_families$family[order(counts_families$n)])

family_plot <- ggplot(counts_families, aes(family, n)) + geom_bar(stat = "identity", fill = "forestgreen", width = 0.25) + coord_flip()
family_plot

#number of species per loiczid

#we keep rowname as a column
library(data.table)
setDT(species_hs_cc_50_random_spread, keep.rownames = TRUE)[]

#we rename the column to loiczid
species_hs_cc_50_random_spread <- species_hs_cc_50_random_spread %>% 
  rename(loiczid =rn)

#we sum the rows
spp_byloiczid <- rowSums (species_hs_cc_50_random_spread[ , 2:51])
spp_byloiczid



