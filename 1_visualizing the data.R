#--------------------------------------------------------------------------------------------------------#
#LOADING REQUIRED LIBRARIES 
#--------------------------------------------------------------------------------------------------------#
library (tidyverse) #for data cleaning and grouping
library (ggplot2) #for ploting
library (skimr) #for summarizing 
library(VIM) #impute missing values
#install.packages ("ggcorrplot")
library(ggcorrplot)#for ploting correlation graphs
#--------------------------------------------------------------------------------------------------------#
#READING THE TWO DIFFERENT DATASETS-
#--------------------------------------------------------------------------------------------------------#

##1) presence-absence matrix

species_hs_cc_50_random_spread <- read.csv('data/presence_absence_matrix_cc_50_random.csv')

#rename the colnames
names(species_hs_cc_50_random_spread) <- gsub(x = names(species_hs_cc_50_random_spread),
                                              pattern = "\\.",
                                              replacement = "-")
#verify is changes have been correctly done
as.tibble (species_hs_cc_50_random_spread)

#deleting the X column that is created when we read the csv
species_hs_cc_50_random_spread <- species_hs_cc_50_random_spread[-1]
as.tibble (species_hs_cc_50_random_spread)
str (species_hs_cc_50_random_spread)

##2)Traits matrix
#we assign am_sid to row names while reading
traits_sp_cc_50_random<-read.csv("data/traits_matrix_cc_50_random.csv", header=T, row.names = "am_sid")

head (traits_sp_cc_50_random)
str (traits_sp_cc_50_random)

#deleting the X column that is created when we read the csv
traits_sp_cc_50_random <- traits_sp_cc_50_random[-1]
as.tibble (traits_sp_cc_50_random)
head (traits_sp_cc_50_random)

#Loo = Asymptotic maximum length (cm) 
#K = individual growth rate (yr^-1) brody growth coefficient
#Winfinity = asymptotic maximum mass (gr) 
#tmax = maximum age (yr) 
#tm = age at maturation (yr)
#M = natural mortality rate (yr^-1)
#Lm = length at maturation (cm) 

#--------------------------------------------------------------------------------------------------------#
#UNDERSTANDING AND VISUALIZING THE DATA
#--------------------------------------------------------------------------------------------------------#

#count up by families
counts_families <- traits_sp_cc_50_random %>% 
  group_by(family) %>% 
  count()

counts_families$family <- factor(counts_families$family, levels = counts_families$family[order(counts_families$n)])

family_plot <- ggplot(counts_families, aes(family, n)) + geom_bar(stat = "identity", fill = "forestgreen", width = 0.25) + coord_flip()
family_plot

#count number of species by loiczid
spp_byloiczid <- species_hs_cc_50_random_spread %>% 
  mutate (spp_byloiczid =rowSums (species_hs_cc_50_random_spread[ , 2:51]))
spp_byloiczid

#--------------------------------------------------------------------------------------------------------#
#MISSING VALUES
#--------------------------------------------------------------------------------------------------------#
#what is the proportion of missing data for each variable?
pctmiss <- colSums(is.na(traits_sp_cc_50_random))/nrow(traits_sp_cc_50_random)
                  round(pctmiss, 2)
pctmiss

#Feeding type and importance have a high percentage of missing values, bodyshape and pricecateg has some missing value too

#let's visualize the distribution of the data using histograms

#we create new variable per family
myctophidae <- subset(traits_sp_cc_50_random, family =="Myctophidae", ) 
myctophidae

ggplot(myctophidae, aes(x = depthrangedeep)) +
  geom_histogram(fill = "cornflowerblue", 
                 color = "white") 
  
#let's calculate the mean and the median
depthrange_by_family<- traits_sp_cc_50_random %>% 
  group_by(family)  %>%  
  summarise (meandeep = mean (depthrangedeep)) %>% 
  summarise (medianLoo = median (depthrangedeep))


#we can impute missing values using the 5 nearest neighbors. There are other packages that can be used:
#VIM, mice, Amelia and missForest
#install.packages ("VIM")
library(VIM)
traits_noNA <- kNN(traits_sp_cc_50_random, k=5)

pctmiss_traits_noNA <- colSums(is.na(traits_noNA))/nrow(traits_noNA)
round(pctmiss_traits_noNA, 2)

#or we can remove NAs
traits_doppedna <- traits_sp_cc_50_random %>% 
  drop_na()

unique (traits_sp_cc_50_random$am_sid) %>% length () #50
unique (traits_doppedna$am_sid) %>% length () #25

#we create a df with the species id 
spp_noNA <- select (traits_doppedna, am_sid)

list_spp_noNA <- c ("Fis-23277", "Fis-23106", "Fis-23149", "Fis-30362", "Fis-23185", "Fis-22715", "Fis-29313", "Fis-29338", "Fis-22752", "Fis-25410", "Fis-25577", "Fis-24789", "Fis-25172", "Fis-23173", "Fis-22852", "Fis-23580", "Fis-22747", "
Fis-6123", "Fis-22768", "Fis-23826", "Fis-23995", "Fis-23407", "Fis-30758", "Fis-29753", "Fis-30477")

#we select only the values from those 25 species from the spread df
head (species_hs_cc_50_random_spread)
species_hs_cc_random_spread_noNA <-species_hs_cc_50_random_spread[ , names(species_hs_cc_50_random_spread) %in% list_spp_noNA, drop=F]

str (species_hs_cc_random_spread_noNA)

#--------------------------------------------------------------------------------------------------------#
#CORRELATION
#--------------------------------------------------------------------------------------------------------#

# select numeric variables
df <- dplyr::select_if(traits_sp_cc_50_random, is.numeric)
# select meanfull traits
df <- df %>% 
  select(Loo, K, Winfinity, tmax, tm, M, Lm, depthrangedeep)


# calulate the correlations
r <- cor(df, use="complete.obs")
round(r,2)

#plot
ggcorrplot(r, 
           hc.order = TRUE, 
           type = "lower",
           lab = TRUE)

#--------------------------------------------------------------------------------------------------------#
#DATA CLEANING 
#--------------------------------------------------------------------------------------------------------#
traits_sp_cc_50_random<-read.csv("data/traits_matrix_cc_50_random.csv", header=T)
#deleting the X column that is created when we read the csv
traits_sp_cc_50_random <- traits_sp_cc_50_random[-1]
as.tibble (traits_sp_cc_50_random)
head (traits_sp_cc_50_random)


traits <- traits_sp_cc_50_random %>% 
  select(-class, -order, -family, - EnvTemp, -IUCN_Code, -Resilience, -pd50, -vulnerability, -importance, -pricecateg)
head(traits)
str (traits)

table (traits$bodyshape)
traits[,"bodyshape"]<-as.character(traits [,"bodyshape"])
#traits[which(traits[,"bodyshape"]=="other (see remarks)"),"bodyshape"]<-"other"
traits[which(traits[,"bodyshape"]=="eel-like"),"bodyshape"]<-"eel_like"
traits[which(traits[,"bodyshape"]=="elongated"),"bodyshape"]<-"elongated"
traits[which(traits[,"bodyshape"]=="fusiform / normal"),"bodyshape"]<-"fusiform_normal"
traits[which(traits[,"bodyshape"]=="short and / or deep"),"bodyshape"]<-"short_andor_deep"
traits[which(traits[,"bodyshape"]=="other"),"bodyshape"]<-"other"

traits[,"bodyshape"]<-factor(traits[,"bodyshape"], levels=c("eel_like", "elongated", "fusiform_normal", "short_andor_deep", "other" ), ordered = TRUE )
table (traits$bodyshape)

#--------------------------------------------------------------------------------------------------------#
#GROUPING
#--------------------------------------------------------------------------------------------------------#

#install.packages ("skimr")
library (skimr)

skim (traits)

spp_list<-read.csv("data/spp_list.csv", header=T)
#deleting the X column that is created when we read the csv
spp_list <- spp_list[-1]

traits <- left_join(traits, spp_list, by="am_sid")
traits <- traits %>% 
  select(am_sid, sciname, Loo, K, Winfinity, tmax, tm, M, Lm, depthrangedeep, demerspelag, FeedingType, bodyshape)

head(traits)

summary (traits$Loo) #Loo = Asymptotic maximum length (cm) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#6.148  18.729  30.635  74.966  60.103 899.551 
summary (traits$K) #K = individual growth rate (yr^-1) brody growth coefficient
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0525  0.1830  0.3502  0.4337  0.5990  2.0514
summary (traits$Winfinity) #Winfinity = asymptotic maximum mass (gr) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2      60     290  112140    2789 5159534 
summary (traits$tmax)#tmax = maximum age (yr) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.260   4.474   8.287  10.935  14.464  54.552
summary (traits$tm)#tm = age at maturation (yr)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5989  1.3167  2.3217  3.4758  4.1229 21.5705 
summary (traits$M)#M = natural mortality rate (yr^-1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0732  0.3266  0.5762  0.7746  1.0073  2.8745 
summary (traits$Lm)#Lm = length at maturation (cm) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.421  10.997  18.087  41.460  33.325 491.172 
summary (traits$depthrangedeep)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   12.0   296.5   500.0   939.4  1192.5  4938.0       3 

traits_grouped <- traits %>% 
  mutate (Loo_class= cut(Loo, c(0,25,50,100,1500), labels=c("short", "medium", "long", "very_long"))) %>% 
  mutate (K_class= cut(K, c(0,0.5,1,1.5,2.5), labels=c("low_growth", "medium_growth", "high_growth", "very_high_growth"))) %>%
  mutate (W_class= cut(Winfinity, c(0,100,1000,10000,6000000), labels=c("light_weigth", "medium_weigth", "heavy_weigth", "very_heavy_weigth"))) %>%
  mutate (age_class= cut(tmax, c(0,5,15,25,55), labels=c("short_life", "medium_life", "long_life", "verylong_life"))) %>% 
  mutate (mature_class= cut(tm, c(0,1,5,10,25), labels=c("fast_mature", "medium_mature", "long_mature", "verylong_mature"))) %>%
  mutate (mortality_class= cut(M, c(0,0.5,1,2,3), labels=c("low_mortal", "medium_mortal", "high_mortal", "veryhigh_mortal"))) %>%
  mutate (maturelength_class=cut(Lm, c(0,10,20,100,500), labels=c("short_maturelength", "medium_maturelength", "long_maturelength", "verylong_maturelength"))) %>%
  mutate (deeprange_class=cut(depthrangedeep, c(0,250,500,2500,5000), labels=c("shallow", "medium_deep", "deep", "very_deep")))

functional_groups <- traits_grouped %>% 
  group_by(Loo_class,K_class, W_class, age_class, mature_class, mortality_class, maturelength_class) %>% 
  summarise(n=n())


#--------------------------------------------------------------------------------------------------------#
#FUNCTIONAL DIVERSITY CALCULATIONS
#--------------------------------------------------------------------------------------------------------#

# Before applying functional diversity functions we need to name the rows of presence-abscence matrix with the loiczid value
species_hs_cc_50_random_spread <- column_to_rownames(species_hs_cc_50_random_spread, var = "loiczid")

#and check that species names are the same in the two matrices
#sum(row.names(TRAITS DATA SET NAME) %in% colnames(species_hs_cc_50_random_spread)) == ncol(species_hs_cc_50_random_spread)