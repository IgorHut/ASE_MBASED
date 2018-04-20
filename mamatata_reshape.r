library(tidyverse)
load("mamatata.rda")
names(mamatata)
mamatata$geneName

mamatata %>% mutate_if(is.factor, as.character) -> mamatata
nchar(mamatata$geneName[1])

mamatata_sep <- filter(mamatata, nchar(mamatata$geneName) > 15)
head(mamatata_sep2)

mamatata_sep = separate(mamatata_sep, geneName, 
                         into = c("geneName1", "geneName2"), sep = ';')

head(mamatata_sep)

mamatata_sep1 <- select(mamatata_sep, - geneName2)
head(mamatata_sep1)
names(mamatata_sep1)[names(mamatata_sep1) == "geneName1"] <- 'geneName'
head(mamatata_sep1)

mamatata_sep2 <- select(mamatata_sep, - geneName1)
head(mamatata_sep2)
names(mamatata_sep2)[names(mamatata_sep2) == "geneName2"] <- 'geneName'
head(mamatata_sep2)

mamatata_0 <- filter(mamatata, nchar(mamatata$geneName) == 15)
head(mamatata_0)

mamatata_final <- rbind(mamatata_0, mamatata_sep1, mamatata_sep2)
head(mamatata_final)

nrow(mamatat_final)
nrow(mamatata)
