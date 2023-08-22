maintitle <- ""

titletab1 <- "F/T"
titletab2 <- "F/M"
titletab3 <- "T/M"
titletab4 <- "KEGG Pathways Heatmap"
titletab5 <- "KEGG Pathway"
titletab6 <- "Settings"
titletab7 <- "Metadata"
titleForEmtyDataPopup <- "Empty selection!"
textForEmptyDataPopup <- "Please select some function or taxon."

titleForDimErrorPopup <- "Wrong Heatmap dimensions!"
textForDimErrorPopup <- "Please select other function or taxon or go one level higher."
titleForNRowErrorPopup <- "Wrong number of rows!"
textForNRowErrorPopup <- "Number of rows suppose to be integer."
titleForNonMatchKostat <- "Error with annotation files!"
textForNonMatchKostat <- "There is inconsistency within KEGG Orthology Annotation file."
titleForNonMatchZerosKostat <- "There is no genes in selected KEGG  pathway."
textForNonMatchZerosKostat <- 'Please find another pathway or add more metagenomes.'
titleForSavingMetadata <- "Saving your metadata..."
textForSavingMetadata <- "Your changes will be applied after you restart the app. Saving process may take a few seconds. Please, restart the app!"
textForDownloadingMetadata <- "KEGG Map is saved in your working directory."
ColNameSelectorTitle <- "Please select a column of metadata to use as metagenome names"
metagenomeone <-"Select multiple metagenomes"
metagenome1selected <- c(1,2)
metagenometwo <- "Select metagenome"
metagenome2selected <- 1

plot1Title <-'F vs. T heatmap'
plot2Title <-'F vs. M heatmap'
plot3Title <-'T vs. M heatmap'
plot4Title <-'KEGG pathway heatmap'


taxone <- "Taxonomy Selection level"
tax1n <- c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain", "root" = "toplevel")
tax1selected <- "genus"
taxthree <- "Select taxons"
taxtwo <- "Taxonomy Aggregation level"
tax2n <- c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain")
tax2selected <- "usp"

funcone <- "Function Selection level"
func1n <- c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4", "root" = "toplevel")
func1selected <- "FUN4"
functhree <- "Select function"
functwo <- "Function Aggregation level"
func2n <- c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4")
func2selected <- "ufun"

truncone <- "Select Truncation Function"
pruneone <- "Select Pruning Function"
prunec1n <- c('sum'='sum','max'='max','min'='min','sd'='sd')
prunec1selected <- 'sd'

pathwayone <- "Input Pathway ID"

pathwayfile<- '../../pathview.RData'
#colname for sample name pickup.
#Color Palette for Heatmaps
#Settings
labelColPal <- "Select color palette for heatmaps"
currentPalette <- "Blues"
colName <- 'Metagenome.Name'#"Group"
if(file.exists("Settings.Rdata")){
  load("Settings.Rdata")
}

# if(file.exists('~/Documents/Projects/R/MFCMAP/R/pathview.Rdata')){
#   load('~/Documents/Projects/R/MFCMAP/R/pathview.Rdata')
# }
if(file.exists(pathwayfile)){
  load(pathwayfile)
}
