#################################
# Libraries used
#################################

library(shiny)
library(ggplot2)
library(vegan)
library(Rtsne)
library(limma)
library(edgeR)

################################
# Modified PERMANOVA fucntion
################################

source("parwise.adonis2_hr.r")

################################
# Additional handling functions
################################
normalize_TMM <- function(obj) {
  dge <- DGEList(counts = obj)
  dge <- calcNormFactors(dge)
  obj_norm <- cpm(dge, log = FALSE)
  return(obj_norm)
}

# Normalization function using DESeq2's median of ratios
normalize_DESeq2 <- function(obj) {
  row_medians <- apply(obj, 1, median, na.rm = TRUE)
  med_ratio <- median(row_medians, na.rm = TRUE) / row_medians
  obj_norm <- t(t(obj) * med_ratio)
  return(obj_norm)
}

# Normalization function using RPKM/FPKM method
normalize_RPKM <- function(obj) {
  total_reads <- rowSums(obj, na.rm = TRUE)
  effective_lengths <- colSums(obj > 0, na.rm = TRUE)
  obj_norm <- 1e9 * t(t(obj) / (total_reads / 1e6 * effective_lengths))
  return(obj_norm)
}


# Transformation function obtained from ASAR's source code
returnAppropriateObj <- function(obj, norm, log, method) {
  if (class(obj) != 'matrix' | any(dim(obj) == 0)) {
    return(matrix(NA, ncol = 1, nrow = 1))
  }

  if (log) {
    obj <- log2(obj + 1)
  }

  if (norm) {
    switch(method,
           'TMM' = {
             library(edgeR)
             dge <- DGEList(counts = obj)
             dge <- calcNormFactors(dge)
             obj <- cpm(dge, log = FALSE)
           },
           'DESeq2' = {
             row_medians <- apply(obj, 1, median, na.rm = TRUE)
             med_ratio <- median(row_medians, na.rm = TRUE) / row_medians
             obj <- t(t(obj) * med_ratio)
           },
           'RPKM' = {
             total_reads <- rowSums(obj, na.rm = TRUE)
             effective_lengths <- colSums(obj > 0, na.rm = TRUE)
             obj <- 1e9 * t(t(obj) / (total_reads / 1e6 * effective_lengths))
           }
    )
  }

  return(obj)
}


# Data handling function with normalization options
c_matrix <- function(a, normalize = "none") {
  ag_vec <- c("usp", "species", "genus", "family", "order", "class","phylum", "domain", "ufun", "FUN2", "FUN3","FUN4","md5")
  ag_vec <- ag_vec[ag_vec != a]
  ab_tfmat <- as.data.frame(funtaxall)
  ad_tfvec <- names(ab_tfmat) %in% ag_vec
  ad_tfmat <- ab_tfmat[,!ad_tfvec]
  ag_amat <- aggregate(. ~ ad_tfmat[,a], as.data.frame(ad_tfmat), FUN = sum)
  ag_amat <- ag_amat[, 3:ncol(ag_amat)]
  ag_amat <- as.matrix(ag_amat)

  if (normalize == "TMM") {
    libSizes <- colSums(ag_amat)
    ag_amat <- cpm(ag_amat, lib.size = libSizes)  # Apply TMM normalization
  } else if (normalize() == "DESeq2") {
    dge <- DGEList(counts = ag_amat)
    dge <- calcNormFactors(dge)  # Apply DESeq2 normalization
    ag_amat <- cpm(dge)
  } else if (normalize == "RPKM") {
    totalReads <- sum(ag_amat)
    ag_amat <- (ag_amat * 1e6) / totalReads  # Apply RPKM normalization
  } else {
    ag_amat <- returnAppropriateObj(ag_amat, norm = FALSE, log = TRUE)  # No normalization
  }

  ag_amat <- t(ag_amat)
  ag_amat <- as.data.frame(ag_amat)
  return(ag_amat)
}

#########################################
# Data loading
#########################################

load("pathview.Rdata")


#########################################
# AMVAR's Shiny UI
#########################################

ui <-fluidPage(
  titlePanel(h2("AMVAR: Automated Metagenome Visual Analysis in R", align = "center")),
  sidebarLayout(
    sidebarPanel(
      selectInput("taxfun", "Select the taxon or function to be used for aggregation:",
                  c("species", "genus", "family", "order", "class","phylum", "domain", "ufun", "FUN2", "FUN3","FUN4"), selected = "species"),
      selectInput("distance","Select the distance or dissimilarity index: ",
                  c("bray", "euclidean", "manhattan", "canberra","kulczynski","jaccard","gower","altGower","horn","mountford","raup","binomial","chao"), selected = "bray"),
      p("When euclidian distance is selected, Principal Components Analysis is computed."),
      selectInput("metadata1", "Select the first level of comparison:" ,
                  c("MG.RAST.ID","Metagenome.Name", "MGN", "Biome", "Feature", "Material", "Location", "Country", "Sequence.Type", "Sequence.Method", "Name", "EID", "MFCID", "Source","Origin"),
                  selected = "Metagenome.Name"),
      selectInput("metadata2", "Select the second level of comparison:" ,
                  c("MG.RAST.ID","Metagenome.Name", "MGN", "Biome", "Feature", "Material", "Location", "Country", "Sequence.Type", "Sequence.Method", "Name", "EID", "MFCID", "Source","Origin"),
                  selected = "Source"),
      radioButtons(inputId = "var3", label = "Select the file type", choices = list("png","pdf")),
      radioButtons("normalize", "Normalization method:",
                   choices = c("TMM", "DESeq2", "RPKM", "none"),
                   selected = "none"),
    ),
    mainPanel(
      tabsetPanel(type = "tab",
                  tabPanel("PCOA",plotOutput("pcoa"), downloadButton(outputId = "down1", label = " Download")),
                  tabPanel("T-SNE",plotOutput("tsne"), downloadButton(outputId = "down2", label = " Download")),
                  tabPanel("PERMANOVA",tableOutput("pnova"),p("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"))
      )
    )
  )
)

#######################################
# AMVAR's Shiny Computation Server
#######################################


server <-function(input,output){

  ######################################
  # Data Inputs
  ######################################
  select <- reactive({ input$taxfun })
  normalize <- reactive({
    input$normalize
  })
  # Added normalization reactive

  or_mat <- reactive({
    res <- c_matrix(select())
    # Perform normalization based on selected method
    if (normalize == "TMM") {
      # Perform TMM normalization
      res <- normalize_TMM(res)
    } else if (normalize == "DESeq2") {
      # Perform DESeq2's median of ratios normalization
      res <- normalize_DESeq2(res)
    } else if (normalize == "RPKM") {
      # Perform RPKM/FPKM normalization
      res <- normalize_RPKM(res)
    } else {
      # If no normalization is selected, use raw data as is
      # No normalization needed
    }
    return(res)
  })

  dis <- reactive({ input$distance })

  #iex <- reactive({iris[,1:4]})

  ################################################
  # METHODS COMPUTATIONS
  ################################################

  # PCOA computation
  pmat <- reactive({
    dmat <- vegdist(or_mat(), method = toString(dis()))
    pco.mat <- cmdscale(dmat, k = 2, eig = TRUE)
    return(pco.mat)
  })

  # Pairwise PERMANOVA computation
  pnova <- reactive({
    sel <- input$metadata1
    factors <- mdt[,sel]
    pnova <- pairwise.adonis(or_mat(), factors, sim.method = toString(dis()))
    return(data.frame(pnova))
  })

  # T-SNE computation
  tsne_out <- reactive({
    set.seed(42)
    dmat <- vegdist(or_mat(), method = toString(dis()))
    tsne.mat <- Rtsne(dmat, perplexity = 0.01)
    tsne.o <- as.data.frame(tsne.mat$Y)
    return(tsne.o)
  })

  #################################
  #Server Outputs
  #################################

  # PCOA
  output$pcoa <- renderPlot({
    sel1 <- input$metadata1
    level1 <- as.factor(mdt[,sel1])
    sel2 <- input$metadata2
    level2 <- as.factor(mdt[,sel2])
    points <- pmat()$points
    points <- as.data.frame(points)
    var1 <- format((round(((pmat()$GOF[1])*100),2 )), nsmall =2)
    var2 <-  format((round(((pmat()$GOF[2])*100),2 )), nsmall =2)
    ggplot(points,aes(x=points$V1,y=points$V2))+
      geom_point(aes(col = level1, shape=level2),size=3,alpha=.8)+
      labs(x = paste("Principal Coordinate 1","(",var1,"%", ")"), y = paste("Principal Coordinate 2","(",var2,"%", ")"))+
      ggtitle("Principal Coordinates Analysis")+
      theme(plot.title = element_text(hjust = 0.5, size=14, face='bold'), axis.text=element_text(size=12),axis.title=element_text(size=14))+
      guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))+
      scale_colour_discrete(name = as.character(sel1))+
      scale_shape_discrete(name = as.character(sel2))
  })

  # PERMANOVA
  output$pnova <- renderTable(pnova(), digits = 3, align = "c")

  # T-SNE
  output$tsne <- renderPlot({
    sel1 <- input$metadata1
    level1 <- as.factor(mdt[,sel1])
    sel2 <- input$metadata2
    level2 <- as.factor(mdt[,sel2])
    sel <- input$metadata1
    col <- mdt[,sel]
    ggplot(tsne_out(),aes(x=tsne_out()[,1],y=tsne_out()[,2]))+
      geom_point(aes(col = level1, shape=level2), size=3,alpha=.8)+
      labs(x = "t-SNE 1", y = "t-SNE 2")+
      ggtitle("T-distributed Stochastic Neighbor Embedding (t-SNE)")+
      theme(plot.title = element_text(hjust = 0.5, size=14, face='bold'), axis.text=element_text(size=12),axis.title=element_text(size=14))+
      guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))+scale_colour_discrete(name = as.character(sel1))+scale_shape_discrete(name = as.character(sel2))
  })


  ################################################
  # GRAPHS DOWNLOAD
  ################################################
  # Download PCOA
  plotInput1 <- reactive({
    sel1 <- input$metadata1
    level1 <- as.factor(mdt[,sel1])
    sel2 <- input$metadata2
    level2 <- as.factor(mdt[,sel2])
    points <- pmat()$points
    points <- as.data.frame(points)
    var1 <- format((round(((pmat()$GOF[1])*100),2 )), nsmall =2)
    var2 <-  format((round(((pmat()$GOF[2])*100),2 )), nsmall =2)
    ggplot(points,aes(x=points$V1,y=points$V2))+
      geom_point(aes(col = level1, shape=level2),size=3,alpha=.8)+
      labs(x = paste("Principal Coordinate 1","(",var1,"%", ")"), y = paste("Principal Coordinate 2","(",var2,"%", ")"))+
      ggtitle("Principal Coordinates Analysis") + theme(plot.title = element_text(hjust = 0.5, size=14, face='bold'), axis.text=element_text(size=12),axis.title=element_text(size=14))+
      guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) + scale_colour_discrete(name = as.character(sel1)) + scale_shape_discrete(name = as.character(sel2))
  })

  output$down1 <- downloadHandler(
    #Specify the file name
    filename = function(){
      paste("pcoa",input$var3, sep = ".")
    },
    content = function(file){
      ggsave(file, plot = plotInput1(), scale= 0.9, dpi = 400, height = 10, width = 20, units = "cm")
    })

  # Download T-SNE
  plotInput2 <- reactive({
    sel1 <- input$metadata1
    level1 <- as.factor(mdt[,sel1])
    sel2 <- input$metadata2
    level2 <- as.factor(mdt[,sel2])
    sel <- input$metadata1
    col <- mdt[,sel]
    ggplot(tsne_out(),aes(x=tsne_out()[,1],y=tsne_out()[,2]))+
      geom_point(aes(col = level1, shape=level2), size=3,alpha=.8)+
      labs(x = "t-SNE 1", y = "t-SNE 2")+
      ggtitle("T-distributed Stochastic Neighbor Embedding (t-SNE)")+
      theme(plot.title = element_text(hjust = 0.5, size=14, face='bold'), axis.text=element_text(size=12),axis.title=element_text(size=14))+
      guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))+scale_colour_discrete(name = as.character(sel1))+scale_shape_discrete(name = as.character(sel2))
  })

  output$down2 <- downloadHandler(
    #Specify the file name
    filename = function(){
      paste("tsne",input$var3, sep = ".")
    },
    content = function(file2){
      ggsave(file2, plot = plotInput2(), scale= 0.9, dpi = 400, height = 10, width = 20, units = "cm")
    })
}

shinyApp(ui = ui, server = server)
