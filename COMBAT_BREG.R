##Load packages 
library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)

##  Read the object
combat <- readRDS("20210623_1_Final_All_bcells_Combat (2).rds")
## extract metadata
combat_metadata <- combat@meta.data

## detailed cell annotations are in combat@metadata$name.unique
## To maintain uniformity change the cell types and add an annotations column 

cell_ids_detailed <- c("B.cyc","B.int.1.early.act ", "B.int.1.early.act/sw", "B.int.1.IFN.resp", "B.int.2.early.act.IFN.resp", "B.int.2.early.act/sw", "B.int.2.IFN.resp", "B.int.2.unsw", "B.mitohi", "B.NAIVE", "B.NAIVE.CD1c", "B.NAIVE.IFN.resp", "B.NAIVE.IgDlo", "B.SW.MEM", "B.SW.MEM.IFN.resp", "B.UNSW.MEM", "B.TRANSIT.CD10", "PB", "PB.cyc", "PB.IFN.resp", "PB.mitohi")
cell_ids_broad <- c("memory", "activated", "memory", "activated", "activated", "activated", "activated", "activated", "memory", "naive", "naive", "naive", "naive", "memory", "memory", "memory", "transitional", "Ab- secreting", "Ab- secreting", "Ab- secreting", "Ab- secreting")
names(cell_ids_broad)= cell_ids_detailed
cell_ids_broad
cell_annotations = combat@meta.data$name.unique
cell_annotations_broad = cell_ids_broad[cell_annotations]

##Adding new metadata column and making a new updated seurat object 
combat_updated<- AddMetaData(combat, metadata = cell_annotations_broad, col.name = "broad.annotations") ## saved as RDS
metadata_combatupdated <- combat_updated@meta.data
View(metadata_combatupdated) 

#Genes of interest
breggenes_only <- c("IL10", "MME","GZMB", "PDCD1", "CD274","SDC1", "ENTPD1", "NT5E", "EBI3","IDO1", "IL12A", "IL32", "TGFB1","THBS1", "CD1D", "CD80", "CD86", "CR1", "FASLG", "FOXP3")

## Look at avg gene expression in dataset

DotPlot(combat_updated, assay = "RNA", features = breggenes_only, cols = "RdBu", split.by = "broad.annotations",)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1,size =9) )

###Featureplot is used to visualize the expression of specific features (genes) on a plot of cells (within cells).Feature plots allow for a detailed examination of the expression patterns and heterogeneity of genes within a dataset. 

pdf("umap_allbreggenes_combat.pdf", width = 12, height = 12)  # Specify the filename and dimensions of the PDF

FeaturePlot(
  combat_updated,
  features = breggenes_only,
  reduction = "umap_allBcells"
)
?FeaturePlot
dev.off()  # Close the PDF device and save the file


##  
raw_data_combat <- combat_updated@assays$RNA@counts
genes_combat <- row.names(raw_data_combat)
setdiff(genes_combat,breggenes_only)
genes_found <- intersect(genes_combat,breggenes_only)

##  "CD1D"   "FASLG"  "IL10"   "CR1"    "SDC1"   "PDCD1"  "CD80"   "CD86"   "MME"    "IL12A"  "NT5E"   "FOXP3"  "IDO1"  
## "CD274"  "ENTPD1" "GZMB"   "THBS1"  "IL32"   "EBI3"   "TGFB1" 

exp_subset_combat<- raw_data_combat[genes_found,]

## Generate a count matrix for each gene of interest across all cells

total_counts_combat <- rowSums(exp_subset_combat)
total_counts_matrix <- as.matrix(total_counts_combat)

## Lowly expressed genes 

lowly_expressed_genes_combat <- names(which(total_counts_combat< quantile(total_counts_combat, 0.2)))
lowly_expressed_genes

## Higly expressed genes 

highly_expressed_genes_com <- names(which(total_counts_combat> quantile(total_counts_combat, 0.9))) ## "ENTPD1" "IL32"
highly_expressed_genes <- names(which(total_counts > quantile(total_counts, 0.8))) ## "CD86"   "NT5E"   "ENTPD1" "IL32"

## Plotting the lowly and highly expressed genes 
# extracting the UMAP coordinates of all cells
UMAP_combat <- combat_updated@reductions$umap_allBcells@cell.embeddings # has the umap coordinates of all droplets
rownames(UMAP_combat)<- rownames(combat_updated@meta.data)
View(UMAP_combat)

## find out the list of non-expressors

list_expressors_combat =NULL

for(g in c(1:length(genes_found))){
  w = which(exp_subset_combat[genes_found[g],]!=0)
  list_expressors_combat= c(list_expressors_combat,list(colnames(exp_subset_combat)[w]))
}
names(list_expressors_combat)= genes_found

pdf("lowly_expressed_genes_combat", )

for(g in c(1:length(genes_found))) {
  main= genes_found[g]
  plot(UMAP_combat[,1], UMAP_combat[,2], col= "grey", bg= "grey",  cex = 0.5, main = main, xlab = "UMAP1",ylab = "UMAP2")
  non_zero_expressors_combat = list_expressors_combat[[genes_found[g]]]
  points(UMAP_combat[non_zero_expressors_combat,1],UMAP_combat[non_zero_expressors_combat,2],col = "darkgreen",bg = "darkred", cex = 0.9 )
}
dev.off()

##  Make a matrix of gene expression and then construct a coexpression matrix 

genelist_new <- breggenes_only
genelist_new <- intersect(row.names(exp_subset), genelist_new)

# Creating null matrices
combat_coexpression= matrix(data = 0, nrow = length(genelist_new), ncol = length(genelist_new), dimnames = c(list(genelist_new), list(genelist_new)))
combat_coexpression_jaccard = matrix(data = 0, nrow = length(genelist_new), ncol = length(genelist_new), dimnames = c(list(genelist_new), list(genelist_new)))


# Iterate the loop to make  matrix for gene expression 

for(g1 in c(1:length(genelist_new))){
  w1= which(exp_subset_combat[genelist_new[g1],]!=0)
  for(g2 in c(1: length(genelist_new))){
    w2= which(exp_subset_combat[genelist_new[g2],]!=0)
    combat_coexpression[g1,g2]= length(intersect(w1,w2))
    if (length(w1) > 0 & length(w2) > 0) {
      combat_coexpression_jaccard[g1, g2] <- length(intersect(w1, w2)) / length(union(w1, w2))
    }
  }
}

saveRDS(combat_coexpression, "coexpression_matrix_combat.RDS")
saveRDS(combat_coexpression_jaccard, "coexpression_jaccardmatrix_combat.RDS")



## Plotting the jaccard indices as heatmap
combat_coexpression_jaccard
diag(combat_coexpression_jaccard) <- NA
heatmap.2(combat_coexpression_jaccard,
          trace = "none",  # Remove row and column names
          col = colorRampPalette(c("orange", "red","darkred")),  # Define the color gradient
          main = "Coexpression Jaccard Heatmap",
          key = TRUE,  # Display the color scale legend
)

## Plotting  a UMAP of all Breg genes overlaying clusters

add.alpha <- function(col, alpha=1){
  if(missing(col)) # the function checks for missing vectors in col first and if there is missing values it stops
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, # else the function converts each of the vectors into a RGB color with  col2rgb 
        #and then dvides by 255 to make them between 0 1nd 1 then the apply function takes each column in the RGB matrix
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  # function takes rgb values in column and coverts it to a color with transparency 
}


cols1 =  c(add.alpha (brewer.pal(6, "Dark2"), alpha = 0.95),add.alpha (brewer.pal(8, "Paired"), alpha = 0.95), add.alpha (brewer.pal(8, "Set1"), alpha = 0.95))

# add trabsperancy to vectors but change transparency to 0.5
cols =  add.alpha (cols1, alpha = 0.5)
#extact the rownames of metadata- droplet ids
id <- rownames(combat_updated@meta.data)
# Create a color palette for each  of the droplets , i.e along the length of the id 
# initialize with grey color 
col_cell <- rep("grey",length(id))

# Assign dropletids as names to each element in the color cell vector 

names(col_cell) <- id
genelist <- c( "IL10", "MME","GZMB", "PDCD1", "CD274","SDC1", "ENTPD1", "NT5E","IDO1", "EB13", "IL32", "TGFB1","THBS1", "CD1D", "CD80", "CD86", "CR1", "FASLG", "FOXP3")
raw_combat <- combat_updated@assays$RNA@counts
genes <- rownames(raw_combat )
setdiff(genes, breggenes_only)
genelist <- intersect(genes,breggenes_only)
m_raw =t(as.data.frame(raw_combat[genelist, ]))

# getting the umap coordinates for each droplet
# extracting the UMAP coordinates of all cells

UMAP_combat <- combat_updated@reductions$umap_allBcells@cell.embeddings
UMAPcoord_combat <- as.data.frame(UMAP_combat)

for(i in c(1:length(genelist))){
  w = which(m_raw[,genelist[i]]!=0)
  col_cell[w] = cols1[i]  # cols1[i] retrieves the color value at index i from the cols1 vector.
}

# applying a function to calculate the number of non zero values forneach droplet ( each row)

multiple = apply(m_raw, 1, function(x){length(which(x!=0))})

# find entries / indices in multiple with values greater than one~ expressing more than one gene
multiple = which(multiple>1)
# color them black 
col_cell[multiple] = "black"

fileout1="C:/Users/x/Desktop/10x/umapgbregene_combat_file_new.pdf"
w <- 3.5
pdf(file=fileout1, height <-w*1.5, width=w*3)
par(mfrow= c(1,2), mar= c(4,4,2,2))
size <- 0.05
nonzero <- which(col_cell != "grey")
plot(UMAPcoord_combat[nonzero,], pch = 21, xlab = "UMAP1", ylab = "UMAP2", main = "Expression of Breg genes", cex = size, col = col_cell[nonzero], bg = col_cell[nonzero])
#w = which(col_cell != "grey")

#points(xy[w,], pch = 21, cex = size, col = col_cell[w], bg = col_cell[w])

cex = 0.9 # controls size and text of plots 
plot(c(1, 2), c(-2.5, 2), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1, 
     cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = "", 
     ylim = c(-2.5,2),xlim = c(0.5, 2), axes = F)

col_leg = c("grey", cols1[c(1:length(genelist))],"black")
leg = c("no Breg genes", paste(genelist, "+", sep = ""), "multiple Breg genes")
legend("topleft", leg, pch = 22,cex= 0.6, bty="n", pt.cex = 1.5, pt.bg = col_leg, 
       col = col_leg, pt.lwd = 1, text.font = 1)
dev.off()

### Plotting percent expression of the genes of interest 
#### Plotting the percent expression fo rgenes of interest 

genes_of_interest <- breggenes_only
rna_data <- GetAssayData(combat_updated, assay = "RNA")
goi_v <- genes_of_interest
rna_data <- subset(rna_data, )

for(i in 1:length(goi_v)){
  goi<-goi_v[i]
  rna_gene <- rna_data[goi,]
  rna_gene_df <- as.data.frame(rna_gene)
  matrices_goi <- cbind(rna_gene_df, combat_updated@meta.data$scRNASeq_sample_ID, combat_updated@meta.data$Source,combat_updated@meta.data$broad.annotations)
  matrices_goi <- matrices_goi %>% rownames_to_column(var= "cell_id")
  colnames(matrices_goi)<- c("cell_id", paste0("counts_", goi),"patient_id","covid_status", "celltype")
  matrices_goi <-matrices_goi %>% filter(covid_status!="Flu")
  matrices_goi <- matrices_goi %>% filter(celltype!= "NA") 
  matrices_goi$covid_status[matrices_goi$covid_status == "COVID_HCW_MILD"] <- "COVID_MILD"
  
  ## getting the total counts per patient per cell type
  total_counts_goi <- matrices_goi %>%
    group_by(patient_id, celltype) %>% count()
  
  ## Append total counts to the matrix of gene 
  matrice_test<- matrices_goi%>% left_join(total_counts_goi, by= c("patient_id","celltype"))
  # store them in a new column
  data.table::setnames(matrice_test, "n", "percelltypecount")
  # calculate the percentages
  tmp_goi<- matrice_test %>% group_by(patient_id,celltype) %>%  summarise(across(starts_with("counts_"), ~sum(.x!=0)))
  tmp_goi<- matrice_test %>% left_join(tmp_goi, by= c("patient_id","celltype"))
  colnames(tmp_goi)[7] <- paste0("total_", goi)
  tmp_goi$percentage <- (tmp_goi[, paste0("total_", goi)] / tmp_goi$percelltypecount) * 100
  
  # Group and summarize
  test_goi <- tmp_goi %>% group_by(patient_id, celltype, percentage, covid_status) %>% summarise() 
  # Create the plot
  plot_goi <- ggplot(test_goi, aes(x = factor(celltype), y = percentage, fill = covid_status)) +
    geom_boxplot() +
    labs(x = "Cell Type", y = "Percent Expression", title = paste( "% of cells expressing", goi))
  
  # Save plot as PDF
  pdf(paste0(goi, "_plot.pdf"))
  print(plot_goi)
  dev.off() 
  
}
##ANOVA
genes_of_interest <- breggenes_only
rna_data <- GetAssayData(combat_updated, assay = "RNA")
goi_v <- genes_of_interest
rna_data <- subset(rna_data, )

# List to store the ANOVA results for all genes
all_anova_results <- list()

# Loop over each gene of interest

for (i in 1:length(goi_v)){
  goi<-goi_v[i]
  rna_gene <- rna_data[goi,]
  rna_gene_df <- as.data.frame(rna_gene)
  matrices_goi <- cbind(rna_gene_df, combat_updated@meta.data$scRNASeq_sample_ID, combat_updated@meta.data$Source, combat_updated@meta.data$broad.annotations)
  matrices_goi <- matrices_goi %>% rownames_to_column(var = "cell_id")
  colnames(matrices_goi) <- c("cell_id", paste0("counts_", goi), "patient_id", "covid_status", "celltype")
  matrices_goi <- matrices_goi %>% filter(covid_status != "Flu")
  matrices_goi <- matrices_goi %>% filter(!is.na(celltype))
  matrices_goi$covid_status[matrices_goi$covid_status == "COVID_HCW_MILD"] <- "COVID_MILD"
  
  # Getting the total counts per patient per cell type
  total_counts_goi <- matrices_goi %>%
    group_by(patient_id, celltype) %>% count()
  
  # Append total counts to the matrix of gene
  matrice_test <- matrices_goi %>% left_join(total_counts_goi, by = c("patient_id", "celltype"))
  data.table::setnames(matrice_test, "n", "percelltypecount")
  
  # Calculate the percentages
  tmp_goi <- matrice_test %>% group_by(patient_id, celltype) %>% summarise(across(starts_with("counts_"), ~sum(.x != 0)))
  tmp_goi <- matrice_test %>% left_join(tmp_goi, by = c("patient_id", "celltype"))
  colnames(tmp_goi)[7] <- paste0("total_", goi)
  tmp_goi$percentage <- (tmp_goi[, paste0("total_", goi)] / tmp_goi$percelltypecount) * 100
  
  # Group and summarize
  test_goi <- tmp_goi %>% group_by(patient_id, celltype, percentage, covid_status) %>% summarise()
  
  # Perform one-way ANOVA and store the results for each cell type separately
  anova_result <- perform_anova(goi, test_goi)
  if (!is.null(anova_result)) {
    all_anova_results[[goi]] <- anova_result
  }
}

# View the ANOVA results for all genes and each cell type separately
all_anova_results

#
output_directory <-  "C:/Users/x/Desktop/10x/combat/results"
##
# Loop over each gene and save the ANOVA results for each gene as a separate CSV file
for (goi in names(all_anova_results)) {
  gene_anova_results <- all_anova_results[[goi]]
  
  # Create a file name for each gene's CSV file
  file_name <- paste0(goi, "_anova_results.csv")
  
  # Set the full path of the file
  file_path <- file.path(output_directory, file_name)
  
  # Convert the list to a data frame
  gene_anova_df <- do.call(rbind, gene_anova_results)
  
  # Save the data frame as a CSV file
  write.csv(gene_anova_df, file = file_path, row.names = FALSE)
}
#### trying to plot
# Load the purrr library
library(purrr)

# Flatten the nested list to a single-level list
flat_anova_results <- flatten(all_anova_results)

# Create the data frame
p_values_df <- data.frame(
  Gene = unlist(lapply(flat_anova_results, `[[`, "Gene")),
  CellType = unlist(lapply(flat_anova_results, `[[`, "CellType")),
  PValue = unlist(lapply(flat_anova_results, `[[`, "ANOVA_PValue"))
)

# Print the resulting data frame
print(p_values_df)



#### Differential gene expression 
test1$covid_status[test1$covid_status== "COVID_HCW_MILD"]<- "COVID_MILD"
Idents(test1)<-test1@meta.data[["covid_status"]]
test1$celltype.condition<-paste(test1@meta.data[["broad.annotations"]], Idents(test1), sep="_")
vec<-unique(test1@meta.data[["broad.annotations"]])
Idents(test1)<-test1$celltype.condition
w<-unique(test1@meta.data[["covid_status"]])

output<-"C:/Users/x/Desktop/10x/combat"

for(i in 1:length(vec))
{
  ident1<-paste0(vec[i],"_HV")
  ident2<-paste0(vec[i],"_COVID_CRIT")
  
  condition.diffgenes<-FindMarkers(test1, ident.1=ident1, ident.2=ident2, min.pct=0.1, logfc.threshold = 0.25, features =  c("IL10", "MME","GZMB", "PDCD1", "CD274","SDC1", "ENTPD1", "NT5E", "EBI3","IDO1", "IL12A", "IL32", "TGFB1","THBS1", "CD1D", "CD80", "CD86", "CR1", "FASLG", "FOXP3"))
  write.csv(condition.diffgenes, file=paste0(vec[i],".csv"))
}

for(m in 1:length(vec))
{
  for(i in 1:length(w))
  {
    for( j in 1:length(w)-i)
    {
      ident1<-paste0(vec[m],"_", w[i])
      ident2<-paste0(vec[m],"_", w[i+j])
      condition.diffgenes<-FindMarkers(test1, ident.1=ident1, ident.2=ident2, min.pct=0.1, logfc.threshold = 0.25, features =  c("IL10", "MME","GZMB", "PDCD1", "CD274","SDC1", "ENTPD1", "NT5E", "EBI3","IDO1", "IL12A", "IL32", "TGFB1","THBS1", "CD1D", "CD80", "CD86", "CR1", "FASLG", "FOXP3"))
      write.csv(condition.diffgenes, file=paste0(vec[m],"_", w[i], "_vs_", w[i+j],".csv"))
    }
  }
}

## Association between pandisease and combat datasets

## library(psych)

proportion_per_celltypenew <- read_excel("percent-expression-matrix_new.xlsx")
proportion_per_celltypenew <- proportion_per_celltypenew[,-9:12] ## empty columns
#rename columns
colnames(proportion_per_celltypenew)<- c("genes", "ag-presenting",  "GC/homing",  "exhausted", "activated", "memory", "naive","ab-secreting", "naive_COMBAT", "ab-secreting_COMBAT", "memory_COMBAT", "transitional_COMBAT", "activated_COMBAT" )

#Absecreting cell
#Fit-linear model
fitab <- lm(proportion_per_celltypenew$`ab-secreting`~proportion_per_celltypenew$`ab-secreting_COMBAT`)
p_value_ab <- summary.lm(fitab)$coefficients[2, 4]
tmp_ab <- cbind(proportion_per_celltype$`ab-secreting`, proportion_per_celltype$`ab-secreting_COMBAT`)

## Correlation and correction
cortest_ab= corr.test(tmp_ab,adjust="holm")
pval_ab= cortest_ab$p[1,2]
rval_ab= cortest_ab$r[1,2]
## repeat for all cell types

##Plot
ggplot(proportion_per_celltypenew, aes(x = ab-secreting, y = ab-secreting_COMBAT)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "activated B cells in Pan-disease dataset", y = "activated B cells in COMBAT dataset") +
  ggtitle(" ") +
  stat_cor(method = "pearson", label.x = 0.5, label.y = 4, label.sep = "\n",
           vjust = -1.5, hjust = 0.7)