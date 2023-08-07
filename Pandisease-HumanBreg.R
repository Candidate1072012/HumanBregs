


#load packages

library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(ggplot2)

##setting working directory 

## read the object 
bcell_combined <- readRDS("Bcellcombined_obj_10_21.rds")
metadata_bcells <- bcell_combined[[]]

## list of Breg Genes 

breggenes_only <- c("IL10", "MME","GZMB", "PDCD1", "CD274","SDC1", "ENTPD1", "NT5E", "EBI3","IDO1", "IL12A", "IL32", "TGFB1","THBS1", "CD1D", "CD80", "CD86", "CR1", "FASLG", "FOXP3")



###Plotting the Breg genes_genes

## Dotplot to see average expression of B regulatory genes across cell types. The DotPlot function in Seurat is used to visualize the expression of specific features (genes) across different groups or conditions in a dataset.
##It creates a dot plot where each dot represents a cell, and the color and size of the dot indicate the expression level of the feature,hughlights differentially expressed genes;

DotPlot(bcell_combined, assay = "RNA", features = breggenes_only, cols = "RdBu",)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1,size =9))


## Featureplotis used to visualize the expression of specific features (genes) on a plot of cells (within cells).Feature plots allow for a detailed examination of the expression patterns and heterogeneity of genes within a dataset. 
pdf("umap_allbreggenes", width = 12, height = 12 )
FeaturePlot(bcell_combined, features = breggenes_only, reduction = "umap")
dev.off()

##  UMAP of cell annotations.... plot clusters of cells present 
DimPlot(bcell_combined, group.by = "broad_annotations", reduction = 'umap', label = 'TRUE')+
  ggtitle("B cell subsets in the Pan-disease B cell dataset")

## Identifying and plotting lowly expressed genes

# Extract the data for the genes of interest 

raw_data <- bcell_combined@assays$RNA@counts
genes <- row.names(raw_data)
setdiff(genes,breggenes_only)
genes_found <- intersect(genes,breggenes_only)

##  "CD1D"   "FASLG"  "IL10"   "CR1"    "SDC1"   "PDCD1"  "CD80"   "CD86"   "MME"    "IL12A"  "NT5E"   "FOXP3"  "IDO1"  
## "CD274"  "ENTPD1" "GZMB"   "THBS1"  "IL32"   "EBI3"   "TGFB1" 

exp_subset <- raw_data[genes_found,]

## Generate a count matrix of total counts each gene of interest across all cells

total_counts <- rowSums(exp_subset)
total_counts_mat <- as.matrix(total_counts )

## Lowly expressed genes 

lowly_expressed_genes <- names(which(total_counts< quantile(total_counts, 0.4)))
lowly_expressed_genes

## Higly expressed genes 

highly_expressed_genes <- names(which(total_counts > quantile(total_counts, 0.9))) ## SDC1", "TGFB1
highly_expressed_genes <- names(which(total_counts > quantile(total_counts, 0.8))) ## "SDC1"   "ENTPD1" "IL32"   "TGFB1" 

## Plotting the lowly and highly expressed genes 

# extracting the UMAP coordinates of all cells
UMAP <- bcell_combined@reductions$umap@cell.embeddings # has the umap coordinates of all droplets
rownames(UMAP)<- rownames(bcell_combined@meta.data)
View(UMAP)


## 

## find out the list of non-expressors
list_expressors =NULL

for(g in c(1:length(genes_found))){
  w = which(exp_subset[genes_found[g],]!=0)
  list_expressors= c(list_expressors,list(colnames(exp_subset)[w]))
}
names(list_expressors)= genes_found

## Plotting lowly expressed_ non-zero expressors 

pdf("plots_lowlyexpressed_genes_new", )

for(g in c(1:length(genes_found))) {
  main= genes_found[g]
  plot(UMAP[,1], UMAP[,2], col= "grey", bg= "grey",  cex = 0.5, main = main, xlab = "UMAP1",ylab = "UMAP2")
  non_zero_expressors = list_expressors[[genes_found[g]]]
  points(UMAP[non_zero_expressors,1], UMAP[non_zero_expressors,2],col = "darkred",bg = "darkred", cex = 0.9 )
}
dev.off()

##  Make a matrix of gene expression and then construct a coexpression matrix 

genelist_new <- breggenes_only
genelist_new <- intersect(row.names(exp_subset), genelist_new)

# Creating null matrices
mat_coexpression= matrix(data = 0, nrow = length(genelist_new), ncol = length(genelist_new), dimnames = c(list(genelist_new), list(genelist_new)))
mat_coexpression_jaccard = matrix(data = 0, nrow = length(genelist_new), ncol = length(genelist_new), dimnames = c(list(genelist_new), list(genelist_new)))


# Iterate the loop to make  matrix for gene expression 

for(g1 in c(1:length(genelist_new))){
  w1= which(exp_subset[genelist_new[g1],]!=0)
  for(g2 in c(1: length(genelist_new))){
    w2= which(exp_subset[genelist_new[g2],]!=0)
    mat_coexpression[g1,g2]= length(intersect(w1,w2))
    if (length(w1) > 0 & length(w2) > 0) {
      mat_coexpression_jaccard[g1, g2] <- length(intersect(w1, w2)) / length(union(w1, w2))
    }
  }
}

saveRDS(mat_coexpression, "coexpression_matrix_rep_up.RDS")
saveRDS(mat_coexpression_jaccard, "coexpression_jaccardmatrix_rep_up.RDS")



## Plotting the jaccard indexes as heatmap



mat_coexpression_jaccard
diag(mat_coexpression_jaccard) <- NA
heatmap.2(mat_coexpression_jaccard_nulldiagonal,
          trace = "none",  # Remove row and column names
          col = colorRampPalette(c("orange", "red","darkred")),  # Define the color gradient
          main = "Coexpression Jaccard Heatmap",
          key = TRUE,  # Display the color scale legend
)


#### Plotting UMAP overlaying genes 

raw_data <- bcell_combined@assays$RNA@counts
genes <- row.names(raw_data)
setdiff(genes,breggenes_only)
genes_found <- intersect(genes,breggenes_only)
##  "CD1D"   "FASLG"  "IL10"   "CR1"    "SDC1"   "PDCD1"  "CD80"   "CD86"   "MME"    "IL12A"  "NT5E"   "FOXP3"  "IDO1"  
## "CD274"  "ENTPD1" "GZMB"   "THBS1"  "IL32"   "EBI3"   "TGFB1"
exp_subset <- raw_data[genes_found,]
raw_counts <- bcell_combined@assays$RNA@data
m_raw  =t(as.data.frame(raw_counts[genes_found, ]))



#### Plotting the proportion(percent) of cells expressing a gene

add.alpha <- function(col, alpha=1){
  if(missing(col)) # the function checks for missing vectors in col first and if there is missing values it stops
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, # else the function converts each of the vectors into a RGB color with  col2rgb 
        #and then dvides by 255 to make them between 0 1nd 1 then the apply function takes each column in the RGB matrix
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  # function takes rgb values in column and coverts it to a color with transparency 
}
# Use brewer. pal to generate vectors of colors ,example- here 6 is the number of vetcors and Dark2 is the palette , set transparency to 0.95 and apply add alpha function to it to add transperancy to vectors

cols1 =  c(add.alpha (brewer.pal(6, "Dark2"), alpha = 0.95),add.alpha (brewer.pal(8, "Paired"), alpha = 0.95), add.alpha (brewer.pal(8, "Set1"), alpha = 0.95))
# add trabsperancy to vectors but change transparency to 0.5
cols =  add.alpha (cols1, alpha = 0.5)

#extact the rownames of metadata- droplet ids
id <- rownames(bcell_combined@meta.data)
# Create a color palette for each  of the droplets , i.e along the length of the id 
# initialize with grey color 
col_cell <- rep("grey",length(id))
# Assign dropletids as names to each element in the color cell vector 
names(col_cell) <- id

Plot_blood_versus_tissue_levels<-function(){
  fileout1=concat(c(out_dir,"/Seurat_filtering_", batch,"_", analysis,"_1.pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*3, width=w*1.1*5)
  par(mfrow= c(3,5), mar = c(8,5,2,3))
}

pdf("plots_expression_ in_percetage", )
for(g in c(1:length(genes_of_interest))){
  gene_of_interest = genes_of_interest[g]
  w_gene = names(which(count_data[gene_of_interest,]!=0))
  table(cell_type,broad_source)
  list_gene = NULL
  for(i1 in c(1:length(cell_types))){
    w1 = which(cell_type==cell_types[i1])
    list1 = NULL
    for(i2 in c(1:length(broad_sources))){
      w2 = which(broad_source==broad_sources[i2])
      w = intersect(w1,w2)
      disease_sub = disease[w]
      sample_sub = sample[w]
      samples_sub = table(sample_sub)
      samples_sub = names(which(samples_sub>=15))
      disease_sub = NULL
      values = NULL
      for(c in c(1:length(samples_sub))){
        wx = cell_ids[intersect(w, which(sample==samples_sub[c]))]
        n_total = length(wx)
        n_gene = length(intersect(wx, w_gene))
        perc = n_gene*100/n_total
        values = c(values, perc)
        disease_sub = c(disease_sub, disease_sample_map[samples_sub[c]])
      }
      list1 = c(list1, list(c(list(values), list(disease_sub))))
    }
    names(list1) = broad_sources
    list_gene = c(list_gene, list(list1))
    
  }
  names(list_gene) = cell_types
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2"), alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  cbind(diseases,cols1)
  cols_box =  add.alpha (c("blue","black"), alpha = 0.5)
  
  groups = NULL
  col_groups = NULL
  p_values = NULL
  use = NULL
  names = NULL
  for(i in c(1:length(list_gene))){
    m = list_gene[[i]]
    g1 = NULL
    data = NULL
    factor = NULL
    for(g in c(1:length(m))){
      factor = c(factor, list(match(m[[g]][[2]], diseases)))
      g1 = c(g1, list(m[[g]][[1]]))}
    if(min(c(length(g1[[1]]), length(g1[[2]]))) >=3){
      #p = t.test((g1[[1]]), (g1[[2]]), alternative = "less")$p.value
      p = wilcox.test((g1[[1]]), (g1[[2]]), alternative = "less")$p.value # not assuming parametric distribution
      groups = c(groups, list(g1))
      col_groups = c(col_groups, list(factor))
      p_values = c(p_values, p)
      names = c(names, names(list_gene)[i])
    }}
  length(which(p_values<0.05)) ## check
  factors1 = names(list_gene[[1]])
  factors = names
  main = concat(c(gene_of_interest," % of cells "))
  max = max(c(unlist(groups), unlist(groups))*1.35)
  min = 0
  b = (max-min)*0.034
  ylab = "% cells"
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = length(groups)
  max_scale = min(c(max, 100))
  range = max-min
  if(range>55){scale = c(0:100)*20}
  if(range<=55){scale = c(0:100)*10}
  if(range<=40){scale = c(0:100)*5}
  if(range <10){scale = c(0:100)*2.5}
  if(range <5){scale = c(0:100)*1}
  if(range <4){scale = c(0:100)*0.5}
  if(range <1.5){scale = c(0:1000)*0.2}
  if(range <0.5){scale = c(0:100)*0.1}
  if(range <0.1){scale = c(0:100)*0.01}
  if(range <0.01){scale = c(0:100)*0.001}
  cex = 0.9
  Fun<-function(x){x}
  scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
  
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
  mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
  mtext(side = 1, text = gsub("emory","em.", factors), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  width = 0.18
  index = 1
  l = length(groups)
  l1 = length(groups[[1]])
  shift = c(1:l1)
  shift = (mean(shift)-shift)
  shift = shift*0.21/max(shift)
  pches = c(22,21,23,24)
  for(i in c(1:l)){
    for(i1 in c(1:l1)){
      points1=as.numeric(groups[[i]][[i1]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i-shift[i1],width,cols_box[i1],1, cols_box[i1])
      points(rep(i-shift[i1], length(points1)),points1, pch =pches[col_groups[[i]][[i1]]], col=cols_box[i1],bg = cols_box[i1], cex = 0.7)
    }}
  
  for(i in c(1:l)){	
    b = max*0.035
    signif_threshold = 0.05
    if(p_values[i]<signif_threshold){
      pval1 = "*"
      # if(p_values[i] <signif_threshold/10){pval1 = "**"}
      # if(p_values[i] <signif_threshold/100){pval1 = "***"}
      y = max(unlist(groups))
      y = y+2*b
      segments(i-shift,y+b, i+shift,y+b,lwd = 3, col = "darkgrey")
      text(i, y+2*b, labels = pval1, cex = 1.7)
    }
    # Add p-values
    text(i, y + 4 * b, labels = sprintf("p = %.3f", p_values[i]), cex = 1.2)
  }
}

col_leg = "black"
plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE, ylim = c(min, max))
legend("topright", factors1, pch = 23,cex= 0.9, bty="n", pt.bg = "white", col = "white", text.col = cols_box, text.font = 2)
legend("bottomright", diseases, pch = pches,cex= 0.9, bty="n", pt.bg = col_leg, col = col_leg, text.col = col_leg, text.font = 2)
dev.off()

### percent expression in cancer 

bcell_data <- bcell_combined
# Subset the data to include only cancer samples
bcell_data_cancer <- subset(bcell_data, condition == "cancer")
library("Seurat")
library('harmony') 
library(ggplot2)
library(pryr)
library(future) 

library(fgsea)
library(ggplot2)
library(BiocParallel)
library(org.Hs.eg.db)

options(future.globals.maxSize = 30000 * 1024^2)
plan(multisession) 

##Functions


concat = function(v) {
  res = ""
  for (i in 1:length(v)){res = paste0(res,v[i])}
  res
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))
}

Draw_box_plot<-function(box,x,width,c,lwd,line_col){
  segments(x, box[2], x, box[3], col = line_col,lwd =lwd)
  segments(x-(width/2), box[2], x+(width/2), box[2], col = line_col,lwd =lwd)
  segments(x-(width/2), box[3], x+(width/2), box[3], col = line_col,lwd =lwd)
  rect(x-width, box[4], x+width, box[5], col = c,lwd =lwd, border = line_col)
  segments(x-width, box[1], x+width, box[1], col = line_col,lwd=2*lwd)}



file <- "C:/Users/x/Desktop/10xpractice/Cell_type_grouping.csv"
p <- (read.csv(file, head = TRUE))
p <- p[,-1]
colnames(p) <- c("cell_types", "Level.1", "Level.2")
genes_of_interest <-breggenes_only

##### counting the number of cells expressing each gene and determining if they are present more in blood or tissue
count_data= bcell_data_cancer@assays$RNA@data
genes = rownames(count_data)
genes_of_interest[which(genes_of_interest %in% genes == F)]
genes_of_interest = genes_of_interest[which(genes_of_interest %in% genes == T)]
genes_of_interest = sort(unique(genes_of_interest))
# Subset counts of genes of interest
exp_sub = count_data[genes_of_interest,]
w = names(which(rowSums(exp_sub)!=0)) ## exclude all genes that are not picked up (HAVCR1)
exp_sub = exp_sub[w,]
genes_of_interest = w

## subset the metadata by sample type- pBMCor primary tissue 
broad_source = bcell_data_cancer@meta.data$sample.type
#subset the metadata by sample type and dataset  and change all pdac datasets
# to pdac, and all 

disease_type = bcell_data_cancer@meta.data$overall_dataset
disease_type[which(disease_type %in% c("pdac Peng",  "pdac Steele"))] = "pdac"

##subset the metadata by condition  and change all types of cancer ("brain metastasis",  "intestinal metaplasia") to cancer
# to cancer and all inflammatory disease  to  "autoinfammatory"
disease = bcell_data_cancer@meta.data$condition
disease[which(disease %in% c("brain metastasis",  "intestinal metaplasia"))] = "cancer"
disease[which(disease %in% c("ai",  "inflamed tissue"))] = "autoinfammatory"

# extract the perpatient data
sample = bcell_data_cancer@meta.data$orig.ident
## extract metadata of celltypes
cell_type = bcell_data_cancer@meta.data$newid_strict


# Extracting the cell type information
file <- "C:/Users/x/Desktop/10xpractice/Cell_type_grouping.csv"
p <- (read.csv(file, head = TRUE))

cell_type[which(cell_type %in% p[,1]==F)]
cell_type_level1 = p[,"Level.1"]
cell_type_level2 = p[,"Level.2"]
names(cell_type_level1) = p[,1]
names(cell_type_level2) = p[,1]

cells_type_level1 = cell_type_level1[cell_type]
cells_type_level2 = cell_type_level2[cell_type]
cell_type = cells_type_level2

broad_sources = sort(unique(broad_source)) ### only the two sources
diseases = sort(unique(disease)) ## only the disease names
samples = sort(unique(sample)) ## all the patient sample names
cell_types = sort(unique(cell_type)) ## all the cell types names




disease_sample = unique(cbind(sample,disease))  #### makes a dataframe of sample and disease,unique omits repitions from cell
disease_sample_map = disease_sample[,2]
names(disease_sample_map) = disease_sample[,1]


batch="Breg"
Output = "C:/Users/x/Desktop/10xpractice/"

cell_ids <- rownames (bcell_data_cancer@meta.data)


Plot_blood_versus_tissue_levels<-function(){
  fileout1=concat(c(out_dir,"/Seurat_filtering_", batch,"_", analysis,"_1.pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*3, width=w*1.1*5)
  par(mfrow = c(3, 5), mar = c(5, 4, 1, 2))
}
pdf("expression_in_percetage_cancer", )
for(g in c(1:length(genes_of_interest))){
  gene_of_interest = genes_of_interest[g]
  w_gene = names(which(count_data[gene_of_interest,]!=0))
  table(cell_type,broad_source)
  list_gene = NULL
  for(i1 in c(1:length(cell_types))){
    w1 = which(cell_type==cell_types[i1])
    list1 = NULL
    for(i2 in c(1:length(broad_sources))){
      w2 = which(broad_source==broad_sources[i2])
      w = intersect(w1,w2)
      disease_sub = disease[w]
      sample_sub = sample[w]
      samples_sub = table(sample_sub)
      samples_sub = names(which(samples_sub>=15))
      disease_sub = NULL
      values = NULL
      for(c in c(1:length(samples_sub))){
        wx = cell_ids[intersect(w, which(sample==samples_sub[c]))]
        n_total = length(wx)
        n_gene = length(intersect(wx, w_gene))
        perc = n_gene*100/n_total
        values = c(values, perc)
        disease_sub = c(disease_sub, disease_sample_map[samples_sub[c]])
      }
      list1 = c(list1, list(c(list(values), list(disease_sub))))
    }
    names(list1) = broad_sources
    list_gene = c(list_gene, list(list1))
  }
  names(list_gene) = cell_types
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2"), alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  cbind(diseases,cols1)
  cols_box =  add.alpha (c("blue","black"), alpha = 0.5)
  
  groups = NULL
  col_groups = NULL
  p_values = NULL
  use = NULL
  names = NULL
  for(i in c(1:length(list_gene))){
    m = list_gene[[i]]
    g1 = NULL
    data = NULL
    factor = NULL
    for(g in c(1:length(m))){
      factor = c(factor, list(match(m[[g]][[2]], diseases)))
      g1 = c(g1, list(m[[g]][[1]]))}
    if(min(c(length(g1[[1]]), length(g1[[2]]))) >=3){
      #p = t.test((g1[[1]]), (g1[[2]]), alternative = "less")$p.value
      p = wilcox.test((g1[[1]]), (g1[[2]]), alternative = "less")$p.value # not assuming parametric distribution
      groups = c(groups, list(g1))
      col_groups = c(col_groups, list(factor))
      p_values = c(p_values, p)
      names = c(names, names(list_gene)[i])
    }}
  length(which(p_values<0.05)) ## check
  factors1 = names(list_gene[[1]])
  factors = names
  main = concat(c(gene_of_interest," % of cells "))
  max = max(c(unlist(groups), unlist(groups))*1.35)
  min = 0
  b = (max-min)*0.034
  ylab = "% cells"
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = length(groups)
  max_scale = min(c(max, 100))
  range = max-min
  if(range>55){scale = c(0:100)*20}
  if(range<=55){scale = c(0:100)*10}
  if(range<=40){scale = c(0:100)*5}
  if(range <10){scale = c(0:100)*2.5}
  if(range <5){scale = c(0:100)*1}
  if(range <4){scale = c(0:100)*0.5}
  if(range <1.5){scale = c(0:1000)*0.2}
  if(range <0.5){scale = c(0:100)*0.1}
  if(range <0.1){scale = c(0:100)*0.01}
  if(range <0.01){scale = c(0:100)*0.001}
  cex = 0.9
  Fun<-function(x){x}
  scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
  
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
  mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
  mtext(side = 1, text = gsub("emory","em.", factors), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  width = 0.18
  index = 1
  l = length(groups)
  l1 = length(groups[[1]])
  shift = c(1:l1)
  shift = (mean(shift)-shift)
  shift = shift*0.21/max(shift)
  pches = c(22,21,23,24)
  for(i in c(1:l)){
    for(i1 in c(1:l1)){
      points1=as.numeric(groups[[i]][[i1]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i-shift[i1],width,cols_box[i1],1, cols_box[i1])
      points(rep(i-shift[i1], length(points1)),points1, pch =pches[col_groups[[i]][[i1]]], col=cols_box[i1],bg = cols_box[i1], cex = 0.7)
    }}
  
  for(i in c(1:l)){	
    b = max*0.035
    signif_threshold = 0.05
    if(p_values[i]<signif_threshold){
      pval1 = "*"
      # if(p_values[i] <signif_threshold/10){pval1 = "**"}
      # if(p_values[i] <signif_threshold/100){pval1 = "***"}
      y = max(unlist(groups))
      y = y+2*b
      segments(i-shift,y+b, i+shift,y+b,lwd = 3, col = "darkgrey")
      text(i, y+2*b, labels = pval1, cex = 1.7)
    }
  }
  for (i in c(1:l)) {
    b = max * 0.035
    signif_threshold = 0.05
    if (p_values[i] < signif_threshold) {
      pval1 = "*"
      y = max(unlist(groups))
      y = y + 2 * b
      segments(i - shift, y + b, i + shift, y + b, lwd = 3, col = "darkgrey")
      text(i, y + 2 * b, labels = pval1, cex = 1.7)
      # Add p-values
      text(i, y + 4 * b, labels = sprintf("p = %.3f", p_values[i]), cex = 1.2)
    }
  }
  
}
col_leg = "black"
plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE, ylim = c(min, max))
legend("topright", factors1, pch = 23,cex= 0.9, bty="n", pt.bg = "white", col = "white", text.col = cols_box, text.font = 2)
legend("bottomright", diseases, pch = pches,cex= 0.9, bty="n", pt.bg = col_leg, col = col_leg, text.col = col_leg, text.font = 2)

dev.off()


