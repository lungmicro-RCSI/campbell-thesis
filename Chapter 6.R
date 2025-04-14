#Chapter 6 code 

#libraries
library(ggplot2)
library(ggrepel)
library(data.table)
library(dplyr)
library(wesanderson)
library(gtsummary)
library(readxl)
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(vegan)
library(DESeq2)
library(data.table)

###############################################################################
#Nasal Rinse Analysis 
###############################################################################
#set wd 
setwd("/Users/christinacampbell/Desktop/Thesis /Data/merged")

#Taxa table
taxa <- read.csv("/Users/christinacampbell/Desktop/Thesis /Data/merged/tax_covid_sleep.csv", row.names=1)

#OTU Table 
otu<-read.csv("/Users/christinacampbell/Desktop/Thesis /Data/merged/otu_covid_sleep.csv",
              row.names = 1, header = T, check.names=F)

#Load MetaData
meta<-read.delim2("/Users/christinacampbell/Desktop/Thesis /Data/merged/metadata_git.txt", sep="\t") 

#Give Rownames
rownames(meta)<-meta$SampleID

#create Mapping table 
mapping.table=sample_data(meta)

#Create OTU Table 
otu = otu_table(otu, taxa_are_rows=TRUE)

#Create Taxa table 
taxa <-as.matrix(taxa)
TAX <- tax_table(taxa)

#make phyloseq
physeq=phyloseq(otu, TAX, meta)

#only pts who've had a visit 2 
nasal <- subset_samples(nasal, !is.na(sample_data(nasal)$redcap_event_name.y))

#only those prescribed CPAP 
nasal<-subset_samples(nasal, ref_cpap =="1")

#keep only the nasal rinse samples for now
nasal<-subset_samples(nasal, Description == "Nasal.Rinse")

#make relative abundance table 
nasal.rel <- transform_sample_counts(nasal, function(x) x/sum(x))

theme<-theme(panel.background = element_blank(),
             panel.border=element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.line = element_line(colour = "black"),
             axis.text.x=element_text(colour="black"),
             text = element_text(size = 40, colour="black"),
             axis.ticks = element_line(colour = "black", size = 2),  
             axis.ticks.length = unit(0.5, "cm"))

###############################################################################
#Figure 6.2A
###############################################################################
#Calculates Shannon Diversity
sample_data(nasal.rel)$Shannon = diversity(otu_table(nasal.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(nasal.rel))

#Make variable of interest a factor
Shannon$visit <- as.factor(Shannon$visit)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#plot
pdf("Alpha Diversity - V1 v V2 - points connected - CPAP Only - Nasal1.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=visit, y=Shannon, fill=visit)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_point(aes(shape=adherent, color = adherent), size = 5)+
  geom_line(aes(group=subj_id), linetype = 2) +
  scale_fill_manual(values=c( "1" ="#e7c6ff",
                              "2" = "#b9c0ff")) + 
  scale_shape_manual(values = c(16, 17), labels = c("0" = "Non-Adherent", "1" = "Adherent"), name = NULL) +
  scale_color_manual(values = c("1" = "#113ff7", "0" = "black"), labels = c("0" = "Non-Adherent", "1" = "Adherent"), name = NULL) + 
  ylab("Shannon Diversity") + 
  xlab("Visit") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = "none") 
dev.off()

#check stats - Wilcoxan as paired 
# Subset the data by visit
Shannon <- Shannon[rownames(Shannon) != "C.Sleep032.V2.NR", ]
x <- Shannon$Shannon[Shannon$visit == "1"]
y <- Shannon$Shannon[Shannon$visit == "2"]
result <- wilcox.test(x, y, paired = TRUE)

###############################################################################
#Figure 6.2B
###############################################################################

##Create Distance Matrix
vegdist   = vegdist(t(otu_table(nasal.rel)), method="bray")

##Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

##calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)

##Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

##Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(nasal.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~visit,data= newResults, mean)

centroids$visit<-as.factor(centroids$visit)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="visit",suffixes=c("",".centroid"))

newResults$visit<-as.character(newResults$visit)

#plot it 
pdf("Beta Diversity - v1 v v2 - NR - CPAP only.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=visit)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("1" ="#e7c6ff",
                               "2" = "#b9c0ff"))+
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=visit), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=visit))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=visit), parse=TRUE,size=10, show.legend=FALSE) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 40, colour = "black"),
        axis.title.y = element_text(size = 40, colour = "black"),
        axis.text.x = element_text(size = 40, colour = "black"), 
        axis.text.y = element_text(size = 40, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"))
dev.off()

#Create Table for Statistics    
data.adonis <- data.frame(sample_data(nasal.rel))

#Run the Statistics
samplepermanova <- adonis2(vegdist ~ visit, data = data.adonis, strata = data.adonis$SubjID)
print(samplepermanova)

###############################################################################
#Figure 6.3
###############################################################################

#Set alpha for differential Anlaysis, for adjusted p value 
alpha <- 0.2

#Convert To DESEQ, different type of phyloseq object for this program
diagdds <- phyloseq_to_deseq2(nasal, ~SubjID + visit)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds$visit<-as.factor(diagdds$visit)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$visit <- droplevels(diagdds$visit)

#Relevel Data
diagdds$visit <- relevel(diagdds$visit, ref ="1")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]


#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(nasal)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)


#Create name with family and (u.g), creating a new column gs with names thats the highest known class
res$gs <- ifelse(is.na(res$Species),paste0(res$Genus,"_",res$ASV), paste0(res$Genus,"_",res$Species,"_",res$ASV))
res$gs <- ifelse(is.na(res$Genus),paste0(res$Family,"_",res$ASV), res$gs)
res$gs <- ifelse(is.na(res$Family), paste0(res$Order,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Order), paste0(res$Class, "_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Class), paste0(res$Phylum,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Phylum), paste0(res$Kingdom,"_",res$ASV),res$gs)

#Make the full trail the First Column, make this new column the 1st column. 
res$names <- res$ASV
res$Gene.symbol <- res$gs

#Subset the different categories, subset by smokers / COPD 
df.1 = subset_samples(nasal.rel, visit %in% c("1"))
df.2 = subset_samples(nasal.rel, visit  %in% c("2"))

#decide what otu to save 
otu.to.save <-as.character(res$names)

#convert to dataframe
df.1.df <- data.frame(otu_table(df.1))
df.2.df <- data.frame(otu_table(df.2))

#from relative table we should get the mean across the row of the otu table, mean of the relative abundance. 
df.1.meanRA <- rowMeans(df.1.df)
df.2.meanRA <- rowMeans(df.2.df)

#need to subset AND reorder just the otus that we have 
df.1.meanRA.save <- df.1.meanRA[otu.to.save]
df.2.meanRA.save <- df.2.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance.1 <- df.1.meanRA.save
res$abundance.2 <- df.2.meanRA.save

#Keep only the count data
drops <- c("Domain","Phylum","Class","Order","Family","Genus","OTU","gs","Species",
           "ASV","Kingdom","names","row2")
res <- res[ , !(names(res) %in% drops)]

#Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL),
                c("ASV","baseMean", "logFC", "lfcSE", "stat", "pvalue", 
                  "adj.P.Val","Gene.symbol","abundance.1","abundance.2")) 

#Write Tables to TXT file
#write.table(res,file=  "Stratum - Differential - smokers vs COPD.txt", sep="\t",
#            col.names = NA, row.names = TRUE)

#Convert to data.frame
res <- as.data.frame(res)

#make an abundance variable for the size of the dots on the plot
res$abundance.2 <- as.numeric(as.character(res$abundance.2))
res$abundance.1 <- as.numeric(as.character(res$abundance.1))
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha, res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment direction to colour the graph
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Visit 2", "Visit 1")

pdf("Differential - v1 v v2 - paired analysis - CPAP only.pdf",width=35, height=12)
ggplot(sig_results, aes(x=Gene.symbol, y=logFC)) +
  geom_point(aes(fill = enriched, size = abundance),   shape = 21, color = "black")+
  scale_alpha_identity() +
  geom_segment( aes(xend=Gene.symbol, yend=0)) +
  #geom_col(width=0.005, color="black") +
  #scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  #ylim(12,15) +
  scale_size_continuous(name="Relative\nAbundance", range = c(3, 20), ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  scale_fill_manual(values = c("Visit 1"="#AA7DCE","Visit 2"="#F4A5AE"), guide = "none") +
  scale_color_manual(values = c("Visit 1"="#AA7DCE","Visit 2"="#F4A5AE"), guide = "none") +
  xlab("") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 50, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

###############################################################################
#Figure 6.5A
###############################################################################

#subset for adherence 
adherent = subset_samples(nasal, adherent == "1")

#make rel abundance
adherent.rel <- transform_sample_counts(adherent, function(x) x/sum(x))

#Calculates Shannon Diversity
sample_data(adherent.rel)$Shannon = diversity(otu_table(adherent.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(adherent.rel))

#Make variable of interest a factor
Shannon$adherent <- as.factor(Shannon$adherent)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#Make variable of interest a factor
Shannon$visit <- as.factor(Shannon$visit)


pdf("Alpha Diversity - V1 & V2 - adherent group - nasal.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=visit, y=Shannon, fill=visit)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_point(size=2.5)+
  geom_line(aes(group=subj_id), linetype = 2) +
  scale_fill_manual(values=c( "1" ="#e7c6ff",
                              "2" = "#b9c0ff")) + 
  ylab("Shannon Diversity") + 
  xlab("Visit") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme+
  guides(fill = FALSE) 
dev.off()


#check stats - Wilcoxan as paired 
# Subset the data by visit
Shannon <- Shannon[rownames(Shannon) != "C.Sleep032.V2.NR", ]
x <- Shannon$Shannon[Shannon$visit == "1"]
y <- Shannon$Shannon[Shannon$visit == "2"]
result <- wilcox.test(x, y, paired = TRUE)



###############################################################################
#Figure 6.5B
###############################################################################
##Create Distance Matrix
vegdist   = vegdist(t(otu_table(adherent.rel)), method="bray")

##Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

##calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)

##Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

##Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(adherent.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~visit,data= newResults, mean)

centroids$visit<-as.factor(centroids$visit)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="visit",suffixes=c("",".centroid"))

newResults$visit<-as.character(newResults$visit)


#plot it 
pdf("Beta Diversity - v1 v v2 - adherent group - nasal.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=visit)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c( "1" ="#e7c6ff",
                                "2" = "#b9c0ff"))+
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=visit), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=visit))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=visit), parse=TRUE,size=10, show.legend=FALSE) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 40, colour = "black"),
        axis.title.y = element_text(size = 40, colour = "black"),
        axis.text.x = element_text(size = 40, colour = "black"), 
        axis.text.y = element_text(size = 40, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"))
dev.off()

#Create Table for Statistics    
data.adonis <- data.frame(sample_data(adherent.rel))

#Run the Statistics
samplepermanova <- adonis2(vegdist ~ visit, data = data.adonis, strata = data.adonis$SubjID)
print(samplepermanova)

###############################################################################
##Figure 6.6
###############################################################################

#Set alpha for differential Anlaysis, for adjusted p value 
alpha <- 0.2


#Convert To DESEQ, different type of phyloseq object for this program
diagdds <- phyloseq_to_deseq2(adherent, ~SubjID + visit)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
#diagdds = estimateSizeFactors(diagdds)
diagdds$visit<-as.factor(diagdds$visit)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$visit <- droplevels(diagdds$visit)

#Relevel Data
diagdds$visit <- relevel(diagdds$visit, ref ="1")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(nasal)[rownames(res), ], "matrix"))
#res.bal = cbind(as(res.bal, "data.frame"), as(tax_table(Species.bal.table)[rownames(res.bal), ], "matrix"))


#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)
print(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Create name with family and (u.g), creating a new column gs with names thats the highest known class
res$gs <- ifelse(is.na(res$Species),paste0(res$Genus,"_",res$ASV), paste0(res$Genus,"_",res$Species,"_",res$ASV))
res$gs <- ifelse(is.na(res$Genus),paste0(res$Family,"_",res$ASV), res$gs)
res$gs <- ifelse(is.na(res$Family), paste0(res$Order,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Order), paste0(res$Class, "_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Class), paste0(res$Phylum,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Phylum), paste0(res$Kingdom,"_",res$ASV),res$gs)

#Make the full trail the First Column, make this new column the 1st column. 
res$names <- res$ASV
res$Gene.symbol <- res$gs

#Subset the different categories, subset by smokers / COPD 
df.1 = subset_samples(nasal.rel, visit %in% c("1"))
df.2 = subset_samples(nasal.rel, visit  %in% c("2"))

#decide what otu to save 
otu.to.save <-as.character(res$names)

#convert to dataframe
df.1.df <- data.frame(otu_table(df.1))
df.2.df <- data.frame(otu_table(df.2))

#from relative table we should get the mean across the row of the otu table, mean of the relative abundance. 
df.1.meanRA <- rowMeans(df.1.df)
df.2.meanRA <- rowMeans(df.2.df)

#need to subset AND reorder just the otus that we have 
df.1.meanRA.save <- df.1.meanRA[otu.to.save]
df.2.meanRA.save <- df.2.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance.1 <- df.1.meanRA.save
res$abundance.2 <- df.2.meanRA.save

#Keep only the count data
drops <- c("Domain","Phylum","Class","Order","Family","Genus","OTU","gs","Species",
           "ASV","Kingdom","names","row2")
res <- res[ , !(names(res) %in% drops)]

#Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL),
                c("ASV","baseMean", "logFC", "lfcSE", "stat", "pvalue", 
                  "adj.P.Val","Gene.symbol","abundance.1","abundance.2")) 

#Convert to data.frame
res <- as.data.frame(res)

#make an abundance variable for the size of the dots on the plot
res$abundance.2 <- as.numeric(as.character(res$abundance.2))
res$abundance.1 <- as.numeric(as.character(res$abundance.1))
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha, res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment direction to colour the graph
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Visit 2", "Visit 1")

#plot it 
pdf("Differential - v1 v v2 - paired analysis - adherent group.pdf",width=35, height=7)
ggplot(sig_results, aes(x=Gene.symbol, y=logFC)) +
  geom_point(aes(fill = enriched, size = abundance),   shape = 21, color = "black")+
  scale_alpha_identity() +
  geom_segment( aes(xend=Gene.symbol, yend=0)) +
  #geom_col(width=0.005, color="black") +
  #scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  #ylim(12,15) +
  scale_size_continuous(name="         Relative\n         Abundance", range = c(3, 20), ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  scale_fill_manual(values = c("Visit 1"="#e7c6ff","Visit 2"="#b9c0ff"), guide = "none") +
  scale_color_manual(values = c("Visit 1"="#e7c6ff","Visit 2"="#b9c0ff"), guide = "none") +
  xlab("") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 50, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

###############################################################################
##Figure 6.8A
###############################################################################
#subset for visit 2
nasal1<-subset_samples(nasal, Visit == "2")
nasal.rel1<-subset_samples(nasal.rel, Visit == "2")

#Calculates Shannon Diversity
sample_data(nasal.rel1)$Shannon = diversity(otu_table(nasal.rel1), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(nasal.rel1))

#Make variable of interest a factor
Shannon$adherent <- as.factor(Shannon$adherent)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#set factors
Shannon$adherent<-factor(Shannon$adherent, c("Adherent", "Non-Adherent"))

#plot it
pdf("Alpha Diversity - adherence - Nasal.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=adherent, y=Shannon, fill=adherent)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Adherent"="#BAD7F2","Non-Adherent"="#96D09A")) + 
  ylab("Shannon Diversity") + 
  xlab("Adherent") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme+
  guides(fill = FALSE) 
dev.off()

#check stats
sampleshannon <- kruskal.test(Shannon ~ adherent, data = Shannon)
sampleshannon <- sampleshannon$p.value
print(sampleshannon)

###############################################################################
##Figure 6.8B
###############################################################################


##Create Distance Matrix
vegdist   = vegdist(t(otu_table(nasal.rel1)), method="bray")

##Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

##calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)

##Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))
#
##Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(nasal.rel1), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
#

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~adherent,data= newResults, mean)

centroids$adherent<-as.factor(centroids$adherent)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="adherent",suffixes=c("",".centroid"))

newResults$adherent<-as.character(newResults$adherent)


pdf("Beta Diversity - adherent - nasal.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=adherent)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Adherent"="#BAD7F2","Non-Adherent"="#96D09A"), guide="none")+
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=adherent), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=adherent))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=adherent), parse=TRUE,size=10, show.legend=FALSE) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 30))
dev.off()

#Create Table for Statistics    
data.adonis <- data.frame(sample_data(nasal.rel))

#Run the Statistics
samplepermanova <- adonis(vegdist ~ adherent, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]



###############################################################################
##Figure 6.9
###############################################################################

#Set alpha for differential Analysis, for adjusted p value 
alpha <- 0.2


#Convert To DESEQ, different type of phyloseq object for this program
diagdds <- phyloseq_to_deseq2(nasal1, ~adherent)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
#diagdds = estimateSizeFactors(diagdds)
diagdds$adherent<-as.factor(diagdds$adherent)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$adherent<- droplevels(diagdds$adherent)

#Relevel Data
diagdds$adherent <- relevel(diagdds$adherent, ref ="Non-Adherent")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(nasal)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Create name with family and (u.g), creating a new column gs with names thats the highest known class
res$gs <- ifelse(is.na(res$Genus),paste0(res$Family,"_",res$ASV), paste0(res$Genus,"_",res$ASV))
res$gs <- ifelse(is.na(res$Family), paste0(res$Order,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Order), paste0(res$Class, "_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Class), paste0(res$Phylum,"_",res$ASV),res$gs)

#Make the full trail the First Column, make this new column the 1st column. 
res$names <- res$ASV
res$Gene.symbol <- res$gs

#Subset the different categories, subset by smokers / COPD 
df.1 = subset_samples(nasal.rel, adherent %in% c("Non-Adherent"))
df.2 = subset_samples(nasal.rel, adherent  %in% c("Adherent"))

#decide what otu to save 
otu.to.save <-as.character(res$names)

#convert to dataframe
df.1.df <- data.frame(otu_table(df.1))
df.2.df <- data.frame(otu_table(df.2))

#from relative table we should get the mean across the row of the otu table, mean of the relative abundance. 
df.1.meanRA <- rowMeans(df.1.df)
df.2.meanRA <- rowMeans(df.2.df)

#need to subset AND reorder just the otus that we have 
df.1.meanRA.save <- df.1.meanRA[otu.to.save]
df.2.meanRA.save <- df.2.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance.1 <- df.1.meanRA.save
res$abundance.2 <- df.2.meanRA.save

#Keep only the count data
drops <- c("Domain","Phylum","Class","Order","Family","Genus","OTU","gs","Species",
           "ASV","Kingdom","names","row2")
res <- res[ , !(names(res) %in% drops)]

#Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL),
                c("ASV","baseMean", "logFC", "lfcSE", "stat", "pvalue", 
                  "adj.P.Val","Gene.symbol","abundance.1","abundance.2")) 

#Convert to data.frame
res <- as.data.frame(res)

#make an abundance variable for the size of the dots on the plot
res$abundance.2 <- as.numeric(as.character(res$abundance.2))
res$abundance.1 <- as.numeric(as.character(res$abundance.1))
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha, res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment direction to colour the graph
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Adherent", "Non-Adherent")

#plot it 
pdf("Differential - adherence.pdf",width=35, height=12)
ggplot(sig_results, aes(x=Gene.symbol, y=logFC)) +
  geom_point(aes(fill=enriched, size=abundance, color = enriched), shape = 21) +
  geom_segment( aes(xend=Gene.symbol, yend=0)) +
  #geom_col(width=0.005, color="black") +
  #scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  #ylim(12,15) +
  scale_size_continuous(name="Relative\nAbundance", range = c(3, 20), ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  scale_fill_manual(values = c("Non-Adherent"="#96D09A","Adherent"="#BAD7F2"), guide = "none") +
  scale_color_manual(values = c("Non-Adherent"="#96D09A","Adherent"="#BAD7F2"), guide = "none") +
  xlab("") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 50, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

###############################################################################
##Figure 6.17A
###############################################################################


#Set alpha for differential Anlaysis, for adjusted p value 
alpha <- 0.2

#Subset only variables for comparisons
nasal1 = subset_samples(nasal, patient_type1  %in% c("Post_Covid"))

#Subset only variables for comparisons
nasal1 = subset_samples(nasal1, change.ess  %in% c("No","Yes"))

#Convert To DESEQ, different type of phyloseq object for this program
diagdds <- phyloseq_to_deseq2(nasal1, ~change.ess)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds$change.ess<-as.factor(diagdds$change.ess)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$change.ess <- droplevels(diagdds$change.ess)

#Relevel Data
diagdds$change.ess <- relevel(diagdds$change.ess, ref ="No")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(nasal1)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)


#Create name with family and (u.g), creating a new column gs with names thats the highest known class
res$gs <- ifelse(is.na(res$Species),paste0(res$Genus,"_",res$ASV), paste0(res$Genus,"_",res$Species,"_",res$ASV))
res$gs <- ifelse(is.na(res$Genus),paste0(res$Family,"_",res$ASV), res$gs)
res$gs <- ifelse(is.na(res$Family), paste0(res$Order,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Order), paste0(res$Class, "_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Class), paste0(res$Phylum,"_",res$ASV),res$gs)

#Make the full trail the First Column, make this new column the 1st column. 
res$names <- res$ASV
res$Gene.symbol <- res$gs

#Subset the different categories, subset by smokers / COPD 
df.1 = subset_samples(nasal.rel, change.ess %in% c("No"))
df.2 = subset_samples(nasal.rel, change.ess  %in% c("Yes"))

#decide what otu to save 
otu.to.save <-as.character(res$names)

#convert to dataframe
df.1.df <- data.frame(otu_table(df.1))
df.2.df <- data.frame(otu_table(df.2))

#from relative table we should get the mean across the row of the otu table, mean of the relative abundance. 
df.1.meanRA <- rowMeans(df.1.df)
df.2.meanRA <- rowMeans(df.2.df)

#need to subset AND reorder just the otus that we have 
df.1.meanRA.save <- df.1.meanRA[otu.to.save]
df.2.meanRA.save <- df.2.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance.1 <- df.1.meanRA.save
res$abundance.2 <- df.2.meanRA.save

#Keep only the count data
drops <- c("Domain","Phylum","Class","Order","Family","Genus","OTU","gs","Species",
           "ASV","Kingdom","names","row2")
res <- res[ , !(names(res) %in% drops)]

#Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL),
                c("ASV","baseMean", "logFC", "lfcSE", "stat", "pvalue", 
                  "adj.P.Val","Gene.symbol","abundance.1","abundance.2")) 

#Convert to data.frame
res <- as.data.frame(res)

#make an abundance variable for the size of the dots on the plot
res$abundance.2 <- as.numeric(as.character(res$abundance.2))
res$abundance.1 <- as.numeric(as.character(res$abundance.1))
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha, res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment direction to colour the graph
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Sleepiness Resolved", "Persistantly Sleepy")

#plot it
pdf("Differential - change in ESS - adherent.pdf",width=35, height=10)
ggplot(sig_results, aes(x=Gene.symbol, y=logFC)) +
  geom_point(aes(fill = enriched, size = abundance),   shape = 21, color = "black")+
  scale_alpha_identity() +
  geom_segment( aes(xend=Gene.symbol, yend=0)) +
  #geom_col(width=0.005, color="black") +
  #scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  #ylim(12,15) +
  scale_size_continuous(name="         Relative\n         Abundance", range = c(3, 20), ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  scale_fill_manual(values = c("Sleepiness Resolved"="#EFD18D","Persistantly Sleepy"="#D09D9F"), guide = "none") +
  scale_color_manual(values = c("Sleepiness Resolved"="#EFD18D","Persistantly Sleepy"="#D09D9F"), guide = "none") +
  xlab("") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 50, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()


###############################################################################
##Figure 6.17B
###############################################################################
#Set alpha for differential Anlaysis, for adjusted p value 
alpha <- 0.2

#Subset only variables for comparisons
nasal1 = subset_samples(nasal, patient_type1  %in% c("Post_Covid"))

#Subset only variables for comparisons
nasal1 = subset_samples(nasal, change.cfs  %in% c("Yes","No"))
nasal1.rel <- transform_sample_counts(nasal1, function(x) x/sum(x))

#Convert To DESEQ, different type of phyloseq object for this program
diagdds <- phyloseq_to_deseq2(nasal1, ~change.cfs)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds$change.cfs<-as.factor(diagdds$change.cfs)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$change.cfs <- droplevels(diagdds$change.cfs)

#Relevel Data
diagdds$change.cfs <- relevel(diagdds$change.cfs, ref ="No")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(nasal1)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)


#Create name with family and (u.g), creating a new column gs with names thats the highest known class
res$gs <- ifelse(is.na(res$Species),paste0(res$Genus,"_",res$ASV), paste0(res$Genus,"_",res$Species,"_",res$ASV))
res$gs <- ifelse(is.na(res$Genus),paste0(res$Family,"_",res$ASV), res$gs)
res$gs <- ifelse(is.na(res$Family), paste0(res$Order,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Order), paste0(res$Class, "_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Class), paste0(res$Phylum,"_",res$ASV),res$gs)

#Make the full trail the First Column, make this new column the 1st column. 
res$names <- res$ASV
res$Gene.symbol <- res$gs

#Subset the different categories, subset by smokers / COPD 
df.1 = subset_samples(nasal1.rel, change.cfs %in% c("No"))
df.2 = subset_samples(nasal1.rel, change.cfs  %in% c("Yes"))

#decide what otu to save 
otu.to.save <-as.character(res$names)

#convert to dataframe
df.1.df <- data.frame(otu_table(df.1))
df.2.df <- data.frame(otu_table(df.2))

#from relative table we should get the mean across the row of the otu table, mean of the relative abundance. 
df.1.meanRA <- rowMeans(df.1.df)
df.2.meanRA <- rowMeans(df.2.df)

#need to subset AND reorder just the otus that we have 
df.1.meanRA.save <- df.1.meanRA[otu.to.save]
df.2.meanRA.save <- df.2.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance.1 <- df.1.meanRA.save
res$abundance.2 <- df.2.meanRA.save

#Keep only the count data
drops <- c("Domain","Phylum","Class","Order","Family","Genus","OTU","gs","Species",
           "ASV","Kingdom","names","row2")
res <- res[ , !(names(res) %in% drops)]

#Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL),
                c("ASV","baseMean", "logFC", "lfcSE", "stat", "pvalue", 
                  "adj.P.Val","Gene.symbol","abundance.1","abundance.2")) 

#Convert to data.frame
res <- as.data.frame(res)

#make an abundance variable for the size of the dots on the plot
res$abundance.2 <- as.numeric(as.character(res$abundance.2))
res$abundance.1 <- as.numeric(as.character(res$abundance.1))
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha, res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment direction to colour the graph
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Fatigue Resolved", "Persistantly Fatigued")

pdf("Differential - change in CFS - adherent group -just Post-COVID - updated - 1.pdf",width=35, height=10)
ggplot(sig_results, aes(x=Gene.symbol, y=logFC)) +
  geom_point(aes(fill = enriched, size = abundance),   shape = 21, color = "black")+
  scale_alpha_identity() +
  geom_segment( aes(xend=Gene.symbol, yend=0)) +
  #geom_col(width=0.005, color="black") +
  #scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  #ylim(12,15) +
  scale_size_continuous(name="Relative\nAbundance", range = c(3, 20), ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  scale_fill_manual(values = c("Fatigue Resolved"="#c9e4d3","Persistantly Fatigued"="#faedcb"), guide = "none") +
  scale_color_manual(values = c("Fatigue Resolved"="#c9e4d3","Persistantly Fatigued"="#faedcb"), guide = "none") +
  xlab("") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 50, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

###############################################################################
##Oral Wash Analysis 
###############################################################################
#set wd 
setwd("/Users/christinacampbell/Desktop/Thesis /Data/merged")

#Taxa table
taxa <- read.csv("/Users/christinacampbell/Desktop/Thesis /Data/merged/tax_covid_sleep.csv", row.names=1)


#OTU Table 
otu<-read.csv("/Users/christinacampbell/Desktop/Thesis /Data/merged/otu_covid_sleep.csv",
              row.names = 1, header = T, check.names=F)


#Load MetaData
meta<-read.delim2("/Users/christinacampbell/Desktop/Thesis /Data/merged/metadata_git.txt", sep="\t") 

#Give Rownames
rownames(meta)<-meta$SampleID

#create Mapping table 
mapping.table=sample_data(meta)

#Create OTU Table 
otu = otu_table(otu, taxa_are_rows=TRUE)

#Create Taxa table 
taxa <-as.matrix(taxa)
TAX <- tax_table(taxa)

#make phyloseq
physeq=phyloseq(otu, TAX, meta)
#only pts who've had a visit 2 
oral <- subset_samples(oral, !is.na(sample_data(oral)$redcap_event_name.y))

#only those prescribed CPAP 
oral<-subset_samples(oral, ref_cpap =="1")

#keep only the oral wash samples for now
oral<-subset_samples(oral, Description == "Mouth.Wash")

#make relative abundance table 
oral.rel <- transform_sample_counts(oral, function(x) x/sum(x))

theme<-theme(panel.background = element_blank(),
             panel.border=element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.line = element_line(colour = "black"),
             axis.text.x=element_text(colour="black"),
             text = element_text(size = 40, colour="black"),
             axis.ticks = element_line(colour = "black", size = 2),  
             axis.ticks.length = unit(0.5, "cm"))


###############################################################################
##Figure 6.4A
###############################################################################

#Calculates Shannon Diversity
sample_data(oral.rel)$Shannon = diversity(otu_table(oral.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(oral.rel))

#Make variable of interest a factor
Shannon$visit <- as.factor(Shannon$visit)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#plot
pdf("Alpha Diversity - V1 v V2 - points connected - CPAP Only - MW.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=visit, y=Shannon, fill=visit)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_point(aes(shape=adherent, color = adherent), size = 5)+
  geom_line(aes(group=subj_id), linetype = 2) +
  scale_fill_manual(values=c("1" ="#e7c6ff",
                             "2" = "#b9c0ff")) + 
  scale_shape_manual(values = c(16, 17), labels = c("0" = "Non-Adherent", "1" = "Adherent"), name = NULL) +
  scale_color_manual(values = c("1" = "#113ff7", "0" = "black"), labels = c("0" = "Non-Adherent", "1" = "Adherent"), name = NULL) + 
  ylab("Shannon Diversity") + 
  xlab("Visit") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = "none") 
dev.off()


#check stats - Wilcoxan as paired 
# Subset the data by visit
x <- Shannon$Shannon[Shannon$visit == "1"]
y <- Shannon$Shannon[Shannon$visit == "2"]
result <- wilcox.test(x, y, paired = TRUE)



###############################################################################
##Figure 6.4B
###############################################################################

##Create Distance Matrix
vegdist   = vegdist(t(otu_table(oral.rel)), method="bray")

##Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

##calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)

##Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))
#
##Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
#

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~visit,data= newResults, mean)

centroids$visit<-as.factor(centroids$visit)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="visit",suffixes=c("",".centroid"))

newResults$visit<-as.character(newResults$visit)

#set colours

pdf("Beta Diversity - v1 v v2 - MW - CPAP only.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=visit)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("1" ="#e7c6ff",
                               "2" = "#b9c0ff"))+
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=visit), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=visit))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=visit), parse=TRUE,size=10, show.legend=FALSE) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 40, colour = "black"),
        axis.title.y = element_text(size = 40, colour = "black"),
        axis.text.x = element_text(size = 40, colour = "black"), 
        axis.text.y = element_text(size = 40, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"))
dev.off()



#Create Table for Statistics    
data.adonis <- data.frame(sample_data(oral.rel))

#Run the Statistics
samplepermanova <- adonis2(vegdist ~ visit, data = data.adonis, strata = data.adonis$SubjID)
print(samplepermanova)
#p=0.418

##############################################################################
#Figure 6.7A
###############################################################################
#subset for adherence 
adherent = subset_samples(oral, adherent == "1")

#make rel abundance
adherent.rel <- transform_sample_counts(adherent, function(x) x/sum(x))

#Calculates Shannon Diversity
sample_data(adherent.rel)$Shannon = diversity(otu_table(adherent.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(adherent.rel))

#Make variable of interest a factor
Shannon$adherent <- as.factor(Shannon$adherent)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#Make variable of interest a factor
Shannon$visit <- as.factor(Shannon$visit)


pdf("Alpha Diversity - V1 & V2 - adherent group - MW.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=visit, y=Shannon, fill=visit)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_point(size=3)+
  geom_line(aes(group=subj_id), linetype = 2) +
  scale_fill_manual(values=c("1" ="#e7c6ff",
                             "2" = "#b9c0ff")) + 
  ylab("Shannon Diversity") + 
  xlab("Visit") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

getwd()
#check stats - Wilcoxan as paired 
# Subset the data by visit

x <- Shannon$Shannon[Shannon$visit == "1"]
y <- Shannon$Shannon[Shannon$visit == "2"]
result <- wilcox.test(x, y, paired = TRUE)

###############################################################################
##Figure 6.7B
##############################################################################

##Create Distance Matrix
vegdist   = vegdist(t(otu_table(adherent.rel)), method="bray")

##Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

##calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)

##Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))
#
##Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(adherent.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
#

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~visit,data= newResults, mean)

centroids$visit<-as.factor(centroids$visit)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="visit",suffixes=c("",".centroid"))

newResults$visit<-as.character(newResults$visit)

#plot it 
pdf("Beta Diversity - v1 v v2 - adherent group - MW.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=visit)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("1" ="#e7c6ff",
                               "2" = "#b9c0ff"), guide = "none")+
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=visit), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=visit))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=visit), parse=TRUE,size=10, show.legend=FALSE) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 40, colour = "black"),
        axis.title.y = element_text(size = 40, colour = "black"),
        axis.text.x = element_text(size = 40, colour = "black"), 
        axis.text.y = element_text(size = 40, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"))
dev.off()

#Create Table for Statistics    
data.adonis <- data.frame(sample_data(adherent.rel))

#Run the Statistics
samplepermanova <- adonis2(vegdist ~ visit, data = data.adonis, strata = data.adonis$SubjID)
print(samplepermanova)

###############################################################################
##Figure 6.10A
###############################################################################
#subset for visit 2
oral1<-subset_samples(oral, Visit == "2")
oral.rel1<-subset_samples(oral.rel, Visit == "2")

#Calculates Shannon Diversity
sample_data(oral.rel1)$Shannon = diversity(otu_table(oral.rel1), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(oral.rel1))

#Make variable of interest a factor
Shannon$adherent <- as.factor(Shannon$adherent)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#plot it
pdf("Alpha Diversity - adherence - MW.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=adherent, y=Shannon, fill=adherent)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Adherent"="#BAD7F2","Non-Adherent"="#96D09A")) + 
  xlab("") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

#check stats
sampleshannon <- kruskal.test(Shannon ~ adherent, data = Shannon)
sampleshannon <- sampleshannon$p.value
print(sampleshannon)

###############################################################################
##Figure 6.10B
###############################################################################
##Create Distance Matrix
vegdist   = vegdist(t(otu_table(oral.rel1)), method="bray")

##Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

##calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)

##Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))
#
##Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~adherent,data= newResults, mean)

centroids$adherent<-as.factor(centroids$adherent)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="adherent",suffixes=c("",".centroid"))

newResults$adherent<-as.character(newResults$adherent)

#plot it
pdf("Beta Diversity - adherent - MW.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=adherent)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Adherent"="#BAD7F2","Non-Adherent"="#96D09A"))+
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=adherent), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=adherent))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=adherent), parse=TRUE,size=10, show.legend=FALSE) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 40, colour = "black"),
        axis.title.y = element_text(size = 40, colour = "black"),
        axis.text.x = element_text(size = 40, colour = "black"), 
        axis.text.y = element_text(size = 40, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"))
dev.off()



#Create Table for Statistics    
data.adonis <- data.frame(sample_data(oral.rel))

#Run the Statistics
samplepermanova <- adonis(vegdist ~ adherent, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]
print(samplepermanova)
#p=0.469

###############################################################################
##Figure 6.11
###############################################################################


#Set alpha for differential Anlaysis, for adjusted p value 
alpha <- 0.2

#Convert To DESEQ, different type of phyloseq object for this program
diagdds <- phyloseq_to_deseq2(oral1, ~adherent)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds$adherent<-as.factor(diagdds$adherent)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$adherent<- droplevels(diagdds$adherent)

#Relevel Data
diagdds$adherent <- relevel(diagdds$adherent, ref ="Non-Adherent")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral1)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Create name with family and (u.g), creating a new column gs with names thats the highest known class
res$gs <- ifelse(is.na(res$Species),paste0(res$Genus,"_",res$ASV), paste0(res$Genus,"_",res$Species,"_",res$ASV))
res$gs <- ifelse(is.na(res$Genus),paste0(res$Family,"_",res$ASV), res$gs)
res$gs <- ifelse(is.na(res$Family), paste0(res$Order,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Order), paste0(res$Class, "_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Class), paste0(res$Phylum,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Phylum), paste0(res$Kingdom,"_",res$ASV),res$gs)

#Make the full trail the First Column, make this new column the 1st column. 
res$names <- res$ASV
res$Gene.symbol <- res$gs

#Subset the different categories, subset by smokers / COPD 
df.1 = subset_samples(oral.rel1, adherent %in% c("Non-Adherent"))
df.2 = subset_samples(oral.rel1, adherent  %in% c("Adherent"))

#decide what otu to save 
otu.to.save <-as.character(res$names)

#convert to dataframe
df.1.df <- data.frame(otu_table(df.1))
df.2.df <- data.frame(otu_table(df.2))

#from relative table we should get the mean across the row of the otu table, mean of the relative abundance. 
df.1.meanRA <- rowMeans(df.1.df)
df.2.meanRA <- rowMeans(df.2.df)

#need to subset AND reorder just the otus that we have 
df.1.meanRA.save <- df.1.meanRA[otu.to.save]
df.2.meanRA.save <- df.2.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance.1 <- df.1.meanRA.save
res$abundance.2 <- df.2.meanRA.save

#Keep only the count data
drops <- c("Domain","Phylum","Class","Order","Family","Genus","OTU","gs","Species",
           "ASV","Kingdom","names","row2")
res <- res[ , !(names(res) %in% drops)]

#Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL),
                c("ASV","baseMean", "logFC", "lfcSE", "stat", "pvalue", 
                  "adj.P.Val","Gene.symbol","abundance.1","abundance.2")) 

#Convert to data.frame
res <- as.data.frame(res)

#make an abundance variable for the size of the dots on the plot
res$abundance.2 <- as.numeric(as.character(res$abundance.2))
res$abundance.1 <- as.numeric(as.character(res$abundance.1))
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha, res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment direction to colour the graph
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Adherent", "Non-Adherent")

#plot it
pdf("Differential - adherence.pdf",width=35, height=7)
ggplot(sig_results, aes(x=Gene.symbol, y=logFC)) +
  geom_point(aes(fill=enriched, size=abundance, color = "black"), shape = 21) +
  geom_segment( aes(xend=Gene.symbol, yend=0)) +
  #geom_col(width=0.005, color="black") +
  #scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  #ylim(12,15) +
  scale_size_continuous(name="         Relative\n         Abundance", range = c(3, 20), ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  scale_fill_manual(values = c("Non-Adherent"="#96D09A","Adherent"="#BAD7F2"), guide = "none") +
  scale_color_manual(values = c("Non-Adherent"="#96D09A","Adherent"="#BAD7F2"), guide = "none") +
  xlab("") +
  ylim(-5,25) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 50, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()


###############################################################################
##figure 6.18
###############################################################################


#Set alpha for differential Anlaysis, for adjusted p value 
alpha <- 0.2

#Subset only variables for comparisons
oral1 = subset_samples(oral, patient_type1  %in% c("Post_Covid"))


#Subset only variables for comparisons
oral1 = subset_samples(oral, change.cfs  %in% c("Yes","No"))
oral1.rel <- transform_sample_counts(oral1, function(x) x/sum(x))

#Convert To DESEQ, different type of phyloseq object for this program
diagdds <- phyloseq_to_deseq2(oral1, ~change.cfs)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
#diagdds = estimateSizeFactors(diagdds)
diagdds$change.cfs<-as.factor(diagdds$change.cfs)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$change.cfs <- droplevels(diagdds$change.cfs)

#Relevel Data
diagdds$change.cfs <- relevel(diagdds$change.cfs, ref ="No")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral1)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)


#Create name with family and (u.g), creating a new column gs with names thats the highest known class
res$gs <- ifelse(is.na(res$Species),paste0(res$Genus,"_",res$ASV), paste0(res$Genus,"_",res$Species,"_",res$ASV))
res$gs <- ifelse(is.na(res$Genus),paste0(res$Family,"_",res$ASV), res$gs)
res$gs <- ifelse(is.na(res$Family), paste0(res$Order,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Order), paste0(res$Class, "_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Class), paste0(res$Phylum,"_",res$ASV),res$gs)

#Make the full trail the First Column, make this new column the 1st column. 
res$names <- res$ASV
res$Gene.symbol <- res$gs

#Subset the different categories, subset by smokers / COPD 
df.1 = subset_samples(oral1.rel, change.cfs %in% c("No"))
df.2 = subset_samples(oral1.rel, change.cfs  %in% c("Yes"))

#decide what otu to save 
otu.to.save <-as.character(res$names)

#convert to dataframe
df.1.df <- data.frame(otu_table(df.1))
df.2.df <- data.frame(otu_table(df.2))

#from relative table we should get the mean across the row of the otu table, mean of the relative abundance. 
df.1.meanRA <- rowMeans(df.1.df)
df.2.meanRA <- rowMeans(df.2.df)

#need to subset AND reorder just the otus that we have 
df.1.meanRA.save <- df.1.meanRA[otu.to.save]
df.2.meanRA.save <- df.2.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance.1 <- df.1.meanRA.save
res$abundance.2 <- df.2.meanRA.save

#Keep only the count data
drops <- c("Domain","Phylum","Class","Order","Family","Genus","OTU","gs","Species",
           "ASV","Kingdom","names","row2")
res <- res[ , !(names(res) %in% drops)]

#Set Names of Results Table
res <- setNames(cbind(rownames(res), res, row.names = NULL),
                c("ASV","baseMean", "logFC", "lfcSE", "stat", "pvalue", 
                  "adj.P.Val","Gene.symbol","abundance.1","abundance.2")) 

#Convert to data.frame
res <- as.data.frame(res)

#make an abundance variable for the size of the dots on the plot
res$abundance.2 <- as.numeric(as.character(res$abundance.2))
res$abundance.1 <- as.numeric(as.character(res$abundance.1))
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha, res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment direction to colour the graph
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Fatigue Resolved", "Persistantly Fatigued")

pdf("Differential - change in CFS - adherent group -just Post-COVID - updated1.pdf",width=35, height=7)
ggplot(sig_results, aes(x=Gene.symbol, y=logFC)) +
  geom_point(aes(fill = enriched, size = abundance),   shape = 21, color = "black")+
  scale_alpha_identity() +
  geom_segment( aes(xend=Gene.symbol, yend=0)) +
  #geom_col(width=0.005, color="black") +
  #scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  #ylim(12,15) +
  scale_size_continuous(name="         Relative\n         Abundance", range = c(3, 20), ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  scale_fill_manual(values = c("Fatigue Resolved"="#c9e4d3","Persistantly Fatigued"="#faedcb"), guide = "none") +
  scale_color_manual(values = c("Fatigue Resolved"="#c9e4d3","Persistantly Fatigued"="#faedcb"), guide = "none") +
  xlab("") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 50, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

