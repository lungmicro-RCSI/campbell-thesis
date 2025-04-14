#Chapter 5

#load library 
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
library(scales)
library(patchwork)
library(microbiome)
#import data 

#set wd 
setwd("")

#Taxa table
taxa <- read.csv(", row.names=1)


#OTU Table 
otu<-read.csv("",
              row.names = 1, header = T, check.names=F)


#Load MetaData
meta<-read.delim2("", sep="\t") 

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


#remove 0 abundance 
physeq = subset_taxa(physeq, rowSums(otu_table(physeq)) != 0)

#only visit 1
oral<-subset_samples(physeq, visit == "1")

#just oral 
oral<-subset_samples(oral, Description == "Mouth.Wash")

#make rel abundance table 
oral.rel <- transform_sample_counts(oral, function(x) x/sum(x))

#set theme 
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
#Figure 5.1A
###############################################################################

#Count up Reads & make a column of Reads
summary <- as.data.frame(rowSums(t(otu_table(oral))))

#Merge Reads with MetaData
Reads <- as.data.frame(merge( x = summary, y=sample_data(oral), by ="row.names", all.x = TRUE))

#Rename column from rowSums(t(otu_table(OTU)))) to Reads
colnames(Reads)[colnames(Reads)=="rowSums(t(otu_table(OTU))))"] <- "Reads"

colnames(Reads)[2]<-"Reads"

#Figure
pdf("Reads Count by Covid Status.pdf", height = 10, width = 15)
ggplot(Reads, aes(x=patient_type1, y=Reads, 
                  fill=patient_type1)) +
  stat_boxplot(geom = "errorbar", width=0.1) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(shape = 1, position=position_jitter(0.2)) +
  scale_y_continuous(name="Reads", trans="log10", breaks=trans_breaks('log10', function(x)10^x), 
                     labels=trans_format('log10', math_format(10^.x))) +
  scale_fill_manual(values=c("Post_COVID" = "#F0A4B0", "Incidental_COVID" = "#8ed4ac","Never_COVID" = "#F7D068")) + 
  guides(fill=FALSE) +
  xlab("COVID Status") +
  ylab("Reads") +
  theme
dev.off()

#check stats - overall 
sampleshannon <- kruskal.test(Reads ~ patient_type1, data = Reads)
sampleshannon <- sampleshannon$p.value

#Check stats for post vs incidental 
shannon1 <- subset(Reads, patient_type1 == "Post_COVID") 
shannon2 <- subset(Reads, patient_type1 == "Incidental_COVID") 
shannon3<-rbind(shannon1, shannon2)

sampleshannon2 <- kruskal.test(Reads ~ patient_type1, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value

#Check stats for post vs never
shannon4 <- subset(Reads, patient_type1 == "Post_COVID") 
shannon5 <- subset(Reads, patient_type1 == "Never_COVID") 
shannon6<-rbind(shannon4, shannon5)

sampleshannon3 <- kruskal.test(Reads ~ patient_type1, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value

#Check stats for Incidental vs never 
shannon7 <- subset(Reads, patient_type1 == "Incidental_COVID") 
shannon8 <- subset(Reads, patient_type1 == "Never_COVID") 
shannon9<-rbind(shannon7, shannon8)

sampleshannon4 <- kruskal.test(Reads ~ patient_type1, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value

###############################################################################
##Figure 5.1B
###############################################################################
#Calculates Shannon Diversity
sample_data(oral)$Shannon = diversity(otu_table(oral.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(oral))

#Make variable of interest a factor
Shannon$patient_type1 <- as.factor(Shannon$patient_type1)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#PLOT IT
pdf("Alpha Diversity - Whole Group by Covid Status - oral.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=patient_type1, y=Shannon, fill=patient_type1)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Post_COVID" = "#F0A4B0", "Incidental_COVID" = "#8ed4ac","Never_COVID" = "#F7D068")) + 
  ylab("Shannon Diversity") + 
  xlab("COVID Status") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

#check stats - overall 
sampleshannon <- kruskal.test(Shannon ~ patient_type1, data = Shannon)
sampleshannon <- sampleshannon$p.value


#Check stats for post vs incidental 
shannon1 <- subset(Shannon, patient_type1 == "Post_COVID") 
shannon2 <- subset(Shannon, patient_type1 == "Incidental_COVID") 
shannon3<-rbind(shannon1, shannon2)

sampleshannon2 <- kruskal.test(Shannon ~ patient_type1, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value

#Check stats for post vs never
shannon4 <- subset(Shannon, patient_type1 == "Post_COVID") 
shannon5 <- subset(Shannon, patient_type1 == "Never_COVID") 
shannon6<-rbind(shannon4, shannon5)

sampleshannon3 <- kruskal.test(Shannon ~ patient_type1, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value

#Check stats for Incidental vs never 
shannon7 <- subset(Shannon, patient_type1 == "Incidental_COVID") 
shannon8 <- subset(Shannon, patient_type1 == "Never_COVID") 
shannon9<-rbind(shannon7, shannon8)

sampleshannon4 <- kruskal.test(Shannon ~ patient_type1, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value

###############################################################################
##Figure 5.2
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
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"


##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~patient_type1,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="patient_type1",suffixes=c("",".centroid"))

#Figure
pdf("Beta Diversity - Whole Group - by pt type .pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, colour=patient_type1)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Post_COVID" = "#F0A4B0", "Incidental_COVID" = "#8ed4ac","Never_COVID" = "#F7D068")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=patient_type1), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=patient_type1))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=patient_type1), parse=TRUE,size=5, show.legend=FALSE) +
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

#Run the Statistics everyone 
samplepermanova <- adonis(vegdist ~ patient_type1, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]

#Check stats for Incidental vs never 
#subset samples
subset1<-subset_samples(oral.rel, patient_type1 %in% c("Incidental_COVID", "Never_COVID"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ patient_type1, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Check stats for post vs never 
#subset samples
subset1<-subset_samples(oral.rel, patient_type1 %in% c("Post_COVID", "Never_COVID"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ patient_type1, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Check stats for incidental vs post
#subset samples
subset1<-subset_samples(oral.rel, patient_type1 %in% c("Post_COVID", "Incidental_COVID"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ patient_type1, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

###############################################################################
##Figure 5.3
###############################################################################

##Never Covid

#subset for never covid
never.covid = subset_samples(nasal.rel, patient_type1  %in% c("Never_Covid"))

# Extract OTU table
otu <- otu_table(never.covid)
otu<-as.data.frame(otu)

# Calculate row means of taxa abundance
otu$MeanRA <- rowMeans(otu, na.rm = TRUE)

#sort my mean abundance 
otu_sorted <- otu %>%
  arrange(desc(MeanRA))

#take top 50
top10 <- otu_sorted %>%
  head(50)

#get genus names from tax table
tax<-tax_table(never.covid)
tax<-as.data.frame(tax)

#merge tax and top10
table<-merge(top10, tax, by="row.names", all.x=TRUE, all.y=FALSE)


#make a list of Genus 
genus<-list(table$Genus)

#make the list a vector 
genus <- unlist(genus)

#convert physeq to df
data<- data.table(psmelt(never.covid))

#subset so we only keep abundance data for the genus we want 
data1 <- data %>%
  filter(Genus %in% genus)


#cobvert Genus to factor and order

#data1$Genus <- factor(data1$Genus, levels = rev(genus))
data1$Genus<-as.factor(data1$Genus)

#log transform Abudance for ordering taxa
data1$logAbundance<-log(data1$Abundance)

#remove infinite values 
data1 <- data1 %>%
  filter(is.finite(logAbundance))


#order by Median 
data1 <- data1 %>%
  mutate(Genus = reorder(Genus, logAbundance, FUN = median))


#plot

never.covid<-ggplot(data1, aes(x=Abundance, y=Genus, fill = Genus)) +
  stat_boxplot(geom ='errorbar', width=0.1, color = "black")+
  geom_boxplot(outlier.shape = NA, width=0.5, fill = "#F7D068", color = "black")+
  geom_jitter(shape=1, position=position_jitter(0.2), size=0.7)+
  ylab("") + 
  xlab("") +
  ggtitle("         Never Covid") +
  scale_x_continuous(name="Relative Abundance",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  guides(fill = "none") 


##Incidental Covid

#subset for incidental 
incidental = subset_samples(nasal.rel, patient_type1  %in% c("Incidental_Covid"))

# Extract OTU table
otu1 <- otu_table(incidental)
otu1<-as.data.frame(otu1)

# Calculate row means of taxa abundance
otu1$MeanRA <- rowMeans(otu1, na.rm = TRUE)

#sort my mean abundance 
otu_sorted1 <- otu1 %>%
  arrange(desc(MeanRA))


#take top 50
top10.1 <- otu_sorted1 %>%
  head(50)

#get genus names from tax table
tax1<-tax_table(incidental)
tax1<-as.data.frame(tax1)

#merge tax and top50
table1<-merge(top10.1, tax1, by="row.names", all.x=TRUE, all.y=FALSE)


#make a list of Genus 
genus1<-list(table1$Genus)

#make the list a vector 
genus1 <- unlist(genus1)

#convert physeq to df
df<- data.table(psmelt(incidental))

#subset so we only keep abundance data for the genus we want 
df <- df %>%
  filter(Genus %in% genus1)


#cobvert Genus to factor
df$Genus<-as.factor(df$Genus)

#log transform Abudance for plot 
df$logAbundance<-log10(df$Abundance)

#remove infinite values 
df <- df %>%
  filter(is.finite(logAbundance))


#order by Median 
df <- df %>%
  mutate(Genus = reorder(Genus, logAbundance, FUN = median))


#plot

incidental<-ggplot(df, aes(x=Abundance, y=Genus, fill = Genus)) +
  stat_boxplot(geom ='errorbar', width=0.1, color = "black")+
  geom_boxplot(outlier.shape = NA, width=0.5, fill = "#8ed4ac", color = "black")+
  geom_jitter(shape=1, position=position_jitter(0.2), size=0.7)+
  ylab("") + 
  xlab("") +
  ggtitle("   Incidental COVID") +
  scale_x_continuous(name="Relative Abundance",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  guides(fill = "none") 


##PostCOVID

#subset for post
post = subset_samples(nasal.rel, patient_type1  %in% c("Post_Covid"))

# Extract OTU table
otu2 <- otu_table(post)
otu2<-as.data.frame(otu2)

# Calculate row means of taxa abundance
otu2$MeanRA <- rowMeans(otu2, na.rm = TRUE)

#sort my mean abundance 
otu_sorted2 <- otu2 %>%
  arrange(desc(MeanRA))


#take top 50
top10.2 <- otu_sorted2 %>%
  head(50)

#get genus names from tax table
tax2<-tax_table(post)
tax2<-as.data.frame(tax2)

#merge tax and top50
table2<-merge(top10.2, tax2, by="row.names", all.x=TRUE, all.y=FALSE)


#make a list of Genus 
genus2<-list(table2$Genus)

#make the list a vector 
genus2 <- unlist(genus2)

#convert physeq to df
df1<- data.table(psmelt(post))

#subset so we only keep abundance data for the genus we want 
df1 <- df1 %>%
  filter(Genus %in% genus2)


#convert Genus to factor
df1$Genus<-as.factor(df1$Genus)

#log transform Abudance for plot 
df1$logAbundance<-log10(df1$Abundance)

#remove infinite values 
df1 <- df1 %>%
  filter(is.finite(logAbundance))


#order by Median 
df1 <- df1 %>%
  mutate(Genus = reorder(Genus, logAbundance, FUN = median))


#plot

post<-ggplot(df1, aes(x=Abundance, y=Genus, fill = Genus)) +
  stat_boxplot(geom ='errorbar', width=0.1, color = "black")+
  geom_boxplot(outlier.shape = NA, width=0.5, fill = "#F0A4B0", color = "black")+
  geom_jitter(shape=1, position=position_jitter(0.2), size=0.7)+
  ylab("") + 
  xlab("") +
  ggtitle("        Post-COVID") +
  scale_x_continuous(name="Relative Abundance",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  guides(fill = "none") 

plot<-never.covid+incidental+post

pdf("Top_50_Taxa_by_COVID_Status.pdf", width = 10, height = 7)
print(plot)
dev.off()


###############################################################################
##Figure 5.4
###############################################################################
#Comparison A

#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, patient_type1  %in% c("Never_COVID","Incidental_COVID"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~patient_type1)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$patient_type1<- droplevels(diagdds$patient_type1)

#Relevel Data
diagdds$patient_type1 <- relevel(diagdds$patient_type1, ref ="Never_COVID")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]


#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#make column with names
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
df.1 = subset_samples(oral.rel, patient_type1  %in% c("Never_COVID"))
df.2 = subset_samples(oral.rel, patient_type1  %in% c("Incidental_COVID"))

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
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, NA))


# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Incidental_COVID", "Never_COVID")

#compa never vs incidenceal 
compA<-sig_results

#make a column for group 
compA$group <- "A"

#Comparison B

#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, patient_type1  %in% c("Post_COVID","Incidental_COVID"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~patient_type1)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$patient_type1<- droplevels(diagdds$patient_type1)

#Relevel Data
diagdds$patient_type1 <- relevel(diagdds$patient_type1, ref ="Incidental_COVID")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]


#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
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
df.1 = subset_samples(oral.rel, patient_type1  %in% c("Incidental_COVID"))
df.2 = subset_samples(oral.rel, patient_type1  %in% c("Post_COVID"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha,  res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Post_COVID", "Incidental_COVID")


#compB
compB<-sig_results

#make a column for group 
compB$group <- "B"

#Comparison C 

#Set alpha for differential Anlaysis, for adjusted p value 
alpha <- 0.2

#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, patient_type1  %in% c("Post_COVID","Never_COVID"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~patient_type1)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$patient_type1<- droplevels(diagdds$patient_type1)

#Relevel Data
diagdds$patient_type1 <- relevel(diagdds$patient_type1, ref ="Never_COVID")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#make column with names
res$gs <- ifelse(is.na(res$Species),paste0(res$Genus,"_",res$ASV), paste0(res$Genus,"_",res$Species,"_",res$ASV))
res$gs <- ifelse(is.na(res$Genus),paste0(res$Family,"_",res$ASV), res$gs)
res$gs <- ifelse(is.na(res$Family), paste0(res$Order,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Order), paste0(res$Class, "_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Class), paste0(res$Phylum,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Phylum), paste0(res$Kingdom,"_",res$ASV),res$gs)


#Make the full trail the First Column, make this new column the 1st column. 
res$names <- res$ASV
res$Gene.symbol <- res$gs

#Subset the different categories, subset by patient type 
df.1 = subset_samples(oral.rel, patient_type1  %in% c("Never_COVID"))
df.2 = subset_samples(oral.rel, patient_type1  %in% c("Post_COVID"))

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

#Convert to data.frame
res <- as.data.frame(res)

#make an abundance variable for the size of the dots on the plot
res$abundance.2 <- as.numeric(as.character(res$abundance.2))
res$abundance.1 <- as.numeric(as.character(res$abundance.1))
res$abundance <- ifelse(res$logFC>=1, res$abundance.2, 
                        ifelse(res$logFC<=-1, res$abundance.1, NA))


# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Post_COVID", "Never_COVID")

#compc 
compC<-sig_results

#make a column for group 
compC$group <- "C"

#merge df
df <- rbind(compA, compB, compC)

#set order by abundance
df <- df %>%
  arrange(factor(group, levels = c("A", "B", "C")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - COVID Status.pdf",width=45, height=20) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 15)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x") +  
  scale_fill_manual(values = c("Incidental_COVID"= "#8ed4ac", "Never_COVID" = "#F7D068", "Post_COVID" = "#F0A4B0", "notsig" = "white"), guide="none") +
  xlab("LogFC") +
  ylab("") +
  xlim(-30, 30) + 
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
##Figure 5.5A & 5.7A
##############################################################################

#Count up Reads & make a column of Reads
summary <- as.data.frame(rowSums(t(otu_table(oral))))

#Merge Reads with MetaData
Reads <- as.data.frame(merge( x = summary, y=sample_data(oral), by ="row.names", all.x = TRUE))

#Rename column from rowSums(t(otu_table(OTU)))) to Reads
colnames(Reads)[colnames(Reads)=="rowSums(t(otu_table(OTU))))"] <- "Reads"

colnames(Reads)[2]<-"Reads"

typeof(Reads$Description)

#Create figure
pdf("Reads Count by OSA.pdf", height = 10, width = 15)
ggplot(Reads, aes(x=OSA, y=Reads, 
                  fill=OSA)) +
  stat_boxplot(geom = "errorbar", width=0.1) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(shape = 1, position=position_jitter(0.2)) +
  scale_y_continuous(name="Reads", trans="log10", breaks=trans_breaks('log10', function(x)10^x), 
                     labels=trans_format('log10', math_format(10^.x))) +
  scale_fill_manual(values=c("OSA"="#D95980",
                             "Normal"="#63AAC0")) +
  guides(fill=FALSE) +
  xlab("OSA Severity") +
  ylab("Reads") +
  theme
dev.off()

#check stats - overall 
sampleshannon <- kruskal.test(Reads ~ OSA, data = Reads)
sampleshannon <- sampleshannon$p.value

#create figure
pdf("Reads Count by OSA.pdf", height = 10, width = 15)
ggplot(Reads, aes(x=severity, y=Reads, 
                  fill=severity)) +
  stat_boxplot(geom = "errorbar", width=0.1) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(shape = 1, position=position_jitter(0.2)) +
  scale_y_continuous(name="Reads", trans="log10", breaks=trans_breaks('log10', function(x)10^x), 
                     labels=trans_format('log10', math_format(10^.x))) +
  scale_fill_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75")) +
  guides(fill=FALSE) +
  xlab("OSA Severity") +
  ylab("Reads") +
  theme
dev.off()


#Check stats for Normal vs Mild 
shannon1 <- subset(Reads, severity == "Normal") 
shannon2 <- subset(Reads, severity == "Mild") 
shannon3<-rbind(shannon1, shannon2)
droplevels(shannon3$severity)

sampleshannon2 <- kruskal.test(Reads ~ severity, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value
print(sampleshannon2)

#Check stats for Normal vs Moderate 
shannon4 <- subset(Reads, severity == "Normal") 
shannon5 <- subset(Reads, severity == "Moderate") 
shannon6<-rbind(shannon4, shannon5)
droplevels(shannon6$severity)


sampleshannon3 <- kruskal.test(Reads ~ severity, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value

#Check stats for Normal vs Severe 
shannon7 <- subset(Reads, severity == "Normal") 
shannon8 <- subset(Reads, severity == "Severe") 
shannon9<-rbind(shannon7, shannon8)
droplevels(shannon9$severity)


sampleshannon4 <- kruskal.test(Reads ~ severity, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value
print(sampleshannon3)

#Check stats for Mild vs Mod
shannon10 <- subset(Reads, severity == "Mild") 
shannon11 <- subset(Reads, severity == "Moderate") 
shannon12<-rbind(shannon10, shannon11)
droplevels(shannon12$severity)

sampleshannon5 <- kruskal.test(Reads ~ severity, data = shannon12)
sampleshannon5 <- sampleshannon5$p.value
print(sampleshannon5)


#Check stats for Mild vs Severe
shannon13 <- subset(Reads, severity == "Mild") 
shannon14 <- subset(Reads, severity == "Severe") 
shannon15<-rbind(shannon13, shannon14)
droplevels(shannon15$severity)

sampleshannon6 <- kruskal.test(Reads ~ severity, data = shannon15)
sampleshannon6 <- sampleshannon6$p.value
print(sampleshannon6)


#Check stats for Moderate vs Severe
shannon16 <- subset(Reads, severity == "Moderate") 
shannon17 <- subset(Reads, severity == "Severe") 
shannon18<-rbind(shannon16, shannon17)
droplevels(shannon18$severity)


sampleshannon7 <- kruskal.test(Reads ~ severity, data = shannon18)
sampleshannon7 <- sampleshannon7$p.value
print(sampleshannon7)


###############################################################################
##Figure 5.5B & 5.7A
###############################################################################

#Calculates Shannon Diversity
sample_data(oral)$Shannon = diversity(otu_table(oral.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(oral))

#Make variable of interest a factor
Shannon$OSA <- as.factor(Shannon$OSA)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#plot
pdf("Alpha Diversity - Whole Group - OSA AHI 5 - oral.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=OSA, y=Shannon, fill=OSA)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("OSA"="#D95980","Normal"="#63AAC0"), guide = "none") + 
  ylab("Shannon Diversity") + 
  xlab("OSA") +
  theme 
dev.off()

#check stats
sampleshannon <- kruskal.test(Shannon ~ OSA, data = Shannon)
sampleshannon <- sampleshannon$p.value

#plot
pdf("Alpha Diversity - Whole Group - OSA Severity - oral.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=severity, y=Shannon, fill=severity)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75")) + 
  ylab("Shannon Diversity") + 
  xlab("Severity") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

#check stats
sampleshannon <- kruskal.test(Shannon ~ severity, data = Shannon)
sampleshannon <- sampleshannon$p.value


#Check stats for Normal vs Mild 
shannon1 <- subset(Shannon, severity == "Normal") 
shannon2 <- subset(Shannon, severity == "Mild") 
shannon3<-rbind(shannon1, shannon2)

sampleshannon2 <- kruskal.test(Shannon ~ severity, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value

#Check stats for Normal vs Moderate 
shannon4 <- subset(Shannon, severity == "Normal") 
shannon5 <- subset(Shannon, severity == "Moderate") 
shannon6<-rbind(shannon4, shannon5)

sampleshannon3 <- kruskal.test(Shannon ~ severity, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value

#Check stats for Normal vs Severe 
shannon7 <- subset(Shannon, severity == "Normal") 
shannon8 <- subset(Shannon, severity == "Severe") 
shannon9<-rbind(shannon7, shannon8)

sampleshannon4 <- kruskal.test(Shannon ~ severity, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value


#Check stats for Mild vs Mod
shannon10 <- subset(Shannon, severity == "Mild") 
shannon11 <- subset(Shannon, severity == "Moderate") 
shannon12<-rbind(shannon10, shannon11)

sampleshannon5 <- kruskal.test(Shannon ~ severity, data = shannon12)
sampleshannon5 <- sampleshannon5$p.value

#Check stats for Mild vs Severe
shannon13 <- subset(Shannon, severity == "Mild") 
shannon14 <- subset(Shannon, severity == "Severe") 
shannon15<-rbind(shannon13, shannon14)

sampleshannon6 <- kruskal.test(Shannon ~ severity, data = shannon15)
sampleshannon6 <- sampleshannon6$p.value


#Check stats for Moderate vs Severe
shannon16 <- subset(Shannon, severity == "Moderate") 
shannon17 <- subset(Shannon, severity == "Severe") 
shannon18<-rbind(shannon16, shannon17)

sampleshannon7 <- kruskal.test(Shannon ~ severity, data = shannon18)
sampleshannon7 <- sampleshannon7$p.value


###############################################################################
##Figure 5.6
###############################################################################

##Create Distance Matrix
vegdist   = vegdist(t(otu_table(oral.rel)), method="bray")

##Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

##calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)

##Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

##Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~OSA,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="OSA",suffixes=c("",".centroid"))


#Fig 4.9
pdf("Beta Diversity - Whole group - OSA AHI 5 - scaled by AHI.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=OSA)) +
  geom_point(aes(size=ahi, fill=OSA, color=OSA),alpha=0.5, shape=21) +
  scale_size_continuous(range = c(5, 20), name= "         AHI") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("OSA"="#D95980","Normal"="#63AAC0"), guide=FALSE) + 
  scale_fill_manual(values=c("OSA"="#D95980","Normal"="#63AAC0"), guide=FALSE) +
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=OSA), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=OSA), size=0.5)+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=OSA), parse=TRUE,size=10, show.legend=FALSE) +
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
samplepermanova <- adonis(vegdist ~ OSA, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]

###############################################################################
##Figure 5.8
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
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~severity,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="severity",suffixes=c("",".centroid"))

pdf("Beta Diversity - Whole group - OSA Severity.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=severity)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75"), guide = "none") + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=severity), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=severity))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=severity), parse=TRUE,size=5, show.legend=FALSE) +
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
samplepermanova <- adonis(vegdist ~ severity, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]

#Normal vs mild 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Mild"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Moderate
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Moderate 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Moderate vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Moderate", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]


###############################################################################
##Figure 5.9
###############################################################################

#Comparison A
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Normal","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Severe", "Normal")

#compA normal vs mild
compA<-sig_results
compA$group <- "A"

#Comparison B 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Mild","Moderate"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Moderate"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results1 <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results1 = sig_results1[order(sig_results1$logFC, na.last = NA), ]
sig_results1$Gene.symbol <- factor(sig_results1$Gene.symbol, levels = sig_results1$Gene.symbol[order(sig_results1$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results1$enriched<-ifelse(sig_results1$logFC>0, "Moderate", "Mild")

#compB norm vs moderate 
compB<-sig_results1
compB$group <- "B"

#Comparason C

Comp1.OTU.Rel.Table = subset_samples(oral.rel, severity  %in% c("Mild","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Severe", "Mild")

#compC 
compC<-sig_results

compC$group<-"C"


#Comparison D
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Moderate","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Moderate")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Moderate"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Severe", "Moderate")

#compD
compD<-sig_results
compD$group <- "D"

#combine variables 
df <- rbind(compA, compB, compC, compD)

#reodrder 
df <- df %>%
  arrange(factor(group, levels = c("A", "B", "C", "D")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - normal comps - xlims.pdf",width=30, height=7) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 10)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75", "Severe OSA" = "#832c75"), guide="none") +
xlab("LogFC") +
  ylab("") + 
  xlim(-30, 25) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

###############################################################################
##RDI Analysis 
##############################################################################
#load metadata 
sleepy.nr <- read.delim("")
rownames(sleepy.nr) <- sleep$sample.id

#make OSA variables using RDI 
#make a variable for OSA severity 
sleepy.nr$severity <-ifelse(sleepy.nr$rdi <=5, "Normal", 
                            ifelse(sleepy.nr$rdi <15, "Mild",
                                   ifelse(sleepy.nr$rdi <30, "Moderate", "Severe")))


#Make a variable for OSA 
sleepy.nr$OSA<-ifelse(sleepy.nr$rdi >=5, "OSA", "No OSA")
sleepy.nr$OSA<-as.factor(as.character(sleepy.nr$OSA))

#convert to metadata format
meta<-sample_data(sleepy.nr)

#make phyloseq
oral=phyloseq(otu, TAX, meta)

#only visit 1
oral<-subset_samples(oral, visit == "1")

#just oral 
oral<-subset_samples(oral, Description == "Mouth.Wash")

#remove 0 abundance 
oral = subset_taxa(oral, rowSums(otu_table(oral)) != 0)

#make rel abundance table 
oral.rel <- transform_sample_counts(oral, function(x) x/sum(x))

###############################################################################
##Figure 5.10A 
##############################################################################

#Count up Reads & make a column of Reads
summary <- as.data.frame(rowSums(t(otu_table(oral))))


#Merge Reads with MetaData
Reads <- as.data.frame(merge( x = summary, y=sample_data(oral), by ="row.names", all.x = TRUE))

#Rename column from rowSums(t(otu_table(OTU)))) to Reads
colnames(Reads)[colnames(Reads)=="rowSums(t(otu_table(OTU))))"] <- "Reads"

colnames(Reads)[2]<-"Reads"

typeof(Reads$Description)

#Create figure
pdf("Reads Count by OSA RDI.pdf", height = 10, width = 15)
ggplot(Reads, aes(x=OSA, y=Reads, 
                  fill=OSA)) +
  stat_boxplot(geom = "errorbar", width=0.1) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(shape = 1, position=position_jitter(0.2)) +
  scale_y_continuous(name="Reads", trans="log10", breaks=trans_breaks('log10', function(x)10^x), 
                     labels=trans_format('log10', math_format(10^.x))) +
  scale_fill_manual(values=c("OSA"="#D95980",
                             "Normal"="#63AAC0")) +
  guides(fill=FALSE) +
  xlab("OSA Severity") +
  ylab("Reads") +
  theme
dev.off()

#check stats - overall 
sampleshannon <- kruskal.test(Reads ~ OSA, data = Reads)
sampleshannon <- sampleshannon$p.value


###############################################################################
##Figure 5.10B & 5.12A
###############################################################################

#Calculates Shannon Diversity
sample_data(oral)$Shannon = diversity(otu_table(oral.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(oral))

#Make variable of interest a factor
Shannon$OSA <- as.factor(Shannon$OSA)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#Fig 4.8B 
pdf("Alpha Diversity - Whole Group - OSA RDI 5 - oral.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=OSA, y=Shannon, fill=OSA)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("OSA"="#D95980","Normal"="#63AAC0"), guide = "none") + 
  ylab("Shannon Diversity") + 
  xlab("OSA") +
  theme 
dev.off()

#check stats
sampleshannon <- kruskal.test(Shannon ~ OSA, data = Shannon)
sampleshannon <- sampleshannon$p.value

#Plot it 
pdf("Alpha Diversity - Whole Group - OSA RDI Severity - oral.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=severity, y=Shannon, fill=severity)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75")) + 
  ylab("Shannon Diversity") + 
  xlab("Severity") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

#check stats
sampleshannon <- kruskal.test(Shannon ~ severity, data = Shannon)
sampleshannon <- sampleshannon$p.value


#Check stats for Normal vs Mild 
shannon1 <- subset(Shannon, severity == "Normal") 
shannon2 <- subset(Shannon, severity == "Mild") 
shannon3<-rbind(shannon1, shannon2)

sampleshannon2 <- kruskal.test(Shannon ~ severity, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value

#Check stats for Normal vs Moderate 
shannon4 <- subset(Shannon, severity == "Normal") 
shannon5 <- subset(Shannon, severity == "Moderate") 
shannon6<-rbind(shannon4, shannon5)

sampleshannon3 <- kruskal.test(Shannon ~ severity, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value

#Check stats for Normal vs Severe 
shannon7 <- subset(Shannon, severity == "Normal") 
shannon8 <- subset(Shannon, severity == "Severe") 
shannon9<-rbind(shannon7, shannon8)

sampleshannon4 <- kruskal.test(Shannon ~ severity, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value


#Check stats for Mild vs Mod
shannon10 <- subset(Shannon, severity == "Mild") 
shannon11 <- subset(Shannon, severity == "Moderate") 
shannon12<-rbind(shannon10, shannon11)

sampleshannon5 <- kruskal.test(Shannon ~ severity, data = shannon12)
sampleshannon5 <- sampleshannon5$p.value

#Check stats for Mild vs Severe
shannon13 <- subset(Shannon, severity == "Mild") 
shannon14 <- subset(Shannon, severity == "Severe") 
shannon15<-rbind(shannon13, shannon14)

sampleshannon6 <- kruskal.test(Shannon ~ severity, data = shannon15)
sampleshannon6 <- sampleshannon6$p.value


#Check stats for Moderate vs Severe
shannon16 <- subset(Shannon, severity == "Moderate") 
shannon17 <- subset(Shannon, severity == "Severe") 
shannon18<-rbind(shannon16, shannon17)

sampleshannon7 <- kruskal.test(Shannon ~ severity, data = shannon18)
sampleshannon7 <- sampleshannon7$p.value


###############################################################################
##Figure 5.11
###############################################################################

##Create Distance Matrix
vegdist   = vegdist(t(otu_table(oral.rel)), method="bray")

##Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

##calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)

##Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

##Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~OSA,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="OSA",suffixes=c("",".centroid"))


#Plot it 
pdf("Beta Diversity - Whole group - OSA RDI 5 - scaled by RDI.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=OSA)) +
  geom_point(aes(size=rdi, fill=OSA, color=OSA),alpha=0.5, shape=21) +
  scale_size_continuous(range = c(5, 20), name= "         RDI") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("OSA"="#D95980","Normal"="#63AAC0"), guide=FALSE) + 
  scale_fill_manual(values=c("OSA"="#D95980","Normal"="#63AAC0"), guide=FALSE) +
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=OSA), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=OSA), size=0.5)+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=OSA), parse=TRUE,size=10, show.legend=FALSE) +
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
samplepermanova <- adonis(vegdist ~ OSA, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]



###############################################################################
##Figure 5.12B
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
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~severity,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="severity",suffixes=c("",".centroid"))

#plot it
pdf("Beta Diversity - Whole group - OSA RDI Severity.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=severity)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75"), guide = "none") + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=severity), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=severity))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=severity), parse=TRUE,size=5, show.legend=FALSE) +
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
samplepermanova <- adonis(vegdist ~ severity, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]

#Normal vs mild 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Mild"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Moderate
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Moderate 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Moderate vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Moderate", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

###############################################################################
##Figure 4.13
###############################################################################

#Comparison A
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Normal","Mild"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Mild"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Mild", "Normal")

#compA 
compA<-sig_results
compA$group <- "A"

#Comparason B 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Normal","Moderate"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Moderate"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results1 <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results1 = sig_results1[order(sig_results1$logFC, na.last = NA), ]
sig_results1$Gene.symbol <- factor(sig_results1$Gene.symbol, levels = sig_results1$Gene.symbol[order(sig_results1$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results1$enriched<-ifelse(sig_results1$logFC>0, "Moderate", "Normal")

#compB 
compB<-sig_results1
compB$group <- "B"

#Comparason C
Comp1.OTU.Rel.Table = subset_samples(oral.rel, severity  %in% c("Normal","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Severe", "Normal")

#compC 
compC<-sig_results

compC$group<-"C"

#combine variables 
df <- rbind(compA, compB, compC)

#reodrder 
df <- df %>%
  arrange(factor(group, levels = c("A", "B", "C")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - normal comps - xlims - RDI.pdf",width=30, height=7) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 10)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75", "Severe OSA" = "#832c75"), guide="none") +
xlab("LogFC") +
  ylab("") + 
  xlim(-30, 25) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()


#Comparison D
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Mild","Moderate"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Moderate"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Moderate", "Mild")

#compD
compD<-sig_results
compD$group <- "D"

#Comparison E
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Mild","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results1 <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results1 = sig_results1[order(sig_results1$logFC, na.last = NA), ]
sig_results1$Gene.symbol <- factor(sig_results1$Gene.symbol, levels = sig_results1$Gene.symbol[order(sig_results1$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results1$enriched<-ifelse(sig_results1$logFC>0, "Severe", "Mild")

#compB norm vs moderate 
compE<-sig_results1
compE$group <- "E"

#Comparason F
Comp1.OTU.Rel.Table = subset_samples(oral.rel, severity  %in% c("Moderate","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Moderate")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Moderate"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Severe", "Moderate")

#compF
compF<-sig_results

compF$group<-"F"

#combine variables 
df <- rbind(compD, compE, compF)

#reodrder 
df <- df %>%
  arrange(factor(group, levels = c("D", "E", "F")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - normal comps - xlims - RDI.pdf",width=30, height=7) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 10)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75", "Severe OSA" = "#832c75"), guide="none") +
xlab("LogFC") +
  ylab("") + 
  xlim(-30, 25) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

###############################################################################
##ODI Analysis 
##############################################################################
#load metadata 
sleepy.nr <- read.delim("")
rownames(sleepy.nr) <- sleep$sample.id

#make OSA variables using ODI
#make a variable for OSA severity 
sleepy.nr$severity <-ifelse(sleepy.nr$odi <=5, "Normal", 
                            ifelse(sleepy.nr$odi <15, "Mild",
                                   ifelse(sleepy.nr$odi <30, "Moderate", "Severe")))


#Make a variable for OSA 
sleepy.nr$OSA<-ifelse(sleepy.nr$odi >=5, "OSA", "No OSA")
sleepy.nr$OSA<-as.factor(as.character(sleepy.nr$OSA))

#convert to metadata format
meta<-sample_data(sleepy.nr)

#make phyloseq
oral=phyloseq(otu, TAX, meta)

#only visit 1
oral<-subset_samples(oral, visit == "1")

#just oral 
oral<-subset_samples(oral, Description == "Mouth.Wash")

#remove 0 abundance 
oral = subset_taxa(oral, rowSums(otu_table(oral)) != 0)

#make rel abundance table 
oral.rel <- transform_sample_counts(oral, function(x) x/sum(x))


###############################################################################
##Figure 5.14A
###############################################################################


#plot it 
pdf("Alpha Diversity - Whole Group - OSA ODI Severity - oral.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=severity, y=Shannon, fill=severity)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75")) + 
  ylab("Shannon Diversity") + 
  xlab("Severity") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

#check stats
sampleshannon <- kruskal.test(Shannon ~ severity, data = Shannon)
sampleshannon <- sampleshannon$p.value


#Check stats for Normal vs Mild 
shannon1 <- subset(Shannon, severity == "Normal") 
shannon2 <- subset(Shannon, severity == "Mild") 
shannon3<-rbind(shannon1, shannon2)

sampleshannon2 <- kruskal.test(Shannon ~ severity, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value

#Check stats for Normal vs Moderate 
shannon4 <- subset(Shannon, severity == "Normal") 
shannon5 <- subset(Shannon, severity == "Moderate") 
shannon6<-rbind(shannon4, shannon5)

sampleshannon3 <- kruskal.test(Shannon ~ severity, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value

#Check stats for Normal vs Severe 
shannon7 <- subset(Shannon, severity == "Normal") 
shannon8 <- subset(Shannon, severity == "Severe") 
shannon9<-rbind(shannon7, shannon8)

sampleshannon4 <- kruskal.test(Shannon ~ severity, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value


#Check stats for Mild vs Mod
shannon10 <- subset(Shannon, severity == "Mild") 
shannon11 <- subset(Shannon, severity == "Moderate") 
shannon12<-rbind(shannon10, shannon11)

sampleshannon5 <- kruskal.test(Shannon ~ severity, data = shannon12)
sampleshannon5 <- sampleshannon5$p.value

#Check stats for Mild vs Severe
shannon13 <- subset(Shannon, severity == "Mild") 
shannon14 <- subset(Shannon, severity == "Severe") 
shannon15<-rbind(shannon13, shannon14)

sampleshannon6 <- kruskal.test(Shannon ~ severity, data = shannon15)
sampleshannon6 <- sampleshannon6$p.value


#Check stats for Moderate vs Severe
shannon16 <- subset(Shannon, severity == "Moderate") 
shannon17 <- subset(Shannon, severity == "Severe") 
shannon18<-rbind(shannon16, shannon17)

sampleshannon7 <- kruskal.test(Shannon ~ severity, data = shannon18)
sampleshannon7 <- sampleshannon7$p.value


###############################################################################
##Figure 5.14B
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
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~severity,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="severity",suffixes=c("",".centroid"))

#plot it 
pdf("Beta Diversity - Whole group - OSA ODI Severity.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=severity)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75"), guide = "none") + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=severity), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=severity))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=severity), parse=TRUE,size=5, show.legend=FALSE) +
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
samplepermanova <- adonis(vegdist ~ severity, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]

#Normal vs mild 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Mild"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Moderate
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Moderate 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Moderate vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Moderate", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

###############################################################################
##Figure 5.15
###############################################################################

#Comparison A
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Normal","Moderate"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Moderate"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Mild", "Normal")

#compA 
compA<-sig_results
compA$group <- "A"

#Comparason B 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Normal","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results1 <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results1 = sig_results1[order(sig_results1$logFC, na.last = NA), ]
sig_results1$Gene.symbol <- factor(sig_results1$Gene.symbol, levels = sig_results1$Gene.symbol[order(sig_results1$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results1$enriched<-ifelse(sig_results1$logFC>0, "Moderate", "Normal")

#compB 
compB<-sig_results1
compB$group <- "B"


#combine variables 
df <- rbind(compA, compB)

#reodrder 
df <- df %>%
  arrange(factor(group, levels = c("A", "B")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - normal comps - xlims - ODI.pdf",width=30, height=7) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 10)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75", "Severe OSA" = "#832c75"), guide="none") +
xlab("LogFC") +
  ylab("") + 
  xlim(-30, 25) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

#Comparison C
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Mild","Moderate"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Moderate"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Moderate", "Mild")

#compC
compC<-sig_results
compC$group <- "C"

#Comparison D 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Mild","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results1 <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results1 = sig_results1[order(sig_results1$logFC, na.last = NA), ]
sig_results1$Gene.symbol <- factor(sig_results1$Gene.symbol, levels = sig_results1$Gene.symbol[order(sig_results1$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results1$enriched<-ifelse(sig_results1$logFC>0, "Severe", "Mild")

#compD 
compD<-sig_results1
compD$group <- "D"

#Comparison E
Comp1.OTU.Rel.Table = subset_samples(oral.rel, severity  %in% c("Moderate","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Moderate")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Moderate"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Severe", "Moderate")

#compE 
compE<-sig_results

compE$group<-"E"

#combine variables 
df <- rbind(compC, compD, compE)

#reorder
df <- df %>%
  arrange(factor(group, levels = c("C", "D", "E")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - normal comps - xlims - ODI.pdf",width=30, height=7) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 10)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75", "Severe OSA" = "#832c75"), guide="none") +
xlab("LogFC") +
  ylab("") + 
  xlim(-30, 25) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()



###############################################################################
##Figure 5.16A
###############################################################################

#Calculates Shannon Diversity
sample_data(oral.rel)$Shannon = diversity(otu_table(oral.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(oral.rel))

#Make variable of interest a factor
Shannon$t90x <- as.factor(Shannon$t90x)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#make figure 
pdf("Alpha Diversity - Everyone - Hypoxia - MW.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=t90x, y=Shannon, fill=t90x)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Light"="#bae1ff","Mild"="#00C9B7","Moderate/Severe"="#F7C7FF"))+
  ylab("Shannon Diversity") + 
  xlab("Hypoxia") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

#check stats - overall 
sampleshannon <- kruskal.test(Shannon ~ t90x, data = Shannon)
sampleshannon <- sampleshannon$p.value

#Check stats for Light vs Mild 
shannon1 <- subset(Shannon, t90x == "Light") 
shannon2 <- subset(Shannon, t90x == "Mild") 
shannon3<-rbind(shannon1, shannon2)

sampleshannon2 <- kruskal.test(Shannon ~ t90x, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value

#Check stats for Light vs Mod/Severe
shannon4 <- subset(Shannon, t90x == "Light") 
shannon5 <- subset(Shannon, t90x == "Moderate/Severe") 
shannon6<-rbind(shannon4, shannon5)

sampleshannon3 <- kruskal.test(Shannon ~ t90x, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value

#Check stats for Mild vs Mod/Severe
shannon7 <- subset(Shannon, t90x == "Mild") 
shannon8 <- subset(Shannon, t90x == "Moderate/Severe") 
shannon9<-rbind(shannon7, shannon8)

sampleshannon4 <- kruskal.test(Shannon ~ t90x, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value


###############################################################################
##Figure 5.16B
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
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~t90x,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="t90x",suffixes=c("",".centroid"))

#plot it 
pdf("Beta Diversity - Whole group - Hypoxia - MW.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=t90x)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Light"="#bae1ff","Mild"="#00C9B7","Moderate/Severe"="#F7C7FF"), guide="none")+
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=t90x), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=t90x))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=t90x), parse=TRUE,size=5, show.legend=FALSE) +
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
samplepermanova <- adonis(vegdist ~ t90x, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]
print(samplepermanova)

#Light vs mild 
#subset samples
subset1<-subset_samples(oral.rel, t90x %in% c("Light", "Mild"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ t90x, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]
print(samplepermanova1)

#Light vs Mod/Severe 
#subset samples
subset1<-subset_samples(oral.rel, t90x %in% c("Light", "Moderate/Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ t90x, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]
print(samplepermanova1)

#Mild vs Mod/Severe
#subset samples
subset1<-subset_samples(oral.rel, t90x %in% c("Mild", "Moderate/Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ t90x, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]
print(samplepermanova1)

###############################################################################
##Figure 5.17
###############################################################################

#Comparison A 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, t90x  %in% c("Light","Mild"))
Comp1.OTU.Rel.Table = subset_samples(oral.rel, t90x  %in% c("Light","Mild"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~t90x)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$t90x <- droplevels(diagdds$t90x)

#Relevel Data
diagdds$t90x <- relevel(diagdds$t90x, ref ="Light")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]


#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#make column with names
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
df.1 = subset_samples(oral.rel, t90x  %in% c("Light"))
df.2 = subset_samples(oral.rel, t90x  %in% c("Mild"))

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

#compA - light v mild 
compA<-sig_results

#make a column for group 
compA$group <- "A" 

#Comparison B
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, t90x  %in% c("Light","Moderate/Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~t90x)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$t90x <- droplevels(diagdds$t90x)

#Relevel Data
diagdds$t90x <- relevel(diagdds$t90x, ref ="Light")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]


#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))


#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)
print(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#make column with names
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
df.1 = subset_samples(oral.rel, t90x  %in% c("Light"))
df.2 = subset_samples(oral.rel, t90x  %in% c("Moderate/Severe"))

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

#compB - light/mod/severe
compB<-sig_results

#make a column for group 
compB$group <- "B" 

#Comparison C
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, t90x  %in% c("Mild","Moderate/Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~t90x)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
#diagdds = estimateSizeFactors(diagdds)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$t90x <- droplevels(diagdds$t90x)

#Relevel Data
diagdds$t90x <- relevel(diagdds$t90x, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#make column with names
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
df.1 = subset_samples(oral.rel, t90x  %in% c("Mild"))
df.2 = subset_samples(oral.rel, t90x  %in% c("Moderate/Severe"))

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

#compC - mild/mod/severe
compC<-sig_results

#make a column for group 
compC$group <- "C" 

#combine df
df <- rbind(compA, compB, compC)

#reorder taxa 
df <- df %>%
  arrange(factor(group, levels = c("A", "B", "C")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - Hypoxia - whole group - NR.pdf",width=50, height=20) 
ggplot(df, aes(x = logFC, y = Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 15)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Light"="#bae1ff","Mild"="#00C9B7","Moderate/Severe"="#F7C7FF"), guide="none") +
  xlab("LogFC") +
  ylab("") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 60, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()


###############################################################################
#Figure 5.18A
###############################################################################
#subset samples 
oral<- subset(oral, OSA == "OSA")
oral.rel<-subset(oral.rel, OSA == "OSA")

#Calculates Shannon Diversity
sample_data(oral)$Shannon = diversity(otu_table(oral.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(oral))

#Make variable of interest a factor
Shannon$patient_type1 <- as.factor(Shannon$patient_type1)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#PLOT IT
pdf("Alpha Diversity - OSA Group by Covid Status - oral.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=patient_type1, y=Shannon, fill=patient_type1)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Post_COVID" = "#F0A4B0", "Incidental_COVID" = "#8ed4ac","Never_COVID" = "#F7D068")) + 
  ylab("Shannon Diversity") + 
  xlab("Patient_Group") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

#check stats - overall 
sampleshannon <- kruskal.test(Shannon ~ patient_type1, data = Shannon)
sampleshannon <- sampleshannon$p.value
print(sampleshannon)

#Check stats for post vs incidental 
shannon1 <- subset(Shannon, patient_type1 == "Post_COVID") 
shannon2 <- subset(Shannon, patient_type1 == "Incidental_COVID") 
shannon3<-rbind(shannon1, shannon2)
droplevels(shannon3$patient_type1)

sampleshannon2 <- kruskal.test(Shannon ~ patient_type1, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value
print(sampleshannon2)

#Check stats for post vs never
shannon4 <- subset(Shannon, patient_type1 == "Post_COVID") 
shannon5 <- subset(Shannon, patient_type1 == "Never_COVID") 
shannon6<-rbind(shannon4, shannon5)
droplevels(shannon6$patient_type1)

sampleshannon3 <- kruskal.test(Shannon ~ patient_type1, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value
print(sampleshannon3)

#Check stats for Incidental vs never 
shannon7 <- subset(Shannon, patient_type1 == "Incidental_COVID") 
shannon8 <- subset(Shannon, patient_type1 == "Never_COVID") 
shannon9<-rbind(shannon7, shannon8)

sampleshannon4 <- kruskal.test(Shannon ~ patient_type1, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value
print(sampleshannon4)

###############################################################################
#5.18B
###############################################################################

##Create Distance Matrix
vegdist   = vegdist(t(otu_table(oral.rel)), method="bray")

##Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)

##calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)

##Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

##Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~patient_type1,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="patient_type1",suffixes=c("",".centroid"))

#Plot it 
pdf("Beta Diversity - OSA Group - by pt type .pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, colour=patient_type1)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Post_COVID" = "#F0A4B0", "Incidental_COVID" = "#8ed4ac","Never_COVID" = "#F7D068")) + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=patient_type1), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=patient_type1))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=patient_type1), parse=TRUE,size=5, show.legend=FALSE) +
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


#Run the Statistics everyone 
samplepermanova <- adonis(vegdist ~ patient_type1, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]
print(samplepermanova)

#Check stats for Incidental vs never 
#subset samples
subset1<-subset_samples(oral.rel, patient_type1 %in% c("Incidental_COVID", "Never_COVID"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ patient_type1, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]
print(samplepermanova1)

#Check stats for post vs never 
#subset samples
subset1<-subset_samples(oral.rel, patient_type1 %in% c("Post_COVID", "Never_COVID"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ patient_type1, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]
print(samplepermanova1)

#Check stats for incidental vs post
#subset samples
subset1<-subset_samples(oral.rel, patient_type1 %in% c("Post_COVID", "Incidental_COVID"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ patient_type1, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]
print(samplepermanova1)

###############################################################################
#Figure 5.19
###############################################################################

#Comparison A
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, patient_type1  %in% c("Never_COVID","Incidental_COVID"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~patient_type1)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$patient_type1<- droplevels(diagdds$patient_type1)

#Relevel Data
diagdds$patient_type1 <- relevel(diagdds$patient_type1, ref ="Never_COVID")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))


#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)
print(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#make column with names
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
df.1 = subset_samples(oral.rel, patient_type1  %in% c("Never_COVID"))
df.2 = subset_samples(oral.rel, patient_type1  %in% c("Incidental_COVID"))

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
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, NA))


# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

#make group
compA<-sig_results

#make a column for group 
compA$group <- "A"

#Comparison B
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, patient_type1  %in% c("Post_COVID","Incidental_COVID"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~patient_type1)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$patient_type1<- droplevels(diagdds$patient_type1)

#Relevel Data
diagdds$patient_type1 <- relevel(diagdds$patient_type1, ref ="Incidental_COVID")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)


#taxa names 
res$gs <- ifelse(is.na(res$Genus),paste0(res$Family,"_",res$ASV), paste0(res$Genus,"_",res$ASV))
res$gs <- ifelse(is.na(res$Family), paste0(res$Order,"_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Order), paste0(res$Class, "_",res$ASV),res$gs)
res$gs <- ifelse(is.na(res$Class), paste0(res$Phylum,"_",res$ASV),res$gs)

#Make the full trail the First Column, make this new column the 1st column. 
res$names <- res$ASV
res$Gene.symbol <- res$gs

#Subset the different categories, subset by smokers / COPD 
df.1 = subset_samples(oral.rel, patient_type1  %in% c("Incidental_COVID"))
df.2 = subset_samples(oral.rel, patient_type1  %in% c("Post_COVID"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha,  res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

#compB 
compB<-sig_results

#make a column for group 
compB$group <- "B"

#Comparison C 
#Set alpha for differential Anlaysis, for adjusted p value 
alpha <- 0.2

#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, patient_type1  %in% c("Post_COVID","Never_COVID"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~patient_type1)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$patient_type1<- droplevels(diagdds$patient_type1)

#Relevel Data
diagdds$patient_type1 <- relevel(diagdds$patient_type1, ref ="Never_COVID")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)
print(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#make column with names
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
df.1 = subset_samples(oral.rel, patient_type1  %in% c("Never_COVID"))
df.2 = subset_samples(oral.rel, patient_type1  %in% c("Post_COVID"))

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
res$abundance <- ifelse(res$logFC>=1, res$abundance.2, 
                        ifelse(res$logFC<=-1, res$abundance.1, NA))


# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < 0.2)

#compc 
compC<-sig_results

#make a column for group 
compC$group <- "C"

#merge df
df <- rbind(compA, compB, compC)

# Arrange the data by group
df <- df %>%
  arrange(factor(group, levels = c("A", "B", "C")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot4 - OSA Group - COVID Status.pdf",width=40, height=15) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 15)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x") +  
  scale_fill_manual(values = c("Incidental_COVID"= "#8ed4ac", "Never_COVID" = "#F7D068", "Post_COVID" = "#F0A4B0", "Post_Covid" = "#F0A4B0", "notsig" = "white"), guide="none") +
  xlab("LogFC") +
  ylab("") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 45, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

###############################################################################
##Figure 5.20A 
##############################################################################
#subset samples 
oral<- subset(oral, patient_type1 == "Post_COVID")
oral.rel<-subset(oral.rel, OSA == "Post_COVID")


#Calculates Shannon Diversity
sample_data(oral)$Shannon = diversity(otu_table(oral.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(oral))

#Make variable of interest a factor
Shannon$OSA <- as.factor(Shannon$OSA)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))


#Plot it 
pdf("Alpha Diversity - Post-COVID Group - OSA Severity - oral.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=severity, y=Shannon, fill=severity)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75")) + 
  ylab("Shannon Diversity") + 
  xlab("Severity") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

#check stats
sampleshannon <- kruskal.test(Shannon ~ severity, data = Shannon)
sampleshannon <- sampleshannon$p.value


#Check stats for Normal vs Mild 
shannon1 <- subset(Shannon, severity == "Normal") 
shannon2 <- subset(Shannon, severity == "Mild") 
shannon3<-rbind(shannon1, shannon2)

sampleshannon2 <- kruskal.test(Shannon ~ severity, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value

#Check stats for Normal vs Moderate 
shannon4 <- subset(Shannon, severity == "Normal") 
shannon5 <- subset(Shannon, severity == "Moderate") 
shannon6<-rbind(shannon4, shannon5)

sampleshannon3 <- kruskal.test(Shannon ~ severity, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value

#Check stats for Normal vs Severe 
shannon7 <- subset(Shannon, severity == "Normal") 
shannon8 <- subset(Shannon, severity == "Severe") 
shannon9<-rbind(shannon7, shannon8)

sampleshannon4 <- kruskal.test(Shannon ~ severity, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value


#Check stats for Mild vs Mod
shannon10 <- subset(Shannon, severity == "Mild") 
shannon11 <- subset(Shannon, severity == "Moderate") 
shannon12<-rbind(shannon10, shannon11)

sampleshannon5 <- kruskal.test(Shannon ~ severity, data = shannon12)
sampleshannon5 <- sampleshannon5$p.value

#Check stats for Mild vs Severe
shannon13 <- subset(Shannon, severity == "Mild") 
shannon14 <- subset(Shannon, severity == "Severe") 
shannon15<-rbind(shannon13, shannon14)

sampleshannon6 <- kruskal.test(Shannon ~ severity, data = shannon15)
sampleshannon6 <- sampleshannon6$p.value


#Check stats for Moderate vs Severe
shannon16 <- subset(Shannon, severity == "Moderate") 
shannon17 <- subset(Shannon, severity == "Severe") 
shannon18<-rbind(shannon16, shannon17)

sampleshannon7 <- kruskal.test(Shannon ~ severity, data = shannon18)
sampleshannon7 <- sampleshannon7$p.value

###############################################################################
##Figure 5.20B
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
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~severity,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="severity",suffixes=c("",".centroid"))

pdf("Beta Diversity - Post-COVID group - OSA Severity.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=severity)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75"), guide = "none") + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=severity), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=severity))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=severity), parse=TRUE,size=5, show.legend=FALSE) +
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
samplepermanova <- adonis(vegdist ~ severity, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]

#Normal vs mild 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Mild"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Moderate
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Moderate 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Moderate vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Moderate", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

###############################################################################
##Figure 5.21
###############################################################################

#Comparison A
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Normal","Mild"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Mild"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Mild", "Normal")

#compA 
compA<-sig_results
compA$group <- "A"

#Comparason B 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Normal","Moderate"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Moderate"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results1 <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results1 = sig_results1[order(sig_results1$logFC, na.last = NA), ]
sig_results1$Gene.symbol <- factor(sig_results1$Gene.symbol, levels = sig_results1$Gene.symbol[order(sig_results1$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results1$enriched<-ifelse(sig_results1$logFC>0, "Moderate", "Normal")

#compB 
compB<-sig_results1
compB$group <- "B"

#Comparason C
Comp1.OTU.Rel.Table = subset_samples(oral.rel, severity  %in% c("Normal","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Severe", "Normal")

#compC 
compC<-sig_results

compC$group<-"C"

#combine variables 
df <- rbind(compA, compB, compC)

#reodrder 
df <- df %>%
  arrange(factor(group, levels = c("A", "B", "C")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - normal comps - xlims - Post-COVID Group.pdf",width=30, height=7) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 10)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75", "Severe OSA" = "#832c75"), guide="none") +
xlab("LogFC") +
  ylab("") + 
  xlim(-30, 25) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()


#Comparison D
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Mild","Moderate"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Moderate"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Moderate", "Mild")

#compD 
compD<-sig_results
compD$group <- "D"

#Comparason E 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Mild","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results1 <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results1 = sig_results1[order(sig_results1$logFC, na.last = NA), ]
sig_results1$Gene.symbol <- factor(sig_results1$Gene.symbol, levels = sig_results1$Gene.symbol[order(sig_results1$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results1$enriched<-ifelse(sig_results1$logFC>0, "Severe", "Mild")

#compE
compE<-sig_results1
compE$group <- "E"

#ComparasonF
Comp1.OTU.Rel.Table = subset_samples(oral.rel, severity  %in% c("Moderate","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Moderate")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Moderate"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Severe", "Moderate")

#compF
compF<-sig_results

compF$group<-"F"

#combine variables 
df <- rbind(compD, compE, compF)

#reodrder 
df <- df %>%
  arrange(factor(group, levels = c("D", "E", "F")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - severity comps - Post-COVID.pdf",width=30, height=7) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 10)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75", "Severe OSA" = "#832c75"), guide="none") +
xlab("LogFC") +
  ylab("") + 
  xlim(-30, 25) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

###############################################################################
##RDI Analysis 
##############################################################################
#load metadata 
sleepy.nr <- read.delim("")
rownames(sleepy.nr) <- sleep$sample.id

#make OSA variables using RDI 
#make a variable for OSA severity 
sleepy.nr$severity <-ifelse(sleepy.nr$rdi <=5, "Normal", 
                            ifelse(sleepy.nr$rdi <15, "Mild",
                                   ifelse(sleepy.nr$rdi <30, "Moderate", "Severe")))


#Make a variable for OSA 
sleepy.nr$OSA<-ifelse(sleepy.nr$rdi >=5, "OSA", "No OSA")
sleepy.nr$OSA<-as.factor(as.character(sleepy.nr$OSA))

#convert to metadata format
meta<-sample_data(sleepy.nr)

#make phyloseq
oral=phyloseq(otu, TAX, meta, phy)

#only visit 1
oral<-subset_samples(oral, visit == "1")

#just oral 
oral<-subset_samples(oral, Description == "Mouth.Wash")

#just POST-COVID
oral<-subset_samples(oral, patient_type1 == "Post-COVID")

#remove 0 abundance 
oral = subset_taxa(oral, rowSums(otu_table(oral)) != 0)

#make rel abundance table 
oral.rel <- transform_sample_counts(oral, function(x) x/sum(x))


###############################################################################
##Figure 5.22A
###############################################################################

#Calculates Shannon Diversity
sample_data(oral)$Shannon = diversity(otu_table(oral.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(oral))

#Make variable of interest a factor
Shannon$OSA <- as.factor(Shannon$OSA)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#Plot it 
pdf("Alpha Diversity - Post-COVID Group - OSA RDI Severity - oral.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=severity, y=Shannon, fill=severity)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75")) + 
  ylab("Shannon Diversity") + 
  xlab("Severity") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

#check stats
sampleshannon <- kruskal.test(Shannon ~ severity, data = Shannon)
sampleshannon <- sampleshannon$p.value


#Check stats for Normal vs Mild 
shannon1 <- subset(Shannon, severity == "Normal") 
shannon2 <- subset(Shannon, severity == "Mild") 
shannon3<-rbind(shannon1, shannon2)

sampleshannon2 <- kruskal.test(Shannon ~ severity, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value

#Check stats for Normal vs Moderate 
shannon4 <- subset(Shannon, severity == "Normal") 
shannon5 <- subset(Shannon, severity == "Moderate") 
shannon6<-rbind(shannon4, shannon5)

sampleshannon3 <- kruskal.test(Shannon ~ severity, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value

#Check stats for Normal vs Severe 
shannon7 <- subset(Shannon, severity == "Normal") 
shannon8 <- subset(Shannon, severity == "Severe") 
shannon9<-rbind(shannon7, shannon8)

sampleshannon4 <- kruskal.test(Shannon ~ severity, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value


#Check stats for Mild vs Mod
shannon10 <- subset(Shannon, severity == "Mild") 
shannon11 <- subset(Shannon, severity == "Moderate") 
shannon12<-rbind(shannon10, shannon11)

sampleshannon5 <- kruskal.test(Shannon ~ severity, data = shannon12)
sampleshannon5 <- sampleshannon5$p.value

#Check stats for Mild vs Severe
shannon13 <- subset(Shannon, severity == "Mild") 
shannon14 <- subset(Shannon, severity == "Severe") 
shannon15<-rbind(shannon13, shannon14)

sampleshannon6 <- kruskal.test(Shannon ~ severity, data = shannon15)
sampleshannon6 <- sampleshannon6$p.value


#Check stats for Moderate vs Severe
shannon16 <- subset(Shannon, severity == "Moderate") 
shannon17 <- subset(Shannon, severity == "Severe") 
shannon18<-rbind(shannon16, shannon17)

sampleshannon7 <- kruskal.test(Shannon ~ severity, data = shannon18)
sampleshannon7 <- sampleshannon7$p.value

###############################################################################
##Figure 5.22B
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
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~severity,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="severity",suffixes=c("",".centroid"))

#plot it
pdf("Beta Diversity - Post-COVID group - OSA RDI Severity.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=severity)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75"), guide = "none") + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=severity), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=severity))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=severity), parse=TRUE,size=5, show.legend=FALSE) +
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
samplepermanova <- adonis(vegdist ~ severity, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]

#Normal vs mild 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Mild"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Moderate
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Moderate 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Moderate vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Moderate", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

###############################################################################
##Figure 5.23
###############################################################################

#Comparison A
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Normal","Mild"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Mild"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Mild", "Normal")

#compA normal vs mild
compA<-sig_results
compA$group <- "A"

#Comparason B 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Normal","Moderate"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Moderate"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results1 <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results1 = sig_results1[order(sig_results1$logFC, na.last = NA), ]
sig_results1$Gene.symbol <- factor(sig_results1$Gene.symbol, levels = sig_results1$Gene.symbol[order(sig_results1$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results1$enriched<-ifelse(sig_results1$logFC>0, "Moderate", "Normal")

#compB norm vs moderate 
compB<-sig_results1
compB$group <- "B"

#combine variables 
df <- rbind(compA, compB)

#reodrder 
df <- df %>%
  arrange(factor(group, levels = c("A", "B")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - normal comps - Post-COVID - RDI.pdf",width=30, height=7) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 10)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75", "Severe OSA" = "#832c75"), guide="none") +
xlab("LogFC") +
  ylab("") + 
  xlim(-30, 25) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

###############################################################################
##Figure 5.24
###############################################################################

#Comparison A
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Mild","Moderate"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Moderate"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Moderate", "Mild")

#compA 
compA<-sig_results
compA$group <- "A"

#Comparason B 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Mild","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results1 <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results1 = sig_results1[order(sig_results1$logFC, na.last = NA), ]
sig_results1$Gene.symbol <- factor(sig_results1$Gene.symbol, levels = sig_results1$Gene.symbol[order(sig_results1$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results1$enriched<-ifelse(sig_results1$logFC>0, "Severe", "Mild")

#compB 
compB<-sig_results1
compB$group <- "B"

#Comparason C
Comp1.OTU.Rel.Table = subset_samples(oral.rel, severity  %in% c("Moderate","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Moderate")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Moderate"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Severe", "Moderate")

#compC 
compC<-sig_results

compC$group<-"C"

#combine variables 
df <- rbind(compA, compB, compC)

#reodrder 
df <- df %>%
  arrange(factor(group, levels = c("A", "B", "C")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - normal comps - Post-COVID - RDI.pdf",width=30, height=7) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 10)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75", "Severe OSA" = "#832c75"), guide="none") +
xlab("LogFC") +
  ylab("") + 
  xlim(-30, 25) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

###############################################################################
##ODI Analysis 
##############################################################################
#load metadata 
sleepy.nr <- read.delim("")
rownames(sleepy.nr) <- sleep$sample.id

#make OSA variables using ODI
#make a variable for OSA severity 
sleepy.nr$severity <-ifelse(sleepy.nr$odi <=5, "Normal", 
                            ifelse(sleepy.nr$odi <15, "Mild",
                                   ifelse(sleepy.nr$odi <30, "Moderate", "Severe")))


#Make a variable for OSA 
sleepy.nr$OSA<-ifelse(sleepy.nr$odi >=5, "OSA", "No OSA")
sleepy.nr$OSA<-as.factor(as.character(sleepy.nr$OSA))

#convert to metadata format
meta<-sample_data(sleepy.nr)

#make phyloseq
oral=phyloseq(otu, TAX, meta)

#only visit 1
oral<-subset_samples(oral, visit == "1")

#just oral 
oral<-subset_samples(oral, Description == "Mouth.Wash")

#just POST-COVID
oral<-subset_samples(oral, patient_type1 == "Post-COVID")

#remove 0 abundance 
oral = subset_taxa(oral, rowSums(otu_table(oral)) != 0)


#make rel abundance table 
oral.rel <- transform_sample_counts(oral, function(x) x/sum(x))


###############################################################################
##Figure 5.25A
###############################################################################

#Calculates Shannon Diversity
sample_data(oral)$Shannon = diversity(otu_table(oral.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(oral))

#Make variable of interest a factor
Shannon$OSA <- as.factor(Shannon$OSA)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#plot it 
pdf("Alpha Diversity - Post-COVID Group - OSA ODI Severity - oral.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=severity, y=Shannon, fill=severity)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75")) + 
  ylab("Shannon Diversity") + 
  xlab("Severity") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

#check stats
sampleshannon <- kruskal.test(Shannon ~ severity, data = Shannon)
sampleshannon <- sampleshannon$p.value


#Check stats for Normal vs Mild 
shannon1 <- subset(Shannon, severity == "Normal") 
shannon2 <- subset(Shannon, severity == "Mild") 
shannon3<-rbind(shannon1, shannon2)

sampleshannon2 <- kruskal.test(Shannon ~ severity, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value

#Check stats for Normal vs Moderate 
shannon4 <- subset(Shannon, severity == "Normal") 
shannon5 <- subset(Shannon, severity == "Moderate") 
shannon6<-rbind(shannon4, shannon5)

sampleshannon3 <- kruskal.test(Shannon ~ severity, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value

#Check stats for Normal vs Severe 
shannon7 <- subset(Shannon, severity == "Normal") 
shannon8 <- subset(Shannon, severity == "Severe") 
shannon9<-rbind(shannon7, shannon8)

sampleshannon4 <- kruskal.test(Shannon ~ severity, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value


#Check stats for Mild vs Mod
shannon10 <- subset(Shannon, severity == "Mild") 
shannon11 <- subset(Shannon, severity == "Moderate") 
shannon12<-rbind(shannon10, shannon11)

sampleshannon5 <- kruskal.test(Shannon ~ severity, data = shannon12)
sampleshannon5 <- sampleshannon5$p.value

#Check stats for Mild vs Severe
shannon13 <- subset(Shannon, severity == "Mild") 
shannon14 <- subset(Shannon, severity == "Severe") 
shannon15<-rbind(shannon13, shannon14)

sampleshannon6 <- kruskal.test(Shannon ~ severity, data = shannon15)
sampleshannon6 <- sampleshannon6$p.value


#Check stats for Moderate vs Severe
shannon16 <- subset(Shannon, severity == "Moderate") 
shannon17 <- subset(Shannon, severity == "Severe") 
shannon18<-rbind(shannon16, shannon17)

sampleshannon7 <- kruskal.test(Shannon ~ severity, data = shannon18)
sampleshannon7 <- sampleshannon7$p.value


###############################################################################
##Figure 5.25B
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
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~severity,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="severity",suffixes=c("",".centroid"))

#plot it 
pdf("Beta Diversity - Post-COVID group - OSA ODI Severity.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=severity)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75"), guide = "none") + 
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=severity), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=severity))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=severity), parse=TRUE,size=5, show.legend=FALSE) +
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
samplepermanova <- adonis(vegdist ~ severity, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]

#Normal vs mild 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Mild"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Moderate
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Normal vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Normal", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Moderate 
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Moderate"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Mild vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Mild", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

#Moderate vs Severe
#subset samples
subset1<-subset_samples(oral.rel, severity %in% c("Moderate", "Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ severity, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]

###############################################################################
##Figure 5.26
###############################################################################

#Comparison A
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Normal","Mild"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Mild"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Mild", "Normal")

#compA 
compA<-sig_results
compA$group <- "A"

#Comparason B 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Normal","Moderate"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Moderate"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results1 <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results1 = sig_results1[order(sig_results1$logFC, na.last = NA), ]
sig_results1$Gene.symbol <- factor(sig_results1$Gene.symbol, levels = sig_results1$Gene.symbol[order(sig_results1$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results1$enriched<-ifelse(sig_results1$logFC>0, "Moderate", "Normal")

#compB 
compB<-sig_results1
compB$group <- "B"

#Comparason C
Comp1.OTU.Rel.Table = subset_samples(oral.rel, severity  %in% c("Normal","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Normal")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Normal"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Severe", "Normal")

#compC 
compC<-sig_results

compC$group<-"C"

#combine variables 
df <- rbind(compA, compB, compC)

#reorder
df <- df %>%
  arrange(factor(group, levels = c("A", "B", "C")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - normal comps - Post-COVID - ODI.pdf",width=30, height=7) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 10)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75", "Severe OSA" = "#832c75"), guide="none") +
xlab("LogFC") +
  ylab("") + 
  xlim(-30, 25) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

#Comparison D
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Mild","Moderate"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)


#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Moderate"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results = sig_results[order(sig_results$logFC, na.last = NA), ]
sig_results$Gene.symbol <- factor(sig_results$Gene.symbol, levels = sig_results$Gene.symbol[order(sig_results$logFC)])

#make a variable for enrichment for colouring the chart 
sig_results$enriched<-ifelse(sig_results$logFC>0, "Moderate", "Mild")

#compD 
compD<-sig_results
compD$group <- "D"

#Comparason E 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, severity  %in% c("Mild","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))
#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Mild"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
res$abundance <- ifelse(res$logFC>=1 & res$adj.P.Val < alpha,  res$abundance.2, 
                        ifelse(res$logFC<=-1 & res$adj.P.Val < alpha, res$abundance.1, 0))

# Subset significant results (adjust p-value threshold as needed)
sig_results1 <- subset(res, adj.P.Val < alpha)

# Reorder Results based on logFC
sig_results1 = sig_results1[order(sig_results1$logFC, na.last = NA), ]
sig_results1$Gene.symbol <- factor(sig_results1$Gene.symbol, levels = sig_results1$Gene.symbol[order(sig_results1$logFC)])


#make a variable for enrichment for colouring the chart 
sig_results1$enriched<-ifelse(sig_results1$logFC>0, "Severe", "Mild")

#compE 
compE<-sig_results1
compE$group <- "E"

#Comparason F
Comp1.OTU.Rel.Table = subset_samples(oral.rel, severity  %in% c("Moderate","Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~severity)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$severity <- droplevels(diagdds$severity)

#Relevel Data
diagdds$severity <- relevel(diagdds$severity, ref ="Moderate")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]



#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#Replace any no genus annotation as NA so we can get rid of them later
res[res=="Bacteria_unclassified"]<-NA

#make column with names
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
df.1 = subset_samples(oral.rel, severity  %in% c("Moderate"))
df.2 = subset_samples(oral.rel, severity  %in% c("Severe"))

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
sig_results$enriched<-ifelse(sig_results$logFC>=0, "Severe", "Moderate")

#compF
comf<-sig_results

compF$group<-"F"

#combine variables 
df <- rbind(compD, compE, compF)

#reorder
df <- df %>%
  arrange(factor(group, levels = c("D", "E", "F")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - SEVERITY comps - Post-COVID - 0DI.pdf",width=30, height=7) 
ggplot(df, aes(x = logFC, y =Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 10)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Normal"="#63AAC0","Mild"="#28602b","Moderate"="#F99B45", "Severe" = "#832c75", "Severe OSA" = "#832c75"), guide="none") +
xlab("LogFC") +
  ylab("") + 
  xlim(-30, 25) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()

###############################################################################
##Figure 5.27A
###############################################################################

#Calculates Shannon Diversity
sample_data(oral.rel)$Shannon = diversity(otu_table(oral.rel), index = "shannon", MARGIN = 2, base = exp(1))

#Convert to data frame for ggplot
Shannon = data.frame(sample_data(oral.rel))

#Make variable of interest a factor
Shannon$t90x <- as.factor(Shannon$t90x)

#Make Sure Shannon is Numeric
Shannon$Shannon <- as.numeric(as.character(Shannon$Shannon))

#make figure 
pdf("Alpha Diversity - Post-COVID - Hypoxia - MW.pdf", height = 10, width = 15, useDingbats=FALSE)
ggplot(Shannon, aes(x=t90x, y=Shannon, fill=t90x)) + 
  stat_boxplot(geom ='errorbar', width=0.1)+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(shape=1, position=position_jitter(0.2))+
  scale_fill_manual(values=c("Light"="#bae1ff","Mild"="#00C9B7","Moderate/Severe"="#F7C7FF"))+
  ylab("Shannon Diversity") + 
  xlab("Hypoxia") +
  #ggtitle("Alpha Diversity by COPD Status") +
  theme +
  guides(fill = FALSE) 
dev.off()

#check stats - overall 
sampleshannon <- kruskal.test(Shannon ~ t90x, data = Shannon)
sampleshannon <- sampleshannon$p.value

#Check stats for Light vs Mild 
shannon1 <- subset(Shannon, t90x == "Light") 
shannon2 <- subset(Shannon, t90x == "Mild") 
shannon3<-rbind(shannon1, shannon2)

sampleshannon2 <- kruskal.test(Shannon ~ t90x, data = shannon3)
sampleshannon2 <- sampleshannon2$p.value

#Check stats for Light vs Mod/Severe
shannon4 <- subset(Shannon, t90x == "Light") 
shannon5 <- subset(Shannon, t90x == "Moderate/Severe") 
shannon6<-rbind(shannon4, shannon5)

sampleshannon3 <- kruskal.test(Shannon ~ t90x, data = shannon6)
sampleshannon3 <- sampleshannon3$p.value

#Check stats for Mild vs Mod/Severe
shannon7 <- subset(Shannon, t90x == "Mild") 
shannon8 <- subset(Shannon, t90x == "Moderate/Severe") 
shannon9<-rbind(shannon7, shannon8)

sampleshannon4 <- kruskal.test(Shannon ~ t90x, data = shannon9)
sampleshannon4 <- sampleshannon4$p.value


###############################################################################
##Figure 5.27B
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
newResults <- merge(x = CmdScale, y = sample_data(oral.rel), by = "row.names", all.x = TRUE)

##Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

##Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~t90x,data= newResults, mean)

##Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="t90x",suffixes=c("",".centroid"))

#plot it 
pdf("Beta Diversity - Post-COVID group - Hypoxia - MW.pdf", height = 10, width = 10, useDingbats=FALSE)
ggplot(newResults, aes(PC1, PC2, color=t90x)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("Light"="#bae1ff","Mild"="#00C9B7","Moderate/Severe"="#F7C7FF"), guide="none")+
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=t90x), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=t90x))+ 
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=t90x), parse=TRUE,size=5, show.legend=FALSE) +
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
samplepermanova <- adonis(vegdist ~ t90x, data.adonis)
samplepermanova <- as.data.frame(samplepermanova$aov.tab)
samplepermanova <- samplepermanova$'Pr(>F)'[1]
print(samplepermanova)

#Light vs mild 
#subset samples
subset1<-subset_samples(oral.rel, t90x %in% c("Light", "Mild"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ t90x, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]
print(samplepermanova1)

#Light vs Mod/Severe 
#subset samples
subset1<-subset_samples(oral.rel, t90x %in% c("Light", "Moderate/Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ t90x, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]
print(samplepermanova1)

#Mild vs Mod/Severe
#subset samples
subset1<-subset_samples(oral.rel, t90x %in% c("Mild", "Moderate/Severe"))

#calculate distance matrix
vegdist1   = vegdist(t(otu_table(subset1)), method="bray")

#Create Table for Statistics    
data.adonis1 <- data.frame(sample_data(subset1))

#Run the Statistics
samplepermanova1 <- adonis(vegdist1 ~ t90x, data.adonis1)
samplepermanova1 <- as.data.frame(samplepermanova1$aov.tab)
samplepermanova1 <- samplepermanova1$'Pr(>F)'[1]
print(samplepermanova1)

###############################################################################
##Figure 5.28
###############################################################################

#Comparison A 
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, t90x  %in% c("Light","Mild"))
Comp1.OTU.Rel.Table = subset_samples(oral.rel, t90x  %in% c("Light","Mild"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~t90x)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Dispersion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$t90x <- droplevels(diagdds$t90x)

#Relevel Data
diagdds$t90x <- relevel(diagdds$t90x, ref ="Light")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]


#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#make column with names
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
df.1 = subset_samples(oral.rel, t90x  %in% c("Light"))
df.2 = subset_samples(oral.rel, t90x  %in% c("Mild"))

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

#compA - light v mild 
compA<-sig_results

#make a column for group 
compA$group <- "A" 

#Comparison B
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, t90x  %in% c("Light","Moderate/Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~t90x)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$t90x <- droplevels(diagdds$t90x)

#Relevel Data
diagdds$t90x <- relevel(diagdds$t90x, ref ="Light")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]


#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))


#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)
print(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#make column with names
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
df.1 = subset_samples(oral.rel, t90x  %in% c("Light"))
df.2 = subset_samples(oral.rel, t90x  %in% c("Moderate/Severe"))

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

#compB - light/mod/severe
compB<-sig_results

#make a column for group 
compB$group <- "B" 

#Comparison C
#Subset only variables for comparisons
Comp1.OTU.Table = subset_samples(oral, t90x  %in% c("Mild","Moderate/Severe"))

#Convert To DESEQ
diagdds <- phyloseq_to_deseq2(Comp1.OTU.Table, ~t90x)

#Calculate geometric means prior to estimate size factor, estimate means
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
#diagdds = estimateSizeFactors(diagdds)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Make sure all unwanted levels are removed from dataset, removes any 'hidden' non.smokers 
diagdds$t90x <- droplevels(diagdds$t90x)

#Relevel Data
diagdds$t90x <- relevel(diagdds$t90x, ref ="Mild")

#Run the differential Analysis
diagdds<- DESeq(diagdds)

#output the table of differential analysis
res <- results(diagdds)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]


#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(oral)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus)

#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#convert to character, changes the row names to column ASV
res$ASV <- rownames(res)

#make column with names
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
df.1 = subset_samples(oral.rel, t90x  %in% c("Mild"))
df.2 = subset_samples(oral.rel, t90x  %in% c("Moderate/Severe"))

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

#compC - mild/mod/severe
compC<-sig_results

#make a column for group 
compC$group <- "C" 

#combine df
df <- rbind(compA, compB, compC)

#reorder taxa 
df <- df %>%
  arrange(factor(group, levels = c("A", "B", "C")), Gene.symbol) %>%  
  mutate(Gene.symbol = factor(Gene.symbol, levels = rev(unique(Gene.symbol))))  

#plot it 
pdf("Combined bubble plot - Hypoxia - Post-COVID group - NR.pdf",width=50, height=20) 
ggplot(df, aes(x = logFC, y = Gene.symbol)) +
  geom_point(aes(fill=enriched, size=abundance), color="black", shape = 21)  +
  scale_size_continuous(name="Relative\nAbundance", range = c(5, 15)) +
  geom_segment( aes(yend=Gene.symbol, xend=0)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
  facet_wrap(~ group, scales = "free_x", ncol=5) +  
  scale_fill_manual(values = c("Light"="#bae1ff","Mild"="#00C9B7","Moderate/Severe"="#F7C7FF"), guide="none") +
  xlab("LogFC") +
  ylab("") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(color = "grey85", linetype = "dotted", size = 0.5),  
        panel.grid.minor = element_line(color = "grey90", linetype = "dotted", size = 0.3), 
        strip.background=element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 60, colour="black"),
        axis.ticks = element_line(colour = "black", size = 2),  
        axis.ticks.length = unit(0.5, "cm"), 
        panel.spacing = unit(2, "lines"))
dev.off()
