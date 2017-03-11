#####PART 1: DATA MANIPULATION ON BASH
all python scripts used in this project can be found on this repository


###Running DIAMOND against all of UniProt(Swiss-Prot and TrEMBL) locally
wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
gunzip uniprot_*

#merging Swiss-Prot and TrEMBL
cat uniprot_sprot.fasta uniprot_trembl.fasta > uniprot_all.fasta

#preparing UniProt database in DIAMOND format
diamond makedb --in uniprot_all.fasta -d uniprot

#Performing sequence alignment of all samples
ls ../../DiamondUniprot/*fna | while read line
	do sample=`echo $line | cut -f 4 -d "/"`; diamond blastx -d uniprot -q $line -o ${sample}.m8 --sensitive --top 1
	done

#creating a database for the file that maps taxon IDs to accession numbers for fast parsing:
sqlite3 acces2taxidGB.sqlite
.separator "\t"
CREATE TABLE tab (unsused1, accnr, taxid, unused2);
.import nucl_gb.accession2taxid tab
CREATE INDEX accnr ON tab (accnr);

#extracting fungal reads from DIAMOND results:
ls *fna.m8 | cut -f 1 -d "_" | while read line
	do ./allUniprotIfFungi.py speclistUni.sqlite taxonomy-ancestor%3A4751.list ${line}_merd_hiQual.fna.m8 > ../../FungalReads/fungalFasta3/${line}.txt
	done

#extracting the sequences of the fungal reads
ls | while read line
	do sample=`echo $line | cut -f 1 -d "."`; extractFastaSeq.py ../../DiamondUniprot/${sample}*.fna ${sample}.txt > ${sample}.fasta
	done
cat *.fasta > allFungi.fasta

#compressing for uploading to Kaiju
gzip < allFungi.fasta > allFungi.fasta.gz

#translating for uploading to KEGG's GhostKOALA
translator.py allFungi.fasta > uniReads.faa

##Processing KEGG's and Kaiju's output

#extracting the fungal reads from Kaiju output
./kaijushorter.py taxonomy-ancestor%3A4751.list kaiju.out > kaijuReads.txt

#getting the intersect of the reads where both KEGG and Kaiju determine to be fungal
#keggTranslatedCleanup.py gets only the ORF from a read with the highest score and keeps it only if it is fungal
grep -F -f <(keggTranslatedCleanup.py keggReads.txt | cut -f 1 -d "_") kaijuReads.txt > keggKaijuReads.txt


###Performing sequence alignment of samples to NCBI/Nt on Milou cluster at UPPMAX
paralleliser.sh


#contents of paralleliser.sh
#!/bin/bash -l

ls /proj/b2016289/Analyses/DIAMOND/2.trimmed/*.fna | while read myline
	do myfile=`echo $myline | awk ' BEGIN {FS="/"} {print $NF} '`; export myfile; export myline; sbatch /proj/b2016289/Scripts/blastNt.sh
	done

	
#contents of blastNt.sh
#!/bin/bash -l

#SBATCH -A b2016289
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 7-00:00:00
#SBATCH -J blastNt

module load bioinfo-tools
module load blast/2.4.0+

blastn -query $myline -db /sw/data/uppnex/blast_databases/nt -out $myfile -outfmt 6 -max_target_seqs 1

#creating a database for the file that maps taxon IDs to accession numbers for fast parsing:
sqlite3 speclistUni.sqlite
.separator "\t"
CREATE TABLE tab (uniId, domain, taxid);
.import ready4DB.txt tab
CREATE INDEX uniId ON tab (uniId);

#extracting fungal reads from blast results
ls *blastn | while read line
	do sample=`echo $line | cut -f 1 -d "_"`; ./allNCBIFungi.py fungiTaxonomy.txt acces2taxidGB.sqlite ${sample}_merd_hiQual.fna.blastn > ../FungalReads/fungalFasta2/${line}
	done
cd ../FungalReads/fungalFasta2/

#extracting the sequences of the fungal reads
ls | while read line
	do sample=`echo $line | cut -f 1 -d "_"`; extractFastaSeq.py ../../DiamondUniprot/${sample}*.fna ${sample}* > ${sample}.fasta
	done
cat *.fasta > allFungi.fasta

##Processing KEGG's and Kaiju's output

#extracting the fungal reads from Kaiju output
kaijushorter.py fungiTaxonomy.txt kaiju.out > kaijuFungi.txt

#getting the intersect of the reads where both KEGG and Kaiju determine to be fungal
grep -F -f <(./keggTranslatedCleanup.py kegg.readsNt) kaijuFungi.txt > KaijuKeggReads.txt

###Joining reads from Uniprot and NCBI/Nt
cat keggKaijuReads.txt KaijuKeggReads.txt | sort | uniq > kegkaiuniNTall.txt

#extracting sequences of fungal reads
extractFastaSeq.py ../../DiamondUniprot/everything.fna kegkaiuniNTall.txt > filteredFungi.fasta

#Assembling the reads
megahit -r filteredFungi.fasta -o assembly --k-min 27 --k-max 97 --k-step 2 2> thelog.log

#translating the contigs for KEGG
translator.py assembly/final.contigs.fa > transcontigs.faa

#getting the contigs which have a KEGG Orthology Identifier (from keggAnn.txt file)
grep -wF -f <(keggTranslatedCleanup.py user.out.top) keggAnn.txt | awk '$2 != ""{print}' > annContigs.txt

'''
#getting the KEGG Orthology group from the KEGG Orthology identifier for each contig
keggortho2whatever.py annContigs.txt ko > matrixKo.txt
'''

#getting the Pathway from the KEGG Orthology identifier for each contig
keggortho2whatever.py annContigs.txt path > matrixPath.txt

###mapping the reads to the conrigs
bwa index final.contigs.fa
bwa mem assembly/final.contigs.fa filteredFungi.fasta > reads.sam

#usefull command to get the percentage of reads that map back
echo "scale=4;`grep -v ^\@ reads.sam | awk '$3 != "*"{print $1}' | sort | uniq | wc -l`/`grep -v ^\@ reads.sam | awk '{print $1}' | sort | uniq | wc -l`" | bc

#splitting the fungal reads per sample
ls ../../DiamondUniprot/*hiQual.fna | while read line
	do sample=`echo $line | cut -f 4 -d "/" | cut -f 1 -d "_"`; grep -wF -f kegkaiuniNTall.txt $line | cut -f 2 -d ">" | awk '{print $1}' > ${sample}.txt
	done

#creating a table of read counts per sample
maketablenew.py reads.sam G* W* > ../../AnalysesR/theTable.txt

#based on theTable.txt file, a table of read counts per pathway was made
./makeextendedtables.py matrixPath.txt ../../AnalysesR/theTable.txt path > ../../AnalysesR/pathcount.txt

###scanning pfam

#extracting the contigs and frames of contigs for which we had annotation from KEGG
cat annContigs.txt | awk '{print $1}' > correctFrames.txt
extractFastaSeq.py transcontigs.faa correctFrames.txt > correctFrames.faa

#formating the file for pfam
cat correctFrames.faa | tr "*" "X" > pfamReady.faa

#performing the pfam scan
./pfam_scan.pl -fasta pfamReady.faa -dir theDir/ -outfile newRun

#extracting contig, frame and pfam protein family info from pfam results
awk 'substr($1,1,1) != "#" && substr($1,1,1) != ""{print}' newRun | awk '{split($1, a, "_"); print a[1]"_"a[2]"\t"$7}' > going2R.txt

#making the table with read counts for each pfam protein family
makeextendedtables.py going2R.txt ../../../AnalysesR/theTable.txt ko > ../../../AnalysesR/pfamatrix.



#####PART 2: DATA MANIPULATION ON R

#setting the environment
library(DESeq2)
library(ggplot2)
library(psych)
library(genefilter)
library(reshape2)
library(scales)
library("Rgraphviz")
library(DiagrammeR)
library(gridExtra)
library(VennDiagram)

#plotting the fungal reads per sample (number of reads taken from bash)
fungalreads <- c(344, 309, 1179, 431, 475, 1936, 816, 425, 444, 1018)
fungaldataframe <- as.data.frame(fungalreads, c("G1", "G2", "G3", "G4", "G5", "W1", "W2", "W3", "W4", "W5"))
type <- as.factor(c(rep("Grass",5),rep("Wheat",5)))
fungaldataframe$type <- type
fungplot <- ggplot(fungaldataframe, aes(x=row.names(fungaldataframe), y=fungalreads, fill=type)) + 
			geom_bar(stat = "identity") + theme_bw(22) + scale_fill_brewer(palette="Dark2") + scale_y_continuous(labels = comma) + 
			labs(x="Sample", y="Read count", subtitle="Fungal reads per sample", title="b") + theme(legend.position = "none")

#plotting the total number of reads per sample (number of reads taken from bash)
totalreadspersample <- c(1517382, 1452487, 1593155, 1791749, 819006, 981838, 1670596, 371522, 446138, 759296)
totaldatfram <- as.data.frame(totalreadspersample, c("G1", "G2", "G3", "G4", "G5", "W1", "W2", "W3", "W4", "W5"))
totaldatfram$type <- type
totplot <- ggplot(totaldatfram, aes(x=row.names(totaldatfram), y=totalreadspersample, fill=type)) + geom_bar(stat = "identity") + 
			theme_bw(22) + scale_fill_brewer(palette="Dark2") + scale_y_continuous(labels = comma) + 
			labs(x="Sample", y="Read count", subtitle="Total reads per sample", title="a") + theme(legend.position = "none")

			
			
#creating venn diagrams of KEGG and Kaiju fungal reads
NTkegg <- read.table("keggReads.txt")
NTkaiju <- read.table("kaijuFungi.txt")
lelist <- c(NTkegg, NTkaiju)
names(lelist) <- c("KEGG", "Kaiju")

unikegg <- read.table("keggOnlyreads.txt")
unikaiju <- read.table("kaijuReads.txt")
lalist <- c(unikegg, unikaiju)
names(lalist) <- c("KEGG", "Kaiju")

Enti <- read.table("KaijuKeggReads.txt")
uniprall <- read.table("keggKaijuReads.txt")
lilist <- c(uniprall, Enti)
names(lilist) <- c("UniProt", "Nucleotide")

venn.plot1 <- venn.diagram(lelist, NULL, fill=c("darkmagenta", "darkblue"), alpha=c(0.5,0.5), sub = "Nucleotide", sub.cex = 2, main.pos = c(0,1), main = "b", main.cex = 2.5, cex = 2, cat.cex = 2)
venn.plot2 <- venn.diagram(lalist, NULL, fill=c("darkmagenta", "darkblue"), alpha=c(0.5,0.5), sub = "UniProt", sub.cex = 2, main.pos = c(0,1), main = "a", main.cex = 2.5, cex = 2, cat.cex = 2)
venn.plot3 <- venn.diagram(lilist, NULL, fill=c("darkmagenta", "darkblue"), alpha=c(0.5,0.5), main = "c", main.pos = c(0,1), main.cex = 2.5, cex = 2, cat.cex = 2)

lay <- rbind(c(1,1,2,2),c(1,1,2,2),c(NA, 3, 3,NA),c(NA,3,3, NA))
grid.arrange(gTree(children=venn.plot2), gTree(children=venn.plot1), gTree(children=venn.plot3), layout_matrix = lay)



#loading and manipulating the pathway - read count table
pathmatrix <- read.table("pathcount.txt", sep = "\t")
pathmatrix <- pathmatrix[,-1]
pathmatrixcast <- dcast(pathmatrix, V2 ~ V4, fun.aggregate = sum, value.var = 'V3')

#loading the data to DESeq2 package enviroment to perform transformations
fordeseq <- t(pathmatrixcast[,-1])
colnames(fordeseq) <- pathmatrixcast$V2
coldata <- data.frame(row.names=pathmatrixcast[,1], type=as.factor(c(rep("Grass",5),rep("Wheat",5))))
dds <- DESeqDataSetFromMatrix(countData = fordeseq, colData = coldata, design = ~ type)
rld <- rlog(dds)

#making a custom deseq plotPCA function to extract PC3
plotPCA2.DESeqTransform = function(object, intgroup="condition", ntop=500, returnData=FALSE)
{
    # calculate the variance for each gene
    rv <- rowVars(assay(object))
    
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(object)[select,]))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
    
    # add the intgroup factors together to create a new grouping factor
    group <- if (length(intgroup) > 1) {
        factor(apply( intgroup.df, 1, paste, collapse=" : "))
    } else {
        colData(object)[[intgroup]]
    }
    
    # assembly the data for the plot
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], group=group, intgroup.df, name=colnames(object), percentVar)
    
    if (returnData) {
        #attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }
}

#drawing PCA plots
myfuncttest <- plotPCA2.DESeqTransform(rld, intgroup = "type", returnData = TRUE)
p1a <- ggplot(myfuncttest, aes(PC1, PC2, colour=type)) + geom_text(aes(label=colnames(fordeseq)),hjust=0, vjust=0, size = 6) + 
		theme_bw(22) + theme(legend.position="none") + scale_color_brewer(palette="Dark2") + labs(subtitle = "a")
p3a <- ggplot(myfuncttest, aes(PC2, PC3, colour=type)) + geom_text(aes(label=colnames(fordeseq)),hjust=0, vjust=0, size = 6) + 
		theme_bw(22) + theme(legend.position="none") + scale_color_brewer(palette="Dark2")+ labs(subtitle = "b")
grid.arrange(p1a, p3a, ncol = 2)



##reducing pathways
#filtering pathways by number of samples with no reads
rm(columnsmorethanthree)
counter1 <- 0
columnsums <- colSums(as.matrix(pathmatrixcast[,-1]) != 0)
for(i in columnsums){
    counter1 <- counter1 + 1
    if(i > 3){
        if(exists("columnsmorethanthree")){
            columnsmorethanthree <- c(columnsmorethanthree, columnsums[counter1])
        } else{
            columnsmorethanthree <- columnsums[counter1]
        }
    }
}
rm(pathfiltered)
for(i in names(columnsmorethanthree)){
    if (exists("pathfiltered")){
        pathfiltered <- c(pathfiltered, pathmatrixcast[i])
    } else {
        pathfiltered <- pathmatrixcast[i]
        }
}
pathfiltered <- as.data.frame(pathfiltered)
names(columnsmorethanthree)

#filtering pathways by variance (genefilter package)
forgenefilterpath <- new("ExpressionSet", exprs=as.matrix(t(pathmatrixcast[,-1])))
a1 <- varFilter(forgenefilterpath)
varfiltpath <- t(exprs(a1))

#further filtering variance filtered pathway set by excluding the ones with more than 3 samples with more reads
rm(columnsmorethanthree)
counter1 <- 0
columnsums <- colSums(as.matrix(varfiltpath) != 0)
for(i in columnsums){
    counter1 <- counter1 + 1
    if(i > 3){
        if(exists("columnsmorethanthree")){
            columnsmorethanthree <- c(columnsmorethanthree, columnsums[counter1])
        } else{
            columnsmorethanthree <- columnsums[counter1]
        }
    }
}
rm(pathfiltered)
for(i in names(columnsmorethanthree)){
    if (exists("pathfiltered")){
        pathfiltered <- c(pathfiltered, pathmatrixcast[i])
    } else {
        pathfiltered <- pathmatrixcast[i]
        }
}
pathfiltered <- as.data.frame(pathfiltered)
names(columnsmorethanthree)

#the two ways of filtering lead to the same subset
#performing wilcoxon Mann-Whitney statistical test
wilcoxonMannWhitney <- apply(pathfiltered, 2, function(x) wilcox.test(x ~ stype))
rm(pvalue)
for(i in wilcoxonMannWhitney){
    if(exists("pvalue")){
        pvalue <- c(pvalue, i$p.value)
    } else {
        pvalue <- i$p.value
    }
}

#correcting for multiple testing using Benjamini-Hochberg method
p.adjust(pvalue, "BH")



##pathway loading handling
#making a function that extracts the transformed loadings from DESeq2 PCA
loadingsfromRlog1.DESeqTransform = function(object, intgroup="condition", ntop=500, returnData=FALSE)
{
    # calculate the variance for each gene
    rv <- rowVars(assay(object))
    
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(object)[select,]))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
    
    # add the intgroup factors together to create a new grouping factor
    group <- if (length(intgroup) > 1) {
        factor(apply( intgroup.df, 1, paste, collapse=" : "))
    } else {
        colData(object)[[intgroup]]
    }
    
    # assembly the data for the plot
    d1 <- data.frame(pca$rotation[,1])
    d2 <- data.frame(pca$rotation[,2])
    colnames(d1) <- "PC1"
    colnames(d2) <- "PC2"
    d1["Pathway"] <- rownames(d1)
    d2["Pathway"] <- rownames(d2)
    d <- merge(d1, d2)
    
    if (returnData) {
        #attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }
}

#extracting the loadings
loadingstransf <- loadingsfromRlog1.DESeqTransform(rld, intgroup = "type", returnData = TRUE)

#subsetting keeping the loadings of the pahways after filtering
reducedLoadings <- loadingstransf[loadingstransf$Pathway %in% names(columnsmorethanthree), ]

#exporting the pathway loading groups
write.table(data.frame(reducedLoadings[reducedLoadings$PC2 < 0,]$Pathway, reducedLoadings[reducedLoadings$PC2 < 0,]$PC2), file = "WheatLoadings", sep = "\t", quote = FALSE)
write.table(data.frame(reducedLoadings[reducedLoadings$PC2 > 0.05,]$Pathway, reducedLoadings[reducedLoadings$PC2 > 0.05,]$PC2), file = "GrassLoadings", sep = "\t", quote = FALSE)



##Factor Analysis of Pathways
reducedPathways <- trabsvarfiltpath[, names(as.data.frame(trabsvarfiltpath)) %in% names(columnsmorethanthree)]

#determining number of factors
fa.parallel(reducedPathways)
vss(reducedPathways, rotate = "promax")

#factor analysis function calling and diagram drawing
factan2 <- fa(cor(reducedPathways), nfactors = 2, rotate = "promax")
fa.sort(factan2)
#fa.diagram(factan2)
fa.graph(factan2, out.file = "FactorAnalysisDot.txt")



##applying filters (zero counts in samples and variance filtering) to pfam protein families

#loading the pfam protein family - read count table
pfammatrix <- read.table("pfamatrix.txt", sep = "\t")
pfammatrix <- pfammatrix[,-1]
pfammatrixcast <- dcast(pfammatrix, V2 ~ V4, fun.aggregate = sum, value.var = 'V3')

#filtering by zero counts
rm(columnsmorethanthree)
counter1 <- 0
columnsums <- colSums(as.matrix(pfammatrixcast[,-1]) != 0)
for(i in columnsums){
    counter1 <- counter1 + 1
    if(i > 3){
        if(exists("columnsmorethanthree")){
            columnsmorethanthree <- c(columnsmorethanthree, columnsums[counter1])
        } else{
            columnsmorethanthree <- columnsums[counter1]
        }
    }
}
rm(pfamfiltered)
for(i in names(columnsmorethanthree)){
    if (exists("pfamfiltered")){
        pfamfiltered <- c(pfamfiltered, pfammatrixcast[i])
    } else {
        pfamfiltered <- pfammatrixcast[i]
        }
}
pfamfiltered <- as.data.frame(pfamfiltered)
names(columnsmorethanthree)

#variance filtering for pfam
forgenefilterpfam <- new("ExpressionSet", exprs=as.matrix(t(pfammatrixcast[,-1])))
a11 <- varFilter(forgenefilterpfam)
varfiltpfam <- t(exprs(a11))

rm(columnsmorethanthree)
counter1 <- 0
columnsums <- colSums(as.matrix(varfiltpfam) != 0)
for(i in columnsums){
    counter1 <- counter1 + 1
    if(i > 3){
        if(exists("columnsmorethanthree")){
            columnsmorethanthree <- c(columnsmorethanthree, columnsums[counter1])
        } else{
            columnsmorethanthree <- columnsums[counter1]
        }
    }
}
names(columnsmorethanthree)


##calculation of ratios
#calculating proportion of reads in Carbon Metabolism pathways
CarbonMetabolism <- c(
"Glycolysis / Gluconeogenesis",
"Citrate cycle (TCA cycle)",
"Pentose phosphate pathway",
"Pentose and glucuronate interconversions",
"Fructose and mannose metabolism",
"Galactose metabolism",
"Ascorbate and aldarate metabolism",
"Starch and sucrose metabolism",
"Amino sugar and nucleotide sugar metabolism",
"Pyruvate metabolism",
"Glyoxylate and dicarboxylate metabolism",
"Propanoate metabolism",
"Butanoate metabolism",
"C5-Branched dibasic acid metabolism",
"Inositol phosphate metabolism")
pathcolsums <- colSums(as.matrix(pathmatrixcast[,-1]))
rm(carboncounts)
for(i in which(names(pathcolsums) %in% CarbonMetabolism)){
    if (exists("carboncounts")){
        carboncounts <- c(carboncounts, sum(pathmatrixcast[, i+1]))
    } else {
        carboncounts <- sum(pathmatrixcast[, i+1])
        }
}
sum(carboncounts)/sum(pathcolsums)

#calculating propostion of reads in Metabolic pathways
metabpathways <- c(
"Citrate cycle (TCA cycle)",
"Pentose phosphate pathway",
"Pentose and glucuronate interconversions",
"Fructose and mannose metabolism",
"Galactose metabolism",
"Ascorbate and aldarate metabolism",
"Starch and sucrose metabolism",
"Amino sugar and nucleotide sugar metabolism",
"Pyruvate metabolism",
"Glyoxylate and dicarboxylate metabolism",
"Propanoate metabolism",
"Butanoate metabolism",
"C5-Branched dibasic acid metabolism",
"Inositol phosphate metabolism	Enzymes",
"Compounds with biological roles",
"Oxidative phosphorylation",
"Photosynthesis",
"Photosynthesis - antenna proteins",
"Carbon fixation in photosynthetic organisms",
"Carbon fixation pathways in prokaryotes",
"Methane metabolism",
"Nitrogen metabolism",
"Sulfur metabolism",
"Photosynthesis proteins",
"Fatty acid biosynthesis",
"Fatty acid elongation",
"Fatty acid degradation",
"Synthesis and degradation of ketone bodies",
"Cutin, suberine and wax biosynthesis",
"Steroid biosynthesis",
"Primary bile acid biosynthesis",
"Secondary bile acid biosynthesis",
"Steroid hormone biosynthesis",
"Glycerolipid metabolism",
"Glycerophospholipid metabolism",
"Ether lipid metabolism",
"Sphingolipid metabolism",
"Arachidonic acid metabolism",
"Linoleic acid metabolism",
"alpha-Linolenic acid metabolism",
"Biosynthesis of unsaturated fatty acids",
"Lipid biosynthesis proteins",
"Purine metabolism",
"Pyrimidine metabolism",
"Alanine, aspartate and glutamate metabolism",
"Glycine, serine and threonine metabolism",
"Cysteine and methionine metabolism",
"Valine, leucine and isoleucine degradation",
"Valine, leucine and isoleucine biosynthesis",
"Lysine biosynthesis",
"Lysine degradation",
"Arginine biosynthesis",
"Arginine and proline metabolism",
"Histidine metabolism",
"Tyrosine metabolism",
"Phenylalanine metabolism",
"Tryptophan metabolism",
"Phenylalanine, tyrosine and tryptophan biosynthesis",
"beta-Alanine metabolism",
"Taurine and hypotaurine metabolism",
"Phosphonate and phosphinate metabolism",
"Selenocompound metabolism",
"Cyanoamino acid metabolism",
"D-Glutamine and D-glutamate metabolism",
"D-Arginine and D-ornithine metabolism",
"D-Alanine metabolism",
"Glutathione metabolism",
"N-Glycan biosynthesis",
"Various types of N-glycan biosynthesis",
"Mucin type O-glycan biosynthesis",
"Mannose type O-glycan biosynthesis New!",
"Other types of O-glycan biosynthesis",
"Glycosaminoglycan biosynthesis - CS/DS",
"Glycosaminoglycan biosynthesis - HS/Hep",
"Glycosaminoglycan biosynthesis - KS",
"Glycosaminoglycan degradation",
"Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
"Glycosphingolipid biosynthesis - lacto and neolacto series",
"Glycosphingolipid biosynthesis - globo and isoglobo series",
"Glycosphingolipid biosynthesis - ganglio series",
"Lipopolysaccharide biosynthesis",
"Peptidoglycan biosynthesis",
"Other glycan degradation",
"Thiamine metabolism",
"Riboflavin metabolism",
"Vitamin B6 metabolism",
"Nicotinate and nicotinamide metabolism",
"Pantothenate and CoA biosynthesis",
"Biotin metabolism",
"Lipoic acid metabolism",
"Folate biosynthesis",
"One carbon pool by folate",
"Retinol metabolism",
"Porphyrin and chlorophyll metabolism",
"Ubiquinone and other terpenoid-quinone biosynthesis",
"Terpenoid backbone biosynthesis",
"Monoterpenoid biosynthesis",
"Sesquiterpenoid and triterpenoid biosynthesis",
"Diterpenoid biosynthesis",
"Carotenoid biosynthesis",
"Brassinosteroid biosynthesis",
"Insect hormone biosynthesis",
"Zeatin biosynthesis",
"Limonene and pinene degradation",
"Geraniol degradation",
"Type I polyketide structures",
"Biosynthesis of 12-, 14- and 16-membered macrolides",
"Biosynthesis of ansamycins",
"Biosynthesis of enediyne antibiotics",
"Biosynthesis of type II polyketide backbone",
"Biosynthesis of type II polyketide products",
"Tetracycline biosynthesis",
"Polyketide sugar unit biosynthesis",
"Nonribosomal peptide structures",
"Biosynthesis of siderophore group nonribosomal peptides",
"Biosynthesis of vancomycin group antibiotics",
"Polyketide biosynthesis proteins",
"Phenylpropanoid biosynthesis",
"Stilbenoid, diarylheptanoid and gingerol biosynthesis",
"Flavonoid biosynthesis",
"Flavone and flavonol biosynthesis",
"Anthocyanin biosynthesis",
"Isoflavonoid biosynthesis",
"Indole alkaloid biosynthesis",
"Indole diterpene alkaloid biosynthesis",
"Isoquinoline alkaloid biosynthesis",
"Tropane, piperidine and pyridine alkaloid biosynthesis",
"Acridone alkaloid biosynthesis",
"Caffeine metabolism",
"Betalain biosynthesis",
"Glucosinolate biosynthesis",
"Benzoxazinoid biosynthesis",
"Penicillin and cephalosporin biosynthesis",
"Carbapenem biosynthesis",
"Monobactam biosynthesis",
"Clavulanic acid biosynthesis",
"Streptomycin biosynthesis",
"Neomycin, kanamycin and gentamicin biosynthesis",
"Acarbose and validamycin biosynthesis",
"Puromycin biosynthesis",
"Novobiocin biosynthesis",
"Staurosporine biosynthesis",
"Aflatoxin biosynthesis",
"Benzoate degradation",
"Aminobenzoate degradation",
"Fluorobenzoate degradation",
"Chloroalkane and chloroalkene degradation",
"Chlorocyclohexane and chlorobenzene degradation",
"Toluene degradation",
"Xylene degradation",
"Nitrotoluene degradation",
"Ethylbenzene degradation",
"Styrene degradation",
"Atrazine degradation",
"Caprolactam degradation",
"DDT degradation",
"Bisphenol degradation",
"Dioxin degradation",
"Naphthalene degradation",
"Polycyclic aromatic hydrocarbon degradation",
"Furfural degradation",
"Steroid degradation",
"Metabolism of xenobiotics by cytochrome P450",
"Drug metabolism - cytochrome P450",
"Drug metabolism - other enzymes",
"Overview of biosynthetic pathways",
"Biosynthesis of plant secondary metabolites",
"Biosynthesis of phenylpropanoids",
"Biosynthesis of terpenoids and steroids",
"Biosynthesis of alkaloids derived from shikimate pathway",
"Biosynthesis of alkaloids derived from ornithine, lysine and nicotinic acid",
"Biosynthesis of alkaloids derived from histidine and purine",
"Biosynthesis of alkaloids derived from terpenoid and polyketide",
"Biosynthesis of plant hormones")
rm(metabcounts)
for(i in which(names(pathcolsums) %in% metabpathways)){
    if (exists("metabcounts")){
        metabcounts <- c(metabcounts, sum(pathmatrixcast[, i+1]))
    } else {
        metabcounts <- sum(pathmatrixcast[, i+1])
        }
}
sum(metabcounts)/sum(pathcolsums)

#calculating the percentage of reads in glycoside hydrolace protein families
(439 + 213 + 124 )/sum(colSums(pfammatrixcast[,-1]))
