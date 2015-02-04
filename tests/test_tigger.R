# LOAD TIGGER

library(tigger)


# LOAD YOUR DATA

# clip_tab
load("C:/Users/Daniel Gadala-Maria/Documents/Kleinstein/Datasets/gmc2.rdata")
#dat = subset(combSeqData, donor == "Twin_147"  & FUNCTIONAL == "TRUE" & gene_identifier == "IGH|IGL|IGK")
dataset = subset(dat.gmc2, TIME %in% c("-1h", "-2d", "-8d") & FUNCTIONAL == "T")
dataset = dataset[which(!duplicated(dataset$SEQUENCE_GAP)),] # REMOVE DUPLICATES
genes = alakazam::getGene(dataset$V_CALL)
family = substr(genes, 1, 5)
dataset = subset(dataset, family %in% c("IGHV1"))
df = dataset[,colnames(dataset)[1:24]]
#save(df, file = "C:/Users/Daniel Gadala-Maria/repos/tiggerpackage/inst/extdata/df.rdata")

germline_db_file = "C:/Users/Daniel Gadala-Maria/Documents/Kleinstein/Datasets/IMGT/IMGT Variable 2014-12-22.fasta"
germline_db = readGermlineDb(germline_db_file, strip_down_name = TRUE)


# RUN TIGGER

result = runTigger(dataset, germline_db)


# CHECK OUT THE RESULT

# Novel sequences
# novel_sequences will contain only new alleles that made it into the genotype
novel_sequences = novelSummary(result)
plotNovelLines(result$novel[names(novel_sequences)])
plotNovelBars(result$novel[names(novel_sequences)])
plotJunctionBars(result$novel[names(novel_sequences)])

# Genotype
print(result$genotype)

# Corrected allele calls
V_CALL_GENOTYPED = result$new_calls
dat = cbind(dat, V_CALL_GENOTYPED)
head(dat)




