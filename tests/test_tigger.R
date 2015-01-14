# LOAD TIGGER

library(tigger)


# LOAD YOUR DATA

# clip_tab
load("C:/Users/Daniel Gadala-Maria/Documents/Kleinstein/Datasets/Bolen Twins//twinStudy_runs5to17_allSeqs_withFRdata.Rd")
ddply(combSeqData, "donor", nrow)
dat = subset(combSeqData, donor == "Twin_147"  & FUNCTIONAL == "TRUE" & gene_identifier == "IGH|IGL|IGK")
#dat = subset(dat.gmc2, TIME %in% c("-1h", "-2d", "-8d") & FUNCTIONAL == "T")
dat = dat[which(!duplicated(dat$SEQUENCE_GAP)),] # REMOVE DUPLICATES
rm(combSeqData)
germline_db_file = "C:/Users/Daniel Gadala-Maria/Documents/Kleinstein/Datasets/IMGT/IMGT Variable 2014-12-22.fasta"
germline_db = readGermlineDb(germline_db_file, strip_down_name = TRUE)


# RUN TIGGER

result = tigger(dat, germline_db)


# CHECK OUT THE RESULT

print(names(novel))
plotNovelLines(result$novel)
plotNovelBars(result$novel)
plotJunctionBars(result$novel) 
print(result$genotype)
V_CALL_GENOTYPED = result$new_calls
dat = cbind(dat, V_CALL_GENOTYPED)
head(dat)










