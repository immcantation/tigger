sample_db <- file.path("..", "tests-data", "sample_db.rda")
load(sample_db)

germline_ighv <- file.path("..", "tests-data", "germline_ighv.rda")
load(germline_ighv)

context("Evidence")

test_that("Helper functions", {
    expect_true(hasNonImgtGaps("AT.ATT"))
    expect_null(getMutatedAA("...", "..."))
    expect_error(getMutatedAA("...", "NNN"), regexp = "Unexpected N in novel_imgt")
    expect_warning(expect_null(getMutatedAA("...", ".A.")), 
                   regexp="Non IMGT gaps found in novel_imgt")
    expect_equal(getMutatedAA("ATCAAATTC", "CGCAAAGTC"),
                 c("1I>R","3F>V"))
    expect_equal(getMutatedAA("...","ATG"), c("1X>M"))
    expect_warning(expect_equal(getMutatedAA("ATG",".A."), c("1M>X")))
})


test_that("generateEvidence", {
    
    # Find novel alleles and return relevant data
    novel_df <- findNovelAlleles(sample_db, germline_ighv)
    geno <- inferGenotype(sample_db, find_unmutated = TRUE,
                          germline_db = germline_ighv, 
                          novel_df = novel_df)
    # Save the genotype sequences to a vector
    genotype_seqs <- genotypeFasta(geno, germline_ighv, novel_df)
    # Visualize the genotype and sequence counts
    sample_db <- reassignAlleles(sample_db, genotype_seqs)
    ev <- generateEvidence(geno, 
                     novel_df, 
                     c(germline_ighv[!names(germline_ighv) %in% names(genotype_seqs)], 
                       genotype_seqs),
                     germline_ighv,
                     sample_db)
    
    # Iterative, with 1 iteration, should match
    novel_df_i <- itigger(sample_db, germline_ighv, max.iter = 1)
    ev_i <- novel_df_i$summary

    expected_ev <- data.frame(
    'POLYMORPHISM_CALL'='IGHV1-8*02_G234T',
    'GENE'='IGHV1-8',
    'TOTAL'=837,
    'NOTE_GT'='',
    'ALLELE'='02_G234T',
    'COUNTS'=370,
    'GERMLINE_CALL'='IGHV1-8*02',
    'NOTE'='Novel allele found!',
    'NT_SUBSTITUTIONS'='234G>T',
    'NOVEL_IMGT'='CAGGTGCAGCTGGTGCAGTCTGGGGCT---GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTC------------ACCAGCTATGATATCAACTGGGTGCGACAGGCCACTGGACAAGGGCTTGAGTGGATGGGATGGATGAACCCTAAC------AGTGGTAACACAGGCTATGCACAGAAGTTCCAG---GGCAGAGTCACCATTACCAGGAACACCTCCATAAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGG',
    'NOVEL_IMGT_COUNT'=657,
    'NOVEL_IMGT_UNIQUE_J'=6,
    'NOVEL_IMGT_UNIQUE_CDR3'=626,
    'PERFECT_MATCH_COUNT'=661,
    'PERFECT_MATCH_FREQ'=0.729580573951435,
    'GERMLINE_CALL_COUNT'=906,
    'GERMLINE_CALL_PERC'=5.2,
    'MUT_MIN'=1,
    'MUT_MAX'=10,
    'MUT_PASS_COUNT'=760,
    'GERMLINE_IMGT'='CAGGTGCAGCTGGTGCAGTCTGGGGCT---GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTC------------ACCAGCTATGATATCAACTGGGTGCGACAGGCCACTGGACAAGGGCTTGAGTGGATGGGATGGATGAACCCTAAC------AGTGGTAACACAGGCTATGCACAGAAGTTCCAG---GGCAGAGTCACCATGACCAGGAACACCTCCATAAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGG',
    'GERMLINE_IMGT_COUNT'=0,
    'POS_MIN'=1,
    'POS_MAX'=312,
    'Y_INTERCEPT'=0.125,
    'Y_INTERCEPT_PASS'=1,
    'SNP_PASS'=754,
    'UNMUTATED_COUNT'=661,
    'UNMUTATED_FREQ'=0.729580573951435,
    'UNMUTATED_SNP_J_GENE_LENGTH_COUNT'=83,
    'SNP_MIN_SEQS_J_MAX_PASS'=1,
    'ALPHA'=0.05,
    'MIN_SEQS'=50,
    'J_MAX'=0.15,
    'MIN_FRAC'=0.75,
    'SEQUENCES'=864,
    'CLOSEST_REFERENCE'='IGHV1-8*02',
    'NT_DIFF'=1,
    'AA_DIFF'=1,
    'AA_SUBSTITUTIONS'='78M>I',
    'UNMUTATED_SEQUENCES'=370,
    'UNMUTATED_FREQUENCY'=0.428240740740741,
    'ALLELIC_PERCENTAGE'=44.205495818399,
    'UNIQUE_JS'=14,
    'UNIQUE_CDR3S'=747,
    'CLOSEST_REFERENCE_IMGT'='CAGGTGCAGCTGGTGCAGTCTGGGGCT---GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTC------------ACCAGCTATGATATCAACTGGGTGCGACAGGCCACTGGACAAGGGCTTGAGTGGATGGGATGGATGAACCCTAAC------AGTGGTAACACAGGCTATGCACAGAAGTTCCAG---GGCAGAGTCACCATGACCAGGAACACCTCCATAAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGG',
    stringsAsFactors = F
    )
    expect_equivalent(data.frame(ev[,colnames(expected_ev)], stringsAsFactors = F),
                      expected_ev)    
    
    expect_equivalent(data.frame(ev_i[,colnames(expected_ev)], stringsAsFactors = F),
                      expected_ev, tolerance=0.001)    

})

