sample_db <- file.path("..", "tests-data", "sample_db.rda")
load(sample_db)

germline_ighv <- file.path("..", "tests-data", "germline_ighv.rda")
load(germline_ighv)

context("Evidence")

#ensure older version of sample() used
R_v <- paste(version$major, version$minor,sep=".")
w <- getOption("warn")
options(warn = -1)
if ( numeric_version(R_v) >= numeric_version("3.6.0") ) {
    RNGkind(sample.kind="Round")   
}
options(warn = w)

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
    skip_on_cran()
    # Find novel alleles and return relevant data
    novel_df <- findNovelAlleles(sample_db, 
                                 v_call="V_CALL",
                                 j_call="J_CALL",
                                 seq="SEQUENCE_IMGT",
                                 junction="JUNCTION",
                                 junction_length="JUNCTION_LENGTH",
                                 germline_ighv)
    geno <- inferGenotype(sample_db,
                          v_call="V_CALL",
                          seq = "SEQUENCE_IMGT",
                          germline_db = germline_ighv, 
                          novel = novel_df,
                          find_unmutated = TRUE)
    # Save the genotype sequences to a vector
    genotype_db <- genotypeFasta(geno, germline_db=germline_ighv, novel=novel_df)
    # Visualize the genotype and sequence counts
    sample_db <- reassignAlleles(sample_db,
                                 v_call="V_CALL", 
                                 seq="SEQUENCE_IMGT",
                                 genotype_db=genotype_db)
    # ev <- generateEvidence(geno, 
    #                  novel_df, 
    #                  c(germline_ighv[!names(germline_ighv) %in% names(genotype_seqs)], 
    #                    genotype_seqs),
    #                  germline_ighv,
    #                  sample_db)
    ev <- generateEvidence(data=sample_db, 
                           novel=novel_df, 
                           genotype=geno,
                           genotype_db=genotype_db,
                           germline_db=germline_ighv,
                           j_call="J_CALL", junction="JUNCTION")

    # Iterative, with 1 iteration, should match
    # novel_df_i <- itigger(sample_db, germline_ighv, max.iter = 1)
    # ev_i <- novel_df_i$summary

    expected_ev <- data.frame(
    'polymorphism_call'='IGHV1-8*02_G234T',
    'gene'='IGHV1-8',
    'total'=837,
    'allele'='02_G234T',
    'counts'=370,
    'germline_call'='IGHV1-8*02',
    'note'='Novel allele found!',
    'nt_substitutions'='234G>T',
    'novel_imgt'='CAGGTGCAGCTGGTGCAGTCTGGGGCT...GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTC............ACCAGCTATGATATCAACTGGGTGCGACAGGCCACTGGACAAGGGCTTGAGTGGATGGGATGGATGAACCCTAAC......AGTGGTAACACAGGCTATGCACAGAAGTTCCAG...GGCAGAGTCACCATTACCAGGAACACCTCCATAAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGG',
    'novel_imgt_count'=657,
    'novel_imgt_unique_j'=6,
    'novel_imgt_unique_cdr3'=626,
    'perfect_match_count'=661,
    'perfect_match_freq'=0.729580573951435,
    'germline_call_count'=906,
    'germline_call_freq'=0.052,
    'mut_min'=1,
    'mut_max'=10,
    'mut_pass_count'=760,
    'germline_imgt'='CAGGTGCAGCTGGTGCAGTCTGGGGCT...GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTC............ACCAGCTATGATATCAACTGGGTGCGACAGGCCACTGGACAAGGGCTTGAGTGGATGGGATGGATGAACCCTAAC......AGTGGTAACACAGGCTATGCACAGAAGTTCCAG...GGCAGAGTCACCATGACCAGGAACACCTCCATAAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGG',
    'germline_imgt_count'=0,
    'pos_min'=1,
    'pos_max'=312,
    'y_intercept'=0.125,
    'y_intercept_pass'=1,
    'snp_pass'=754,
    'unmutated_count'=661,
    'unmutated_freq'=0.729580573951435,
    'unmutated_snp_j_gene_length_count'=83,
    'snp_min_seqs_j_max_pass'=1,
    'alpha'=0.05,
    'min_seqs'=50,
    'j_max'=0.15,
    'min_frac'=0.75,
    'sequences'=864,
    'closest_reference'='IGHV1-8*02',
    'nt_diff'=1,
    'aa_diff'=1,
    'aa_substitutions'='78M>I',
    'unmutated_sequences'=370,
    'unmutated_frequency'=0.428240740740741,
    'allelic_percentage'=44.205495818399,
    'unique_js'=14,
    'unique_cdr3s'=747,
    'closest_reference_imgt'='CAGGTGCAGCTGGTGCAGTCTGGGGCT...GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTC............ACCAGCTATGATATCAACTGGGTGCGACAGGCCACTGGACAAGGGCTTGAGTGGATGGGATGGATGAACCCTAAC......AGTGGTAACACAGGCTATGCACAGAAGTTCCAG...GGCAGAGTCACCATGACCAGGAACACCTCCATAAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGG',
    stringsAsFactors = F
    )
    expect_equivalent(data.frame(ev[,colnames(expected_ev)], stringsAsFactors = F),
                      expected_ev)    
    
    # expect_equivalent(data.frame(ev_i[,colnames(expected_ev)], stringsAsFactors = F),
                      # expected_ev, tolerance=0.001)    

})

