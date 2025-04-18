sample_db <- file.path("..", "tests-data", "sample_db.rda")
load(sample_db)

airr_db <- file.path("..", "tests-data", "airr_db.rda")
load(airr_db)

SampleDbPosRangeMax <- file.path("..", "tests-data", "SampleDbPosRangeMax.rda")
load(SampleDbPosRangeMax)


germline_ighv <- file.path("..", "tests-data", "germline_ighv.rda")
load(germline_ighv)

context("Core functions")

#ensure older version of sample() used
R_v <- paste(version$major, version$minor,sep=".")
w <- getOption("warn")
options(warn = -1)
if ( numeric_version(R_v) >= numeric_version("3.6.0") ) {
    RNGkind(sample.kind="Round")   
}
options(warn = w)

test_that("Test findNovelAlleles",{ 
    expect_error(novel_df <- findNovelAlleles(sample_db, germline_ighv))
    novel_df <- findNovelAlleles(sample_db, germline_ighv,
                            v_call="V_CALL", j_call="J_CALL",
                            seq="SEQUENCE_IMGT",
                            junction = "JUNCTION",
                            junction_length = "JUNCTION_LENGTH")
    
    expect_equal(selectNovel(novel_df)$polymorphism_call, "IGHV1-8*02_G234T")
        
    novel_df_airr <- findNovelAlleles(airr_db, germline_ighv,
                                      v_call="v_call", j_call="j_call",
                                      seq = "sequence_alignment",
                                      junction="junction",
                                      junction_length = "junction_length")
    expect_equivalent(novel_df, novel_df_airr)
    
    geno <- inferGenotype(sample_db,
                          v_call="V_CALL",
                          seq = "SEQUENCE_IMGT",
                          germline_db = germline_ighv, 
                          novel = novel_df,
                          find_unmutated = TRUE)
    
    geno_airr <- inferGenotype(airr_db,
                          v_call="v_call",
                          seq = "sequence_alignment",
                          germline_db = germline_ighv, 
                          novel = novel_df_airr,
                          find_unmutated = TRUE)
    expect_equivalent(geno, geno_airr)
    
    geno_bay <- inferGenotypeBayesian(sample_db,
                                      germline_db = germline_ighv,
                                      novel = novel_df,
                                      v_call="V_CALL", 
                                      seq="SEQUENCE_IMGT")
    geno_bay_airr <- inferGenotypeBayesian(airr_db,
                                      germline_db = germline_ighv,
                                      novel = novel_df_airr,
                                      v_call="v_call", 
                                      seq="sequence_alignment")
    expect_equivalent(geno_bay, geno_bay_airr)
    
})


test_that("Test findNovelAlleles - pos_range_max",{ 
    nv_pos_range_318 <- findNovelAlleles(SampleDbPosRangeMax, 
                                         germline_db = SampleGermlineIGHV,
                                         pos_range=315:318)
    nv_pos_range_318_vend <- findNovelAlleles(SampleDbPosRangeMax, 
                                              germline_db = SampleGermlineIGHV,
                                              pos_range=315:318,
                                              pos_range_max="v_germline_end")   
    # Finds false positive, position 318.Not found when using pos_range_max="v_germline_end"
    expect_equal(selectNovel(nv_pos_range_318)[['pos_max']],318)
    expect_equal(nrow(selectNovel(nv_pos_range_318_vend)),0)
})

test_that("Test sortAlleles",{ 
    alleles = c("IGHV1-69D*01","IGHV1-69*01","IGHV1-2*01","IGHV1-69-2*01",
                "IGHV2-5*01","IGHV1-NL1*01", "IGHV1-2*01,IGHV1-2*05", 
                "IGHV1-2", "IGHV1-2*02", "IGHV1-69*02",
                "IGHV1S10*01", "IGHV1S1*01",
                "IGHV1-2*02_G234T"
    )
    sorted_alleles <- sortAlleles(alleles)
    expect_equal(sorted_alleles,
                 c( "IGHV1-2", "IGHV1-2*01", "IGHV1-2*01,IGHV1-2*05",
                    "IGHV1-2*02", "IGHV1-2*02_G234T", "IGHV1-69*01",
                    "IGHV1-69D*01", "IGHV1-69*02", "IGHV1-69-2*01",
                    "IGHV1-NL1*01", "IGHV1S1*01", "IGHV1S10*01",
                    "IGHV2-5*01"))
    
    expect_equal(
        sortAlleles(c("TRAV38-2/DV8", "TRAV38-1")),
        c("TRAV38-1","TRAV38-2/DV8"))
})

test_that("subsampleDb",{ 
    
    db <- data.frame(
        "v_call"=c("IGHV1-2*01","IGHV1-2*02","IGHV1-2*01,IGHV1-2*02","IGHV1-2*03"),
        "sample"=c("S1","S1","S2","S3")
    )
    
    set.seed(5)
    expect_equivalent(db[1,], subsampleDb(db, min_n = 1, max_n=1, mode="gene"))
    
    set.seed(5)
    expect_equivalent(db[c(1,3,4),], 
                 subsampleDb(db, min_n = 1, max_n=1, mode="allele"))
    
    set.seed(5)
    expect_equivalent(db[c(1,2,3,4),], 
                 subsampleDb(db, min_n = 1, max_n=1, mode="allele", group = "sample"))
    
    set.seed(5)
    expect_equivalent(db[1,], 
                 subsampleDb(db, min_n = 1, max_n=1, mode="family"))
    
    set.seed(5)
    expect_equivalent(db[c(1,3,4),], 
                 subsampleDb(db, min_n = 1, max_n=1, mode="family", group="sample"))    
    
    set.seed(5)
    expect_equivalent(db[c(1, 3),], subsampleDb(db, min_n = 1, max_n=2))

})


test_that("selectNovel keep_alleles keeps or removes alleles leading to the same novel sequence",{ 
    
    nv <- data.frame(
        list(
            "germline_call"=c("IGHV1-69*13", "IGHV1-69*14"),
            "polymorphism_call"=c("IGHV1-69*13_G244A","IGHV1-69*14_G54A"),
            "note"=c("Novel allele found!. Same as: IGHV1-69*14_G54A",
                     "Novel allele found!. Same as: IGHV1-69*13_G244A"),
            "novel_imgt"=c("CAGGTCCAGCTGGTGCAGTCTGGGGCT...GAGGTGAAGAAGCCTGGGTCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGAGGCACCTTC............AGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGAGGGATCATCCCTATC......TTTGGTACAGCAAACTACGCACAGAAGTTCCAG...GGCAGAGTCACGATTACCGCGGACAAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGA",
                           "CAGGTCCAGCTGGTGCAGTCTGGGGCT...GAGGTGAAGAAGCCTGGGTCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGAGGCACCTTC............AGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGAGGGATCATCCCTATC......TTTGGTACAGCAAACTACGCACAGAAGTTCCAG...GGCAGAGTCACGATTACCGCGGACAAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGA")
        ),
        stringsAsFactors = F
    )
    
    expect_equal(nrow(selectNovel(nv, keep_alleles = F)),1)
    expect_equal(nrow(selectNovel(nv, keep_alleles = T)),2)
})

test_that("Test genotypeFasta",{ 
    gt <- data.frame(
        "gene"=c("IGHV1-2", "IGHV3-23", "IGHV3-23D", "IGHV3-64D"),
        "alleles"=c("04,05","01","01","09"),
        "counts"=c("7,4","11","11","1"),
        "total"=c(11,11,11,1),
        "note"=c("","","",""),
        stringsAsFactors = F
    )
    # Dummy data. Sequence is not needed.
    germline_db <- c(
        "IGHV1-2*04"="A",
        "IGHV1-2*05"="A",
        "IGHV3-23*01"="C",
        "IGHV3-23D*01"="T",
        "IGHV3-64D*09"="G",
        "IGHV7-81*01"="T"
    )
    gtfa <- genotypeFasta(gt, germline_db)
    expect_equal(gtfa, germline_db[1:5])
    
    expect_error(genotypeFasta(gt, germline_db[-1]),
                 regexp="IGHV1-2\\*04")
    
})