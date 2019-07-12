sample_db <- file.path("..", "tests-data", "sample_db.rda")
load(sample_db)

germline_ighv <- file.path("..", "tests-data", "germline_ighv.rda")
load(germline_ighv)

context("Core functions")

test_that("Test findNovelAlleles",{ 
    novel_df <- findNovelAlleles(sample_db, germline_ighv)
    expect_equal(novel_df$NOTE[5], "Novel allele found!")
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
        "V_CALL"=c("IGHV1-2*01","IGHV1-2*02","IGHV1-2*01,IGHV1-2*02","IGHV1-2*03"),
        "SAMPLE"=c("S1","S1","S2","S3")
    )
    
    set.seed(5)
    expect_equivalent(db[1,], subsampleDb(db, min_n = 1, max_n=1, mode="gene"))
    
    set.seed(5)
    expect_equivalent(db[c(1,3,4),], 
                 subsampleDb(db, min_n = 1, max_n=1, mode="allele"))
    
    set.seed(5)
    expect_equivalent(db[c(1,2,3,4),], 
                 subsampleDb(db, min_n = 1, max_n=1, mode="allele", group = "SAMPLE"))
    
    set.seed(5)
    expect_equivalent(db[1,], 
                 subsampleDb(db, min_n = 1, max_n=1, mode="family"))
    
    set.seed(5)
    expect_equivalent(db[c(1,3,4),], 
                 subsampleDb(db, min_n = 1, max_n=1, mode="family", group="SAMPLE"))    
    
    set.seed(5)
    expect_equivalent(db[c(1, 3),], subsampleDb(db, min_n = 1, max_n=2))
    

})
