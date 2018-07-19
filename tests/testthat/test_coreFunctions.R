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
    
})


test_that("Helper functions", {
    expect_true(hasNonImgtGaps("AT.ATT"))
    expect_null(getMutatedAA("...", "..."))
    expect_error(getMutatedAA("...", "NNN"), regexp = "Unexpected N in novel_imgt")
    expect_warning(expect_null(getMutatedAA("...", ".A.")), regexp="Non IMGT gaps found in novel_imgt")
    expect_equal(getMutatedAA("ATCAAATTC", "CGCAAAGTC"),
                 c("1I>R","3F>V"))
    expect_equal(getMutatedAA("...","ATG"), c("1X>M"))
    expect_warning(expect_equal(getMutatedAA("ATG",".A."), c("1M>X")))
})
# Run 3 iterations manually and compare results using itigger
# This test takes a long time, and is set to FALSE
if (FALSE) {
    test_that("Test itigger",{ 
    sample_db$FAMILY <- getFamily(sample_db$V_CALL, first=TRUE)
    sample_db <- sample_db %>%
        mutate(V_CALL_ORIGINAL=V_CALL) %>%
        select(V_CALL, FAMILY, V_CALL_ORIGINAL) %>%
        mutate(
            V_CALL=strsplit(as.character(V_CALL),",")
        ) %>%
        tidyr::unnest(V_CALL) %>%
        # filter(V_CALL==V_CALL_ORIGINAL) %>%
        group_by(FAMILY, V_CALL) %>%
        dplyr::summarise(V_CALL_FREQ=n()) %>%
        group_by(FAMILY) %>%
        mutate(V_CALL_FREQ=V_CALL_FREQ/sum(V_CALL_FREQ)) %>%
        arrange(desc(V_CALL_FREQ)) %>%
        select(-V_CALL_FREQ) %>%
        sample_n(1) %>%
        ungroup() %>%
        right_join(sample_db %>% select(-V_CALL))
    germline_subset <- germline_ighv[unique(sample_db$V_CALL)]
    
    ## Three iterations manually
    # 1
    nv1 <- findNovelAlleles(sample_db, germline_subset, v_call="V_CALL")
    gt1 <- inferGenotype(sample_db, "V_CALL", germline_db = germline_subset, novel_df = nv1)
    germdb1 <- genotypeFasta(gt1, germline_subset, nv1)
    sample_db1 <- sample_db
    sample_db1[['V_CALL_GENOTYPED']] <-  reassignAlleles(sample_db, germdb1, 
                                                         v_call="V_CALL")[['V_CALL_GENOTYPED']]
    germdb1 <- c(germdb1,
                 germline_subset[names(germline_subset) %in% names(germdb1) == F])
    
    setdiff(germdb1, germline_subset)
    setdiff(names(germdb1), names(germline_subset))
    
    # 2
    nv2 <- findNovelAlleles(sample_db1, germdb1, v_call="V_CALL_GENOTYPED")
    gt2 <- inferGenotype(sample_db1, "V_CALL_GENOTYPED", germline_db = germdb1, novel_df = nv2)
    germdb2 <- genotypeFasta(gt2, germdb1, nv2)
    sample_db2 <- sample_db1
    sample_db2[['V_CALL_GENOTYPED']] <-  reassignAlleles(sample_db1, germdb2, 
                                                         v_call="V_CALL_GENOTYPED")[['V_CALL_GENOTYPED']] 
    germdb2 <- c(germdb2,
                 germdb1[names(germdb1) %in% names(germdb2) == F])
    
    setdiff(germdb2, germdb1)
    setdiff(names(germdb2), names(germdb1))
    
    # 3
    nv3 <- findNovelAlleles(sample_db2, germdb2, v_call="V_CALL_GENOTYPED")
    gt3 <- inferGenotype(sample_db2, "V_CALL_GENOTYPED", germline_db = germdb2, novel_df = nv3)
    germdb3 <- genotypeFasta(gt3, germdb2, nv3)
    sample_db3 <- sample_db2
    sample_db3[['V_CALL_GENOTYPED']] <-  reassignAlleles(sample_db2, germdb3, 
                                                         v_call="V_CALL_GENOTYPED")[['V_CALL_GENOTYPED']] 
    # No new germlines in 3rd iteration
    germdb3 <- c(germdb3,
                 germdb2[names(germdb2) %in% names(germdb3) == F])
    setdiff(germdb3, germdb2)
    setdiff(names(germdb3), names(germdb2))
    
    
    # This should stop at 3rd iteration
    inv <- itigger(sample_db, germline_subset, fields=NULL, nproc=1, max.iter = 4)
    
    
    expect_equivalent(inv$nv[inv$nv$ITERATION=="1",colnames(nv1)], nv1)
    expect_equivalent(inv$nv[inv$nv$ITERATION=="2",colnames(nv2)], nv2)
    expect_equivalent(inv$nv[inv$nv$ITERATION=="3",colnames(nv3)], nv3)
    
    expect_equivalent(sort(names(inv$germline[["1"]])),sort(names(germdb2)))
    expect_equivalent(sort(names(inv$germline[["1"]])),sort(names(germdb2)))
    })
    
    test_that("Test itigger PGP1 MiSeq",{ 
        
        load(file.path("..", "tests-data", "pgp1.RData"))
        load(file.path("..", "tests-data", "pgp1_germlines.RData"))
        
        allele_count <- sapply(names(pgp1_germlines), function(g) {
            sum(grepl(pgp1_germlines[g], pgp1$SEQUENCE_IMGT, fixed = T))
        })
        allele_count <- data.frame(
            allele=names(allele_count),
            count=allele_count,
            family=alakazam::getFamily(names(allele_count)))
        germ_alleles <- allele_count %>%
            dplyr::arrange(desc(count)) %>%
            dplyr::filter(family != "IGHV1") %>%
            dplyr::group_by(family) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup() %>%
            select(allele) %>%
            rbind(.,"IGHV1-18*01")
            
        germlines <- pgp1_germlines[germ_alleles[['allele']]]
        
        tigger_results <- itigger(pgp1, germlines[], max.iter=Inf)
        names(pgp1_germlines)[match(tigger_results$new_germlines, cleanSeqs(pgp1_germlines))]
        
        tigger_results_50_25 <- itigger(pgp1, germlines[], max.iter=Inf, germline_min=50, min_seqs=25)
        names(pgp1_germlines)[match(tigger_results_50_25$new_germlines, cleanSeqs(pgp1_germlines))]
    })
}

    