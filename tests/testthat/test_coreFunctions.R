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
    expect_false("genotyped_alleles" %in% colnames(geno_bay_airr))

    geno_bay_gt <- inferGenotypeBayesian(airr_db,
                                      germline_db = germline_ighv,
                                      novel = novel_df_airr,
                                      v_call="v_call",
                                      seq="sequence_alignment",
                                      genotyped_alleles=TRUE)
    expect_true("genotyped_alleles" %in% colnames(geno_bay_gt))
    
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

test_that("findNovelAlleles pos_range_max column is dereferenced correctly", {
    # airr_db and germline_ighv are loaded at the top of this file.
    # pos_range 310:320 overlaps with the v_germ_length boundary (values 310-320),
    # so the POSITION filter has observable effect: sequences whose v_germ_length
    # is shorter than a passing position are excluded.
    nv_null <- findNovelAlleles(airr_db, germline_ighv,
                                v_call="v_call", j_call="j_call",
                                seq="sequence_alignment",
                                junction="junction",
                                junction_length="junction_length",
                                pos_range=310:320,
                                pos_range_max=NULL)
    nv_with_max <- findNovelAlleles(airr_db, germline_ighv,
                                    v_call="v_call", j_call="j_call",
                                    seq="sequence_alignment",
                                    junction="junction",
                                    junction_length="junction_length",
                                    pos_range=310:320,
                                    pos_range_max="v_germ_length")

    # Both calls must complete without error
    expect_s3_class(nv_null, "data.frame")
    expect_s3_class(nv_with_max, "data.frame")

    # Applying pos_range_max must change results, proving the filter is applied
    expect_false(identical(nrow(nv_null), nrow(nv_with_max)))
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

    gt$genotyped_alleles <- c("04", "01", "01", "09")
    gtfa_unseen <- genotypeFasta(gt, germline_db, include_unseen=TRUE)
    expect_equal(gtfa_unseen, germline_db[c(6, 1, 3, 4, 5)])
    
    expect_error(genotypeFasta(gt, germline_db[-1]),
                 regexp="IGHV1-2\\*04")
    
})

test_that("inferGenotype uses locus-specific fractional gene_cutoff", {
    db <- data.frame(
        v_call=c(rep("IGHV1-1*01", 10), rep("IGKV1-1*01", 10)),
        sequence_alignment=rep("AAAA", 20),
        locus=c(rep("IGH", 10), rep("IGK", 10)),
        stringsAsFactors=FALSE
    )

    expect_warning(
        geno <- inferGenotype(db, find_unmutated=FALSE, gene_cutoff=0.75),
        regexp="Mixed loci detected"
    )
    expect_equal(sort(geno$gene), c("IGHV1-1", "IGKV1-1"))

    db$locus <- NULL
    expect_warning(
        geno_from_call <- inferGenotype(db, find_unmutated=FALSE, gene_cutoff=0.75),
        regexp="Mixed loci detected"
    )
    expect_equal(sort(geno_from_call$gene), c("IGHV1-1", "IGKV1-1"))
})

test_that("reassignAlleles supports segment output, overwrite, and trimming", {
    db_v <- data.frame(
        v_call="IGHV1-1*01",
        sequence_alignment="AAAT",
        stringsAsFactors=FALSE
    )
    genotype_v <- c("IGHV1-1*01"="AAAA", "IGHV1-1*02"="AAAT")
    reassigned_v <- reassignAlleles(db_v, genotype_v,
                                    v_call="v_call",
                                    seq="sequence_alignment")
    expect_equal(reassigned_v$v_call_genotyped, "IGHV1-1*02")

    overwritten_v <- reassignAlleles(db_v, genotype_v,
                                     v_call="v_call",
                                     seq="sequence_alignment",
                                     overwrite=TRUE)
    expect_equal(overwritten_v$v_call, "IGHV1-1*02")
    expect_false("v_call_genotyped" %in% colnames(overwritten_v))

    db_d <- data.frame(
        d_call="IGHD1-1*01",
        sequence_alignment="CCCA",
        stringsAsFactors=FALSE
    )
    genotype_d <- c("IGHD1-1*01"="CCCC", "IGHD1-1*02"="CCCA")
    reassigned_d <- reassignAlleles(db_d, genotype_d,
                                    v_call="d_call",
                                    seq="sequence_alignment")
    expect_equal(reassigned_d$d_call_genotyped, "IGHD1-1*02")

    db_j <- data.frame(
        j_call="IGHJ1*01",
        sequence_alignment="TTTA",
        stringsAsFactors=FALSE
    )
    genotype_j <- c("IGHJ1*01"="TTTT", "IGHJ1*02"="TTTA")
    reassigned_j <- reassignAlleles(db_j, genotype_j,
                                    v_call="j_call",
                                    seq="sequence_alignment")
    expect_equal(reassigned_j$j_call_genotyped, "IGHJ1*02")

    db_trim <- data.frame(
        v_call="IGHV1-1*01",
        sequence_alignment="GGAAAT",
        v_germline_start=3,
        v_germline_end=6,
        stringsAsFactors=FALSE
    )
    genotype_trim <- c("IGHV1-1*01"="GGAAAA", "IGHV1-1*02"="CCAAAT")
    reassigned_trim <- reassignAlleles(db_trim, genotype_trim,
                                       v_call="v_call",
                                       seq="sequence_alignment",
                                       trim_seq=TRUE)
    expect_equal(reassigned_trim$v_call_genotyped, "IGHV1-1*02")

    expect_error(
        reassignAlleles(db_v, genotype_v, trim_seq=TRUE),
        regexp="missing columns"
    )
})

test_that("reassignAlleles Rcpp mismatch path matches fallback results", {
    skip_if_not(
        exists("seqMismatchCountRcpp", envir=asNamespace("alakazam"), inherits=FALSE) &&
            exists("seqMismatchMatrixRcpp", envir=asNamespace("alakazam"), inherits=FALSE),
        "Alakazam Rcpp mismatch functions are not installed"
    )

    samples <- c("ACGT", "ACNT", "AC-T", "acgt", "AC.T")
    germlines <- c(g1="ACGA", g2="ACGT", g3="TCGT")
    fallback <- sapply(germlines, function(x) {
        sapply(getMutatedPositions(samples, x, ignored_regex="[\\.N-]",
                                   match_instead=FALSE), length)
    })
    rcpp <- get("seqMismatchMatrixRcpp", envir=asNamespace("alakazam"))(
        samples, germlines, ignore=c(".", "N", "-"))
    expect_equal(unname(rcpp), unname(fallback))

    data_trim <- data.frame(
        sequence_alignment=c("GGAAAT", "TTACGT"),
        v_germline_start=c(3, 3),
        v_germline_end=c(6, 6),
        stringsAsFactors=FALSE
    )
    samples_trim <- substr(data_trim$sequence_alignment,
                           data_trim$v_germline_start,
                           data_trim$v_germline_end)
    germlines_trim <- c(g1="GGAAAA", g2="TTACGT", g3="CCAAAT")
    fallback_trim <- sapply(germlines_trim, function(x) {
        ref <- substr(rep(x, nrow(data_trim)), data_trim$v_germline_start,
                      data_trim$v_germline_end)
        sapply(getMutatedPositions(samples_trim, ref,
                                   ignored_regex="[\\.N-]",
                                   match_instead=FALSE), length)
    })
    rcpp_trim <- sapply(germlines_trim, function(x) {
        ref <- substr(rep(x, nrow(data_trim)), data_trim$v_germline_start,
                      data_trim$v_germline_end)
        get("seqMismatchCountRcpp", envir=asNamespace("alakazam"))(
            samples_trim, ref, ignore=c(".", "N", "-"))
    })
    expect_equal(unname(rcpp_trim), unname(fallback_trim))

    db <- data.frame(
        v_call=c("IGHV1-1*01", "IGHV1-1*01", "IGHV1-1*01"),
        sequence_alignment=c("AAAT", "AAGT", "AAAC"),
        stringsAsFactors=FALSE
    )
    genotype_db <- c("IGHV1-1*01"="AAAA", "IGHV1-1*02"="AAAT")
    old_opt <- getOption("tigger.use_alakazam_rcpp_mismatch")
    on.exit(options(tigger.use_alakazam_rcpp_mismatch=old_opt), add=TRUE)

    options(tigger.use_alakazam_rcpp_mismatch=FALSE)
    fallback_db <- reassignAlleles(db, genotype_db,
                                   v_call="v_call",
                                   seq="sequence_alignment")
    options(tigger.use_alakazam_rcpp_mismatch=TRUE)
    rcpp_db <- reassignAlleles(db, genotype_db,
                               v_call="v_call",
                               seq="sequence_alignment")
    expect_equal(rcpp_db$v_call_genotyped, fallback_db$v_call_genotyped)

    db_trim <- data.frame(
        v_call=c("IGHV1-1*01", "IGHV1-1*01"),
        sequence_alignment=c("GGAAAT", "TTAAAA"),
        v_germline_start=c(3, 3),
        v_germline_end=c(6, 6),
        stringsAsFactors=FALSE
    )
    genotype_trim <- c("IGHV1-1*01"="GGAAAA", "IGHV1-1*02"="CCAAAT")
    options(tigger.use_alakazam_rcpp_mismatch=FALSE)
    fallback_trim_db <- reassignAlleles(db_trim, genotype_trim,
                                        v_call="v_call",
                                        seq="sequence_alignment",
                                        trim_seq=TRUE)
    options(tigger.use_alakazam_rcpp_mismatch=TRUE)
    rcpp_trim_db <- reassignAlleles(db_trim, genotype_trim,
                                    v_call="v_call",
                                    seq="sequence_alignment",
                                    trim_seq=TRUE)
    expect_equal(rcpp_trim_db$v_call_genotyped,
                 fallback_trim_db$v_call_genotyped)
})
