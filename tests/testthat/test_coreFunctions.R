sample_db <- file.path("..", "tests-data", "sample_db.rda")
load(sample_db)

germline_ighv <- file.path("..", "tests-data", "germline_ighv.rda")
load(germline_ighv)

test_that("Test findNovelAlleles",{ 
    novel_df <- findNovelAlleles(sample_db, germline_ighv)
    expect_equal(novel_df$NOTE[5], "Novel allele found!")
})
