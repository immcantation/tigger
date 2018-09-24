sample_db <- file.path("..", "tests-data", "sample_db.rda")
load(sample_db)

germline_ighv <- file.path("..", "tests-data", "germline_ighv.rda")
load(germline_ighv)

context("Evidence")
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