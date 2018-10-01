context("itigger")

sample_db <- file.path("..", "tests-data", "sample_db.rda")
load(sample_db)

germline_ighv <- file.path("..", "tests-data", "germline_ighv.rda")
load(germline_ighv)

# Run 3 iterations manually and compare results using itigger
# This test takes a long time, skip on CRAN
# if (FALSE) {
    test_that("itigger and findNovelAlleles find same alleles",{ 
        skip_on_cran()
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
            dplyr::sample_n(1) %>%
            ungroup() %>%
            dplyr::right_join(sample_db %>% select(-V_CALL))
        germline_subset <- germline_ighv[unique(sample_db$V_CALL)]
        
        ## Three iterations manually
        # 1
        nv1 <- findNovelAlleles(sample_db, germline_subset, v_call="V_CALL")
        gt1 <- inferGenotype(sample_db, germline_db = germline_subset, novel = nv1, v_call="V_CALL")
        germdb1 <- genotypeFasta(gt1, germline_subset, nv1)
        sample_db1 <- sample_db
        sample_db1[['V_CALL_GENOTYPED']] <-  reassignAlleles(sample_db, germdb1, 
                                                             v_call="V_CALL")[['V_CALL_GENOTYPED']]
        germdb1 <- c(germdb1,
                     germline_subset[names(germline_subset) %in% names(germdb1) == F])
        
        
        # Evidence should match evidence from itigger 1 iteration
        ev <- generateEvidence(data=sample_db1, 
                               novel=nv1, 
                               genotype=gt1,
                               genotype_db=genotypeFasta(gt1, germline_subset, nv1),
                               germline_db=germline_subset)
        inv <- itigger(sample_db, germline_subset, fields=NULL, nproc=1, max.iter = 1)
        ev_obs <- inv$summary[inv$summary$ITERATION=="1",colnames(ev)]
        expect_equivalent(ev_obs, ev1)
        
        # setdiff(germdb1, germline_subset)
        # setdiff(names(germdb1), names(germline_subset))
        
        # 2
        nv2 <- findNovelAlleles(sample_db1, germdb1, v_call="V_CALL_GENOTYPED")
        gt2 <- inferGenotype(sample_db1, germline_db = germdb1, novel = nv2, v_call="V_CALL_GENOTYPED")
        germdb2 <- genotypeFasta(gt2, germdb1, nv2)
        sample_db2 <- sample_db1
        sample_db2[['V_CALL_GENOTYPED']] <-  reassignAlleles(sample_db1, germdb2, 
                                                             v_call="V_CALL_GENOTYPED")[['V_CALL_GENOTYPED']] 
        germdb2 <- c(germdb2,
                     germdb1[names(germdb1) %in% names(germdb2) == F])
        
        # setdiff(germdb2, germdb1)
        # setdiff(names(germdb2), names(germdb1))
        
        # 3
        nv3 <- findNovelAlleles(sample_db2, germdb2, v_call="V_CALL_GENOTYPED")
        gt3 <- inferGenotype(sample_db2, germline_db = germdb2, novel = nv3, v_call="V_CALL_GENOTYPED")
        germdb3 <- genotypeFasta(gt3, germdb2, nv3)
        sample_db3 <- sample_db2
        sample_db3[['V_CALL_GENOTYPED']] <-  reassignAlleles(sample_db2, germdb3, 
                                                             v_call="V_CALL_GENOTYPED")[['V_CALL_GENOTYPED']] 
        # No new germlines in 3rd iteration
        germdb3 <- c(germdb3,
                     germdb2[names(germdb2) %in% names(germdb3) == F])
        # setdiff(germdb3, germdb2)
        # setdiff(names(germdb3), names(germdb2))
        
        # Iterative novel
        # This should stop at 3rd iteration
        inv <- itigger(sample_db, germline_subset, fields=NULL, nproc=1, max.iter = 4)
        
        expect_equivalent(inv$novel[inv$novel$ITERATION=="1",colnames(nv1)], nv1)
        expect_equivalent(inv$novel[inv$novel$ITERATION=="2",colnames(nv2)], nv2)
        expect_equivalent(inv$novel[inv$novel$ITERATION=="3",colnames(nv3)], nv3)
        
        expect_equivalent(sort(names(inv$germline[["1"]])),sort(names(germdb2)))
        expect_equivalent(sort(names(inv$germline[["1"]])),sort(names(germdb2)))
        
    })
# }
