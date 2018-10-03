Version 0.3.0 October 3, 2018
-------------------------------------------------------------------------------

Bug Fixes:

+ Fixed a bug in `reassignAlleles` occuring with single match genotypes.
+ Fixed `selectNovel` improperly removing all identical novel alleles, rather 
  than keeping a single entry.
+ `genotypeFasta` will now retain IMGT-numbering spacers as `.` characters
  instead of converting them to `-` characters.
+ Fixed a bug in `findNovelAlleles` causing overly aggressive minimum sequence 
  threshold filtering.
+ Fixed a bug in the grouping behavior of `getPopularMutationCount`.
  
New Features:

+ Added a Bayesian approach to genotype inferrence as the 
  `inferGenotypeBayesian` function.
+ Added the function `generateEvidence` to build a complete evidence table
  from the results of `findNovelAlleles`, `inferGenotype`, 
  `inferGenotypeBayesian`, and `reassignAlleles`.
+ Added multiple new evidence columns to the output of `findNovelAlleles`
  and adjusted the definitions/names of some existing columns.
+ Added behavior to the `keep_gene` argument of `reassignAlleles` to provide
  options for maintaining reassignments at the gene (previous `TRUE` behavior), 
  family, or repertoire level.
+ Improved tie resolution in `findNovelAlleles`.
  
Backwards Incompatible Refactors:

+ Renamed sample data from `germline_ighv`, `sample_db`, `genotype` and 
  `novel_df` to `GermlineIGHV`, `SampleDb`, `SampleGenotype` and `SampleNovel`,
  respectively.
+ Renamed the `novel_df` argument to `novel` in `selectNovel`, `inferGenotype`,
  and `genotypeFasta`.
+ Renamed the `novel_df_row` argument to `novel_row` in `plotNovel`.
+ Argument order in `inferGenotype` was alter for clarity.
+ Changed the return behavior of `reassignAlleles` so that it returns the 
  input data.frame with the `V_CALL_GENOTYPED` column appended or overwritten.
+ `cleanSeqs` will no longer replace `.` characters with `-`.


Version 0.2.11 September 21, 2017
-------------------------------------------------------------------------------

+ Improved memory utilization in `findNovelAlleles`.

  
Version 0.2.10 July 1, 2017
-------------------------------------------------------------------------------

+ Bugfix wherein `inferGenotype` would break when performing check for alleles
  that could not be distinguished.


Version 0.2.9.999 May 16, 2017
-------------------------------------------------------------------------------

+ Bugfix wherein `inferGenotype` would break if all sequences submitted were
  from a single gene and `find_unmutated` was set to `TRUE`.


Version 0.2.9: March 24, 2017
-------------------------------------------------------------------------------

+ License changed to Creative Commons Attribution-ShareAlike 4.0 International
(CC BY-SA 4.0).


Version 0.2.8: August 26, 2016
-------------------------------------------------------------------------------

+ Bugfix following recent update of alakazam (0.2.5) to import selectively.
+ Removed unneeded dependency on shazam package (not needed as of 0.2.5.999).


Version 0.2.7:  July 24, 2016
-------------------------------------------------------------------------------

+ More updates to work with the latest version of dplyr (0.5.0).
+ Bugfix in findNovelAlleles when allele passed germline_min but not min_seqs.
+ Fixed vignette typo and updated findUnmutatedCalls man page.


Version 0.2.6:  July 01, 2016
-------------------------------------------------------------------------------

+ Updated code to work with the latest version of dplyr (0.5.0).


Version 0.2.5.999:  June 10, 2016
-------------------------------------------------------------------------------

+ Fixed a bug werein `findNovelAlleles()` was not running in parallel, even 
  when `nproc` > 1.
+ Changed default to `nproc=1` in `findNovelAlleles()`.

Version 0.2.5:  June 07, 2016
-------------------------------------------------------------------------------

+ Initial CRAN release.
