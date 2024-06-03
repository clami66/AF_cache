# Summary of changes

I made some changes to your code (and AF) to make pairing work better (all data in /proj/berzelius-2021-29/users/x_arnel/herpes/bin). I try to summarize (and also for myself to remember)

## format_alphafold_jobs_AE.py
* Added option to only run pairs in a list --pairs and some other options

## taxid2virusid.py
* A script that I use to modify the pickle file to replace virusIDs with hostIDs (sometimes multiple ones)

## mmseqs_2_uniprot_parallel_test.py
* I added so that also sequences starting with "UniRef100_" are included I added some specie information in the saved names in the a3m files (so that I could trace how paring worked)
* I added one option "--shortname" so that all sequence names are changed to 7 letters otherwise they do not pass the pairing in AF2M
* I added an option "--duplicate" so that identical sequences with multiple specie-IDs are included the output. (this is quite special for my problem so probably not a good idea to include in a default version).

## Changes to AlphaFold (in ../apps/)
* Added an loggin of number of paired sequences msa_pairing.py
* Added an alternative RexExp patters in msa_identifiers.py (not necessary if using short sequences)
* Added a flag "includeidentical" in run_alphafold.py This is passed on to pipeline.py / make_msa_features and produces a warning for duplicated sequences.
* Various other debug outputs

## getnumpairs.sh
* Extraxt number of pairs used in pairing (also exists in the log-files, but rather quick to do anyhow.)