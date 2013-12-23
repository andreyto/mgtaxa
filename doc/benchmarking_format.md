Format of Benchmarking results
==============================
Introduction
------------
Benchmarking results have two dimensions:

1. Taxonomic exclusion rank. To simulate increasing degree of taxonomic 
   novelty, models from progressively higher ranks of sister clades of a 
   given test fragment are excluded from making the predictions. `no_rank` 
   means that no exclusion has been made.
2. Taxonomic prediction rank. For a given exclusion rank and a given test 
   fragment, the prediction is made, and the performance metrics are 
   evaluated at main taxonomic ranks in the predicted lineage. `no_rank` 
   means at the rank attached to the model. For RefSeq microbial genomes, 
   that would typically be a strain. In general, model can be attached to 
   any rank in the hierarchy. In that case, it will contribute to the 
   performance metrics only at that rank (in no-exclusion mode) and above.

Files
-----
There are multiple subdirectories corresponding to different test fragment 
lengths. 

In each subdirectory, there are multiple `.csv` tab delimited files. Prefix 
in file name marks the type of hierarchical placement of the performance 
metrics:

- `lin.`    "Main" taxonomic ranks, starting from species
- `bot.`    The lowest taxonomic node attached to the model. Note that 
  multiple models can be attached to the same taxonomic node, resulting at 
  aggregated performance metrics at this level.
- `mod.`    Model, regardless of the assigned taxonomy. This is the most 
  specific prediction level.

For each prefix above, the next file name component marks the type of 
performance metric:

 - `aggr`    The most high level metric. For each combination of exclusion 
   and prediction levels, the clade-level metric is reported averaged over 
   all clades. Separate tables for sensitivity, specificity and accuracy are 
   included, computed as in [PMID: 17179938]. Table titles are preceded with `#`.
 - `conf.raw`    Confusion table in a relational format. Each row corresponds 
   to a unique combination of exclusion level, prediction level, true ("test") 
   label and predicted label. `cnt` field contains a count of times this 
   combination has been observed. Label is expressed by two fields redundantly:
  - `id_lev_{test,pred}` Unique ID. For taxonomic labels, this is NCBI taxonomy 
    ID.
  - `name_{test,pred}`  Descriptive name. Might not be unique.

