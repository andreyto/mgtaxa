#Introduction

`mgt-icm-classifier` is a unified driver program for ICM-based classification, model training and benchmarking.

In the following text, `<MGT_HOME>` means the directory where MGTAXA was installed.
When trying the provided examples of the commands, you should replace
`<MGT_HOME>` with the actual path.

It is assumed that your are executing the commands from under BASH shell.

Other examples instruct you to "source" MGTAXA environment file before executing
the commands: `source <MGT_HOME>/etc/mgtaxa.shrc`. After that, a set of
environment variables is defined in your current shell, including MGT_HOME
variable that can be referenced as `$MGT_HOME`. PATH is also modified so
that MGTAXA executables can be found without having to specify absolute paths.

**Your current working directory must be on a file system shared across the
nodes of your compute cluster if you are using a distributed cluster execution
mode.**

Run `<MGT_HOME>/bin/mgt-icm-classifier --help` to get a description of individual 
arguments.

#Some common use cases

##To classify with the default models pre-built from NCBI RefSeq

```
<MGT_HOME>/bin/mgt-icm-classifier --mode predict \
        --inp-seq contigs.fna \
        --inp-seq-attrib weights.csv --pred-min-len-samp 1000 \
        --pred-out-dir results --run-mode batchDep \
        --batch-backend makeflow --makeflow-options '-T sge' \
        --workflow-run 1 --lrm-user-options "-b n -P 0116"
```

In this example:

-   Multi-FASTA file contains **Nucleotide** sequences to assign taxonomy to
    (`--inp-seq contigs.fna`). It is required that the defline of your input 
    FASTA file contains a unique ID for each sequence. The ID is defined as 
    the part between the '>' and the first blank character. Only first 32
    characters of the ID will be used and should be unique among **all** input
    sequences.
    You can also use quoted shell glob pattern such as `'some_path/*.fna'`. 
    The pattern must be quoted so it is not expanded by the shell where you are
    running the command. The files might be compressed with gzip or bzip2,
    in which case they must have extensions `.gz` or `.bz2` correspondingly.

-   Distributed batch mode will be used for execution (`--run-mode batchDep`)
    You can switch to local sequential execution in a single process with
    `--run-mode inproc`. The default model library is too large to wait for
    completion of the run in `inproc` mode, but you can use this mode to quickly
    check that your parameters are accepted and then kill the run with
    `Ctrl-C`.

-   Makeflow workflow engine will be used to submit batch jobs in a proper
    order (`--batch-backend makeflow`)

-   The local resource manager (LRM, a.k.a. batch queuing system) used by Makeflow 
    will be SGE (`--makeflow-options '-T sge'`); notice the necessary quoting. 

-   Options specific to your local SGE installation are provided in 
    `--lrm-user-options` and should be quoted as shown above. These options are 
    those that you would need in order to submit some shell script for execution 
    with `qsub`. In the example above, `-b n` tells SGE not to treat the input
    script as an executable program; `-P 0116` is SGE project name required
    on some SGE installations.

-   The generated Makeflow workflow will be executed immediately 
    (`--workflow-run 1`). This means that the command will block and exit
    only after the processing has finished (or failed). This is what you
    typically want if you need to run only one `mgt-icm-classifier` command
    or a simple linear sequence of such commands in a shell script. Note
    that you can also submit this command itself for batch execution with
    the usual `qsub` mechanism, so that you would not have to run it on the submit
    node. This would require that the compute nodes in your cluster are allowed
    to submit new jobs, because this is what the master script would be doing
    when executing the workflow.

-   All input sequences with length less than 1000 bp (`--pred-min-len-samp 1000`) 
    will be ignored.

-   The optional weight file (`--inp-seq-attrib weights.csv`), if provided, will be 
    used when generating the aggregated clade abundance tables. If no weight file is
    provided, each sequence is given weight of 1. See `--help` for the file format.
    There is a sample script to generate the weight file from 454 Newbler assembly 
    output as the number of reads per contig:
    ```
    source <MGT_HOME>/etc/mgtaxa.shrc
    python $MGT_HOME/bin/454_contig_read_cnt.py < 454ReadStatus.txt > weights.csv
    ```

-   The results will be in 'results' directory (`--pred-out-dir results`).
    
    The results directory will contain the following files (default names are used
    here; names can be changed through additional command line options shown 
    by `--help`):

    -   Per-sequence predictions stored in a tab-delimited file [`pred-taxa.csv`].
        
        The fields are:
        -    id :             Sequence ID 
        -    len :            Length of sequence
        -    taxid :          NCBI taxonomy ID of the assigned model
        -    name :           NCBI node name corresponding to `taxid`
        -    rank :           NCBI name of the taxonomic rank for `taxid` 
                              (or `no_rank`)
        -    weight :         Weight as taken from the optional `--inp-seq-attrib`
        -    idscore :        ID of the model assigned to this sequence
        -    namescore :      Name of the model corresponding to idscore
        -    taxid_species :  NCBI taxonomy ID of the species rank in lineage
        -    name_species :   NCBI name of the species rank in lineage
        -    score_species :  Confidence score for the species assignment (0-1.)

        The last three fields are repeated for the higher order ranks in the 
        assigned taxonomic lineage using corresponding rank names as suffixes.

        For the ranks below the lowest assigned rank (from the `rank` field), the 
        lineage-specific fields will contain the value `null`.

        Note that we have to address a situation when NCBI lineage does not have
        a node defined for a specific rank (e.g. `order` follows `genus`, and 
        there is no `family`). We want to maintain a property that aggregating
        counts across different ranks results in the same total number for 
        sequences that were already assigned at some lower rank. Therefore, for 
        missing lineage ranks, we artifically fill them from the nodes at the 
        nearest defined rank either below or above.

    -   Summary counts for each assigned model [`pred-taxa.summ.csv`]. This is an 
        aggregation of `pred-taxa.csv` over `idscore` field, with `score_*` fields
        averaged, `weight` summed, and `cnt` field representing count. `id` field
        should be ignored.

    -   Aggregated counts at different ranks of the taxonomic hierarchy
        [`stats/stats.csv`]. There are multiple tables in one tab-delmited file.
        Tables are separated by `#`-prefixed comment header lines. `weight` and
        `len` fields are summed, `score` is averaged.

    -   Dynamic multi-level pie chart vizualization of the assigned taxonomic
        distribution implemented as a html input file [`stats/stats.html`] for 
        Krona (PMC3190407) JavaScript library. Krona library itself is also
        exported into the results directory, so that `stats.html` can be opened
        in a Web browser from the local file system using Open File dialog (Firefox 
        or Chrome, but not IE 8).

    -   PDF file with aggregated counts shown as bar plots [`stats/stats.pdf`]. 
        The deviation from a standard bar plot is that the **area** of each bar 
        is proportional to the count with both height and width scaled as square
        root of the count.

    -   Compact representation of per-sequence assignments is HDF5 format 
        [`pred_taxa`]

    -   Representation of summaries in SQLite format [`pred-taxa.sqlite`]

The application leaves various intermediate files (lots of them) in the starting
directory after it finishes. You might want to run the program from a temporary
working directory that you can delete wholesale afterwards.

### How to control and monitor Makeflow execution

Makeflow is an external software package (http://www3.nd.edu/~ccl/software/makeflow/) 
that is bundled with MGTAXA. In addition to its main workflow execution program 
`makeflow`, it provides utilities for monitoring the progress of the running workflow 
(`makeflow_monitor`) and to produce a report about a finished workflow 
(`makeflow_log_parser`).

To get access to these utilities, you can either source MGTAXA environment file
once in your current shell `source <MGT_HOME>/etc/mgtaxa.shrc` or prepend a call
to each executable with MGTAXA wrapper script, such as
`<MGT_HOME>/bin/mgt-wrapper makeflow --help`. These two options apply generally to
running any MGTAXA programs, with an exception of `mgt-icm-classifier`, which sources
the MGTAXA envirnoment file automatically. Below, we assume that you have
sourced the environment file, and are using Makeflow utilities without specifying the
full paths.

When user has selected `--batch-backend makeflow` to run computations, MGTAXA 
accepts options for `makeflow` via `--makeflow-options` value string that must be
either in single or double quotes if it contains spaces. MGTAXA will parse some 
of the `makeflow` options and/or set default values, and pass other options 
verbatim to `makeflow` without parsing. You can run `makeflow --help` to get the 
full list of options and also study the Makeflow manual. By default, MGTAXA tells 
Makeflow to resubmit each failed batch task up to 3 times with an interval of 30 seconds.
This provides a reasonable resilience against occasional compute node failures. 
The user can modify these `makeflow` options via `--makeflow-options` variable.
One exception is that `makeflow` option (`-B|--batch-options`) is ignored and
`--lrm-user-options` of `mgt-icm-classifier` program is used in its place.

Note that if you have a large multi-core machine that does not have any LRM
(such as SGE) installed, you can easily use it for the parallel execution of
`mgt-icm-classifier` pipeline. You would have to pass 
`--makeflow-options "-T local"` (or omit `--makeflow-options` altogether
because `-T local` is the default backend in Makeflow) and omit
`--lrm-user-options`. The Makeflow will then run locally, starting as many
tasks at the same time as there are cores on the machine.

It is also possible to run Makeflow with `-T mpi-queue` backend if your cluster's
LRM is configured to schedule efficiently only large MPI jobs. A thorough
example of using that mode can be found in our other package
[PGP](https://bitbucket.org/andreyto/proteogenomics).


Once you have started the `makeflow` execution, you can run 
`makeflow_monitor workflow.mkf.makeflowlog` from a Linux terminal to get
a dynamic view of the workflow progress (note that you need to supply a path to the
worflow log file, not to the workflow file itself).

You can also monitor your jobs or processes outside of Makeflow in the usual
ways (such as with `qstat` or `top`).

Once the workflow has finished, you can use 
`makeflow_log_parser workflow.mkf.makeflowlog` to get summary statistics on
the resources used for computations.

If the `makeflow` finishes successfully, it prints at the standard output
`Nothing else to do` and has exit code 0. Otherwise, it will print
`Makeflow failed` and have non-zero exit code. In that case the standard
output will also have messages about specific commands that failed. To
debug, you can copy-paste from the standard output the command line of 
the failed rule into a local terminal and run it. When executing the workflow 
under a batch system such as SGE, `makeflow` discards standard output and 
error streams of the jobs, so running the failed job commands locally is the most
direct way of getting the more descriptive error messages.

When `mgt-icm-classifier` is told to use Makeflow backend, it will create
a workflow file (default name is `workflow.mkf`) as well as a BASH script
(default name is `workflow.mkf.sh`) that can be used to properly execute 
`makeflow` with the workflow file as its input.

If `--workflow-run 1` is provided, the script will be executed immediately.
You often might want, however, to build higher-order pipelines out of
different invocations of `mgt-icm-classifier` or any other non-MGTAXA commands.
Although you can easily do this by writing shell scripts that call
`mgt-icm-classifier --workflow-run 1` for immediate execution, the most robust 
and scalable approach is to write these pipelines as Makeflow workflow files
as well.

To support that use case, `--workflow-run` option can be omitted or set to `0`.

You can look at hand-edited Makeflow file in the source tree directory that 
we use to train classification models from a local copy of NCBI RefSeq and 
benchmark them in a single invocation of `makeflow`:
`bin/onetime/build-ref-db/build.mkf`

The shell script that executes that workflow is: 
`bin/onetime/build-ref-db/run_build.sh`

The typical pattern in that Makeflow file looks like this:

```
REF_NCBI_ICM_DB.mkf: $REF_NCBI_SEQ_DB 
    $MGT_ICM_CLASSIFIER --mode train \
    --db-seq $REF_NCBI_SEQ_DB \
    --db-imm $REF_NCBI_ICM_DB \
    --run-mode batchDep \
    --batch-backend makeflow \
    --workflow-file REF_NCBI_ICM_DB.mkf

$REF_NCBI_ICM_DB: REF_NCBI_ICM_DB.mkf
    LOCAL MAKEFLOW REF_NCBI_ICM_DB.mkf
```

Above, the first make rule generates the workflow file `REF_NCBI_ICM_DB.mkf` but does 
not run it. The second rule specifies `REF_NCBI_ICM_DB.mkf` as its dependency, so
that it will run only after the file is generated. The command of the second rule
starts with keywords `LOCAL MAKEFLOW` which will result in the execution of
`REF_NCBI_ICM_DB.mkf` as a sub-Makeflow.

If any of the sub-Makeflows fails, you can address the reasons and launch
the script that runs the top-level workflow again - the Makeflow should pick
up the processing after the last step that succeeded. Note that if you have
edited the master workflow file (`build.mkf`), you typically need to delete
the corresponding Makeflow log file (`build.mkf.makeflowlog`), which will
result in the workflow re-executing from the very beginning. Makeflow uses
the log to figure out where to continue processing after a restart. If the
signatures of the rules have changed, the log is no longer valid and likely
to confuse the Makeflow.

##To train new models and use them for classification

Training of the models consists of two steps:

###Build a sequence database in internal representation from input FASTA files

```
<MGT_HOME>/bin/mgt-icm-classifier --mode make-ref-seqdb \
    --inp-train-seq models.fasta.gz \
    --inp-train-seq-format generic \
    --inp-train-model-descr models.json \
    --db-seq seq.db \
    --run-mode batchDep \
    --batch-backend makeflow --makeflow-options '-T sge' \
    --workflow-run 1 --lrm-user-options "-b n -P 0116"
```

In this example:

-   Multi-FASTA file contains **Nucleotide** sequences to use for training the
    models (`--inp-train-seq models.fasta.gz`). Other requirements described 
    above for the `--inp-seq` option also apply here. In addition, an option
    `--inp-train-seq-list` can be used to specify a text file that contains
    a list of files or shell glob patterns (see `--help`). `--inp-train-seq`
    and `--inp-train-seq-list` can be used both at the same time.

-   The training sequences are of a "generic" origin (`--inp-train-seq-format generic`).
    This is opposed to `--inp-train-seq-format ncbi`, in which case the sequences
    would be assumed to come from the NCBI databases and have NCBI GI ID in the 
    defline, which MGTAXA would then try to match with the taxonomic node using
    NCBI-provided GI-to-taxa index.

-   When `--inp-train-seq-format generic` is used, a JSON file has to be provided
    that groups sequence IDs from `--inp-train-seq models.fasta.gz` into models and
    maps the models to the NCBI taxonomy.
    
    An example of such file is:
    
    ```
    [
        {
            "ids_seq": ["contig1","contig2","contig3"], 
            "id": "mod1", 
            "name":"My model 1", 
            "taxid": 766
        }, 
        {
            "ids_seq": ["contig4","contig5"], 
            "id": "mod2", 
            "name": "My model 2", 
            "taxid": 767892
        }
    ]
    ```

    Above:

    -   `ids_seq` contains a list of sequence IDs from `models.fasta.gz`. The
        corresponding model will be trained on the union of these sequences.
        Sequences in `ids_seq` list do not have to be adjacent in the `models.fasta.gz`
        file.
    -   `id` is a unique ID you wish to give to your model (32 characters
        long or less, no spaces)
    -   `name` is a descriptive name you wish to give to your model (may not
        be unique but it helps if it is)
    -   `taxid` is an existing taxonomic ID from the 
        [NCBI Taxonomy DB](http://www.ncbi.nlm.nih.gov/taxonomy). It is entirely
        up to you how you map your training sequences to the NCBI taxonomy.
        For example, you might have establsihed from a protein-level annotation
        that your contigs 1 to 3 belong to organism(s) from the order Rickettsiales
        (taxid 766). You then make that order-level assignment to your model.
        If you are not interested in the taxonomy of your model at all and only
        going to use it for recruiting other sequences, you can set taxid to 1
        (root node of the taxonomy).

-   The output path of this command is given by (`--db-seq seq.db`). You will need 
    to use that path when you run model training.

###Train models based on a sequence DB built by a previous `--mode make-ref-seqdb` call.

```
<MGT_HOME>/bin/mgt-icm-classifier --mode train \
    --db-seq seq.db \
    --db-imm imm.db \
    --train-min-len-samp 5000 \
    --run-mode batchDep \
    --batch-backend makeflow --makeflow-options '-T sge' \
    --workflow-run 1 --lrm-user-options "-b n -P 0116"
```

In this example:

-   Database of model sequences built before is used as input (`--db-seq seq.db`)

-   Models will be placed at the output path (`--db-imm imm.db`). That path you
    will need to use when classifying agains these models.

-   Any models that have a total length of reference sequence less than
    `--train-min-len-samp 5000` will be skipped during training. The default value
    is 2000. It is not recommended to go below that.

Once you have trained your models, you can use them to classify other sequences
with or without combining your new set of models with other sets.

###Classify using a combination of our models and default models.

```
source <MGT_HOME>/etc/mgtaxa.shrc
mgt-icm-classifier --mode predict \
    --db-imm imm.db \
    --db-imm $MGT_DATA/refseq-icm \
    --db-imm $MGT_DATA/ref-extra-icm \
    --inp-seq contigs.fna \
    --inp-seq-attrib weights.csv --pred-min-len-samp 1000 \
    --pred-out-dir results --run-mode batchDep \
    --batch-backend makeflow --makeflow-options '-T sge' \
    --workflow-run 1 --lrm-user-options "-b n -P 0116"
```

In this example:

-   Compared to the first example of classifying sequences that used a default set
    of models, we classify against a combination of RefSeq-based models
    (`--db-imm $MGT_DATA/refseq-icm`), models built for a few other genomes
    such as SAR86 (`--db-imm $MGT_DATA/ref-extra-icm`) and models that we
    have built ourselves in the previous example (`--db-imm imm.db`).

-   We dereference the environment variable `$MGT_DATA` to refer to the location of
    the default models. That variable is defined in the current shell when the 
    command `source <MGT_HOME>/etc/mgtaxa.shrc` is executed.

