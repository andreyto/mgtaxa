"""Classes that assign classification labels to various ranks in the taxonomy of our known sequence.
We adopt a hierarchical approach to classification in order to place a sequence onto a known
taxonomic tree. A recursive procedure is applied to each sample: if a node at some rank has
been assigned to it at a previous step, then a classifier trained for only sub-nodes of this
node at the next lower rank is applied to assign a sub-node. E.g. if a node is assigned
a given order, we apply a classifier trained to assign families under this order.
In addition, a one-class classifier is trained for each sub-node (family in our example).
If its prediction disagrees with the multi-class prediction, the sample is assigned to the
"reject" group.
The code in this module:
    - groups available training sequence under each node of each classification ranks
    - assigns classification labels based on ranks with addition of special label for
    sequences that do not have a required rank assigned to them (e.g. a sequence has
    order and genus but no family aasigned)
    - splits sequence into balanced testing and training sets, restricted by maximum
    and minimum sequence lengths and excluding from training all genera selected for testing
    - the testing sequence length in turn reflects the expected length distribution
    of prediction sequence fragments
    - after the testing sequence is set aside, the splitting procedure is again applied
    to training sequence in order to create datasets for n-fold parameter optimization.

Final output:
A set of tables that describe a hierarchy of classifiers and their training/testing sequence.
Structure of these tables allow straightforward generation of classification training/validation
jobs for each classifier independently, as well as testing/prediction with the entire hierarchical
algorithm.

Implementation:
We describe first a set of operations and then determine how to map them into
SQL or in-memory tree representation.
Terminology: in order to avoid confusion between taxonomic 'class' rank and
classification 'class', we replace the later with 'label' in the following text.

1. For each node at a given level (rank from our predefined set of classification
taxonomic ranks), select all next lower level sub-nodes (e.g. all 'family'
nodes under a given 'order' node). Direct subnodes with missing
next lower level rank form one separate group "unknown". Each group is assigned its
taxid as a class label, except for "unknown" group which is assigned label 1. 1 is an actual
taxid of NCBI root node, but that node will never appear as a classification target,
so it is safe to re-use it for arbitrary label. Label for the "other" node should be
made a parameter just in case.
Special case: it is possible that there are no subnodes with the next lower level
rank (e.g. we eliminated the only two existing 'order' nodes from the viral
hierarchy in order to deal with 'family' directly under viral 'superkingdom'.
And there are no 'phylum' or 'class' nodes for viruses.). The implementation
should use the next-next-lower level rank and proceed.
Special case: only one sub-node. We should replace it with its next level sub-nodes.

Discussion on the meaning of "unknown" group:
Two diferent origins are possible:
    a) The "real" classification hierarchy does not actually have a given rank defined
    for these sequences. Example are viral 'order' (only two orders are defined). Although
    we dealt with the viral case at the earlier stage by eliminating the "order" rank
    completely, such a case can be present at other branches of the tree. In that case
    it appears to be OK for us to lump all such sequences in a separate group, train
    to classify it along with "real" labels, and then apply the standard recursive procedure
    to that group. For example, if all viral families are separable from each other, then
    any non-intersecting sets of them must be also separable, and thus it is safe to
    classify first at the 'order' level (two defined orders and one "unknown" order
    with the rest of the families), and then at the 'family' level under each of the
    three 'order' labels.
    b) The required rank was not assigned to a given subnode for some reason (e.g. lack
    of information about its phenotype), but in reality it belongs to some of the defined
    nodes at this rank. E.g. for a given family node, no 'order' was defined (it is attached
    directly to some 'class' node), but given full information about the constituent
    species, it would be attached to one of the existing 'order' nodes. In such case creating
    a separate training label is wrong - such label will overlap with other labels and lead
    to misclassifications (given our central assumption that sequence composition reflects
    taxonomic hierarchy).
"""

def exclude_rRNA(self):
    db.createTableAs(
    """select * from seq_hdr where hdr rlike '[[:<:]]rRNA[[:>:]]'"""
    )
    