# BioseqSuite
## Description
**BioseqSuite** is a tool designed to analyze any biological sequences i.e. DNA, RNA and proteins. It consists of 3 main functions: fastq_tool, run_dna_rna_tools and protein_tool.

### fastq_tool
**fastq_tool** is a tool designed to analize fastq sequences. More on this [here](https://knowledge.illumina.com/software/general/software-general-reference_material-list/000002211). It takes as input the dictionary, where **_key_** is the name of sequence (str) and **_value_** is a tuple of type (sequence, quality). This function returns the dictionary, but with filtered sequences according to the additional parameters that were set in the input. This function takes the next list of **parameters**:
- `seqs` - fastq sequences.
- `gc_bounts` - GC composition interval (in percent) for filtering (default is (0, 100)), i.e all reads passed. If you pass **a single number** to the argument, it is assumed to be an *upper* bound.
- `length_bounds` - filtering length interval (default is (0, 2**32)).
- `quality_threshold` - threshold value of average read quality for filtering, default is 0 (phred33 scale).
