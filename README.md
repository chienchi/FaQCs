illuminaFastqQC: Quality Control of Next Generateion Sequenceing Data
=======

-------------
PREREQUISITES
-------------

1. The main program is developed in Perl v 5.8.8.
2. Parallel::ForkManager module from CPAN   
   (http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.9/lib/Parallel/ForkManager.pm)
3. R for ploting                 
   (http://www.r-project.org/)
                             
4. Jellyfish for kmer counting   (Optional) 
   (http://www.cbcb.umd.edu/software/jellyfish/) 

-----------
BASIC USAGE
-----------

* Trimming by quality 5 and filtering reads with any ambiguous base or low complexity.
perl illumina_fastq_QC.pl -p 'reads1.fastq reads2.fastq' -d out_directory

* Quailty check only on subsamples of input, no trimming and filtering. 
perl illumina_fastq_QC.pl -p 'reads1.fastq reads2.fastq' -d out_directory -qc_only 

---------------
VERSION HISTORY
---------------

======== Version 1.1
New features and changes in illumina_fastq_qc  version 1.1 with respect to version 1.0:
- add -qc_only option for quick quality check without trimming and filtering
- add -discard option to output discarded reads

======== Version 1.0
Stable function release.
Features:
- trim bidirection
- minimium length filtering after trim
- "N" base filter
- low complexity filter
- average read quality filter
- autocheck quality encoding and quality encoding coversion
- multi-threads  (required Parallel::ForkManager)
- input paired end reads aware

