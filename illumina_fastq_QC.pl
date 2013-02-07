#!/usr/bin/env perl
######   Required  ########################################################
#          1. Parallel::ForkManager module from CPAN                      #
#          2. String::Approx  module from CPAN                            #
#          3. R for ploting                                               #
#          4. Jellyfish for kmer counting                                 #
#             (http://www.cbcb.umd.edu/software/jellyfish/)               #
###########################################################################
######   Inputs   #########################################################
#          1. Fastq files either paired-end or unpaired reads or both     #
#              Can input multiple library fastq files but only output     #
#              concatenate trimmed fastq files                            #
#          2. Output directory                                            #
#          3. Other options                                               #
###########################################################################
######   Output   #########################################################
#          1.  Two Paired-ends files if input paired-end reads            #
#          2.  One unpaired reads file                                    #
#          3.  trimming statistical text file                             #
#          4.  quality report pdf file                                    #
###########################################################################
######   Functions     ####################################################
#       1. trim bidirection. 2. add -min_L flag                           #
#              3. phredScore= ord(Q)-$ascii for diffent quality           #
#                 encoding and conversion                                 #
#              4. -n "N" base filter   5. Low complexity filter           #
#              6. multi-threads  (required Parallel::ForkManager)         #
#              7. output ascii conversion 8. stats report                 #
#              9. input paired end reads 10. average read quality filter  #
#       11. Filter artifact                                               #
# AUTHOR: CHIEN-CHI LO                                                    #                                                
# Copyright (c) 2012 LANS, LLC All rights reserved                        #
# All right reserved. This program is free software; you can redistribute #
# it and/or modify it under the same terms as Perl itself.                #
#                                                                         #
# LAST REVISED: September 2012                                            #
###########################################################################
use strict;
use File::Basename;
use Getopt::Long;
use Parallel::ForkManager;
use String::Approx;

my $version=1.2;
my $debug=0;

sub Usage {
print <<"END";
     Usage: perl $0 [options] [-u unpaired.fastq] -p 'reads1.fastq reads2.fastq' -d out_directory
     Version $version
     Input File: (can use more than once)
            -u            Unpaired reads
            
            -p            Paired reads in two files and separate by space in quote
     Trim:
            -q            Targets # as quality level (default 5) for trimming
                          NOT A HARD cutoff! (see bwa's bwa_trim_read() function in bwaseqio.c)
     Filters:
            -min_L        Trimmed sequence length will have at least minimum length (default:50)

            -avg_q        Average quality cutoff (default:0, no filtering)
            
            -n            <INT> Filter Sequence has more than this number of base "N" will be discarded. 
                          (default: 0, zero tolerance)

            -lc           Low complexity filter ratio, Maximum fraction of mono-/di-nucleotide sequence  (default: 0.85)
            
            -adapter      <bool> Filter reads with illumina adapter/primers (default: no)
                          -rate   Mismatch ratio of adapters' length (default: 0.2, allow 20% mismatches)
            					
            -artifactFile  <file>    additional artifact (adapters/primers/contaminations) reference file in fasta format 
     Q_Format:
            -ascii        Encoding type: 33 or 64 or autoCheck (default)
                          Type of ASCII encoding: 33 (standard) or 64 (illumina 1.3+)

            -out_ascii    Output encoding. (default: 33)
     Output:
            -prefix       Output file prefix. (default: QC)

            -stats        Statistical numbers output file (default: prefix.stats.txt)

            -d            Output directory.
     Options:
            -t            # of CPUs to run the script (default:2 )

            -split_size   Split the input file into several sub files by sequence number (default: 1000000) 

            -qc_only      <bool> no Filters (-q 0 -avg_q 0 -n 1000 -lc 1 -min_L 0)

            -kmer_rarefaction_on     <bool>   
                          Turn on the kmer calculation. Turn on will slow down ~10 times. (default:Calculation is off.)
                          (meaningless if -subset is too small)
                          -m  <INT>     kmer for rarefaction curve (range:[2,31], default 31)

            -subset       <INT>   Use this nubmer x split_size for qc_only and kmer_rarefaction  
                                  (default: 10,  10x1000000 SE reads, 20x1000000 PE reads)

            -discard      <bool> Output discarded reads to prefix.discard.fastq (default: 0, not output)

            -debug        keep intermediate files
END
exit;
}

# magic number of quality score
my $highest_illumina_score=41;
my $lowest_illumina_score=0;
# Options Variable initialization
my $thread=2;
my $opt_q=5;
my $opt_min_L=50;
my $opt_avg_cutoff=0;
my $ascii;
my $N_num_cutoff=0;
my $out_offset=33;
my $low_complexity_cutoff_ratio=0.85;
my $subfile_size=1000000;
#my $subfile_size=100000;
my $kmer_rarefaction_on=0; 
my $kmer=31; 
my $prefix="QC";
my $plots_file;
my $stats_output;
my $trimmed_reads1_fastq_file;
my $trimmed_reads2_fastq_file;
my $trimmed_unpaired_fastq_file;
my $trimmed_discard_fastq_file;
my @paired_files;
my @unpaired_files;
my $outDir;
my $output_discard;
my $qc_only;
my $stringent_cutoff=0;
my $filter_adapter;
my $filterAdapterMismatchRate=0.2;
my $artifactFile;
my $subsample_num=10;

# Options
GetOptions("q=i"          => \$opt_q,
           "min_L=i"      => \$opt_min_L,
           "avg_q=f"      => \$opt_avg_cutoff,
           "p=s"          => \@paired_files,
           "u=s"          => \@unpaired_files,
           "ascii=i"      => \$ascii,
           "n=i"          => \$N_num_cutoff,
           'stringent_q=i'=> \$stringent_cutoff,
           "lc=f"         => \$low_complexity_cutoff_ratio,
           "out_ascii=i"  => \$out_offset,
           "t|threads=i"  => \$thread,
           'kmer_rarefaction_on' => \$kmer_rarefaction_on,
           'm=i'          => \$kmer,
           'split_size=i' => \$subfile_size,
           'prefix=s'     => \$prefix,
           'd=s'          => \$outDir,
           'stats=s'      => \$stats_output,
           'discard'      => \$output_discard,
           'qc_only'      => \$qc_only,
           'subset=i'     => \$subsample_num,
           'debug'        => \$debug,
           'adapter'      => \$filter_adapter,
           'rate=f'       => \$filterAdapterMismatchRate,
           'artifactFile=s'  => \$artifactFile,
           'R1=s'         => \$trimmed_reads1_fastq_file,   # for galaxy impelementation
           'R2=s'         => \$trimmed_reads2_fastq_file,   # for galaxy impelementation
           'Ru=s'         => \$trimmed_unpaired_fastq_file, # for galaxy impelementation
           'Rd=s'         => \$trimmed_discard_fastq_file,  # for galaxy impelementation
           'QRpdf=s'      => \$plots_file,                  # for galaxy impelementation
           "help|?"       => sub{Usage()} );

die "Missing input files.\n", &Usage() unless @unpaired_files or @paired_files;
die "Missing output directory.\n",&Usage() unless $outDir;

###### Output file initialization #####
# temp files for plotting
my $quality_matrix="$outDir/$prefix.quality.matrix";
my $avg_quality_histogram="$outDir/$prefix.for_qual_histogram.txt";
my $base_matrix="$outDir/$prefix.base.matrix";
my $nuc_composition_file="$outDir/$prefix.base_content.txt";
my $length_histogram="$outDir/$prefix.length_count.txt";

# output files
$trimmed_discard_fastq_file="$outDir/$prefix.discard.fastq" if (!$trimmed_discard_fastq_file);
$plots_file="$outDir/${prefix}_qc_report.pdf" if (!$plots_file);
$trimmed_unpaired_fastq_file="$outDir/$prefix.unpaired.trimmed.fastq" if (!$trimmed_unpaired_fastq_file);
$trimmed_reads1_fastq_file="$outDir/$prefix.1.trimmed.fastq" if (!$trimmed_reads1_fastq_file);
$trimmed_reads2_fastq_file="$outDir/$prefix.2.trimmed.fastq" if (!$trimmed_reads2_fastq_file);
$stats_output="$outDir/$prefix.stats.txt" if (!$stats_output);
   
# temp files for kmer counting plot   
my $kmer_rarefaction_file="$outDir/$prefix.Kmercount.txt"; #
my $kmer_files="$outDir/$prefix.KmerFiles.txt"; #
my $kmer_histogram_file="$outDir/$prefix.kmerH.txt"; #
#######################################

######   Output check  ################
if (! -e $outDir)
{
  mkdir $outDir;
}

if (-e $trimmed_reads1_fastq_file)
{
   print "The output $trimmed_reads1_fastq_file file exists and will be overwritten.\n";
   system ("rm $trimmed_reads1_fastq_file");
}
if (-e $trimmed_reads2_fastq_file)
{
   print "The output $trimmed_reads2_fastq_file file exists and will be overwritten.\n";
   system ("rm $trimmed_reads2_fastq_file");
}
if (-e $trimmed_unpaired_fastq_file)
{
   print "The output $trimmed_unpaired_fastq_file file exists and will be overwritten.\n";
   system ("rm $trimmed_unpaired_fastq_file");
}
if (-e $plots_file)
{
   print "The output $plots_file file exists and will be overwritten.\n";
   system ("rm $plots_file");
}
if (-e $trimmed_discard_fastq_file)
{
   print "The output $trimmed_discard_fastq_file file exists and will be overwritten.\n";
   system ("rm $trimmed_discard_fastq_file");
}
system ("rm $kmer_files") if (-e $kmer_files);
#######################################

if ($qc_only)
{
    # set all Filters to pass
    $opt_q=0;
    $opt_min_L=0;
    $opt_avg_cutoff=0;
    $N_num_cutoff=1000;
    $low_complexity_cutoff_ratio=1;
}


my ( $total_count,$total_count_1, $total_count_2, $total_num, $trimmed_num,$total_raw_seq_len, $total_trimmed_seq_len);
my ( $trim_seq_len_std, $trim_seq_len_avg, $max, $min, $median, $numOfReadsWithAdapter);
my ( $paired_seq_num, $total_paired_bases );
my ( @split_files, @split_files_2);
my %position;
my %AverageQ;
my %base_position;
my %base_content;
my %len_hash;
my ( $i_file_name, $i_path, $i_suffix );
  
  foreach my $input (@unpaired_files,@paired_files){
     print "Processing $input file\n";
     #print $STATS_fh "Processing $input file\n";
     my ($reads1_file,$reads2_file) = split /\s+/,$input;
  
     # check file 
     if(&file_check($reads1_file)<0) { die "The file $reads1_file doesn't exist or empty.\n";}
     if(&file_check($reads2_file)<0 and $reads2_file) { die "The file $reads2_file doesn't exist or empty.\n";}

     # check quality offset
     if (! $ascii){$ascii = &checkQualityFormat($reads1_file)}

    #split
    ($total_count_1,@split_files) = &split_fastq($reads1_file,$outDir,$subfile_size);
    ($total_count_2,@split_files_2) = &split_fastq($reads2_file,$outDir,$subfile_size) if ($reads2_file);
     $total_count += $total_count_1 + $total_count_2;
     my $random_num_ref = &random_subsample(scalar(@split_files),$subsample_num);
     $subsample_num -= scalar(@split_files);
     $subsample_num -= scalar(@split_files) if ($reads2_file);

    my $pm = new Parallel::ForkManager($thread);
 # data structure retrieval and handling from multi-threading 
   # $pm->run_on_start(
   #   sub { my ($pid,$ident)=@_;
   #    #print "** $ident started, pid: $pid\n";
   #   }
   # );

    $pm -> run_on_finish ( # called BEFORE the first call to start()
      sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $nums_ref) = @_;
      
        # retrieve data structure from child
        if (defined($nums_ref)) {  # children are not forced to send anything
          my ($total_score, $avg_score);
          $total_num += $nums_ref->{raw_seq_num};
          $trimmed_num += $nums_ref->{trim_seq_num};
          $total_raw_seq_len += $nums_ref->{total_raw_seq_len};
          $total_trimmed_seq_len += $nums_ref->{total_trim_seq_len};
          my $processed_num= $nums_ref->{raw_seq_num};
          $trim_seq_len_std = $nums_ref->{trim_seq_len_std};
          $trim_seq_len_avg = $nums_ref->{trim_seq_len_avg};
          $max = $nums_ref->{max};
          $min = $nums_ref->{min};
          $median = $nums_ref->{median};
          $paired_seq_num += $nums_ref->{paired_seq_num};
          $total_paired_bases +=  $nums_ref->{total_paired_bases};
          $numOfReadsWithAdapter += $nums_ref->{numOfReadsWithAdapter};
          my %temp_avgQ = %{$nums_ref->{ReadAvgQ}};
          map {$AverageQ{$_}->{bases} += $temp_avgQ{$_}->{basesNum};
               $AverageQ{$_}->{reads} += $temp_avgQ{$_}->{readsNum}; 
              } keys %temp_avgQ; 
          my %temp_position= %{$nums_ref->{qual}};
          my %temp_base_position= %{$nums_ref->{base}};
          foreach my $pos (1..($total_raw_seq_len/$total_num))
          {
                for my $score ($lowest_illumina_score..$highest_illumina_score)
                { 
                    $position{$pos}->{$score} += $temp_position{$pos}->{$score};
                    $total_score +=  $score * $temp_position{$pos}->{$score}; 
                }
                for my $nuc ("A","T","C","G","N")
                {
                    $base_position{$pos}->{$nuc} += $temp_base_position{$pos}->{$nuc};
                }
          }
          my %tmp_base_content = %{$nums_ref->{Base_content}};
          for my $nuc ("A","T","C","G","N","GC")
          {
              while (my ($key, $value)= each %{$tmp_base_content{$nuc}})
              {
                 $base_content{$nuc}->{$key} += $value;
              }
          } 
          while (my ($key, $value)= each %{$nums_ref->{ReadLen}} )
          {
		         $len_hash{$key} += $value;   
          }
          if ( $random_num_ref->{$ident}){
              open (KMEROUT, ">>$kmer_files") or die "$!\n";
              print KMEROUT $nums_ref->{trim_file_name_1},"\t",$nums_ref->{trim_seq_num_1},"\n";
              print KMEROUT $nums_ref->{trim_file_name_2},"\t" if ($nums_ref->{trim_file_name_2});
              print KMEROUT $nums_ref->{trim_seq_num_2},"\n" if ($nums_ref->{trim_file_name_2});
              close KMEROUT; 
          }

          #print $STATS_fh " Processed $total_num/$total_count\n";
          print "Processed $total_num/$total_count\n";
          printf (" Post Trimming Length(Mean, Std, Median, Max, Min) of %d reads with Overall quality %.2f\n",$processed_num, $total_score/$nums_ref->{total_trim_seq_len});
          printf (" (%.2f, %.2f, %.1f, %d, %d)\n",$trim_seq_len_avg,$trim_seq_len_std,$median,$max,$min);
          #unlink $split_files[$ident];
          #unlink $split_files_2[$ident] if ($split_files_2[$ident]);
         
        } else {  # problems occuring during storage or retrieval will throw a warning
          print qq|No message received from child process $pid! on $ident\n|;
        }
      }
    );

  foreach my $i(0..$#split_files)
  {
    next if ($qc_only and ! $random_num_ref->{$i});
    $pm->start($i) and next;
    my $hash_ref = &qc_process($split_files[$i],$split_files_2[$i]);
    &run_kmercount($hash_ref,$kmer,1) if ( $random_num_ref->{$i} and $kmer_rarefaction_on );
    &run_kmercount($hash_ref,$kmer,2) if ( $random_num_ref->{$i} and $kmer_rarefaction_on and $split_files_2[$i]);
    $pm->finish(0, $hash_ref);
   
  }
  $pm->wait_all_children;
  # clean up
  foreach my $i(0..$#split_files)
  {
    unlink $split_files[$i];
    unlink $split_files_2[$i] if ($split_files_2[$i]);
  }
} #end foreach $input


# concatenate each thread's trimmed reads and files clean up.
if (! $qc_only)
{
    if (@unpaired_files) 
    {
      foreach my $input (@unpaired_files){
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$input", qr/\.[^.]*/ );
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$i_file_name", qr/\.[^.]*/ ) if ($i_suffix =~ /gz$/);
    
        if (system("cat $outDir/${i_file_name}_?????_trim.fastq >> $trimmed_unpaired_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim.fastq`;
        if ($output_discard){
            if (system("cat $outDir/${i_file_name}_?????_trim_discard.fastq >> $trimmed_discard_fastq_file")) { die "cat failed: $!" }
            `rm $outDir/${i_file_name}_?????_trim_discard.fastq`;
        }
      }
    }
    
    if (@paired_files)
    {
      foreach my $input (@paired_files){
        my ($reads1_file,$reads2_file) = split /\s+/,$input;
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$reads1_file", qr/\.[^.]*/ );
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$i_file_name", qr/\.[^.]*/ ) if ($i_suffix =~ /gz$/);
        if (system("cat $outDir/${i_file_name}_?????_trim.fastq >> $trimmed_reads1_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim.fastq`;
        if ($output_discard){ 
           if (system("cat $outDir/${i_file_name}_?????_trim_discard.fastq >> $trimmed_discard_fastq_file")) { die "cat failed: $!" }
           `rm $outDir/${i_file_name}_?????_trim_discard.fastq`;
        }
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$reads2_file", qr/\.[^.]*/ );
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$i_file_name", qr/\.[^.]*/ ) if ($i_suffix =~ /gz$/);
        if (system("cat $outDir/${i_file_name}_?????_trim.fastq >> $trimmed_reads2_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim.fastq`;
        if (system("cat $outDir/${i_file_name}_?????_trim_unpaired.fastq >> $trimmed_unpaired_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim_unpaired.fastq`;
      }
    }  
}

&Kmer_rarefaction() if ($kmer_rarefaction_on); 
&print_quality_report_files();
&print_final_stats();
&plot_by_R();
unless ($debug){
unlink $nuc_composition_file;
unlink $quality_matrix;
unlink $base_matrix;
unlink $avg_quality_histogram;
unlink $length_histogram;
unlink $kmer_files;
    if ($kmer_rarefaction_on)
    {
#        unlink $kmer_rarefaction_file;
#        unlink $kmer_histogram_file;
    }
}
exit(0);

sub run_kmercount 
{
    my $hash_ref =shift;
    my $kmer =shift;
    my $mate = shift;
    my $file= $hash_ref->{"trim_file_name_$mate"};
    my $cmd = "jellyfish count -C -s 512M -c 3 -m $kmer -o $file$prefix $file";
    if (system($cmd)){
      die "$!\n";
    }
    
}  


sub Kmer_rarefaction {                                                               
    my $old_merge;                                                                     
    my $count=0;                                                                       
    my $new_merge="$outDir/${prefix}tmpKmerMerge${$}$count";                           
    my (@kmercount,$distinct_kmer,$total_kmer);                          
    #opendir (DIR, $outDir);                                                            
    #my @files=grep { /${prefix}_0/ } readdir(DIR);                                     
    #closedir DIR;  
    open (KMER, ">$kmer_rarefaction_file") or die "$!\n";                              
    open (FILES, "$kmer_files") or die "$!\n";
    while (<FILES>)
    #for my $i(0..$#files)                                                              
    {                                          
      chomp;                                        
      #my $subset_kmer_file="$outDir/$files[$i]";                         
      my ($subset_kmer_file, $seq_num) = split /\t/, $_; 
      $subset_kmer_file = $subset_kmer_file.$prefix."_0";
      #my ($subset_file_name) = $subset_kmer_file =~ /(\S+)${prefix}_0/; 
      #$sequence_num += $subfile_size_hash{$subset_file_name};                          
      if ($count==0)                                                                   
      {                                                                                
         @kmercount=`jellyfish stats $subset_kmer_file | grep -A 1 'Distinct' | awk '{print \$2}'`;                                                                         
         ($distinct_kmer) = $kmercount[0] =~ /(\d+)/;                                  
         ($total_kmer) = $kmercount[1] =~ /(\d+)/;                                     
         print KMER $seq_num . "\t". $distinct_kmer. "\t". $total_kmer."\n";      
         $old_merge=$subset_kmer_file;                                                 
         $count++;                                                                     
      }                                                                                
      else                                                                             
      {                                                                                
         `jellyfish merge -o $new_merge $old_merge $subset_kmer_file `;                
         @kmercount=`jellyfish stats $new_merge | grep -A 1 'Distinct' | awk '{print \$2}'`;                                                                                
         ($distinct_kmer) = $kmercount[0] =~ /(\d+)/;                                  
         ($total_kmer) = $kmercount[1] =~ /(\d+)/;                                     
         unlink $old_merge;                                                            
         unlink $subset_kmer_file;                                                     
         $old_merge=$new_merge;                                                        
         $new_merge="$outDir/${prefix}tmpKmerMerge${$}$count";                         
         print KMER $seq_num . "\t". $distinct_kmer. "\t". $total_kmer."\n";      
         $count++;                                                                     
      }                                                                                
    }                                     
    close FILES;                                             
    close KMER;                                                                        
    `jellyfish histo -o $kmer_histogram_file $old_merge`;                              
    unlink $old_merge;                                                                 
} 

sub print_final_stats{
    open (my $fh, ">$stats_output") or die "$!\t$stats_output\n";
    # stats
    print $fh "\nBefore Trimming\n";
    print $fh "Reads #:\t$total_num\n";
    print $fh "Total bases:\t$total_raw_seq_len\n";
    printf $fh ("Reads Length:\t%.2f\n",$total_raw_seq_len/$total_num);
    
    print $fh "\nAfter Trimming\n";
    printf $fh ("Reads #:\t\%d (%.2f %%)\n",$trimmed_num, $trimmed_num/$total_num*100);
    printf $fh ("Total bases:\t\%d (%.2f %%)\n",$total_trimmed_seq_len,$total_trimmed_seq_len/$total_raw_seq_len*100);
    if ($trimmed_num)
    {
      printf $fh ("Mean Reads Length:\t%.2f\n",$total_trimmed_seq_len/$trimmed_num); 
    }
    else
    {
      printf $fh "Mean Reads Length:\t0\n";
    }
    
    if (@paired_files){
      printf $fh ("  Paired Reads #:\t\%d (%.2f %%)\n",$paired_seq_num, $paired_seq_num/$trimmed_num*100);
      printf $fh ("  Paired total bases:\t\%d (%.2f %%)\n",$total_paired_bases,$total_paired_bases/$total_trimmed_seq_len*100);
      printf $fh ("  Unpaired Reads #:\t\%d (%.2f %%)\n", $trimmed_num - $paired_seq_num, ($trimmed_num - $paired_seq_num)/$trimmed_num*100);
      printf $fh ("  Unpaired total bases:\t\%d (%.2f %%)\n", $total_trimmed_seq_len - $total_paired_bases , ($total_trimmed_seq_len - $total_paired_bases)/$total_trimmed_seq_len*100);
    }
    
    if ($filter_adapter){
        printf $fh ("\nReads with Adapters/Primers #:\t\%d (%.2f %%)", $numOfReadsWithAdapter , ($numOfReadsWithAdapter)/$total_num*100);
    }
    printf $fh ("\nDiscarded reads #:\t\%d (%.2f %%)\n", $total_num - $trimmed_num , ($total_num - $trimmed_num)/$total_num*100);
    printf $fh ("Trimmed bases:\t\%d (%.2f %%)\n", $total_raw_seq_len - $total_trimmed_seq_len, ($total_raw_seq_len - $total_trimmed_seq_len)/$total_raw_seq_len*100);
    
    return (0);
}


sub plot_by_R
{
    open (R,">$outDir/tmp$$.R");
    print R <<RSCRIPT; 
def.par <- par(no.readonly = TRUE) # get default parameters

pdf(file = \"$plots_file\",width = 10, height = 8)
#trimmed summary
txt<-c(
"Before Trimming",
paste("Reads #:", prettyNum(\"$total_num\",big.mark=",")),
paste("Total bases:", prettyNum(\"$total_raw_seq_len\",big.mark=",")),
paste("Reads Length:",sprintf("%.2f",$total_raw_seq_len/$total_num)),
" ",
"After Trimming",
paste("Reads #:", prettyNum(\"$trimmed_num\",big.mark=","),sprintf("(%.2f %%)", $trimmed_num/$total_num*100)),
paste("Total bases:",prettyNum(\"$total_trimmed_seq_len\",big.mark=","),sprintf("(%.2f %%)",$total_trimmed_seq_len/$total_raw_seq_len*100))
)



if($trimmed_num>0)
{
txt<-c(txt,
paste("Mean Reads Length:", sprintf("%.2f",$total_trimmed_seq_len/$trimmed_num))
)
}

if($paired_seq_num >0){
     txt <- c(txt, 
 paste("  Paired Reads #: ", prettyNum(\"$paired_seq_num\",big.mark=","), sprintf("(%.2f %%)", $paired_seq_num/$trimmed_num*100)),
 paste("  Paired total bases: ", prettyNum(\"$total_paired_bases\",big.mark=","),sprintf("(%.2f %%)",$total_paired_bases/$total_trimmed_seq_len*100)),
 paste("  Unpaired Reads: ", prettyNum($trimmed_num - $paired_seq_num,big.mark=","),sprintf("(%.2f %%)", ($trimmed_num - $paired_seq_num)/$trimmed_num*100)),
 paste("  Unpaired total bases: ",prettyNum($total_trimmed_seq_len - $total_paired_bases,big.mark=","),sprintf("(%.2f %%)" , ($total_trimmed_seq_len - $total_paired_bases)/$total_trimmed_seq_len*100))
 )
}

txt<-  c( txt, " ")
if ($filter_adapter)
{
   txt<-  c( txt,
   paste("Reads with Adapters/Primers #: ", prettyNum($numOfReadsWithAdapter,big.mark=","), sprintf("(%.2f %%)", ($numOfReadsWithAdapter)/$total_num*100))
   )
}
    txt<-  c( txt,
   paste("Discarded reads #: ", prettyNum($total_num - $trimmed_num,big.mark=","), sprintf("(%.2f %%)", ($total_num - $trimmed_num)/$total_num*100)),
   paste("Discarded bases: ",prettyNum($total_raw_seq_len - $total_trimmed_seq_len,big.mark=","),sprintf("(%.2f %%)",  ($total_raw_seq_len - $total_trimmed_seq_len)/$total_raw_seq_len*100))
  )


plot(1:80,xaxt=\'n\',yaxt=\'n\',type=\'n\',ylab=\'\',xlab=\'\')
#title(paste(\"$prefix\",\"QC report\"),sub = 'DOE Joint Genome Institute/Los Alamos National Laboratory', adj = 0.5, col.sub='darkblue',font.sub=2,cex.sub=0.8)
for (i in seq(along=txt))
{
 text(1,(85-i*5),txt[i],adj=0,font=2)
}
 

#lenght histogram
lengthfile<-read.table(file=\"$length_histogram\")
lengthList<-as.numeric(lengthfile\$V1)
lengthCount<-as.numeric(lengthfile\$V2)
lenAvg<-sum(lengthList * lengthCount)/sum(lengthCount)
lenStd<-sqrt(sum(((lengthList - lenAvg)**2)*lengthCount)/sum(lengthCount))
lenMax<-max(lengthList[lengthCount>0])
lenMin<-min(lengthList[lengthCount>0])
barplot(lengthCount/1000000,names.arg=lengthList,xlab=\"Length\",ylab=\"Count (millions)\",main=\"Reads Length Histogram\",cex.names=0.8)
legend.txt<-c(paste(\"Mean\",sprintf (\"%.2f\",lenAvg),\"±\",sprintf (\"%.2f\",lenStd)),paste(\"Max\",lenMax),paste(\"Min\",lenMin))
legend('topleft',legend.txt,bty='n')

#readGC plot
baseP<-read.table(file=\"$nuc_composition_file\")
Apercent<-baseP\$V2[which(baseP\$V1==\"A\")]
ApercentCount<-baseP\$V3[which(baseP\$V1==\"A\")]
Tpercent<-baseP\$V2[which(baseP\$V1==\"T\")]
TpercentCount<-baseP\$V3[which(baseP\$V1==\"T\")]
Cpercent<-baseP\$V2[which(baseP\$V1==\"C\")]
CpercentCount<-baseP\$V3[which(baseP\$V1==\"C\")]
Gpercent<-baseP\$V2[which(baseP\$V1==\"G\")]
GpercentCount<-baseP\$V3[which(baseP\$V1==\"G\")]
#Npercent<-baseP\$V2[which(baseP\$V1==\"N\")]
#NpercentCount<-baseP\$V3[which(baseP\$V1==\"N\")]
GCpercent<-baseP\$V2[which(baseP\$V1==\"GC\")]
GCpercentCount<-baseP\$V3[which(baseP\$V1==\"GC\")]
aAvg<-sum(Apercent * ApercentCount)/sum(ApercentCount)
aStd<-sqrt(sum(((Apercent - aAvg)**2)*ApercentCount)/sum(ApercentCount))
tAvg<-sum(Tpercent * TpercentCount)/sum(TpercentCount)
tStd<-sqrt(sum(((Tpercent - tAvg)**2)*TpercentCount)/sum(TpercentCount))
cAvg<-sum(Cpercent * CpercentCount)/sum(CpercentCount)
cStd<-sqrt(sum(((Cpercent - cAvg)**2)*CpercentCount)/sum(CpercentCount))
gAvg<-sum(Gpercent * GpercentCount)/sum(GpercentCount)
gStd<-sqrt(sum(((Gpercent - gAvg)**2)*GpercentCount)/sum(GpercentCount))
#nAvg<-sum(Npercent * NpercentCount)/sum(NpercentCount)
#nStd<-sqrt(sum(((Npercent - nAvg)**2)*NpercentCount)/sum(NpercentCount))
gcAvg<-sum(GCpercent * GCpercentCount)/sum(GCpercentCount)
gcStd<-sqrt(sum(((GCpercent - gcAvg)**2)*GCpercentCount)/sum(GCpercentCount))
GCaggregate<-tapply(GCpercentCount,list(cut(GCpercent,breaks=c(seq(0,100,1)))),FUN=sum)
Aaggregate<-tapply(ApercentCount,list(cut(Apercent,breaks=c(seq(0,100,1)))),FUN=sum)
Taggregate<-tapply(TpercentCount,list(cut(Tpercent,breaks=c(seq(0,100,1)))),FUN=sum)
Caggregate<-tapply(CpercentCount,list(cut(Cpercent,breaks=c(seq(0,100,1)))),FUN=sum)
Gaggregate<-tapply(GpercentCount,list(cut(Gpercent,breaks=c(seq(0,100,1)))),FUN=sum)


par(fig=c(0,0.75,0,1),mar=c(5,4,4,2),xpd=FALSE,cex.main=1.2)
plot(GCaggregate/1000000,xlim=c(0,100),type=\"h\",lwd=4, main=\"Reads GC content\",xlab=\"GC (%)\",ylab=\"Number of reads (millions)\",lend=2)
legend.txt<-c(paste(\"GC\",sprintf (\"%.2f%%\",gcAvg),\"±\",sprintf (\"%.2f\",gcStd)))
legend('topright',legend.txt,bty='n')

par(fig=c(0.75,1,0.75,1), mar=c(3, 2, 2, 2),new=TRUE,cex.main=1)
legend.txt<-c(paste(\"A\",sprintf (\"%.2f%%\",aAvg),\"±\",sprintf (\"%.2f\",aStd)))
plot(Aaggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab="",ylab="",)

par(fig=c(0.75,1,0.5,0.75),mar=c(3, 2, 2, 2),new=TRUE)
legend.txt<-c(paste(\"T\",sprintf (\"%.2f%%\",tAvg),\"±\",sprintf (\"%.2f\",tStd)))
plot(Taggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab="",ylab="")

par(fig=c(0.75,1,0.25,0.5),mar=c(3, 2, 2, 2),new=TRUE)
legend.txt<-c(paste(\"C\",sprintf (\"%.2f%%\",cAvg),\"±\",sprintf (\"%.2f\",cStd)))
plot(Caggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab="",ylab="")

par(fig=c(0.75,1,0,0.25),mar=c(3, 2, 2, 2),new=TRUE)
legend.txt<-c(paste(\"G\",sprintf (\"%.2f%%\",gAvg),\"±\",sprintf (\"%.2f\",gStd)))
plot(Gaggregate/1000000,xlim=c(0,50),type=\"h\",lwd=2,main=legend.txt,xlab="",ylab="")
par(def.par)#- reset to default


#ATCG composition per base ATCG plot
baseM<-read.table(file=\"$base_matrix\")
aBase<-baseM\$V1
tBase<-baseM\$V2
cBase<-baseM\$V3
gBase<-baseM\$V4
nBase<-baseM\$V5

aPer<-(aBase/rowSums(baseM))*100
tPer<-(tBase/rowSums(baseM))*100
cPer<-(cBase/rowSums(baseM))*100
gPer<-(gBase/rowSums(baseM))*100

xpos<-seq(1,length(aBase),1)
plot(xpos,aPer,col=\'green3\',type=\'l\',xaxt=\'n\',xlab=\'Base\',ylab=\'Base content (%)\',ylim=c(1,100))
lines(xpos,tPer,col=\'red\')
lines(xpos,cPer,col=\'blue\')
lines(xpos,gPer,col=\'black\')
axis(1,at=xpos,labels=xpos)
legend(\'topright\',c(\'A\',\'T\',\'C\',\'G\'),col=c(\'green3\',\'red\',\'blue\',\'black\'),box.col=0,lwd=1)
title(\"Nucleotide Content Per Cycle\")

if ($N_num_cutoff >0){
#N composition per Base plot
plot(xpos,nBase/$trimmed_num*1000000,col=\'red\',type=\'l\',xaxt=\'n\',xlab=\'Position\',ylab=\'N Base count per million reads\',ylim=c(0,max(nBase/$trimmed_num*1000000)))
axis(1,at=xpos,labels=xpos)
title(\"N Nucleotide Content Per Cycle\")
}


if(file.exists(\"$kmer_rarefaction_file\")){

kmerfile<-read.table(file=\"$kmer_rarefaction_file\")
sampling_size<-sum(kmerfile\$V1)
sampling<-""
if(sampling_size<$trimmed_num){
    sampling<-paste(\"(Sampling\",format(sampling_size/1000000,digit=3),\"M Reads)\")
}
cumSeqNum<-cumsum(kmerfile\$V1);
plot(cumSeqNum/1000000,kmerfile\$V3/1000000,xlab=paste(\"Number of Sequence (million)\",sampling), ylab=\"Number of Distinct Kmer (million,k=$kmer)\",type=\'l\',lty=2)
lines(cumSeqNum/1000000,kmerfile\$V2/1000000,col='blue',lwd=2)
title(\"Kmer Rarefaction Curve\")
y<-kmerfile\$V2/1000000
x<-cumSeqNum/1000000
lres<-lm(y~x)
# y=ax+b
a<-format(coef(lres)[[2]], digits = 2)
b<-format(coef(lres)[[1]], digits = 2)

}

if(file.exists(\"$kmer_histogram_file\")){

sampling<-""
if(sampling_size<$trimmed_num){
    sampling<-paste(\"(Sampling\",format(sampling_size/1000000,digit=3),\"M Reads)\")
}
kmerHfile<-read.table(file=\"$kmer_histogram_file\")
barplot(kmerHfile\$V2/1000000,names.arg=kmerHfile\$V1,xlab=\"Kmer Frequency\",log=\'y\',ylab=\"Number of Kmers (millions,k=$kmer)\",main=paste(\"Kmer Frequency Histogram\",sampling),col=\'darkcyan\',border=\'darkcyan\')
total_kmer<-sum(kmerHfile\$V2)
legend('topleft',paste(\"Total: \",total_kmer),bty='n')
par(fig=c(0.5,0.9,0.5,1), new=TRUE)
barplot(kmerHfile\$V2/1000000,xlim=c(0,100),names.arg=kmerHfile\$V1,log=\'y\',col=\'darkcyan\',border=\'darkcyan\')
par(fig=c(0,1,0,1))
par(def.par)#- reset to default

}


# read avg quality count barplot 
Qhist_file<-read.table(file=\"$avg_quality_histogram\",header=TRUE)
par(mar=c(5,4,4,4))
cumulate<-cumsum(Qhist_file\$readsNum)
plot(Qhist_file\$Score,Qhist_file\$readsNum/1000000,type=\'h\',xlim=c(max(Qhist_file\$Score),min(Qhist_file\$Score)),xlab=\"Avg Score\", ylab=\"Reads Number (millions)\",lwd=12,lend=2)
title('Reads Average Quality Histogram')
par(new=TRUE)
plot(Qhist_file\$Score,cumulate/sum(Qhist_file\$readsNum)*100,type='l',xlim=c(max(Qhist_file\$Score),min(Qhist_file\$Score)),yaxt='n',xaxt='n',ylab="",xlab="",col='blue',lwd=3)
axis(4,col='blue',col.ticks='blue',col.axis='blue')
mtext(side=4,'Cumulative Percentage',line=2,col='blue')
Qover20Reads<-sum(as.numeric(Qhist_file\$readsNum[Qhist_file\$Score>=20]))
Qover20ReadsPer<-sprintf(\"%.2f%%\",Qover20Reads/sum(Qhist_file\$readsNum)*100)
Qover20Bases<-sum(as.numeric(Qhist_file\$readsBases[Qhist_file\$Score>=20]))
Qover20AvgLen<-sprintf(\"%.2f\",Qover20Bases/Qover20Reads)
mtext(side=3,paste(\"Number of Q>=20 reads:\",prettyNum(Qover20Reads,big.mark=\",\"),\"(\",Qover20ReadsPer,\")\",\", mean Length:\",Qover20AvgLen),adj=0,cex=0.8,line=0.5)
par(def.par)#- reset to default


# read in matrix file for the following three plots
z<-as.matrix(read.table(file=\"$quality_matrix\"));
x<-1:nrow(z)
y<-1:ncol(z)
y<-y-1

#quality boxplot per base
is.wholenumber <- function(x, tol = .Machine\$double.eps^0.5)  abs(x - round(x)) < tol
plot(1:length(x),x,type=\'n\',xlab=\"Position\",ylab=\"Quality score\", ylim=c(0,max(y)+1),xaxt=\'n\')
axis(1,at=x,labels=TRUE)
title(\"Quality Boxplot Per Cycle\")

for (i in 1:length(x)) {
  total<-sum(z[i,])
  qAvg<-sum(y*z[i,])/total
  if (is.wholenumber(total/2))
  {
     med<-( min(y[cumsum((z[i,]))>=total/2]) + min(y[cumsum((z[i,]))>=total/2+1]) )/2
  }
  else
  {
     med<-min(y[cumsum((z[i,]))>=ceiling(total/2)])
  }

  if (is.wholenumber(total/4))
  {
     Q1<-( min(y[cumsum((z[i,]))>=total/4]) + min(y[cumsum((z[i,]))>=total/4+1]) )/2
  }
  else
  {
     Q1<-min(y[cumsum((z[i,]))>=round(total/4)])
  }

  if (is.wholenumber(total/4*3))
  {
     Q3<-( min(y[cumsum((z[i,]))>=total/4*3]) + min(y[cumsum((z[i,]))>=total/4*3+1]) )/2
  }
  else
  {
     Q3<-min(y[cumsum((z[i,]))>=round(total/4*3)])
  }
  maxi<-max(y[z[i,]>0])
  mini<-min(y[z[i,]>0])
  #if (Q1 == 'Inf') {Q1 = maxi}
  if (Q3 == \'Inf\') {Q3 = maxi}
  IntQ<-Q3-Q1
  mini<-max(mini,Q1-1.5*IntQ)
  maxi<-min(maxi,Q3+1.5*IntQ)
  rect(i-0.4,Q1,i+0.4,Q3,col=\'bisque\')
  lines(c(i,i),c(Q3,maxi),lty=2)
  lines(c(i,i),c(mini,Q1),lty=2)
  lines(c(i-0.4,i+0.4),c(mini,mini))
  lines(c(i-0.4,i+0.4),c(maxi,maxi))
  lines(c(i-0.4,i+0.4),c(med,med))
  #points(i,qAvg,col=\'red\')
  reads_num<-prettyNum($trimmed_num,big.mark=",")
  reads_base<-prettyNum($total_trimmed_seq_len,big.mark=",")
  abline(h=20, col = \"gray60\")
  legend(\"bottomleft\",c(paste(\"# Reads: \",reads_num),paste(\"# Bases:\",reads_base)))
## for outliers
#points()
}

#quality 3D plot
persp(x,y,z/1000000,theta = 50, phi = 30, expand = 0.7, col = \"#0000ff22\",ntick=10,ticktype=\"detailed\",xlab=\'Position\',ylab=\'Score\',zlab=\"\",r=6,shade=0.75)
mtext(side=2, \"Frequency (millions)\",line=2)
title(\"Quality 3D plot. (Position vs. Score vs. Frequency)\")

#Quality count bar plot
col<-colSums(z)
less30columnNum<-length(col)-$highest_illumina_score+30-1
atleast30columnNum<-$highest_illumina_score-30+1
color<-c(rep(\'blue\',less30columnNum),rep(\'darkgreen\',atleast30columnNum))
over30per<-sprintf(\"%.2f%%\",sum(col[(less30columnNum+1):length(col)])/sum(col)*100)
countInM<-col/1000000
plot(seq($lowest_illumina_score,$highest_illumina_score,1),countInM,col=color,type='h',ylab=\"Total (million)\",xlab=\"Q score\",lwd=12,lend=2,bty='n')
abline(v=29.5,col='darkgreen')
text(30,(max(countInM)-min(countInM))*0.9,labels=\">=Q30\",cex=0.8,adj=0,col=\'darkgreen\')
text(30,(max(countInM)-min(countInM))*0.85,labels=over30per,cex=0.8,adj=0,col=\'darkgreen\')
title(\"Quality report\")

tmp<-dev.off()

quit()

RSCRIPT

    close R;
    system ("R --vanilla --silent --quiet < $outDir/tmp$$.R 1>/dev/null");
    unless ($debug){system ("rm $outDir/tmp$$.*");}
    system ("rm $outDir/Rplots.pdf") if (-e "$outDir/Rplots.pdf");
}

sub print_quality_report_files{
    open (OUT, ">$quality_matrix");
    open (BASE, ">$base_matrix");
    foreach my $pos2(sort {$a<=>$b} keys %position){
        my $q_string;
        my $n_string;
        for ($lowest_illumina_score..$highest_illumina_score)
        {
            if ($position{$pos2}->{$_}){
              $q_string .=  $position{$pos2}->{$_}."\t";
            }
            else
            {
              $q_string .= "0\t";
            }
        }
        for my $base ("A","T","C","G","N")
        { 
            if ($base_position{$pos2}->{$base}){
              $n_string .= $base_position{$pos2}->{$base}."\t";
            }
            else
            {
              $n_string .= "0\t";
            }
        }
        $q_string =~ s/\t$/\n/;
        $n_string =~ s/\t$/\n/;
        print OUT $q_string;
        print BASE $n_string;
    }
    close OUT;
    close BASE;

    open (OUT2, ">$avg_quality_histogram");
    my $Key=$highest_illumina_score;
    my $h_print_string;
    while ($Key >= 0)
    {
        if(defined $AverageQ{$Key}->{reads}){
           $h_print_string .= "$Key\t".
           $AverageQ{$Key}->{reads}. "\t".$AverageQ{$Key}->{bases}."\n";
        }else{
           $h_print_string .= "$Key\t".
              "0\t".
              "0\n";
        }
        --$Key;
    }
    print OUT2 "Score\treadsNum\treadsBases\n";
    print OUT2 $h_print_string; 
    close OUT2;

    open (OUT3, ">$nuc_composition_file");
    for my $nuc ("A","T","C","G","N","GC")
    {
      foreach my $key ( sort {$a<=>$b} keys %{$base_content{$nuc}})
      {
             print OUT3 "$nuc\t$key\t${$base_content{$nuc}}{$key}\n";
      }
    } 
    close OUT3;
 
    open (OUT4,">$length_histogram");
	my @len_list= sort {$a<=>$b} keys %len_hash;
    for my $key (1..$len_list[-1])
    {
	    if ($len_hash{$key}) 
		{
            print OUT4 "$key\t$len_hash{$key}\n";
		}
		else
		{
			print OUT4 "$key\t0\n";
		}
    }
    close OUT4;
}

sub file_check 
{
    #check file exist and non zero size
    my $file=shift;
    my $exist=-1;
    if (-e $file) {$exist=1};
    if (-z $file) {$exist=-1};
    return $exist;
}

sub qc_process {
  my ($input, $input2) = @_;
  my ($h1,$s,$s_trimmed,$h2,$q, $q_trimmed); my $len=0; my $trim_len=0;
  my ($r2_h1,$r2_s,$r2_s_trimmed,$r2_h2,$r2_q,$r2_q_trimmed); my $r2_len=0; my $r2_trim_len=0;
  my %stats;
  my $seq_r;
  my $avg_q;
  my ($pos5,$pos3);
  my $numOfReadsWithAdapter;
  my ($raw_seq_num,$total_raw_seq_len);
  my ($trim_seq_num_1,$trim_seq_num_2,$trim_seq_len, $total_trim_seq_len,@trim_seq_len);
  my ($paired_seq_num,$total_paired_bases); 
  my (%tmp1,%tmp2);
  my ($drop_1,$drop_2)=(1,1);
  my ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$input", qr/\.[^.]*/ );
  my $trim_output_1="$outDir/${i_file_name}_trim.fastq";
  my $trim_output_2;
  my $trim_output_unpaired;
  my $trim_output_discard="$outDir/${i_file_name}_trim_discard.fastq";
  open(IN,"$input") or die "$input\t$!";
  if (! $qc_only)
  {
     open(OUT,"> $trim_output_1");
     open(DISCARD,">$trim_output_discard") if ($output_discard);
  }
  if ($input2) # paired mate
  {
     open(IN2,"$input2") or die "$input2\t$!";
     my ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$input2", qr/\.[^.]*/ );
     $trim_output_2="$outDir/${i_file_name}_trim.fastq";
     $trim_output_unpaired="$outDir/${i_file_name}_trim_unpaired.fastq";
     if (! $qc_only)
     {
         open(OUT2,"> $trim_output_2");
         open(UNPAIR,"> $trim_output_unpaired");
     }
  }

  while ($h1 = <IN>) {  # read first header
        $drop_1=0;
        $raw_seq_num++;
        $s = <IN>;  # read sequence
        chomp $s;
        $h2 = <IN>;  # read second header
        $q = <IN>;  # read quality scores
        chomp $q;
        $len = length($q);
        $total_raw_seq_len += $len;
        if ($filter_adapter)
        {
            $drop_1=&filter_adapter($s,$filterAdapterMismatchRate);
            $numOfReadsWithAdapter++ if ($drop_1);
            $drop_1 =0 if ($qc_only);
        }
        $drop_1=1 if ($len < $opt_min_L);
        if ($drop_1==0){
           
           ($s_trimmed,$q_trimmed,$pos5,$pos3)= &bwa_trim ($len,$s,$q);
           $trim_len=length($q_trimmed); 
           $drop_1=1 if ($trim_len < $opt_min_L);
        }
        if ($drop_1==0)
        {
            ($seq_r,$drop_1)=&get_base_and_quality_info($s_trimmed,$q_trimmed,$trim_len,$pos5,$pos3,\%stats);
            %stats=%{$seq_r};
        }
        if ($drop_1==0){  # pass all filters...
            $q_trimmed=&quality_encoding_coversion($q_trimmed,$ascii,$out_offset) if ($ascii != $out_offset);
            $trim_seq_num_1++;
            push @trim_seq_len, $trim_len;
            $total_trim_seq_len += $trim_len;
        }
        
        if ($input2) {
           $drop_2=0;
           $r2_h1 = <IN2>;
           $raw_seq_num++;
           $r2_s = <IN2>;  # mate read sequence
           chomp $r2_s;
           $r2_h2 = <IN2>;  # mate read second header
           $r2_q = <IN2>;  # mate read quality scores
           chomp $r2_q;
           $r2_len = length($r2_q);
           $total_raw_seq_len += $r2_len;
           if ($filter_adapter)
           {
               $drop_2=&filter_adapter($r2_s,$filterAdapterMismatchRate);
               $numOfReadsWithAdapter++ if ($drop_2);
               $drop_2 =0 if ($qc_only);
           }
           $drop_2=1 if ($r2_len < $opt_min_L);
           if ($drop_2==0){
              ($r2_s_trimmed, $r2_q_trimmed, $pos5, $pos3)=&bwa_trim ($r2_len,$r2_s,$r2_q);
               $r2_trim_len=length($r2_q_trimmed);
               $drop_2=1 if ($r2_trim_len < $opt_min_L);
           }
           if ($drop_2==0)
           {
               ($seq_r,$drop_2)=&get_base_and_quality_info($r2_s_trimmed,$r2_q_trimmed,$r2_trim_len,$pos5,$pos3,\%stats);
               %stats=%{$seq_r};
           }
          
           if ($drop_2==0){ # pass all filters
               $r2_q_trimmed=&quality_encoding_coversion($r2_q_trimmed,$ascii,$out_offset) if ($ascii != $out_offset);
               $trim_seq_num_2++;
               push @trim_seq_len, $r2_trim_len;
               $total_trim_seq_len += $r2_trim_len;
           }
        }
        if ($drop_1==0 and $drop_2==0)
        {
            $paired_seq_num = $paired_seq_num + 2;
            $total_paired_bases += $trim_len + $r2_trim_len;  
        }
 
        # output trimmed files
        if (! $qc_only)
        {
            if ($drop_1==0 and $drop_2==0)
            {
	        print OUT $h1.$s_trimmed."\n".$h2.$q_trimmed."\n";
                print OUT2 $r2_h1.$r2_s_trimmed."\n".$r2_h2.$r2_q_trimmed."\n";
            }
            elsif ($drop_1==1 and $drop_2==0)
            {
                print DISCARD $h1.$s."\n".$h2.$q."\n" if ($output_discard);
                print UNPAIR $r2_h1.$r2_s_trimmed."\n".$r2_h2.$r2_q_trimmed."\n";
            }
            elsif ($drop_1==0 and $drop_2==1)
            {
                if ($input2)
                {
                   print UNPAIR $h1.$s_trimmed."\n".$h2.$q_trimmed."\n";
                   print DISCARD $r2_h1.$r2_s."\n".$r2_h2.$r2_q."\n" if ($output_discard);
                }
                else
                {
                   print OUT $h1.$s_trimmed."\n".$h2.$q_trimmed."\n";
                }
            }
            elsif ($output_discard)
            {   
                print DISCARD $h1.$s."\n".$h2.$q."\n";
                if ($input2)
                {
                   print DISCARD $r2_h1.$r2_s."\n".$r2_h2.$r2_q."\n";
                }
            }
        }
  } # end while

  my ($trim_seq_len_std,$trim_seq_len_avg,$max,$min,$median)=&standard_deviation(@trim_seq_len);
  
  close(IN);
  close(OUT);
  close IN2 if ($input2);
  close OUT2 if ($input2);
  close UNPAIR if ($input2);
  close DISCARD;
  $stats{raw_seq_num}=$raw_seq_num;
  $stats{trim_seq_num}=$trim_seq_num_1+$trim_seq_num_2;
  $stats{trim_seq_num_1}=$trim_seq_num_1;
  $stats{trim_seq_num_2}=$trim_seq_num_2;
  $stats{total_raw_seq_len}=$total_raw_seq_len;
  $stats{total_trim_seq_len}=$total_trim_seq_len;
  $stats{trim_seq_len_std}=$trim_seq_len_std;
  $stats{trim_seq_len_avg}=$trim_seq_len_avg;
  $stats{max}=$max;
  $stats{min}=$min;
  $stats{median}=$median;
  $stats{paired_seq_num}=$paired_seq_num;
  $stats{total_paired_bases}=$total_paired_bases;
  $stats{trim_file_name_1}=$trim_output_1;
  $stats{trim_file_name_2}=$trim_output_2;
  $stats{numOfReadsWithAdapter}=$numOfReadsWithAdapter;
  return \%stats;
}

sub filters  # not used
{
    my $s = shift;
    my $q = shift;
    my $drop = 0;
    # "N" base filter
    my $N_num = $s=~ tr/Nn/Nn/;
    if ( $N_num > $N_num_cutoff and $drop==0) {$drop=1;}
    #low complexity filter
    if (&low_complexity_filter($s,$low_complexity_cutoff_ratio) and $drop==0) {$drop=1;}
    
    return $drop;
} 
     
sub bwa_trim2
{
        # bwa trim implementation from both 5' and 3' end
        my ($len,$s,$q) = @_;
        # trim 3' end
        my $pos_3 = $len;
        my $final_pos_3 =$pos_3+1;  
        my $area=0;  
        my $maxArea=0;
        
	while ($pos_3>0 and $area>=0) {
		$area += $opt_q - (ord(substr($q,$pos_3-1,1))-$ascii);
		if ($area > $maxArea) {
			$maxArea = $area;
			$final_pos_3 = $pos_3;
		}
		$pos_3--;
	}
        # trim 5' end
        my $pos_5=1; 
        my $final_pos_5=0;
        $maxArea=0;
        $area=0;
	while ($pos_5<$final_pos_3 and $area>=0) {
		$area += $opt_q - (ord(substr($q,$pos_5-1,1))-$ascii);
		if ($area > $maxArea) {
			$maxArea = $area;
			$final_pos_5 = $pos_5;
		}
		$pos_5++;
	}
        
        $s = substr($s,$final_pos_5,$final_pos_3-$final_pos_5-1);
        $q = substr($q,$final_pos_5,$final_pos_3-$final_pos_5-1);
        return ($s,$q,$final_pos_5,$final_pos_3);
}

sub bwa_trim
{
    # bwa trim implementation from both 5' and 3' end
    # at least scan 5 bases from both end and 2 more bases after the negative area
    my ($len,$s,$q) = @_;
    my $at_least_scan=5;
    my $num_after_neg=2;
      
    # trim 3' end
    my $pos_3 = $len;
    my $final_pos_3 =$pos_3+1;  
    my $area=0;  
    my $maxArea=0;
    while ($at_least_scan) 
    {
        $at_least_scan--;    
        if ($pos_3>$num_after_neg and $area>=0) {$at_least_scan=$num_after_neg;}
	    $area += $opt_q - (ord(substr($q,$pos_3-1,1))-$ascii);
	    if ($area > $maxArea) {
		    $maxArea = $area;
		    $final_pos_3 = $pos_3;
	    }
	    $pos_3--;
    }
    # trim 5' end
    my $pos_5=1; 
    my $final_pos_5=0;
    $maxArea=0;
    $area=0;
    $at_least_scan=5;
    while ($at_least_scan) 
    {
        $at_least_scan--; 
        if ($pos_5<($final_pos_3-$num_after_neg) and $area>=0) {$at_least_scan=$num_after_neg;}
	    $area += $opt_q - (ord(substr($q,$pos_5-1,1))-$ascii);
	    if ($area > $maxArea) {
		    $maxArea = $area;
		    $final_pos_5 = $pos_5;
	    }
	    $pos_5++;
    }
    
    if ($final_pos_3<=$final_pos_5)
    {
        # since at least scan 5 bases, we need to avoid negative length on substr
        # return empty and it will be filtered by read length;
        $s = substr($s,0,0);
        $q = substr($q,0,0);
    }
    else
    {
        $s = substr($s,$final_pos_5,$final_pos_3-$final_pos_5-1);
        $q = substr($q,$final_pos_5,$final_pos_3-$final_pos_5-1);
    }
    return ($s,$q,$final_pos_5,$final_pos_3);
}


sub get_base_and_quality_info
{
     my $s=shift;
     my $q=shift;
     my $len=shift;
     my $start_pos=shift;
     my $end_pos=shift;
     my $stats=shift;
     my $total_q=0;
     my $avg_q=0;
     my $drop=0;
     my ($a_Base,$t_Base,$c_Base,$g_Base,$n_Base)=(0,0,0,0,0);
     my ($a_Base_percent,$t_Base_percent,$c_Base_percent,$g_Base_percent,$n_Base_percent);
     my $base_percent_string;
     my ($previous_base,$dinucleotide,%di_nucleotide);
     my %seq=%{$stats};
     $start_pos++;
     $end_pos--;
     for my $pos($start_pos..$end_pos)
     {
         my $q_digit=ord(substr($q,$pos-$start_pos,1))-$ascii;
         $drop=1 if ($q_digit < $stringent_cutoff);
         $seq{qual}->{$pos}->{$q_digit}++;
         $total_q += $q_digit; 
         my $base=uc(substr($s,$pos-$start_pos,1));
         $seq{base}->{$pos}->{$base}++;
         $a_Base++ if ($base =~ /A/);  
         $t_Base++ if ($base =~ /T/);  
         $c_Base++ if ($base =~ /C/);  
         $g_Base++ if ($base =~ /G/);  
         $n_Base++ if ($base =~ /N/);  
         if ($previous_base and ($previous_base ne $base))
         {
            $dinucleotide = $previous_base.$base;
            $di_nucleotide{$dinucleotide}++;
         }
         $previous_base=$base;
     }
     $avg_q = int ($total_q/$len);
     if (!$qc_only)
     {
         # apply filters: N num, average quality, low complexity
         $drop=1 if ($n_Base>$N_num_cutoff);
         $drop=1 if ($avg_q < $opt_avg_cutoff);
         $drop=1 if (($a_Base/$len) >= $low_complexity_cutoff_ratio);
         $drop=1 if (($t_Base/$len) >= $low_complexity_cutoff_ratio);
         $drop=1 if (($c_Base/$len) >= $low_complexity_cutoff_ratio);
         $drop=1 if (($g_Base/$len) >= $low_complexity_cutoff_ratio);
         if ($drop == 0)
         {
           foreach my $di_nuc (keys %di_nucleotide)
           {
             if (($di_nucleotide{$di_nuc}*2/$len) >= $low_complexity_cutoff_ratio)
             {
              $drop=1;
              delete $di_nucleotide{$di_nuc};
              last; #end for loop
             }
             delete $di_nucleotide{$di_nuc};
           }
         }
     }
     if ($drop == 1)  #substract the position matrix by one because the dropped read
     {
         for my $pos($start_pos..$end_pos)
         {
            my $q_digit=ord(substr($q,$pos-$start_pos,1))-$ascii;
            $seq{qual}->{$pos}->{$q_digit}--;
            my $base=uc(substr($s,$pos-$start_pos,1));
            $seq{base}->{$pos}->{$base}--;
         }
         return (\%seq,$drop);
     }
     else
     {
         $a_Base_percent = sprintf("%.2f",$a_Base/$len*100);
         $t_Base_percent = sprintf("%.2f",$t_Base/$len*100);
         $c_Base_percent = sprintf("%.2f",$c_Base/$len*100);
         $g_Base_percent = sprintf("%.2f",$g_Base/$len*100);
         $n_Base_percent = sprintf("%.2f",$n_Base/$len*100);
         $seq{ReadAvgQ}->{$avg_q}->{readsNum}++;
         $seq{ReadAvgQ}->{$avg_q}->{basesNum}+=$len;
         $seq{ReadLen}->{$len}++;
         $seq{Base_content}->{A}->{$a_Base_percent}++;
         $seq{Base_content}->{T}->{$t_Base_percent}++;
         $seq{Base_content}->{C}->{$c_Base_percent}++;
         $seq{Base_content}->{G}->{$g_Base_percent}++;
         $seq{Base_content}->{N}->{$n_Base_percent}++;
         $seq{Base_content}->{GC}->{$c_Base_percent+$g_Base_percent}++;
 
         return (\%seq,$drop);
     }
} 

sub standard_deviation {
  my(@numbers) = @_;
  #Prevent division by 0 error in case you get junk data
  return undef unless(scalar(@numbers));
  @numbers= sort {$b<=>$a} @numbers;
  my $max=$numbers[0];
  my $min=$numbers[-1];
  my $median;
  if (scalar(@numbers) % 2) {
      $median = $numbers[int(scalar(@numbers)/2)];
  } else {
      $median = ($numbers[scalar(@numbers)/2] + $numbers[scalar(@numbers)/2 - 1]) / 2;
  }

  # Step 1, find the mean of the numbers
  my $total1 = 0;
  foreach my $num (@numbers) {
  $total1 += $num;
  }
  my $mean1 = $total1 / (scalar @numbers);

  # Step 2, find the mean of the squares of the differences
  # between each number and the mean
  my $total2 = 0;
  foreach my $num (@numbers) {
  $total2 += ($mean1-$num)**2;
  }
  my $mean2 = $total2 / (scalar @numbers);

  # Step 3, standard deviation is the square root of the
  # above mean
  my $std_dev = sqrt($mean2);
  return ($std_dev,$mean1,$max,$min,$median);
}


sub low_complexity_filter {  # not used 
	my ($seq,$cut_off_ratio) = @_;
	my $seq_len = length ($seq);
    my @low_complex_array=("A","T","C","G","AT","AC","AG","TA","TC","TG","CA","CT","CG","GA","GT","GC");
    my $filter=0;
    for my $item (@low_complex_array)
    {
      my $item_len = length ($item);
      my $num_low_complex = $seq =~ s/$item/$item/g;
      if (($num_low_complex*$item_len/$seq_len) >= $cut_off_ratio)
      {
        $filter=1;
        last; #end for loop
      }
    }
    return ($filter);
}

sub random_subsample
{
  my $total_num=shift;
  my $output_num=shift;
  my %get;
  if ($output_num <= 0)
  {
     $get{1}=0;
  }
  else
  { 
    while (1)
    {
      $get{int(rand($total_num))}=1;
      my $num_of_random=scalar (keys %get);
      last if ($num_of_random == $output_num or $num_of_random == $total_num);
    }
  }
  return (\%get);
}
sub quality_encoding_coversion 
{  
    # given quality acsii string, input offset and output offset
    my $q_string=shift;
    my $input_offset=shift;
    my $out_offset=shift;
    $q_string=~ s/(\S)/chr(ord($1)-$input_offset+$out_offset)/eg;
    return($q_string);
}


sub checkQualityFormat {
    # $offset_value=&checkQualityFormat($fastq_file)
    # $offset_value = -1 means not proper fastq format.
    my $fastq_file=shift;
    # open the files
    if ($fastq_file =~ /gz$/)
    {
        open FQ, "zcat $fastq_file | " or die $!;
    }
    else
    {
        open FQ, "<", $fastq_file or die $!;
    }

    # initiate
    my @line;
    my $l;
    my $number;
    my $offset;
    # go thorugh the file
    my $first_line=<FQ>;
    if ($first_line !~ /^@/) {$offset=-1; return $offset;}
    OUTER:while(<FQ>){
      
      # if it is the line before the quality line
      if($_ =~ /^\+/){
    
       $l = <FQ>; # get the quality line
       @line = split(//,$l); # divide in chars
       for(my $i = 0; $i <= $#line; $i++){ # for each char
        $number = ord($line[$i]); # get the number represented by the ascii char
      
        # check if it is sanger or illumina/solexa, based on the ASCII image at http://en.wikipedia.org/wiki/FASTQ_format#Encoding
        if($number > 74){ # if solexa/illumina
          $offset=64;
          #die "This file is solexa/illumina format\n"; # print result to terminal and die
          last OUTER; 
        }elsif($number < 59){ # if sanger
          $offset=33;
          #die "This file is sanger format\n"; # print result to terminal and die
          last OUTER;
        }
       }
      }
    }
   return $offset;
}

sub split_fastq {
   my ($input,$outDir,$subfile_size)=@_;

   my $seqThisFile = 0;
   my $fileCount = 0;
   my $name="";
   my $seq="";
   my $q_name="";
   my $qual_seq="";
   my @subfiles;
   my $seq_num=0;
   my $total_seq_length=0;
   my ($file_name, $i_path, $i_suffix_1)=fileparse("$input", qr/\.[^.]*/);
   ($file_name, $i_path, $i_suffix)=fileparse("$file_name", qr/\.[^.]*/) if ($i_suffix_1 =~ /gz$/);
   my $file_ext = sprintf("%05d", $fileCount);
   my $out_file_name= $outDir . "/". $file_name . "_" . $file_ext . ".fastq";
   push @subfiles, $out_file_name;
   open (OUTFILE, ">" . $out_file_name ) or die ("Cannot open file for output: $!");
   if ($i_suffix_1 =~ /gz/)
   {
      open (IN, "zcat $input | ") || die "cannot open file $input:$!";
   }
   else
   {
      open (IN, $input) || die "cannot open file $input:$!";
   }
   while (<IN>)
   {
          $name = $_;
          $seq=<IN>;
          $seq =~ s/\n//g;
          while ($seq !~ /\+/)
          {
             $seq .= <IN>;
             $seq =~ s/\n//g;
          }
          my $q_name_pos=index($seq,"+");
          $q_name = substr($seq,$q_name_pos);
          $seq = substr($seq, 0, $q_name_pos);
          my $seq_len = length $seq;
          $qual_seq=<IN>;
          $qual_seq =~ s/\n//g;
          my $qual_seq_len = length $qual_seq;
          while ( $qual_seq_len < $seq_len )
          {
              last if ( $qual_seq_len == $seq_len);
              $qual_seq .= <IN>;
              $qual_seq =~ s/\n//g;
              $qual_seq_len = length $qual_seq;
          }
          $total_seq_length += $seq_len;
          print (OUTFILE "$name$seq\n$q_name\n$qual_seq\n");
          $seq_num++;
          $seqThisFile++;

          if ($seqThisFile == $subfile_size)
          {
              $fileCount++;
	      $seqThisFile = 0;
	      close (OUTFILE) or die( "Cannot close file : $!");
	      $file_ext = sprintf("%05d", $fileCount);
	      $out_file_name= $outDir . "/". $file_name . "_" . $file_ext . ".fastq";
              open (OUTFILE, ">" .  $out_file_name ) or die ("Cannot open file for output: $!") if (! eof);
              push @subfiles, $out_file_name if (! eof);
          }

   }
   my $average_len = $total_seq_length/$seq_num;
   if ( $average_len < $opt_min_L) { print "The input average length $average_len < minimum cutoff length(opt_min_L) $opt_min_L\n."; exit;}
   close (IN)  or die( "Cannot close file : $!");
   close (OUTFILE) or die( "Cannot close file : $!") if (! eof OUTFILE);
   return ($seq_num,@subfiles);
}
   
sub filter_adapter
{
    my ($s,$mismatchRate)=@_;
    my $match_flag=0;
    $mismatchRate = $mismatchRate*100;
    my @adapterSeqs=(
		"GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
		"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
		"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
		"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
		"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
		"CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG",
		"CTGTCTCTTATACACATCT",
		"AGATGTGTATAAGAGACAG",
		"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
		"GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
	);
    
    @adapterSeqs=(@adapterSeqs,&read_artifactFile) if ($artifactFile);

    foreach my $adapter (@adapterSeqs)
    {
        $match_flag = String::Approx::amatch($adapter, ["i", "S ${mismatchRate}%"], $s);
        if ($match_flag)
        {
           return 1;
        }
    }
    return 0;
}

sub read_artifactFile 
{
     open (ARTIFACT,$artifactFile);
     my $seq;
     my @seqs;
     while(<ARTIFACT>)
     {
        chomp;
        if (/^>/)
        {
           if ($seq)
           { 
               push @seqs,$seq;
           }
           $seq="";
        }
        else
        {
           $seq .= $_;
        }
     }
     if ($seq)
     { 
        push @seqs,$seq;
     }
     close ARTIFACT;
     return @seqs;
}

1;


