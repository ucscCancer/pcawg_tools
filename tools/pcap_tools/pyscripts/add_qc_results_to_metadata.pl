#!/usr/bin/perl
use strict;
use JSON;
use Time::Piece;

my $analysisF = shift;
#my $qcF=shift;
my $bwa_output_dir = shift; #also includes "mode" tumor or normal? BAMStats bwa output bams are $mode/rgname.stats $mode/out_rgname.bam
my $stats_output_dir = shift;
my $download_timing = shift;
my $bammarkduplicates_metrics_file =  shift;
my $merged_timing = shift;

#`rsync -av $analysisF.old $analysisF`;
#`rm $analysisF.new`; 
 
open(IN,"<$analysisF");
open(OUT,">$analysisF.new");

#need to add the following:
#1) timing files:
	#a) docker_root/uuid_download_timing.txt
	#b) bwa_output_path_bwa_timing.txt
	#c) align_rg_bam_path*_merge_timing.txt
	#d) align_rg_bam_path*_qc_timing.txt
#2) bammarkduplciates metrics

my $DONE=undef;

while(my $line=<IN>)
{
  chomp($line);
  if($DONE || ($line !~ /<\/ANALYSIS_ATTRIBUTES>/i && $line !~ /<\/ANALYSIS>/i))
  {
	if($line =~ /bamsort\'+/)
	{
		$line =~ s/\'//g;
	}
	print OUT "$line\n";
	next;
  }

  # QC
  print OUT "      <ANALYSIS_ATTRIBUTE>
        <TAG>qc_metrics</TAG>
        <VALUE>" . &getQcResult($stats_output_dir) . "</VALUE>
      </ANALYSIS_ATTRIBUTE>\n";
  
  # Runtime
  print OUT "      <ANALYSIS_ATTRIBUTE>
          <TAG>timing_metrics</TAG>
          <VALUE>" . &getRuntimeInfo($bwa_output_dir,$stats_output_dir,$download_timing,$merged_timing) . "</VALUE>
      </ANALYSIS_ATTRIBUTE>\n";
  
  # Markduplicates metrics
  print OUT "      <ANALYSIS_ATTRIBUTE>
          <TAG>markduplicates_metrics</TAG>
          <VALUE>" . &getMarkduplicatesMetrics($bammarkduplicates_metrics_file) . "</VALUE>
      </ANALYSIS_ATTRIBUTE>\n";

  print OUT "$line\n";
  $DONE=1;
}
close(IN);
close(OUT);
`mv $analysisF $analysisF.old`;
`rsync -av $analysisF.new $analysisF`;

#original for multiple rg bams
sub getQcResult {
  my $dir = shift;
  # detect all the QC report files by checking file name pattern

  opendir(DIR, "$dir");

  my @qc_result_files = grep { /\.stats$/ } readdir(DIR);

  close(DIR);

  my $ret = { "qc_metrics" => [] };

  foreach (@qc_result_files) {

    open (QC, "< $dir/$_");

    my @header = split /\t/, <QC>;
    my @data = split /\t/, <QC>;
    chomp ((@header, @data));

    close (QC);

    my $qc_metrics = {};
    $qc_metrics->{$_} = shift @data for (@header);

    push @{ $ret->{qc_metrics} }, {"read_group_id" => $qc_metrics->{readgroup}, "metrics" => $qc_metrics};
  }

  return to_json $ret;
}

#sub getQcResult 
#{
#  # detect all the QC report files by checking file name pattern
#  my $qcF = shift;
#
#  my $ret = { "qc_metrics" => [] };
#  open (QC, "<$qcF");
#  my @header = split /\t/, <QC>;
#  chomp(@header);
#  while(my $line = <QC>)
#  {
#	chomp($line);
#    my @data = split /\t/, $line;
#    my $qc_metrics = {};
#    $qc_metrics->{$_} = shift @data for (@header);
#    push @{ $ret->{qc_metrics} }, ["read_group_id" => $qc_metrics->{readgroup}, "metrics" => $qc_metrics];
#  }    
#  close (QC);
#
#  return to_json $ret;
#}

sub read_timing {
  my ($file) = @_;
  open T, "<$file" or return "not_collected"; # very quick workaround to deal with no download_timing file generated due to skip gtdownload option. Brian, please handle it as you see it appropriate
  my $start = <T>;
  my $stop = <T>;
  chomp $start;
  chomp $stop;
  my $delta = $stop - $start;
  close T;
  return($delta);
}

sub getRuntimeInfo {
  my ($bwa_output_dir,$stats_output_dir,$download_timing_file,$merged_timing_file) = @_;
  # detect all the timing files by checking file name pattern, read QC data
  # to pull back the read group and associate with timing

  opendir(DIR, $stats_output_dir);
  #TODO: fix this naming
  my @qc_result_files = grep { /\.stats$/ } readdir(DIR);

  close(DIR);

  my $ret = { "timing_metrics" => [] };

  foreach (@qc_result_files) {
    my $rg = $_;
    # find the index number so we can match with timing info
  #TODO: fix this naming
    #$_ =~ /out_(\d+)\.bam\.stats\.txt/;
    #my $i = $1;

    open (QC, "< $stats_output_dir/$rg");

    my @header = split /\t/, <QC>;
    my @data = split /\t/, <QC>;
    chomp ((@header, @data));

    close (QC);

    my $qc_metrics = {};
    $qc_metrics->{$_} = shift @data for (@header);

    my $read_group = $qc_metrics->{readgroup};
    my @rg_name = split(/\./,$rg);
    pop(@rg_name);
    $read_group = join(".",@rg_name);
    # now go ahead and read that index file for timing
  #TODO: fix this naming
    my $download_timing = read_timing($download_timing_file);
    my $bwa_timing = read_timing("$bwa_output_dir/out_$read_group".".bam_bwa_timing.txt");
    my $qc_timing = read_timing("$bwa_output_dir/out_$read_group".".bam_qc_timing.txt");
    my $merge_timing = read_timing($merged_timing_file);

    # fill in the data structure
    push @{ $ret->{timing_metrics} }, { "read_group_id" => $read_group, "metrics" => { "download_timing_seconds" => $download_timing, "bwa_timing_seconds" => $bwa_timing, "qc_timing_seconds" => $qc_timing, "merge_timing_seconds" => $merge_timing } };

  }

  # and return hash
  return to_json $ret;

}

sub getMarkduplicatesMetrics {
  my $metrics_file_path = shift;
  my $dup_metrics = `cat $metrics_file_path`;
  my @rows = split /\n/, $dup_metrics;
  
  my @header = ();
  my @data = ();  
  my $data_row = 0;
  foreach (@rows) {
    last if (/^## HISTOGRAM/); # ignore everything with this and after
    next if (/^#/ || /^\s*$/);
    
    $data_row++;
    do {@header = split /\t/; next} if ($data_row == 1); # header line

    push @data, $_;
  }

  my $ret = {"markduplicates_metrics" => [], "note" => "The metrics are produced by bammarkduplicates tool of the Biobambam package"};
  foreach (@data) {
    my $metrics = {};
    my @fields = split /\t/;
    
    $metrics->{lc($_)} = shift @fields for (@header);
    delete $metrics->{'estimated_library_size'}; # this is irrelevant
    
    push @{ $ret->{"markduplicates_metrics"} }, {"library" => $metrics->{'library'}, "metrics" => $metrics};
    #push @{ $ret->{"markduplicates_metrics"} }, {"metrics" => $metrics};
  }

  return to_json $ret;
}
