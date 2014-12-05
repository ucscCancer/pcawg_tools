#!/usr/bin/perl


use strict;
use Cwd;

die qq(
Bad numbr of inputs

) if(!@ARGV);

my $options ="";
my $normal="";
my $command="";
my $tumor="";
my $tumorbam = "";
my $output="";
my $snp="";
my $indel="";
my $working_dir = cwd();
my $log = '';

foreach my $input (@ARGV) 
{
	my @tmp = split "::", $input;
	if($tmp[0] eq "COMMAND") 
	{
		$command = $tmp[1];
	} 
	elsif($tmp[0] eq "NORMAL") 
	{
		$normal = $tmp[1];
	} 
	elsif($tmp[0] eq "TUMOR") 
	{
		$tumor = $tmp[1];
	}
	elsif($tmp[0] eq "TUMORBAM")
	{
		$tumorbam = $tmp[1];
	}
	elsif($tmp[0] eq "OPTION") 
	{
		my @p = split(/\s+/,$tmp[1]);
		if ($p[0] =~ m/validation|strand-filter|output-vcf/ && $p[1] == 0) {
			next;
		}
		$options = "$options ${tmp[1]}";
	}
	elsif($tmp[0] eq "OUTPUT") 
	{
		$output = $tmp[1];
	}
	elsif($tmp[0] eq "SNP") 
	{
		$snp = $tmp[1];
	}
	elsif($tmp[0] eq "INDEL") 
	{
		$indel = $tmp[1];
	}

	elsif($tmp[0] eq "LOG")
        {
                $log = $tmp[1];
        }
	
	else 
	{
		die("Unknown Input: $input\n");
	}
}

## VCF OUTPUT 
if ($output ne '') {
	system ("$command $normal $tumor $options --output-snp $working_dir/out.snp --output-indel $working_dir/out.indel 2>$log");
	my $indels = "$working_dir/out.indel.vcf";
	my $snps = "$working_dir/out.snp.vcf";
	
	system("grep -v '^\#' $indels | grep -v '^chrom position' >> $snps");
	my @chr_ord = chromosome_order($tumorbam);

	vs2vcf($snps, $output,\@chr_ord);
}
## SNP/INDEL OUTPUT
else {
	system  ("$command $normal $tumor $options --output-snp $snp --output-indel $indel 2>$log");
}

sub vs2vcf 
{

	#
	# G l o b a l     v a r i a b l e s 
	#
	my $version = '0.1';

	#
	# Read in file
	#
	my $input = shift;
	my $output = shift;
	my $chr_ord = shift;
	open(IN, $input) or die "Can't open $input': $!\n";
	open(OUT, ">$output") or die "Can't create $output': $!\n";
	my %output;

	while ( <IN> )
	{
		if ( /^#/ )
		{
			print OUT;
			next;
		}
		chomp;
		my $line = $_;

		my @flds = split ( "\t", $line );
		my $ref = $flds[3];
		my $alt = $flds[4];
		#
		# Deletion of bases
		#
		if ( $alt =~ /^\-/ )
		{
			($flds[3], $flds[4]) = ($ref.substr($alt,1), $ref);
		}

		#
		# Insertion of bases
		#
		if ( $alt =~ /^\+/ )
		{
			$flds[4] = $ref.substr($alt,1);
		}
		## insert dot for reference positions.
		if ($flds[4] eq '') {
			$flds[4] = '.';
		}
		print OUT join( "\t", @flds),"\n" unless defined $chr_ord;
		$output{$flds[0]}{$flds[1]} = join( "\t", @flds)."\n" if defined $chr_ord;
	}
	close(IN);
	# if chromosome order given return in sorted order
	if(defined $chr_ord) 
	{
		for my $chrom (@{ $chr_ord }) 
		{
			for my $pos (sort {$a<=>$b} keys %{ $output{$chrom} }) 
			{
				print OUT $output{$chrom}{$pos};
			}
		}
	}
	close(OUT);
}


sub chromosome_order 
{
	my $input = shift;
	# calculate flagstats
	my $COMM = "samtools view -H $input | grep '^\@SQ'";
	my @SQ = `$COMM`;
	chomp @SQ;
	for(my $i = 0; $i <= $#SQ; $i++) 
	{
		$SQ[$i] =~ s/^\@SQ\tSN:(.*?)\tLN:\d+$/$1/;
	} 
	return(@SQ);
}


