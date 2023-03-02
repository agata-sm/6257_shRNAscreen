#!c:/perl/bin/perl.exe

# script to parse count table
# part of shRNA screen processing pipeline
# proj 6257

# add gene names, format columns and colnames

# 22xii2022

use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use List::Util 'all';

my $script_name="6257_proc_cnt_table.pl";


if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/infile\n";
	print "--library: /path/to/library info file\n";
	print "--outfile: /path/to/outfile\n";
	print "--setup: one of: all - all shRNAs retained; 0rm - rows with 0 counts across ALL samples removed; 0rm_noAlt - rows with 0 counts and shRNAs with annotated off-targets removed\n";

}else{

	my $parameters=join(' ', @ARGV);

	#commandline parsing for parameters
	GetOptions(
		'infile=s'		=>	\(my $infile),
		'outfile=s'		=>	\(my $outfile),
		'library=s'		=>	\(my $library),
		'setup=s'		=>	\(my $setup)
	) or die "Error in command line arguments";


	my %geneID_probe;
	my %probes_alt;
	open(INFILE_LIB, "<","$library") or die "Cannot open input file $library: $!";
	
	while(<INFILE_LIB>){

		chomp $_;

		unless ($_=~m/cloneId/) {

			my @line=split /\t/;

			my $probe_id=$line[0];
			my $gene_id=$line[5];

			#print "$probe_id\t$gene_id\n";
			$geneID_probe{$probe_id}=$gene_id;

			if ($line[10] ne qw /noAltHuman/){
				$probes_alt{$probe_id}=$line[10];
				#print "@line\n";
			}
		}

	}
	close(INFILE_LIB);


	#print Dumper %geneID_probe;	


	open(INFILE, "<","$infile") or die "Cannot open input file $infile: $!";
	open(OUTFILE, ">","$outfile") or die "Cannot open output file $outfile: $!";


	while(<INFILE>){

		chomp $_;

		if  ($_ =~m/^Geneid/){

			my @line=split /\t/;
	        my $probe_id=$line[0];
			my @other=splice @line,0,6;
			#print "@line\n"; #SMPL.mapped.filt_mapq255_NM1.bowtie2.bam
			s/\.mapped\.filt_mapq255_NM\d+\.bowtie2\.bam$// foreach @line;
			#print "@line\n"; #SMPL
			my $samples=join("\t",@line);
			my $header="ProbeID\tGeneID\t$samples";
			print OUTFILE "$header\n";
		}


		unless ( ($_=~m/^#/) | ($_ =~m/^Geneid/) ) {

			my @line=split /\t/;
			my $probe_id=$line[0];
			my @other=splice @line,0,6;
			#print "@line\n";

			my $gene_id=$geneID_probe{$probe_id};
			my $samples=join("\t",@line);


			if($setup eq qw /0rm/){

				unless (all { $_ == 0 } @line ){
					print OUTFILE "$probe_id\t$gene_id\t$samples\n";
				}

			}elsif($setup eq qw /0rm_noAlt/){

				unless (all { $_ == 0 } @line ){
					
					unless (exists $probes_alt{$probe_id} ){
						print OUTFILE "$probe_id\t$gene_id\t$samples\n";
					}
				}
			}			
			else{
				print OUTFILE "$probe_id\t$gene_id\t$samples\n";
			}


		}
	}
	close(INFILE);
	close(OUTFILE);

}
exit;

