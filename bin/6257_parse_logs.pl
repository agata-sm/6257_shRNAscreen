#!c:/perl/bin/perl.exe

# script to parse logs from cutadapt and alignment filtering
# part of shRNA screen processing pipeline
# proj 6257

# logs in current directory are parsed for stats

# 21xii2022

use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use File::Basename;
use Data::Dumper;

my $script_name="6257_parse_logs.pl";


if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--outdir: /path/to/outdir where outfile log_stats.txt will be saved\n";

}else{

	my $parameters=join(' ', @ARGV);

	#commandline parsing for parameters
	GetOptions(
		'outdir=s'		=>	\(my $outdir)
	) or die "Error in command line arguments";

	my $outfile="$outdir\/log_stats.txt";

	my @cutadapt_fwd_logs = glob( './*cutadapt_trim_fwd.log' );
	my @cutadapt_rc_logs = glob( './*cutadapt_trim_rc.log' );
	my @alnfilt_logs = glob( './*read_stats.log' );

	print "\nfiles included in this log\n";
	print "@cutadapt_fwd_logs\n";
	print "@cutadapt_rc_logs\n";
	print "@alnfilt_logs\n\n";

	#HoH with stats and values
	my %smpl_stat_value;

	foreach my $file (@cutadapt_fwd_logs){
	    #print basename($file), "\n";
	    my $filename=basename($file);
	    if ($filename=~m/(\S+).cutadapt_trim_fwd.log/){
	    	my $smpl_id=$1;
	    	#print "$smpl_id\n";

	    	open(INFILE, "<","$file") or die "Cannot open input file $file: $!";

			while(<INFILE>){

				chomp $_;

				if($_=~m/^Total read pairs processed:\s+(.+)$/) { #Total read pairs processed:         11,275,495
					#print "$_\n";
					my $value=$1;
					$value=~s/,//g;
					#print "$value\n";
					$smpl_stat_value{$smpl_id}{"total_read_pairs"}=$value;
				}

				if($_=~m/^Pairs written \(passing filters\):\s+(.+)$/) { #Pairs written (passing filters):     4,506,668 (40.0%)
					#print "$_\n";
					my $strng=$1;
					$strng=~m/^(\S+).+/;
					my $value=$1;
					$value=~s/,//g;
					#print "$value\n";
					$smpl_stat_value{$smpl_id}{"fwd_passing_read_pairs"}=$value;
				}
			}
		close(INFILE);

		}
	}

	foreach my $file (@cutadapt_rc_logs){
	    #print basename($file), "\n";
	    my $filename=basename($file);

	    if ($filename=~m/(\S+).cutadapt_trim_rc.log/){
	    	my $smpl_id=$1;
	    	#print "$smpl_id\n";

	    	open(INFILE, "<","$file") or die "Cannot open input file $file: $!";

			while(<INFILE>){

				chomp $_;

				if($_=~m/^Pairs written \(passing filters\):\s+(.+)$/) { #Pairs written (passing filters):     4,506,668 (40.0%)
					#print "$_\n";
					my $strng=$1;
					$strng=~m/^(\S+).+/;
					my $value=$1;
					$value=~s/,//g;
					#print "$value\n";
					$smpl_stat_value{$smpl_id}{"rc_passing_read_pairs"}=$value;
				}
			}
		close(INFILE);
	    }
	}


	foreach my $file (@alnfilt_logs){
	    #print basename($file), "\n";
	    my $filename=basename($file);

	    if ($filename=~m/(\S+).read_stats.log/){
	    	my $smpl_id=$1;
	    	#print "$smpl_id\n";

	    	open(INFILE, "<","$file") or die "Cannot open input file $file: $!";

			while(<INFILE>){

				chomp $_;

				if($_=~m/^aligned read pairs:\s+(\d+)$/) { #Pairs written (passing filters):     4,506,668 (40.0%)
					#print "$_\n";
					my $value=$1;
					#print "$value\n";
					$smpl_stat_value{$smpl_id}{"aln_read_pairs_all"}=$value;
				}


				if($_=~m/^aligned read pairs with MAPQ 255:\s+(\d+)$/) { #Pairs written (passing filters):     4,506,668 (40.0%)
					#print "$_\n";
					my $value=$1;
					#print "$value\n";
					$smpl_stat_value{$smpl_id}{"aln_read_pairs_MAPQ255"}=$value;
				}


				if($_=~m/^aligned read pairs with MAPQ 255 filtered for NM \d+:\s+(\d+)$/) { #Pairs written (passing filters):     4,506,668 (40.0%)
					#print "$_\n";
					my $value=$1;
					#print "$value\n";
					$smpl_stat_value{$smpl_id}{"aln_read_pairs_MAPQ255_NM"}=$value;
				}

			}
			close(INFILE);
	    }
	}


	#print Dumper %smpl_stat_value;	



	open(OUTFILE, ">","$outfile") or die "Cannot open input file $outfile: $!";
	my $header="sample\tall_read_pairs\ttrimmed_read_pairs\talignments\talignments_mapq255\talignments_mapq255_mismatch";
	print OUTFILE "$header\n";

	#traverse the hash
	for my $sample (sort keys %smpl_stat_value){
    	#print "$sample: \n";
    
    	#for my $ele (keys %{$smpl_stat_value{$sample}}) {
        	
        	#print "  $ele: " . $smpl_stat_value{$sample}->{$ele} . "\n";


   		 	my $all_rps=$smpl_stat_value{$sample}{"total_read_pairs"};
   		 	#print "ALL proc rps: $all_rps\n";

   		 	my $fwd_rps=$smpl_stat_value{$sample}{"fwd_passing_read_pairs"};
   		 	#print "FWD rps: $fwd_rps\n";

   		 	my $rc_rps=$smpl_stat_value{$sample}{"rc_passing_read_pairs"};
   		 	#print "RC rps: $rc_rps\n";


   		 	my $tot_passing_rps=$fwd_rps+$rc_rps;
   		 	my $rounded_pass = sprintf("%.3f", $tot_passing_rps/$all_rps);

   		 	#print "TOT passing: $tot_passing_rps ($rounded_pass)\n";

   		 	my $all_aln=$smpl_stat_value{$sample}{"aln_read_pairs_all"};
   		 	my $mapq_aln=$smpl_stat_value{$sample}{"aln_read_pairs_MAPQ255"};
   		 	my $nm_aln=$smpl_stat_value{$sample}{"aln_read_pairs_MAPQ255_NM"};

   		 	my $rounded_frac_mapq = sprintf("%.3f", $mapq_aln/$all_aln);
   		 	my $rounded_frac_nm = sprintf("%.3f", $nm_aln/$all_aln);

   		 	#print "aln passing MAPQ $rounded_frac_mapq \t aln passing NM $rounded_frac_nm\n";



   		 #}
   		#print OUTFILE "$sample\t$all_rps\t$tot_passing_rps ($rounded_pass)\t$all_aln\t$mapq_aln ($rounded_frac_mapq)\t$nm_aln ($rounded_frac_nm)\n";
   		print OUTFILE "$sample\t$all_rps\t$tot_passing_rps\t$all_aln\t$mapq_aln\t$nm_aln\n";

   	}

close(OUTFILE);
### to do: sort hash, format as table (for easy reporting in R), add filtering stats

}
exit;
