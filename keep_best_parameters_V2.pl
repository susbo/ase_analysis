#!/usr/bin/perl
# Used to keep only the best alignment when the same fragment is aligned to two (e.g. parental) genomes.
# This script was developed by Susanne Bornelöv and Nima Rafati
#
#use strict;
use Getopt::Long;
use Pod::Usage;

my $ncpu = 2;
#&GetOptions('sample=s' =>\$sample, 'fai=s' =>\$genome, 'help' => \$help);

#pod2usage(1) if $help;
#pod2usage(-exitval => 0, -verbose => 2) if $usage;

if( !GetOptions( "help" => \$help,
		 "sample=s" =>\$sample,
		 "fai=s" =>\$genome,
		 "ncpu=i" =>\$ncpu)){
			pod2usage( { -message => "Failed to parse command line.\n",
	                 -verbose => 1,
        	         -exitval => 1 } );
}

if ($help){
 #|| @ARGV == "") {
        pod2usage( {-message => "$header\n",
                    -verbose => 99,
                    -exitval => 0 } );
}


#print $ARGV, '\n', $sample, '\n', $genome, '\n'; <STDIN>;
#if ($genome eq "" || $sample eq "" )
#{
#        exit print $usage;
#}

my $usage ="keep_best.pl sample_name reference.fa
written by susanne.bornelov\@gmail.com
modified/improved by nimarafati\@gmail.com\n";
my $verbose = 0;
#if ($ARGV[0] eq "" || $ARGV[1] eq ""){print $usage;exit;}
# Keep temporary files on SSD disk for faster read/write
my $tmp_prefix = "/tmp_fast";
#foreach my $file (@files) {
	chomp $file;

	$file = $sample;
	$reference = $genome;
	$file1= $file."-AA";
	$file2=$file."-BB";
	my $name1 = "$file-gsnap-snpT-AA";
	my $name2 = "$file-gsnap-snpT-BB";
	print "$name\n";

#	$name = "$file/$file.Reordered.sort.bami" ;
	my $shortName = $1;

	# Filter secondary alignment and unmapped read/mate
	print "\@CMD:samtools view -@ $ncpu  -F 268 $file1/$name1.Reordered.sort.bam | samtools view -bS -t $reference - >$file1/$name1.Reordered.sortP.bam\n"; 
	system "samtools view -@ $ncpu  -F 268 $file1/$name1.Reordered.sort.bam | samtools view -bS -t $reference - >$file1/$name1.Reordered.sortP.bam ";
	print "\@CMD:samtools view -@ $ncpu -F 268 $file2/$name2.Reordered.sort.bam | samtools view -bS -t $reference - >$file2/$name2.Reordered.sortP.bam\n";
	system "samtools view -@ $ncpu  -F 268 $file2/$name2.Reordered.sort.bam | samtools view -bS -t $reference - >$file2/$name2.Reordered.sortP.bam ";
	# Sort files by name
	print "\@CMD:samtools sort -@ $ncpu  -m 6G -n $file1/$name1.Reordered.sortP.bam -o $file1/$name1.Reordered.sortPName.bam\n";
	system "samtools sort -@ $ncpu  -m 6G -n $file1/$name1.Reordered.sortP.bam -o $file1/$name1.Reordered.sortPName.bam ";
	print "\@CMD:samtools sort -@ $ncpu  -m 6G -n $file2/$name2.Reordered.sortP.bam -o $file2/$name2.Reordered.sortPName.bam\n";
	system "samtools sort -@ $ncpu  -m 6G -n $file2/$name2.Reordered.sortP.bam -o $file2/$name2.Reordered.sortPName.bam ";

	my $id;
	print "\@CMD:samtools view -@ $ncpu $file1/$name1.Reordered.sortPName.bam|\n";
	open PAR1, "samtools view -@ $ncpu $file1/$name1.Reordered.sortPName.bam|";
	print "\@CMD:samtools view -@ $ncpu $file2/$name2.Reordered.sortPName.bam|\n";
	open PAR2, "samtools view -@ $ncpu $file2/$name2.Reordered.sortPName.bam|";

	# Output files, named by each sample
	open FIRST,">$file1/$name1.Reordered.sortPName.Best.sam";
	open SECOND,">$file2/$name2.Reordered.sortPName.Best.sam";
	open THIRD, ">$file2/$name2.Reordered.sortPName.Tie.sam";

	my $par1 = <PAR1>;
	my $par2 = <PAR2>;

	my $single = 0;	# Aligned only to one parent
	my $both = 0;	# Align to both
	my $tie = 0;	# Tie in alignment score

	my @single;
	my @both;
	my @tie;

	# Read and compare lines until EOF
	while (42) {
		# Print all remaining if one file has reached EOF
		if ($par1 eq "0") {
			print SECOND $par2;
			while ($par2 = <PAR2>) {
				print SECOND $par2;
			}
			last;
		} elsif ($par2 eq "0") {
			print FIRST $par1;
			while ($par1 = <PAR1>) {
				print FIRST $par1;
			}
			last;
		}

		# Split active line
		my @line1 = split "\t",$par1;
		my @line2 = split "\t",$par2;

		# Check if different names. If true, print the one "lowest" name and read next line.
		# This check makes sure that the same fragment is read on both files.
		if ($line1[0] ne $line2[0]) {
			$single++;	# statistics

			if ($line1[0] le $line2[0]) {
				$single[0]++;
				print FIRST $par1;
				next unless $par1 = <PAR1>//0;
				@line1 = split "\t",$par1;
			} else {
				$single[1]++;
				print SECOND $par2;
				next unless $par2 = <PAR2>//0;
				@line2 = split "\t",$par2;
			}
		}
		# If we have the same fragment in both files, print the best alignment.
		else {
			$both++;	# statistics

			# Save all rows corresponding to the correct fragment
			my @F = ($par1);
			my @S = ($par2);
			
			# Fragment name
			my $id = $line1[0];

			# Save all scores corresponding to the correct fragment (separately for First/Second file)
#			my $Fscore = $line1[4];
#			my $FN = 1;
#			my $Sscore = $line2[4];
#			my $SN = 1;

			# Compuate alignment score based on cigar
			my $Fcigar = 0;
			my $Scigar = 0;

			# Initialize with the first alignments of each fragment (First/Second)
			$Fcigar += cigar(@line1);
			$Scigar += cigar(@line2);

			# Read all qualities from first file corresponding to the correct fragment name ($id)
			if ($par1 = <PAR1>//0) {
				@line1 = split "\t",$par1;
				while ($line1[0] eq $id) {
					push @F,$par1;
#					$Fscore += $line1[4];
#					$FN++;
					$Fcigar += cigar(@line1);
					last unless $par1 = <PAR1>//0;
					@line1 = split "\t",$par1;
				}
			}

			# Read all qualities from second file corresponding to the correct fragment name ($id)
		       if ($par2 = <PAR2>//0) {
	         		@line2 = split "\t",$par2;
   	      			while ($line2[0] eq $id) {
			      	        push @S,$par2;
#					$Sscore += $line2[4];
#					$SN++;
					$Scigar += cigar(@line2);
		  	 	        last unless $par2 = <PAR2>//0;
		 		         @line2 = split "\t",$par2;
	         }
			}

			# Normalize by the number of aligned fragments per file
#			$Fscore /= $FN;
#			$Sscore /= $SN;

			if ($verbose) {
				if ($Fcigar != $Scigar) {
					print "CIGAR F: $Fcigar\tS: $Scigar ";
				}
			}

			# Break tie by randomization, assuming nothing can be worse than -10k
#			if ($Fcigar == $Scigar) {
#				if (rand()>0.5) { $Fcigar = -10000; $tie[0]++;}
#				else { $Scigar = -10000; $tie[1]++;}
#				$tie++;
			if ($Fcigar>$Scigar) {
				$both[0]++;
			} else {
				$both[1]++;
			}

			# Print all lines corresponding to the best alignment
			if ($Fcigar > $Scigar) {
				print FIRST @F;
			} elsif($fcigar < $Scigar) {
				print SECOND @S;
			} elsif ($Fcigar == $Scigar){
#				print "HERE\n 1:@F\n2:@S";<STDIN>;
				print THIRD @F;
			}		
#print "Single: @single\tBoth: @both\tTie: @tie\n" if $both[0]%10000==0;
		}			
	}
#	print "Single: $single\tBoth: $both\tTie: $tie\n";
	close PAR1;
	close PAR2;
	close FIRST;
	close SECOND;
	close THIRD;

#	exit();	# Uncomment to break after first file
#}

sub cigar {
	print join "\t",@_ if $verbose;
	# This should be improved since the column order may differ. However, the current implementation works with GSNAP.
	foreach my $ar (@_){
		if ($ar=~ m/MD.*/){
			$MD=$ar;
		}
	}
	
	my ($qual,$cigar) = ($_[10],$_[5]);

	if ($MD =~ /MD:Z:(.+)$/) {
		$MD = $1;
	} else {
		print "ERROR!\n";
		print "MD=$MD\n";
		print "qual=$qual\n";
		exit();
	}
	my $score = 0;
	my @qualities = split "",$qual;

	# If there is an insertion or a soft clip according to the cigar string
	# then the quality array has to be edited to delete the 
	# positions corresponding to the insertion.
	# Otherwise, the quality scores will be shifted.
	if ($cigar =~ /[IS]/) {
		my @cigar = split "",$cigar;
		my $i = 0;
		my $pos = 0;
		while ($i < @cigar) {
			my $num = $cigar[$i];
			while ($cigar[$i+1] =~ /\d/) {
				$num .= $cigar[$i+1];
				$i++;
			}
			if ($cigar[$i+1] eq "M") {
				$pos += $num;
			} elsif ($cigar[$i+1] eq "I" || $cigar[$i+1] eq "S") {
				print "".(join "",@qualities)."\n" if $verbose;
				splice @qualities,$pos,$num;
				print "".(join "",@qualities)."\n" if $verbose;
			}
			$i+=2;
		}
	}

	my @MD = split "",$MD;
	my $pos = 0;
	my $i = 0;

	print "$MD\n" if $verbose;

	# Read the MD tag and calculate alignment score.
	while ($i < @MD) {
		if ($MD[$i] eq "N" || $MD[$i] eq "n") { # "N" in reference. Do not count.
			print "N[pos=$pos|i=$i] " if $verbose;
			$pos++;
		}
		elsif ($MD[$i] =~ /[acgt]/) {	# SNP position using GSNAP. Count as correct match (+6).
			print "m[pos=$pos|i=$i] " if $verbose;
			$score += 6;
			$pos++;
		}
		elsif ($MD[$i] =~ /[ACGT]/) {	# Mismatch. Give a negative score corresponding to base quality.

			# Double-check if not acceptable?
#			my $coord = $pos+$_[3];
#			my @base = `samtools faidx reference/Galgal.$parent.fa $_[2]:$coord-$coord`; # must compure real position
#			print "samtools faidx reference/Galgal.$parent.fa $_[2]:$coord-$coord\n";
##>1:100-100
##C
#			chomp($base[1]);
#			my $strand = "+";
#			if ($_[1] & 16) { $strand = "-"; } # Check if reverse strand bit is set
#			print "Reference base: $base[1] ($strand) Expected: $MD[$i]\n";

			print "X[pos=$pos|i=$i] " if $verbose;
			$score -= ord($qualities[$pos])-33-1;	# add 6-(q+5), or equivalently 1-q
			$pos++;
		} elsif ($MD[$i] eq "^") {	# Deletion. Do not count.
			print "D[pos=$pos|i=$i] " if $verbose;
			while (defined($MD[$i+1]) && $MD[$i+1] =~ /[ACGTN]/) {
				$i++;
			}
		} elsif ($MD[$i] =~ /\d/) {	# Match. Give a +6 score corresponding to the length of the match.
			print "M[pos=$pos|i=$i] " if $verbose;
			my $num = $MD[$i];
			while (defined($MD[$i+1]) && $MD[$i+1] =~ /\d/) {
				$i++;
				$num .= $MD[$i];
			}
			$score += 6*$num;
			$pos += $num;
		} else {	# Unknown situation.
			print "ERROR\n";
		}
		$i++;
		print "score=$score\n" if $verbose;
	}
	print "\n" if $verbose;
	return $score;
}

__END__

=head1 NAME

keep_best_parameters.pl

=head1 DESCRIPTION

This script takes a sample name and genome_size or genome.fai.
Please note that the structure and naming of the file are important. Each individual should have two directories each with specific namming:
sample-AA should contain a bam file with following name "sample-AA/sample-gsnap-snpT-AA.Reordered.sort.bam"
sample-BB should contain a bam file with following name "sample-BB/sample-gsnap-snpT-BB.Reordered.sort.bam"

For instance sample "Sample_B16-T13-G":

Sample_B16-T13-G-AA/Sample_B16-T13-G-gsnap-snpT-AA.Reordered.sort.bam 
Sample_B16-T13-G-BB/Sample_B16-T13-G-gsnap-snpT-BB.Reordered.sort.bam

This script was developed by Susane Bornelöv and Nima Rafati
susanne.bornelov@gmail.com and nimarafati@gmail.com

=head1 SYNOPSIS

	keep_best_parameters.pl -sample -fai
	keep_best_parameters.pl -help

=head1 OPTIONS

=over  3

=item  B<-sample>

Sample name

=item  B<-fai>

genome.fai

=item  B<-cpu>

n cpu (default 2)

=item B<-h> or B<--help>

Display this helpful text.

=back

