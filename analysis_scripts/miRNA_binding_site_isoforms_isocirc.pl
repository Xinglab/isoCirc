use strict;
my $dir = '/home/gaoy/genome_index/hg19_chromosomes/';
my $circ = '/mnt/isilon/xing_lab/aspera/isocirc_datasets/corr_long_reads/isoCirc_output/pb_3/isoCirc.out';
my $mir = '/mnt/isilon/xing_lab/aspera/isocirc_datasets/scripts/mature.fa';
#my @bases = ('A','G','C','T');
#my $mer6Sort_ref = &mer6Sort(\@bases);
#my @sorted_mer6 = @{$mer6Sort_ref};
my %mir_mers;
my $mir_seed_length = 6;
my %circRNA_info;
my %circRNA_seq;
my %titleNum;
my $pre_mir;
my $tag = 0;
open MIR, "<", $mir or die "cannot open $mir: $!";
while (<MIR>) {
	s/\r\n//;
	chomp;
	if (/^>(\S+)/ and index($_, 'Homo sapiens') > 1) {
		$pre_mir = $1;
		$tag = 1;
	} elsif (/^>/) {
		$tag = 0;
	} elsif($tag == 1) {
		my $seed = substr($_, 1, $mir_seed_length);
		my $RC_seed = &comp_rev_RNA($seed);
		push @{$mir_mers{$RC_seed}}, $pre_mir;
	}
}
print "Total 6mer seed number:",scalar(keys %mir_mers), "\n";

open CIRC, "<", $circ or die "cannot open $circ: $!";
while (<CIRC>) {
	chomp;
	my @line = split /\t/;
	if ($line[0] eq '#isoformID') {
		$titleNum{$line[$_]} = $_ for 0 .. $#line;
	} else {
		my $circRNA = $line[$titleNum{'chrom'}].':'.$line[$titleNum{'startCoor0base'}].'-'.$line[$titleNum{'endCoor'}];
		my $strand = substr($line[$titleNum{'canoBSJMotif'}], 0, 1);
		$circRNA_info{$line[$titleNum{'chrom'}]}{$circRNA} = [$line[$titleNum{'startCoor0base'}], $line[$titleNum{'endCoor'}], $strand, $line[$titleNum{'#isoformID'}], $line[$titleNum{'blockStarts'}], $line[$titleNum{'blockSize'}]];
	}
}

#print "rRNA\tcategory";
#for my $mer6(@sorted_mer6) {
#	print "\t$mer6";
#}
#print "\n";
my $n = 0;
my %mer6Dist_cat;
while (my ($chr, $circRNA_ref) = each %circRNA_info) {
	#next if $n > 0;
	next if length($chr) > 5 or $chr eq 'chrM';
	$n ++;
	open FA, "<", $dir.$chr.'.fa' or die "cannot open $dir$chr.fa: $!";
	my $seq;
	my $n = 0;
	while (<FA>) {
		chomp;
		$n++;
		next if $n==1;
		$seq .= $_;
	}
	print "$chr loaded: ", length($seq), "bp\n" if defined $seq;
	print scalar(keys %{$circRNA_ref}), " circRNAs in $chr:";
	#print "$_," for keys %{$circRNA_ref};
	print "\n";
	while (my ($circRNA, $info_ref) = each %{$circRNA_ref}) {
		my ($start0, $end0, $strand, $category, $block_starts, $block_sizes) = @{$info_ref};
		#next if $end - $start < 20;
		my @block_start_sort = split ',', $block_starts;
		my @block_size_sort = split ',', $block_sizes;
		for my $i (0 .. $#block_start_sort) {
			my $start = $block_start_sort[$i];
			my $length = $block_size_sort[$i];
			if($start0+$start+$length <= $end0 and $start >= 0 and $length > 0){
				#my $mer6Dist_ref = &mer6Finder(substr($seq, $start, $end-$start), $strand);
				my $block_seq = substr($seq, $start0+$start, $length);
				$block_seq = "\U$block_seq";
				#print "$rRNA\t$category";
				if ($strand eq '+') {
					$circRNA_seq{$circRNA} .= $block_seq;
				} elsif ($strand eq '-') {
					$circRNA_seq{$circRNA} = &comp_rev($block_seq).$circRNA_seq{$circRNA};
				} else {
					die "no strand: @{$info_ref}";
				}
				#print "\n";
			} else {
				die "error: @{$info_ref}";
			}
		}
		
	}
}

my %circRNA_mer_freq;
while (my ($circRNA, $circ_seq) = each %circRNA_seq) {
	my ($chr, undef) = split ':', $circRNA;
	#for my $mer (@sorted_mers) {
	while (my ($mer, $mir_ID_ref) = each %mir_mers) {
		my $index = -1;
		while (1) {
			$index = index($circ_seq.substr($circ_seq, 0, length($mir_seed_length)-1), $mer, $index+1);
			if ($index >= 0) {
				push @{$circRNA_mer_freq{$circRNA}{$mer}}, $index unless (exists $circRNA_mer_freq{$circRNA} and exists $circRNA_mer_freq{$circRNA}{$mer} and scalar(@{$circRNA_mer_freq{$circRNA}{$mer}}) >= 1 and $index - $circRNA_mer_freq{$circRNA}{$mer}[-1] < length($mer));
			}
			last if $index < 0;
		}
	}
	print ">${circRNA}_", length($circ_seq),"_$circRNA_info{$chr}{$circRNA}[3]_$circRNA_info{$chr}{$circRNA}[2]\n^$circ_seq\n";
	if (scalar(keys %{$circRNA_mer_freq{$circRNA}}) >= 1) {
		while (my ($mer, $pos_ref) = each %{$circRNA_mer_freq{$circRNA}}) {
			if (@{$pos_ref} < 3) {
				print "<";
				next;
			}
			print "<$mer(";
			print "$_/" for @{$mir_mers{$mer}};
			print "):";
			print "$_;" for @{$pos_ref};
			print ",";
		}
	} else {
		print "<0";
	}
	print "\n";
}


sub mer6Sort{
	my @single_base = @{$_[0]};
	my @mer6s;
	my @bi_base;
	for my $base1 (@single_base) {
		for my $base2 (@single_base) {
			push @bi_base, "$base1$base2";
		}
	}
	for my $base1 (@bi_base) {
		for my $base2 (@bi_base) {
			for my $base3 (@bi_base) {
				push @mer6s, "$base1$base2$base3";
			}
		}
	}

	if(@mer6s != 4**6){
		print scalar(@mer6s);
		die "bug for 6 mer categories";
	} else {
		\@mer6s;
	}
}

sub comp_rev {
	my $seq = reverse($_[0]);
	$seq =~ tr/ATCG/TAGC/;
	$seq;
}

sub comp_rev_RNA {
	my $seq = reverse($_[0]);
	$seq =~ tr/AUCG/TAGC/;
	$seq;
}