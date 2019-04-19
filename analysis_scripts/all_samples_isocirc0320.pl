use strict;
my $dir = shift;	#'/mnt/isilon/xing_lab/aspera/isocirc_datasets/corr_long_reads/isoCirc_output/'
my $pe = '/mnt/isilon/xing_lab/aspera/isocirc_datasets/CIRI_results/sample_BSJ_reads_PE_all.txt';
#my $out = '/mnt/isilon/xing_lab/aspera/isocirc_datasets/CIRI_results/sample_BSJ_reads_PE_all4overlap.txt';
my $num = 200;
my %titleNum;
my @libraries;
my %circ_library;
my %circ_bsj;
my %lib_total_read;
my %lib_total_circ;
opendir DIR, $dir or die "cannot opendir $dir: $!";
for my $folder (readdir DIR) {
	#print "$folder\n";
	if (-d $dir.$folder and length($folder)>=3 and -f $dir.$folder.'/isoCirc.out') {
		push @libraries, $folder;
		open CIRC, "<", $dir.$folder.'/isoCirc.out' or die "cannot open $dir$folder/isoCirc.out: $!";
		#print "$folder:";
		while (<CIRC>) {
			s/\r\n//;
			chomp;
			my @line = split /\t/;
			if ($line[0] eq '#isoformID') {
				$titleNum{$line[$_]} = $_ for 0 .. $#line;
			} elsif ($line[0] ne '#isoCirc'){
				$circ_library{ $line[$titleNum{'chrom'}].':'.($line[$titleNum{'startCoor0base'}]+1).'|'.$line[$titleNum{'endCoor'}] }{$folder} += $line[$titleNum{'readCount'}];
				$lib_total_read{$folder} += $line[$titleNum{'readCount'}];
				$lib_total_circ{$folder}{ $line[$titleNum{'chrom'}].':'.($line[$titleNum{'startCoor0base'}]+1).'|'.$line[$titleNum{'endCoor'}] } ++;
			}
		}
		#print scalar(keys %{$lib_total_circ{$folder}}),"\n";
	}
}

open PE, "<", $pe or die "cannot open $pe: $!";
while (<PE>) {
	chomp;
	my @line = split /\t/;
	next if $line[0] eq 'circRNA_ID' or @line < 6;
	$circ_bsj{$line[0]}{'CIRI_RNaseR'} = $line[3]+$line[4]+$line[6];
	$circ_bsj{$line[0]}{'CIRI_polyA'} = $line[1]+$line[2]+$line[5];
	#if ($line[0] =~ /(chr\d+):(\d+)\|(\d+)/ ) {
	#	$circ_loc{$line[0]} = [$1, $2, $3];
	#}
	$lib_total_read{'CIRI_RNaseR'} += $line[3]+$line[4]+$line[6];
	$lib_total_read{'CIRI_polyA'} += $line[1]+$line[2]+$line[5];
	$lib_total_circ{'CIRI_RNaseR'}{ $line[0] } ++ if $line[3]+$line[4]+$line[6] > 0;
	$lib_total_circ{'CIRI_polyA'}{ $line[0] } ++ if $line[1]+$line[2]+$line[5] > 0;
	print "CIRI\t",$line[3]+$line[4]+$line[6],"\t$line[0]\n" if $line[3]+$line[4]+$line[6] > 3*($line[1]+$line[2]+$line[5]);
}

my @sort_lib = sort {$a cmp $b} @libraries;
#print "library";
for my $lib (@sort_lib, 'CIRI_RNaseR', 'CIRI_polyA') {
	#print "\t$lib";
}

my @sort_circ = sort {$a cmp $b} keys %circ_library;
for my $circ (@sort_circ){

	#print "\n$circ";
	for my $lib (@sort_lib) {
		if (exists $circ_library{$circ}{$lib}){
			#print "\t$circ_library{$circ}{$lib}";
			print "$lib\t$circ_library{$circ}{$lib}\t$circ\n";
		} else {
			#print "\t0";
		}
	}
	for my $lib ('CIRI_RNaseR', 'CIRI_polyA') {
		if (exists $circ_bsj{$circ}{$lib}){
			#print "$lib\t$circ_bsj{$circ}{$lib}\t$circ\n";
			#print "\t$circ_bsj{$circ}{$lib}";
		} else {
			#print "\t0";
		}
	}
	#print "\n";
}
#print "\ntotal_reads";

for my $lib (@sort_lib, 'CIRI_RNaseR', 'CIRI_polyA') {
	#print "\t$lib_total_read{$lib}";
}
#print "\ntotal_circ";
for my $lib (@sort_lib, 'CIRI_RNaseR', 'CIRI_polyA') {
	#print "\t", scalar(keys %{$lib_total_circ{$lib}});
}
#print "\n";