#/usr/bin/perl -w
use strict;

my ($Len_seq, $No_seq, $id, %hash,$con,@arr, @brr,$i,$j,$k1,$k2,$k,$sum,$start,$Island_len,$line1,$end,$Sea_len,$num,%id,$j0,$j1,$j2, $No_pairs,%cs);
my (@seq1,@seq2,$seq1,$seq2,$count_snp,@nseq1,@nseq2,$count_indel,$PIC,$len);

open (FH, "$ARGV[0]") or die;   ## Input file: cluster result of multi-cp sequences in the foemat of fq
open (CR, ">ClusterResult.1line") or die;

$No_seq=0;

while (<FH>) {
	chomp;
	if (/^>/) {
		$No_seq ++;
		if ($No_seq == 1) {
			print CR "$_\n";
		}else{
			print CR "\n"."$_\n";		
		}
	}else{
		print CR "$_";
	}
}
print CR "\n";

close FH;
close CR;

open (CR, "ClusterResult.1line") or die;
open (TMP, ">tmp") or die;

while (<CR>) {
	chomp;
	if (/^>/) {
		$id=$_;
	}else{
		$Len_seq=length ($_);
		$hash{$id}=$_;
		$con=$_;
		@arr=split("");
		push (@brr, @arr);
	}
}
close FH;

for ($j=0;$j<$Len_seq;$j++) {
	$k1=$j;
	print TMP "$k1";
	for ($i=1;$i<$No_seq;$i++) {
		$k2=$j+$Len_seq*$i;
		if ($brr[$k1] eq $brr[$k2]) {
			print TMP "\t0";
		}else{
			print TMP "\t1";
		}
	}
	print TMP "\n";
	$k1=$k2;
}
close TMP;

open (TMP, "tmp") or die;
open (ST, ">star") or die;
while (<TMP>) {
	chomp;
	@arr = split "\t";
	$sum = 0;
	for ($i=1;$i<$No_seq;$i++) {
		$sum=$sum+$arr[$i];
	}
	if ($sum == 0){
		print ST "1";
	}else{
		print ST "0";
	}
}
print ST "\n";

close TMP;
close ST;

open (ST, "star") or die;
open (CI, ">Con_Island.txt") or die;

while (<ST>) {
	chomp;
	@arr=split("");
}

$i=0;
$k=0;
$Island_len=$ARGV[1];
for ($j=0;$j<$Len_seq;$j++) {
	chomp;
	if ($arr[$j] eq 1) {
		$i++;
	}else{
		if ($i >= $ARGV[1]) {
			$k++;
			$start=$j-$i+1;
			print CI "$k\t"."$start\t"."$j\t"."$i\n";
		}
		$i=0;
	}
}

close ST;
close CI;

open (CI, "Con_Island.txt") or die;
open (CS, ">Con_Sea.txt") or die;
open (CIS, ">Con_Island.seq.txt") or die;

$i=0;
while (<CI>) {
	chomp;
	$i++;
	@arr=split("\t");
	if ($i eq 1) {
		$line1=$arr[1]-1;
		$end=$arr[2]+1;
		print CS "$i\t"."$end\t";
		print CIS "$_\t".substr($con,$arr[1]-1,$arr[3])."\n";
	}else{
		$start=$arr[1]-1;
		$Sea_len=$start-$end+1;
		print CS "$start\t"."$Sea_len\n";
		$end=$arr[2]+1;
		print CS "$i\t"."$end\t";
		print CIS "$_\t".substr($con,$arr[1]-1,$arr[3])."\n";
	}
}
$Sea_len=$Len_seq-$end+1+$line1;
print CS "$line1\t"."$Sea_len\n";
$num = $i;

close CI;
close CS;
close CIS;

open (CS, "Con_Sea.txt") or die;
open (CSS, ">Con_Sea.seq.txt") or die;

while (<CS>) {
	chomp;
	@arr=split("\t");
	$cs{$arr[0]}=$_;
	foreach my $key (keys %hash) {
		if ($arr[2] >= $arr[1]) {
			print CSS "$arr[0]."."$key\t".substr($hash{$key},$arr[1]-1,$arr[3])."\n";
		}else{
			print CSS "$arr[0]."."$key\t".substr($hash{$key},$arr[1]-1,$Len_seq-$arr[1]+2).substr($hash{$key},0,$line1)."\n";
		}
	}
}
close CS;
close CSS;

open (CSS, "Con_Sea.seq.txt") or die;
open (OUT, ">Con_Sea.PIC.txt") or die;

$i=0;
while (<CSS>) {
	chomp;
	if ($i<$No_seq) {
		@arr=split("\t");
		@brr=split(/>/,$arr[0]);
		$id{$i}=$brr[1];
		$i++;
	}
}

$No_pairs=$i*($i-1)/2;

print OUT "Con_Sea\t\t\t\t";
for ($k1=0;$k1<$No_seq-1;$k1++) {
	for ($k2=$k1+1;$k2<$No_seq;$k2++) {
		print OUT "$id{$k1}"." vs "."$id{$k2}\t\t\t";
	}
}

print OUT "\n"."No.\t"."Start site\t"."End site\t"."Length";
print OUT ("\tSNP"."\tIndel"."\tPIC") x $No_pairs;
print OUT "\n";

close CSS;

open (CSS, "Con_Sea.seq.txt") or die;
$i=0;
while (<CSS>) {
	chomp;
	@arr=split("\t");
	$hash{$i}=$arr[1];
	$i++;
}
close CSS;

for ($j0=0;$j0<$num;$j0++) {
	$id=$j0+1;
	print OUT $cs{$id};
	for ($k1=0;$k1<$No_seq-1;$k1++) {
		for ($k2=$k1+1;$k2<$No_seq;$k2++) {
			$j1=$j0*$No_seq+$k1;
			$j2=$j0*$No_seq+$k2;
			$seq1=$hash{$j1};
			$seq2=$hash{$j2};
			$PIC=&SNP_gbq+&Indel_gbq;
			print OUT "\t".&SNP_gbq."\t".&Indel_gbq."\t".$PIC;
		}
	}
	print OUT "\n";
}

close OUT;

sub SNP_gbq {
	@seq1 = split ("",$seq1);
	@seq2 = split ("",$seq2);
	$count_snp = 0;
	$len = @seq1;
	for ($i=0;$i<$len;$i++) {
		if ($seq1[$i] ne "-" and $seq2[$i] ne "-" and $seq1[$i] ne $seq2[$i]) {
			$count_snp ++;
		}
	}
	$count_snp;
}

sub Indel_gbq {
	@seq1 = split ("",$seq1);
	@seq2 = split ("",$seq2);
	$count_indel = 0;
	$len = @seq1;
	$i=0;
	$j=0;
	while ($i<$len) {
		if ($seq1[$i] eq "-" and $seq2[$i] eq "-") {
			$i ++;
		}else{
			$nseq1[$j]=$seq1[$i];
			$nseq2[$j]=$seq2[$i];
			$j++;$i++
		}
	}
	for ($k=0;$k<=$j;$k++) {
		if ($nseq1[$k] ne "-" and $nseq1[$k+1] eq "-") {
			$count_indel ++;
		}
	}
	for ($k=0;$k<=$j;$k++) {
		if ($nseq2[$k] ne "-" and $nseq2[$k+1] eq "-") {
			$count_indel ++;
		}
	}
	$count_indel;
}
