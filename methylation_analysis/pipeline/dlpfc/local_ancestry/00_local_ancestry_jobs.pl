use warnings;

mkdir "out_local";
open(OUT, ">local_ancestry_jobs.txt");
open(IN, "../out/vmr.bed");
while(<IN>){
	chomp;
	@tokens=split(/\t/,$_);
	$out = $tokens[0]."_".$tokens[1]."_".$tokens[2];
	$tokens[0] =~ s/chr//;
	$start = $tokens[1] - 100000;
	$end = $tokens[2] + 100000;
	print OUT "perl local_ancestry.pl --chr $tokens[0] --start $start --end $end --out ./out_local/$out\n";
}
close(IN);
close(OUT);