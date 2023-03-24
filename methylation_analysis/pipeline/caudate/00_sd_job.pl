use warnings;

mkdir "out";
open(OUT, ">sd_jobs.txt");
foreach $i (1..22){
	print OUT "Rscript sd.R chr$i\n";
}
close(OUT);