use warnings;

open(OUT, ">vmr_level_jobs.txt");
foreach $i (1..22){
	print OUT "Rscript vmr_level.R chr$i\n";
}
close(OUT);