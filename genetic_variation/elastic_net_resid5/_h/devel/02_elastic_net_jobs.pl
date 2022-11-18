use warnings;

open(OUT, ">elastic_net_jobs.txt");
chdir "out";
foreach $file (<p*>){
	print OUT "Rscript elastic_net.R ./out/$file\n";
}
close(OUT);