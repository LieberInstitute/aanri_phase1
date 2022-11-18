use warnings;

open(OUT, ">pred_jobs.txt");
chdir "out";
foreach $file (<p*>){
	if($file =~ /\d$/){
	print OUT "Rscript pred.R ./out/$file\n";
	}
}
close(OUT);