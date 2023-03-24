use warnings;

$indir = "/dcs04/lieber/statsgen/shizhong/wgbs/finemap/hippo/pheno";
$outdir = "/dcs04/lieber/statsgen/shizhong/AANRI/VMR3/99/hippo/aa2/out";

open(OUT, ">residual_var_jobs.txt");

chdir $indir;
opendir($dh, "./");
@dirs = grep { !/^\./ && /chr/ && -d "$_"} readdir($dh);
closedir $dh;

foreach $d (@dirs){
	chdir $d;	
	opendir($dh, "./");
	@dirs2 = grep { !/^\./ && -d "$_"} readdir($dh);
	closedir $dh;
	foreach $d2 (@dirs2){
		chdir $d2;
		`mkdir -p $outdir/$d/$d2`;
		foreach $file (<p*>){
			$f1 = "$indir/$d/$d2/$file";
			$f2 = "$indir/$d/$d2/cgnames_$file";
			$f3 = "$outdir/$d/$d2/$file";
			print OUT "Rscript residual_var.R $f1 $f2 $d $f3\n";
		}
		chdir "..";
	}
	chdir "..";
}

close(OUT);