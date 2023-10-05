use warnings;

# trait meta
open(IN, "/dcs04/lieber/statsgen/shizhong/marcc/shizhong/database/ldscore/gwas/GWAS_statistics/GWAS_statistics2/trait_names_keys");
while(<IN>){
	chomp;
	@tokens=split(/\t/,$_);
	$filename=$tokens[0];
	$filename =~ s/\.gz/\.out\.results/;
	$type{$filename}=$tokens[1];
	$trait{$filename}=$tokens[2];
	print "$filename\t$trait{$filename}\n";
}
close(IN);

opendir($dh, "./");
@outdirs = grep { !/^\./ && /_2$/ && -d "$_"} readdir($dh);
closedir $dh;

open(OUT, ">ldsc_results.txt");
print OUT "group\ttrait\ttype\tProp._SNPs	Prop._h2	Prop._h2_std_error	Enrichment	Enrichment_std_error	Enrichment_p	Coefficient	Coefficient_std_error	Coefficient_z-score\n";
foreach $dir (@outdirs){
	$group=$dir;
	$group =~ s/out_//;
	$group =~ s/_2//;
	chdir $dir;
	foreach $file(<*results>){
		print OUT "$group\t$trait{$file}\t$type{$file}\t";
		open(IN,$file);
		<IN>;
		$line=<IN>;
		chomp $line;
		@tokens=split(/\t/,$line);
		shift @tokens;
		print OUT join("\t", @tokens), "\n";
		close(IN);
	}
	chdir "..";
}
close(OUT);