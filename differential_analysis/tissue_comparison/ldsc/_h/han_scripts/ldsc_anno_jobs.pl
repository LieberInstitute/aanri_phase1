use warnings;
use File::Basename;

$make_annot="/dcs04/lieber/statsgen/shizhong/ldsc/make_annot.py";
$referenceDir="/dcs04/lieber/statsgen/shizhong/marcc/jiyun/ldsc/referencefiles/1000G_EUR_Phase3_plink/";

open(OUT, ">ldsc_anno_jobs.txt");
foreach $bedfile (<bedfiles/*.bed>){
	$out = basename($bedfile);
	$out =~ s/\.bed//;
	$outdir="out_".$out;
	mkdir $outdir;
	foreach $i(1..22){
		$bimfile=$referenceDir."1000G.EUR.QC.".$i.".bim";
		print OUT "python $make_annot ";
		print OUT "--bed-file $bedfile ";
		print OUT "--bimfile $bimfile ";
		print OUT "--annot-file ./${outdir}/chr.${i}.annot.gz\n";
	}
}
close(OUT);