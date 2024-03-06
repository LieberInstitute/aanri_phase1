use warnings;

mkdir "gwas";
chdir "gwas";

# all AIMS from CEU and AFR comparison
open(IN, "/dcl02/lieber/shan/shizhong/database/1KG/QC/CEU/freq_diff_CEU_AFR/AFR_aims.map");
while(<IN>){
	chomp;
	@tokens=split(' ',$_);
	$bar=$tokens[0]."-".$tokens[3];
	if($tokens[1] ne "rs116332642" and $tokens[1] ne "rs369717105"){ # remove duplicate snps in gwas
	$snp{$bar}=$tokens[1];
	}
}
close(IN);

# get AIMS overlap with gwas
$gwas="/dcl01/lieber/ImagingGenetics/Share/Imputation_TopMed/LIBD_Brain.topmed_102920.maf_0.005.info_0.1.call_0.8.cleaned_brnum.sex_imputed.ea_aa.geno_0.1.maf_0.05.hwe_0.000001";
$bimfile=$gwas.".bim";
open(IN, $bimfile);
open(OUT, ">aims");
open(OUT2, ">aims2");
open(OUT3, ">aims3");
while(<IN>){
	chomp;
	@tokens=split(' ',$_);
	$bar=$tokens[0]."-".$tokens[3];
	if(defined $snp{$bar}){
		print OUT "$tokens[1]\n";
		print OUT2 "$tokens[1]\t$snp{$bar}\n";
		print OUT3 "$snp{$bar}\n";
	}
}
close(IN);
close(OUT);
close(OUT2);
close(OUT3);

# get gwas data of our samples
`plink --bfile $gwas --extract aims --make-bed --out temp`;
`plink --bfile temp --update-name aims2 --recode --out aims_gwas`;
`rm temp*`;

# get gwas data of 1kg samples
$ceu="/dcl02/lieber/shan/shizhong/database/1KG/QC/CEU/freq_diff_CEU_AFR/CEU_aims";
$AFR="/dcl02/lieber/shan/shizhong/database/1KG/QC/CEU/freq_diff_CEU_AFR/AFR_aims";
`plink --file $ceu --extract aims3 --recode --out aims_ceu`;
`plink --file $AFR --extract aims3 --recode --out aims_AFR`;



