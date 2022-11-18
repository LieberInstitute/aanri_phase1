use warnings;

$outdir = "gwas";
mkdir $outdir;
chdir $outdir;

# read AA brnum
open(IN, "../out/p1");
$line = <IN>;
chomp $line;
@tokens = split(' ',$line);
shift @tokens;
shift @tokens;
shift @tokens;
shift @tokens;
foreach (@tokens){
	$tag{$_} = 1;
}
close(IN);

# ID dict
$infile = "/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/phenotype/pheno_PC";
open(IN, $infile);
open(OUT1, ">id_map");
open(OUT2, ">id_subset");
<IN>;
while(<IN>){
	@tokens = split(' ',$_);
	if(defined $tag{$tokens[1]}){
		print OUT1 "$tokens[0]\t$tokens[1]\n";
		print OUT2 "$tokens[0]\n";
	}
}
close(IN);
close(OUT1);
close(OUT2);

# subset gwas by chr
foreach $i (1..22){
	$pfile = "/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/topmed/merge_H650_1M_2.5M_5M/AA/all/plink/AA_chr".$i;
	`plink2 --pfile $pfile --keep id_subset --maf 0.01 --hwe 1e-5  --make-pgen --out chr$i`;
	# update sample names to brnum
	`plink2 --pfile chr$i --update-ids id_map --make-pgen --out chr${i}_brnum`;
}
