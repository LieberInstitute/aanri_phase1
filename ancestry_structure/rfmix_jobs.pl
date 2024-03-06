use warnings;

open(OUT, ">rfmix_jobs.txt");
foreach $i (1..22){
	print OUT "/dcs04/lieber/statsgen/shizhong/software/RFMix/rfmix/rfmix ";
	print OUT "-f /dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/topmed/merge_H650_1M_2.5M_5M/EA_AA/all/vcf/chr$i.ea.aa.vcf.gz ";
	print OUT "-r /dcs04/lieber/statsgen/shizhong/database/1KG/GRCh38_phased_vcf/refmix_ref/out/1kGP_high_coverage_Illumina.chr$i.filtered.SNV_INDEL_SV_phased_panel.snpsOnly.eur.afr.vcf.gz ";
	print OUT "-m /dcs04/lieber/statsgen/shizhong/database/1KG/GRCh38_phased_vcf/refmix_ref/out/samples_id2 ";
	print OUT "-g /dcs04/lieber/statsgen/shizhong/database/1KG/GRCh38_phased_vcf/refmix_ref/out/genetic_map38 ";
	print OUT "-o chr$i ";
	print OUT "--chromosome=chr$i\n";
}
close(OUT);