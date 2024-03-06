use warnings;

# gene strand information
$file="/dcs04/lieber/statsgen/shizhong/database/gene_meta_hg38.txt";
open(IN, $file);
<IN>;
while(<IN>){
	chomp;
	@tokens=split(/\t/,$_);
	$strand{$tokens[0]} = $tokens[7];
}
close(IN);

$file="/dcs04/lieber/statsgen/jbenjami/projects/aanri_phase1/differential_analysis/tissue_comparison/summary_table/_m/BrainSeq_ancestry_4features_4regions.txt.gz";
open(IN, "gunzip -c $file |");
<IN>;
while(<IN>){
	chomp;
	@tokens=split(/\t/,$_);
	$tokens[0] =~ s/ //;
	if($tokens[9] ne "Gene"){
		next;
	}
	@temp = split(/\./,$tokens[1]);
	if(not defined $strand{$temp[0]}){
		next;
	}
	$grp=$tokens[0]."_".$tokens[9];
	$feature = $tokens[1];
	$sign{$grp}{$feature}=$tokens[8];
	$chr{$grp}{$feature}=$tokens[4];
	if($strand{$temp[0]} == 1){
		$p1{$grp}{$feature}=$tokens[5] - 500000;
		if($p1{$grp}{$feature} < 0){
			$p1{$grp}{$feature} = 0;
		}
		$p2{$grp}{$feature}=$tokens[6] + 50000;
	}
	else{
		$p1{$grp}{$feature}=$tokens[6] + 500000;		
		$p2{$grp}{$feature}=$tokens[5] - 50000;
		if($p2{$grp}{$feature} < 0){
			$p2{$grp}{$feature} = 0;
		}
	}
}
close(IN);
	
@grps = keys %sign;
foreach $grp (@grps){
	$out = $grp.".bed";
	$out1 = $grp."_up.bed";
	$out2 = $grp."_down.bed";
	@fs = keys %{$sign{$grp}};
	open(OUT1, ">hg38/$out1");
	open(OUT2, ">hg38/$out2");
	open(OUT, ">hg38/$out");
	foreach $f (@fs){
		 if($chr{$grp}{$f} ne "chrX" and $chr{$grp}{$f} ne "chrY"){
			if($sign{$grp}{$f} > 0){
				if($p1{$grp}{$f} < $p2{$grp}{$f}){
					print OUT1 "$chr{$grp}{$f}\t$p1{$grp}{$f}\t$p2{$grp}{$f}\n";
				}
				else{
					print OUT1 "$chr{$grp}{$f}\t$p2{$grp}{$f}\t$p1{$grp}{$f}\n";
				}
			}
			else{
				if($p1{$grp}{$f} < $p2{$grp}{$f}){
					print OUT2 "$chr{$grp}{$f}\t$p1{$grp}{$f}\t$p2{$grp}{$f}\n";
				}
				else{
					print OUT2 "$chr{$grp}{$f}\t$p2{$grp}{$f}\t$p1{$grp}{$f}\n";
				
				}
			}
			if($p1{$grp}{$f} < $p2{$grp}{$f}){
				print OUT "$chr{$grp}{$f}\t$p1{$grp}{$f}\t$p2{$grp}{$f}\n";
			}
			else{
				print OUT "$chr{$grp}{$f}\t$p2{$grp}{$f}\t$p1{$grp}{$f}\n";
			}
		 }
	}
	close(OUT1);
	close(OUT2);
	close(OUT);
}
