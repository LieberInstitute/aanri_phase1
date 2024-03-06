use warnings;
use File::Basename;

$liftover="/dcs04/lieber/statsgen/shizhong/scripts/liftover/liftover_bed.pl";
foreach $file(<./hg38/*.bed>){
	$out = basename($file);
	`perl $liftover --in $file --38to19`;
	`mv lifted.bed $out`;
}

