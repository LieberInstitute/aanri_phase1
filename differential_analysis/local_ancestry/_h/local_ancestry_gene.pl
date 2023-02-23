use warnings;
use Getopt::Long;

GetOptions ("chr=i"   => \$chr,    
            "start=i" => \$start,   
	    "end=i"   => \$end,			
            "out=s"   => \$out)   
or die("Error in command line arguments\n");

# extract data
$infile = "/dcs04/lieber/statsgen/shizhong/AANRI/local_ancestry/chr".$chr.".fb.tsv";
open(IN, $infile);
open(OUT, ">$out");
<IN>;
$line=<IN>;
print OUT "$line";
while(<IN>){
    @tokens=split(/\t/,$_);
    if($tokens[1] > $start and $tokens[1] < $end){
	print OUT "$_";
    }
    if($tokens[1] > $end){
	last;
    }
}
close(IN);
close(OUT);

# average across SNPs
`Rscript ../_h/local_ancestry.R $out ${out}_AFR`;

