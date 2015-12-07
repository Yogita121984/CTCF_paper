#!/usr/bin/perl -w
open(RH2,"GSE27463_RefSeq.reg.threecolumn");
while($line=<RH2>)
{
chomp($line);
@split_line=split(/\t/,$line);
$gname=$split_line[0];
$rpkm_0=$split_line[1];
$rpkm_160=$split_line[2];
#$fold=2;
$fold_try=$rpkm_0 * 2;
#$gname=$split_line[4];

if($rpkm_160 >= $fold_try)
{
 #   print "$chr\t$start\t$end\t$strand\t$gname\t$rpkm_160\n";
print "$gname\t$rpkm_0\t$rpkm_160\n";

}
}
exit;
