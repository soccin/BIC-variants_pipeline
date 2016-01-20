# Donavan Cheng
# Annotate with gene intervals list with exon coordinates, cytoBand

# perl ~/pubdata/hg19/scripts/annotate_IMPACT_gene_intervals_list.pl > /home/chengd1/LabPipeline/CopyNumberWork/gene_intervals.list.annotated

#use lib "/home/chengd1/Analysis/software/perl/lib";
#use lib qw(/ifs/data/zeng/dmp/resources/perllib);
use Tie::IxHash;

#open IN, "</home/chengd1/Analysis/pubdata/hg19/hg19_refFlat.canonical_coding_exons.2013Nov21st.interval_list";
#open IN, "</home/chengd1/Analysis/pubdata/hg19/hg19_refFlat.canonical_coding_exons.2013Nov27th.interval_list";
open IN, "</home/chengd1/pubdata/hg19/hg19_refFlat.canonical_coding_exons.2013Dec12th.interval_list";
my $line="";
tie(my %exonHash,'Tie::IxHash');
while(chomp($line=<IN>)){
    my @f = split('\t',$line);
    $f[0]=~s/chr//g;
    $exonHash{$f[0]}{$f[1]}=$line;
}
close IN;

open IN, "</dmp/data/pubdata/cytoband/VERSIONS/hg19_20090614/cytoBand.txt";
while(chomp($line=<IN>)){
    my @f = split('\t',$line);
    $f[0]=~s/chr//g;
    $cytHash{$f[0]}{$f[1]}=$line;
}
close IN;

my $count=1;
#open IN, "</dmp/data/mskdata/interval-lists/VERSIONS/cv3/gene_intervals.list";
open IN, "</dmp/data/mskdata/interval-lists/VERSIONS/cv1/gene_intervals.list";
while(chomp($line=<IN>)){
    my @f = split("\:",$line);
    my $chr = $f[0];
    my @g = split("\-",$f[1]);
    my $start = $g[0];
    my $stop = $g[1];
    
    my @hits=();
    print STDERR "$count\n";
    $count++;

    foreach my $key (sort {$a <=> $b} keys %{$exonHash{$chr}}){
	my @h = split('\t',$exonHash{$chr}{$key});
	$h[0]=~s/chr//g;
	if($start > $h[2]){next;}
	if($stop < $h[1]){last;}
	push @hits, $h[4].":".$h[0].":".$h[1].":".$h[2];
    }

    my $chit="-";
    foreach my $key (sort {$a <=> $b} keys %{$cytHash{$chr}}){
	my @h = split('\t',$cytHash{$chr}{$key});
	$h[0]=~s/chr//g;
	if($start > $h[2]){next;}
	if($stop < $h[1]){last;}
	$chit=$h[0].$h[3];
    }
    
    print "$line\t$chit\t".join(",",@hits)."\n";
}
close IN;
