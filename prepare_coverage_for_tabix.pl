use strict;

my $covfile=$ARGV[0];


open(COV,"$covfile") || die "$!";

my $header=<COV>;

$header=~s/Depth_for_//g;
$header=~s/^Locus\s+Total_Depth\s+Average_Depth_sample\s+//g;
$header="#chr\tpos\t".$header;
print "$header";
while(<COV>)
{
    my @line=split(/\t/);
    my $newline='';
    if($line[0]=~/^(\S+):(\d+)/)
    {
        $newline="$1\t$2\t";
    }else
    {
        print "NO MATCH\n";
    }
    $newline.= join ( "\t", @line[ 3 .. $#line ]);
    print "$newline";
}