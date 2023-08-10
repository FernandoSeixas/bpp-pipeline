#load libraries ----
use strict; # help you to keep track of variable names and scope


##### arguments /////////////////////////
my $fas = $ARGV[0];    # fasta file to read
my $out = $ARGV[1];    # final file to write to
my $loc = $ARGV[2];    # locus code
my $minlen = $ARGV[3]; # minimum alignment length


my $icount = 0;
my $seqlen = 0;
my %gphocs;

## read fas and filter individuals with too many missing sites
my $name;

if ($fas =~ /.gz$/) { open(file, "gunzip -c $fas |") || die "can't open pipe to $fas"; }
else { open(file, $fas) or die "can't open $fas $!\n"; }
while (my $line = <file>) {
    chomp($line);
    ##  header lines 
    if ($line =~ /^\>/) {
        $line =~ s/\>//g;
        $name = $line;
        $icount++;
        }
    # sequence lines
    else {
        $gphocs{$name} = $line;
        $seqlen = length($line);
        }
}


## Write to file in gphocs format
if ($seqlen >= $minlen) {
    open (my $fh, '>', $out) or die "can't open $out $!\n"; 
    # print alignment header
    print $fh "$icount $seqlen\n";
    # print name and sequence
    foreach my $k (sort keys %gphocs) {
        print $fh "$loc\^$k\t$gphocs{$k}\n";
    }
    close $fh;
}

