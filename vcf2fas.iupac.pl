#load libraries ----
#use strict; # help you to keep track of variable names and scope


## arguments ----
my $vcf = $ARGV[0];
my $out = $ARGV[1]; # output file name


## define iupac codes ----
my %iupac = (
    "AA" => "A",
    "CC" => "C",
    "GG" => "G",
    "TT" => "T",
    "AG" => "R", "GA" => "R",
    "CT" => "Y", "TC" => "Y",
    "GC" => "S", "CG" => "S",
    "AT" => "W", "TA" => "W",
    "GT" => "K", "TG" => "K",
    "AC" => "M", "CA" => "M",
    "NN" => "N"
);


## read vcf file ---
my %fasta;
my %idmis;
my $allsites=0;
my $missing=0;
my %samples2numbers;
my $lastcol;
my $seqlen=0;

if ($vcf =~ /.gz$/) { open(file, "gunzip -c $vcf |") || die "can't open pipe to $vcf"; }
else { open(file, $vcf) or die "can't open $vcf $!\n"; }
while (my $line = <file>) {
    chomp($line);
    ## ignore comment lines
    if ($line =~ /^\#\#/) { next; }
    ## get header line (with sample names)
    elsif ($line =~ /^\#CHROM/) {
        my @cols = split("\t", $line);
        chomp @cols;
        $lastcol = scalar @cols - 1;
        # create key for individual in hash; 
        for my $x (9 .. $lastcol) {
            $samples2numbers{$x} = $cols[$x];
            $fasta{$cols[$x]} = "";
            $idmis{$cols[$x]} = 0;
        }
    }
    # get individuals genotypes
    else {
        my @cols = split("\t", $line);
        chomp @cols;
        my $CHR = $cols[0];
        my $POS = $cols[1];
        my $REF = $cols[3];
        my $ALT = $cols[4];
        my @ALTERN = split(",", $ALT);
        my $scalt = scalar @ALTERN;
        if ($scalt == 1) {$ALT1 = $ALTERN[0];}
        if ($scalt > 1 & $scal <= 2 ) {$ALT1 = $ALTERN[0]; $ALT2 = $ALTERN[1];}
        if ($scalt > 2) {next;}
        # each individual genotype
        $allsites++;
        my @GTline;
        for my $x (9 .. $lastcol) {
        my $sitemiss=0;
            my @GENO = split('\:', $cols[$x]);
            my $GT = $GENO[0];
            my $gen;
            # homozygous
            if    ( $GT eq "0/0" | $GT eq "0|0" ) {$gen = $REF . $REF;}
            elsif ( $GT eq "1/1" | $GT eq "1|1" ) {$gen = $ALT1 . $ALT1;}
            elsif ( $GT eq "2/2" | $GT eq "2|2" ) {$gen = $ALT2 . $ALT2;}
            elsif ( $GT eq "0") {$gen = $REF . $REF;}
            elsif ( $GT eq "1") {$gen = $ALT1 . $ALT1;}
            elsif ( $GT eq "2") {$gen = $ALT2 . $ALT2;}
            # heterozygous
            elsif ( $GT eq "0/1" | $GT eq "1/0" | $GT eq "0|1" | $GT eq "1|0" ) {$gen = $REF . $ALT1;}
            elsif ( $GT eq "0/2" | $GT eq "2/0" | $GT eq "0|2" | $GT eq "2|0" ) {$gen = $REF . $ALT2;}
            elsif ( $GT eq "1/2" | $GT eq "2/1" | $GT eq "1|2" | $GT eq "2|1" ) {$gen = $ALT1 . $ALT2;}
            # missing
            elsif ($GT eq "./." | $GT eq ".\/." ) {$gen = "NN"; $idmis{$samples2numbers{$x}}++; $missing++;}
            elsif ($GT eq ".")    {$gen = "NN"; $idmis{$samples2numbers{$x}}++; $missing++;}
            else { $gen == "NN"; }
            my $code = $iupac{$gen};
            push(@GTline, $code);
        }
        # 
        # my $ns = grep(/\-/, @GTline);
        # if ($ns <= $mis) {
        $seqlen++;
        for (my $i = 9; $i <= $lastcol; $i++) {
            my $index = $i - 9;
            my $ind = $samples2numbers{$i};
            $fasta{$ind} .= $GTline[$index];
        }
        # }
    }
}

## Transform to fasta format
open (my $fh, '>', $out) or die "can't open $out $!\n";
foreach my $k (sort keys %fasta) {
    print $fh ">$k\n";
    print $fh "$fasta{$k}\n";
}
close fh;
