#!/usr/bin/perl
use strict;
use Getopt::Long;

my $read1;
my $read2;
my $OUT1;
my $OUT2;
my $result;


$result = GetOptions("read1=s"=>\$read1,
		     "read2=s"=>\$read2,
		     "out1=s"=>\$OUT1,
		     "out2=s"=>\$OUT2);
#my $syncs = shift;

#open(IN,"<$syncs") || die "Error opening $syncs\n";


open(OUTFILES,">syncs.out") || die "Error opening syncs.out\n";
open(COUNTS,">syncs.counts") || die "Error opening syncs.counts\n";


my @pairedCleanFastqFiles = ();
my @cleanFastqFiles = ();
my @totalCleanReads = ();

# while(<IN>) {
#     my $file = $_;
#     chomp($file);
#     push(@pairedCleanFastqFiles,$file);
# }

push(@pairedCleanFastqFiles,$read1);
push(@pairedCleanFastqFiles,$read2);

foreach my $pairedFastq (@pairedCleanFastqFiles) {
    if ($pairedFastq =~ $read1) {
        my $mate = $read2;
#        $mate =~ s/_1_clean/_2_clean/;
#        my $sync1 = $pairedFastq;
#        $sync1 =~ s/\.fastq/\_sync\.fastq/;
#        my $sync2 = $mate;
#        $sync2 =~ s/\.fastq/_sync\.fastq/;
        
        open(IN1,"<$pairedFastq") || die "Error opening $pairedFastq\n";
        open(IN2,"<$mate") || die "Error opening $mate\n";
        
        my %in1 = ();
        my %in2 = ();
        
        my %final = ();
        
        my $id;
        my $i=1;
        while(<IN1>) {
            my $line = $_;
            chomp($line);
            if ($i==1) {
                $id = $line;
                $id =~ s/\/1//g;
                $id =~ s/1:N:0://g;
                $in1{$id} = 1;
            }
            elsif ( ( ($i-1) % 4) == 0) {
                $id = $line;
                $id =~ s/\/1//g;
                $id =~ s/1:N:0://g;
                $in1{$id} = 1;
            }
            $i++;
        }
        close(IN1);
        $id=0;
        $i=1;
        while(<IN2>) {
            my $line = $_;
            chomp($line);
            if ($i==1) {
                $id = $line;
                $id =~ s/\/2//g;
                $id =~ s/2:N:0://g;
                $in2{$id} = 1;
            }
            elsif ( ( ($i-1) % 4) == 0) {
                $id = $line;
                $id =~ s/\/2//g;
                $id =~ s/2:N:0://g;
                $in2{$id} = 1;
            }
        }
        close(IN2);
        foreach my $key (keys %in1) {
            if (defined($in2{$key})) {
                $final{$key} = 1;
            }
        }
        open(IN1,"<$pairedFastq") || die "Error opening $pairedFastq\n";
        open(IN2,"<$mate") || die "Error opening $mate\n";
        
        open(OUT1,">$OUT1") || die "Error opening $OUT1\n";
        open(OUT2,">$OUT2") || die "Error opening $OUT2\n";
        $id = 0;
        $i=1;
        while(<IN1>) {
            my $line = $_;
            chomp($line);
            if ($i==1) {
                $id = $line;
                $id =~ s/\/1//g;
                $id =~ s/1:N:0://g;
                if(defined($final{$id})) {
                    print OUT1 "$line\n";
                    my $next = <IN1>;
                    chomp($next);
                    print OUT1 "$next\n";
                    $next = <IN1>;
                    chomp($next);
                    print OUT1 "$next\n";
                    $next = <IN1>;
                    chomp($next);
                    print OUT1 "$next\n";
                }
            }
            elsif ( ( ($i-1) % 4) == 0) {
                $id = $line;
                $id =~ s/\/1//g;
                $id =~ s/1:N:0://g;
                if(defined($final{$id})) {
                    print OUT1 "$line\n";
                    my $next = <IN1>;
                    chomp($next);
                    print OUT1 "$next\n";
                    $next = <IN1>;
                    chomp($next);
                    print OUT1 "$next\n";
                    $next = <IN1>;
                    chomp($next);
                    print OUT1 "$next\n";
                }
            }
        }
        
        close(IN1);
        close(OUT1);
        $id = 0;
        $i=1;
        while(<IN2>) {
            my $line = $_;
            chomp($line);
            if ($i==1) {
                $id = $line;
                $id =~ s/\/2//g;
                $id =~ s/2:N:0://g;
                if(defined($final{$id})) {
                    print OUT2 "$line\n";
                    my $next = <IN2>;
                    chomp($next);
                    print OUT2 "$next\n";
                    $next = <IN2>;
                    chomp($next);
                    print OUT2 "$next\n";
                    $next = <IN2>;
                    chomp($next);
                    print OUT2 "$next\n";
                }
            }
            elsif ( ( ($i-1) % 4) == 0) {
                $id = $line;
                $id =~ s/\/2//g;
                $id =~ s/2:N:0://g;
                if(defined($final{$id})) {
                    print OUT2 "$line\n";
                    my $next = <IN2>;
                    chomp($next);
                    print OUT2 "$next\n";
                    $next = <IN2>;
                    chomp($next);
                    print OUT2 "$next\n";
                    $next = <IN2>;
                    chomp($next);
                    print OUT2 "$next\n";
                }
            }
        }
        close(IN2);
        close(OUT2);

        push(@cleanFastqFiles,$OUT1);
        push(@cleanFastqFiles,$OUT2);
        my $linecount = keys %final;
#        my $clean_lines = `wc -l $sync1`;
#        $clean_lines =~ s/^\s*|\s$sync1\n$//g;
#        push(@totalCleanReads,($clean_lines/4));
        push(@totalCleanReads,$linecount);
#        $clean_lines = `wc -l $sync2`;
#        $clean_lines =~ s/^\s*|\s$sync2\n$//g;
#        push(@totalCleanReads,($clean_lines/4));
        push(@totalCleanReads,$linecount);
    }
    else {
        next;
    }
}

foreach(@cleanFastqFiles) {
    print OUTFILES "$_\n";
}
close(OUTFILES);

foreach(@totalCleanReads) {
    print COUNTS "$_\n";
}
close(COUNTS);
