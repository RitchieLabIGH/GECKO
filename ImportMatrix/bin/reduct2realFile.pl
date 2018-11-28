#!/usr/bin/env perl
use warnings;
use strict;


###################################################################################
# Usage : perl script.pl file_reduct.matrix file_realcounts.matrix > newfile.matrix
#
#

my $reductfile = $ARGV[0];
my $initfile = $ARGV[1];

my $line;
my @content;

my %dict;
open(IN, $reductfile) or die();
my $header1 = <IN>;
my $header2 = <IN>;
while($line=<IN>){
    @content = split("\t",$line);
    $dict{$content[0]} = 1;
}
close(IN);

open(IN, $initfile) or die();
$header1 = <IN>;
$header2 = <IN>;
print $header1,$header2;
while($line=<IN>){
    @content = split("\t",$line);
    if(exists($dict{$content[0]})){
        print $line;
    }
}
close(IN);


