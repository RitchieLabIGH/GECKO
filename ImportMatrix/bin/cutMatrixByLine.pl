#!/usr/bin/env perl
use warnings;
use strict;


my $fileMatrix = $ARGV[0];
my $numLines   = $ARGV[1];


my $line;
my $header;
my $nline = 0;
my $nfile = 1;

open(IN, $fileMatrix) or die();
$line=<IN>;$header = $line;
$line=<IN>;$header .= $line;

my $fileout = fillZero("submatrix", $nfile);
$fileout.="$nfile";
$fileout.=".matrix";
open(OUT, "> $fileout") or die();
print OUT $header;

while($line = <IN>){
    $nline++;
    print OUT $line;
    if($nline==$numLines){
        close(OUT);
        $nfile++;
        $fileout = fillZero("submatrix", $nfile);
        $fileout.="$nfile";
        $fileout.=".matrix";
        open(OUT, "> $fileout") or die();
        print OUT $header;
        $nline = 0;
    }
}
close(IN);
close(OUT);






sub fillZero{
    my $name = $_[0];
    my $number = $_[1];

    for(my $start=10; $start<1000000000; $start=$start*10){
        if($number<$start){
            $name.="0";
        }
    }

    return $name;

}