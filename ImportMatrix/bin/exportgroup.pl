#!/usr/bin/env perl
use warnings;
use strict;


my $fileMatrix = $ARGV[0];


open(IN, $fileMatrix) or die();
my $line1 = <IN>;
my $line2 = <IN>;
close(IN);

chomp($line1);
chomp($line2);
my @c1 = split("\t", $line1);
my @c2 = split("\t", $line2);

for(my $i = 1; $i<scalar(@c1); $i++){
    print $c1[$i],"\t",$c1[$i],"\t",$c2[$i],"\n";
}



