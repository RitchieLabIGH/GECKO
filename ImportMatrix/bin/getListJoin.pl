#!/usr/bin/env perl
use warnings;
use strict;

my @listOfFiles =();
foreach my $elt (@ARGV){
    push(@listOfFiles, $elt);
}

my @listOfFilesSorted = sort(@listOfFiles);
for(my $i = 0; $i<scalar(@listOfFilesSorted); $i=$i+2){
    if($i<(scalar(@listOfFilesSorted)-1)){
        my $file1 = $listOfFilesSorted[$i];
        my $file2 = $listOfFilesSorted[$i+1];
        $file1=~s/subgroup//;$file1=~s/.conf.matrix//;
        $file2=~s/subgroup//;$file2=~s/.conf.matrix//;
        my $newfile = "subgroup".$file1.$file2.".conf.matrix";
        print $newfile,",",$listOfFilesSorted[$i],",",$listOfFilesSorted[$i+1],"\n";
    }
    else{
        my $file1 = $listOfFilesSorted[$i];
        $file1=~s/subgroup//;$file1=~s/.conf.matrix//;
        my $newfile = "subgroup".$file1.$file1.".conf.matrix";
        print $newfile,",",$listOfFilesSorted[$i],",",$listOfFilesSorted[$i],"\n";
    }
}
