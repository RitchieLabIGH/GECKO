#!/usr/bin/env perl
use warnings;
use strict;


my $state = 0;
my $line;
for(my $i = 0; $i<scalar(@ARGV); $i++){
    open(IN, $ARGV[$i]) or die();
    if($state==0){
        $line=<IN>;print $line;
        $line=<IN>;print $line;
        $state=1;
    }
    else{
        $line=<IN>;
        $line=<IN>;
    }
    while($line=<IN>){
        print $line;
    }
    close(IN);
}