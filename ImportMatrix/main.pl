#!/usr/bin/env perl
use warnings;
use strict;


###########################################################
# main script for dealing with fastq files to create matrix
#
#
###########################################################


my %parameters;
$parameters{'reads'} = "/home/aubin/SeqOthello/example/fq/*{1,2}.fastq";
$parameters{'saveTrimmed'} = "";
$parameters{'notrim'} = "";
$parameters{'outdir'} = "./results";
$parameters{'singleEnd'} = "";
$parameters{'resume'} = "";
$parameters{'bam'} = 0;
$parameters{'kmersize'} = 30;


if(scalar(@ARGV)==0){
    help("");
    die();
}

############################################################
#### Reading parameters
$parameters{'function'} = $ARGV[0];
my $p = 1;
while( $p<scalar(@ARGV) ){
    if($ARGV[$p] eq '-h'){
        help($parameters{'function'});
        die();
    }

    if($ARGV[$p] eq '-resume'){
        $parameters{'resume'} = $ARGV[$p];
    }

    if($ARGV[$p] eq '-B'){
        $parameters{'bam'} = 1;
    }
    
    if($ARGV[$p] eq '--threshold'){
        $p++;
        $parameters{'threshold'} = $ARGV[$p];
    }    

    if($ARGV[$p] eq '--outdir'){
        $p++;
        $parameters{'outdir'} = $ARGV[$p];
    }

    if($ARGV[$p] eq '--reads'){
        $p++;
        $parameters{'reads'} = $ARGV[$p];
    }

    if($ARGV[$p] eq '--notrim'){
        $parameters{'notrim'} = $ARGV[$p];
    }

    if($ARGV[$p] eq '--saveTrimmed'){
        $parameters{'saveTrimmed'} = $ARGV[$p];
    }

    if($ARGV[$p] eq '--singleEnd'){
        $parameters{'singleEnd'} = '--singleEnd';
    }

    if($ARGV[$p] eq '--groupconfig'){
        $p++;
        $parameters{'groupconfig'} = $ARGV[$p];
    }

	if($ARGV[$p] eq '--kmersize'){
        $p++;
        $parameters{'kmersize'} = $ARGV[$p];
    }

    if($ARGV[$p] eq '--matrix'){
        $p++;
        $parameters{'matrix'} = $ARGV[$p];
    }

    if($ARGV[$p] eq '--matrixDiscrete'){
        $p++;
        $parameters{'matrixDiscrete'} = $ARGV[$p];
    }
    if($ARGV[$p] eq '--matrixRaw'){
        $p++;
        $parameters{'matrixRaw'} = $ARGV[$p];
    }

    $p++;
}



############################################################
#### decomposition
if($parameters{'function'} eq "decomposition"){
    if($parameters{'bam'}==0){
        my $command = "./nextflow run 01_decomposition.nf --reads '".$parameters{'reads'}."' --outdir \'".$parameters{'outdir'}."\' --kmersize ".$parameters{'kmersize'};
        $command.=" ".$parameters{'notrim'}." ".$parameters{'saveTrimmed'}." ".$parameters{'singleEnd'}." ".$parameters{'resume'};
        print $command,"\n";
        system($command) == 0
            or die "system failed: $?";
    }
    else{
        my $command = "./nextflow run 01_decomposition_bam.nf --reads '".$parameters{'reads'}."' --outdir \'".$parameters{'outdir'}."\' --kmersize ".$parameters{'kmersize'};
        $command.=" ".$parameters{'notrim'}." ".$parameters{'saveTrimmed'}." ".$parameters{'resume'};
        print $command,"\n";
        system($command) == 0
            or die "system failed: $?";
    }
}





############################################################
#### importation
if($parameters{'function'} eq "importation"){
    if(exists $parameters{'groupconfig'}){

        # initial importation
        my $command = "./nextflow run 02_importation.nf --groupconfig \'".$parameters{'groupconfig'}."\' --outdir \'".$parameters{'outdir'}."\'"." ".$parameters{'resume'};
        system($command) == 0
            or die "system failed: $?";
        my $directoryMatrix = $parameters{'outdir'}."/rawimport/matrix";
        my $listingparameters = $directoryMatrix."/list_param.txt";

        # recursive join
        my @list_matrix = glob("$directoryMatrix/*matrix");

        while(scalar(@list_matrix)>1){

            open(OUT, "> $listingparameters") or die();
            #get only the name without $directoryMatrix
            for(my $i = 0; $i<scalar(@list_matrix); $i++){
                my @content = split(/\//,$list_matrix[$i]);
                $list_matrix[$i] = $content[scalar(@content)-1];
            }

            #sort alpha
            my @list_matrix_Sorted = sort(@list_matrix);
            for(my $i = 0; $i<scalar(@list_matrix_Sorted); $i=$i+2){
                if($i<(scalar(@list_matrix_Sorted)-1)){
                    my $file1 = $list_matrix_Sorted[$i];
                    my $file2 = $list_matrix_Sorted[$i+1];
                    $file1=~s/subgroup//;$file1=~s/.conf.matrix//;
                    $file2=~s/subgroup//;$file2=~s/.conf.matrix//;
                    my $newfile = "subgroup".$file1.$file2.".conf.matrix";
                    my $testpath = $directoryMatrix."/".$list_matrix_Sorted[$i];
                    if ($testpath=~/^\//){
                        print OUT $newfile," ",$directoryMatrix."/".$list_matrix_Sorted[$i]," ",$directoryMatrix."/".$list_matrix_Sorted[$i+1],"\n";
                        #path absolu
                    }
                    else{
                        print OUT $newfile," ","../../../".$directoryMatrix."/".$list_matrix_Sorted[$i]," ","../../../".$directoryMatrix."/".$list_matrix_Sorted[$i+1],"\n";
                        #path relatif
                    }
                }
                else{
                    #move the file into new format
                    my $file1 = $list_matrix_Sorted[$i];
                    $file1=~s/subgroup//;$file1=~s/.conf.matrix//;
                    my $newfile = "subgroup".$file1.$file1.".conf.matrix";
                    $command = "cp $directoryMatrix"."/".$list_matrix_Sorted[$i]." $directoryMatrix"."/".$newfile;
                    system($command);
                }
            }
            close(OUT);

            $command = "./nextflow run 03_Join.nf --listParams \'".$listingparameters."\' --outdir \'".$parameters{'outdir'}."\'";
            system($command) == 0
                or die "system failed: $?";

            foreach my $elt (@list_matrix){
                $command = "rm $directoryMatrix"."/".$elt;
                system($command);
            }

            @list_matrix = glob("$directoryMatrix/*matrix");

        }


        my $newfile = $directoryMatrix."/RAWmatrix.matrix";
        $command = "mv ".$list_matrix[0]." $newfile";
        system($command);





    }
    else{
        help($parameters{'function'});
        die();
    }


}
############################################################
#### Anova Filter
if($parameters{'function'} eq "anova"){
    if(exists $parameters{'matrix'}){
        # initial importation
        my $command = "./nextflow run 04_anova.nf --matrix \'".$parameters{'matrix'}."\' --outdir \'".$parameters{'outdir'}."\'"." --threshold \'".$parameters{'threshold'}."\'"." ".$parameters{'resume'};
        system($command) == 0
            or die "system failed: $?";

    }
}

############################################################
#### Discretisation
if($parameters{'function'} eq "discretization"){
    if(exists $parameters{'matrix'}){
        # initial importation
        my $command = "./nextflow run 04_discretization.nf --matrix \'".$parameters{'matrix'}."\' --outdir \'".$parameters{'outdir'}."\'"." ".$parameters{'resume'};
        system($command) == 0
            or die "system failed: $?";

    }
}




############################################################
#### filter
if($parameters{'function'} eq "filter"){
    if(exists $parameters{'matrix'}){
        # search for most little group
        my %groups;
        open(IN, $parameters{'matrix'}) or die;
        my $line1 = <IN>;
        my $line2 = <IN>;chomp($line2);
        close(IN);
        my @contentline = split("\t",$line2);
        for(my $place = 1; $place<scalar(@contentline); $place++){
            if(exists $groups{$contentline[$place]}){
                $groups{$contentline[$place]}+=1;
            }
            else{
                $groups{$contentline[$place]} = 1;
            }
        }
        my @keys = sort { $groups{$a} <=> $groups{$b} } keys(%groups);
        my $littlevalue = $groups{$keys[0]};
        my $selectedvalue = int($littlevalue/3);


        # initial importation
	my $command = "./nextflow run 05_filter_CutByClass.nf --matrix \'".$parameters{'matrix'}."\' --filter $selectedvalue --outdir \'".$parameters{'outdir'}."\'"." ".$parameters{'resume'};
        system($command) == 0
            or die "system failed: $?";
    }
}



############################################################
#### get Real Counts
if($parameters{'function'} eq "real"){
     if(exists $parameters{'matrixDiscrete'} and exists $parameters{'matrixRaw'}){
        my $command = "./nextflow run 06_toreal.nf --matrixDiscrete \'".$parameters{'matrixDiscrete'}."\' --matrixReal \'".$parameters{'matrixRaw'}."\' --outdir \'".$parameters{'outdir'}."\'"." ".$parameters{'resume'};
        system($command) == 0
            or die "system failed: $?";
     }

}


sub help{
    my $function = $_[0];
}
