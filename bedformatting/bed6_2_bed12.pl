#!/usr/bin/perl -w

# <short description> 
#
# bed6_2_bed12.pl   copyright by Ruediger Lehmann (2020-07)
#
#
#-------------------------------------------------------------------------------
# Bug-Report
#
# <list of known bugs>
#
#-------------------------------------------------------------------------------

################################################################################
# GLOBAL VARIABLES, MODULES and LIBRARIES
################################################################################
use strict;

use IO::File;
use File::Basename;
use Getopt::Std;
use Cwd;

#use lib "<path of used library>";
#use <name of module>;

################################################################################
my( $file_bed6, $file_bed12, $newdir, $hash_size_bed6, $hash_size_bed12 );
my( %bed6, %bed12 );
#my( @<arrays> );
my( $ProgName, $DirName, $VERSION );
#-------------------------------------------------------------------------------
# last edition <year>-<month>-<day>

if( defined(readlink(__FILE__)) ){
	$ProgName = basename(readlink(__FILE__));
	$DirName  = dirname(readlink(__FILE__));
} else {
	$ProgName = basename(__FILE__);
	$DirName  = dirname(__FILE__);
}

$VERSION  = 0.1;

=pod

=head1 NAME

bed6_2_bed12.pl v0.1

This program changes a bed-File from bed6-format to bed12-format

=cut

################################################################################
#-------------------------------- Options --------------------------------------
#-------------------------------------------------------------------------------

use vars qw( $opt_i $opt_o );

getopts("i:o:");



$newdir = cwd();


#-------------------------------------------------------------------------------

if( $opt_i ){
	$file_bed6 = $opt_i;
} else {
    print "\n\t\tSorry, you have forgotten the name of the bed6-File \n\n";
    system "perldoc $0";
    exit 0;
}

#-------------------------------------------------------------------------------

if( $opt_o ){
	$file_bed12 = $opt_o;
} else {
    print "\n\t\tSorry, you have forgotten the name of the bed12-File \n\n";
    system "perldoc $0";
    exit 0;
}

#-------------------------------------------------------------------------------


################################################################################
#
#	    M A I N
#
################################################################################

$hash_size_bed6 = 0;


#-------------------------------------------------------------------------------

print "\n\n\t\t\t***************** READING BED6 FILE ***************** \n\n\n";
read_bed6_file( \%bed6, $file_bed6 );

$hash_size_bed6 = keys %bed6;
#printf "TEST2\t%d\t%d\n", $hash_size_bed6, $bed6{$hash_size_bed6-1}{start};

#-------------------------------------------------------------------------------

print "\n\n\t\t\t***************** CONDENSING TABLE ***************** \n\n\n";
condense_overlap( \%bed6, \%bed12, $hash_size_bed6 );



#-------------------------------------------------------------------------------


print "\n\n\t\t\t***************** WRITTING BED12 FILE ***************** \n\n\n";
$hash_size_bed12 = keys %bed12;
write_bed12_file( \%bed12, $file_bed12, $hash_size_bed12 );


#-------------------------------------------------------------------------------


################################################################################
################################################################################
# function

sub read_bed6_file{
	my( $rBed6, $filename ) = @_;
    my( $i );
    
	open( IN, "< $filename" ) || die "Cannot open $filename file!\n";

    print "\t-----------   reading bed6 file   -----------\n\n\n";

	$i = 0;
    
	while( <IN> ){                                           
	    chomp;

#print "----------- TEST 1 -------------\n";
        
		# Aufspalten der Zeile an den gesetzten Tabs in "Chromosom - Position1 - Position2 - Annotation"
		if( /^(\w+)\t(\d+)\t(\d+)\t(.+)$/ ){
			
#    print "$1\t$2\t$3\t$4\n";
			$$rBed6{$i}{chr}   = $1;
			$$rBed6{$i}{start} = $2;
			$$rBed6{$i}{stop}  = $3;
			$$rBed6{$i}{anno}  = $4;
		} 
#   	Manche Zeilen enthalten keine 4. Spalte, d.h. keine Annotation
		elsif ( /^(\w+)\t(\d+)\t(\d+)$/ ){

			$$rBed6{$i}{chr}   = $1;
			$$rBed6{$i}{start} = $2;
			$$rBed6{$i}{stop}  = $3;
			$$rBed6{$i}{anno}  = "";
		} 
		else {
		}

		# Inkrementieren des Zählers für den Hash
		$i++;
        
    }
#printf "TEST1\t%d\n", $$rBed6{200000}{start};

print "$i Eintraege in alter Tabelle.\n";
	close IN ;

}

#-------------------------------------------------------------------------------

sub condense_overlap{
	my( $rBed6, $rBed12, $hash_size ) = @_;
    my( $i, $j, $file_log, $comment );
    

    print "\t-----------   condensing table   -----------\n\n\n";

    $j = 0;
	
	$file_log = "bed6_2_bed12.log";
	
	$comment = "ANNOTATION";
	
	open( LOG, "> $file_log")  || die "Cannot open $file_log file!\n";

    $$rBed12{$j}{chr}   = $$rBed6{0}{chr};
    $$rBed12{$j}{start} = $$rBed6{0}{start};
    $$rBed12{$j}{stop}  = $$rBed6{0}{stop};
    $$rBed12{$j}{anno}  = $$rBed6{0}{anno};
    
	for ($i = 1; $i < $hash_size; $i++){
	
		if ($$rBed6{$i}{chr} eq $$rBed12{$j}{chr}){                                                      # Wir vergleichen die Positionen innerhalb eines Chromosoms
			
			if ($$rBed6{$i}{start} >= $$rBed12{$j}{start} && $$rBed6{$i}{stop} <= $$rBed12{$j}{stop}){   # der neue Bereich wird vom alten umschlossen
			
				if ($$rBed12{$j}{anno} ne $$rBed6{$i}{anno}){                                             # Es gibt Unterschiede in der Annotation bei überlappenden Bereichen
				
					printf LOG "%6d\t%s\t%d\t%d\t%s\n", $i, $$rBed12{$j}{chr}, $$rBed12{$j}{start}, $$rBed12{$j}{stop}, $$rBed12{$j}{anno};
					printf LOG "-------------- %s\n", $comment;
					printf LOG "%6d\t%s\t%d\t%d\t%s\n\n\n", $i+1, $$rBed6{$i}{chr}, $$rBed6{$i}{start}, $$rBed6{$i}{stop}, $$rBed6{$i}{anno};
					
					if ( length($$rBed6{$i}{anno}) > length($$rBed12{$j}{anno}) ){                       # Die neue Annotation ist umfangreicher
						$$rBed12{$j}{anno}  = $$rBed6{$i}{anno};
					}
				} else {
					printf LOG "%6d\t%s\t%d\t%d\t%s\n", $i, $$rBed12{$j}{chr}, $$rBed12{$j}{start}, $$rBed12{$j}{stop}, $$rBed12{$j}{anno};
					printf LOG "--------------\n";
					printf LOG "%6d\t%s\t%d\t%d\t%s\n\n\n", $i+1, $$rBed6{$i}{chr}, $$rBed6{$i}{start}, $$rBed6{$i}{stop}, $$rBed6{$i}{anno};
				}
			} elsif ($$rBed6{$i}{start} < $$rBed12{$j}{stop} && $$rBed6{$i}{stop} > $$rBed12{$j}{stop}){ # der neue Bereich überlappt sich mit dem alten
			
				if ($$rBed12{$j}{anno} ne $$rBed6{$i}{anno}){                                             # Es gibt Unterschiede in der Annotation bei überlappenden Bereichen
				
					printf LOG "%6d\t%s\t%d\t%d\t%s\n", $i, $$rBed12{$j}{chr}, $$rBed12{$j}{start}, $$rBed12{$j}{stop}, $$rBed12{$j}{anno};
					printf LOG "-------------- %s\n", $comment;
					printf LOG "%6d\t%s\t%d\t%d\t%s\n\n\n", $i+1, $$rBed6{$i}{chr}, $$rBed6{$i}{start}, $$rBed6{$i}{stop}, $$rBed6{$i}{anno};
					
					if ( length($$rBed6{$i}{anno}) > length($$rBed12{$j}{anno}) ){                       # Die neue Annotation ist umfangreicher
						$$rBed12{$j}{anno}  = $$rBed6{$i}{anno};
					}
				} else {
					printf LOG "%6d\t%s\t%d\t%d\t%s\n", $i, $$rBed12{$j}{chr}, $$rBed12{$j}{start}, $$rBed12{$j}{stop}, $$rBed12{$j}{anno};
					printf LOG "--------------\n";
					printf LOG "%6d\t%s\t%d\t%d\t%s\n\n\n", $i+1, $$rBed6{$i}{chr}, $$rBed6{$i}{start}, $$rBed6{$i}{stop}, $$rBed6{$i}{anno};
				}
				
				$$rBed12{$j}{stop}  = $$rBed6{$i}{stop};                                                 # Der Bereich wird erweitert
				
			} elsif ($$rBed6{$i}{start} > $$rBed12{$j}{stop}) {                                           # Es gibt keine Überlappung
			
				$j++;
				$$rBed12{$j}{chr}   = $$rBed6{$i}{chr};
				$$rBed12{$j}{start} = $$rBed6{$i}{start};
				$$rBed12{$j}{stop}  = $$rBed6{$i}{stop};
				$$rBed12{$j}{anno}  = $$rBed6{$i}{anno};
				
			}
			
		} else {		                                                                                 # Es beginnen die Eintraege eines neuen Chromosoms

			$j++;
			$$rBed12{$j}{chr}   = $$rBed6{$i}{chr};
			$$rBed12{$j}{start} = $$rBed6{$i}{start};
			$$rBed12{$j}{stop}  = $$rBed6{$i}{stop};
			$$rBed12{$j}{anno}  = $$rBed6{$i}{anno};
		}
		
	}

printf "Neue Tabelle:\t$j Eintraege\n";
    
	close LOG;

}

#-------------------------------------------------------------------------------

sub write_bed12_file{
	my( $rBed12, $filename, $hash_size ) = @_;
    my( $i );
    
 	open( OUT, "> $filename" ) || die "Cannot open $filename file!\n";

 	for ($i = 0; $i < $hash_size; $i++){
	
		if ($$rBed12{$i}{anno} ne ""){
			printf OUT "%s\t%d\t%d\t%s\n", $$rBed12{$i}{chr}, $$rBed12{$i}{start}, $$rBed12{$i}{stop}, $$rBed12{$i}{anno};
        } else {
			printf OUT "%s\t%d\t%d\t-\n",  $$rBed12{$i}{chr}, $$rBed12{$i}{start}, $$rBed12{$i}{stop};
        }		
	}

	close OUT;
}

################################################################################
################################################################################
################################################################################
################################################################################


=head1 SYNOPSIS

description of Synopsis

=head1 OUTPUT
	
output description

=cut
