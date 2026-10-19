#!/usr/bin/perl

if ($#ARGV != 2){
    die << "EOT" ;
This program is aimed at transforming an XSCALE output file
(with FRIEDEL'S_LAW=FALSE) in a mtz file (CCP4 reflections
file)

The program should be used like this:
 xds2mtz.pl infile outfile label

the label is the label prepended to the column names:
F_label(+), SIGF_label(+), F_label(-) & SIGF_label(-) ... etc

EOT
    ;

}

open (IN,"$ARGV[0]");
#open (OUT,">${ARGV[1]}.hkl");
$label = $ARGV[2];

# Extracting the information from the header of XSCALE output 

while (($line = <IN>) !~ /\!END_OF_HEADER/){
    if ($line =~ /\!UNIT_CELL_CONSTANTS=[\s]+([0-9]+.[0-9]*)[\s]+([0-9]+.[0-9]*)[\s]+([0-9]+.[0-9]*)[\s]+([0-9]+.[0-9]*)[\s]+([0-9]+.[0-9]*)[\s]+([0-9]+.[0-9]*)/){
        $CELL = "$1 $2 $3 $4 $5 $6";
        print "CELL $CELL\n";
    }
    if ($line =~ /\!SPACE_GROUP_NUMBER=[\s]+([0-9]+)/){
        $SYM = $1;
        print "SPACE GROUP NUMBER $SYM\n";
    }
};

#while (($line = <IN>) !~ /END_OF_DATA/){
#    chomp $line;
#
#    print OUT ("$line\n");
#}

close (IN);
#close (OUT);

# Generating the command file for XDSCONV

open (COM,">XDSCONV.INP");
print COM << "EOT";

SPACE_GROUP_NUMBER= $SYM
UNIT_CELL_CONSTANTS= $CELL

!name, type and resolution window
INPUT_FILE= ${ARGV[0]} XDS_ASCII 100 1.0

!name, target, freidel's law validity
 OUTPUT_FILE= ${ARGV[0]}_conv FALL FRIEDEL'S_LAW= FALSE
!

EOT
;

# Run xdsconv and add M/ISYM column
system "xdsconv > xdsconv.log";

open (IN,"${ARGV[0]}_conv");
open (OUT,">${ARGV[0]}_conv2");

while (<IN>){
    chomp;
    split;
    $Fp = 1;
    $Fm = 1;
    $NL = "";
    $NL .= sprintf("%d,%d,%d,",$_[0],$_[1],$_[2]);
    if (($_[3] == 0) && ($_[4] == -1)) {
        $NL .= ",,";
        $Fp = 0;
    }
    else {
        $NL .= sprintf("%10.4e,%10.4e,",$_[3],$_[4]);
    }
    if (($_[5] == 0) && ($_[6] == -1)) {
        $NL .= ",,";
        $Fm = 0;
    }
    else {
        $NL .= sprintf("%10.4e,%10.4e,",$_[5],$_[6]);
    }
    $NL .= sprintf("%10.4e,%10.4e,",$_[7],$_[8]);
    if (($Fp == 0) || ($Fm == 0)) {
        $NL .= ",,";
    }
    else {
        $NL .= sprintf("%10.4e,%10.4e,",$_[9],$_[10]);
    }

    $toto = (($Fp == 0)? (($Fm == 0 )? 0 : 2 ):(($Fm == 0 )? 1 : 0 ));
    $NL .= sprintf "%1d/",($toto) ;
    print OUT ("$NL\n");
#    print "$_\n M/ISYM = $toto, F+ = $Fp and F- = $Fm \n\n" if ($toto != 0);
}

close (IN);
close (OUT);
system "mv ${ARGV[0]}_conv2 ${ARGV[0]}_conv";

# Run f2mtz 

system << "EOT";

f2mtz hklin ${ARGV[0]}_conv hklout ${ARGV[1]} > f2mtz.log << EOF-f2mtz

TITLE transform $label (XDS) in CCP4 reflection file
PNAME GCM
DNAME $label
CELL $CELL
SYMMETRY $SYM

LABOUT H    K    L     F_${label}(+)  SIGF_${label}(+) F_${label}(-)  SIGF_${label}(-)   -
                     F_${label}     SIGF_${label}    DANO_${label}  SIGDANO_${label}   ISYM_${label}
CTYPE  H    H    H     G              L                G              L                  -
                     F              Q                D              Q                  Y
#FORMAT '(3(F5.0),8(1X,F11.9),2X,F1.0)'

EOF-f2mtz

EOT
;


