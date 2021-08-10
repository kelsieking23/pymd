#!/usr/bin/perl

#
# plot_ss_xpm.pl - a program that reads in an ss.xpm file to generate individual plots
# of secondary structure propensity per residue
#   * Probabilities instead of raw counts
#

use strict;

# input defined as first command-line argument
my $input = $ARGV[0];

# open the input
open(IN, "<$input");
my @in = <IN>;
close(IN);

# determine number of residues and frames
my $nres = 0;
my $nframes = 0;
for (my $i=0; $i<scalar(@in); $i++) {
    if ($in[$i] =~ /static char/) {
        my $res_line = $in[$i+1];
        my @info = split(" ", $res_line);

        $nframes = $info[0];
        my @nframes = split('', $nframes);
        shift(@nframes);    # get rid of the "
        $nframes = join('', @nframes);

        $nres = $info[1];
    }
}

print "I found $nres residues.\n";
print "There are $nframes frames.\n";

# initialize hashes for later output writing
my %a_helix;
for (my $a=1; $a<=$nres; $a++) {
    $a_helix{$a} = 0;
}

my %three_helix;
for (my $th=1; $th<=$nres; $th++) {
    $three_helix{$th} = 0;
}

my %pi_helix;
for (my $p=1; $p<=$nres; $p++) {
    $pi_helix{$p} = 0;
}

my %total_helix;
for (my $h=1; $h<=$nres; $h++) {
    $total_helix{$h} = 0;
}

my %sheet;
for (my $s=1; $s<=$nres; $s++) {
    $sheet{$s} = 0;
}

my %coil;
for (my $c=1; $c<=$nres; $c++) {
    $coil{$c} = 0;
}

my %turn;
for (my $t=1; $t<=$nres; $t++) {
    $turn{$t} = 0;
}

my %bridge;
for (my $b=1; $b<=$nres; $b++) {
    $bridge{$b} = 0;
}

my %bend;
for (my $b=1; $b<=$nres; $b++) {
    $bend{$b} = 0;
}

# The last nres lines have the secondary structure information
# So, if there are 40 residues that were analyzed, there are 40 lines
# of continuous secondary structure codes

# clean up the output - up to 18 lines of comments, etc.
splice(@in, 0, 18);

# remove any "x-axis" or "y-axis" lines
for (my $n=0; $n<scalar(@in); $n++) {
    if (($in[$n] =~ /x-axis/) || ($in[$n] =~ /y-axis/)) {
            shift(@in);
            $n--;
    }
}

# DEBUG
# open(OUT_DUMMY, ">>test");
# printf OUT_DUMMY "@in\n";
# close(OUT_DUMMY);

# There should now be $nres lines left in the file
# The SS codes for the last residue are written first (top-down in .xpm file)
#   * Element 0 is residue $nres, element 1 is $nres-1, etc.

for (my $i=$nres; $i>=1; $i--) {
    # There will be $nframes+2 elements in @line (extra two are " at beginning
    # and end of the line)
    # Establish a conversion factor and split the input lines
    my $j = $nres - $i;
    my @line = split('', $in[$j]);

    # for each type of ss, write to hashes

    for (my $k=1; $k<=($nframes+1); $k++) {
        if ($line[$k] =~ /H/) {
            $a_helix{$i}++;
            $total_helix{$i}++;
        } elsif ($line[$k] =~ /G/) {
            $three_helix{$i}++;
            $total_helix{$i}++;
        } elsif ($line[$k] =~ /I/) {
            $pi_helix{$i}++;
            $total_helix{$i}++;
        } elsif ($line[$k] =~ /E/) {
            $sheet{$i}++;
        } elsif ($line[$k] =~ /~/) {
            $coil{$i}++;
        } elsif ($line[$k] =~ /T/) {
            $turn{$i}++;
        } elsif ($line[$k] =~ /B/) {
            $bridge{$i}++;
        } elsif ($line[$k] =~ /S/) {
            $bend{$i}++;
        }
    }
}

# open a single output file for writing
# each secondary structure type is a separate data series

open(OUT, ">>summary_SS.xvg") or die "Cannot open \"summary_SS.dat\"\n";
printf OUT "# Probability of various secondary structure elements, by residue.\n";
printf OUT "@\ttitle\t\"Secondary Structure Content\"\n";
printf OUT "@\txaxis\tlabel \"Residue\"\n";
printf OUT "@\tyaxis\tlabel \"Probability\"\n";
printf OUT "\@TYPE xy\n";
printf OUT "@ s0 legend \"\\f{Symbol}a\\f{Times}-Helix\"\n";
printf OUT "@ s1 legend \"\\f{Symbol}p\\f{Times}-Helix\"\n";
printf OUT "@ s2 legend \"3\\s10\\N-Helix\"\n";
printf OUT "@ s3 legend \"\\f{Symbol}b\\f{Times}-Strand\"\n";
printf OUT "@ s4 legend \"\\f{Symbol}b\\f{Times}-Bend\"\n";
printf OUT "@ s5 legend \"\\f{Symbol}b\\f{Times}-Turn\"\n";
printf OUT "@ s6 legend \"\\f{Symbol}b\\f{Times}-Bridge\"\n";
printf OUT "@ s7 legend \"Random Coil\"\n";

for (my $o=1; $o<=$nres; $o++) {
    printf OUT "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", $o, ($a_helix{$o}/$nframes), ($pi_helix{$o}/$nframes), ($three_helix{$o}/$nframes), ($sheet{$o}/$nframes), ($bend{$o}/$nframes), ($turn{$o}/$nframes), ($bridge{$o}/$nframes), ($coil{$o}/$nframes);
}

close(OUT);

exit;
