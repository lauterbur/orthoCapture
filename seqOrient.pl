#!/usr/local/bin/perl -w

# Copyright 2013, Naoki Takebayashi <ntakebayashi@alaska.edu>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# Version: 20130612

# POSIX::tmpnam() deprecated, replace with File::Temp - MEL April 2017

my $usage =
    "Usage: $0 [-h] [-r refSeqNumber] inputfile\n".
    "  -h: help\n" .
    "  -r integer: Specified sequence is used as the reference\n" .
    "\n" .
    "This program read in the sequence file, which may contain sequences " .
    "with opposite orientation.  By default, it will use the 1st sequence as the " .
    "reference. It makes the complement of the sequences with revseq of " .
    "EMBOSS and see if the complement aligns better with the reference seq. ".
    "If so, the complement will be used.  It will print out the fasta file " .
    "with the corrected orientation to the STDOUT.  For the pairwise ".
    "alignment, matcher of EMBOSS is used.  The scores of the alignments ".
    "are printed to STDERR.\n";

use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Factory::EMBOSS;
use IO::File;
#use POSIX qw(tmpnam);
use File::Temp qw/ tmpnam /;
use Getopt::Std;

getopts('r:h') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

my $infile = shift;

my $outFH = Bio::SeqIO->newFh(-format => 'fasta');

# read in the data
$in = Bio::SeqIO->new(-file => $infile, '-format' => 'fasta');

my @seqArray =();
while(my $seq = $in->next_seq()) {
    push @seqArray, $seq;
}


# compare
my @result = ();

my $refIndex = 1;
if (defined($opt_r)) {
    if ($opt_r <= scalar(@seqArray)) {
	$refIndex = $opt_r;
    } else {
	warn "WARN: $opt_r is greater than the number of sequences=", 
	scalar(@seqArray), 
	".  Using the sequence $refIndex for the reference.\n";
    }
}

my $ref = splice (@seqArray, $refIndex-1, 1);
push @result, $ref;
foreach my $i (@seqArray) {
    my ($ret1, $ret2) = compOrientation($ref, $i);
    push @result, $ret2;
}

foreach my $i (@result) {
    print $outFH $i;
}

exit(0);

sub compOrientation {
    my ($seq1, $seq2) = @_;

    my $id1 = $seq1->display_id();
    my $id2 = $seq2->display_id();

    # align with matcher
    my $align1 = embossMatcher($seq1, $seq2);

    # check the complement align
    my $comp = Complement($seq2);
    my $align2 = embossMatcher($seq1, $comp);

#    print $align1->length, "\n";
#    print $align1->no_residues, "\n";
#    print $align1->score, "\n";
#    print $align1->percentage_identity, "\n";

    my $score1 = $align1->score;
    my $score2 = $align2->score;
    if ($score1 < $score2) {
	warn "$id1 - $id2: score reg=$score1, comp=$score2, " .
	    "complement of $id2 is used\n";
	return ($seq1, $comp);
    } else {
	warn "$id1 - $id2: score reg=$score1, comp=$score2\n";
	return ($seq1, $seq2)
    }
}

sub embossMatcher {
    my ($seq1, $seq2) = @_;
    my $factory = new Bio::Factory::EMBOSS;
    my $prog = $factory->program('matcher');
    
    my $tempOutfile;
    my $fh;
    do {$tempOutfile = tmpnam()} 
    until $fh = IO::File->new($tempOutfile, O_RDWR|O_CREAT|O_EXCL);
    $fh->close;

    $prog->run({ -asequence => $seq1,
		 -bsequence => $seq2,
		 -aformat      => "pair",
		 -alternatives => 1,
		 -outfile     => $tempOutfile});

    my $alignio_fmt = "emboss";
    my $align_io = new Bio::AlignIO(-format => $alignio_fmt,
				    -file   => $tempOutfile);
#    $out = Bio::AlignIO->new(-format => 'pfam');
#    while ( my $aln = $align_io->next_aln() ) { $out->write_aln($aln); };

    unlink $tempOutfile || die "ERROR: Unable to unlink $tempOutfile\n";
    return($align_io->next_aln());
}

sub Complement {
    my $seq = shift;

    my $tempOutfile;
    my $fh;
    do {$tempOutfile = tmpnam()} 
    until $fh = IO::File->new($tempOutfile, O_RDWR|O_CREAT|O_EXCL);
    $fh->close;

    my $factory = new Bio::Factory::EMBOSS;
    my $revseq = $factory->program('revseq');

    my %input = (-sequence => $seq,
		 -outseq => $tempOutfile );

    $revseq->run(\%input);
    my $seqio = Bio::SeqIO->new (-file => $tempOutfile);

    unlink $tempOutfile || die "ERROR: Unable to unlink $tempOutfile\n";

    return($seqio->next_seq())
}
