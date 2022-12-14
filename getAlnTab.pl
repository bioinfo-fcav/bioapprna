#!/usr/bin/perl
#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2018  Universidade Estadual Paulista - UNESP
#
#  Universidade Estadual Paulista - UNESP
#  Laboratório de Bioinformática
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#
# $Id$

=head1 NAME

=head1 SYNOPSIS

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Arguments:

        -h/--help   Help
        -l/--level  Log level [Default: FATAL] 
            OFF
            FATAL
            ERROR
            WARN
            INFO
            DEBUG
            TRACE
            ALL

=head1 AUTHOR

Daniel Guariz Pinheiro E<lt>dgpinheiro@gmail.comE<gt>

Copyright (c) 2018 Universidade Estadual Paulista - UNESP

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Readonly;
use Getopt::Long;

use File::Basename;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $indir, $aligner, $stats_file);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help"      => sub { &Usage(); },
            "l|level=s"     => \$level,
            "i|indir=s"     => \$indir,
            "a|aligner=s"   => \$aligner,
            "s|stats=s"     => \$stats_file
    ) or &Usage();


if ($level) {
    my %LEVEL = (   
    'OFF'   =>$OFF,
    'FATAL' =>$FATAL,
    'ERROR' =>$ERROR,
    'WARN'  =>$WARN,
    'INFO'  =>$INFO,
    'DEBUG' =>$DEBUG,
    'TRACE' =>$TRACE,
    'ALL'   =>$ALL);
    $LOGGER->logdie("Wrong log level ($level). Choose one of: ".join(', ', keys %LEVEL)) unless (exists $LEVEL{$level});
    Log::Log4perl->easy_init($LEVEL{$level});
}

$LOGGER->logdie("Missing input directory") unless ($indir);
$LOGGER->logdie("Wrong input directory ($indir)") unless (-d $indir);

my %valid=( 'tophat'    =>  1,
            'star'      =>  1);

$LOGGER->logdie("Missing aligner") unless ($aligner);
$LOGGER->logdie("Wrong aligner ($aligner)") unless (exists $valid{$aligner});

$LOGGER->logdie("Missing previous stats file") unless ($stats_file);
$LOGGER->logdie("Wrong previous stats file ($stats_file)") unless (-e $stats_file);

my %count;

my @sample;
foreach my $if (sort { $a cmp $b } glob("$indir/*.$aligner.out.txt")) {
    my ($name) = basename($if, '.'.$aligner.'.out.txt');
    #print $name,"\t",$if,"\n";
    push(@sample, $name);

    open(IN, "<", $if) or $LOGGER->logdie($!);
    while(<IN>) {
        chomp;
        if ($_=~/^\s+(\d+) were paired;/) {
            $count{$name}->{'pe'}->{'total'} = $1;
        } elsif ($_=~/^\s+(\d+) aligned concordantly 0 times/) {
            $count{$name}->{'pe'}->{'nonconcordant'} = $1;
        } elsif ($_=~/^\s+(\d+) aligned concordantly exactly 1 time/) {
            $count{$name}->{'pe'}->{'concordant_unique'} = $1;
        } elsif ($_=~/^\s+(\d+) aligned concordantly >1 times/) {
            $count{$name}->{'pe'}->{'concordant_multiple'} = $1;
        } elsif ($_=~/^\s+(\d+) were single-end reads;/) {
            $count{$name}->{'se'}->{'total'} = $1;
        } elsif ($_=~/^\s+(\d+) single-end reads mapped uniquely/) {
            $count{$name}->{'se'}->{'unique'} = $1;
        } elsif ($_=~/^\s+(\d+) single-end reads mapped >1 times/) {
            $count{$name}->{'se'}->{'multiple'} = $1;
        }
    }
    close(IN);

}

open(STATS, "<", $stats_file) or $LOGGER->logdie($!);
my $stats_header_line=<STATS>;
chomp($stats_header_line);
my @stats_header=split(/\t/, $stats_header_line);

while(<STATS>) {
    chomp;
    my %data;
    @data{@stats_header} = split(/\t/, $_);
    $count{ $data{'Sample'} }->{ 'pe' }->{ 'input' } = sprintf("%.0f", $data{'Free of contamination PE R1'});
    $count{ $data{'Sample'} }->{ 'se' }->{ 'input' } = sprintf("%.0f", $data{'Free of contamination SE R1'} + $data{'Free of contamination SE R2'});
}

close(STATS);

my @header = ('pe.input', 'pe.total', 'pe.nonconcordant', 'pe.concordant_unique', 'pe.concordant_multiple', 'se.input', 'se.total', 'se.unique', 'se.multiple');

print join("\t", 'Sample', @header),"\n";

foreach my $name (@sample) {
    my @value;
    foreach my $h (@header) {
        my ($layout, $type) = split(/\./, $h);
        push(@value, $count{$name}->{$layout}->{$type}||0);
    }
    print join("\t", $name,@value),"\n";
}

# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -i      --indir     Input directory with the output files from SAM_nameSorted_to_uniq_count_stats.pl ("./output/align_out_info")
        -a      --aligner   Aligner (tophat or star)
        -s      --statsfile Previous stats file from "stats_preprocess.sh"

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

