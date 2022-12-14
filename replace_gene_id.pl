#!/bin/env perl

use strict;
use warnings;


my $gff = $ARGV[0];

die "Missing GFF" unless ($gff);
die "Wrong GFF ($gff)" unless (-e $gff);

my $gtf = $ARGV[1];

die "Missing GTF" unless ($gtf);
die "Wrong GTF ($gtf)" unless (-e $gtf);

my %gene_id;
my %rna_id;
open(GFF, "<", $gff) or die $!;
while(<GFF>) {
	chomp;
	next if ($_=~/^#/);
	my @F=split(/\t/, $_);
	#NZ_CP058354.1   RefSeq  gene    5624    7792    .       +       .       ID=gene-HXT67_RS00030_1;Dbxref=GeneID:64065474;Name=HXT67_RS00030;gbkey=Gene;gene_biotype=protein_coding;locus_tag=HXT67_RS00030
	
	if ($F[2] eq 'gene') {
		if ($F[8]=~/ID=([^;]+)/) {
			my $ID=$1;
			if ($F[8]=~/GeneID:(\d+)/) {
				my $GeneID=$1;
				$gene_id{$ID}=$GeneID;
			} else {
				$gene_id{$ID}=$ID;
			}
		} else {
			die "Not found an ID for the gene ($F[8])";
		}
	} elsif ($F[2] eq 'RNA') {
		if ($F[8]=~/ID=([^;]+)/) {
			my $ID=$1;
			if ($F[8]=~/transcript_id=([^;]+)/) {
				my $transcript_id=$1;
				$rna_id{$ID}=$transcript_id;
			} else {
				$rna_id{$ID}=$ID;
			}
		} else {
			die "Not found an ID for the RNA ($F[8])";
		}
	}
}
close(GFF);

#foreach my $ID (keys %gene_id) {
#	print $ID,"\t",$gene_id{$ID},"\n";
#}

open(GTF, "<", $gtf) or die $!;
while(<GTF>) {
	chomp;
	if ($_=~/gene_id "([^"]+)"/) {
		my $ID=$1;
		my $newID=$gene_id{$ID};
		$_=~s/gene_id "[^"]+"/gene_id "$newID"/;
		print $_,"\n";
	} else {
		die "Not foun gene_id attribute in GTF line ($_)";
	}
	if ($_=~/transcript_id "([^"]+)"/) {
		my $ID=$1;
		my $newID=$rna_id{$ID};
		$_=~s/transcript_id "[^"]+"/transcript_id "$newID"/;
		print $_,"\n";
	} else {
		die "Not foun transcript_id attribute in GTF line ($_)";
	}
}
close(GTF);

