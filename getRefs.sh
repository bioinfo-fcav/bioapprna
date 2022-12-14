#!/bin/bash

mkdir -p ./ref

if [ ! -e "./ref/genome.fa" ]; then
	wget "https://i5k.nal.usda.gov/sites/default/files/data/Arthropoda/tricas-%28Tribolium_castaneum%29/GCF_000002335.3/1.Genome%20Assembly/Tcas5.2-GCF_000002335.3/Scaffolds/GCF_000002335.3_Tcas5.2_genomic_RefSeqIDs.fna.gz" -O ./ref/GCF_000002335.3_Tcas5.2_genomic_RefSeqIDs.fna.gz
	gunzip ./ref/GCF_000002335.3_Tcas5.2_genomic_RefSeqIDs.fna.gz
	sed 's/^>\(\S\+\).*/>\1/' ./ref/GCF_000002335.3_Tcas5.2_genomic_RefSeqIDs.fna > ./ref/genome.fa
fi
exit
if [ ! -e "./ref/genome.gtf" ]; then
	wget "https://i5k.nal.usda.gov/sites/default/files/data/Arthropoda/tricas-%28Tribolium_castaneum%29/GCF_000002335.3/3.Additional%20Gene%20Sets%20and%20Annotation%20Projects/NCBI/NCBI_Tribolium_castaneum_Annotation_Release_103/GCF_000002335.3_Tcas5.2_genomic.gff.gz" -O ./ref/GCF_000002335.3_Tcas5.2_genomic.gff.gz
	gunzip ./ref/GCF_000002335.3_Tcas5.2_genomic.gff.gz
	fixNCBIgff.sh ./ref/GCF_000002335.3_Tcas5.2_genomic.gff ./ref/genome.gff
	gffread ./ref/genome.gff -g ./ref/genome.fa -T -o ./ref/original_genome.gtf
	gffread ./ref/genome.gff -g ./ref/genome.fa -w ./ref/transcritome.fa
	./replace_gene_id.pl ./ref/genome.gff ./ref/original_genome.gtf > ./ref/genome.gtf
fi

if [ ! -e "./ref/gene_info" ]; then
	wget "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz" -O ./ref/gene_info.gz
	gunzip ./ref/gene_info.gz
fi

head -1  ./ref/gene_info > ./ref/gene_info_Tcas.txt
grep -P "^7070\t" ./ref/gene_info >> ./ref/gene_info_Tcas.txt

