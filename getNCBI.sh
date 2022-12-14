#!/bin/bash

#HASH / ASSOCIATIVE ARRAY / DICTIONARY / KEY=VALUE
declare -A CORRIDA=(	["SRR15082126"]="RNAI_B1"
			["SRR15082127"]="RNAI_B2"
			["SRR15082128"]="RNAI_B3"
			["SRR15082129"]="CTRL_B1"
			["SRR15082130"]="CTRL_B2"
			["SRR15082131"]="CTRL_B3"
) 

mkdir -p ./raw
mkdir -p ./input

for srr in ${!CORRIDA[@]}; do 
	echo "Baixando: ${srr}"
	wget "https://sra-pub-run-odp.s3.amazonaws.com/sra/${srr}/${srr}" -O ./raw/${srr}.sra 2>/dev/null
	echo "Extraindo fastq: ${srr}"
	fastq-dump --split-3 ./raw/${srr}.sra --outdir ./raw
	
	echo "Removendo sra: ${srr}"
	rm -f ./raw/${srr}.sra

	echo "Obtendo uma amostra (primeiras leituras) de 1.000.000 de leituras: ${srr}"
	head -4000000 ./raw/${srr}_1.fastq > ./raw/${srr}_1mi_1.fastq
	head -4000000 ./raw/${srr}_2.fastq > ./raw/${srr}_1mi_2.fastq

	echo "Criando atalhos no diretório de entrada: ${srr} ${CORRIDA[${srr}]}"
	ln -s ../raw/${srr}_1mi_1.fastq ./input/${CORRIDA[${srr}]}_R1.fastq
	ln -s ../raw/${srr}_1mi_2.fastq ./input/${CORRIDA[${srr}]}_R2.fastq
done

