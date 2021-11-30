#!/bin/bash
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
#  Copyright (C) 2012  Universidade de São Paulo
#
#  Universidade de São Paulo
#  Laboratório de Biologia do Desenvolvimento de Abelhas
#  Núcleo de Bioinformática (LBDA-BioInfo)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://zulu.fmrp.usp.br/bioinfo 
#

#indir="./sim"
indir=$1
#outdir="./outsim"
outdir=$2
#refgtf="./ref/star_index/CanFam3.1.gtf"
refgtf=$3
#refstaridx="./ref/star_index"
refstaridx=$4

num_threads=20

#Diretório deste script
THIS_SCRIPT_HOME=$(dirname $0)


if [ ! ${indir} ]; then
	echo "Missing input directory"
	exit
else
	if [ ! -d ${indir} ]; then
		echo "Wrong input directory"
		exit
	fi
fi

if [ ! ${outdir} ]; then
	echo "Missing output directory"
	exit
else
	if [ ! -d ${outdir} ]; then
		echo "Wrong output directory"
		exit
	fi
fi

if [ ! ${refgtf} ]; then
	echo "Missing reference GTF file"
	exit
else
	if [ ! -e ${r2file} ]; then
		echo "Wrong reference GTF file (${refgtf})"
		exit
	else
		if [[ ! ${refgtf} =~ .gtf$ ]]; then
			echo "Wrong name of reference GTF file"
			exit
		fi
	fi
fi

if [ ! ${refstaridx} ]; then
	echo "Missing reference STAR index directory"
	exit
else
	if [ ! -d ${refstaridx} ]; then
		echo "Wrong reference STAR index directory"
		exit
	else
		if [ ! -e ${refstaridx}/SAindex ]; then
			echo "Cannot found STAR index files"
			exit
		fi
		if [ ! -e ${refstaridx}/genome.fa ]; then
			echo "Genome file named genome.fa must be inside STAR index directory (${refstaridx})"
			exit
		fi
	fi
fi

# Obtendo o nome das amostras baseado na convenção definida para o nome dos arquivos .fastq
# referentes a cada uma das amostras/réplicas biológicas
samps=()
samps=(`find ${indir} -name '*_R[12].fastq' -exec bash -c 'basename {} .fastq | sed "s/_R[12]$//"' \; | sort -u | xargs`)

if [ ${#samps[@]} == 0 ]; then
	echo "Warnings! Not found any ${indir}/*.fastq file to process!"

	samps=(`find ${indir} -name '*_R[12].fastq.gz' -exec bash -c 'basename {} .fastq.gz | sed "s/_R[12]$//"' \; | sort -u | xargs`)
	if [ ${#samps[@]} != 0 ]; then
		echo "There are ${#samps[@]} samples already processed!"
	fi
else
	for s in ${samps[@]}; do
		echo "Running preprocess.sh & alnasm.sh for ${s}"
		
		time ${THIS_SCRIPT_HOME}/alnasm.sh "${indir}/${s}_R1.fastq" "${indir}/${s}_R2.fastq" "${outdir}" "${refgtf}" "${refstaridx}"

	done
fi

mkdir -p ${outdir}/cufflinks

ready=(`find ${outdir}/cufflinks/* -type d -exec bash -c 'basename {}' 2> /dev/null \; | sort -u | xargs`)
if [ ${#ready[@]} != 0 ]; then
	echo "There are ${#ready[@]} samples assembled ready to proceed with gene expression analysis!"
else
	echo "Not found any assembled sample!"
	exit
fi

mkdir -p ${outdir}/cuffmerge
#mkdir -p ${outdir}/stringmerge

find ${outdir}/cufflinks/ -name transcripts.gtf > ${outdir}/cuffmerge/assembly_GTF_list.txt
#find ${outdir}/stringtie/ -name transcripts.gtf > ${outdir}/stringmerge/assembly_GTF_list.txt

if [ ! -s ${outdir}/cuffmerge/assembly_GTF_list.txt ]; then
#if [ ! -s ${outdir}/stringmerge/assembly_GTF_list.txt ]; then
	echo "Not found any transcriptome assembly"
	exit
fi

if [ ! -e ${outdir}/cuffmerge/merged.gtf ]; then
#if [ ! -e ${outdir}/stringmerge/merged.gtf ]; then
	cuffmerge	-o ${outdir}/cuffmerge \
			--ref-gtf ${refgtf} \
			--ref-sequence ${refstaridx}/genome.fa \
			--num-threads ${num_threads} \
			--min-isoform-fraction 0.20 \
			${outdir}/cuffmerge/assembly_GTF_list.txt

#	stringtie --merge 	-o ${outdir}/stringmerge/merged.gtf \
#				-G ${refgtf} \
#				-c 1 \
#				-F 1 \
#				-T 1 \
#				-g 50 \
#				-p ${num_threads} \
#				${outdir}/stringmerge/assembly_GTF_list.txt
fi

samps=(`find ${outdir}/star_out -name 'Aligned.out.sorted.bam' -exec bash -c 'dirname {}' \; | sort | xargs `)

declare -A GROUP

for s in ${samps[@]}; do
	bn=`basename ${s}`

	if [ ! -e ${outdir}/cuffquant/${bn}/abundances.cxb ]; then
		echo "Running cuffquant for ${bn}"
		cuffquant 	--output-dir ${outdir}/cuffquant/${bn} \
				--num-threads ${num_threads} \
				--frag-bias-correct ${refstaridx}/genome.fa \
				--multi-read-correct \
				--library-type fr-firststrand \
				--no-update-check \
				${outdir}/cuffmerge/merged.gtf \
				${outdir}/star_out/${bn}/Aligned.out.sorted.bam

#		cuffquant 	--output-dir ${outdir}/cuffquant/${bn} \
#				--num-threads ${num_threads} \
#				--frag-bias-correct ${refstaridx}/genome.fa \
#				--multi-read-correct \
#				--library-type fr-firststrand \
#				--no-update-check \
#				${outdir}/stringmerge/merged.gtf \
#				${outdir}/star_out/${bn}/Aligned.out.sorted.bam
	else
		echo "Using quantification for ${bn}"	
	fi
	
	grp=`echo ${bn} | sed 's/_B.*_T.*$//'`
	if [ "${GROUP[${grp}]}" ]; then
		GROUP[${grp}]="${GROUP[${grp}]} ${outdir}/cuffquant/${bn}/abundances.cxb"
	else
		GROUP[${grp}]="${outdir}/cuffquant/${bn}/abundances.cxb"
	fi
done

echo "Collecting samples for sample groups: $(IFS=, ; echo "${!GROUP[*]}") ..."
CXBS=""
for K in "${!GROUP[@]}"; do 
	GROUP[$K]=`echo "${GROUP[$K]}" | sed 's/ /,/g'`
	echo -e "\tSamples for $K: ${GROUP[$K]}"
	if [ "${CXBS}" ]; then
		CXBS="${CXBS} ${GROUP[$K]}"
	else
		CXBS="${GROUP[$K]}"
	fi
done

if [ ! -e ${outdir}/cuffdiff/gene_exp.diff ] || [ ! -s ${outdir}/cuffdiff/gene_exp.diff ]; then

	echo "Running cuffdiff among sample groups: $(IFS=, ; echo "${!GROUP[*]}") ..."

	cuffdiff	--output-dir ${outdir}/cuffdiff \
			--labels $(IFS=, ; echo "${!GROUP[*]}") \
			--num-threads ${num_threads} \
			--library-type fr-firststrand \
			--dispersion-method per-condition \
			--library-norm-method geometric \
			--total-hits-norm \
			--min-reps-for-js-test 2 \
			--no-update-check \
			${outdir}/cuffmerge/merged.gtf \
			${CXBS}

#	cuffdiff	--output-dir ${outdir}/cuffdiff \
#			--labels $(IFS=, ; echo "${!GROUP[*]}") \
#			--num-threads ${num_threads} \
#			--library-type fr-firststrand \
#			--dispersion-method per-condition \
#			--library-norm-method geometric \
#			--total-hits-norm \
#			--min-reps-for-js-test 2 \
#			--no-update-check \
#			${outdir}/stringmerge/merged.gtf \
#			${CXBS}
fi

if [ ! -e ${outdir}/cuffnorm/genes.count_table ] || [ ! -s ${outdir}/cuffnorm/genes.count_table ]; then
	echo "Running cuffnorm for sample groups: $(IFS=, ; echo "${!GROUP[*]}") ..."
	
	cuffnorm	--output-dir ${outdir}/cuffnorm \
			--labels $(IFS=, ; echo "${!GROUP[*]}") \
			--num-threads ${num_threads} \
			--library-type fr-firststrand \
			--library-norm-method geometric \
			--total-hits-norm \
			--no-update-check \
			${outdir}/cuffmerge/merged.gtf \
			${CXBS}

#	cuffnorm	--output-dir ${outdir}/cuffnorm \
#			--labels $(IFS=, ; echo "${!GROUP[*]}") \
#			--num-threads ${num_threads} \
#			--library-type fr-firststrand \
#			--library-norm-method geometric \
#			--total-hits-norm \
#			--no-update-check \
#			${outdir}/stringmerge/merged.gtf \
#			${CXBS}
fi

echo "... CONGRATULATIONS !!!"
