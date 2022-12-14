#!/bin/bash

# definição do número de threads
num_threads=$1

# input - diretório contendo os arquivos de entrada no formato .fastq
input=$2
#Diretorio que contem os arquivos processados que servirão para alinhamento com o bowtie2

refassembly=$3
# Diretorio que contem a montagem de novo Trinity.fa

reftrans=$4
# Arquivo que contem o transcriptoma de referencia mais proximo para ser utilizado como base, para se fazer a distribuição dos transcritos

if [ ! ${num_threads} ]
then   
        echo "Missing number of threads"
        exit
else   
        if [[ ! ${num_threads} =~ ^[0-9]+$ ]]
        then   
                echo "Wrong number of threads ${num_threads}"
                exit
        fi
fi

if [ ! ${input} ]
then   
	echo "Missing input directory (<input-dir>/*.fastq files)"
        exit
else   
        if [ ! -d ${input} ]
        then   
		echo "Wrong input directory (<input-dir>/*.fastq files): ${input}"
                exit
        fi
fi

if [ ! ${refassembly} ]
then
        echo "Missing Trinity / PASA output directory"
        exit
else
        if [ ! -d ${refassembly} ]
        then
                echo "Wrong Trinity / PASA output directory ${refassembly}"
                exit
        fi
fi

if [ ! ${reftrans} ]
then
        echo "Missing reference transcriptome fasta file"
        exit
else
        if [ ! -d ${refs} ]
        then
                echo "Wrong reference transcriptome fasta file ${reftrans}"
                exit
        fi
fi

# output - diretório para armazenar o resultado do processo de montagem
output=$5


if [ ! ${output} ]
then   
        echo "Missing output directory"
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "Wrong output directory ${output}"
                exit
        fi
fi

readlen=$6

if [ ! ${readlen} ]
then
	readlen=150
	echo "Warning: Read length not defined. Using default (${readlen})"
fi

prefix=$7

if [ ! ${prefix} ]
then
	prefix='Trinity'
fi


reffasta=`find ${refassembly} -type f -maxdepth 1 -name "*${prefix}*.fasta" 2>/dev/null`
refgene_trans_map=`find ${refassembly} -type f -maxdepth 1 -name "*${prefix}*.gene_trans_map" 2>/dev/null`

if [ ! ${reffasta} ]; then
	echo "Error: Could'nt find reference assembly fasta file using prefix (${prefix})"
	exit
fi

if [ ! ${refgene_trans_map} ]; then
	echo "Error: Could'nt find reference assembly gene_trans_map file using prefix (${prefix})"
	exit
fi

### Arquivos e diretórios de saída (output) 

basedir_out="${output}/Detonate"
mkdir -p ${basedir_out}

bowtie2_index="${basedir_out}/bowtie2_index"
mkdir -p ${bowtie2_index}

bowtie2_index="${bowtie2_index}/transcriptome"

bowtie2="${basedir_out}/bowtie2"
mkdir -p ${bowtie2}

if [ ! -e "${bowtie2_index}.1.bt2" ]; then 
	bowtie2-build 	-f ${reffasta} ${bowtie2_index} \
	 		 > ${bowtie2_index}.bowtie2-build.log.out.txt \
			2> ${bowtie2_index}.bowtie2-build.log.err.txt
fi




if [ ! -e "${basedir_out}/All.aligned.sorted.bam" ]; then

	for r1 in `ls ${input}/*_1.fastq`; do

		r2=`echo ${r1} | sed 's/_1.fastq/_2.fastq/'`

	        if [ ! -e "${r2}" ]; then
			echo "Read 2 (${r2}) paired with Read 1 ($r1) not found."
		       exit
		fi

		name=`basename ${r1} | sed 's/\..*//' `

	        echo -e "Read alignment using bowtie2 ${name}: ${r1} & ${r2} ...\n"

		bowtie2 -p ${num_threads} \
			--dpad 0 \
			--gbar 99999999 \
			--mp 1,1 \
			--np 1 \
			--score-min L,0,-0.1 \
			--no-mixed \
			--no-discordant \
			-t -I 1 -X 1000 \
			--end-to-end \
			--fr \
			--very-sensitive \
			-x ${bowtie2_index} \
			-1 ${r1} \
			-2 ${r2}  \
			 2>  ${bowtie2}/${name}.log.txt \
			| samtools view -S -b -o ${bowtie2}/${name}.bam -


		if [ ! -e "${basedir_out}/bowtie2_header.sam" ]; then	
			echo "Samtools View header"
			samtools view -H ${bowtie2}/${name}.bam > ${basedir_out}/bowtie2_header.sam
		fi
       
	done

	if [ ! -e "${basedir_out}/All.aligned.bam" ]; then
		bamfiles=()

		bamfiles=( $( find ${bowtie2} -name '*.bam' ) )

		echo " Merging bamfiles"
		samtools merge -f ${basedir_out}/All.aligned.bam ${bamfiles[*]}
		samtools merge --threads ${num_threads} -f -h ${basedir_out}/bowtie2_header.sam ${basedir_out}/All.aligned.bam ${bamfiles[*]}
	fi
	
	echo " Sorting ${basedir_out}/All.aligned.bam "
	
	samtools sort -n --threads ${num_threads} -o ${basedir_out}/All.aligned.sorted.bam ${basedir_out}/All.aligned.bam

	if [ -e "${basedir_out}/All.aligned.sorted.bam" ]; then
		echo "   Removing All.aligned.bam ..."
		rm -f ${basedir_out}/All.aligned.bam
		echo "   Removing ${bowtie2} ..."
		rm -fr ${bowtie2}
	fi
	
	if [ ! -e "${basedir_out}/All.selected.aligned.sorted.bam" ]; then
		echo "	Exclude alignments with reads with less than 25 bp ..."
		samtools view ${basedir_out}/All.aligned.sorted.bam | perl -F"\t" -lane ' if ( length($F[9]) < 25 ) { print $F[0]; }' | nsort -u > ${basedir_out}/shorts25.txt

		samtools view -h ${basedir_out}/All.aligned.sorted.bam | grep -v -wf ${basedir_out}/shorts25.txt  > ${basedir_out}/All.selected.aligned.sorted.sam

		samtools view -bS ${basedir_out}/All.selected.aligned.sorted.sam > ${basedir_out}/All.selected.aligned.sorted.bam

		if [ -e "${basedir_out}/All.selected.aligned.sorted.bam" ]; then
			echo "	Removing All.aligned.sorted.bam ..."
			rm -f ${basedir_out}/All.aligned.sorted.bam
		fi

	fi		
fi


if [ ! -e "${basedir_out}/transcript-len-distrib.txt" ]; then

	echo "Evaluating the reference transcript length distribution ..."

	rsem-eval-estimate-transcript-length-distribution ${reftrans} ${basedir_out}/transcript-len-distrib.txt
fi

if [ ! -e "${basedir_out}/PE-detonate.score" ]; then


	echo "Calculating score ..."

	rsem-eval-calculate-score 	--forward-prob 0.5 \
					--transcript-to-gene-map ${refgene_trans_map} \
					--tag "" \
        	               	        --num-threads ${num_threads} \
					--transcript-length-parameters ${basedir_out}/transcript-len-distrib.txt \
					--paired-end \
               		                --bam ${basedir_out}/All.selected.aligned.sorted.bam \
					${reffasta} ${basedir_out}/PE-detonate ${readlen}
fi
