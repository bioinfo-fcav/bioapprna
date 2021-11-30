#!/bin/bash

# Recebendo valores das variáveis a partir dos argumentos 
# da linha de comando
r1file=$1
r2file=$2
outdir=$3

# Ex.: refgtf="./ref/star_index/CanFam3.1.gtf"
refgtf=$4

# Ex.: refstaridx="./ref/star_index"
refstaridx=$5

# Número de threads (subprocessos) que serão lançados em paralelo
num_threads=20

THIS_SCRIPT_HOME=$(dirname $0)

# Execução do pré-processamento
${THIS_SCRIPT_HOME}/preprocess.sh "${r1file}" "${r2file}" "${outdir}" "${refgtf}"

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

# Capturando os nomes dos arquivos, ou seja, sem os diretórios, e, também, 
# sem a extensão .fastq
bn1=`basename ${r1file} .fastq`
bn2=`basename ${r2file} .fastq`
# Substituindo no nome do arquivo R1 o trecho "_R1" por nada para ser
# o nome base independente de ser R1 ou R2, por exemplo, no caso dos arquivos
# de LOG para atropos e prinseq
bn=`echo ${bn1} | sed 's/_R1//'`

outdir_star_pe="${outdir}/star_out_pe/${bn}"
outdir_star_se="${outdir}/star_out_se/${bn}"
outdir_star="${outdir}/star_out/${bn}"
mkdir -p ${outdir_star_pe}
mkdir -p ${outdir_star_se}
mkdir -p ${outdir_star}

if [ ! -e "${outdir_star}/Aligned.out.sorted.bam" ]; then
	
	echo "Running STAR alignment ..."
	
	if [ ! -e "${outdir_star_pe}/Aligned.out.bam" ]; then
		echo "	Running STAR alignment [PE] ..."

		# Alinhamento com STAR (paired-end)

		STAR  	--runMode alignReads \
			--limitBAMsortRAM 15000000000 \
			--runThreadN ${num_threads} \
			--alignEndsType EndToEnd \
			--sjdbInsertSave All \
			--outSAMattributes All \
			--outSAMprimaryFlag AllBestScore \
			--outSAMtype BAM Unsorted \
			--outSAMorder paired \
			--sjdbOverhang 149\
			--sjdbGTFfile ${refgtf} \
			--genomeDir ${refstaridx} \
			--outFileNamePrefix ${outdir_star_pe}/ \
			--readFilesIn  ${outdir}/processed/prinseq/${bn1}.atropos.prinseq.fastq \
				       ${outdir}/processed/prinseq/${bn2}.atropos.prinseq.fastq \
			--sjdbScore 8 \
			--outSAMunmapped Within \
			--chimOutType WithinBAM \
			-outFilterMultimapNmax 20 \
			--alignSJoverhangMin 8 \
			--alignSJDBoverhangMin 1 \
			--outFilterMismatchNoverReadLmax 0.05 \
			--alignIntronMin 40 \
			--alignIntronMax 72170 \
			--alignMatesGapMax 1000000 \
			--outFilterScoreMinOverLread 0 \
			--outFilterMatchNminOverLread 0 \
			--outFilterMatchNmin 0  \
			--scoreGapNoncan -4 \
			--seedSearchStartLmax 10 \
			--outFilterMismatchNmax 999

		# conversão da saída SJ.out.tab gerada pelo STAR para o formato BED reconhecido pelo STAR para
		# o parâmetro "--sjdbFileChrStartEnd"

		perl -F"\t" -lane ' if (  ($F[4]>0) && ($F[5]>=0) && ($F[6]>=5) && ($F[7]>=0) && ($F[8]>10)  ) { print join("\t", @F[0,1,2],(($F[3]==1) ? "+" : (($F[3]==2) ? "-" : ".")) ); }' ${outdir_star_pe}/SJ.out.tab > ${outdir_star_pe}/junctions.junc

	# Alinhamento com STAR (single-end)

	fi

	if [ ! -e "${outdir_star_se}/Aligned.out.bam" ]; then
		
		echo "Running STAR alignment [SE] ..."
		
		STAR    --runMode alignReads \
			--limitBAMsortRAM 15000000000 \
			--runThreadN ${num_threads} \
			--alignEndsType EndToEnd \
			--sjdbInsertSave All \
			--outSAMattributes All \
			--outSAMprimaryFlag AllBestScore \
			--outSAMtype BAM Unsorted \
			--outSAMorder paired \
			--sjdbOverhang 149 \
			--sjdbGTFfile ${refgtf} \
			--genomeDir ${refstaridx} \
			--outFileNamePrefix ${outdir_star_se}/ \
			--readFilesIn ${outdir}/processed/prinseq/${bn1}_singletons.atropos.prinseq.fastq,${outdir}/processed/prinseq/${bn2}_singletons.atropos.prinseq.fastq \
			--sjdbScore 8 \
			--outSAMunmapped Within \
			--chimOutType WithinBAM \
			-outFilterMultimapNmax 20 \
			--alignSJoverhangMin 8 \
			--alignSJDBoverhangMin 1 \
			--outFilterMismatchNoverReadLmax 0.05 \
			--alignIntronMin 40 \
			--alignIntronMax 72170 \
			--alignMatesGapMax 1000000 \
			--outFilterScoreMinOverLread 0 \
			--outFilterMatchNminOverLread 0 \
			--outFilterMatchNmin 0  \
			--scoreGapNoncan -4 \
			--seedSearchStartLmax 10 \
			--outFilterMismatchNmax 999
	fi

	# Combinar resultados do alinhamento com reads paired-end e alinhamento com reads single-end (singletons)	
	samtools merge -@ ${num_threads} -f -n 	${outdir_star}/Aligned.out.bam \
						${outdir_star_pe}/Aligned.out.bam \
						${outdir_star_se}/Aligned.out.bam

	# Ordenando o resultado do alinhamento por coordenadas genômicas
	# - exigência para executar o cufflinks
	samtools sort -@ ${num_threads} -o	${outdir_star}/Aligned.out.sorted.bam \
						${outdir_star}/Aligned.out.bam

	rm -f ${outdir_star}/Aligned.out.bam
fi

###
# Montagem dos transcritomas
###

mkdir -p ${outdir}/cufflinks/${bn}

if [ ! -e "${outdir}/cufflinks/${bn}/transcripts.gtf" ]; then

	echo "Running Assembly [cufflinks] ..."

	cufflinks 	--output-dir ${outdir}/cufflinks/${bn} \
			--num-threads ${num_threads} \
			--GTF ${refgtf}\
			--frag-bias-correct ${refstaridx}/genome.fa \
			--multi-read-correct \
			--library-type fr-firststrand \
			--frag-len-mean 300 \
			--frag-len-std-dev 50 \
			--total-hits-norm \
			--min-isoform-fraction 0.20 \
			--small-anchor-fraction 0.10 \
			--min-frags-per-transfrag 10 \
			--max-intron-length 72170 \
			--min-intron-length 40 \
			--no-update-check \
			${outdir_star}/Aligned.out.sorted.bam
		
	#mkdir -p ${outdir}/stringtie/${bn}
	#
	#stringtie 	${outdir_star}/Aligned.out.sorted.bam \
	#		-o ${outdir}/stringtie/${bn}/transcripts.gtf \
	#		-p ${num_threads} \
	#		-G ${refgtf} \
	#		--rf \
	#		-B \
	#		-t \
	#		-c 2 \
	#		-C ${outdir}/stringtie/${bn}/cov_refs.gtf \
	#		-A ${outdir}/stringtie/${bn}/gene_abund.tab
fi

