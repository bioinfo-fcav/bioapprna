#!/bin/bash

# ALINHADOR tophat OU star
aligner=${1}

if [ ! ${aligner} ]; then
	echo "[ERROR] Missing aligner (tophat or star)." 1>&2
	exit
fi

if [ "${aligner}" != "tophat" ] &&
   [ "${aligner}" != "star" ]; then

	echo "[ERROR] Aligner must be \"tophat\" or \"star\" (${aligner})." 1>&2
  	exit
fi

indir=${2}

# SE ${indir} NÃO EXISTE, OU SEJA, SE NÃO FOI PASSADO ARGUMENTO 1 NA LINHA DE COMANDO
if [ ! ${indir} ]; then
	echo "[ERROR] Missing input directory." 1>&2
	exit
fi

# SE ${indir} NÃO É DIRETÓRIO
if [ ! -d ${indir} ]; then
	echo "[ERROR] Wrong input directory (${indir})." 1>&2
	exit
fi

outdir=${3}

# SE ${outdir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 2 NA LINHA DE COMANDO
if [ ! ${outdir} ]; then
	echo "[ERROR] Missing output directory." 1>&2
	exit
fi

# SE ${outdir} NÃO É DIRETÓRIO
if [ ! -d ${outdir} ]; then
	echo "[ERROR] Wrong output directory (${outdir})." 1>&2
	exit
fi

# Número de CORES para o processamento
# ATENÇÃO: Não exceder o limite da máquina
THREADS=${4}

if [ ! ${THREADS} ]; then
	echo "[ERROR] Missing number of threads." 1>&2
	exit
fi

refgtf=${5}
# SE ${refgtf} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 3 NA LINHA DE COMANDO
if [ ! ${refgtf} ]; then
	echo "[ERROR] Missing GTF file." 1>&2
	exit
fi

if [ ! -e "${refgtf}" ]; then
	echo "[ERROR] Not found GTF file (${refgtf})." 1>&2
	exit
fi
absrefgtf=`readlink -f ${refgtf}`

refseq=${6}
# SE ${refseq} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 4 NA LINHA DE COMANDO
if [ ! ${refseq} ]; then
	echo "[ERROR] Missing GENOME fasta file." 1>&2
	exit
fi

if [ ! -e "${refseq}" ]; then
	echo "Not found GENOME fasta file (${refseq})." 1>&2
	exit
fi

# Opção cufflinks/stringtie
assembler=${7}

if [ ! ${assembler} ]; then
	echo "[ERROR] Missing assembler (cufflinks, stringtie, reference)." 1>&2
	exit
fi

if [ "${assembler}" != "cufflinks" ] &&
   [ "${assembler}" != "stringtie" ] &&
   [ "${assembler}" != "reference" ]; then

	echo "[ERROR] Assembler must be \"cufflinks\" or \"stringtie\" (${assembler}). You can also select \"reference\" to not assemble." 1>&2
  	exit
fi

design_file=${8}

if [ ! ${design_file} ]; then
	echo "[ERROR] Missing design file." 1>&2
	exit
fi

if [ ! -e ${design_file} ]; then
	echo "[ERROR] Wrong design file (${design_file})." 1>&2
	exit
fi

design_name=$(basename ${design_file} .txt | sed 's/^[^_]\+_//')

echo "Processing design ($design_name) from file: $design_file ...";

# Contaminantes
contaminants=${9}

if [ ! ${contaminants} ]; then
	echo "[ERROR] Missing contaminant info. Please set \"NA\" here for execution without contaminants." 1>&2
	exit
fi

if [ "${contaminants}" == "NA" ]; then
	contaminants=""
fi

if [ ! "$(ls -A ${indir}/*.fastq* 2> /dev/null)" ]; then
	echo "[ERROR] Input dir if empty of FASTQ data." 1>&2
	exit
fi

./preprocess5.sh "${indir}" "${outdir}" "${THREADS}" ${contaminants}

# gene_info

gene_info=${10}

echo -e "Starting Transcriptome Assembly ..."

# Criação de estrutura de diretórios
	
curdir=`pwd`
	
refseq_abs_path=$(readlink -f ${refseq})

if [ "${aligner}" == "tophat" ]; then

	mkdir -p ${outdir}/tophat_index
	mkdir -p ${outdir}/tophat_out_pe
	mkdir -p ${outdir}/tophat_out_se
	mkdir -p ${outdir}/tophat_out_final
	
	if [ ! -e "${outdir}/tophat_index/genome.fa" ]; then
		cd ${outdir}/tophat_index
		ln -s ${refseq_abs_path} genome.fa
		cd ${curdir}
	fi
		
	if [ ! -e "${outdir}/tophat_index/genome.1.bt2" ]; then
		echo -e "Indexing genome with TopHat2 ..."
		cd ${outdir}/tophat_index
		bowtie2-build 	--threads ${THREADS} \
				genome.fa genome > bowtie2.out.txt 2> bowtie2.err.txt
		cd ${curdir}
	fi
else
	# CASO CONTRÁRIO SERÁ star

	mkdir -p ${outdir}/star_index
	mkdir -p ${outdir}/star_out_pe
	mkdir -p ${outdir}/star_out_se
	mkdir -p ${outdir}/star_out_final
	
	if [ ! -e "${outdir}/star_index/genome.fa" ]; then
		cd ${outdir}/star_index
		ln -s ${refseq_abs_path} genome.fa
		cd ${curdir}
	fi
		
	if [ ! -e "${outdir}/star_index/SAindex" ]; then


		cd ${outdir}/star_index
		
		echo -e "Indexing genome with STAR ..."

		STAR 	--runThreadN ${THREADS} \
	       		--runMode genomeGenerate \
			--genomeFastaFiles genome.fa \
			--genomeDir ./ \
			--sjdbGTFfile ${absrefgtf} \
			--genomeSAindexNbases 12 \
			--sjdbOverhang 100 \
			 > STAR.genomeGenerate.log.out.txt \
			2> STAR.genomeGenerate.log.err.txt

		cd ${curdir}
	fi

fi

if [ "${assembler}" == "cufflinks" ]; then
	mkdir -p ${outdir}/${design_name}/${aligner}_cufflinks
	mkdir -p ${outdir}/${design_name}${aligner}_cuffmerge
elif [ "${assembler}" == "stringtie" ]; then
	mkdir -p ${outdir}/${design_name}/${aligner}_stringtie
	mkdir -p ${outdir}/${design_name}/${aligner}_stringmerge
elif [ "${assembler}" == "reference" ]; then
	echo "[WARNING] Using reference transcriptome (${refgtf})"
else
	echo "[ERROR] Unexpected error!" 1>&2
	exit
fi

mkdir -p ${outdir}/${design_name}/${aligner}_${assembler}_cuffcompare
mkdir -p ${outdir}/${design_name}/${aligner}_${assembler}_cuffquant
mkdir -p ${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm
mkdir -p ${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff

###
# ALIGNMENT STEPS
###

for r1 in `find -L ${outdir}/ -name '*.prinseq.cleaned_1.fastq'`; do 
	
	r2=`echo ${r1} | sed 's/prinseq.cleaned_1.fastq/prinseq.cleaned_2.fastq/'`
	
	if [ ! -e ${r2} ]; then
		echo "[ERROR] Not found R2 (${r2})." 1>&2
		exit
	fi
	
	echo -e "\tFound R1 ($(basename ${r1})) & R2 ($(basename ${r2})) ..."
	
	r1_singletons=`echo ${r1} | sed 's/prinseq.cleaned_1.fastq/prinseq.cleaned_1_singletons.fastq/'`
	r2_singletons=`echo ${r2} | sed 's/prinseq.cleaned_2.fastq/prinseq.cleaned_2_singletons.fastq/'`
	
	if [ ! -e ${r1_singletons} ]; then
		echo "[ERROR] Not found R1 singletons (${r1_singletons})." 1>&2
		exit
	fi

	if   [ ! -e ${r2_singletons} ]; then
		echo "[ERROR] Not found R2 singletons (${r2_singletons})." 1>&2
		exit
	fi
	
	name=`basename ${r1} .fastq | sed 's/.atropos_final.prinseq.cleaned_1//'`
	
	mkdir -p ${outdir}/align_out_final/${name}
	
	if [ "${aligner}" == "tophat" ]; then

		if [ ! -e "${outdir}/tophat_out_pe/${name}/accepted_hits.bam" ]; then

			echo -e "\tTopHat2 alignment (${name}) paired-end reads X genome ..." 

			tophat2 --min-anchor 8 \
				--min-intron-length 30 \
				--max-intron-length 10000 \
				--max-multihits 20 \
				--transcriptome-max-hits 10 \
				--prefilter-multihits \
				--num-threads ${THREADS} \
				--GTF ${refgtf} \
				--transcriptome-index ${outdir}/tophat_index/transcriptome \
				--mate-inner-dist 0 \
				--mate-std-dev 60 \
				--coverage-search \
				--microexon-search \
				--b2-very-sensitive \
				--library-type fr-firststrand \
				--output-dir ${outdir}/tophat_out_pe/${name} \
				--no-sort-bam \
				${outdir}/tophat_index/genome \
				${r1} \
				${r2}	 > ${outdir}/tophat_out_pe/${name}.log.out.txt \
					2> ${outdir}/tophat_out_pe/${name}.log.err.txt
		else
			echo -e "\tFound Tophat2 output for PE (${name})..."
		fi		
		
		if [ ! -e "${outdir}/tophat_out_se/${name}/accepted_hits.bam" ]; then
			
			mkdir -p ${outdir}/tophat_out_se/${name}

			cat ${r1_singletons} ${r2_singletons} > ${outdir}/tophat_out_se/${name}/singletons.fastq
			
			if [ -s "${outdir}/tophat_out_se/${name}/singletons.fastq" ]; then

				echo -e "\tTopHat2 alignment (${name}) singleton reads X genome ..." 

				tophat2 --min-anchor 8 \
					--min-intron-length 30 \
					--max-intron-length 10000 \
					--max-multihits 20 \
					--transcriptome-max-hits 10 \
					--prefilter-multihits \
					--num-threads ${THREADS} \
					--GTF ${refgtf} \
					--transcriptome-index ${outdir}/tophat_index/transcriptome \
					--coverage-search \
					--microexon-search \
					--b2-very-sensitive \
					--library-type fr-firststrand \
					--output-dir ${outdir}/tophat_out_se/${name} \
					--no-sort-bam \
					${outdir}/tophat_index/genome \
					${outdir}/tophat_out_se/${name}/singletons.fastq \
					  > ${outdir}/tophat_out_se/${name}.log.out.txt \
					 2> ${outdir}/tophat_out_se/${name}.log.err.txt
				
				 # Considerar a implementação do TopHat-Recondition https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1058-x
			fi
		else
			echo -e "\tFound Tophat2 output for SE (${name})..."
		fi
		
		if [ ! -e "${outdir}/tophat_out_final/${name}/accepted_hits.bam" ]; then
			
			mkdir -p ${outdir}/tophat_out_final/${name}
		
			if [ -s "${outdir}/tophat_out_pe/${name}/accepted_hits.bam" ]; then

				if [ -s "${outdir}/tophat_out_se/${name}/accepted_hits.bam" ]; then

					echo -e "\tMerging TopHat2 results ..."
		
					samtools view -H ${outdir}/tophat_out_pe/${name}/accepted_hits.bam > ${outdir}/tophat_out_final/${name}/Header.txt
					samtools merge 	-n --threads ${THREADS} \
							-h ${outdir}/tophat_out_final/${name}/Header.txt \
							${outdir}/tophat_out_final/${name}/accepted_hits.bam \
							${outdir}/tophat_out_pe/${name}/accepted_hits.bam \
							${outdir}/tophat_out_se/${name}/accepted_hits.bam \
						 > ${outdir}/tophat_out_final/${name}.log.out.txt \
						2> ${outdir}/tophat_out_final/${name}.log.err.txt

				else
					pe_result_abs_path=$(readlink -f ${outdir}/tophat_out_pe/${name}/accepted_hits.bam)
					cd ${outdir}/tophat_out_final/${name}/
					ln -s ${pe_result_abs_path} accepted_hits.bam
					cd ${curdir}
				fi				
			else
				if [ -s "${outdir}/tophat_out_se/${name}/accepted_hits.bam" ]; then
					se_result_abs_path=$(readlink -f ${outdir}/tophat_out_se/${name}/accepted_hits.bam)
					cd ${outdir}/tophat_out_final/${name}/
					ln -s ${se_result_abs_path} accepted_hits.bam
					cd ${curdir}
				else
					echo -e "[ERROR] Not found any alignment for PE or SE reads." 1>&2
				fi
			fi
		else
			echo -e "\tFound Tophat2 output final (${name})..."
		fi
					
		if [ ! -e "${outdir}/tophat_out_final/${name}/Aligned.sorted.bam" ]; then
			echo -e "\tSorting alignments (${name})..."
			samtools sort 	--threads ${THREADS} \
					-o ${outdir}/tophat_out_final/${name}/Aligned.sorted.bam \
			       		${outdir}/tophat_out_final/${name}/accepted_hits.bam \
					 > ${outdir}/tophat_out_final/${name}/Aligned.sorted.out.txt \
					2> ${outdir}/tophat_out_final/${name}/Aligned.sorted.err.txt
		fi
		
		# SEMPRE VAMOS REMOVER O LINK SIMBÓLICO PARA QUE AO ESCOLHER UM OUTRO 
		# ALINHADOR ELE SEJA SUBSTITUÍDO
		#if [ ! -e "${outdir}/align_out_final/${name}/Aligned.out.bam" ]; then
			rm -f ${outdir}/align_out_final/${name}/Aligned.out.bam
			rm -f ${outdir}/align_out_final/${name}/Aligned.sorted.bam
			if [ -e "${outdir}/tophat_out_final/${name}/accepted_hits.bam" ]; then
				align_final_out=`readlink -f ${outdir}/tophat_out_final/${name}/accepted_hits.bam`
				align_sorted_out=`readlink -f ${outdir}/tophat_out_final/${name}/Aligned.sorted.bam`
				cd ${outdir}/align_out_final/${name}
				ln -s ${align_final_out} Aligned.out.bam
				ln -s ${align_sorted_out} Aligned.sorted.bam
				cd ${curdir}
			else
				echo "[ERROR] Not found Tophat final output (${outdir}/tophat_out_final/${name}/accepted_hits.bam)" 2>&1
				exit
			fi				
		#fi
	else
		# SE NÃO FOR tophat ENTÃO star
		
		if [ ! -e "${outdir}/star_out_pe/${name}/Aligned.out.bam" ]; then

			echo -e "\tSTAR alignment (${name}) paired-end reads X genome ..." 

			mkdir -p ${outdir}/star_out_pe/${name}/
			# Para a execução do cufflinks é necessário: --outSAMstrandField intronMotif e --outFilterIntronMotifs RemoveNoncanonical
			STAR 	--runThreadN ${THREADS} \
				--genomeDir ${outdir}/star_index/ \
				--readFilesIn ${r1} ${r2} \
				--outSAMstrandField intronMotif \
				--outFilterIntronMotifs RemoveNoncanonical \
				--sjdbGTFfile ${refgtf} \
				--outFilterMultimapNmax 20 \
				--outFileNamePrefix ${outdir}/star_out_pe/${name}/ \
				--outSAMtype BAM Unsorted \
				--outFilterType BySJout \
				--outSJfilterReads Unique \
				--alignSJoverhangMin 8 \
				--alignSJDBoverhangMin 1 \
				--outFilterMismatchNmax 999 \
				--outFilterMismatchNoverReadLmax 0.02 \
				--alignIntronMin 30 \
				--alignIntronMax 10000 \
				--alignMatesGapMax 10000 \
			 > ${outdir}/star_out_pe/${name}.log.out.txt \
			2> ${outdir}/star_out_pe/${name}.log.err.txt
		else
			echo -e "\tFound STAR output for PE (${name})..."
		fi

		if [ ! -e "${outdir}/star_out_se/${name}/Aligned.out.bam" ]; then
			
			mkdir -p ${outdir}/star_out_se/${name}

			cat ${r1_singletons} ${r2_singletons} > ${outdir}/star_out_se/${name}/singletons.fastq
			
			if [ -s "${outdir}/star_out_se/${name}/singletons.fastq" ]; then

				echo -e "\tSTAR alignment (${name}) singleton reads X genome ..." 
			
				STAR 	--runThreadN ${THREADS} \
					--genomeDir ${outdir}/star_index/ \
					--readFilesIn ${outdir}/star_out_se/${name}/singletons.fastq \
					--outSAMstrandField intronMotif \
					--outFilterIntronMotifs RemoveNoncanonical \
					--sjdbGTFfile ${refgtf} \
					--outFilterMultimapNmax 20 \
					--outFileNamePrefix ${outdir}/star_out_se/${name}/ \
					--outSAMtype BAM Unsorted \
					--outFilterType BySJout \
					--outSJfilterReads Unique \
					--alignSJoverhangMin 8 \
					--alignSJDBoverhangMin 1 \
					--outFilterMismatchNmax 999 \
					--outFilterMismatchNoverReadLmax 0.02 \
					--alignIntronMin 30 \
					--alignIntronMax 10000 \
					--alignMatesGapMax 10000 \
				 > ${outdir}/star_out_se/${name}.log.out.txt \
				2> ${outdir}/star_out_se/${name}.log.err.txt

			fi
		else
			echo -e "\tFound STAR output for SE (${name})..."
		fi
		
		if [ ! -e "${outdir}/star_out_final/${name}/Aligned.out.bam" ]; then
			
			mkdir -p ${outdir}/star_out_final/${name}
		
			if [ -s "${outdir}/star_out_pe/${name}/Aligned.out.bam" ]; then

				if [ -s "${outdir}/star_out_se/${name}/Aligned.out.bam" ]; then

					echo -e "\tMerging STAR results ..."
		
					samtools view -H ${outdir}/star_out_pe/${name}/Aligned.out.bam > ${outdir}/star_out_final/${name}/Header.txt
					samtools merge 	-n --threads ${THREADS} \
							-h ${outdir}/star_out_final/${name}/Header.txt \
							${outdir}/star_out_final/${name}/Aligned.out.bam \
							${outdir}/star_out_pe/${name}/Aligned.out.bam \
							${outdir}/star_out_se/${name}/Aligned.out.bam \
						 > ${outdir}/star_out_final/${name}.log.out.txt \
						2> ${outdir}/star_out_final/${name}.log.err.txt
						
						samtools sort   -n --threads ${THREADS} \
								${outdir}/star_out_final/${name}/Aligned.out.bam \
								-o ${outdir}/star_out_final/${name}/Aligned.named.out.bam
						
						rm -f ${outdir}/star_out_final/${name}/Aligned.out.bam
						
						mv 	${outdir}/star_out_final/${name}/Aligned.named.out.bam \
							${outdir}/star_out_final/${name}/Aligned.out.bam
						
				else
					pe_result_abs_path=$(readlink -f ${outdir}/star_out_pe/${name}/Aligned.out.bam)
					cd ${outdir}/star_out_final/${name}/
					ln -s ${pe_result_abs_path} Aligned.out.bam
					cd ${curdir}
				fi
			else
				if [ -s "${outdir}/star_out_se/${name}/Aligned.out.bam" ]; then
					se_result_abs_path=$(readlink -f ${outdir}/star_out_se/${name}/Aligned.out.bam)
					cd ${outdir}/star_out_final/${name}/
					ln -s ${se_result_abs_path} Aligned.out.bam
					cd ${curdir}
				else
					echo -e "[ERROR] Not found any alignment for PE or SE reads." 1>&2
				fi
			fi
		else
			echo -e "\tFound STAR output final (${name})..."
		fi
		
		if [ ! -e "${outdir}/star_out_final/${name}/Aligned.sorted.bam" ]; then

			echo -e "\tSorting alignments (${name})..."

			samtools sort 	--threads ${THREADS} \
					-o ${outdir}/star_out_final/${name}/Aligned.sorted.bam \
			       		${outdir}/star_out_final/${name}/Aligned.out.bam \
					 > ${outdir}/star_out_final/${name}/Aligned.sorted.out.txt \
					2> ${outdir}/star_out_final/${name}/Aligned.sorted.err.txt
		fi
		
		# SEMPRE VAMOS REMOVER O LINK SIMBÓLICO PARA QUE AO ESCOLHER UM OUTRO 
		# ALINHADOR ELE SEJA SUBSTITUÍDO
		#if [ ! -e "${outdir}/align_out_final/${name}/Aligned.out.bam" ]; then
			rm -f ${outdir}/align_out_final/${name}/Aligned.out.bam
			rm -f ${outdir}/align_out_final/${name}/Aligned.sorted.bam
			if [ -e "${outdir}/star_out_final/${name}/Aligned.out.bam" ]; then
				align_final_out=`readlink -f ${outdir}/star_out_final/${name}/Aligned.out.bam`
				align_sorted_out=`readlink -f ${outdir}/star_out_final/${name}/Aligned.sorted.bam`
				cd ${outdir}/align_out_final/${name}
				ln -s ${align_final_out} Aligned.out.bam
				ln -s ${align_sorted_out} Aligned.sorted.bam
				cd ${curdir}
			else
				echo "[ERROR] Not found STAR final output (${outdir}/star_out_final/${name}/Aligned.out.bam)" 2>&1
				exit
			fi
		#fi

	fi
	
	mkdir -p ${outdir}/align_out_info/
	
	if [ -e "${outdir}/align_out_final/${name}/Aligned.out.bam" ]; then
	
		if [ ! -e "${outdir}/align_out_info/${name}.${aligner}.out.txt" ]; then

			echo -e "\tGet alignment information (${name}) [${aligner}] ..."
			
			SAM_nameSorted_to_uniq_count_stats.pl ${outdir}/align_out_final/${name}/Aligned.out.bam > ${outdir}/align_out_info/${name}.${aligner}.out.txt 2> ${outdir}/align_out_info/${name}.${aligner}.err.txt
			
		fi


###
# ASSEMBLY STEPS
###

		if [ "${assembler}" == "cufflinks" ]; then
			if [ ! -e "${outdir}/${aligner}_cufflinks/${name}/transcripts.gtf" ]; then

				mkdir -p ${outdir}/${aligner}_cufflinks/${name}

				echo -e "\tAssembly transcriptome (${name}) [cufflinks]\n"

				cufflinks	--output-dir ${outdir}/${aligner}_cufflinks/${name} \
						--num-threads ${THREADS} \
						--GTF ${refgtf} \
						--frag-bias-correct ${refseq} \
						--multi-read-correct \
						--library-type fr-firststrand \
						--frag-len-mean 160 \
						--frag-len-std-dev 60 \
						--total-hits-norm \
						--min-isoform-fraction 0.25 \
						--pre-mrna-fraction 0.15 \
						--min-frags-per-transfrag 10 \
						--junc-alpha 0.001 \
						--small-anchor-fraction 0.08 \
						--overhang-tolerance 8 \
						--min-intron-length 30 \
						--max-intron-length 10000 \
						--trim-3-avgcov-thresh 5 \
						--trim-3-dropoff-frac 0.1 \
						--max-multiread-fraction 0.75 \
						--overlap-radius 25 \
						--3-overhang-tolerance 600 \
						--intron-overhang-tolerance 50 \
						${outdir}/align_out_final/${name}/Aligned.sorted.bam \
					 > ${outdir}/${aligner}_cufflinks/${name}/cufflinks.out.txt \
					2> ${outdir}/${aligner}_cufflinks/${name}/cufflinks.err.txt
			fi

		elif [ "${assembler}" == "stringtie" ]; then
			if [ ! -e "${outdir}/${aligner}_stringtie/${name}/transcripts.gtf" ]; then		
				mkdir -p ${outdir}/${aligner}_stringtie/${name}
					
				echo -e "\tAssembly transcriptome (${name}) [stringtie]\n"

				stringtie	${outdir}/align_out_final/${name}/Aligned.sorted.bam \
						--rf \
						-G ${refgtf} \
						-f 0.25 \
						-m 200 \
						-o ${outdir}/${aligner}_stringtie/${name}/transcripts.gtf \
						-a 8 \
						-j 2 \
						-c 4 \
						-v \
						-g 25 \
						-C ${outdir}/${aligner}_stringtie/${name}/coverages.txt \
						-M 0.9 \
						-p ${THREADS} \
						-A ${outdir}/${aligner}_stringtie/${name}/abundances.txt \
						-e \
						-B \
					 > ${outdir}/${aligner}_stringtie/${name}/stringtie.out.txt \
					2> ${outdir}/${aligner}_stringtie/${name}/stringtie.err.txt				
			fi

		fi


	else
		echo -e "[ERROR] Not found alignment data (${outdir}/align_out_final/${name}/Aligned.out.bam)"	1>&2
	fi

done

rm -f ${outdir}/${design_name}/assembly_GTF_list.txt

transcriptomeref=""

if [ "${assembler}" != "reference" ]; then

	while IFS=$'\t' read -r -a smparray; do
		name=${smparray[0]}
		groupname=${smparray[1]}
		
		transc="${outdir}/${aligner}_${assembler}/${name}/transcripts.gtf"
		
		if [ -e ${transc} ]; then
			echo -e "\tProcessing transcriptome (${name}) biological group (${groupname}) file ${transc} ..."
			echo ${transc} >> ${outdir}/${design_name}/assembly_GTF_list.txt
		else
			echo "[ERROR] Missing transcriptome file ${transc}" 1>&2
			exit		
		fi

	done < "${design_file}"

	if [ "${assembler}" == "cufflinks" ]; then

		if [ ! -e "${outdir}/${design_name}/${aligner}_cuffmerge/merged.gtf" ]; then

			echo -e "\tMerging transcriptomes (${outdir}/${design_name}/assembly_GTF_list.txt) in a transcriptome reference [cuffmerge]"
		
			cuffmerge 	-o ${outdir}/${design_name}/${aligner}_cuffmerge \
					--ref-gtf ${refgtf} \
					--ref-sequence ${refseq} \
					--min-isoform-fraction 0.25 \
					--num-threads ${THREADS} \
					${outdir}/${design_name}/assembly_GTF_list.txt \
				 > ${outdir}/${design_name}/${aligner}_cuffmerge/cuffmerge.out.txt \
				2> ${outdir}/${design_name}/${aligner}_cuffmerge/cuffmerge.err.txt
		fi
		
		transcriptomeref="${outdir}/${design_name}/${aligner}_cuffmerge/merged.gtf"

	elif [ "${assembler}" == "stringtie" ]; then

		if [ ! -e "${outdir}/${design_name}/${aligner}_stringmerge/merged.gtf" ]; then

			echo -e "\tMerging transcriptomes (${outdir}/${design_name}/assembly_GTF_list.txt) in a transcriptome reference [stringtie]\n"
			
			stringtie 	--merge	\
					--rf \
					-G ${refgtf} \
					-o ${outdir}/${design_name}/${aligner}_stringmerge/merged.gtf \
					-m 200 \
					-c 4 \
					-F 4 \
					-T 4 \
					-f 0.25 \
					-e \
					-g 100 \
					${outdir}/${design_name}/assembly_GTF_list.txt \
				 > ${outdir}/${design_name}/${aligner}_stringmerge/stringmerge.out.txt \
				2> ${outdir}/${design_name}/${aligner}_stringmerge/stringmerge.err.txt

		fi
			
		transcriptomeref="${outdir}/${design_name}/${aligner}_stringmerge/merged.gtf"
		
	fi
	
	if [ ! -e "${outdir}/${design_name}/${aligner}_${assembler}_cuffcompare/cuffcmp.combined.gtf" ]; then
		
		echo -e "\tRunning cuffcompare with ${aligner} & ${assembler} transcriptome reference (${transcriptomeref})..."

		cuffcompare 	-r ${refgtf} \
				-s ${refseq} \
				-o ${outdir}/${design_name}/${aligner}_${assembler}_cuffcompare/cuffcmp \
				${transcriptomeref} \
			 > ${outdir}/${design_name}/${aligner}_${assembler}_cuffcompare/cuffcmp.out.txt \
			2> ${outdir}/${design_name}/${aligner}_${assembler}_cuffcompare/cuffcmp.err.txt
		
		# ANOTAÇÃO DOS TRANSCRITOS "TCONS" 
		perl -F"\t" -lane 'INIT { print join("\t","transcript_id","nearest_ref","class_code"); } my ($transcript_id)=$F[8]=~/transcript_id \"([^\"]+)\"/; my ($nearest_ref)=$F[8]=~/nearest_ref \"([^\"]+)\"/; $nearest_ref=~s/^rna-//; $nearest_ref=~s/_[1-9]+$//; $nearest_ref=~s/^rna_gene-//; my ($class_code)=$F[8]=~/class_code \"([^\"]+)\"/; print $transcript_id,"\t",$nearest_ref||"","\t",$class_code||"";' ${outdir}/${design_name}/${aligner}_${assembler}_cuffcompare/cuffcmp.combined.gtf | awk 'NR == 1; NR > 1 {print $0 | "sort -u"}' > ${outdir}/${design_name}/${aligner}_${assembler}_cuffcompare/TCONS.nearest_ref.txt
	fi
			
	transcriptomeref="${outdir}/${design_name}/${aligner}_${assembler}_cuffcompare/cuffcmp.combined.gtf"
else
	transcriptomeref="${refgtf}"

	cd ${outdir}/${design_name}/${aligner}_${assembler}_cuffcompare/
	if [ ! -e ./cuffcmp.combined.gtf ]; then
		ln -s ${absrefgtf} ./cuffcmp.combined.gtf
	fi
	cd ${curdir}
fi

# LISTA DE VALORES NÃO REDUNDANTES (NOME DO GRUPO BIOLÓGICO)
# Ex.: (CONTROL TEST)
biogroup_label=()

declare -A BIOLABEL;

while IFS=$'\t' read -r -a smparray; do
	name=${smparray[0]}
	groupname=${smparray[1]}

	bamfile="${outdir}/align_out_final/${name}/Aligned.sorted.bam"

	if [ ! -e ${bamfile} ]; then
		echo "[ERROR] Missing ${bamfile}." 1>&2
		exit
	fi


	if [ ! -e "${outdir}/${design_name}/${aligner}_${assembler}_cuffquant/${name}/abundances.cxb" ]; then

		echo -e "\tRunning cuffquant using sample ${name} as using ${aligner} & ${assembler} (${transcriptomeref}) ..."
		mkdir -p ${outdir}/${design_name}/${aligner}_${assembler}_cuffquant/${name}
			
		cuffquant 	--output-dir ${outdir}/${design_name}/${aligner}_${assembler}_cuffquant/${name} \
				--frag-bias-correct ${refseq} \
				--multi-read-correct \
				--num-threads ${THREADS} \
				--library-type fr-firststrand \
				--frag-len-mean 160 \
				--frag-len-std-dev 60 \
				--max-bundle-frags 9999999 \
				--max-frag-multihits 20 \
				${transcriptomeref} \
				${bamfile} \
			 > ${outdir}/${design_name}/${aligner}_${assembler}_cuffquant/${name}/cuffquant.log.out.txt \
			2> ${outdir}/${design_name}/${aligner}_${assembler}_cuffquant/${name}/cuffquant.log.err.txt

	fi
	
	biogroup_label=($(printf "%s\n" ${biogroup_label[@]} ${groupname} | sort -u ))
	
	BIOLABEL[${name}]="${groupname}"

done < "${design_file}"

biogroup_files=()

echo -e "\tCollecting Expression Data from cuffquant output (*.cxb) ..."


for label in ${biogroup_label[@]}; do
	echo -e "\t\tCollecting .cxb files for ${label} ..."
	group=()

	for k in ${!BIOLABEL[@]}; do
		if [ ${BIOLABEL[${k}]} == ${label} ]; then
			cxbfile="${outdir}/${design_name}/${aligner}_${assembler}_cuffquant/${k}/abundances.cxb"
			echo -e "\t\t\tFound ${cxbfile}"
			group=(${group[@]} "${cxbfile}")
		fi			
	done
	biogroup_files=(${biogroup_files[@]} $(IFS=, ; echo "${group[*]}") )
done

echo -e "Starting Gene Expression Analysis ..."
echo -e "\t\tLabels.: " $(IFS=, ; echo "${biogroup_label[*]}")
echo -e "\t\tFiles..: " ${biogroup_files[*]}

if [ ! -e "${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/isoforms.count_table" ]; then

	echo -e "\t\t\tGenerating abundance matrices (cuffnorm) ..."

	cuffnorm 	--output-dir ${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm \
			--labels $(IFS=, ; echo "${biogroup_label[*]}") \
			--num-threads ${THREADS} \
			--library-type fr-firststrand \
			--library-norm-method geometric \
			--output-format simple-table \
			${transcriptomeref} \
			${biogroup_files[*]} \
			 > ${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/cuffnorm.log.out.txt \
			2> ${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/cuffnorm.log.err.txt

fi

if [ ! -e "${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/isoforms.raw_count_table.txt" ]; then
	de-normalize-cuffnorm.R --in=${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/isoforms.count_table \
				--st=${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/samples.table \
				--out=${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/isoforms.raw_count_table.txt
	 > ${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/de-normalize-cuffnorm.isoforms.out.txt \
	2> ${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/de-normalize-cuffnorm.isoforms.err.txt

fi

if [ ! -e "${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/genes.raw_count_table.txt" ]; then
	de-normalize-cuffnorm.R --in=${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/genes.count_table \
				--st=${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/samples.table \
				--out=${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/genes.raw_count_table.txt
	 > ${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/de-normalize-cuffnorm.genes.out.txt \
	2> ${outdir}/${design_name}/${aligner}_${assembler}_cuffnorm/de-normalize-cuffnorm.genes.err.txt

fi

if [ ! -e "${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/isoform_exp.diff" ]; then
	echo -e "\t\t\tAnalysing differential expression (cuffdiff) ..."

	cuffdiff 	--output-dir ${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff \
			--labels $(IFS=, ; echo "${biogroup_label[*]}") \
			--frag-bias-correct ${refseq} \
			--multi-read-correct \
			--num-threads ${THREADS} \
			--library-type fr-firststrand \
			--frag-len-mean 160 \
			--frag-len-std-dev 60 \
			--max-bundle-frags 9999999 \
			--max-frag-multihits 20 \
			--total-hits-norm \
			--min-reps-for-js-test 2 \
			--library-norm-method geometric \
			--dispersion-method per-condition \
			--min-alignment-count 10 \
			${transcriptomeref} \
			${biogroup_files[*]} \
			 > ${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/cuffdiff.log.out.txt \
			2> ${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/cuffdiff.log.err.txt
fi

if [ ! -e "${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/isoform_exp.diff.annot.txt" ]; then

	if [ -e "${outdir}/${design_name}/${aligner}_${assembler}_cuffcompare/TCONS.nearest_ref.txt" ]; then

		echo "Annotating isoform_exp.diff ..."

		mergeR.R	--x=${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/isoform_exp.diff \
				--by.x="test_id" \
				--y=${outdir}/${design_name}/${aligner}_${assembler}_cuffcompare/TCONS.nearest_ref.txt \
				--by.y="transcript_id" \
				--all.x \
				--print.out.label \
				--out=${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/isoform_exp.diff.annot.txt \
				 > ${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/mergeR.isoform.log.out.txt \
				2> ${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/mergeR.isoform.log.err.txt
	fi
fi

if [ ! -e "${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/gene_exp.diff.annot.txt" ]; then

	if [ ${gene_info} ]; then

		if [ ! -e ${gene_info} ]; then
			echo "[ERROR] Wrong gene_info file (${gene_info})." 1>&2
			exit
		fi
		

		echo "Annotating gene_exp.diff ..."

		splitteR.R	--x="${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/gene_exp.diff" \
				--col.x="gene" \
				--by.x="," \
				--out="${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/gene_exp.diff.spplitted.txt" \
				 > ${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/mergeR.gene.log.out.txt \
				2> ${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/mergeR.gene.log.err.txt
		
		mergeR.R	--x=${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/gene_exp.diff.spplitted.txt \
				--by.x="gene_id" \
				--y=${gene_info} \
				--by.y="GeneID" \
				--all.x \
				--print.out.label \
				--out=${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/gene_exp.diff.annot.txt \
				 > ${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/mergeR.gene.log.out.txt \
				2> ${outdir}/${design_name}/${aligner}_${assembler}_cuffdiff/mergeR.gene.log.err.txt
	
	fi
fi

echo "Obtaining count's information about preprocess and alignments ..."

mkdir -p ${outdir}/backup
./stats_preprocess.sh ${indir} ${outdir} ${outdir}/backup > ${outdir}/STATS.txt
./getAlnTab.pl -a star -s ${outdir}/STATS.txt -i ${outdir}/align_out_info &> ${outdir}/ALIGN_STATS.txt

