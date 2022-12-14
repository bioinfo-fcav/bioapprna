#!/bin/bash

# input - diretório contendo os arquivos de entrada no formato .fastq
input=$1

# validação do parâmetro "input"
if [ ! ${input} ]
then   
        echo "[ERROR] Missing input directory." 1>&2
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "[ERROR] Wrong input directory (${input})." 1>&2
                exit
        fi
fi

# output - diretório para armazenar o resultado do processo de montagem
output=$2

# validação do parâmetro "output"
if [ ! ${output} ]
then   
        echo "[ERROR] Missing output directory." 1>&2
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "[ERROR] Wrong output directory (${output})." 1>&2
                exit
        fi
fi

# Número de CORES para o processamento
# ATENÇÃO: Não exceder o limite da máquina
THREADS=$3

if [ ! ${THREADS} ]; then
	echo "[ERROR] Missing number of threads." 1>&2
	exit
fi

# Quantidade de memória para o processamento com Jellyfish
# ATENÇÃO: Não exceder o limite da máquina
MEM=$4

if [ ! ${MEM} ]; then
	echo "[ERROR] Missing memory." 1>&2
	exit
fi

reftranscriptomefa=$5

if [ ! ${reftranscriptomefa} ]; then
	echo "[ERROR] Missing reference transcriptome or NA." 1>&2
	exit
else
	if [ ${reftranscriptomefa} == "NA" ]; then
		echo "[WARNING] Not using a transcriptome reference for Detonate evaluation" 1>&2
	else
		if [ ! -e ${reftranscriptomefa} ]
		then   
			echo "[ERROR] Wrong reference transcriptome fasta file (${reftranscriptomefa})." 1>&2
			exit
		fi
	fi		
fi

refproteomefa=$6

if [ ! ${refproteomefa} ]; then
	echo "[ERROR] Missing reference proteome or NA." 1>&2
	exit
else
	if [ ${refproteomefa} == "NA" ]; then
		echo "[WARNING] Not using a proteome reference for Transrate evaluation" 1>&2
	else
		if [ ! -e ${refproteomefa} ]
		then   
			echo "[ERROR] Wrong reference proteome fasta file (${refproteomefa})." 1>&2
			exit
		fi
	fi		
fi

refpfamhmm=$7

if [ ! ${refpfamhmm} ]; then
	echo "[ERROR] Missing reference Pfam HMM file or NA." 1>&2
	exit
else
	if [ ${refpfamhmm} == "NA" ]; then
		echo "[WARNING] Not using a Pfam reference for TransDecoder evaluation" 1>&2
	else
		if [ ! -e ${refpfamhmm} ]
		then   
			echo "[ERROR] Wrong reference Pfam hmm file (${refpfamhmm})." 1>&2
			exit
		fi
	fi		
fi



reference_sample=$8

if [ ${reference_sample} ]; then
	reference_sample_param="--reference_sample ${reference_sample}"
else 
	reference_sample_param=""
fi

###
# Arquivos e diretórios de saída (output) 
#

basedir_out="${output}/"

renamed_out="${basedir_out}/renamed"

trinity_out="${basedir_out}/trinity_assembled"

mkdir -p ${renamed_out}
mkdir -p ${trinity_out}

left=()
left_singleton=()

right=()
right_singleton=()

echo "Performing renaming step ..."

for fastq in `ls ${input}/*.fastq`; do
	# obtendo nome do arquivo 
	fastqbn=`basename ${fastq}`;
	if [[ ! $fastqbn =~ \.bad_ ]]; then
		renamed_fastq="${renamed_out}/${fastqbn}"
		if [ ! -e ${renamed_fastq} ]; then
			echo -e "\tRenaming ${fastqbn} ..."
			if [[ ${fastqbn} =~ _1[\._] ]]; then
				awk '{ if (NR%4==1) { if ($1!~/\/1$/) { print $1"/1" } else { print $1 } } else if (NR%4==3) { print "+" } else { print $1 } }' ${fastq} > ${renamed_fastq}
			elif [[ ${fastqbn} =~ _2[\._]  ]]; then
				awk '{ if (NR%4==1) { if ($1!~/\/2$/) { print $1"/2" } else { print $1 } } else if (NR%4==3) { print "+" } else { print $1 } }' ${fastq} > ${renamed_fastq}
			else 
				echo "Warning: ${fastqbn} discarded!"
			fi
		fi
		
		if [[ ${fastqbn} =~ _1[\._] ]]; then
			if [[ ${fastqbn} =~ singletons ]]; then
				if [ -s ${renamed_fastq} ]; then
					left_singleton=($(printf "%s\n" ${left_singleton[@]} ${renamed_fastq} | sort -u ))
				fi
			else
				left=($(printf "%s\n" ${left[@]} ${renamed_fastq}  | sort -u ))
			fi
		elif [[ ${fastqbn} =~ _2[\._] ]]; then
			if [[ ${fastqbn} =~ singleton ]]; then
				if [ -s ${renamed_fastq} ]; then
					right_singleton=($(printf "%s\n" ${right_singleton[@]} ${renamed_fastq}  | sort -u ))
				fi
			else
				right=($(printf "%s\n" ${right[@]} ${renamed_fastq}  | sort -u ))
			fi
		else
			echo "Warning: ${fastqbn} discarded!"
		fi
	fi
done


#for l in ${left[@]}; do
#	echo -e "L: ${l}";
#done
#
#for r in ${right[@]}; do
#	echo -e "R: ${r}";
#done
#
#for ls in ${left_singleton[@]}; do
#	echo -e "LS: ${ls}";
#done
#
#for rs in ${right_singleton[@]}; do
#	echo -e "RS: ${rs}";
#done

trinity_out_bn=`echo ${trinity_out} | sed 's/\/+$//'`

echo "Checking existance of ${trinity_out_bn}.Trinity.fasta"

if [ ! -e ${trinity_out_bn}.Trinity.fasta ]; then
	
	echo -e "Assembling step (Trinity) ..."
	
	rm -fr ${trinity_out}
	mkdir -p ${trinity_out}

	Trinity --output ${trinity_out} \
		--seqType fq \
		--max_memory ${MEM} \
		--CPU ${THREADS} \
		--SS_lib_type RF \
		--min_kmer_cov 3 \
		--min_glue 2 \
		--max_chrysalis_cluster_size 30 \
		--group_pairs_distance 600 \
		--path_reinforcement_distance 30 \
		--min_per_id_same_path 95 \
		--max_diffs_same_path 5 \
		--max_internal_gap_same_path 10 \
		--min_contig_length 300 \
		--normalize_max_read_cov 200 \
		--left $(IFS=, ; echo "${left[*]},${left_singleton[*]}") \
		--right $(IFS=, ; echo "${right[*]},${right_singleton[*]}") \
		 > ${trinity_out}/Trinity.log.out.txt \
		2> ${trinity_out}/Trinity.log.err.txt
else	
	echo "Assembly found: ${trinity_out_bn}.Trinity.fasta"
fi

trinity_fasta="${trinity_out_bn}.Trinity.fasta"
trinity_trans_map="${trinity_out_bn}.Trinity.fasta.gene_trans_map"


if [ ${reftranscriptomefa} != 'NA' ]; then

	echo "Calculating max read length ..."

	#maxreadlength=`cat ${input}/*_1.fastq | perl -lane ' INIT { our $max=0; } if ($.%4==2) { if (length($_)>$max) { $max=length($_);} } END {print $max;}'`
	
	#echo "Found maximum read length: ${maxreadlength}"

	echo "Calling Detonate.sh ..."

	./Detonate.sh ${THREADS} ${input} ${output} ${reftranscriptomefa} ${output} ${maxreadlength}
fi

if [ ${refproteomefa} != 'NA' ]; then

	echo "Calling Transrate.sh ..."

	./Transrate.sh ${THREADS} ${renamed_out} ${trinity_out_bn}.Trinity.fasta ${refproteomefa} ${output} > Transrate.txt

	trinity_fasta=`find ${output} -type f -name 'good.*.fasta' -print 2>/dev/null | grep -v 'single_component_bad'`
	trinity_trans_map=`find ${output} -type f -name 'good.*.gene_trans_map' -print 2>/dev/null`
fi


abundance_out="${output}/Abundance"

mkdir -p ${abundance_out}

echo -e "id\tname\tgroup" > ${abundance_out}/groups.txt

rm -f ${abundance_out}/samples.txt
rm -f ${abundance_out}/quant_files.txt

for l in ${left[@]}; do
	#echo ${l}
	repname=`basename ${l} | sed 's/\..*$//'`
	#echo ${repname}
	condname=`echo ${repname} | sed 's/_B[0-9]\+//'`
	#echo ${condname}
	r=`echo ${l} | sed 's/_1.fastq/_2.fastq/'`
	right=(${right[@]} ${r})

	echo -e "${condname}\t${abundance_out}/${repname}\t${l}\t${r}" >> ${abundance_out}/samples.txt
	echo -e "${condname}\t${repname}\t${l}\t${r}" >> ${abundance_out}/samples_DE.txt
	echo -e "${abundance_out}/${repname}/quant.sf" >> ${abundance_out}/quant_files.txt
done

echo -e "Using assembly:\t${trinity_fasta}\n\t\t${trinity_trans_map}\n"

if [ ! ${TRINITY_HOME} ]; then
	echo "[ERROR] Missing TRINITY_HOME environmental variable." 1>&2
	exit
fi

if [ ! -e "${abundance_out}/abund.gene.counts.matrix" ]; then

	echo "Calling align_and_estimate_abundance.pl ..."

	${TRINITY_HOME}/util/align_and_estimate_abundance.pl 	--transcripts ${trinity_fasta} \
								--SS_lib_type RF \
								--est_method salmon \
								--samples_file ${abundance_out}/samples.txt \
								--gene_trans_map ${trinity_trans_map} \
								--prep_reference \
								--thread_count ${THREADS} \
								--seqType fq \
								--output_dir ${abundance_out} \
				 > ${abundance_out}/align_and_estimate_abundance.log.out.txt \
				2> ${abundance_out}/align_and_estimate_abundance.log.err.txt

	echo "Constructing abundance matrix ..."


	${TRINITY_HOME}/util/abundance_estimates_to_matrix.pl	--est_method salmon \
								--gene_trans_map ${trinity_trans_map} \
								--name_sample_by_basedir \
								--cross_sample_norm none \
								--quant_files ${abundance_out}/quant_files.txt \
								--out_prefix ${abundance_out}/abund \
				 > ${abundance_out}/abundance_estimates_to_matrix.log.out.txt \
				2> ${abundance_out}/abundance_estimates_to_matrix.log.err.txt

fi


de_out="${output}/DE"
mkdir -p ${de_out}

de_results=(`find output/DE/ -name '*.DE_results' -print`)

if [ "${#de_results[@]}" -eq "0"  ]; then

	echo "Performing Differential Expression analysis with DESeq2 ..."

	${TRINITY_HOME}/Analysis/DifferentialExpression/run_DE_analysis.pl	--matrix ${abundance_out}/abund.isoform.counts.matrix \
										--method DESeq2 \
										--samples_file ${abundance_out}/samples_DE.txt \
										--min_reps_min_cpm 2,1 \
										--output ${de_out} \
										${reference_sample_param} \
					 > ${de_out}/run_DE_analysis.log.out.txt \
					2> ${de_out}/run_DE_analysis.log.err.txt
fi

td_out="${output}/TransDecoder"
mkdir -p ${td_out}

#TransDecoder.LongOrfs	-t ${trinity_fasta} \
#			--gene_trans_map ${trinity_trans_map} \
#			-m 90 \
#			--output_dir ${td_out} \
#			 > ${td_out}/TransDecoder.LongOrfs.log.out.txt \
#			2> ${td_out}/TransDecoder.LongOrfs.log.err.txt

proteome_hits_param=""
if [ ${refproteomefa} != 'NA' ]; then
	
	proteomerefidx="${td_out}/proteome"
	
	if [ ! -e "${proteomerefidx}.pin" ]; then
		
		echo "Make blast index for reference proteome (${refproteomefa}) ..."

		makeblastdb	-in ${refproteomefa} \
				-dbtype prot \
				-out ${proteomerefidx}	\
				 > ${td_out}/makeblastdb.log.out.txt \
				2> ${td_out}/makeblastdb.log.err.txt
	fi
	
	if [ ! -e "${td_out}/proteome.blastp.outfmt6" ]; then
		
		echo "Running blastp: predicted ORFs x reference proteome (${refproteomefa}) ..."

		blastp	-query ${td_out}/longest_orfs.pep \
			-db ${proteomerefidx} \
			-num_threads ${THREADS} \
			-max_target_seqs 1 \
			-outfmt 6 \
			-qcov_hsp_perc 90 \
			-evalue 1e-5 \
			 > ${td_out}/proteome.blastp.outfmt6 \
			2> ${td_out}/blastp.log.err.txt
	fi

	perl -F"\t" -lane 'if ($F[2]>=25) { print $_; }' ${td_out}/proteome.blastp.outfmt6 > ${td_out}/proteome.blastp.selected.outfmt6

	proteome_hits_param="--retain_blastp_hits ${td_out}/proteome.blastp.selected.outfmt6 "

fi


pfam_hits_param=""

if [ ${refpfamhmm} != 'NA' ]; then

	if [ ! -e "${refpfamhmm}.h3m" ]; then
		
		echo "Preparing Pfam database (hmmpress) ..."

		hmmpress ${refpfamhmm} > ${refpfamhmm}.hmmpress.log.out.txt 2> ${refpfamhmm}.hmmpress.log.err.txt
	fi


	# https://github.com/TransDecoder/TransDecoder/issues/94
	# As per this site, hmmsearch is significantly faster and parallelizes better than hmmscan (at least when searching ~70,000 transcripts against Pfam-A)


	if [ ! -e "${td_out}/pfam.domtblout" ]; then

		echo "Running hmmsearch predicted ORFs x Pfam ..."

		hmmsearch 	--cpu ${THREADS} -o /dev/null --domtblout ${td_out}/hmmsearch.tmp \
			${refpfamhmm} ${td_out}/longest_orfs.pep \
			 > ${td_out}/hmmsearch.log.out.txt \
			2> ${td_out}/hmmsearch.log.err.txt

		awk 'BEGIN {OFS=FS=" "} NR<=3{print}; NR>3{tmp=$1; $1=$4; $4=tmp; tmp=$2; $2=$5; $5=tmp; print}' ${td_out}/hmmsearch.tmp > ${td_out}/pfam.domtblout
				

	fi
		
	pfam_hits_param="--retain_pfam_hits ${td_out}/pfam.domtblout" 
fi


transdecoder_input_bn=`basename ${trinity_fasta}`

if [ ! -e "${td_out}/${transdecoder_input_bn}.transdecoder.pep" ]; then
	TransDecoder.Predict	-t ${trinity_fasta} \
				--retain_long_orfs_mode dynamic \
				${pfam_hits_param} \
				${proteome_hits_param} \
				--single_best_only \
				--output_dir ${td_out} \
				-T 1000 \
				 > ${td_out}/TransDecoder.Predict.log.out.txt \
				2> ${td_out}/TransDecoder.Predict.log.err.txt


	mv ./${transdecoder_input_bn}.transdecoder.* ${td_out}/
fi

