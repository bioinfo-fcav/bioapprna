#!/bin/bash

# definição do número de threads
num_threads=$1

# input - diretório contendo os arquivos de entrada no formato .fastq
input=$2

# reftrinity - Arquivo da montagem de novo Trinity.fa
reftrinity=$3

# refprotfaa - Arquivo que contem o proteoma de referencia mais proximo para ser utilizado como base
refprotfaa=$4

# output - diretório para armazenar o resultado do processo de montagem
output=$5

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
        echo "Missing input directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input directory ${input}"
                exit
        fi
fi

if [ ! ${reftrinity} ]
then
        echo "Missing assembled fasta file"
        exit
else
        if [ ! -e ${reftrinity} ]
        then
                echo "Wrong assembled fasta file  ${reftrinity}"
                exit
        fi
fi


if [ ! ${refprotfaa} ]
then
        echo "Missing reference protein fasta file"
        exit
else
        if [ ! -e ${refprotfaa} ]
        then
                echo "Wrong reference protein fasta file  ${refprotfaa}"
                exit
        fi
fi


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


mkdir -p ${output}/Transrate

basedir_out="${output}/Transrate"


leftreads=()
rightreads=()

leftreads=( $( find ${input} -name '*_1.fastq' ) )
for r in ${leftreads[@]}; do
	r=`echo ${r} | sed 's/_1\./_2\./'`
	rightreads=(${rightreads[@]} ${r})
done


bn=`basename ${reftrinity} .fasta`

if [ ! -e "${basedir_out}/${bn}/good.${bn}.fasta" ]; then

	echo -e "Runnig Transrate ..."

	transrate --reference ${refprotfaa} \
	          --assembly ${reftrinity} \
	          --left $(IFS=, ; echo "${leftreads[*]}") \
	          --right $(IFS=, ; echo "${rightreads[*]}") \
	          --threads ${num_threads} \
		  --output ${basedir_out}
	
	grep '^>' ${basedir_out}/${bn}/good.${bn}.fasta | sed 's/^>//' > ${basedir_out}/${bn}/good.${bn}.isoform.txt
	cat ${basedir_out}/${bn}/good.${bn}.isoform.txt | sed 's/_i.*//' > ${basedir_out}/${bn}/good.${bn}.gene.txt

	paste ${basedir_out}/${bn}/good.${bn}.gene.txt ${basedir_out}/${bn}/good.${bn}.isoform.txt > ${basedir_out}/${bn}/good.${bn}.fasta.gene_trans_map

fi
