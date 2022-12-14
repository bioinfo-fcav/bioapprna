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
#  Copyright (C) 2019  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade Estadual Paulista "Júlio de Mesquita Filho" (UNESP)
#  Faculdade de Ciências Agrárias e Veterinárias (FCAV)
#  Laboratório de Bioinformática (LB)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://www.fcav.unesp.br 
#

indir=$1

# SE ${indir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO NA LINHA DE COMANDO
if [ ! ${indir} ]; then
	echo "[ERROR] Missing input directory." 1>&2
	exit
fi

# SE ${indir} NÃO É DIRETÓRIO
if [ ! -d ${indir} ]; then
	echo "[ERROR] Wrong input directory (${indir})." 1>&2
	exit
fi

outdir=$2

# SE ${outdir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO NA LINHA DE COMANDO
if [ ! ${outdir} ]; then
	echo "[ERROR] Missing output directory." 1>&2
	exit
fi

# SE ${outdir} NÃO É DIRETÓRIO
if [ ! -d ${outdir} ]; then
	echo "[ERROR] Wrong output directory (${outdir})." 1>&2
	exit
fi

backupdir=$3

if [ ${backupdir} ]; then
	if [ ! -d ${backupdir} ]; then
		echo "[ERROR] Wrong backup directory (${backupdir})." 1>&2
		exit
	fi
fi

echo -e "Sample\tRaw R1\tRaw R2\tFree of adapter R1\tFree of adapter R2\tHigh quality PE R1\tHigh quality PE R2\tHigh quality SE R1\tHigh quality SE R2\tFree of contamination PE R1\tFree of contamination PE R2\tFree of contamination SE R1\tFree of contamination SE R2"
for infiler1 in `ls ${indir}/*_R1.fastq ${indir}/*_R1.fastq.gz 2>/dev/null`; do
	
	infiler2=$(echo ${infiler1} | sed 's/_R1\.fastq/_R2.fastq/')
	
	SMP=`basename $( basename ${infiler1} .gz ) _R1.fastq`;
	
	if [ ${backupdir} ]; then
		if [ -e ${backupdir}/${SMP}.txt ]; then
			source "${backupdir}/${SMP}.txt"
		fi
	fi

	if [ ! ${R1} ]; then
	
		if [[ ${infiler1} =~ .gz$ ]]; then
			R1=`echo "scale=2; $(gunzip -c ${indir}/${SMP}_R1.fastq.gz | wc -l)/4" | bc`;
		else
			R1=`echo "scale=2; $(cat ${indir}/${SMP}_R1.fastq | wc -l)/4" | bc`;
		fi
		
		if [ ${backupdir} ]; then
			echo "R1=${R1}" >> ${backupdir}/${SMP}.txt
		fi
	fi

	if [ ! ${R2} ]; then

		if [[ ${infiler2} =~ .gz$ ]]; then
		
			R2=`echo "scale=2; $(gunzip -c ${indir}/${SMP}_R2.fastq.gz | wc -l)/4" | bc`;
		else
			R2=`echo "scale=2; $(cat ${indir}/${SMP}_R2.fastq | wc -l)/4" | bc`;
		fi

		if [ ${backupdir} ]; then
			echo "R2=${R2}" >> ${backupdir}/${SMP}.txt
		fi
	fi

	if [ ! ${ATROPOS_PE1} ]; then
		ATROPOS_PE1=`echo "scale=2; $(cat ${outdir}/processed/atropos/${SMP}_R1.atropos_final.fastq | wc -l)/4" | bc`;
		if [ ${backupdir} ]; then
			echo "ATROPOS_PE1=${ATROPOS_PE1}" >> ${backupdir}/${SMP}.txt
		fi
	fi		
	if [ ! ${ATROPOS_PE2} ]; then
		ATROPOS_PE2=`echo "scale=2; $(cat ${outdir}/processed/atropos/${SMP}_R2.atropos_final.fastq | wc -l)/4" | bc`;
		if [ ${backupdir} ]; then
			echo "ATROPOS_PE2=${ATROPOS_PE2}" >> ${backupdir}/${SMP}.txt
		fi
	fi
	
	if [ ! ${PRINSEQ_PE1} ]; then
		PRINSEQ_PE1=`echo "scale=2; $(cat ${outdir}/processed/prinseq/${SMP}.atropos_final.prinseq_1.fastq | wc -l)/4" | bc`;
		if [ ${backupdir} ]; then
			echo "PRINSEQ_PE1=${PRINSEQ_PE1}" >> ${backupdir}/${SMP}.txt
		fi
	fi
	if [ ! ${PRINSEQ_PE2} ]; then
		PRINSEQ_PE2=`echo "scale=2; $(cat ${outdir}/processed/prinseq/${SMP}.atropos_final.prinseq_2.fastq | wc -l)/4" | bc`;
		if [ ${backupdir} ]; then
			echo "PRINSEQ_PE2=${PRINSEQ_PE2}" >> ${backupdir}/${SMP}.txt
		fi
	fi
	if [ ! ${PRINSEQ_SE1} ]; then
		PRINSEQ_SE1=`echo "scale=2; $(cat ${outdir}/processed/prinseq/${SMP}.atropos_final.prinseq_1_singletons.fastq | wc -l)/4" | bc`;
		if [ ${backupdir} ]; then
			echo "PRINSEQ_SE1=${PRINSEQ_SE1}" >> ${backupdir}/${SMP}.txt
		fi
	fi
	if [ ! ${PRINSEQ_SE2} ]; then
		PRINSEQ_SE2=`echo "scale=2; $(cat ${outdir}/processed/prinseq/${SMP}.atropos_final.prinseq_2_singletons.fastq | wc -l)/4" | bc`;
		if [ ${backupdir} ]; then
			echo "PRINSEQ_SE2=${PRINSEQ_SE2}" >> ${backupdir}/${SMP}.txt
		fi
	fi

	if [ ! -L ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_1.fastq ]; then

		if [ ! -e ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_1.fastq ]; then
			echo "[ERROR] Missing ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_1.fastq" 1>&2
#			exit
		fi
		if [ ! -e ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_2.fastq ]; then
			echo "[ERROR] Missing ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_2.fastq" 1>&2
#			exit
		fi
		if [ ! -e ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_1_singletons.fastq ]; then
			echo "[ERROR] Missing ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_1_singletons.fastq" 1>&2
#			exit
		fi
		if [ ! -e ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_2_singletons.fastq ]; then
			echo "[ERROR] Missing ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_2_singletons.fastq" 1>&2
#			exit
		fi

		if [ ! ${CLEANED_PE1} ]; then
			CLEANED_PE1=`echo "scale=2; $(cat ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_1.fastq | wc -l)/4" | bc`;
			if [ ${backupdir} ]; then
				echo "CLEANED_PE1=${CLEANED_PE1}" >> ${backupdir}/${SMP}.txt
			fi
		fi
		if [ ! ${CLEANED_PE2} ]; then
			CLEANED_PE2=`echo "scale=2; $(cat ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_2.fastq | wc -l)/4" | bc`;
			if [ ${backupdir} ]; then
				echo "CLEANED_PE2=${CLEANED_PE2}" >> ${backupdir}/${SMP}.txt
			fi
		fi
		if [ ! ${CLEANED_SE1} ]; then
			CLEANED_SE1=`echo "scale=2; $(cat ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_1_singletons.fastq | wc -l)/4" | bc`;
			if [ ${backupdir} ]; then
				echo "CLEANED_SE1=${CLEANED_SE1}" >> ${backupdir}/${SMP}.txt
			fi
		fi
		if [ ! ${CLEANED_SE2} ]; then
			CLEANED_SE2=`echo "scale=2; $(cat ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_2_singletons.fastq | wc -l)/4" | bc`;
			if [ ${backupdir} ]; then
				echo "CLEANED_SE2=${CLEANED_SE2}" >> ${backupdir}/${SMP}.txt
			fi
		fi
	else
		if [ ! ${CLEANED_PE1} ]; then
			if [ ! -L ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_1.fastq ]; then
				echo "[ERROR] Missing ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_2.fastq" 1>&2
#				exit
			fi
			CLEANED_PE1=${PRINSEQ_PE1}
			if [ ${backupdir} ]; then
				echo "CLEANED_PE1=${CLEANED_PE1}" >> ${backupdir}/${SMP}.txt
			fi
		fi
		if [ ! ${CLEANED_PE2} ]; then
			if [ ! -L ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_2.fastq ]; then
				echo "[ERROR] Missing ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_2.fastq" 1>&2
#				exit
			fi
			CLEANED_PE2=${PRINSEQ_PE2}
			if [ ${backupdir} ]; then
				echo "CLEANED_PE2=${CLEANED_PE2}" >> ${backupdir}/${SMP}.txt
			fi
		fi
		if [ ! ${CLEANED_SE1} ]; then
			if [ ! -L ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_1_singletons.fastq ]; then
				echo "[ERROR] Missing ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_1_singletons.fastq" 1>&2
#				exit
			fi
			CLEANED_SE1=${PRINSEQ_SE1}
			if [ ${backupdir} ]; then
				echo "CLEANED_SE1=${CLEANED_SE1}" >> ${backupdir}/${SMP}.txt
			fi
		fi			
		if [ ! ${CLEANED_SE2} ]; then
			if [ ! -L ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_2_singletons.fastq ]; then
				echo "[ERROR] Missing ${outdir}/processed/cleaned/${SMP}.atropos_final.prinseq.cleaned_2_singletons.fastq" 1>&2
#				exit
			fi
			CLEANED_SE2=${PRINSEQ_SE2}
			if [ ${backupdir} ]; then
				echo "CLEANED_SE2=${CLEANED_SE2}" >> ${backupdir}/${SMP}.txt
			fi
		fi
	fi

	echo -e "${SMP}\t${R1}\t${R2}\t${ATROPOS_PE1}\t${ATROPOS_PE2}\t${PRINSEQ_PE1}\t${PRINSEQ_PE2}\t${PRINSEQ_SE1}\t${PRINSEQ_SE2}\t${CLEANED_PE1}\t${CLEANED_PE2}\t${CLEANED_SE1}\t${CLEANED_SE2}"
	
	unset R1
	unset R2
	unset ATROPOS_PE1
	unset ATROPOS_PE2
	unset PRINSEQ_PE1
	unset PRINSEQ_PE2
	unset PRINSEQ_SE1
	unset PRINSEQ_SE2
	unset CLEANED_PE1
	unset CLEANED_PE2
	unset CLEANED_SE1
	unset CLEANED_SE2
done
