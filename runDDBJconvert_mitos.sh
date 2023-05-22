#!/bin/sh

# DDBJ file converter for mitochondrial genomes ver. 2.0 (2023.04.14)

if [ $# != 4 ]; then
  echo "argument error!"
  echo "usage: <this> <general_information.tsv> <genome_file.fasta> <annotation.gff> <output_name>"
  exit 1
fi

scriptdir=`dirname ${0}`

#pre-processing for sequence and gff files
mkdir ${4}_DDBJfileMittemp
sed -e 's/\ .*//' $2 > ${4}_DDBJfileMittemp/renamedseq.fasta
seqkit fx2tab ${4}_DDBJfileMittemp/renamedseq.fasta > ${4}_DDBJfileMittemp/renamedseq.tab
seqkit locate -i -r -p '"N{10,}"' --gtf -P $2 > ${4}_DDBJfileMittemp/gaps.gtf
cat ${4}_DDBJfileMittemp/gaps.gtf $3 > ${4}_DDBJfileMittemp/annot_pre.gff

#run scripts
python3 ${scriptdir}/DDBJconverter_mitos.py $1 ${4}_DDBJfileMittemp/renamedseq.tab ${4}_DDBJfileMittemp/annot_pre.gff > ${4}.ann.txt




