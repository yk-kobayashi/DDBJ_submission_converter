#!/bin/sh

# DDBJ file converter for prokaryotic genomes ver. 2.0 (2023.04.14)

if [ $# != 5 ]; then
  echo "argument error!"
  echo "usage: <this> <general_information.tsv> <genome_file.fasta> <annotation.gtf> <functional_annotation.csv> <output_name>"
  exit 1
fi

scriptdir=`dirname ${0}`

#pre-processing for sequence and gff files
mkdir ${5}_DDBJfiletemp
sed -e 's/\ .*//' $2 > ${5}_DDBJfiletemp/renamedseq.fasta
seqkit fx2tab ${5}_DDBJfiletemp/renamedseq.fasta > ${5}_DDBJfiletemp/renamedseq.tab
seqkit locate -i -r -p '"N{10,}"' --gtf -P $2 > ${5}_DDBJfiletemp/gaps.gtf
cat ${5}_DDBJfiletemp/gaps.gtf $3 > ${5}_DDBJfiletemp/annot_pre.gff

#main scripts
python3 ${scriptdir}/DDBJconverter_prokaryote.py $1 ${5}_DDBJfiletemp/renamedseq.tab ${5}_DDBJfiletemp/annot_pre.gff $4 > ${5}.ann.txt
#python3 outputDDBJseqfile.py $1 ${2}_DDBJfiletemp/renamedseq.tab > ${2}_DDBJseqfile.txt
