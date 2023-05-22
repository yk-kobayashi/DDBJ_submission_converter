#!/bin/sh

# DDBJ file converter for eukaryotic genomes ver. 2.1 (2023.04.14)

if [ $# != 5 ]; then
  echo "argument error!"
  echo "usage: <this> <general_information.tsv> <genome_file.fasta> <annotation.gtf> <functional_annotation.csv> <output_name>"
  exit 1
fi

scriptdir=`dirname ${0}`

#pre-processing for sequence and gff files
mkdir ${5}_DDBJfiletemp
sed -e 's/\ .*//' $2 | seqkit fx2tab > ${5}_DDBJfiletemp/renamedseq.tab
seqkit locate -i -r -p '"N{10,}"' --gtf -P $2 > ${5}_DDBJfiletemp/gaps.gtf
grep -v ^"#" $3 | grep -e "\tgene\t" -e "\ttranscript\t" -e "\tstart_codon\t" -e "\tCDS\t" -e "\tstop_codon\t" -e "\ttRNA\t" -e "\trRNA\t"  > ${5}_DDBJfiletemp/annot_simple.gtf
cat ${5}_DDBJfiletemp/gaps.gtf ${5}_DDBJfiletemp/annot_simple.gtf > ${5}_DDBJfiletemp/annot_pre.gtf

#main scripts
python3 ${scriptdir}/DDBJconverter_eukaryote.py $1 ${5}_DDBJfiletemp/renamedseq.tab ${5}_DDBJfiletemp/annot_pre.gtf $4 > ${5}.ann.txt
