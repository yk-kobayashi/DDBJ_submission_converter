#!/bin/sh

if [ $# != 3 ]; then
  echo "argument error!"
  echo "usage: <this> <general_information.tsv> <genome_file.fasta> <annotation.gff>"
  exit 1
fi

#pre-processing for sequence and gff files
mkdir DDBJfiletemp
sed -e 's/\ .*//' $2 > DDBJfiletemp/renamedseq.fasta
seqkit fx2tab DDBJfiletemp/renamedseq.fasta > DDBJfiletemp/renamedseq.tab
seqkit locate -i -r -p '"N{10,}"' --gtf -P $2 > DDBJfiletemp/gaps.gtf
cat DDBJfiletemp/gaps.gtf $3 > DDBJfiletemp/annot_pre.gff

#main scripts
python3 DDBJconverter_mitos.py $1 DDBJfiletemp/renamedseq.tab DDBJfiletemp/annot_pre.gff > DDBJannotfile.txt
python3 outputDDBJseqfile.py DDBJfiletemp/renamedseq.tab > DDBJseqfile.txt




