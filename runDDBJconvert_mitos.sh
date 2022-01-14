#!/bin/sh

if [ $# != 3 ]; then
  echo "argument error!"
  echo "usage: <this> <general_information.tsv> <genome_file.fasta> <annotation.gff>"
  exit 1
fi

#pre-processing for sequence and gff files
mkdir DDBJfileMittemp
sed -e 's/\ .*//' $2 > DDBJfileMittemp/renamedseq.fasta
seqkit fx2tab DDBJfileMittemp/renamedseq.fasta > DDBJfileMittemp/renamedseq.tab
seqkit locate -i -r -p '"N{10,}"' --gtf -P $2 > DDBJfileMittemp/gaps.gtf
cat DDBJfileMittemp/gaps.gtf $3 > DDBJfileMittemp/annot_pre.gff

#main scripts
python3 DDBJconverter_mitos.py $1 DDBJfileMittemp/renamedseq.tab DDBJfileMittemp/annot_pre.gff > DDBJannotfile_mit.txt
python3 outputDDBJseqfile.py $1 DDBJfileMittemp/renamedseq.tab > DDBJseqfile_mit.txt




