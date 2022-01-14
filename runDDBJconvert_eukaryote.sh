#!/bin/sh

if [ $# != 4 ]; then
  echo "argument error!"
  echo "usage: <this> <general_information.tsv> <genome_file.fasta> <annotation.gff> <functional_annotation.csv>"
  exit 1
fi

#pre-processing for sequence and gff files
mkdir DDBJfiletemp
sed -e 's/\ .*//' $2 > DDBJfiletemp/renamedseq.fasta
seqkit fx2tab DDBJfiletemp/renamedseq.fasta > DDBJfiletemp/renamedseq.tab
seqkit locate -i -r -p '"N{10,}"' --gtf -P $2 > DDBJfiletemp/gaps.gtf
grep -v ^"#" $3 | grep -e "\tgene\t" -e "\ttranscript\t" -e "\tstart_codon\t" -e "\tCDS\t" -e "\tstop_codon\t" -e "\ttRNA\t" -e "\trRNA\t"  > DDBJfiletemp/annot_simple.gff
cat DDBJfiletemp/gaps.gtf DDBJfiletemp/annot_simple.gff > DDBJfiletemp/annot_pre.gff

#main scripts
python3 DDBJconverter_eukaryote.py $1 DDBJfiletemp/renamedseq.tab DDBJfiletemp/annot_pre.gff $4 > DDBJannotfile.txt
python3 outputDDBJseqfile.py $1 DDBJfiletemp/renamedseq.tab > DDBJseqfile.txt
