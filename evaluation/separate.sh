#!/bin/bash

infile=$1 # fasta file of annotated genes
name=${infile:0:4}
declare -A arr
arr['Mmus']=10090
arr['Rnor']=10116
arr['Ppan']=9598
arr['Ptro']=9598
arr['Cfam']=9615

echo "GeneID,Alignment %ID,Query Length, Query Aligned Length, Target Length, Target Aligned Length" > ${name}_metrics.csv
while read gene;
do
	echo $gene
        echo $name
	gawk -i inplace '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $infile # unwrap fasta file
	grep -i -A1 ">"$gene $infile > ${name}_${gene}.temp # find gene entry
        awk -v gene=$gene -v val=${arr[${name}]} 'BEGIN{IGNORECASE=1} $1==val && $3==gene {print $5}' ../gene_info | tr "|" "\n" | sed 's/^/>/' > ${gene}_synonyms.txt # get synonyms for that gene name
        grep -A1 -iFf ${gene}_synonyms.txt $infile >> ${name}_${gene}.temp # add those entries
	header=`head -1 ${name}_${gene}.temp`
	echo ${header##*/} | cut -d_ -f1-2 > ${name}_${gene}_T.fa
#	cat ${name}_${gene}.temp | awk '{print length, $0}' | sort -nr | head -1 | awk '{print $2}' >> ${name}_${gene}_T.fa # gets longest annotated exon
	cat ${name}_${gene}.temp > ${name}_${gene}_T.fa # gets longest annotated exon
	grep -A1 $gene ${name}_extended.fa > ${name}_${gene}.temp
#	header1=`head -1 ${name}_${gene}.temp`
#	echo '>'${header1##*/} | cut -d_ -f1-2 >> ${name}_${gene}.fa
#	cat ${name}_${gene}.temp | awk '{print length, $0}' | sort -nr | head -1 | awk '{print $2}'>> ${name}_${gene}.fa
	cat ${name}_${gene}.temp > ${name}_${gene}_Q.fa
#	cat ${name}_${gene}.temp | while read L; do if [[ $L =~ ^'>' ]]; then echo $L; else echo $L | rev | tr "ATGC" "TACG"; fi ; done >> ${name}_${gene}_Q.fa
		# https://www.biostars.org/p/189325/#277090
#	cat ${name}_${gene}.temp >> ${name}_${gene}_Q.fa
	sed -i '/\-\-/d' ${name}_${gene}_T.fa
	sed -i '/\-\-/d' ${name}_${gene}_Q.fa
	exonerate -Q dna -T dna -m affine:local -q "${name}_${gene}_Q.fa" -t "${name}_${gene}_T.fa" -n 1 --ryo "%qid,%pi,%ql,%qal,%tl,%tal\n" --verbose 0 --showalignment no --showvulgar no > ${name}_${gene}_exonerate.out
#	echo ${gene}
	head -1 ${name}_${gene}_exonerate.out >> ${name}_metrics.csv
done < ../exon_lists/homo_1500_exons_succeeded_names.txt
