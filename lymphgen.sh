uniq shortVariants3.snpeff.sort.annotate.vcf > shortVariants3.snpeff.sort.annotate2.vcf

java -Xms10g -Xmx40g -jar /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/picard.jar LiftoverVcf I=shortVariants3.snpeff.sort.annotate2.vcf O=shortVariants3.snpeff.sort.annotate2.hg19.vcf CHAIN=/mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg38ToHg19.over.chain REJECT=annotation.reject.vcf R=/mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg19.fa WRITE_ORIGINAL_POSITION=true

java -Xms10g -Xmx30g -jar ../snpEff/SnpSift.jar extractFields shortVariants3.snpeff.sort.annotate2.hg19.vcf ANN[0].GENE ANN[0].IMPACT ANN[0].EFFECT ANN[0].HGVS_C ANN[0].HGVS_P > shortVariants3.snpeff.sort.annotate2.hg19.extraFields.tsv

R CMD BATCH ../script/vcfr_heme_lymphgen.R

tr -d \\r <shortVariants3.snpeff.sort.annotate2.hg19.tsv> shortVariants3.snpeff.sort.annotate2.hg19.tr.tsv
paste -d "\t" shortVariants3.snpeff.sort.annotate2.hg19.tr.tsv shortVariants3.snpeff.sort.annotate2.hg19.extraFields.tsv | tr "," "/" | tr "\t" "," > shortVariants3.snpeff.sort.annotate2.hg19.csv

awk -F"," -v "OFS=," '{print $38":"$39"_"$40"_"$37,$0}' shortVariants3.snpeff.sort.annotate2.hg19.csv > shortVariants3.snpeff.sort.annotate2.hg19_2.csv
join -t "," --header -a 1 -1 1 -2 1 -e "." -o auto <(head -n +1 shortVariants.csv && tail -n +2 shortVariants.csv | sort -k 1,1 -t ',') <(head -n +1 shortVariants3.snpeff.sort.annotate2.hg19_2.csv && tail -n +2 shortVariants3.snpeff.sort.annotate2.hg19_2.csv | sort -k 1,1 -t ',' | uniq ) > shortVariants3.snpeff.sort.annotate2.hg19_3.csv

R CMD BATCH ../script/lymphgen.R

