xlsx2csv --all -i -d tab *.xlsx ./

paste -d '\t' <(cut -f 2 SNV.csv) <(cut -f 3 SNV.csv) <(cut -f 4 SNV.csv | cut -d '>' -f 1) <(cut -f 4 SNV.csv | cut -d '>' -f 2) <(cut -f 1 SNV.csv) <(cut -f 5- SNV.csv) | sed -e '1d' > annotation.tsv
paste -d '\t' <(cut -f 2 Indels.csv) <(cut -f 3 Indels.csv) <(cut -f 4 Indels.csv | cut -d '>' -f 1) <(cut -f 4 Indels.csv | cut -d '>' -f 2) <(cut -f 1 Indels.csv) <(cut -f 5,6,10,11,12 Indels.csv) <(cut -f 8,13,15 Indels.csv) | sed -e '1d' > annotation2.tsv

cat annotation.tsv annotation2.tsv | awk -F"\t" -v "OFS=\t" '{print "chr"$1,$2,".",$3,$4,".",".",".","GT:DP", "0|1:100"}'  > annotation3.tsv
cat /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/vcf_head2.tsv annotation3.tsv > annotation.vcf

cut -f 10 annotation.tsv annotation2.tsv | cut -d '.' -f 1 > gene_list.txt

R CMD BATCH ../script/Biomart.R

tr -d \\r <annotation_NM.tsv> annotation_NM.tr.tsv
sed -e '1d' annotation_NM.tr.tsv > annotation_NM.tr2.tsv
cut -f 4 annotation_NM.tr2.tsv > ENST_list.txt

cat annotation.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > annotation_sort.vcf
bgzip annotation_sort.vcf
tabix -p vcf annotation_sort.vcf.gz
bcftools norm -f /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg19.fa -o annotation_sort.leftalign.vcf annotation_sort.vcf.gz
awk -F"\t" -v "OFS=\t" '{gsub("c","C",$4);gsub("g","G",$4);gsub("t","T",$4);gsub("a","A",$4);print $0}' annotation_sort.leftalign.vcf > annotation_sort.leftalign2.vcf


java -jar /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/snpEff/snpEff.jar -onlyTr ENST_list.txt GRCh37.87 annotation_sort.leftalign2.vcf > annotation.snpeff.vcf
cat annotation.snpeff.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > annotation.snpeff.sort.vcf

java -Xms10g -Xmx30g -jar ../snpEff/SnpSift.jar annotate ../database/clinvar_*.hg19.recode2.vcf.gz annotation.snpeff.sort.vcf | 
java -jar ../snpEff/SnpSift.jar annotate ../database/tommo-54kjpn-GRCh37-af2.sort.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/tommo-38kjpn-GRCh37-af2.sort.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/tommo-14kjpn-GRCh37-af2.sort.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/tommo-8.3kjpn-20200831-af2.sort.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.exomes.r2.1.1.sites.recode2.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.genomes.r2.1.1.sites.all.recode2.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.genomes.v3.1.2.sites.all.recode2.hg19.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.exomes.v4.1.sites.all.recode2.hg19.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.genomes.v4.1.sites.all.recode2.hg19.norm.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate ../database/BBJ_RIKEN_TMM_20200309_all2.sort.norm.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate ../database/GCF_000001405.25.recode2.norm2.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/CosmicCodingMuts.normal.v99.recode2.variants.uniq.sort.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/CosmicNonCodingVariants.normal.v99.recode2.variants.uniq.sort.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate ../database/all.frequency.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/MGeND.recode2.norm2.vcf.gz > annotation.snpeff.sort.annotate.vcf
java -Xms10g -Xmx30g -jar ../snpEff/SnpSift.jar extractFields annotation.snpeff.sort.annotate.vcf ANN[0].GENE ANN[0].IMPACT ANN[0].EFFECT ANN[0].HGVS_C ANN[0].HGVS_P > annotation.snpeff.sort.annotate.extraFields.tsv

R CMD BATCH ../script/vcfr_Genmine.R
tr -d \\r <annotation.snpeff.sort.annotate.tsv> annotation.snpeff.sort.annotate.tr.tsv
paste -d "\t" annotation.snpeff.sort.annotate.tr.tsv annotation.snpeff.sort.annotate.extraFields.tsv | tr "," "/" | tr "\t" "," > annotation.snpeff.sort.annotate.csv

paste -d "," <(awk -F"," -v "OFS=," '{print $2":"$3"_"$5"_"$6,$4}' annotation.snpeff.sort.annotate.csv) <(cut -d , -f 7- annotation.snpeff.sort.annotate.csv) > annotation.snpeff.sort.annotate2.csv
cat <(paste -d '\t' <(cut -f 2 SNV.csv) <(cut -f 3 SNV.csv) <(cut -f 4 SNV.csv | cut -d '>' -f 1) <(cut -f 4 SNV.csv | cut -d '>' -f 2) <(cut -f 1 SNV.csv) <(cut -f 5- SNV.csv) | head -n 1) annotation.tsv annotation2.tsv > annotation4.tsv
awk -F"\t" -v "OFS=," '{print "chr"$1":"$2"_"$3"_"$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' annotation4.tsv | tr "\t" "," > annotation5.csv
join --header -t "," -a 1 -1 1 -2 1 -e "." -o auto <(head -n +1 annotation5.csv && tail -n +2 annotation5.csv | sort -k 1,1 -t ',') <(head -n +1 annotation.snpeff.sort.annotate2.csv && tail -n +2 annotation.snpeff.sort.annotate2.csv | sort -k 1,1 -t ',') > annotation.snpeff.sort.annotate3.csv

awk -F, 'NR == 1 || $6 == 1' annotation.snpeff.sort.annotate3.csv > annotation.snpeff.sort.annotate4.csv

rm annotation.snpeff.sort.annotate.csv annotation.snpeff.sort.annotate.extraFields.tsv annotation.snpeff.sort.annotate.tr.tsv annotation.snpeff.sort.annotate.tsv annotation.snpeff.sort.annotate.vcf annotation.snpeff.sort.annotate2.csv annotation.snpeff.sort.vcf annotation.snpeff.vcf annotation.tsv annotation.vcf annotation2.tsv annotation3.tsv annotation4.tsv annotation5.csv annotation_NM.tr.tsv annotation_NM.tr2.tsv annotation_sort.leftalign.vcf annotation_sort.leftalign2.vcf annotation_sort.vcf.gz annotation_sort.vcf.gz.tbi snpEff_genes.txt 