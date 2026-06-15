cat *.xml | cut -f 6 | grep -v "<ag-class>" > Book.tsv

awk -F'\t' '{
  lines[NR % 12] = $0;
}
$1 ~ /<type>snv<\/type>/ || /<type>insertion<\/type>/ || /<type>deletion<\/type>/ || /<type>deletion-insertion<\/type>/ {
  buffer = ""
  for (i = NR; i >= NR-7; i--) {
    buffer = buffer lines[i % 12] "\t";
  }
  found = 4;
  output = buffer;
}
found && found-- {
  output = output $0 "\t";
  if (found == 0) {
    print output;
  }
}' Book.tsv | sed 's/<[^>]*>//g' | sed -e "s/&gt;/>/g" > Book2.tsv

cut -f 7 Book2.tsv > gene_list.txt
cut -f 2-12 Book2.tsv > Book3.tsv

R CMD BATCH ../script/Biomart.R

tr -d \\r <annotation_NM.tsv> annotation_NM.tr.tsv
sed -e '1d' annotation_NM.tr.tsv > annotation_NM.tr2.tsv
join -t "$(printf '\011')" -a 1 -1 6 -2 2 -e "." -o auto <(cat Book3.tsv | sort -k 6,6 -t$'\t') <(cat annotation_NM.tr2.tsv | sort -k 2,2 -t$'\t') > Book5.tsv
awk -F"\t" -v "OFS=\t" '{gsub(/[0-9]/,"",$9);gsub(/-/,"",$9);gsub(/+/,"",$9);gsub(/_/,"",$9);gsub(/c./,"",$9);print $0}' Book5.tsv > Book6.tsv
paste -d '\t' <(cut -f 6 Book6.tsv | cut -d : -f 1) <(cut -f 6 Book6.tsv | cut -d : -f 2) <(cut -f 4,5 Book6.tsv) | awk -F"\t" -v "OFS=\t" '{print $1,$2,".",$4,$3,".",".",".","GT:DP", "0|1:100"}' > Book7.tsv
cat /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/vcf_head2.tsv Book7.tsv > annotation.vcf
java -Xms10g -Xmx40g -jar /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/picard.jar LiftoverVcf I=annotation.vcf O=annotation.hg19.vcf CHAIN=/mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg38ToHg19.over.chain REJECT=annotation.reject.vcf R=/mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg19.fa WRITE_ORIGINAL_POSITION=true

bgzip annotation.hg19.vcf
tabix -p vcf annotation.hg19.vcf.gz
bcftools norm -f /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg19.fa -o annotation.hg19.leftalign.vcf annotation.hg19.vcf.gz
awk -F"\t" -v "OFS=\t" '{gsub("c","C",$4);gsub("g","G",$4);gsub("t","T",$4);gsub("a","A",$4);print $0}' annotation.hg19.leftalign.vcf > annotation.hg19.leftalign2.vcf

cut -f 4 annotation_NM.tr2.tsv > ENST_list.txt

# 対応するENST numberのないものがあるか
blank_number=`awk -F"\t" -v "OFS=\t" '$1=="" {print $0}' ENST_list.txt | wc -l`
if [[ $blank_number >0 ]] ; then
	echo "There are some blanks in ENST_list.txt"
fi

java -jar /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/snpEff/snpEff.jar -onlyTr ENST_list.txt GRCh37.87 annotation.hg19.leftalign2.vcf > annotation.snpeff.vcf
cat annotation.snpeff.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > annotation.snpeff.sort.vcf

java -Xms10g -Xmx30g -jar ../snpEff/SnpSift.jar annotate ../database/clinvar_*.hg19.recode2.vcf.gz annotation.snpeff.sort.vcf | 
java -jar ../snpEff/SnpSift.jar annotate ../database/tommo-54kjpn-GRCh37-af2.sort.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/tommo-38kjpn-GRCh37-af2.sort.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/tommo-14kjpn-GRCh37-af2.sort.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/tommo-8.3kjpn-20200831-af2.sort.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.exomes.r2.1.1.sites.recode2.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.genomes.r2.1.1.sites.all.recode2.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.genomes.v3.1.2.sites.all.recode2.hg19.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.exomes.v4.0.sites.all.recode2.hg19.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.genomes.v4.0.sites.all.recode2.hg19.norm.vcf.gz |
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

paste -d "," <(awk -F"," -v "OFS=," '{print $10":"$11"_"$5"_"$6,$2,$3,$4,$9}' annotation.snpeff.sort.annotate.csv) <(cut -d , -f 14- annotation.snpeff.sort.annotate.csv) > annotation.snpeff.sort.annotate2.csv
cat ../database/GeneMine.head.csv <(awk -F"\t" -v "OFS=," '{print $6"_"$5"_"$4,$2,$3,$7,$8,$9,$10,$11,$13,$1,$14}' Book5.tsv) > Book8.csv
join --header -t "," -a 1 -1 1 -2 1 -e "." -o auto <(head -n +1 Book8.csv && tail -n +2 Book8.csv | sort -k 1,1 -t ',') <(head -n +1 annotation.snpeff.sort.annotate2.csv && tail -n +2 annotation.snpeff.sort.annotate2.csv | sort -k 1,1 -t ',') > annotation.snpeff.sort.annotate3.csv


awk '/<values>/,/<\/values>/{if(/<item>/)print}' *.xml | sed -e 's/<item>//g' -e 's/<\/item>//g' | cut -f 7 > mutation_list.txt
R CMD BATCH ../script/signature.R


rm Book.tsv Book2.tsv Book3.tsv Book5.tsv Book6.tsv Book7.tsv Book8.csv annotation.hg19.leftalign.vcf annotation.hg19.leftalign2.vcf annotation.hg19.vcf.gz annotation.hg19.vcf.gz.tbi annotation.hg19.vcf.idx annotation.reject.vcf annotation.snpeff.sort.annotate.csv annotation.snpeff.sort.annotate.extraFields.tsv annotation.snpeff.sort.annotate.tr.tsv annotation.snpeff.sort.annotate.tsv annotation.snpeff.sort.annotate.vcf annotation.snpeff.sort.annotate2.csv annotation.snpeff.sort.vcf annotation.snpeff.vcf annotation.vcf annotation_NM.tr.tsv annotation_NM.tsv gene_list.txt mutation_list.txt snpEff_summary.html vcfr_Genmine.Rout
