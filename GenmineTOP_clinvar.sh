
# ファイル名の抽出と処理
for file in /mnt/c/database/clinvar_*.hg38.recode2.vcf.gz; do
  # * の部分（数字列）を抽出
  date=$(basename $file | sed -E 's/clinvar_(.*)\.hg38.recode2\.vcf\.gz/\1/')

  # ファイルが存在するか確認
  if ! curl --head --silent --fail "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_${date}.vcf.gz" > /dev/null; then
    echo "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_${date}.vcf.gz が存在しません。ダウンロードを実行します。"
    # ダウンロードを実行
    latest_file=$(curl -s https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/ | grep -oP 'clinvar_\d+\.vcf\.gz' | sort -r | head -n 1)
    wget -c "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/$latest_file" -P /mnt/c/database
    rm /mnt/c/database/clinvar_${date}.hg38.recode2.vcf.gz /mnt/c/database/clinvar_${date}.hg38.recode2.vcf.gz.tbi
    
    filename_base=$(echo "$latest_file" | sed 's/\.vcf\.gz//')
    vcftools --gzvcf /mnt/c/database/$latest_file --recode --out /mnt/c/database/$filename_base.hg38 --recode-INFO AF_ESP --recode-INFO AF_EXAC --recode-INFO AF_TGP --recode-INFO CLNSIG --recode-INFO CLNSIGCONF --recode-INFO SCI --recode-INFO ONC
    sed '9,15d' /mnt/c/database/$filename_base.hg38.recode.vcf | sed '11,22d' | sed '12,22d' | sed '13,14d' > /mnt/c/database/$filename_base.hg38.recode2.vcf
    bgzip -@ 24 /mnt/c/database/$filename_base.hg38.recode2.vcf
    tabix -p vcf /mnt/c/database/$filename_base.hg38.recode2.vcf.gz
    rm /mnt/c/database/$filename_base.vcf.gz /mnt/c/database/$filename_base.hg38.recode.vcf /mnt/c/database/$filename_base.hg38.log
    
  else
    echo "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_${date}.vcf.gz は既に存在しています。"
  fi
done


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

cat /mnt/c/database/vcf_head3.tsv Book7.tsv > shortVariants.vcf

cat shortVariants.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > shortVariants.sort.vcf

awk 'BEGIN {OFS="\t"} 
     # ヘッダー行はそのまま出力
     /^#/ {print; next} 
     {
         # 元の情報をまとめる
         orig_info = "ORIG_CHROM=" $1 ";ORIG_POS=" $2 ";ORIG_REF=" $4 ";ORIG_ALT=" $5;
         # INFOフィールドが"."の場合とそうでない場合で追加の仕方を分ける
         if($8 == ".") {
             $8 = orig_info;
         } else {
             $8 = $8 ";" orig_info;
         }
         print
     }' shortVariants.sort.vcf > shortVariants.sort2.vcf

bgzip shortVariants.sort2.vcf
tabix -p vcf shortVariants.sort2.vcf.gz
bcftools norm -f /mnt/c/database/hg38.fa -o shortVariants2.vcf shortVariants.sort2.vcf.gz -c s
awk -F"\t" -v "OFS=\t" '{gsub("c","C",$4);gsub("g","G",$4);gsub("t","T",$4);gsub("a","A",$4);gsub("c","C",$5);gsub("g","G",$5);gsub("t","T",$5);gsub("a","A",$5);print $0}' shortVariants2.vcf > shortVariants3.vcf
bgzip -d shortVariants.sort2.vcf.gz

cut -f 5 annotation_NM.tr2.tsv | sort | uniq > shortVariants.refseq.tsv

java -jar ../snpEff/snpEff.jar -onlyTr shortVariants.refseq.tsv GRCh38.p13.RefSeq shortVariants3.vcf > shortVariants3.snpeff.vcf

cat shortVariants3.snpeff.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > shortVariants3.snpeff.sort.vcf

java -Xms10g -Xmx30g -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/clinvar_*.hg38.recode2.vcf.gz shortVariants3.snpeff.sort.vcf | 
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/tommo-60kjpn-GRCh38-af.norm.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/tommo-54kjpn-GRCh38-af.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/tommo-38kjpn-GRCh38-af.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/tommo-14kjpn-GRCh38-af.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/gnomad.exomes.v4.1.sites.all.hg38.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/gnomad.genomes.v4.1.sites.all.hg38.norm.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/gem_j_wga.all.grch38.recode3.norm.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/jga_snp.grch38.recode3.norm.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/jga_wes.grch38.recode3.norm.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/jga_wgs.recode2.norm.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/ncbn.jpn.all.grch38.recode3.norm.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/GCF_000001405.40.recode2.norm3.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/Cosmic_CompleteTargetedScreensMutant_Normal_v101_GRCh38.recode2.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/Cosmic_GenomeScreensMutant_Normal_v101_GRCh38.recode2.norm.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate /mnt/c/database/Cosmic_NonCodingVariants_Normal_v101_GRCh38.recode2.norm.vcf.gz > shortVariants3.snpeff.sort.annotate.vcf
java -Xms10g -Xmx30g -jar ../snpEff/SnpSift.jar extractFields shortVariants3.snpeff.sort.annotate.vcf ANN[0].GENE ANN[0].IMPACT ANN[0].EFFECT ANN[0].HGVS_C ANN[0].HGVS_P > shortVariants3.snpeff.sort.annotate.extraFields.tsv

R CMD BATCH ../script/vcfr_heme.R

tr -d \\r <shortVariants3.snpeff.sort.annotate.tsv> shortVariants3.snpeff.sort.annotate.tr.tsv
paste -d "\t" shortVariants3.snpeff.sort.annotate.tr.tsv shortVariants3.snpeff.sort.annotate.extraFields.tsv | tr "," "/" | tr "\t" "," > shortVariants3.snpeff.sort.annotate.csv

cat /mnt/c/database/GeneMine.head.csv <(awk -F"\t" -v "OFS=," '{print $6"_"$5"_"$4,$2,$3,$7,$8,$9,$10,$11,$13,$1,$14}' Book5.tsv) > Book8.csv
awk -F"," -v "OFS=," '{print $9":"$10"_"$11"_"$12,$0}' shortVariants3.snpeff.sort.annotate.csv > shortVariants3.snpeff.sort.annotate2.csv
join --header -t "," -a 1 -1 1 -2 1 -e "." -o auto <(head -n +1 Book8.csv && tail -n +2 Book8.csv | sort -k 1,1 -t ',') <(head -n +1 shortVariants3.snpeff.sort.annotate2.csv && tail -n +2 shortVariants3.snpeff.sort.annotate2.csv | sort -k 1,1 -t ',') > annotation.snpeff.sort.annotate3.csv

join -t "," --header -a 1 -1 4 -2 1 -e "." -o auto <(head -n +1 annotation.snpeff.sort.annotate3.csv && tail -n +2 annotation.snpeff.sort.annotate3.csv | sort -k 4,4 -t ',') /mnt/c/database/cancerGeneList3.csv > annotation.snpeff.sort.annotate4.csv


# awk '/<values>/,/<\/values>/{if(/<item>/)print}' *.xml | sed -e 's/<item>//g' -e 's/<\/item>//g' | cut -f 7 > mutation_list.txt
# R CMD BATCH ../script/signature.R


rm Book.tsv Book2.tsv Book3.tsv Book5.tsv Book6.tsv Book7.tsv Book8.csv annotation_NM.tr.tsv annotation_NM.tr2.tsv annotation_NM.tsv shortVariants.sort.vcf shortVariants.sort2.vcf shortVariants.sort2.vcf.gz.tbi shortVariants.vcf shortVariants2.vcf shortVariants3.snpeff.sort.annotate.csv shortVariants3.snpeff.sort.annotate.extraFields.tsv shortVariants3.snpeff.sort.annotate.tr.tsv shortVariants3.snpeff.sort.annotate.tsv shortVariants3.snpeff.sort.annotate.vcf shortVariants3.snpeff.sort.annotate2.csv shortVariants3.snpeff.sort.vcf shortVariants3.snpeff.vcf shortVariants3.vcf vcfr_heme.Rout gene_list.txt annotation.snpeff.sort.annotate3.csv

