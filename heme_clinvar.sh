

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




# 入力 JSON ファイルと出力 TSV ファイルのパスを指定
json_file=$(ls *.json)
tsv_file="shortVariants.tsv"

# ヘッダー行を作成（必要な列名を適宜指定）
header="itemId\tchromosome\tposition\treferenceAllele\talternateAllele\talternateAlleleFrequency\ttotalReadDepth\talternateAlleleReadDepth\ttranscriptId\tgeneSymbol\tcdsChange\taminoAcidsChange\tcalculatedEffects\tgeneID\tnumber_exonIntron\tnumber_number\tnumber_numberTotal\tvalidated\tdetectedInPaired\tcosmicCountThreshold\tcytoband\tdatabase_cosmicAll\tdatabase_cosmicHeme\tdatabase_clinVar\tdatabase_tommo\tdatabase_gnomAD\tdatabase_hgvd\tdatabase_igsr1000Genomes\tdatabase_dbSNP\tmnvReadNumber_refRef\tmnvReadNumber_altAlt\tmnvReadNumber_refAlt\tmnvReadNumber_altRef\tmnvReadNumber_refRefRef\tmnvReadNumber_refRefAlt\tmnvReadNumber_refAltRef\tmnvReadNumber_refAltAlt\tmnvReadNumber_altRefRef\tmnvReadNumber_altRefAlt\tmnvReadNumber_altAltRef\tmnvReadNumber_altAltAlt\tfunction_fastTrack\tfunction_mitelman"

# ヘッダー行を出力ファイルに書き出す
echo -e "$header" > "$tsv_file"

# jq で各 variant をフラットに展開して TSV 形式に変換
jq -r '
  .variants.shortVariants[] |
  # transcripts 配列の最初の要素（存在しない場合は {}）を$t として取得
  (if (.transcripts | length) > 0 then .transcripts[0] else {} end) as $t |
  # transcripts 内の number オブジェクトを $num として取得（存在しなければ {}）
  (if ($t.number | type) == "object" then $t.number else {} end) as $num |
  [
    .itemId,
    .chromosome,
    .position,
    .referenceAllele,
    .alternateAllele,
    .alternateAlleleFrequency,
    .totalReadDepth,
    .alternateAlleleReadDepth,
    $t.transcriptId,
    $t.geneSymbol,
    $t.cdsChange,
    $t.aminoAcidsChange,
    ($t.calculatedEffects | join(",")),
    ($t.geneID | join(",")),
    $num.exonIntron,
    $num.number,
    $num.numberTotal,
    .validated,
    .detectedInPaired,
    .cosmicCountThreshold,
    .cytoband,
    .database.cosmicAll,
    .database.cosmicHeme,
    (.database.clinVar | join(",")),
    .database.tommo,
    .database.gnomAD,
    .database.hgvd,
    .database.igsr1000Genomes,
    (.database.dbSNP | join(",")),
    .mnvReadNumber.refRef,
    .mnvReadNumber.altAlt,
    .mnvReadNumber.refAlt,
    .mnvReadNumber.altRef,
    .mnvReadNumber.refRefRef,
    .mnvReadNumber.refRefAlt,
    .mnvReadNumber.refAltRef,
    .mnvReadNumber.refAltAlt,
    .mnvReadNumber.altRefRef,
    .mnvReadNumber.altRefAlt,
    .mnvReadNumber.altAltRef,
    .mnvReadNumber.altAltAlt,
    .function.fastTrack,
    (.function.mitelman | join(","))
  ] | @tsv
' "$json_file" >> "$tsv_file"


cat /mnt/c/database/vcf_head3.tsv <(sed -e '1d' shortVariants.tsv | awk -F"\t" -v "OFS=\t" '{print "chr"$2,$3,".",$4,$5,".",".",".","GT:DP", "0|1:100"}') > shortVariants.vcf

cat shortVariants.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | uniq > shortVariants.sort.vcf

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


awk -F"\t" -v "OFS=\t" '{print $9}' shortVariants.tsv | sed -e '1d' | sort | uniq > shortVariants.refseq.tsv

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

awk -F"\t" -v "OFS=\t" '{print "chr"$2":"$3"_"$4"_"$5,$0}' shortVariants.tsv | tr "," "/" | tr "\t" "," > shortVariants.csv

awk -F"," -v "OFS=," '{print $10":"$11"_"$12"_"$13,$0}' shortVariants3.snpeff.sort.annotate.csv > shortVariants3.snpeff.sort.annotate2.csv
join -t "," --header -a 1 -1 1 -2 1 -e "." -o auto <(head -n +1 shortVariants.csv && tail -n +2 shortVariants.csv | sort -k 1,1 -t ',') <(head -n +1 shortVariants3.snpeff.sort.annotate2.csv && tail -n +2 shortVariants3.snpeff.sort.annotate2.csv | sort -k 1,1 -t ',' | uniq ) > shortVariants2.csv


join -t "," --header -a 1 -1 11 -2 1 -e "." -o auto <(head -n +1 shortVariants2.csv && tail -n +2 shortVariants2.csv | sort -k 11,11 -t ',') /mnt/c/database/cancerGeneList3.csv > shortVariants3.csv

(head -n +1 shortVariants3.csv && tail -n +2 shortVariants3.csv | sort -k 3,3 -t ',') > shortVariants4.csv





