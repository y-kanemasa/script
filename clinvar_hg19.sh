# databaseで行う

#!/usr/bin/env bash
set -euo pipefail

shopt -s nullglob

files=(clinvar_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].vcf.gz)
old_files=(clinvar_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].hg19.recode2.vcf.gz)

if ((${#files[@]} == 1)); then
  latest_file="${files[0]}"
  filename_base=$(echo "$latest_file" | sed 's/\.vcf\.gz//')
  vcftools --gzvcf ../database/$latest_file --recode --out ../database/$filename_base.hg19 --recode-INFO AF_ESP --recode-INFO AF_EXAC --recode-INFO AF_TGP --recode-INFO CLNSIG --recode-INFO CLNSIGCONF --recode-INFO SCI --recode-INFO ONC
  sed '9,15d' ../database/$filename_base.hg19.recode.vcf | sed '11,22d' | sed '12,22d' | sed '13,14d' > ../database/$filename_base.hg19.recode2.vcf
  bgzip -@ 24 ../database/$filename_base.hg19.recode2.vcf
  tabix -p vcf ../database/$filename_base.hg19.recode2.vcf.gz
  rm ../database/$filename_base.vcf.gz ../database/$filename_base.hg19.recode.vcf ../database/$filename_base.hg19.log ../database/$old_files ../database/$old_files.tbi

else
  echo "clinvar_YYYYMMDD.vcf.gz はありません"
fi






