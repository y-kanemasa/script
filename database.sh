vcftools --gzvcf gnomad.genomes.v3.1.2.sites.chr${id}.vcf.bgz --recode --out gnomad.genomes.v3.1.2.sites.chr${id} --recode-INFO AF --recode-INFO AF_raw --recode-INFO AF_popmax --recode-INFO faf95_popmax
sed '7,8d' gnomad.genomes.v3.1.2.sites.chr${id}.recode.vcf | sed '8d' | sed '9,458d' | sed '10,373d' | sed '11,100d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_g3ge",$0);print $0}' | awk -F"\t" -v "OFS=\t" '{gsub("faf95","faf95_g3ge",$0);print $0}' > gnomad.genomes.v3.1.2.sites.chr${id}.recode2.vcf
java -Xms10g -Xmx40g -jar picard.jar LiftoverVcf I=gnomad.genomes.v3.1.2.sites.chr${id}.recode2.vcf O=gnomad.genomes.v3.1.2.sites.chr${id}.recode2.hg19.vcf CHAIN=hg38ToHg19.over.chain.gz REJECT=gnomad.genomes.v3.1.2.sites.chr${id}.recode2.reject.vcf R=./hg19.fa

java -Xms10g -Xmx40g -jar picard.jar MergeVcfs I=list.txt O=gnomad.genomes.v3.1.2.sites.all.recode2.hg19.vcf


vcftools --gzvcf gnomad.exomes.r2.1.1.sites.vcf.bgz --recode --out gnomad.exomes.r2.1.1.sites --recode-INFO AF --recode-INFO AF_raw --recode-INFO AF_popmax --recode-INFO faf95

sed '7,8d' gnomad.exomes.r2.1.1.sites.recode.vcf | sed '8,66d' | sed '9,736d' | sed '10,57d' | sed '11,39d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_g2ex",$0);print $0}' | awk -F"\t" -v "OFS=\t" '{gsub("faf95","faf95_g2ex",$0);print $0}' > gnomad.exomes.r2.1.1.sites.recode2.vcf


vcftools --gzvcf gnomad.genomes.r2.1.1.sites.${id}.vcf.bgz --recode --out gnomad.genomes.r2.1.1.sites.${id} --recode-INFO AF --recode-INFO AF_raw --recode-INFO AF_popmax --recode-INFO faf95
sed '7,8d' gnomad.genomes.r2.1.1.sites.${id}.recode.vcf | sed '8,58d' | sed '9,464d' | sed '10,46d' | sed '11,28d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_g2ge",$0);print $0}' | awk -F"\t" -v "OFS=\t" '{gsub("faf95","faf95_g2ge",$0);print $0}' > gnomad.genomes.r2.1.1.sites.${id}.recode2.vcf

java -Xms10g -Xmx40g -jar picard.jar MergeVcfs I=list.txt O=gnomad.genomes.r2.1.1.sites.all.recode2.vcf


vcftools --gzvcf gnomad.exomes.v4.1.sites.chr${id}.vcf.bgz --recode --out gnomad.exomes.v4.1.sites.chr${id} --recode-INFO AF --recode-INFO AF_raw --recode-INFO AF_grpmax
sed '7,8d' gnomad.exomes.v4.1.sites.chr${id}.recode.vcf | sed '8,228d' | sed '9,35d' | sed '10,169d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_g4ex",$0);print $0}' > gnomad.exomes.v4.1.sites.chr${id}.recode2.vcf
java -Xms10g -Xmx40g -jar ../picard.jar LiftoverVcf I=gnomad.exomes.v4.1.sites.chr${id}.recode2.vcf O=gnomad.exomes.v4.1.sites.chr${id}.recode2.hg19.vcf CHAIN=../hg38ToHg19.over.chain.gz REJECT=gnomad.exomes.v4.1.sites.chr${id}.recode2.reject.vcf R=../hg19.fa

ls *.recode2.hg19.vcf > list.txt
java -Xms10g -Xmx40g -jar ../picard.jar MergeVcfs I=list.txt O=gnomad.exomes.v4.1.sites.all.recode2.hg19.vcf

bcftools norm --check-ref w --threads 24 -f ../hg19.fa gnomad.exomes.v4.1.sites.all.recode2.hg19.vcf -o gnomad.exomes.v4.1.sites.all.recode2.hg19.norm.vcf


vcftools --gzvcf gnomad.genomes.v4.1.sites.chr${id}.vcf.bgz --recode --out gnomad.genomes.v4.1.sites.chr${id} --recode-INFO AF --recode-INFO AF_raw --recode-INFO AF_grpmax
sed '7,8d' gnomad.genomes.v4.1.sites.chr${id}.recode.vcf | sed '8,116d' | sed '9,35d' | sed '10,108d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_g4ge",$0);print $0}' > gnomad.genomes.v4.1.sites.chr${id}.recode2.vcf
java -Xms10g -Xmx40g -jar ../picard.jar LiftoverVcf I=gnomad.genomes.v4.1.sites.chr${id}.recode2.vcf O=gnomad.genomes.v4.1.sites.chr${id}.recode2.hg19.vcf CHAIN=../hg38ToHg19.over.chain.gz REJECT=gnomad.genomes.v4.1.sites.chr${id}.recode2.reject.vcf R=../hg19.fa

ls *.recode2.hg19.vcf > list.txt
java -Xms10g -Xmx40g -jar ../picard.jar MergeVcfs I=list.txt O=gnomad.genomes.v4.1.sites.all.recode2.hg19.vcf

bcftools norm --check-ref w --threads 24 -f ../hg19.fa gnomad.genomes.v4.1.sites.all.recode2.hg19.vcf -o gnomad.genomes.v4.1.sites.all.recode2.hg19.norm.vcf




vcftools --gzvcf tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.vcf.gz --recode --out tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome --recode-INFO AF




vcftools --gzvcf tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.vcf.gz --recode --out tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3 --recode-INFO AF






cat <(grep '^#' tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.vcf) <(grep -v '^#' tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.vcf | awk -F"\t" -v "OFS=\t" '{gsub("W","N",$4);gsub("Y","N",$4);gsub("H","N",$4);gsub("R","N",$4);gsub("S","N",$4);gsub("B","N",$4);gsub("M","N",$4);gsub("V","N",$4);gsub("K","N",$4);gsub("D","N",$4);gsub("W","N",$5);gsub("Y","N",$5);gsub("H","N",$5);gsub("R","N",$5);gsub("S","N",$5);gsub("B","N",$5);gsub("M","N",$5);gsub("V","N",$5);gsub("K","N",$5);gsub("D","N",$5);print $0}') > tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome2.vcf
java -Xms10g -Xmx40g -jar /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/picard.jar LiftoverVcf I=tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome2.vcf O=tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.hg19.vcf CHAIN=/mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg38ToHg19.over.chain REJECT=annotation.reject.vcf R=/mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg19.fa WRITE_ORIGINAL_POSITION=true
vcftools --vcf tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.hg19.vcf --recode --out tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.hg19 --recode-INFO AF
sed '21,40d' tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.hg19.recode.vcf | sed '21,109d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_tommo60",$0);print $0}' > tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.hg19.recode2.vcf



java -Xms10g -Xmx40g -jar /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/picard.jar LiftoverVcf I=tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.vcf O=tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.hg19.vcf CHAIN=/mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg38ToHg19.over.chain REJECT=annotation.reject.vcf R=/mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg19.fa WRITE_ORIGINAL_POSITION=true
vcftools --vcf tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.hg19.vcf --recode --out tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.hg19 --recode-INFO AF
sed '17,31d' tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.hg19.recode.vcf | sed '18,88d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_tommo60",$0);print $0}' > tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.hg19.recode2.vcf


java -Xms10g -Xmx40g -jar /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/picard.jar MergeVcfs I=tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.hg19.recode2.vcf I=tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.hg19.recode2.vcf O=tommo-60kjpn-hg19-af.vcf




vcftools --vcf tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.vcf --recode --out tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome --recode-INFO AF
sed '21,40d' tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.recode.vcf | sed '22,82d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_tommo60",$0);print $0}' > tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.recode2.vcf


vcftools --vcf tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.vcf --recode --out tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3 --recode-INFO AF
sed '21,40d' tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.recode.vcf | sed '22,82d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_tommo60",$0);print $0}' > tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.recode2.vcf

java -Xms10g -Xmx40g -jar /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/picard.jar SortVcf I=tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.recode2.vcf I=tommo-60kjpn-20240904-GRCh38-snvindel-af-chrX_PAR3.recode2.vcf O=tommo-60kjpn-af.vcf

java -Xms10g -Xmx40g -jar /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/picard.jar LiftoverVcf I=tommo-60kjpn-af.vcf O=tommo-60kjpn-af.hg19.vcf CHAIN=/mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg38ToHg19.over.chain REJECT=tommo-60kjpn-af.reject.vcf R=/mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg19.fa WRITE_ORIGINAL_POSITION=true



sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-autosome.vcf > tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-autosome.4.2.vcf
vcftools --vcf tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-autosome.4.2.vcf --recode --out tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-autosome.4.2 --recode-INFO AF
sed '20,21d' tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode.vcf | sed '21,109d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_tommo54",$0);print $0}' > tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode2.vcf
cat <(grep '^#' tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode2.vcf) <(grep -v '^#' tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode2.vcf | awk -F"\t" -v "OFS=\t" '{gsub("W","N",$4);gsub("Y","N",$4);gsub("H","N",$4);gsub("R","N",$4);gsub("S","N",$4);gsub("B","N",$4);gsub("M","N",$4);gsub("V","N",$4);gsub("K","N",$4);gsub("D","N",$4);gsub("W","N",$5);gsub("Y","N",$5);gsub("H","N",$5);gsub("R","N",$5);gsub("S","N",$5);gsub("B","N",$5);gsub("M","N",$5);gsub("V","N",$5);gsub("K","N",$5);gsub("D","N",$5);print $0}') > tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode3.vcf


sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.vcf > tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.vcf
vcftools --vcf tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.vcf --recode --out tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2 --recode-INFO AF
sed '20,21d' tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.recode.vcf | sed '21,109d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_tommo54",$0);print $0}' > tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.recode2.vcf

java -Xms10g -Xmx40g -jar ../picard.jar MergeVcfs I=tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode3.vcf I=tommo-54kjpn-20230626-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.recode2.vcf O=tommo-54kjpn-GRCh37-af.vcf

# jacquard expand tommo-54kjpn-GRCh37-af.vcf tommo-54kjpn-GRCh37-af.tsv

# Rでexoand
 library(tidyverse)
 library(vcfR)
 vcf <- read.vcfR("tommo-54kjpn-GRCh37-af.vcf")
 vcf_df <- vcfR2tidy(vcf)
 write.table(vcf_df$fix, "tommo-54kjpn-GRCh37-af2.tsv", sep = "\t", row.names = F, na = ".", quote=F)

tr -d \\r <tommo-54kjpn-GRCh37-af2.tsv> tommo-54kjpn-GRCh37-af.tsv


awk -F"\t" -v "OFS=\t" '{print $5}' tommo-54kjpn-GRCh37-af.tsv | awk -F"," -v "OFS=\t" '{print $1}' > tommo-54kjpn-GRCh37-af_1.1.tsv
awk -F"\t" -v "OFS=\t" '{print $9}' tommo-54kjpn-GRCh37-af.tsv | awk -F"," -v "OFS=\t" '{print $1}' > tommo-54kjpn-GRCh37-af_1.2.tsv

paste <(cut -f 1,2,3,4,6,7 tommo-54kjpn-GRCh37-af.tsv) tommo-54kjpn-GRCh37-af_1.1.tsv tommo-54kjpn-GRCh37-af_1.2.tsv | sed '1d' | awk -F"\t" -v "OFS=\t" '{print $1,$2,$3,$4,$7,$5,$6,"AF_tommo54="$8}' > tommo-54kjpn-GRCh37-af_1.tsv

grep ^# tommo-54kjpn-GRCh37-af.vcf > tommo-54kjpn-GRCh37-af.head.vcf
cat tommo-54kjpn-GRCh37-af.head.vcf <(cat tommo-54kjpn-GRCh37-af_1.tsv tommo-54kjpn-GRCh37-af_2.tsv tommo-54kjpn-GRCh37-af_3.tsv tommo-54kjpn-GRCh37-af_4.tsv tommo-54kjpn-GRCh37-af_5.tsv tommo-54kjpn-GRCh37-af_6.tsv tommo-54kjpn-GRCh37-af_7.tsv | awk -F"\t" -v "OFS=\t" '$5!="" {print $0}') > tommo-54kjpn-GRCh37-af2.vcf

cat tommo-54kjpn-GRCh37-af2.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > tommo-54kjpn-GRCh37-af2.sort.vcf
bcftools norm --check-ref w --threads 24 -f /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hs37d5.fa tommo-54kjpn-GRCh37-af2.sort.vcf -o tommo-54kjpn-GRCh37-af2.sort.norm.vcf




sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' tommo-38kjpn-20220630-GRCh37_lifted_from_GRCh38-af-autosome.vcf > tommo-38kjpn-20220630-GRCh37_lifted_from_GRCh38-af-autosome.4.2.vcf
vcftools --vcf tommo-38kjpn-20220630-GRCh37_lifted_from_GRCh38-af-autosome.4.2.vcf --recode --out tommo-38kjpn-20220630-GRCh37_lifted_from_GRCh38-af-autosome.4.2 --recode-INFO AF
sed '20,21d' tommo-38kjpn-20220630-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode.vcf | sed '21,118d' | sed '88,93d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_tommo38",$0);print $0}' > tommo-38kjpn-20220630-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode2.vcf
cat <(grep '^#' tommo-38kjpn-20220630-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode2.vcf) <(grep -v '^#' tommo-38kjpn-20220630-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode2.vcf | awk -F"\t" -v "OFS=\t" '{gsub("W","N",$4);gsub("Y","N",$4);gsub("H","N",$4);gsub("R","N",$4);gsub("S","N",$4);gsub("B","N",$4);gsub("M","N",$4);gsub("V","N",$4);gsub("K","N",$4);gsub("D","N",$4);gsub("W","N",$5);gsub("Y","N",$5);gsub("H","N",$5);gsub("R","N",$5);gsub("S","N",$5);gsub("B","N",$5);gsub("M","N",$5);gsub("V","N",$5);gsub("K","N",$5);gsub("D","N",$5);print $0}') > tommo-38kjpn-20220630-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode3.vcf


sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' tommo-38kjpn-20220929-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.vcf > tommo-38kjpn-20220929-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.vcf
vcftools --vcf tommo-38kjpn-20220929-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.vcf --recode --out tommo-38kjpn-20220929-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2 --recode-INFO AF
sed '20,21d' tommo-38kjpn-20220929-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.recode.vcf | sed '21,88d' | sed '88,93d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_tommo38",$0);print $0}' > tommo-38kjpn-20220929-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.recode2.vcf

java -Xms10g -Xmx40g -jar ../picard.jar MergeVcfs I=tommo-38kjpn-20220630-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode3.vcf I=tommo-38kjpn-20220929-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.recode2.vcf O=tommo-38kjpn-GRCh37-af.vcf

jacquard expand tommo-38kjpn-GRCh37-af.vcf tommo-38kjpn-GRCh37-af.tsv

awk -F"\t" -v "OFS=\t" '{print $5}' tommo-38kjpn-GRCh37-af.tsv | awk -F"," -v "OFS=\t" '{print $1}' > tommo-38kjpn-GRCh37-af_1.1.tsv
awk -F"\t" -v "OFS=\t" '{print $9}' tommo-38kjpn-GRCh37-af.tsv | awk -F"," -v "OFS=\t" '{print $1}' > tommo-38kjpn-GRCh37-af_1.2.tsv

paste <(cut -f 1,2,3,4,6,7 tommo-38kjpn-GRCh37-af.tsv) tommo-38kjpn-GRCh37-af_1.1.tsv tommo-38kjpn-GRCh37-af_1.2.tsv | sed '1d' | awk -F"\t" -v "OFS=\t" '{print $1,$2,$3,$4,$7,$5,$6,"AF_tommo38="$8}' > tommo-38kjpn-GRCh37-af_1.tsv

grep ^# tommo-38kjpn-GRCh37-af.vcf > tommo-38kjpn-GRCh37-af.head.vcf
cat tommo-38kjpn-GRCh37-af.head.vcf <(cat tommo-38kjpn-GRCh37-af_1.tsv tommo-38kjpn-GRCh37-af_2.tsv tommo-38kjpn-GRCh37-af_3.tsv tommo-38kjpn-GRCh37-af_4.tsv tommo-38kjpn-GRCh37-af_5.tsv tommo-38kjpn-GRCh37-af_6.tsv tommo-38kjpn-GRCh37-af_7.tsv | awk -F"\t" -v "OFS=\t" '$5!="" {print $0}') > tommo-38kjpn-GRCh37-af2.vcf

cat tommo-38kjpn-GRCh37-af2.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > tommo-38kjpn-GRCh37-af2.sort.vcf
bcftools norm --check-ref w --threads 24 -f /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hs37d5.fa tommo-38kjpn-GRCh37-af2.sort.vcf -o tommo-38kjpn-GRCh37-af2.sort.norm.vcf



sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-autosome.vcf > tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-autosome.4.2.vcf
vcftools --vcf tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-autosome.4.2.vcf --recode --out tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-autosome.4.2 --recode-INFO AF
sed '21d' tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode.vcf | sed '22,78d' | sed '88,93d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_tommo14",$0);print $0}' > tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode2.vcf
cat <(grep '^#' tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode2.vcf) <(grep -v '^#' tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode2.vcf | awk -F"\t" -v "OFS=\t" '{gsub("W","N",$4);gsub("Y","N",$4);gsub("H","N",$4);gsub("R","N",$4);gsub("S","N",$4);gsub("B","N",$4);gsub("M","N",$4);gsub("V","N",$4);gsub("K","N",$4);gsub("D","N",$4);gsub("W","N",$5);gsub("Y","N",$5);gsub("H","N",$5);gsub("R","N",$5);gsub("S","N",$5);gsub("B","N",$5);gsub("M","N",$5);gsub("V","N",$5);gsub("K","N",$5);gsub("D","N",$5);print $0}') > tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode3.vcf

sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.vcf > tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.vcf
vcftools --vcf tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.vcf --recode --out tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2 --recode-INFO AF
sed '21d' tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.recode.vcf | sed '22,78d' | sed '88,93d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_tommo14",$0);print $0}' > tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.recode2.vcf

java -Xms10g -Xmx40g -jar ../picard.jar MergeVcfs I=tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-autosome.4.2.recode3.vcf I=tommo-14KJPN-20211208r2-GRCh37_lifted_from_GRCh38-af-chrX_PAR3.4.2.recode2.vcf O=tommo-14kjpn-GRCh37-af.vcf

jacquard expand tommo-14kjpn-GRCh37-af.vcf tommo-14kjpn-GRCh37-af.tsv

awk -F"\t" -v "OFS=\t" '{print $5}' tommo-14kjpn-GRCh37-af.tsv | awk -F"," -v "OFS=\t" '{print $1}' > tommo-14kjpn-GRCh37-af_1.1.tsv
awk -F"\t" -v "OFS=\t" '{print $9}' tommo-14kjpn-GRCh37-af.tsv | awk -F"," -v "OFS=\t" '{print $1}' > tommo-14kjpn-GRCh37-af_1.2.tsv

paste <(cut -f 1,2,3,4,6,7 tommo-14kjpn-GRCh37-af.tsv) tommo-14kjpn-GRCh37-af_1.1.tsv tommo-14kjpn-GRCh37-af_1.2.tsv | sed '1d' | awk -F"\t" -v "OFS=\t" '{print $1,$2,$3,$4,$7,$5,$6,"AF_tommo14="$8}' > tommo-14kjpn-GRCh37-af_1.tsv

grep ^# tommo-14kjpn-GRCh37-af.vcf > tommo-14kjpn-GRCh37-af.head.vcf
cat tommo-14kjpn-GRCh37-af.head.vcf <(cat tommo-14kjpn-GRCh37-af_1.tsv tommo-14kjpn-GRCh37-af_2.tsv tommo-14kjpn-GRCh37-af_3.tsv tommo-14kjpn-GRCh37-af_4.tsv tommo-14kjpn-GRCh37-af_5.tsv tommo-14kjpn-GRCh37-af_6.tsv tommo-14kjpn-GRCh37-af_7.tsv | awk -F"\t" -v "OFS=\t" '$5!="" {print $0}') > tommo-14kjpn-GRCh37-af2.vcf

cat tommo-14kjpn-GRCh37-af2.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > tommo-14kjpn-GRCh37-af2.sort.vcf
bcftools norm --check-ref w --threads 24 -f ../hs37d5.fa tommo-14kjpn-GRCh37-af2.sort.vcf -o tommo-14kjpn-GRCh37-af2.sort.norm.vcf







cat simple_somatic_mutation.aggregated.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > simple_somatic_mutation.aggregated.sort.vcf
vcftools --vcf simple_somatic_mutation.aggregated.sort.vcf --recode --out simple_somatic_mutation.aggregated.sort --recode-INFO affected_donors
sed '2,3d' simple_somatic_mutation.aggregated.sort.recode.vcf | sed '3,6d' > simple_somatic_mutation.aggregated.sort.recode2.vcf


vcftools --gzvcf clinvar_20221129.vcf.gz --recode --out clinvar_20221129 --recode-INFO AF_ESP --recode-INFO AF_EXAC --recode-INFO AF_TGP --recode-INFO CLNSIG
sed '9,15d' clinvar_20221129.recode.vcf | sed '10,20d' > clinvar_20221129.recode2.vcf


vcftools --gzvcf clinvar_20240206.vcf.gz --recode --out clinvar_20240206 --recode-INFO AF_ESP --recode-INFO AF_EXAC --recode-INFO AF_TGP --recode-INFO CLNSIG --recode-INFO CLNSIGCONF --recode-INFO SCI
sed '9,15d' clinvar_20240206.recode.vcf | sed '11,32d' | sed '12d' > clinvar_20240206.recode2.vcf


vcftools --gzvcf GCF_000001405.25.gz --out GCF_000001405.25 --recode --recode-INFO FREQ --recode-INFO COMMON
sed '7,27d' GCF_000001405.25.recode.vcf | sed '9,16d' > GCF_000001405.25.recode2.vcf

bgzip -@ 24 GCF_000001405.25.recode2.vcf
tabix -p vcf GCF_000001405.25.recode2.vcf.gz
bcftools norm -m-any --check-ref w --thread 24 -f GRCh37.p13.fa GCF_000001405.25.recode2.vcf.gz -o GCF_000001405.25.recode2.norm.vcf
bcftools annotate --rename-chrs GCF_000001405.25_GRCh37.p13_assembly_report.chrnames --threads 24 -o GCF_000001405.25.recode2.norm2.vcf GCF_000001405.25.recode2.norm.vcf
bgzip -@ 24 GCF_000001405.25.recode2.norm2.vcf
tabix -p vcf GCF_000001405.25.recode2.norm2.vcf.gz

jacquard expand GCF_000001405.25.recode2.vcf GCF_000001405.25.tsv

awk -F"\t" -v "OFS=\t" '{print $5}' GCF_000001405.25.tsv | awk -F"," -v "OFS=\t" '{print $1}' > GCF_000001405.25_1.1.tsv
awk -F"\t" -v "OFS=\t" '{print $9}' GCF_000001405.25.tsv | awk -F"," -v "OFS=\t" '{print $1}' > GCF_000001405.25_1.2.tsv



vcftools --vcf Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.vcf --recode --out Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37 --recode-INFO SAMPLE_COUNT
sed '9,17d' Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.recode.vcf | sed '10,12d' | sed '34d' | awk -F"\t" -v "OFS=\t" '{gsub("SAMPLE_COUNT","CNT_target",$0);print $0}' > Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.recode2.vcf
grep "^#" Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.recode2.vcf > Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.recode2.head.vcf
grep -v -e '^\s*#' -e '^\s*$' Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.recode2.vcf > Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.recode2.variants.vcf
cat Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.recode2.variants.vcf | sort | uniq > Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.recode2.variants.uniq.vcf
cat Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.recode2.head.vcf Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.recode2.variants.uniq.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > Cosmic_CompleteTargetedScreensMutant_Normal_v98_GRCh37.recode2.variants.uniq.sort.vcf

vcftools --vcf Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.vcf --recode --out Cosmic_GenomeScreensMutant_Normal_v98_GRCh37 --recode-INFO SAMPLE_COUNT
sed '9,17d' Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.recode.vcf | sed '10,12d' | sed '35d' | awk -F"\t" -v "OFS=\t" '{gsub("SAMPLE_COUNT","CNT_Genome",$0);print $0}' > Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.recode2.vcf
grep "^#" Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.recode2.vcf > Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.recode2.head.vcf
grep -v -e '^\s*#' -e '^\s*$' Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.recode2.vcf > Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.recode2.variants.vcf
cat Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.recode2.variants.vcf | sort | uniq >Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.recode2.variants.uniq.vcf
cat Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.recode2.head.vcf Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.recode2.variants.uniq.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.recode2.variants.uniq.sort.vcf

vcftools --vcf Cosmic_NonCodingVariants_Normal_v98_GRCh37.vcf --recode --out Cosmic_NonCodingVariants_Normal_v98_GRCh37 --recode-INFO SAMPLE_COUNT
sed '9,13d' Cosmic_NonCodingVariants_Normal_v98_GRCh37.recode.vcf | sed '10,12d' | sed '35d' | awk -F"\t" -v "OFS=\t" '{gsub("SAMPLE_COUNT","CNT_NonCoding",$0);print $0}' > Cosmic_NonCodingVariants_Normal_v98_GRCh37.recode2.vcf
grep "^#" Cosmic_NonCodingVariants_Normal_v98_GRCh37.recode2.vcf > Cosmic_NonCodingVariants_Normal_v98_GRCh37.recode2.head.vcf
grep -v -e '^\s*#' -e '^\s*$' Cosmic_NonCodingVariants_Normal_v98_GRCh37.recode2.vcf > Cosmic_NonCodingVariants_Normal_v98_GRCh37.recode2.variants.vcf
cat Cosmic_NonCodingVariants_Normal_v98_GRCh37.recode2.variants.vcf | sort | uniq > Cosmic_NonCodingVariants_Normal_v98_GRCh37.recode2.variants.uniq.vcf
cat Cosmic_NonCodingVariants_Normal_v98_GRCh37.recode2.head.vcf Cosmic_NonCodingVariants_Normal_v98_GRCh37.recode2.variants.uniq.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > Cosmic_NonCodingVariants_Normal_v98_GRCh37.recode2.variants.uniq.sort.vcf



id=v99
vcftools --vcf CosmicCodingMuts.normal.${id}.vcf --recode --out CosmicCodingMuts.normal.${id} --recode-INFO CNT
sed '9,17d' CosmicCodingMuts.normal.${id}.recode.vcf | sed '35d' | awk -F"\t" -v "OFS=\t" '{gsub("CNT","CNT_coding",$0);print $0}' > CosmicCodingMuts.normal.${id}.recode2.vcf
grep "^#" CosmicCodingMuts.normal.${id}.recode2.vcf > CosmicCodingMuts.normal.${id}.recode2.head.vcf
grep -v -e '^\s*#' -e '^\s*$' CosmicCodingMuts.normal.${id}.recode2.vcf > CosmicCodingMuts.normal.${id}.recode2.variants.vcf
cat CosmicCodingMuts.normal.${id}.recode2.variants.vcf | sort | uniq > CosmicCodingMuts.normal.${id}.recode2.variants.uniq.vcf
cat CosmicCodingMuts.normal.${id}.recode2.head.vcf CosmicCodingMuts.normal.${id}.recode2.variants.uniq.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > CosmicCodingMuts.normal.${id}.recode2.variants.uniq.sort.vcf

vcftools --vcf CosmicNonCodingVariants.normal.${id}.vcf --recode --out CosmicNonCodingVariants.normal.${id} --recode-INFO CNT
sed '9,17d' CosmicNonCodingVariants.normal.${id}.recode.vcf | sed '35d' | awk -F"\t" -v "OFS=\t" '{gsub("CNT","CNT_noncoding",$0);print $0}' > CosmicNonCodingVariants.normal.${id}.recode2.vcf
grep "^#" CosmicNonCodingVariants.normal.${id}.recode2.vcf > CosmicNonCodingVariants.normal.${id}.recode2.head.vcf
grep -v -e '^\s*#' -e '^\s*$' CosmicNonCodingVariants.normal.${id}.recode2.vcf > CosmicNonCodingVariants.normal.${id}.recode2.variants.vcf
cat CosmicNonCodingVariants.normal.${id}.recode2.variants.vcf | sort | uniq > CosmicNonCodingVariants.normal.${id}.recode2.variants.uniq.vcf
cat CosmicNonCodingVariants.normal.${id}.recode2.head.vcf CosmicNonCodingVariants.normal.${id}.recode2.variants.uniq.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > CosmicNonCodingVariants.normal.${id}.recode2.variants.uniq.sort.vcf


jacquard expand CosmicCodingMuts.normal.${id}.vcf CosmicCodingMuts.normal.${id}.tsv
awk -F"\t" -v "OFS=\t" '{print $1,$2,$3,$4,$5,$18}' CosmicCodingMuts.normal.${id}.tsv | awk -F"\t" -v "OFS=\t" '$6!="." {print $0}' > CosmicCodingMuts.normal.${id}.OV.tsv
cut -f 6 CosmicCodingMuts.normal.${id}.OV.tsv | awk -F":" -v "OFS=\t" '{print "chr"$1":"$2"_"$3}' | tr "/" "_" > CosmicCodingMuts.normal.${id}.OV_1.tsv
awk -F"\t" -v "OFS=\t" '{print "chr"$1":"$2"_"$4"_"$5}' CosmicCodingMuts.normal.${id}.OV.tsv > CosmicCodingMuts.normal.${id}.OV_2.tsv
paste -d , CosmicCodingMuts.normal.${id}.OV_1.tsv CosmicCodingMuts.normal.${id}.OV_2.tsv > CosmicCodingMuts.normal.${id}.OV2.csv
(head -n +1 CosmicCodingMuts.normal.${id}.OV2.csv && tail -n +2 CosmicCodingMuts.normal.${id}.OV2.csv | sort) | uniq > CosmicCodingMuts.normal.${id}.OV3.csv

jacquard expand CosmicNonCodingVariants.normal.${id}.vcf CosmicNonCodingVariants.normal.${id}.tsv
awk -F"\t" -v "OFS=\t" '{print $1,$2,$3,$4,$5,$18}' CosmicNonCodingVariants.normal.${id}.tsv | awk -F"\t" -v "OFS=\t" '$6!="." {print $0}' > CosmicNonCodingVariants.normal.${id}.OV.tsv
cut -f 6 CosmicNonCodingVariants.normal.${id}.OV.tsv | awk -F":" -v "OFS=\t" '{print "chr"$1":"$2"_"$3}' | tr "/" "_" > CosmicNonCodingVariants.normal.${id}.OV_1.tsv
awk -F"\t" -v "OFS=\t" '{print "chr"$1":"$2"_"$4"_"$5}' CosmicNonCodingVariants.normal.${id}.OV.tsv > CosmicNonCodingVariants.normal.${id}.OV_2.tsv
paste -d , CosmicNonCodingVariants.normal.${id}.OV_1.tsv CosmicNonCodingVariants.normal.${id}.OV_2.tsv > CosmicNonCodingVariants.normal.${id}.OV2.csv
(head -n +1 CosmicNonCodingVariants.normal.${id}.OV2.csv && tail -n +2 CosmicNonCodingVariants.normal.${id}.OV2.csv | sort) | uniq > CosmicNonCodingVariants.normal.${id}.OV3.csv




vcftools --vcf Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37.vcf --recode --out Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37 --recode-INFO TARGETED_SCREEN_SAMPLE_COUNT
sed '9,17d' Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37.recode.vcf | sed '10,12d' | sed '34d' | awk -F"\t" -v "OFS=\t" '{gsub("TARGETED_SCREEN_SAMPLE_COUNT","CNT_target",$0);print $0}' > Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37.recode2.vcf
grep "^#" Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37.recode2.vcf > Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37.recode2.head.vcf
grep -v -e '^\s*#' -e '^\s*$' Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37.recode2.vcf > Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37.recode2.variants.vcf
cat Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37.recode2.variants.vcf | sort | uniq > Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37.recode2.variants.uniq.vcf
cat Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37.recode2.head.vcf Cosmic_CompleteTargetedScreensMutant_Normal_v100_GRCh37.recode2.variants.uniq.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > Cosmic_CompleteTargetedScreensMutant_Normal_v10_GRCh37.recode2.variants.uniq.sort.vcf

vcftools --vcf Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.vcf --recode --out Cosmic_GenomeScreensMutant_Normal_v100_GRCh37 --recode-INFO GENOME_SCREEN_SAMPLE_COUNT
sed '9,17d' Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.recode.vcf | sed '10,12d' | sed '35d' | awk -F"\t" -v "OFS=\t" '{gsub("GENOME_SCREEN_SAMPLE_COUNT","CNT_Genome",$0);print $0}' > Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.recode2.vcf
grep "^#" Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.recode2.vcf > Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.recode2.head.vcf
grep -v -e '^\s*#' -e '^\s*$' Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.recode2.vcf > Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.recode2.variants.vcf
cat Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.recode2.variants.vcf | sort | uniq >Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.recode2.variants.uniq.vcf
cat Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.recode2.head.vcf Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.recode2.variants.uniq.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > Cosmic_GenomeScreensMutant_Normal_v100_GRCh37.recode2.variants.uniq.sort.vcf

vcftools --vcf Cosmic_NonCodingVariants_Normal_v100_GRCh37.vcf --recode --out Cosmic_NonCodingVariants_Normal_v100_GRCh37 --recode-INFO SAMPLE_COUNT
sed '9,13d' Cosmic_NonCodingVariants_Normal_v100_GRCh37.recode.vcf | sed '10,12d' | sed '35d' | awk -F"\t" -v "OFS=\t" '{gsub("SAMPLE_COUNT","CNT_NonCoding",$0);print $0}' > Cosmic_NonCodingVariants_Normal_v100_GRCh37.recode2.vcf
grep "^#" Cosmic_NonCodingVariants_Normal_v100_GRCh37.recode2.vcf > Cosmic_NonCodingVariants_Normal_v100_GRCh37.recode2.head.vcf
grep -v -e '^\s*#' -e '^\s*$' Cosmic_NonCodingVariants_Normal_v100_GRCh37.recode2.vcf > Cosmic_NonCodingVariants_Normal_v100_GRCh37.recode2.variants.vcf
cat Cosmic_NonCodingVariants_Normal_v100_GRCh37.recode2.variants.vcf | sort | uniq > Cosmic_NonCodingVariants_Normal_v100_GRCh37.recode2.variants.uniq.vcf
cat Cosmic_NonCodingVariants_Normal_v100_GRCh37.recode2.head.vcf Cosmic_NonCodingVariants_Normal_v100_GRCh37.recode2.variants.uniq.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > Cosmic_NonCodingVariants_Normal_v100_GRCh37.recode2.variants.uniq.sort.vcf



unpigz -k BBJ_RIKEN_TMM_20200309_autosome.chr${id}.vqsr_99.5_99.0.region_flag.genotype_filtering.hwe_flag.annotated.vcf.gz
vcftools --vcf BBJ_RIKEN_TMM_20200309_autosome.chr${id}.vqsr_99.5_99.0.region_flag.genotype_filtering.hwe_flag.annotated.vcf --recode --out BBJ_RIKEN_TMM_20200309_autosome.chr${id} --recode-INFO AF
vcftools --vcf BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_nonPAR.vqsr_99.5_99.0.region_flag.genotype_filtering.hwe_flag.annotated.vcf --recode --out BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_nonPAR --recode-INFO AF
vcftools --vcf BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_PAR1.vqsr_99.5_99.0.region_flag.genotype_filtering.hwe_flag.annotated.vcf --recode --out BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_PAR1 --recode-INFO AF
vcftools --vcf BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_PAR2.vqsr_99.5_99.0.region_flag.genotype_filtering.hwe_flag.annotated.vcf --recode --out BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_PAR2 --recode-INFO AF


sed '105d' BBJ_RIKEN_TMM_20200309_autosome.chr${id}.recode.vcf | sed '106,131d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_GEMJ",$0);print $0}' > BBJ_RIKEN_TMM_20200309_autosome.chr${id}.recode2.vcf
sed '105d' BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_nonPAR.recode.vcf | sed '106,131d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_GEMJ",$0);print $0}' > BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_nonPAR.recode2.vcf
sed '105d' BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_PAR1.recode.vcf | sed '106,131d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_GEMJ",$0);print $0}' > BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_PAR1.recode2.vcf
sed '105d' BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_PAR2.recode.vcf | sed '106,131d' | awk -F"\t" -v "OFS=\t" '{gsub("AF","AF_GEMJ",$0);print $0}' > BBJ_RIKEN_TMM_20200309_chrX_PAR2.chrX_PAR2.recode2.vcf

ls *.recode2.vcf > list.txt
java -Xms10g -Xmx40g -jar ../picard.jar MergeVcfs I=list.txt O=BBJ_RIKEN_TMM_20200309_all.recode2.vcf

jacquard expand BBJ_RIKEN_TMM_20200309_all.recode2.vcf BBJ_RIKEN_TMM_20200309_all.tsv

awk -F"\t" -v "OFS=\t" '{print $5}' BBJ_RIKEN_TMM_20200309_all.tsv | awk -F"," -v "OFS=\t" '{print $1}' > BBJ_RIKEN_TMM_20200309_all_1.1.tsv
awk -F"\t" -v "OFS=\t" '{print $9}' BBJ_RIKEN_TMM_20200309_all.tsv | awk -F"," -v "OFS=\t" '{print $1}' > BBJ_RIKEN_TMM_20200309_all_1.2.tsv

paste <(cut -f 1,2,3,4,6,7 BBJ_RIKEN_TMM_20200309_all.tsv) BBJ_RIKEN_TMM_20200309_all_1.1.tsv BBJ_RIKEN_TMM_20200309_all_1.2.tsv | sed '1d' | awk -F"\t" -v "OFS=\t" '{print $1,$2,$3,$4,$7,$5,$6,"AF_GEMJ="$8}' > BBJ_RIKEN_TMM_20200309_all_1.tsv

grep ^# BBJ_RIKEN_TMM_20200309_all.recode2.vcf > BBJ_RIKEN_TMM_20200309_all.recode2.head.vcf
cat BBJ_RIKEN_TMM_20200309_all_1.tsv BBJ_RIKEN_TMM_20200309_all_2.tsv BBJ_RIKEN_TMM_20200309_all_3.tsv BBJ_RIKEN_TMM_20200309_all_4.tsv BBJ_RIKEN_TMM_20200309_all_5.tsv BBJ_RIKEN_TMM_20200309_all_6.tsv | awk -F"\t" -v "OFS=\t" '$5!="" {print $0}' | sort | uniq > BBJ_RIKEN_TMM_20200309_all2.tsv
cat BBJ_RIKEN_TMM_20200309_all.recode2.head.vcf BBJ_RIKEN_TMM_20200309_all2.tsv > BBJ_RIKEN_TMM_20200309_all2.vcf
cat BBJ_RIKEN_TMM_20200309_all2.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > BBJ_RIKEN_TMM_20200309_all2.sort.vcf
bcftools norm --check-ref w --threads 24 -f ./hs37d5.fa BBJ_RIKEN_TMM_20200309_all2.sort.vcf -o BBJ_RIKEN_TMM_20200309_all2.sort.norm.vcf



for id in {1..22} X Y
do
sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' chr_${id}_frequency.vcf > chr_${id}_frequency.4.2.vcf
vcftools --vcf chr_${id}_frequency.4.2.vcf --recode --out chr_${id}_frequency --recode-INFO AF_JGA_NGS --recode-INFO AF_JGA_SNP --recode-INFO AF_TOMMO --recode-INFO AF_HGVD --recode-INFO AF_GEM_J_WGA --recode-INFO AF_GNOMAD_EXOMES --recode-INFO AF_GNOMAD_GENOMES
cat chr_${id}_frequency.recode.vcf | sed '28,31d' | sed '29,31d' | sed '30,33d' | sed '31,33d' | sed '32,34d' | sed '33,35d' | sed '34,36d' | sed '35d' > chr_${id}_frequency.recode2.vcf
done


ls *.recode2.vcf > list.txt
java -Xms10g -Xmx40g -jar ../picard.jar MergeVcfs I=list.txt O=all.frequency.vcf
bcftools norm --check-ref w --threads 24 -f /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hs37d5.fa all.frequency.vcf -o all.frequency.norm.vcf





for id in {1..22} X Y; do sed '1d' chr_${id}_frequency.tsv | awk -F"\t" -v "OFS=\t" '{print "chr"$4":"$5"_"$6"_"$7,$1,$2,$11,$15,$22,$26,$30,$34,$38}' > chr_${id}_frequency2.tsv; done
head -n 1 chr_1_frequency.tsv | awk -F"\t" -v "OFS=\t" '{print "chr"$4":"$5"_"$6"_"$7,$1,$2,$11,$15,$22,$26,$30,$34,$38}' > head.tsv
cat *2.tsv | sort > all.frequency.sort.tsv
cat head.tsv all.frequency.sort.tsv | tr "\t" "," | tr "-" "." > all.frequency.sort.csv


sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' MGeND.vcf > MGeND.4.2.vcf
vcftools --vcf MGeND.4.2.vcf --recode --out MGeND --recode-INFO CS
sed '29,32d' MGeND.recode.vcf | sed '30,95d' > MGeND.recode2.vcf
bcftools norm --check-ref w --threads 24 -f /mnt/c/Users/Yusuke\ Kanemasa/Desktop/EPdata/database/hg19.fa MGeND.recode2.vcf -o MGeND.recode2.norm.vcf
awk -F"\t" -v "OFS=\t" '{gsub(" ", "_", $8); print}' MGeND.recode2.norm.vcf > MGeND.recode2.norm2.vcf



cut -f 1,2,4 CosmicMutantExport.tsv > CosmicMutantExport2.tsv
cat CosmicMutantExport2.tsv | sed -e '1d' | sort | uniq -c > CosmicMutantExport3.tsv
sed 's/^[ \t]*//' CosmicMutantExport3.tsv | tr " " "\t" > CosmicMutantExport4.tsv
cut -f 2 CosmicMutantExport4.tsv | cut -d _ -f 1 > CosmicMutantExport5.tsv
paste -d "\t" CosmicMutantExport4.tsv CosmicMutantExport5.tsv > CosmicMutantExport6.tsv

# CosmicMutantExport6.tsvをLibreOfficeでカウント数を逆順、遺伝子名で順に並べ替え　CosmicMutantExport7.csvとして保存

tr -d \\r <CosmicMutantExport7.csv> CosmicMutantExport8.csv
cat CosmicMutantExport8.csv | tr "," "\t" > CosmicMutantExport8.tsv
awk -F"\t" '!colname[$1]++{print}' CosmicMutantExport8.tsv > CosmicMutantExport9.tsv

# ENST変更　GNAS, AKT3, ERRFI1, MDM4, PBRM1, SNCAIP, CHKE2
cat CosmicMutantExport9.tsv | awk -F"\t" -v "OFS=\t" '{gsub("GNAS_ENST00000313949","GNAS_ENST00000371085",$3);gsub("ENST00000313949.7","ENST00000371085.3",$4);gsub("AKT3_ENST00000336199","AKT3",$3);gsub("ENST00000336199.5","ENST00000366539.1",$4);gsub("ERRFI1_ENST00000474874","ERRFI1",$3);gsub("ENST00000474874.1","ENST00000377482.5",$4);gsub("MDM4_ENST00000507825","MDM4",$3);gsub("ENST00000507825.2","ENST00000367182.3",$4);gsub("PBRM1_ENST00000337303","PBRM1_ENST00000296302",$3);gsub("ENST00000337303.4","ENST00000296302.7",$4);gsub("SNCAIP_ENST00000542191","SNCAIP",$3);gsub("ENST00000542191.1","ENST00000261368.8",$4);gsub("CHEK2_ENST00000544772","CHEK2_ENST00000328354",$3);gsub("ENST00000544772.1","ENST00000328354.6",$4);print $0}' | sort > CosmicMutantExport10.tsv

vi CosmicMutantExport10.tsv  # Gene    count   Gene_ENST       ENST    ID



awk -F"\t" -v "OFS=," '{print $1":"$2"_"$3"_"$4,$5,$6,$7,$8,$9,$10}' AlphaMissense_hg19.tsv > AlphaMissense_hg19_2.csv
sed '1,3d' AlphaMissense_hg19_2.csv > AlphaMissense_hg19_3.csv
(head -n +1 AlphaMissense_hg19_3.csv && tail -n +2 AlphaMissense_hg19_3.csv | sort -k 1,1 -t ',') > AlphaMissense_hg19_sort.csv