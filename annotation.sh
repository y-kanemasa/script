

cat short-variant.csv | tr -d '\r' | tr -d '"' > short-variant_tr.csv
vi -c 'set nobomb' -c 'wq!' short-variant_tr.csv

# subcolnal variant があるか否か
subclonal_number=`grep subclonal short-variant_tr.csv | wc -l`
if [[ $subclonal_number >0 ]] ; then
	echo "There are some subcolnal variants in short-variant_tr.csv"
	R CMD BATCH ../script/subclonal.R
	tr -d \\r < short-variant_tr.cut.csv > short-variant_tr.cut2.csv
	tr -d \\r < short-variant_tr.sub.csv > short-variant_tr.sub2.csv
	rm short-variant_tr.csv short-variant_tr.cut.csv short-variant_tr.sub.csv
	mv short-variant_tr.cut2.csv short-variant_tr.csv
fi

awk -F"," -v "OFS=," '{print $1,$2,$3,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' short-variant_tr.csv | awk -F"," -v "OFS=," '{gsub(/[0-9]/,"",$4);gsub(/-/,"",$4);gsub(/+/,"",$4);gsub(/_/,"",$4);print $0}' > short-variant_tr.change.csv
awk -F"," -v "OFS=," '{print $1,$2,$3,$4,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' short-variant_tr.change.csv > short-variant_tr.change2.csv
awk -F"," -v "OFS=," '$14=="+" {print $0}' short-variant_tr.change2.csv > short-variant_tr.change2.plus.csv
awk -F"," -v "OFS=," '$14=="-" {print $0}' short-variant_tr.change2.csv > short-variant_tr.change2.minus.csv
awk -F"," -v "OFS=," '{gsub(">","_",$5);print $0}' short-variant_tr.change2.plus.csv > short-variant_tr.change2.plus.change.csv
awk -F"," -v "OFS=," '{print $5}' short-variant_tr.change2.minus.csv > short-variant_tr.change2.minus.change.csv
awk -F">" -v "OFS=," '{print $1}' short-variant_tr.change2.minus.change.csv > short-variant_tr.change2.minus.change_1.csv
awk -F">" -v "OFS=," '{print $2}' short-variant_tr.change2.minus.change.csv > short-variant_tr.change2.minus.change_2.csv
cat short-variant_tr.change2.minus.change_1.csv | tr A-Z a-z | awk -F"," -v "OFS=," '{gsub("c","G",$1);gsub("g","C",$1);gsub("t","A",$1);gsub("a","T",$1);print $0}' | rev > short-variant_tr.change2.minus.change2_1.csv
cat short-variant_tr.change2.minus.change_2.csv | tr A-Z a-z | awk -F"," -v "OFS=," '{gsub("c","G",$1);gsub("g","C",$1);gsub("t","A",$1);gsub("a","T",$1);print $0}' | rev > short-variant_tr.change2.minus.change2_2.csv
paste -d "_" short-variant_tr.change2.minus.change2_1.csv short-variant_tr.change2.minus.change2_2.csv > short-variant_tr.change2.minus.change2.csv
paste -d , short-variant_tr.change2.minus.csv short-variant_tr.change2.minus.change2.csv | awk -F"," -v "OFS=," '{print $1,$2,$3,$4,$16,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' > short-variant_tr.change2.minus.change3.csv
head -n 1 short-variant_tr.change2.csv > short-variant_tr.change2.head.csv
cat short-variant_tr.change2.head.csv short-variant_tr.change2.plus.change.csv short-variant_tr.change2.minus.change3.csv > short-variant_tr.change2.combined.csv
awk -F"," -v "OFS=," '{print $11,$5,$1,$2,$3,$4,$6,$7,$8,$9,$10,$12,$13,$14,$15}' short-variant_tr.change2.combined.csv > short-variant_tr.change2.combined2.csv
awk -F"," -v "OFS=," '$6 ~ /del/ {print $0}' short-variant_tr.change2.combined2.csv > short-variant_tr.change2.combined2.del.csv
awk -F"," -v "OFS=," '$14=="+" {print $0}' short-variant_tr.change2.combined2.del.csv > short-variant_tr.change2.combined2.del.plus.csv
awk -F"," -v "OFS=," '$14=="-" {print $0}' short-variant_tr.change2.combined2.del.csv > short-variant_tr.change2.combined2.del.minus.csv
awk -F"," -v "OFS=," '{gsub("del","",$6);print $6}' short-variant_tr.change2.combined2.del.plus.csv > short-variant_tr.change2.combined2.del.plus2.csv
awk -F"," -v "OFS=," '{gsub("del","",$6);print $6}' short-variant_tr.change2.combined2.del.minus.csv | tr A-Z a-z | awk -F"," -v "OFS=," '{gsub("c","G",$1);gsub("g","C",$1);gsub("t","A",$1);gsub("a","T",$1);print $0}' | rev > short-variant_tr.change2.combined2.del.minus2.csv
awk -F"," -v "OFS=," '$6 ~ /ins/ {print $0}' short-variant_tr.change2.combined2.csv > short-variant_tr.change2.combined2.ins.csv
awk -F"," -v "OFS=," '$14=="+" {print $0}' short-variant_tr.change2.combined2.ins.csv > short-variant_tr.change2.combined2.ins.plus.csv
awk -F"," -v "OFS=," '$14=="-" {print $0}' short-variant_tr.change2.combined2.ins.csv > short-variant_tr.change2.combined2.ins.minus.csv
awk -F"," -v "OFS=," '{gsub("ins","",$6);print $6}' short-variant_tr.change2.combined2.ins.plus.csv > short-variant_tr.change2.combined2.ins.plus2.csv
awk -F"," -v "OFS=," '{gsub("ins","",$6);print $6}' short-variant_tr.change2.combined2.ins.minus.csv | tr A-Z a-z | awk -F"," -v "OFS=," '{gsub("c","G",$1);gsub("g","C",$1);gsub("t","A",$1);gsub("a","T",$1);print $0}' | rev > short-variant_tr.change2.combined2.ins.minus2.csv
paste -d , short-variant_tr.change2.combined2.del.plus.csv short-variant_tr.change2.combined2.del.plus2.csv > short-variant_tr.change2.combined2.del.plus3.csv
paste -d , short-variant_tr.change2.combined2.del.minus.csv short-variant_tr.change2.combined2.del.minus2.csv > short-variant_tr.change2.combined2.del.minus3.csv
paste -d , short-variant_tr.change2.combined2.ins.plus.csv short-variant_tr.change2.combined2.ins.plus2.csv > short-variant_tr.change2.combined2.ins.plus3.csv
paste -d , short-variant_tr.change2.combined2.ins.minus.csv short-variant_tr.change2.combined2.ins.minus2.csv > short-variant_tr.change2.combined2.ins.minus3.csv
cat short-variant_tr.change2.combined2.del.plus3.csv short-variant_tr.change2.combined2.del.minus3.csv > short-variant_tr.change2.combined2.del2.csv
cat short-variant_tr.change2.combined2.ins.plus3.csv short-variant_tr.change2.combined2.ins.minus3.csv > short-variant_tr.change2.combined2.ins2.csv
awk -F"," -v "OFS=," '$2 ~ /T/ || $2 ~ /C/ || $2 ~ /G/ || $2 ~ /A/ {print $0}' short-variant_tr.change2.combined2.del2.csv > short-variant_tr.change2.combined2.del2.short.csv
awk -F"," -v "OFS=," '$2 !~ /T/ && $2 !~ /C/ && $2 !~ /G/ && $2 !~ /A/ {print $0}' short-variant_tr.change2.combined2.del2.csv > short-variant_tr.change2.combined2.del2.large.csv
cut -d , -f 5 short-variant_tr.change2.combined2.del2.large.csv | sed 's/.*del\([0-9]*\).*/\1/' > short-variant_tr.change2.combined2.del2.large.number.csv
paste -d , short-variant_tr.change2.combined2.del2.large.csv short-variant_tr.change2.combined2.del2.large.number.csv > short-variant_tr.change2.combined2.del2.large2.csv
awk -F"," -v "OFS=," '$6 !~ /del/ && $6 !~ /ins/ {print $0}' short-variant_tr.change2.combined2.csv | awk -F"," -v "OFS=," '{print $6,$1,$2,$3,$4,$5,$7,$8,$9,$10,$11,$12,$13,$14,$15}' | grep ^">" | awk -F"," -v "OFS=," '{print $2,$3,$4,$5,$6,$1,$7,$8,$9,$10,$11,$12,$13,$14,$15}' > short-variant_tr.change2.combined2.splice.csv
awk -F"," -v "OFS=," '$14=="+" {print $0}' short-variant_tr.change2.combined2.splice.csv > short-variant_tr.change2.combined2.splice.plus.csv
awk -F"," -v "OFS=," '$14=="-" {print $0}' short-variant_tr.change2.combined2.splice.csv > short-variant_tr.change2.combined2.splice.minus.csv
awk -F"," -v "OFS=," '{gsub(">","",$6);print $6}' short-variant_tr.change2.combined2.splice.plus.csv > short-variant_tr.change2.combined2.splice.plus2.csv
awk -F"," -v "OFS=," '{gsub(">","",$6);print $6}' short-variant_tr.change2.combined2.splice.minus.csv | tr A-Z a-z | awk -F"," -v "OFS=," '{gsub("c","G",$1);gsub("g","C",$1);gsub("t","A",$1);gsub("a","T",$1);print $0}' | rev > short-variant_tr.change2.combined2.splice.minus2.csv
paste -d , short-variant_tr.change2.combined2.splice.plus.csv short-variant_tr.change2.combined2.splice.plus2.csv > short-variant_tr.change2.combined2.splice.plus3.csv
paste -d , short-variant_tr.change2.combined2.splice.minus.csv short-variant_tr.change2.combined2.splice.minus2.csv > short-variant_tr.change2.combined2.splice.minus3.csv
cat short-variant_tr.change2.combined2.splice.plus3.csv short-variant_tr.change2.combined2.splice.minus3.csv > short-variant_tr.change2.combined2.splice2.csv
awk -F"," -v "OFS=," '{print $5}' short-variant_tr.change2.combined2.splice2.csv | awk '{ sub(">.*$",""); print $0; }' > short-variant_tr.change2.combined2.splice3.csv
splice_number=`wc -l short-variant_tr.change2.combined2.splice3.csv | cut -d " " -f 1`

if [[ $splice_number >0 ]] ; then
	for id in `seq $splice_number`
	do
	calc1=$(cut -d "_" -f 1 short-variant_tr.change2.combined2.splice3.csv | sed -n ${id}P | tee /dev/tty)
	echo $((calc1)) >> short-variant_tr.change2.combined2.splice3.calc1.csv
	calc2=$(cut -d "_" -f 2 short-variant_tr.change2.combined2.splice3.csv | sed -n ${id}P | tee /dev/tty)
	echo $((calc2)) >> short-variant_tr.change2.combined2.splice3.calc2.csv
	done
	paste -d , short-variant_tr.change2.combined2.splice3.calc1.csv short-variant_tr.change2.combined2.splice3.calc2.csv | awk -F"," -v "OFS=," '{print $2-$1}' > short-variant_tr.change2.combined2.splice4.csv
	paste -d , short-variant_tr.change2.combined2.splice2.csv short-variant_tr.change2.combined2.splice4.csv > short-variant_tr.change2.combined2.splice5.csv
	rm short-variant_tr.change2.combined2.splice3.calc1.csv short-variant_tr.change2.combined2.splice3.calc2.csv short-variant_tr.change2.combined2.splice4.csv
else 
	echo "" > short-variant_tr.change2.combined2.splice5.csv
fi

R CMD BATCH ../script/base.R

tr -d \\r <short-variant_tr.change2.combined2.del2.short.join.csv> short-variant_tr.change2.combined2.del2.short.join.tr.csv
awk -F"," -v "OFS=," '{print $1"_"$17$16"_"$17,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' short-variant_tr.change2.combined2.del2.short.join.tr.csv > short-variant_tr.change2.combined2.del2.short.join2.csv
tr -d \\r <short-variant_tr.change2.combined2.del2.large.join.csv> short-variant_tr.change2.combined2.del2.large.join.tr.csv
awk -F"," -v "OFS=," '{print $1"_"$16"_"$17,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' short-variant_tr.change2.combined2.del2.large.join.tr.csv > short-variant_tr.change2.combined2.del2.large.join2.csv
tr -d \\r <short-variant_tr.change2.combined2.ins2.join.csv> short-variant_tr.change2.combined2.ins2.join.tr.csv
awk -F"," -v "OFS=," '{print $1"_"$17"_"$17$16,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' short-variant_tr.change2.combined2.ins2.join.tr.csv > short-variant_tr.change2.combined2.ins2.join2.csv
tr -d \\r <short-variant_tr.change2.combined2.splice5.join.csv > short-variant_tr.change2.combined2.splice5.join.tr.csv
awk -F"," -v "OFS=," '{print $1"_"$17"_"$16,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' short-variant_tr.change2.combined2.splice5.join.tr.csv > short-variant_tr.change2.combined2.splice5.join2.csv
awk -F"," -v "OFS=," '$6 !~ /del/ && $6 !~ /ins/ {print $0}' short-variant_tr.change2.combined2.csv | awk -F"," -v "OFS=," '{print $6,$1,$2,$3,$4,$5,$7,$8,$9,$10,$11,$12,$13,$14,$15}' | grep -v ^">" | awk -F"," -v "OFS=," '{print $2,$3,$4,$5,$6,$1,$7,$8,$9,$10,$11,$12,$13,$14,$15}' | awk -F"," -v "OFS=," '{print $1"_"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' > short-variant_tr.change2.combined3.csv
cat short-variant_tr.change2.combined3.csv short-variant_tr.change2.combined2.del2.short.join2.csv short-variant_tr.change2.combined2.del2.large.join2.csv short-variant_tr.change2.combined2.ins2.join2.csv short-variant_tr.change2.combined2.splice5.join2.csv > short-variant_tr.change2.combined4.csv
cat short-variant_tr.change2.combined4.csv | tr "," "\t" > short-variant_tr.change2.combined4.tsv
awk -F"\t" -v "OFS=\t" '{print $1}' short-variant_tr.change2.combined4.tsv | awk -F":" -v "OFS=\t" '{print $1}' | sed -e '1d' > short-variant_tr.change2.combined4.chr.tsv
awk -F"\t" -v "OFS=\t" '{print $1}' short-variant_tr.change2.combined4.tsv | awk -F":" -v "OFS=\t" '{print $2}' | awk -F"_" -v "OFS=\t" '{print $1}' | sed -e '1d' > short-variant_tr.change2.combined4.pos.tsv
awk -F"\t" -v "OFS=\t" '{print $1}' short-variant_tr.change2.combined4.tsv | awk -F":" -v "OFS=\t" '{print $2}' | awk -F"_" -v "OFS=\t" '{print $2}' | sed -e '1d' > short-variant_tr.change2.combined4.REF.tsv
awk -F"\t" -v "OFS=\t" '{print $1}' short-variant_tr.change2.combined4.tsv | awk -F":" -v "OFS=\t" '{print $2}' | awk -F"_" -v "OFS=\t" '{print $3}' | sed -e '1d' > short-variant_tr.change2.combined4.ALT.tsv
paste short-variant_tr.change2.combined4.chr.tsv short-variant_tr.change2.combined4.pos.tsv short-variant_tr.change2.combined4.REF.tsv short-variant_tr.change2.combined4.ALT.tsv > short-variant_tr.change2.combined5.tsv
awk -F"\t" -v "OFS=\t" '{print $1,$2,".",$3,$4,".",".",".","."}' short-variant_tr.change2.combined5.tsv > short-variant_tr.change2.combined6.tsv
cat ../database/vcf_head2.tsv <(awk -F"\t" -v "OFS=\t" '{print $1, $2, ".", $4, $5, ".", ".", ".", "GT:DP", "0|1:100"}' short-variant_tr.change2.combined6.tsv) > short-variant_tr.change2.combined6.vcf
awk -F"," -v "OFS=," '$0 !~/-/ {print $0}' short-variant_tr.change2.combined6.vcf > short-variant_tr.change2.combined6_2.vcf
awk -F"," -v "OFS=\t" '{print $13,$7}' short-variant_tr.csv | sed -e '1d' | sort > short-variant_tr.refseq.tsv
join -t "$(printf '\011')" -a 1 -1 1 -2 1 -e "." -o auto short-variant_tr.refseq.tsv ../database/Refseq_ENST.tsv > short-variant_tr.refseq.join.tsv

# 対応するENST numberのないものがあるか
blank_number=`cat short-variant_tr.refseq.join.tsv | awk -F"\t" -v "OFS=\t" '$3=="." {print $0}' | wc -l`
if [[ $blank_number >0 ]] ; then
	echo "There are some blanks in annotation_ENST.tsv"
fi

cat short-variant_tr.change2.combined6_2.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > short-variant_tr.change2.combined6_3.vcf
bgzip short-variant_tr.change2.combined6_3.vcf
tabix -p vcf short-variant_tr.change2.combined6_3.vcf.gz
bcftools norm -f ../database/hg19.fa -o short-variant_tr.change2.combined6_4.vcf short-variant_tr.change2.combined6_3.vcf.gz
awk -F"\t" -v "OFS=\t" '{gsub("c","C",$4);gsub("g","G",$4);gsub("t","T",$4);gsub("a","A",$4);gsub("c","C",$5);gsub("g","G",$5);gsub("t","T",$5);gsub("a","A",$5);print $0}' short-variant_tr.change2.combined6_4.vcf > short-variant_tr.change2.combined6_5.vcf
bgzip -d short-variant_tr.change2.combined6_3.vcf.gz

# ATM, ALOX12B, WT1, MUTYH, MYD88, NOTCH1, IGF1R, AKT3, PIK3CB, PARP3, MYCのENSTを変更
cat short-variant_tr.refseq.join.tsv | cut -f 3 | sort | uniq > annotation_ENST.tsv
cat annotation_ENST.tsv | awk -F"\t" -v "OFS=\t" '{gsub("ENST00000675843","ENST00000278616",$1);gsub("ENST00000647874","ENST00000319144",$1);gsub("ENST00000452863","ENST00000332351",$1);gsub("ENST00000355498","ENST00000372115",$1);gsub("ENST00000650905","ENST00000396334",$1);gsub("ENST00000651671","ENST00000277541",$1);gsub("ENST00000650285","ENST00000268035",$1);gsub("ENST00000673466","ENST00000336199",$1);gsub("ENST00000674063","ENST00000477593",$1);gsub("ENST00000398755","ENST00000417220",$1);gsub("ENST00000621592","ENST00000377970",$1);gsub("ENST00000646891","ENST00000288602",$1);print $0}' > annotation_ENST2.tsv
java -jar ../snpEff/snpEff.jar -onlyTr annotation_ENST2.tsv GRCh37.87 short-variant_tr.change2.combined6_5.vcf > short-variant_tr.change2.combined6.snpeff.vcf
cat short-variant_tr.change2.combined6.snpeff.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > short-variant_tr.change2.combined6.snpeff.sort.vcf
java -Xms10g -Xmx30g -jar ../snpEff/SnpSift.jar annotate ../database/clinvar_*.hg19.recode2.vcf.gz short-variant_tr.change2.combined6.snpeff.sort.vcf | 
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
java -jar ../snpEff/SnpSift.jar annotate ../database/MGeND.recode2.norm2.vcf.gz > short-variant_tr.change2.combined6.snpeff.sort.annotate.vcf
java -Xms10g -Xmx30g -jar ../snpEff/SnpSift.jar extractFields short-variant_tr.change2.combined6.snpeff.sort.annotate.vcf ANN[0].GENE ANN[0].IMPACT ANN[0].EFFECT ANN[0].HGVS_C ANN[0].HGVS_P > short-variant_tr.change2.combined6.snpeff.sort.annotate.extraFields.tsv

R CMD BATCH ../script/vcfr.R

tr -d \\r <short-variant_tr.change2.combined6.snpeff.sort.annotate.tsv> short-variant_tr.change2.combined6.snpeff.sort.annotate.tr.tsv
paste -d "\t" short-variant_tr.change2.combined6.snpeff.sort.annotate.tr.tsv short-variant_tr.change2.combined6.snpeff.sort.annotate.extraFields.tsv | tr "," "/" | tr "\t" "," > short-variant_tr.change2.combined6.snpeff.sort.annotate.csv

# 6引く faf95_g2ex, faf95_g2ge, faf95_g3ge_popmax, SwappedAlleles, ReverseComplementedAlleles
paste -d , <(awk -F"," -v "OFS=," '{print $2":"$3"_"$5"_"$6,$4}' short-variant_tr.change2.combined6.snpeff.sort.annotate.csv) <(cut -d , -f 9- short-variant_tr.change2.combined6.snpeff.sort.annotate.csv) | awk -F, '{for (i=1; i<=NF; i++) if (i!=4 && i!=5 && i!=18 && i!=23 && i!=26 && i!=29 && i!=30) printf "%s%s", $i, (i==NF?"\n":",")}' > short-variant_tr.change2.combined6.snpeff.sort.annotate2.csv


# Left alignで修正されたバリアント
left_align_number=`diff <(grep -v ^'#' short-variant_tr.change2.combined6_3.vcf) <(grep -v ^'#' short-variant_tr.change2.combined6_5.vcf) | wc -l`
if [[ $left_align_number >0 ]] ; then
	echo "There are some variants needed to be letf-aligned"
	diff <(grep -v ^'#' short-variant_tr.change2.combined6_3.vcf) <(grep -v ^'#' short-variant_tr.change2.combined6_5.vcf) | grep ^'<' | awk -F"\t" -v "OFS=\t" '{print $1":"$2"_"$4"_"$5}' | cut -d ' ' -f 2 > left_align_before.tsv
	diff <(grep -v ^'#' short-variant_tr.change2.combined6_3.vcf) <(grep -v ^'#' short-variant_tr.change2.combined6_5.vcf) | grep ^'>' | awk -F"\t" -v "OFS=\t" '{print $1":"$2"_"$4"_"$5}' | cut -d ' ' -f 2 > left_align_after.tsv
	paste -d ' ' left_align_before.tsv left_align_after.tsv | awk -F" " -v "OFS= " '{print "s/"$1"/"$2"/"}' > left_align_sed.txt
	sed -f left_align_sed.txt short-variant_tr.change2.combined4.csv > short-variant_tr.change2.combined4_2.csv
	join -t "," --header -a 1 -1 1 -2 1 -e "." -o auto <(head -n +1 short-variant_tr.change2.combined4_2.csv && tail -n +2 short-variant_tr.change2.combined4_2.csv | sort -k 1,1 -t ',') <(head -n +1 short-variant_tr.change2.combined6.snpeff.sort.annotate2.csv && tail -n +2 short-variant_tr.change2.combined6.snpeff.sort.annotate2.csv | sort -k 1,1 -t ',') > short-variant_tr.change2.combined6.snpeff.sort.annotate6.csv
	rm left_align_before.tsv left_align_after.tsv short-variant_tr.change2.combined4_2.csv left_align_sed.txt
else
	join -t "," --header -a 1 -1 1 -2 1 -e "." -o auto <(head -n +1 short-variant_tr.change2.combined4.csv && tail -n +2 short-variant_tr.change2.combined4.csv | sort -k 1,1 -t ',') <(head -n +1 short-variant_tr.change2.combined6.snpeff.sort.annotate2.csv && tail -n +2 short-variant_tr.change2.combined6.snpeff.sort.annotate2.csv | sort -k 1,1 -t ',') > short-variant_tr.change2.combined6.snpeff.sort.annotate6.csv
fi


cut -d , -f 9 short-variant_tr.change2.combined6.snpeff.sort.annotate6.csv | awk -F"\t" -v "OFS=\t" '{gsub("MLL2","KMT2D",$1);gsub("MLL","KMT2A",$1);gsub("PARK2","PRKN",$1);gsub("MYCL1","MYCL",$1);gsub("MRE11A","MRE11",$1);gsub("FAM123B","AMER1",$1);gsub("WHSC1L1","NSD3",$1);gsub("WHSC1","NSD2",$1);print $0}' > short-variant_tr.change2.combined6.snpeff.sort.annotate6.gene.csv
paste -d , short-variant_tr.change2.combined6.snpeff.sort.annotate6.csv short-variant_tr.change2.combined6.snpeff.sort.annotate6.gene.csv > short-variant_tr.change2.combined6.snpeff.sort.annotate7.csv

# 56を変更
join -t "," --header -a 1 -1 56 -2 1 -e "." -o auto <(head -n +1 short-variant_tr.change2.combined6.snpeff.sort.annotate7.csv && tail -n +2 short-variant_tr.change2.combined6.snpeff.sort.annotate7.csv | sort -k 56,56 -t ',') ../database/cancerGeneList3.csv > short-variant_tr.change2.combined6.snpeff.sort.annotate8.csv
(head -n +1 short-variant_tr.change2.combined6.snpeff.sort.annotate8.csv && tail -n +2 short-variant_tr.change2.combined6.snpeff.sort.annotate8.csv | sort -t ',' -k 13,13 -k 10,10) | cut -d , -f 2- > short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort.csv


#gene名がSNPeffとCOSMICで同一か
awk -F"," -v "OFS=," '$2 !~ /-/ {print$0}' short-variant_tr.change2.combined6.snpeff.sort.annotate8.csv > short-variant_tr.change2.combined6.snpeff.sort.annotate9.csv
# 61を変更
diff <(sed '1d' short-variant_tr.change2.combined6.snpeff.sort.annotate9.csv | cut -d , -f 1) <(sed '1d' short-variant_tr.change2.combined6.snpeff.sort.annotate9.csv | cut -d , -f 61)


if [[ -e short-variant_tr.sub2.csv ]]; then
    (head -n +1 short-variant_tr.change2.combined6.snpeff.sort.annotate8.csv && tail -n +2 short-variant_tr.change2.combined6.snpeff.sort.annotate8.csv | sort -t ',' -n -k 3,3) | cut -d , -f 2- > short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort2.csv
    paste -d , short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort2.csv short-variant_tr.sub2.csv > short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort3.csv
    (head -n +1 short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort3.csv && tail -n +2 short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort3.csv | sort -t ',' -k 12,12 -k 9,9) > short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort4.csv
    rm short-variant_tr.sub2.csv short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort.csv short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort2.csv short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort3.csv
    mv short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort4.csv short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort.csv
fi


rm annotation_ENST.tsv annotation_ENST2.tsv base.Rout short-variant_tr.change.csv short-variant_tr.change2.combined.csv short-variant_tr.change2.combined2.csv short-variant_tr.change2.combined2.del.csv short-variant_tr.change2.combined2.del.minus.csv short-variant_tr.change2.combined2.del.minus2.csv short-variant_tr.change2.combined2.del.minus3.csv short-variant_tr.change2.combined2.del.plus.csv short-variant_tr.change2.combined2.del.plus2.csv short-variant_tr.change2.combined2.del.plus3.csv short-variant_tr.change2.combined2.del2.csv short-variant_tr.change2.combined2.del2.large.csv short-variant_tr.change2.combined2.del2.large.join.csv short-variant_tr.change2.combined2.del2.large.join.tr.csv short-variant_tr.change2.combined2.del2.large.join2.csv short-variant_tr.change2.combined2.del2.large.number.csv short-variant_tr.change2.combined2.del2.large2.csv short-variant_tr.change2.combined2.del2.short.csv short-variant_tr.change2.combined2.del2.short.join.csv short-variant_tr.change2.combined2.del2.short.join.tr.csv short-variant_tr.change2.combined2.del2.short.join2.csv short-variant_tr.change2.combined2.ins.csv short-variant_tr.change2.combined2.ins.minus.csv short-variant_tr.change2.combined2.ins.minus2.csv short-variant_tr.change2.combined2.ins.minus3.csv short-variant_tr.change2.combined2.ins.plus.csv short-variant_tr.change2.combined2.ins.plus2.csv short-variant_tr.change2.combined2.ins.plus3.csv short-variant_tr.change2.combined2.ins2.csv short-variant_tr.change2.combined2.ins2.join.csv short-variant_tr.change2.combined2.ins2.join.tr.csv short-variant_tr.change2.combined2.ins2.join2.csv short-variant_tr.change2.combined2.splice.csv short-variant_tr.change2.combined2.splice.minus.csv short-variant_tr.change2.combined2.splice.minus2.csv short-variant_tr.change2.combined2.splice.minus3.csv short-variant_tr.change2.combined2.splice.plus.csv short-variant_tr.change2.combined2.splice.plus2.csv short-variant_tr.change2.combined2.splice.plus3.csv short-variant_tr.change2.combined2.splice2.csv short-variant_tr.change2.combined2.splice3.csv short-variant_tr.change2.combined2.splice5.csv short-variant_tr.change2.combined2.splice5.join.csv short-variant_tr.change2.combined2.splice5.join.tr.csv short-variant_tr.change2.combined2.splice5.join2.csv short-variant_tr.change2.combined3.csv short-variant_tr.change2.combined4.ALT.tsv short-variant_tr.change2.combined4.REF.tsv short-variant_tr.change2.combined4.chr.tsv short-variant_tr.change2.combined4.csv short-variant_tr.change2.combined4.pos.tsv short-variant_tr.change2.combined4.tsv short-variant_tr.change2.combined5.tsv short-variant_tr.change2.combined6.snpeff.sort.annotate.csv short-variant_tr.change2.combined6.snpeff.sort.annotate.extraFields.tsv short-variant_tr.change2.combined6.snpeff.sort.annotate.tsv short-variant_tr.change2.combined6.snpeff.sort.annotate.vcf short-variant_tr.change2.combined6.snpeff.sort.annotate2.csv short-variant_tr.change2.combined6.snpeff.sort.annotate6.csv short-variant_tr.change2.combined6.snpeff.sort.annotate6.gene.csv short-variant_tr.change2.combined6.snpeff.sort.annotate7.csv short-variant_tr.change2.combined6.snpeff.sort.annotate8.csv short-variant_tr.change2.combined6.snpeff.sort.annotate9.csv short-variant_tr.change2.combined6.snpeff.sort.vcf short-variant_tr.change2.combined6.snpeff.vcf short-variant_tr.change2.combined6.tsv short-variant_tr.change2.combined6.vcf short-variant_tr.change2.combined6_2.vcf short-variant_tr.change2.combined6_3.vcf short-variant_tr.change2.combined6_3.vcf.gz.tbi short-variant_tr.change2.combined6_4.vcf short-variant_tr.change2.combined6_5.vcf short-variant_tr.change2.csv short-variant_tr.change2.head.csv short-variant_tr.change2.minus.change.csv short-variant_tr.change2.minus.change2.csv short-variant_tr.change2.minus.change2_1.csv short-variant_tr.change2.minus.change2_2.csv short-variant_tr.change2.minus.change3.csv short-variant_tr.change2.minus.change_1.csv short-variant_tr.change2.minus.change_2.csv short-variant_tr.change2.minus.csv short-variant_tr.change2.plus.change.csv short-variant_tr.change2.plus.csv short-variant_tr.csv short-variant_tr.refseq.tsv snpEff_genes.txt short-variant_tr.change2.combined6.snpeff.sort.annotate.tr.tsv