num=($(awk -F"\t" -v "OFS=\t" '{print $5}' RG_*.xml | grep "<no>" -c ))
for id in `seq 1 ${num}`; do awk -F"\t" -v "OFS=\t" '{print $5}' RG_*.xml | grep -v ^$ | awk -F"\t" -v "OFS=\t" '/<no>'${id}'<\/no>/,/source_id>/{print $0}' | grep -e "<no>" -e "<gene>" -e "<cds_change>" -e "<function>" -e "<change>" -e "<location>" -e "<copy_ratio>" -e "<allele-freq>" -e "<cosmic_status>" -e "<method>" | sed -e "s/&gt;/>/g" | sed 's/<no>//g' | sed 's/<gene>//g' | sed 's/<cds_change>//g' | sed 's/<function>//g' | sed 's/<change>//g' | sed 's/<location>//g' | sed 's/<copy_ratio>//g' | sed 's/<allele-freq>//g' | sed 's/<cosmic_status>//g' | sed 's/<method>//g' | sed 's/,//g' | awk '{sub("</.*", ""); print $0}' > marker${id}.tsv; row=($(wc -l marker${id}.tsv)); if [ $row = 20 ]; then head marker${id}.tsv > marker${id}_2.tsv; elif [ $row = 8 ]; then cat marker${id}.tsv | sed 3i- | sed 5i- > marker${id}_2.tsv; elif [ $row = 0 ]; then rm marker${id}.tsv; else cat marker${id}.tsv > marker${id}_2.tsv; fi; done
paste -s *_2.tsv > Book.tsv
cat ../database/marker2.head.tsv Book.tsv | tr ',' '\t' > Book2.tsv
join -t "$(printf '\011')" --header -a 1 -1 2 -2 1  -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 2.2  <(head -n +1 Book2.tsv && tail -n +2 Book2.tsv | sort -k 2,2) ../strand.uniq.sort.LF.tsv > Book3.tsv
awk -F"\t" -v "OFS=\t"  '{print $1,$2,$3,$3,$4,$5,$6,$7,$8,$9,$10,$11}' Book3.tsv | awk -F"\t" -v "OFS=\t" '{gsub(/[0-9]/,"",$4);gsub(/-/,"",$4);gsub(/+/,"",$4);gsub(/_/,"",$4);gsub(/c./,"",$4);print $0}' | tr "\t" "," > Book4.csv

awk -F"," -v "OFS=," '$4 ~ "del" && $4 !~ "ins" {print $0}' Book4.csv > Book4.del.csv
awk -F"," -v "OFS=," '$12=="+" {print $0}' Book4.del.csv > Book4.del.plus.csv
awk -F"," -v "OFS=," '$12=="-" {print $0}' Book4.del.csv > Book4.del.minus.csv
awk -F"," -v "OFS=," '{gsub("del","",$4);print $4}' Book4.del.plus.csv > Book4.del.plus2.csv
awk -F"," -v "OFS=," '{gsub("del","",$4);print $4}' Book4.del.minus.csv | tr A-Z a-z | awk -F"," -v "OFS=," '{gsub("c","G",$1);gsub("g","C",$1);gsub("t","A",$1);gsub("a","T",$1);print $0}' | rev > Book4.del.minus2.csv

awk -F"," -v "OFS=," '$4 ~ "ins" && $4 !~ "del" {print $0}' Book4.csv > Book4.ins.csv
awk -F"," -v "OFS=," '$12=="+" {print $0}' Book4.ins.csv > Book4.ins.plus.csv
awk -F"," -v "OFS=," '$12=="-" {print $0}' Book4.ins.csv > Book4.ins.minus.csv
awk -F"," -v "OFS=," '{gsub("ins","",$4);print $4}' Book4.ins.plus.csv > Book4.ins.plus2.csv
awk -F"," -v "OFS=," '{gsub("ins","",$4);print $4}' Book4.ins.minus.csv | tr A-Z a-z | awk -F"," -v "OFS=," '{gsub("c","G",$1);gsub("g","C",$1);gsub("t","A",$1);gsub("a","T",$1);print $0}' | rev > Book4.ins.minus2.csv

awk -F"," -v "OFS=," '$4 ~ "delins" {print $0}' Book4.csv > Book4.delins.csv
awk -F"," -v "OFS=," '$12=="+" {print $0}' Book4.delins.csv > Book4.delins.plus.csv
awk -F"," -v "OFS=," '$12=="-" {print $0}' Book4.delins.csv > Book4.delins.minus.csv
awk -F"," -v "OFS=," '{gsub("delins","",$4);print $4}' Book4.delins.plus.csv > Book4.delins.plus2.csv
awk -F"," -v "OFS=," '{gsub("delins","",$4);print $4}' Book4.delins.minus.csv | tr A-Z a-z | awk -F"," -v "OFS=," '{gsub("c","G",$1);gsub("g","C",$1);gsub("t","A",$1);gsub("a","T",$1);print $0}' | rev > Book4.delins.minus2.csv

awk -F"," -v "OFS=," '{print $7}' Book4.ins.plus.csv | awk -F":" -v "OFS=," '{print $1}' > Book4.ins.plus3.csv
awk -F"," -v "OFS=," '{print $7}' Book4.ins.plus.csv | awk -F":" -v "OFS=," '{print $2-1}' > Book4.ins.plus4.csv
paste -d : Book4.ins.plus3.csv Book4.ins.plus4.csv > Book4.ins.plus5.csv
paste -d , Book4.ins.plus.csv Book4.ins.plus2.csv Book4.ins.plus5.csv | awk -F"," -v "OFS=," '{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13,$14}' > Book4.ins.plus6.csv
awk -F"," -v "OFS=," '{print $7}' Book4.ins.minus.csv | awk -F":" -v "OFS=," '{print $1}' > Book4.ins.minus3.csv
awk -F"," -v "OFS=," '{print $7}' Book4.ins.minus.csv | awk -F":" -v "OFS=," '{print $2-1}' > Book4.ins.minus4.csv
paste -d : Book4.ins.minus3.csv Book4.ins.minus4.csv > Book4.ins.minus5.csv
paste -d , Book4.ins.minus.csv Book4.ins.minus2.csv Book4.ins.minus5.csv | awk -F"," -v "OFS=," '{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13,$14}' > Book4.ins.minus6.csv

awk -F"," -v "OFS=," '{print $7}' Book4.del.minus.csv | awk -F":" -v "OFS=," '{print $1}' > Book4.del.minus3.csv
awk -F"," -v "OFS=," '{print $7}' Book4.del.minus.csv | awk -F":" -v "OFS=," '{print $2-1}' > Book4.del.minus4.csv
paste -d : Book4.del.minus3.csv Book4.del.minus4.csv > Book4.del.minus5.csv
paste -d , Book4.del.minus.csv Book4.del.minus2.csv Book4.del.minus5.csv | awk -F"," -v "OFS=," '{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13,$14}' > Book4.del.minus6.csv
awk -F"," -v "OFS=," '{print $7}' Book4.del.plus.csv | awk -F":" -v "OFS=," '{print $1}' > Book4.del.plus3.csv
awk -F"," -v "OFS=," '{print $7}' Book4.del.plus.csv | awk -F":" -v "OFS=," '{print $2-1}' > Book4.del.plus4.csv
paste -d : Book4.del.plus3.csv Book4.del.plus4.csv > Book4.del.plus5.csv
paste -d , Book4.del.plus.csv Book4.del.plus2.csv Book4.del.plus5.csv | awk -F"," -v "OFS=," '{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13,$14}' > Book4.del.plus6.csv

awk -F"," -v "OFS=," '{print $7}' Book4.delins.minus.csv | awk -F":" -v "OFS=," '{print $1}' > Book4.delins.minus3.csv
awk -F"," -v "OFS=," '{print $7}' Book4.delins.minus.csv | awk -F":" -v "OFS=," '{print $2}' > Book4.delins.minus4.csv
paste -d : Book4.delins.minus3.csv Book4.delins.minus4.csv > Book4.delins.minus5.csv
paste -d , Book4.delins.minus.csv Book4.delins.minus2.csv Book4.delins.minus5.csv | awk -F"," -v "OFS=," '{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13,$14}' > Book4.delins.minus6.csv
awk -F"," -v "OFS=," '{print $7}' Book4.delins.plus.csv | awk -F":" -v "OFS=," '{print $1}' > Book4.delins.plus3.csv
awk -F"," -v "OFS=," '{print $7}' Book4.delins.plus.csv | awk -F":" -v "OFS=," '{print $2}' > Book4.delins.plus4.csv
paste -d : Book4.delins.plus3.csv Book4.delins.plus4.csv > Book4.delins.plus5.csv
paste -d , Book4.delins.plus.csv Book4.delins.plus2.csv Book4.delins.plus5.csv | awk -F"," -v "OFS=," '{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13,$14}' > Book4.delins.plus6.csv

cat Book4.del.plus6.csv Book4.del.minus6.csv | awk -F"," -v "OFS=," '{print "chr"$13,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > Book4.del2.csv
cat Book4.ins.plus6.csv Book4.ins.minus6.csv | awk -F"," -v "OFS=," '{print "chr"$13,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > Book4.ins2.csv
cat Book4.delins.plus6.csv Book4.delins.minus6.csv | awk -F"," -v "OFS=," '{print "chr"$13,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > Book4.delins2.csv

R CMD BATCH ../script/base_NCC.R

tr -d \\r <Book4.del2.join.csv> Book4.del2.join.tr.csv
tr -d \\r <Book4.ins2.join.csv> Book4.ins2.join.tr.csv
tr -d \\r <Book4.delins2.join.csv> Book4.delins2.join.tr.csv
awk -F"," -v "OFS=," '{print $1"_"$14$13"_"$14,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12}' Book4.del2.join.tr.csv > Book4.del2.join.tr2.csv
awk -F"," -v "OFS=," '{print $1"_"$14"_"$14$13,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12}' Book4.ins2.join.tr.csv > Book4.ins2.join.tr2.csv
awk -F"," -v "OFS=," '{print $1"_"$14"_"$13,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12}' Book4.delins2.join.tr.csv > Book4.delins2.join.tr2.csv

awk -F"," -v "OFS=," '$4 !~ "del" {print $0}' Book4.csv | awk -F"," -v "OFS=," '$4 !~ "ins" {print $0}' > Book5.csv
awk -F"," -v "OFS=," '$12=="+" {print $0}' Book5.csv > Book5.plus.csv
awk -F"," -v "OFS=," '$12=="-" {print $0}' Book5.csv > Book5.minus.csv
awk -F"," -v "OFS=," '{gsub(">","_",$4);print $0}' Book5.plus.csv > Book5.plus2.csv
awk -F"," -v "OFS=," '{gsub("C>G","G_C",$4);gsub("C>T","G_A",$4);gsub("C>A","G_T",$4);gsub("G>C","C_G",$4);gsub("G>T","C_A",$4);gsub("G>A","C_T",$4);gsub("T>G","A_C",$4);gsub("T>C","A_G",$4);gsub("T>A","A_T",$4);gsub("A>G","T_C",$4);gsub("A>C","T_G",$4);gsub("A>T","T_A",$4);print $0}' Book5.minus.csv > Book5.minus2.csv
head -n 1 Book5.csv > Book5.head.csv
cat Book5.head.csv Book5.plus2.csv Book5.minus2.csv | awk -F"," -v "OFS=," '{print "chr"$7"_"$4,$1,$2,$3,$5,$6,$8,$9,$10,$11,$12}' > Book6.csv
cat Book6.csv Book4.del2.join.tr2.csv Book4.ins2.join.tr2.csv Book4.delins2.join.tr2.csv | uniq > Book7.csv
cat Book7.csv | tr "," "\t" > Book7.tsv
awk -F"\t" -v "OFS=\t" '{print $1}' Book7.tsv | awk -F":" -v "OFS=\t" '{print $1}' | sed -e '1d' > Book7.chr.tsv
awk -F"\t" -v "OFS=\t" '{print $1}' Book7.tsv | awk -F":" -v "OFS=\t" '{print $2}' | awk -F"_" -v "OFS=\t" '{print $1}' | sed -e '1d' > Book7.pos.tsv
awk -F"\t" -v "OFS=\t" '{print $1}' Book7.tsv | awk -F":" -v "OFS=\t" '{print $2}' | awk -F"_" -v "OFS=\t" '{print $2}' | sed -e '1d' > Book7.REF.tsv
awk -F"\t" -v "OFS=\t" '{print $1}' Book7.tsv | awk -F":" -v "OFS=\t" '{print $2}' | awk -F"_" -v "OFS=\t" '{print $3}' | sed -e '1d' > Book7.ALT.tsv
paste Book7.chr.tsv Book7.pos.tsv Book7.REF.tsv Book7.ALT.tsv > Book8.tsv
awk -F"\t" -v "OFS=\t" '{print $1,$2,".",$3,$4,".",".",".","GT:DP","0|1:100"}' Book8.tsv > Book9.tsv
cat ../database/vcf_head2.tsv Book9.tsv > Book9.vcf

cat Book9.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > Book9_2.vcf
bgzip Book9_2.vcf
tabix -p vcf Book9_2.vcf.gz
bcftools norm -f ../database/hg19.fa -o Book9_3.vcf Book9_2.vcf.gz
awk -F"\t" -v "OFS=\t" '{gsub("c","C",$4);gsub("g","G",$4);gsub("t","T",$4);gsub("a","A",$4);gsub("c","C",$5);gsub("g","G",$5);gsub("t","T",$5);gsub("a","A",$5);print $0}' Book9_3.vcf > Book9_4.vcf
bgzip -d Book9_2.vcf.gz

awk -F"," -v "OFS=," '$0 !~/-/ {print $0}' Book9_4.vcf > Book9_5.vcf
java -jar ../snpEff/snpEff.jar hg19 Book9_5.vcf > Book9.snpeff.vcf
cat Book9.snpeff.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > Book9.snpeff.sort.vcf
java -Xms10g -Xmx30g -jar ../snpEff/SnpSift.jar annotate ../database/clinvar_20230604.recode2.vcf.gz Book9.snpeff.sort.vcf | 
java -jar ../snpEff/SnpSift.jar annotate ../database/tommo-38kjpn-GRCh37-af2.sort.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/tommo-14kjpn-GRCh37-af2.sort.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/tommo-8.3kjpn-20200831-af2.sort.norm.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.exomes.r2.1.1.sites.recode2.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.genomes.r2.1.1.sites.all.recode2.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/gnomad.genomes.v3.1.2.sites.all.recode3.hg19.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/BBJ_RIKEN_TMM_20200309_all2.sort.norm.vcf.gz |
java -jar ../snpEff/SnpSift.jar annotate ../database/GCF_000001405.25.recode2.norm2.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/CosmicCodingMuts.normal.v98.recode2.variants.uniq.sort.vcf.gz | 
java -jar ../snpEff/SnpSift.jar annotate ../database/CosmicNonCodingVariants.normal.v98.recode2.variants.uniq.sort.vcf.gz > Book9.snpeff.sort.annotate.vcf
java -Xms10g -Xmx30g -jar ../snpEff/SnpSift.jar extractFields Book9.snpeff.sort.annotate.vcf ANN[0].GENE ANN[0].IMPACT ANN[0].EFFECT ANN[0].HGVS_C ANN[0].HGVS_P > Book9.snpeff.sort.annotate.extraFields.tsv

R CMD BATCH ../script/vcfr_NCC.R

tr -d \\r <Book9.snpeff.sort.annotate.tsv> Book9.snpeff.sort.annotate.tr.tsv
paste -d "\t" Book9.snpeff.sort.annotate.tr.tsv Book9.snpeff.sort.annotate.extraFields.tsv | tr "," "/" | tr "\t" "," > Book9.snpeff.sort.annotate.csv
awk -F"," -v "OFS=," '{print $2":"$3"_"$5"_"$6,$4,$9,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35}' Book9.snpeff.sort.annotate.csv > Book9.snpeff.sort.annotate2.csv
join -t "," --header -a 1 -1 1 -2 1 -e "." -o auto <(head -n +1 Book9.snpeff.sort.annotate2.csv && tail -n +2 Book9.snpeff.sort.annotate2.csv | sort -k 1,1 -t ',') ../database/all.frequency.sort.csv > Book9.snpeff.sort.annotate3.csv
(head -n +1 Book7.csv && tail -n +2 Book7.csv | sort -k 1,1 -t ',') > Book7.sort.csv
join -t "," --header -a 1 -1 1 -2 1 -e "." -o auto Book7.sort.csv Book9.snpeff.sort.annotate3.csv > Book9.snpeff.sort.annotate6.csv
(head -n +1 Book9.snpeff.sort.annotate6.csv && tail -n +2 Book9.snpeff.sort.annotate6.csv | sort -t ',' -k 2,2n) > Book9.snpeff.sort.annotate6.sort.csv
rm Book3.tsv Book4.csv Book4.del.csv Book4.del.minus.csv Book4.del.minus2.csv Book4.del.minus3.csv Book4.del.plus.csv Book4.del.plus2.csv Book4.del.plus3.csv Book4.del2.csv Book4.del2.join.csv Book4.ins.csv Book4.ins.minus.csv Book4.ins.minus2.csv Book4.ins.minus3.csv Book4.ins.minus4.csv Book4.ins.minus5.csv Book4.ins.minus6.csv Book4.ins.plus.csv Book4.ins.plus2.csv Book4.ins.plus3.csv Book4.ins.plus4.csv Book4.ins.plus5.csv Book4.ins.plus6.csv Book4.ins2.csv Book4.ins2.join.csv Book5.csv Book5.head.csv Book5.minus.csv Book5.minus2.csv Book5.plus.csv Book5.plus2.csv Book6.csv Book7.ALT.tsv Book7.REF.tsv Book7.chr.tsv Book7.csv Book7.pos.tsv Book7.tsv Book8.tsv Book9.snpeff.sort.annotate.csv Book9.snpeff.sort.annotate.tsv Book9.snpeff.sort.annotate.vcf Book9.snpeff.sort.annotate2.csv Book9.snpeff.sort.annotate3.csv Book9.snpeff.sort.annotate6.csv Book9.snpeff.sort.vcf Book9.snpeff.vcf Book9.tsv Book9.vcf snpEff_genes.txt Book4.del.minus4.csv Book4.del.minus5.csv Book4.del.minus6.csv Book4.del.plus4.csv Book4.del.plus5.csv Book4.del.plus6.csv Book.tsv Book9_2.vcf Book2.tsv marker*.tsv Book4.del2.join.tr.csv Book4.del2.join.tr2.csv Book4.delins.csv Book4.delins.minus.csv Book4.delins.minus2.csv Book4.delins.minus3.csv Book4.delins.minus4.csv Book4.delins.minus5.csv Book4.delins.minus6.csv Book4.delins.plus.csv Book4.delins.plus2.csv Book4.delins.plus3.csv Book4.delins.plus4.csv Book4.delins.plus5.csv Book4.delins.plus6.csv Book4.delins2.csv Book4.delins2.join.csv Book4.delins2.join.tr.csv Book4.delins2.join.tr2.csv Book4.ins2.join.tr.csv Book4.ins2.join.tr2.csv Book7.sort.csv Book9.snpeff.sort.annotate.extraFields.tsv Book9.snpeff.sort.annotate.tr.tsv Book9_2.vcf.gz.tbi Book9_3.vcf Book9_4.vcf Book9_5.vcf
