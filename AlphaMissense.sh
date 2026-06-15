(head -n +1 short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort2.csv && tail -n +2 short-variant_tr.change2.combined6.snpeff.sort.annotate8.sort2.csv | sort -t ',' -k 1,1) > short-variant_tr.change2.combined6.snpeff.sort.annotate10.csv

join -t "," --header -a 1 -1 1 -2 1 -e "." -o auto short-variant_tr.change2.combined6.snpeff.sort.annotate10.csv ../database/AlphaMissense_hg19_sort.csv > short-variant_Alphamissense.csv



awk -F"," -v "OFS=," '$55 != "." {print$0}' short-variant_Alphamissense.csv > short-variant_Alphamissense2.csv
diff <(sed '1d' short-variant_Alphamissense2.csv | cut -d , -f 11) <(sed '1d' short-variant_Alphamissense2.csv | cut -d , -f 55)