sed -i $'1s/^\uFEFF//' copy-number-alteration.csv
cat copy-number-alteration.csv | tr -d '\r' | tr -d '"' > copy-number-alteration_tr.csv
join -t "," --header -a 1 -1 4 -2 1 -e "-" -o auto <(head -n +1 copy-number-alteration_tr.csv && tail -n +2 copy-number-alteration_tr.csv | sort -k 4,4 -t ',') ../database/cancerGeneList.sort2.csv > copy-number-alteration_tr2.csv
(head -n +1 copy-number-alteration_tr2.csv && tail -n +2 copy-number-alteration_tr2.csv | sort -t ',' -k 2,2n) > copy-number-alteration_tr3.csv
rm copy-number-alteration_tr.csv copy-number-alteration_tr2.csv
