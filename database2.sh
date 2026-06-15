
# ClinVar

unpigz -k clinvar_*.vcf.gz

awk 'BEGIN {
  print "variant\tCLNSIG_other\tONC_other\tCLNSIGCONF_other"
}
!/^##/ {
  split($8, arr, ";"); clnSig=""; oncVal=""; clnSigConf="";
  for(i in arr) {
    if (arr[i] ~ /^CLNSIG=/) { sub(/^CLNSIG=/,"",arr[i]); clnSig=arr[i] }
    else if (arr[i] ~ /^ONC=/) { sub(/^ONC=/,"",arr[i]); oncVal=arr[i] }
    else if (arr[i] ~ /^CLNSIGCONF=/) { sub(/^CLNSIGCONF=/,"",arr[i]); clnSigConf=arr[i] }
  }
  print "chr"$1":"$2"_"$4"_"$5 "\t" clnSig "\t" oncVal "\t" clnSigConf
}' clinvar_*.vcf | sed '2d' > clinvar.tsv

(head -n +1 clinvar.tsv && tail -n +2 clinvar.tsv | sort -t$'\t' -k 1,1) > clinvar.sort.tsv



# 同じアミノ酸置換
#シェルスクリプトを用いて、あるファイルの##から始まる行を削除し、その後"chr", 1列目, ":", ２列目, "_", ４列目, "_", 5列目とつないで１列目としてください。
#もとのファイルの8列目には";"で区切られた文字列が複数あります。この中で、"GENE="、"CNT="で始まる文字列を抽出して、それぞれ新しいファイルの2,3列目としてください。だたし"GENE"、"CNT="は削除してください。さらに、”HGVSP="で始まる文字列はさらに":"で区切られています。この":"で区切られた文字列のうち2番目の文字列を抽出し、新しいファイルの4列目にしてください。さらに"HGVSC="で始まる文字列のうち、"="と"."の間にある文字列を抽出して新しいファイルの5列目にしてください。
#新しいファイルの列名は、variant, gene, CNT, AAchange, ENSTとしてください。

awk '
  BEGIN {
    FS = "\t"; OFS = "\t";
    # ヘッダーを出力
    print "variant","gene","CNT","AAchange","ENST";
  }
  # "##" で始まる行はスキップ
  /^##/ { next }
  {
    # 8列目を";"で分割
    split($8, info_arr, ";");
    
    # 変数を初期化
    gene = "";
    cnt  = "";
    aachange = "";
    enst = "";
    
    # 分割した配列を1つずつチェック
    for(i in info_arr) {
      # GENE= から始まるもの → GENE名抽出
      if (info_arr[i] ~ /^GENE=/) {
        # "GENE=" は5文字なので6文字目以降を取り出す
        gene = substr(info_arr[i], 6);
      }
      # CNT= から始まるもの → CNT値抽出
      else if (info_arr[i] ~ /^CNT=/) {
        cnt = substr(info_arr[i], 5);
      }
      # HGVSP= から始まるもの → ":" 区切りの2番目
      else if (info_arr[i] ~ /^HGVSP=/) {
        tmp_str = substr(info_arr[i], 7);     # "HGVSP=" を除去
        split(tmp_str, tmp_arr, ":");
        # ":" が全く含まれないケースがあった場合には要注意
        # ここでは最低でも2要素あると仮定
        aachange = tmp_arr[2];
      }
      # HGVSC= から始まるもの → "=" から "." の手前
      else if (info_arr[i] ~ /^HGVSC=/) {
        tmp_str = substr(info_arr[i], 7);     # "HGVSC=" を除去
        # 最初の "." 以降を削る (それ以降を空文字に)
        sub(/\..*/, "", tmp_str);
        enst = tmp_str;
      }
    }

    # 1列目(variant)を作成: chr + 1列目 : 2列目 _ 4列目 _ 5列目
    variant = "chr" $1 ":" $2 "_" $4 "_" $5;
    
    # 出力
    print variant, gene, cnt, aachange, enst;
  }
' CosmicCodingMuts.normal.v99.vcf | sed '2d' | awk 'BEGIN {
    FS = OFS = "\t"  
}
{
     # (1) 第4列の末尾が "[A-Za-z]{3}fsTer[0-9]+" なら丸ごと "fs" に置換
    sub(/[A-Za-z]{3}fsTer[0-9]+$/, "fs", $4)
    # (2) 第4列内にある "Ter" を "*" に置換（すべて置換したい場合は gsub）
    gsub(/Ter/, "*", $4)
    print
}' > CosmicCodingMuts.normal.v99.tsv


awk 'BEGIN {
  FS=OFS="\t"
}
{
  # 6列目を「2列目_4列目」として生成
  $6 = $2"_"$4

  # 6列目をキーに3列目を加算する
  sum[$6] += $3

  # 元の行を保存しておく
  lines[NR] = $0
}
END {
  # 全行に対して集計結果を付け加えて出力
  for(i=1; i<=NR; i++) {
    split(lines[i], f, FS)
    # 7列目に合計値を付与（キーは6列目）
    f[7] = sum[f[6]]

    # 4列目が""の場合、その行は出力しない
    if(f[4] != "") {
      print f[1], f[2], f[3], f[4], f[5], f[6], f[7]
    }
  }
}' CosmicCodingMuts.normal.v99.tsv | sed '1s/\<0\>/CNT2/g' > CosmicCodingMuts.normal.v99_2.tsv


awk 'BEGIN { FS=OFS="\t" }
     NR==1 { print; next }       # ヘッダーはそのまま出力
     $2 !~ /_ENST0000/ { print }' CosmicCodingMuts.normal.v99_2.tsv > CosmicCodingMuts.normal.v99_3.tsv


#6列目が共通の文字列である行に対して次の処理を行ってください。これらのグループとなっている行のそれぞれに対して、行を複製して8列目にそのグループの行の1列目を記載してください。たとえば2つの行からなるグループについては、1つの行に対して8列目が異なる2行が作成され、トータルで4行になります。3つの行からなるグループについては、1つの行に対して8列目が異なる3行が作成され、トータルで9行になります。
#6列目に共通の文字列を有しない行に対しては、8列目にその行の1列目の文字列を記載してください。
#8列目の列名は"other_variant"としてください。

awk '
BEGIN {
  FS  = "\t"     # 入力の区切り文字: タブ
  OFS = "\t"     # 出力の区切り文字: タブ
}

NR == 1 {
  # ---- [1] ヘッダー行の処理 ----
  # 先頭行を分割して8列目を"other_variant"に書き換えて出力
  split($0, arr, FS)
  
  # 万一、列数が 8 未満なら拡張して 8 列目を作る
  # （少なくとも arr[8] に "other_variant" を入れる）
  arr[8] = "other_variant"
  
  # 再構築して出力
  $0 = arr[1]
  for (i = 2; i <= length(arr); i++) {
    $0 = $0 OFS arr[i]
  }
  print $0
  
  # この後の処理はスキップして次行へ
  next
}

{
  # ---- [2] データ行の読み込み ----
  # 行全体を記憶 (NR 行目)
  lines[NR] = $0
  
  # 6列目をキーとし、行番号 (NR) をカンマ区切りで連結
  key = $6
  if (key in groupIndices) {
    # 既に値があれば末尾にカンマ＋NRを追加
    groupIndices[key] = groupIndices[key] "," NR
  } else {
    # 初回ならそのまま格納
    groupIndices[key] = NR
  }
}

END {
  # ---- [3] グループごとに展開し、(n × n) 行を出力 ----
  for (k in groupIndices) {
    # まず、そのグループの行番号リストを取り出す
    # 例: "2,5,9" など
    split(groupIndices[k], idxs, ",")
    n = length(idxs)
    
    # グループ内の各行 i(1..n) をベースに出力し、
    # 8列目をグループ内の各行 m(1..n) の 1列目で複製
    for (i = 1; i <= n; i++) {
      # i番目の行を取り出す
      line_i = lines[idxs[i]]
      split(line_i, arrI, FS)   # 分割して配列 arrI に
      
      for (m = 1; m <= n; m++) {
        # m番目の行から1列目を抜き出す
        line_m = lines[idxs[m]]
        split(line_m, arrM, FS)
        
        # arrI[8] を arrM[1] に差し替え
        arrI[8] = arrM[1]
        
        # 再構築して出力
        #   列数が arrI 内で増減している場合に備え、
        #   max_col を計算してループで出力する方法を使う
        out = arrI[1]
        max_col = (length(arrI) > 8 ? length(arrI) : 8)
        for (x = 2; x <= max_col; x++) {
          out = out OFS arrI[x]
        }
        print out
      }
    }
  }
}' CosmicCodingMuts.normal.v99_3.tsv > CosmicCodingMuts.normal.v99_4.tsv


(head -n +1 CosmicCodingMuts.normal.v99_4.tsv && tail -n +2 CosmicCodingMuts.normal.v99_4.tsv | sort -t$'\t' -k 8,8) > CosmicCodingMuts.normal.v99_4.sort.tsv

join -t$'\t' --header -a 1 -1 8 -2 1 -e "." -o auto CosmicCodingMuts.normal.v99_4.sort.tsv clinvar.sort.tsv > CosmicCodingMuts.normal.v99_4.sort.join.tsv

(head -n +1 CosmicCodingMuts.normal.v99_4.sort.join.tsv && tail -n +2 CosmicCodingMuts.normal.v99_4.sort.join.tsv | sort -t$'\t' -k 2,2) | cut -f 2,7,8,9,10,11 | uniq | cut -f 2,3,4,5,6 > CosmicCodingMuts.normal.v99_4.sort.join2.tsv


#シェルスクリプトを用いて、あるファイルで1列目が同一の文字列である行に対して操作を行います。
#3,4,5列目には文字列が入っています。1列目が同一の文字列である行のグループに対して、これらの行の3列目の文字列を"/"で結合した文字列を出力し、新しく6列目として各行に格納してください。同様の処理を4列目、5列目に対して行い、新しく7列目、8列目としてください。
#1列目に共通の文字列を有しない行に対しては、6,7,8列目にその行の3,4,5列の文字列を記載してください。
#最終的に6,7,8列目の列名をCLNSIG_other2, ONC_other2, CLNSIGCONF_other2としてください。

awk '
BEGIN {
  FS  = "\t"   # 入力の区切り文字 (タブ)
  OFS = "\t"   # 出力の区切り文字 (タブ)
}

NR == 1 {
  #
  # -- [1] ヘッダー行の処理 --
  #
  # 入力の1行目がヘッダーと想定し、配列に取り込む
  split($0, headerCols, FS)
  
  # 6,7,8 列目の新しい列名を設定
  headerCols[6] = "CLNSIG_other2"
  headerCols[7] = "ONC_other2"
  headerCols[8] = "CLNSIGCONF_other2"
  
  # 再構築して出力（最大で8列になる）
  max_hcol = (length(headerCols) > 8 ? length(headerCols) : 8)
  out = headerCols[1]
  for(i = 2; i <= max_hcol; i++){
    out = out OFS headerCols[i]
  }
  print out
  
  # 次行以降の処理へ
  next
}

#
# -- [2] データ行の読み込み・グルーピング --
#
{
  # 行全体を lines[NR] に保存 (後で再出力に使う)
  lines[NR] = $0
  
  # 1列目をキーとする
  key = $1
  
  # 同じ1列目をもつ行番号の一覧（カンマ区切り）を groupIndices[key] に管理
  if (key in groupIndices) {
    groupIndices[key] = groupIndices[key] "," NR
  } else {
    groupIndices[key] = NR
  }
  
  # 3列目, 4列目, 5列目を "/結合" で保持していく
  if (key in groupCol3) {
    groupCol3[key] = groupCol3[key] "/" $3
    groupCol4[key] = groupCol4[key] "/" $4
    groupCol5[key] = groupCol5[key] "/" $5
  } else {
    groupCol3[key] = $3
    groupCol4[key] = $4
    groupCol5[key] = $5
  }
}

END {
  #
  # -- [3] グループごとに行を再構築して出力 --
  #
  # groupIndices には "同じ1列目を持つ行番号リスト" が入っている
  for (key in groupIndices) {
    # カンマ区切りの行番号リストを split
    nlist = groupIndices[key]
    split(nlist, idxs, ",")  # idxs[] = {行番号1, 行番号2, ...}
    
    # グループに含まれる 3,4,5 列目を "/結合" した文字列
    merged3 = groupCol3[key]
    merged4 = groupCol4[key]
    merged5 = groupCol5[key]
    
    # グループ内の各行について再出力
    for (i = 1; i <= length(idxs); i++) {
      rowNum = idxs[i]
      # 事前に保存しておいた元の行を分割
      split(lines[rowNum], arr, FS)
      
      # 6,7,8列目として merged3, merged4, merged5 を格納
      arr[6] = merged3
      arr[7] = merged4
      arr[8] = merged5
      
      # 再構築して出力
      max_col = (length(arr) > 8 ? length(arr) : 8)
      out = arr[1]
      for (j = 2; j <= max_col; j++){
        out = out OFS arr[j]
      }
      print out
    }
  }
}' CosmicCodingMuts.normal.v99_4.sort.join2.tsv > CosmicCodingMuts.normal.v99_4.sort.join3.tsv

cut -f 1,2,6,7,8 CosmicCodingMuts.normal.v99_4.sort.join3.tsv > CosmicCodingMuts.normal.v99_4.sort.join4.tsv

(head -n +1 CosmicCodingMuts.normal.v99_4.sort.join4.tsv && tail -n +2 CosmicCodingMuts.normal.v99_4.sort.join4.tsv | sort -t$'\t' -k 1,1) | uniq | tr "\t" "," > CosmicCodingMuts.normal.v99.clinvar_sameAminochange.csv

(head -n +1 CosmicCodingMuts.normal.v99.clinvar_sameAminochange.csv && tail -n +2 CosmicCodingMuts.normal.v99.clinvar_sameAminochange.csv | sort -t ',' -k 1,1) > CosmicCodingMuts.normal.v99.clinvar_sameAminochange2.csv

# CNT2は同じアミノ酸置換変異でのCNT総数
# CLNSIG_other2	ONC_other2 CLNSIGCONF_other2は同じアミノ酸置換変異でのClinVar登録





# アミノ酸部位が同一

cat CosmicCodingMuts.normal.v99.tsv | awk '$4 !~ /fs/' | awk '$4 !~ /%3D/' | awk '$4 !~ /*/' | 
 awk 'BEGIN {
  FS = OFS = "\t"
}
{
  # 第1列が「数字 + 3文字以上アルファベット」で終わるか？
  # (.*[0-9]+) = 途中なんでも + 数字(1文字以上)
  # ([A-Za-z]{3,}) = 3文字以上のアルファベット
  # $ = 文字列末尾
  if (match($4, /(.*[0-9]+)([A-Za-z]{3,})$/, arr)) {
      # arr[1] には「末尾の数字を含む部分」までが入る
      # arr[2] には「末尾アルファベット3文字以上」が入る
      # → アルファベットだけを除外し、数字部分は残す
      $4 = arr[1]
  }
  print
}' > CosmicCodingMuts.normal.v99_5.tsv


awk 'BEGIN {
  FS=OFS="\t"
}
{
  # 6列目を「2列目_4列目」として生成
  $6 = $2"_"$4

  # 6列目をキーに3列目を加算する
  sum[$6] += $3

  # 元の行を保存しておく
  lines[NR] = $0
}
END {
  # 全行に対して集計結果を付け加えて出力
  for(i=1; i<=NR; i++) {
    split(lines[i], f, FS)
    # 7列目に合計値を付与（キーは6列目）
    f[7] = sum[f[6]]

    # 4列目が""の場合、その行は出力しない
    if(f[4] != "") {
      print f[1], f[2], f[3], f[4], f[5], f[6], f[7]
    }
  }
}' CosmicCodingMuts.normal.v99_5.tsv | sed '1s/\<0\>/CNT3/g' > CosmicCodingMuts.normal.v99_6.tsv


awk 'BEGIN { FS=OFS="\t" }
     NR==1 { print; next }       # ヘッダーはそのまま出力
     $2 !~ /_ENST0000/ { print }' CosmicCodingMuts.normal.v99_6.tsv > CosmicCodingMuts.normal.v99_7.tsv


#6列目が共通の文字列である行に対して次の処理を行ってください。これらのグループとなっている行のそれぞれに対して、行を複製して8列目にそのグループの行の1列目を記載してください。たとえば2つの行からなるグループについては、1つの行に対して8列目が異なる2行が作成され、トータルで4行になります。3つの行からなるグループについては、1つの行に対して8列目が異なる3行が作成され、トータルで9行になります。
#6列目に共通の文字列を有しない行に対しては、8列目にその行の1列目の文字列を記載してください。
#8列目の列名は"other_variant"としてください。

awk '
BEGIN {
  FS  = "\t"     # 入力の区切り文字: タブ
  OFS = "\t"     # 出力の区切り文字: タブ
}

NR == 1 {
  # ---- [1] ヘッダー行の処理 ----
  # 先頭行を分割して8列目を"other_variant"に書き換えて出力
  split($0, arr, FS)
  
  # 万一、列数が 8 未満なら拡張して 8 列目を作る
  # （少なくとも arr[8] に "other_variant" を入れる）
  arr[8] = "other_variant"
  
  # 再構築して出力
  $0 = arr[1]
  for (i = 2; i <= length(arr); i++) {
    $0 = $0 OFS arr[i]
  }
  print $0
  
  # この後の処理はスキップして次行へ
  next
}

{
  # ---- [2] データ行の読み込み ----
  # 行全体を記憶 (NR 行目)
  lines[NR] = $0
  
  # 6列目をキーとし、行番号 (NR) をカンマ区切りで連結
  key = $6
  if (key in groupIndices) {
    # 既に値があれば末尾にカンマ＋NRを追加
    groupIndices[key] = groupIndices[key] "," NR
  } else {
    # 初回ならそのまま格納
    groupIndices[key] = NR
  }
}

END {
  # ---- [3] グループごとに展開し、(n × n) 行を出力 ----
  for (k in groupIndices) {
    # まず、そのグループの行番号リストを取り出す
    # 例: "2,5,9" など
    split(groupIndices[k], idxs, ",")
    n = length(idxs)
    
    # グループ内の各行 i(1..n) をベースに出力し、
    # 8列目をグループ内の各行 m(1..n) の 1列目で複製
    for (i = 1; i <= n; i++) {
      # i番目の行を取り出す
      line_i = lines[idxs[i]]
      split(line_i, arrI, FS)   # 分割して配列 arrI に
      
      for (m = 1; m <= n; m++) {
        # m番目の行から1列目を抜き出す
        line_m = lines[idxs[m]]
        split(line_m, arrM, FS)
        
        # arrI[8] を arrM[1] に差し替え
        arrI[8] = arrM[1]
        
        # 再構築して出力
        #   列数が arrI 内で増減している場合に備え、
        #   max_col を計算してループで出力する方法を使う
        out = arrI[1]
        max_col = (length(arrI) > 8 ? length(arrI) : 8)
        for (x = 2; x <= max_col; x++) {
          out = out OFS arrI[x]
        }
        print out
      }
    }
  }
}' CosmicCodingMuts.normal.v99_7.tsv > CosmicCodingMuts.normal.v99_8.tsv

(head -n +1 CosmicCodingMuts.normal.v99_8.tsv && tail -n +2 CosmicCodingMuts.normal.v99_8.tsv | sort -t$'\t' -k 8,8) > CosmicCodingMuts.normal.v99_8.sort.tsv

join -t$'\t' --header -a 1 -1 8 -2 1 -e "." -o auto CosmicCodingMuts.normal.v99_8.sort.tsv clinvar.sort.tsv > CosmicCodingMuts.normal.v99_8.sort.join.tsv

(head -n +1 CosmicCodingMuts.normal.v99_8.sort.join.tsv && tail -n +2 CosmicCodingMuts.normal.v99_8.sort.join.tsv | sort -t$'\t' -k 2,2) | cut -f 2,7,8,9,10,11 | uniq | cut -f 2,3,4,5,6 > CosmicCodingMuts.normal.v99_8.sort.join2.tsv


#シェルスクリプトを用いて、あるファイルで1列目が同一の文字列である行に対して操作を行います。
#3,4,5列目には文字列が入っています。1列目が同一の文字列である行のグループに対して、これらの行の3列目の文字列を"/"で結合した文字列を出力し、新しく6列目として各行に格納してください。同様の処理を4列目、5列目に対して行い、新しく7列目、8列目としてください。
#1列目に共通の文字列を有しない行に対しては、6,7,8列目にその行の3,4,5列の文字列を記載してください。
#最終的に6,7,8列目の列名をCLNSIG_other2, ONC_other2, CLNSIGCONF_other2としてください。

awk '
BEGIN {
  FS  = "\t"   # 入力の区切り文字 (タブ)
  OFS = "\t"   # 出力の区切り文字 (タブ)
}

NR == 1 {
  #
  # -- [1] ヘッダー行の処理 --
  #
  # 入力の1行目がヘッダーと想定し、配列に取り込む
  split($0, headerCols, FS)
  
  # 6,7,8 列目の新しい列名を設定
  headerCols[6] = "CLNSIG_other3"
  headerCols[7] = "ONC_other3"
  headerCols[8] = "CLNSIGCONF_other3"
  
  # 再構築して出力（最大で8列になる）
  max_hcol = (length(headerCols) > 8 ? length(headerCols) : 8)
  out = headerCols[1]
  for(i = 2; i <= max_hcol; i++){
    out = out OFS headerCols[i]
  }
  print out
  
  # 次行以降の処理へ
  next
}

#
# -- [2] データ行の読み込み・グルーピング --
#
{
  # 行全体を lines[NR] に保存 (後で再出力に使う)
  lines[NR] = $0
  
  # 1列目をキーとする
  key = $1
  
  # 同じ1列目をもつ行番号の一覧（カンマ区切り）を groupIndices[key] に管理
  if (key in groupIndices) {
    groupIndices[key] = groupIndices[key] "," NR
  } else {
    groupIndices[key] = NR
  }
  
  # 3列目, 4列目, 5列目を "/結合" で保持していく
  if (key in groupCol3) {
    groupCol3[key] = groupCol3[key] "/" $3
    groupCol4[key] = groupCol4[key] "/" $4
    groupCol5[key] = groupCol5[key] "/" $5
  } else {
    groupCol3[key] = $3
    groupCol4[key] = $4
    groupCol5[key] = $5
  }
}

END {
  #
  # -- [3] グループごとに行を再構築して出力 --
  #
  # groupIndices には "同じ1列目を持つ行番号リスト" が入っている
  for (key in groupIndices) {
    # カンマ区切りの行番号リストを split
    nlist = groupIndices[key]
    split(nlist, idxs, ",")  # idxs[] = {行番号1, 行番号2, ...}
    
    # グループに含まれる 3,4,5 列目を "/結合" した文字列
    merged3 = groupCol3[key]
    merged4 = groupCol4[key]
    merged5 = groupCol5[key]
    
    # グループ内の各行について再出力
    for (i = 1; i <= length(idxs); i++) {
      rowNum = idxs[i]
      # 事前に保存しておいた元の行を分割
      split(lines[rowNum], arr, FS)
      
      # 6,7,8列目として merged3, merged4, merged5 を格納
      arr[6] = merged3
      arr[7] = merged4
      arr[8] = merged5
      
      # 再構築して出力
      max_col = (length(arr) > 8 ? length(arr) : 8)
      out = arr[1]
      for (j = 2; j <= max_col; j++){
        out = out OFS arr[j]
      }
      print out
    }
  }
}' CosmicCodingMuts.normal.v99_8.sort.join2.tsv > CosmicCodingMuts.normal.v99_8.sort.join3.tsv

cut -f 1,2,6,7,8 CosmicCodingMuts.normal.v99_8.sort.join3.tsv > CosmicCodingMuts.normal.v99_8.sort.join4.tsv

(head -n +1 CosmicCodingMuts.normal.v99_8.sort.join4.tsv && tail -n +2 CosmicCodingMuts.normal.v99_8.sort.join4.tsv | sort -t$'\t' -k 1,1) | uniq | tr "\t" "," > CosmicCodingMuts.normal.v99.clinvar_sameAminoposition.csv


# CNT3は同じアミノ酸部位でのCNT総数
# CLNSIG_other3	ONC_other3 CLNSIGCONF_other3は同じアミノ酸部位でのClinVar登録



# 塩基違い

cat CosmicCodingMuts.normal.v99.vcf | sed -e '/Ter/d' | sed -e '/%3D/d' > CosmicCodingMuts.normal.v99.otherbase.tsv

awk '!/^##/ {
  # 第8列($8)を";"区切りで分解し、CNT=で始まる要素を探す
  split($8, arr, ";")
  cntVal = ""
  for(i in arr) {
    if(arr[i] ~ /^CNT=/) {
      sub(/^CNT=/, "", arr[i])  # "CNT="を削除
      cntVal = arr[i]
    }
  }
  # 1列目: chr$1:$2_$4
  # 2列目: cntVal
  print "chr"$1 ":" $2 "_" $4 "\t" cntVal
}' "$1" CosmicCodingMuts.normal.v99.otherbase.tsv > CosmicCodingMuts.normal.v99.otherbase2.tsv

awk 'BEGIN {
  FS=OFS="\t"
}
NF > 0 {
  # 空白行はスキップ(NF>0で非空行のみ処理)
  sum[$1] += $2
}
END {
  # ヘッダー行を出力
  print "variant","CNT4"
  # 1列目ごとに2列目(数値)を合計した結果を出力
  for (v in sum) {
    # 1列目(v)が空でなければ出力
    if(v != "") {
      print v, sum[v]
    }
  }
}' CosmicCodingMuts.normal.v99.otherbase2.tsv > CosmicCodingMuts.normal.v99.otherbase3.tsv

(head -n +1 CosmicCodingMuts.normal.v99.otherbase3.tsv && tail -n +2 CosmicCodingMuts.normal.v99.otherbase3.tsv | sort -t$'\t' -k 1,1) | sed '2d' | tr "\t" "," > CosmicCodingMuts.normal.v99.otherbase3.sort.csv

# CNT4は同じ塩基部位でのCNT総数



awk 'BEGIN {
  print "variant\tCLNSIG_otherbase\tONC_otherbase\tCLNSIGCONF_otherbase"
}
!/^##/ {
  split($8, arr, ";"); clnSig=""; oncVal=""; clnSigConf="";
  for(i in arr) {
    if (arr[i] ~ /^CLNSIG=/) { sub(/^CLNSIG=/,"",arr[i]); clnSig=arr[i] }
    else if (arr[i] ~ /^ONC=/) { sub(/^ONC=/,"",arr[i]); oncVal=arr[i] }
    else if (arr[i] ~ /^CLNSIGCONF=/) { sub(/^CLNSIGCONF=/,"",arr[i]); clnSigConf=arr[i] }
  }
  print "chr"$1":"$2"_"$4 "\t" clnSig "\t" oncVal "\t" clnSigConf
}' clinvar_*.vcf | sed '2d' > clinvar.otherbase.tsv


