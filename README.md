# DDBJ_submission_converter
Scripts to convert genome and annotation data into DDBJ submission files  
ゲノム塩基配列データと遺伝子アノテーションのデータから、DDBJに登録する形式のファイルに変換するためのスクリプト群。  
真核生物用、原核生物用、ミトコンドリア用に分かれています。

入力元ファイルとして、  
・ゲノムの塩基配列が記されたfastaファイル  
・遺伝子領域のアノテーションされたgffファイル  
・遺伝子機能をアノテーションしたcsvファイル(EggNOG mapperの結果のtableそのままでOK。ミトコンドリアゲノムの場合は不要)  
と、基本情報をExcelで記入したtsvファイルを使用し、  
DDBJの登録に必要な配列ファイルとアノテーションファイルを出力します。

遺伝子アノテーションは  
・真核生物はAugustus系列(BRAKERもAugustus依存)  
・原核生物はProkka  
・ミトコンドリアゲノムはMITOS  
でgffファイルを作成したと想定しています。

○依存ツール  
・Python3  
・seqkit  

○使用法  
(1) general_templateのtsvファイル(本体ゲノム用orミトコンドリア用)を編集し、登録者や生物種、解析手法などの情報を記入

(2) シェルスクリプトを起動  
　　2-a) 真核生物ゲノムの場合、  
　　./runDDBJconvert_eukaryote.sh  <(1)の基本情報ファイル> <塩基配列.fasta> <gffファイル(Augustus準拠)> <機能アノテーション.csv>

　　2-b) 原核生物ゲノムの場合、  
　　./runDDBJconvert_prokaryote.sh  <(1)の基本情報ファイル> <塩基配列.fasta> <gffファイル(Prokka準拠)> <機能アノテーション.csv>

　　2-c) ミトコンドリアゲノムの場合、  
　　./runDDBJconvert_mitos.sh  <(1)の基本情報ファイル> <塩基配列.fasta> <gffファイル(MITOS準拠)>

(3) 出力結果  
DDBJに登録する際の配列ファイル(元のゲノムのfastaファイルのエントリ間に//の文字を挟んだだけのもの)が  
DDBJseqfile.txt  
アノテーションファイル(形式の独特な登録用ファイル)が  
DDBJannotfile.txt  
というファイル名で出力されます。

最初に起動したシェルスクリプトが作業用の中途ファイルをDDBJfiletempというフォルダ内に作成し、これを使ってpythonスクリプトを走らせる流れです。  
最終的なDDBJannotfile.txtとDDBJseqfile.txtがうまく生成されていれば、DDBJfiletempのフォルダは削除してしまって大丈夫です。  

※注  
作成時の実行環境ではUME/jParserでエラーの出ないアノテーションファイルを作成できていましたが、登録前には自身のファイルをUME/jParserでチェックしておく必要があるかと思います。
なお、当然ながら元データに問題がある場合はどうにもならないので、その場合は元データを作り直すか、作成したデータをマニュアル補正する必要があります。  
(gapが遺伝子領域に被っていたり、CDS内部にstopコドンが存在していたりする場合など)
