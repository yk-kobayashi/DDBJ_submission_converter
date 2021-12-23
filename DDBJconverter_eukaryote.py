# DDBJに登録するためのアノテーションファイル作成スクリプト(真核生物版、AugustusベースのGFFを想定) version 1.01 (2021.12.23)
# 入力ファイルとして、
# 1.著者や論文、生物種などの情報を記述したtsvファイル(Excelフォームに入力してtsv変換する)
# 2.ゲノムのfastaファイルをtab形式に変換したもの(元のfastaファイルにseqkit fx2tabを使用すれば良い)
# 3.アノテーションのGFFファイル(AUGUSTUS出力の結果を想定)
#    不要な行が入っているとうまく機能しないので、あらかじめ以下の作業をして必要な要素のみを含む作業用GFFに変換しておく必要がある。
#    grep -v ^"#" (GFFファイル) | grep -e "\tgene\t" -e "\ttranscript\t" -e "\tstart_codon\t" -e "\tCDS\t" -e "\tstop_codon\t"-e "\ttRNA\t" -e "\trRNA\t" > (使用するファイル)
# 4.機能アノテーション結果の表ファイル(1列目がtranscript IDと等しいことが必須。記載したいアノテーション情報の入っている列が何列目か、1.のファイルに記載する)
# の4つのファイルを使用する。
# runDDBJconvert_eukaryote.shのシェルスクリプトから起動する場合、fastaファイルとgffファイルを指定すると2-3の変換作業は自動的に実行される

import sys
import linecache

# ファイル読み込み
args = sys.argv
f1_general = args[1]
f2_seq = args[2]
f3_gff = args[3]
f4_func = args[4]

# 一般情報ファイルをディクショナリ形式で読み込む。Excelで作ったことを想定してencoding設定している
general = {}
with open(f1_general, "r", encoding="cp932") as general_info:
    for general_line in general_info:
        generallinecontent = general_line.split("\t")
        general[generallinecontent[0]] = generallinecontent[1]

# 一般情報を書き込む。Excelで入力したフォームは区切り文字の可能性のある文字を含むと""で囲まれるので、含みそうな項目は""を消す処理をする
num_authors = int(general["number_of_authors"])
print("COMMON\tSUBMITTER\t\tab_name\t"+general["submitter1"].replace("\"",""))
if num_authors > 1:
    for author in range(1, num_authors):
        authornum = "submitter" + str(author + 1)
        print("\t\t\tab_name\t" + general[authornum].replace("\"",""))
print("\t\t\tcontact\t"+general["contact_person"])
print("\t\t\temail\t"+general["email"])
print("\t\t\tphone\t"+general["phone"])
print("\t\t\tinstitute\t"+general["institute"].replace("\"",""))
print("\t\t\tdepartment\t"+general["department"].replace("\"",""))
print("\t\t\tcountry\t"+general["country"])
print("\t\t\tstate\t"+general["state"])
print("\t\t\tcity\t"+general["city"])
print("\t\t\tstreet\t"+general["street"].replace("\"",""))
print("\t\t\tzip\t"+general["zip"])
print("\tREFERENCE\t\ttitle\t"+general["reference"].replace("\"",""))
print("\t\t\tab_name\t"+general["submitter1"].replace("\"",""))
if num_authors > 1:
    for author in range(1, num_authors):
        authornum = "submitter" + str(author + 1)
        print("\t\t\tab_name\t" + general[authornum].replace("\"",""))
print("\t\t\tyear\t"+general["year"])
print("\t\t\tstatus\t"+general["status"])
print("\tST_COMMENT\t\ttagset_id\tGenome-Assembly-Data")
print("\t\t\tAssembly Method\t"+general["assembly_method"].replace("\"",""))
print("\t\t\tAssembly Name\t"+general["assembly_name"].replace("\"",""))
print("\t\t\tGenome Coverage\t"+general["genome_coverage"].replace("\"",""))
print("\t\t\tSequencing Technology\t"+general["sequencing_technology"].replace("\"",""))

# 一般情報ファイルをもとに、以下の記述に必要な情報を変数に入れておく
speciesname = general["species_name"]
name_column = int(general["name_position"]) - 1
function_column = int(general["function_position"]) - 1
locusnum = 10

# あらかじめGFFファイルとアノテーションファイルの行数を数えておく
gff_for_count = open(f3_gff)
gffcount = 0
for countgfflines in gff_for_count:
    gffcount += 1
annot_for_count = open(f4_func)
annotcount = 0
for countannotlines in annot_for_count:
    annotcount += 1

# 配列ファイルをもとに、コンティグの情報を書き込む
sequence = open(f2_seq)
for contig in sequence:
    seqcontent = contig.split()
    print(seqcontent[0] + "\tsource\t1.." + str(len(seqcontent[1])) + "\torganism\t" + speciesname)
    print("\t\t\tmol_type\tgenomic DNA")

#　GFFファイルを1行ずつ確認しながらアノテーション情報を書き込んでいく
#　なお、元のGFFファイルが、gene行　→　transcript行　→　start_codonあるいはstop_codon行　→　
#　 →　(この間にinternal, intron, initial, terminalなどがあっても良い)　→　
#　 →　CDS行　→　stop_codon行あるいはstart_codon行、の順に並んでいることが前提。
#　start_codonやstop_codon が存在するのに書き込まれていないと不完全な遺伝子扱いになるかもしれないので注意。

    for gffline in range(gffcount):
        currentline = linecache.getline(f3_gff, gffline+1)
        nextline = linecache.getline(f3_gff, gffline+2)
        gffcontent = currentline.split()
        nextcontent = nextline.split()
        if len(nextcontent) < 8:
            nextcontent = ["-", "-", "-", "-", "-", "-", "-", "-", "-"]
        direction, directstart, directend = "forward", "", ""
        if gffcontent[0] == seqcontent[0]:
            if gffcontent[6] == "-":
                direction, directstart, directend = "complement", "complement(", ")"
            if gffcontent[2] == "location":
                featuretype = "gap"
                transcriptid, locustag = "", ""
                print("\tgap\t" + gffcontent[3] + ".." + gffcontent[4] + "\testimated_length\tknown")
            if gffcontent[2] == "rRNA" or gffcontent[2] == "tRNA": #rRNAとtRNAの場合。Augustusではそもそも書き出さないので、マニュアル追加されていることになるだろうか。
                featuretype = "noncoding"
                transcriptid, locustag = "", ""
                if general["locus_tag_renaming"] == "0":
                    print("\t" + gffcontent[2] + "\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + gffcontent[8])
                if general["locus_tag_renaming"] == "1":
                    print("\t" + gffcontent[2] + "\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + general["locus_tag_prefix"] + "_" + gffcontent[8])
                if general["locus_tag_renaming"] == "2":
                    print("\t" + gffcontent[2] + "\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + general["locus_tag_prefix"] + "_" + str(locusnum).zfill(6))
                    locusnum = locusnum + 10
                print("\t\t\tproduct\t" + gffcontent[8])
            if gffcontent[2] == "gene": #遺伝子の開始。包含transcript情報を初期化。
                featuretype = "gene"
                locustag = gffcontent[8]
                transcriptid = ""
                transcriptrange, transcriptlength, transcriptnames, transcriptproducts, codonstartlist = [], [], [], [], []
            if gffcontent[2] == "transcript": #transcriptの開始。機能を確認し、transcript情報を初期化。
                transcriptid = gffcontent[8]
                genename, product = "", "hypothetical protein"
                for annotline in range(annotcount):
                    currentannot = linecache.getline(f4_func, annotline+1)
                    annotcontent = currentannot.split('\t')
                    if (annotcontent[0] == transcriptid):
                        if annotcontent[name_column] != "-" and annotcontent[name_column] != "-\n" and annotcontent[name_column] != "" and annotcontent[name_column] != "\n":
                            annotcontent[name_column] = annotcontent[name_column].replace("\n","")  #対象列が末尾の場合、改行コードが入るのを防ぐ
                            genename = annotcontent[name_column]
                        if annotcontent[function_column] != "-" and annotcontent[function_column] != "-\n" and annotcontent[name_column] != "" and annotcontent[name_column] != "\n":
                            annotcontent[function_column] = annotcontent[function_column].replace("\n","")  #対象列が末尾の場合、改行コードが入るのを防ぐ
                            annotcontent[function_column] = annotcontent[function_column].replace("  "," ")  #functional annotationの記述に連続スペースが入っている場合があるので、これも除く
                            product = annotcontent[function_column]+", putative"
                transcriptnames = transcriptnames + [genename]
                transcriptproducts = transcriptproducts + [product]
                transcriptstart, transcriptend = gffcontent[3], gffcontent[4]
                cdsrange, exonnum, exonlength, codonstart = "", 0, 0, 1
                startpos, endpos, startjoint, endjoint  = "nt", "nt", "", ""
                join, joinend = "", ""
            if (gffcontent[2] == "start_codon" or gffcontent[2] == "stop_codon") and exonnum == 0: # 5'側の開始/終止コドンを確認する
                startpos = "OK"
            if gffcontent[2] == "CDS": #CDS情報を処理する。不完全型CDSになっている可能性を考慮して確認するため処理がやや複雑。
                exonnum += 1
                exonlength = exonlength + (int(gffcontent[4]) - int(gffcontent[3]) + 1)
                if startpos == "nt" and exonnum == 1 and ((nextcontent[2] != "start_codon" and nextcontent[2] != "stop_codon") or gffcontent[3] != nextcontent[3]): #exonが短くstart/stopと最初のCDSの順番が入れ替わっている場合を考慮している
                    cdsrange ="<"
                    if direction == "forward":
                        codonstart = int(gffcontent[7]) + 1
                if transcriptid not in nextline:
                    gffcontent[4] = ">" + gffcontent[4]
                    if direction == "complement":
                        codonstart = int(gffcontent[7]) + 1
                cdsrange = cdsrange + gffcontent[3] + ".." + gffcontent[4]
                if (transcriptid in nextline) and (gffcontent[4] != transcriptend):
                    cdsrange = cdsrange + ","
            if (transcriptid not in nextline) and (featuretype == "gene"): #transcriptの終わりを判定し、各種情報をストック
                if exonnum > 1:
                    join, joinend = "join(", ")"
                transcriptrange = transcriptrange + [directstart + join + cdsrange + joinend + directend]
                transcriptlength = transcriptlength + [exonlength]
                codonstartlist = codonstartlist + [codonstart]
            if (locustag not in nextline) and (featuretype == "gene"): #geneの終わりを判定し、最長トランスクリプト情報を出力
                maxtranscript = transcriptlength.index(max(transcriptlength))
                if general["locus_tag_renaming"] == "0":
                    print("\tCDS\t" + transcriptrange[maxtranscript] + "\tlocus_tag\t" + locustag)
                if general["locus_tag_renaming"] == "1":
                    print("\tCDS\t" + transcriptrange[maxtranscript] + "\tlocus_tag\t" + general["locus_tag_prefix"] + "_" + locustag)
                if general["locus_tag_renaming"] == "2":
                    print("\tCDS\t" + transcriptrange[maxtranscript] + "\tlocus_tag\t" + general["locus_tag_prefix"] + "_" + str(locusnum).zfill(6))
                    locusnum = locusnum + 10
                if transcriptnames[maxtranscript] != "":
                    print("\t\t\tgene\t"+transcriptnames[maxtranscript])
                print("\t\t\tproduct\t"+transcriptproducts[maxtranscript])
                print("\t\t\ttransl_table\t" + general["transl_table"])
                print("\t\t\tcodon_start\t"+str(codonstartlist[maxtranscript]))



