# DDBJに登録するためのアノテーションファイル作成スクリプト(ミトコンドリア版、MITOSベースのGFFを想定)
# 入力ファイルとして、
# 1.著者や論文、生物種などの情報を記述したtsvファイル(Excelフォームに入力してtsv変換する)
# 2.ゲノムのfastaファイルをtab形式に変換したもの(元のfastaファイルにseqkit fx2tabを使用すれば良い)
# 3.アノテーションのGFFファイル(MITOS出力の結果を想定)
#    CDSやrRNA,tRNAといった書き込む要素の行のみ使用する
# の4つのファイルを使用する。

import sys
import linecache

# ファイル読み込み
args = sys.argv
f1_general = args[1]
f2_seq = args[2]
f3_gff = args[3]

# 一般情報ファイルをディクショナリ形式で読み込む
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

# 一般情報ファイルをもとに、以下の記述に必要な情報を変数に入れておく。
# MITOSはlocus tagを設定しないので、locus tagは新規に設定する前提。
speciesname = general["species_name"]
locustagprefix = general["locus_tag_prefix"]
if general["locus_tag_shared"] == "1":
    locustagprefix = locustagprefix + "_mit"
locusnum = 10

# あらかじめGFFファイルの行数を数えておく
gff_for_count = open(f3_gff)
gffcount = 0
for countgfflines in gff_for_count:
    gffcount += 1

# 配列ファイルをもとに、コンティグの情報を書き込む
sequence = open(f2_seq)
for contig in sequence:
    seqcontent = contig.split()
    print(seqcontent[0] + "\tsource\t1.." + str(len(seqcontent[1])) + "\torganism\t" + speciesname)
    print("\t\t\tmol_type\tgenomic DNA")

#　GFFファイルからアノテーション情報を書き込んでいく
#　CDSやrRNA,tRNAといった書き込む要素の行のみ使用する

    for gffline in range(gffcount):
        currentline = linecache.getline(f3_gff, gffline+1)
        gffcontent = currentline.split()
        direction, directstart, directend = "forward", "", ""
        if gffcontent[0] == seqcontent[0]:
            description = gffcontent[8].split("=")
            if gffcontent[6] == "-":
                direction, directstart, directend = "complement", "complement(", ")"
            if gffcontent[2] == "location": #アセンブリギャップの情報を処理
                print("\tgap\t" + gffcontent[3] + ".." + gffcontent[4] + "\testimated_length\tknown")
            if gffcontent[2] == "rep_origin":
                print("\t" + gffcontent[2] + "\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tgene\t" + description[1])
            if gffcontent[2] == "rRNA" or gffcontent[2] == "tRNA":
                print("\t" + gffcontent[2] + "\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + locustagprefix + str(locusnum).zfill(3))
                locusnum = locusnum + 10
                print("\t\t\tgene\t" + description[1])
            if gffcontent[2] == "gene":
                print("\tCDS\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + locustagprefix + str(locusnum).zfill(3))
                locusnum = locusnum + 10
                print("\t\t\tgene\t" + description[1])
                print("\t\t\ttransl_table\t" + general["transl_table"])
                print("\t\t\tcodon_start\t1")
