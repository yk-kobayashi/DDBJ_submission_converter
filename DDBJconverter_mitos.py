# DDBJに登録するためのアノテーションファイル作成スクリプト(ミトコンドリア版、MITOSベースのGFFを想定) version 2.0 (2023.4.14)
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
print("COMMON\tDBLINK\t\tproject\t"+general["bioprojectID"])
print("\t\t\tbiosample\t"+general["biosampleID"])
print("\t\t\tsequence read archive\t"+general["sequenceReadArchive"])
num_authors = int(general["number_of_authors"])
print("\tSUBMITTER\t\tab_name\t"+general["submitter1"].replace("\"",""))
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
print("\tDATE\t\thold_date\t"+general["hold_date"])
print("\tST_COMMENT\t\ttagset_id\tGenome-Assembly-Data")
print("\t\t\tAssembly Method\t"+general["assembly_method"].replace("\"",""))
print("\t\t\tAssembly Name\t"+general["assembly_name"].replace("\"",""))
print("\t\t\tGenome Coverage\t"+general["genome_coverage"].replace("\"",""))
print("\t\t\tSequencing Technology\t"+general["sequencing_technology"].replace("\"",""))

# 一般情報ファイルをもとに、以下の記述に必要な情報を変数に入れておく。
# MITOSはlocus tagを設定しないので、locus tagは新規に設定する前提。
speciesname = general["species_name"]
strain = general["strain_or_isolate"].split("=")
locustagprefix = general["locus_tag_prefix"]
if general["locus_tag_shared"] == "1":
    locustagprefix = locustagprefix + "_mit"
locusnum = 10
contignameprefix = ""
if general["contigname_numonly"] == "yes":
    contignameprefix = "contig"

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
    print("\t\t\t" + strain[0] + "\t" + strain[1])
    print("\t\t\tcountry\t" + general["geo_loc"])
    print("\t\t\tcollection_date\t" + general["collection_date"])
    print("\t\t\torganelle\tmitochondrion")
    print("\t\t\tmol_type\tgenomic DNA")
    if general["circularseq"] == "yes":
        print("\tTOPOLOGY\t\tcircular\t")

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
#            if gffcontent[2] == "rep_origin": # replication originが正しくアノテーションされていれば記入したいが、mitos結果はOHの断片がそれぞれrep_originとして出てきてしまうので入れない
#                print("\t" + gffcontent[2] + "\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tgene\t" + description[1])
            if gffcontent[2] == "rRNA" or gffcontent[2] == "tRNA":
                print("\t" + gffcontent[2] + "\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + locustagprefix + str(locusnum).zfill(3))
                locusnum = locusnum + 10
                print("\t\t\tgene\t" + description[1])
            if gffcontent[2] == "gene":
                print("\tCDS\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + locustagprefix + str(locusnum).zfill(3))
                locusnum = locusnum + 10
                print("\t\t\tgene\t" + description[1])
                if description[1] == "atp6":
                    print("\t\t\tproduct\tATP synthase F0 subunit 6")
                if description[1] == "atp8":
                    print("\t\t\tproduct\tATP synthase F0 subunit 8")
                if description[1] == "cob" or description[1] == "cytB":
                    print("\t\t\tproduct\tcytochrome b")
                if description[1] == "cox1":
                    print("\t\t\tproduct\tcytochrome c oxidase subunit 1")
                if description[1] == "cox2":
                    print("\t\t\tproduct\tcytochrome c oxidase subunit 2")
                if description[1] == "cox3":
                    print("\t\t\tproduct\tcytochrome c oxidase subunit 3")
                if description[1] == "nad1":
                    print("\t\t\tproduct\tNADH dehydrogenase subunit 1")
                if description[1] == "nad2":
                    print("\t\t\tproduct\tNADH dehydrogenase subunit 2")
                if description[1] == "nad3":
                    print("\t\t\tproduct\tNADH dehydrogenase subunit 3")
                if description[1] == "nad4":
                    print("\t\t\tproduct\tNADH dehydrogenase subunit 4")
                if description[1] == "nad4l":
                    print("\t\t\tproduct\tNADH dehydrogenase subunit 4L")
                if description[1] == "nad5":
                    print("\t\t\tproduct\tNADH dehydrogenase subunit 5")
                if description[1] == "nad6":
                    print("\t\t\tproduct\tNADH dehydrogenase subunit 6")
                print("\t\t\ttransl_table\t" + general["transl_table"])
                print("\t\t\tcodon_start\t1")
