# DDBJに登録するためのアノテーションファイル作成スクリプト(原核生物版、ProkkaベースのGFFを想定) version 2.0 (2023.4.14)
# 入力ファイルとして、
# 1.著者や論文、生物種などの情報を記述したtsvファイル(Excelフォームに入力してtsv変換する)
# 2.ゲノムのfastaファイルをtab形式に変換したもの(元のfastaファイルにseqkit fx2tabを使用すれば良い)
#   runDDBJconvert_prokaryote.shのシェルスクリプトから実行する場合、元のfastaから自動的にタブ区切りファイルが生成される
# 3.アノテーションのGFFファイル(Prokka出力の結果を想定)
#    CDSやrRNA,tRNAといった書き込む要素の行のみ使用する(geneやmRNA行は使用しない)
# 4.機能アノテーション結果の表ファイル(1列目がlocus tagと等しいことが必須。記載したい遺伝子名とアノテーション情報の入っている列が何列目か、1.のファイルに記載する)
# の4つのファイルを使用する。

import sys
import linecache

# ファイル読み込み
args = sys.argv
f1_general = args[1]
f2_seq = args[2]
f3_gff = args[3]
f4_func = args[4]

# 一般情報ファイルをディクショナリ形式で読み込む
general = {}
with open(f1_general, "r", encoding="cp932") as general_info:
    for general_line in general_info:
        generallinecontent = general_line.split("\t")
        general[generallinecontent[0]] = generallinecontent[1]

# 一般情報を書き込む。Excelで入力したフォームは区切り文字の可能性のある文字を含むと""で囲まれるので、含みそうな項目は""を消す処理をする
print("COMMON\tDATATYPE\t\ttype\tWGS")
print("\tKEYWORD\t\tkeyword\tWGS")
print("\t\t\tkeyword\t"+general["dataquality"])
print("\tDBLINK\t\tproject\t"+general["bioprojectID"])
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

# 一般情報ファイルをもとに、以下の記述に必要な情報を変数に入れておく
speciesname = general["species_name"]
strain = general["strain_or_isolate"].split("=")
name_column = int(general["name_position"]) - 1
function_column = int(general["function_position"]) - 1
locusnum = 10
contignameprefix = ""
if general["contigname_numonly"] == "yes":
    contignameprefix = "contig"

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
    print(contignameprefix + seqcontent[0] + "\tsource\t1.." + str(len(seqcontent[1])) + "\torganism\t" + speciesname)
    print("\t\t\t" + strain[0] + "\t" + strain[1])
    print("\t\t\tcountry\t" + general["geo_loc"])
    print("\t\t\tcollection_date\t" + general["collection_date"])
    print("\t\t\tmol_type\tgenomic DNA")
    print("\t\t\tsubmitter_seqid\t@@[entry]@@")
    print("\t\t\tff_definition\t@@[organism]@@ @@[" + strain[0] + "]@@ DNA, @@[submitter_seqid]@@")

#　GFFファイルからアノテーション情報を書き込んでいく
#　CDSやrRNA,tRNAといった書き込む要素の行のみ使用する

    for gffline in range(gffcount):
        currentline = linecache.getline(f3_gff, gffline+1)
        gffcontent = currentline.split()
        direction, directstart, directend = "forward", "", ""
        if gffcontent[0] == seqcontent[0]:
            if gffcontent[6] == "-":
                direction, directstart, directend = "complement", "complement(", ")"
            if gffcontent[2] == "location": #アセンブリギャップの情報を処理
                print("\tassembly_gap\t" + gffcontent[3] + ".." + gffcontent[4] + "\testimated_length\tknown")
                print("\t\t\tgap_type\twithin scaffold")
                print("\t\t\tlinkage_evidence\tunspecified")
            if gffcontent[2] == "rRNA" or gffcontent[2] == "tRNA"
                description = gffcontent[8].split(";")
                descriptionnum = len(description)
                detail = {}
                for category in range(descriptionnum): # GFFの説明カラムにある情報を読み込む
                    eachcategory = description[category].split("=")
                    eachcategory[1] = eachcategory[1].replace("\n", "")
                    detail[eachcategory[0]] = eachcategory[1]
                if general["locus_tag_renaming"] == "0":
                    print("\t" + gffcontent[2] + "\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + detail["ID"])
                if general["locus_tag_renaming"] == "1":
                    print("\t" + gffcontent[2] + "\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + general["locus_tag_prefix"] + "_" + detail["ID"])
                if general["locus_tag_renaming"] == "2":
                    print("\t" + gffcontent[2] + "\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + general["locus_tag_prefix"] + "_" + str(locusnum).zfill(6))
                    locusnum = locusnum + 10
                print("\t\t\tproduct\t" + detail["product"])
            if gffcontent[2] == "CDS":
                description = gffcontent[8].split(";")
                descriptionnum = len(description)
                detail = {}
                for category in range(descriptionnum): # GFFの説明カラムにある情報を読み込む
                    eachcategory = description[category].split("=")
                    eachcategory[1] = eachcategory[1].replace("\n", "")
                    detail[eachcategory[0]] = eachcategory[1]
                if general["locus_tag_renaming"] == "0":
                    print("\tCDS\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + detail["ID"])
                if general["locus_tag_renaming"] == "1":
                    print("\tCDS\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + general["locus_tag_prefix"] + "_" + detail["ID"])
                if general["locus_tag_renaming"] == "2":
                    print("\tCDS\t" + directstart + gffcontent[3] + ".." + gffcontent[4] + directend + "\tlocus_tag\t" + general["locus_tag_prefix"] + "_" + str(locusnum).zfill(6))
                    locusnum = locusnum + 10
                gene, product, genenote = "", "hypothetical protein", ""
                if "gene" in detail:
                    gene = detail["gene"]
                for annotline in range(annotcount): # アノテーションファイルから翻訳産物の情報を抽出する
                    currentannot = linecache.getline(f4_func, annotline+1)
                    annotcontent = currentannot.split('\t')
                    if (annotcontent[0] == detail["ID"]):
                        if annotcontent[name_column] != "-" and annotcontent[name_column] != "-\n" and annotcontent[name_column] != "" and annotcontent[name_column] != "\n":
                            annotcontent[name_column] = annotcontent[name_column].replace("\n","")  #対象列が末尾の場合、改行コードが入るのを防ぐ
                            gene = annotcontent[name_column]
                        if annotcontent[function_column] != "-" and annotcontent[function_column] != "-\n" and annotcontent[function_column] != "" and annotcontent[function_column] != "\n":
                            annotfunction = annotcontent[function_column].replace("\n","")  #対象列が末尾の場合、改行コードが入るのを防ぐ
                            annotfunction = annotfunction.replace("  "," ")  #functional annotationの記述に連続スペースが入っている場合があるので、これも除く
                            if annotfunction[0].isupper() and annotfunction[1].islower(): #略語等以外の先頭大文字を小文字に変換する。略号等かどうかは2文字めが英字小文字かどうかで判定
                                annotfunction = annotfunction[0].lower() + annotfunction[1:]
                            #機能アノテーションが文章だった場合はnoteに、そうでなければproductにする
                            if (annotfunction.startswith("it ") or annotfunction.startswith("the ") or annotfunction.startswith("this ") or annotfunction.startswith("to ")) == True:
                                genenote = annotfunction
                            elif ("functions " in annotfunction) or ("acts " in annotfunction) or ("plays " in annotfunction) or ("catalyzes " in annotfunction) or ("hydrolyzes " in annotfunction) or ("demethylases " in annotfunction) or ("promotes " in annotfunction) or ("splits " in annotfunction) or ("adds " in annotfunction) or ("has " in annotfunction) or ("assists " in annotfunction) or ("mediates " in annotfunction) or ("cleaves " in annotfunction) or ("produces " in annotfunction) or ("stimulates " in annotfunction) or ("allows " in annotfunction) or ("occurs " in annotfunction) or ("shuttles " in annotfunction) or ("cooperates " in annotfunction) or ("methylates " in annotfunction) or ("condensation" in annotfunction) or ("constitutes " in annotfunction) or ("may " in annotfunction) or ("might " in annotfunction) or ("seems to " in annotfunction) or ("binds to " in annotfunction) or ("present " in annotfunction) or ("control of " in annotfunction) or ("regulation " in annotfunction) or ("probably " in annotfunction) or ("annotation" in annotfunction):
                                genenote = annotfunction
                            elif annotfunction.startswith("putative ") == True:
                                product = annotfunction
                            else:
                                product = "putative " + annotfunction
                            #機能アノテーションが説明を含んでいる場合、説明部分はnoteに移す
                            if "It is " in product:
                                productdetail = product.split("It is ")
                                product = productdetail[0]
                                genenote = productdetail[1] + ". " + genenote
                            if "that is " in product:
                                productdetail = product.split("that is ")
                                product = productdetail[0]
                                genenote = productdetail[1] + ". " + genenote
                            if "which is " in product:
                                productdetail = product.split("which is ")
                                product = productdetail[0]
                                genenote = productdetail[1] + ". " + genenote
                            if "predicted to be " in product:
                                productdetail = product.split("predicted to be ")
                                product = productdetail[0]
                                genenote = "predicted to be " + productdetail[1] + ". " + genenote
                            if "belongs to " in product:
                                productdetail = product.split("belongs to ")
                                product = productdetail[0]
                                genenote = "belongs to " + productdetail[1] + ". " + genenote
                            if "functions in " in product:
                                productdetail = product.split("functions in ")
                                product = productdetail[0]
                                genenote = "functions in " + productdetail[1] + ". " + genenote
                            if "involved in " in product:
                                productdetail = product.split("involved in ")
                                product = productdetail[0]
                                genenote = "involved in " + productdetail[1] + ". " + genenote
                            if "required for " in product:
                                productdetail = product.split("required for ")
                                product = productdetail[0]
                                genenote = "required for " + productdetail[1] + ". " + genenote
                            if "participates in " in product:
                                productdetail = product.split("participates in ")
                                product = productdetail[0]
                                genenote = "participates in " + productdetail[1] + ". " + genenote
                            if "implicated in " in product:
                                productdetail = product.split("implicated in ")
                                product = productdetail[0]
                                genenote = "implicated in " + productdetail[1] + ". " + genenote
                            if " which " in product:
                                productdetail = product.split(" which ")
                                product = productdetail[0]
                                genenote = productdetail[1] + ". " + genenote
                            if " that " in product:
                                productdetail = product.split(" that ")
                                product = productdetail[0]
                                genenote = productdetail[1] + ". " + genenote
                            if ". " in product:
                                productdetail = product.split(". ")
                                product = productdetail[0]
                                genenote = productdetail[1] + ". " + genenote
                            if (", " in product) and (", and " not in product):
                                productdetail = product.split(", ")
                                product = productdetail[0]
                                genenote = productdetail[1] + ". " + genenote
                            if " are " in product:
                                productdetail = product.split(" are ")
                                product = productdetail[0]
                                genenote = productdetail[1] + ". " + genenote
                            if " is " in product:
                                productdetail = product.split(" is ")
                                product = productdetail[0]
                                genenote = productdetail[1] + ". " + genenote
                            #product末尾の調整。文末の句読点の除去、ドメイン・モチーフやファミリーのタンパク質productとしての適正化
                            if product.endswith(" ") == True:
                                product = product.rstrip()
                            if product.endswith(".") == True:
                                product = product.rstrip(".")
                            if product.endswith(",") == True:
                                product = product.rstrip(",")
                            if product.endswith("activity") == True:
                                product = product.rstrip("activity")
                            if (product.endswith("proteins") == True) or (product.endswith("homologues") == True):
                                product = product.rstrip("s")
                            if (product.endswith("domain") == True) or (product.endswith("domains") == True) or (product.startswith("putative domain") == True) or ("domain of unknown function" in product):
                                product = product + " containing protein"
                            if (product.endswith("motif") == True) or (product.endswith("motifs") == True):
                                product = product + " containing protein"
                            if (product.endswith("repeat") == True) or (product.endswith("repeats") == True):
                                product = product + " containing protein"
                            if (product.endswith("family") == True) or (product.endswith("like") == True) or (product.endswith("binding") == True) or (product.endswith("finger") == True) or (product.endswith("zipper") == True) or (product.endswith("nuckle") == True) or (product.endswith("membrane") == True) or (product.endswith("carrier") == True):
                                product = product + " protein"
                            #説明をnoteに移した結果として内容がなくなったproductはhypothetical proteinに戻す
                            if product == "putative":
                                product = "hypothetical protein"
                if gene != "":
                    print("\t\t\tgene\t"+gene)
                print("\t\t\tproduct\t"+product)
                if genenote != "":
                    print("\t\t\tnote\t" + genenote)
                print("\t\t\ttransl_table\t" + general["transl_table"])
                print("\t\t\tcodon_start\t"+str(int(gffcontent[7]) + 1))
