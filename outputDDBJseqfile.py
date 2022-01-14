import sys

args = sys.argv
generalfile = args[1]
genomeseq = args[2]

# コンティグ名が数値だけの場合の改名に対応するためだけのために、一般情報ファイルを読み込む。使うのは１箇所だけ。
# Excelで作ったことを想定してencoding設定している
general = {}
with open(generalfile, "r", encoding="cp932") as general_info:
    for general_line in general_info:
        generallinecontent = general_line.split("\t")
        general[generallinecontent[0]] = generallinecontent[1]

contignameprefix = ""
if general["contigname_numonly"] == "yes":
    contignameprefix = "contig"

sequence = open(genomeseq)
for contig in sequence:
    seqcontent = contig.split()
    print(">" + contignameprefix + seqcontent[0])
    print(seqcontent[1])
    print("//")


