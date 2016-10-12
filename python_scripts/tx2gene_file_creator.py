from sys import argv


def main(path, fasta):

    fasta = path + fasta
    output = path + "tx2gene_rrna.csv"
    fasta = open(fasta, 'r')
    output = open(output, 'w')
    output.write("TXNAME,GENEID" + "\n")

    for line in fasta:
        if line[0] == ">":
            output.write(line[1:line.index(" ")])
            output.write(",")
            output.write(line[line.index("gene:") + 5: line.index(" ", line.index("gene:"))] + "\n")

    fasta.close()
    output.close()

if __name__ == "__main__":
    main(argv[1], argv[2])
