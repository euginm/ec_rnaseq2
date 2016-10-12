goa = open("/home/euginm/Documents/ref/18.E_coli_MG1655.goa", 'r')
output = open("/home/euginm/Documents/ref/18.E_coli_MG1655_extracted.goa", 'w')

goa.next()
goa.next()

for line in goa:
    gene_name = ''
    for name in line.split('\t')[10].split('|'):
        if name[0] == 'b' and len(name) == 5:
            gene_name = name
    output.write(gene_name + '\t' + 'GO:' + line.split('GO:')[1][:7] + '\n')

goa.close()
output.close()
