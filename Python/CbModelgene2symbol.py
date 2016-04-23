
import itertools, re, time

##filename = raw_input('Input GEO file name\n')

genesfilename = raw_input('Input COBRA genes file name\n')

NCBIfile = open('Homo_sapiens.gene_info', 'r')

t_0 = time.clock()

gene2symbol = {}
gene2synonym = {}
for line in NCBIfile.readlines():
    cols = line.split('\t')
    geneID = cols[1]
    symbol = cols[2]
    gene2symbol[geneID] = symbol
NCBIfile.close()

Cbgenesfile = open(genesfilename, 'r')
Cbsymbolsfilename = re.sub('genes.txt', 'symbols.txt', genesfilename)
Cbsymbolsfile = open(Cbsymbolsfilename, 'w')

nomatch = 0
count = 0
for line in Cbgenesfile.readlines():
    genes = line.split('\t')
    for gene in genes:
        if gene.isdigit():
            if gene in gene2symbol:
                Cbsymbolsfile.write(gene2symbol[gene])
                Cbsymbolsfile.write('\t')
            else:
                Cbsymbolsfile.write(gene)
                Cbsymbolsfile.write('\t')
                nomatch += 1

        elif gene.find('\n')!=0:
            gene = gene[0:gene.find('\n')]
            if gene in gene2symbol:
                Cbsymbolsfile.write(gene2symbol[gene])
                Cbsymbolsfile.write('\n')
            else:
                Cbsymbolsfile.write(gene)
                Cbsymbolsfile.write('\n')
                nomatch += 1

        else:
            Cbsymbolsfile.write(gene)
            Cbsymbolsfile.write('\t')

##    count += 1
##    if count > 2:
##        break    
Cbsymbolsfile.close()

print nomatch
print time.clock() - t_0
