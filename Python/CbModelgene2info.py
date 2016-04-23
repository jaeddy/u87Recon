
import itertools, re, time

##filename = raw_input('Input GEO file name\n')

genesfilename = raw_input('Input COBRA genes file name\n')

NCBIfile = open('Homo_sapiens.gene_info', 'r')

t_0 = time.clock()

gene2info = {}
count = 0
for line in NCBIfile.readlines():
    cols = line.split('\t')
    if count==0:
        print len(cols)
    geneID = cols[1]
    info = line
    gene2info[geneID] = info
NCBIfile.close()

Cbgenesfile = open(genesfilename, 'r')
Cbinfofilename = re.sub('genes.txt', 'info.txt', genesfilename)
Cbinfofile = open(Cbinfofilename, 'w')

nomatch = 0
count = 0
for line in Cbgenesfile.readlines():
    genes = line.split('\t')
    for gene in genes:
        if gene.isdigit():
            if gene in gene2info:
                Cbinfofile.write(gene2info[gene])
                Cbinfofile.write('\t')
            else:
                Cbinfofile.write(gene)
                Cbinfofile.write('\t')
                nomatch += 1

        elif gene.find('\n')!=0:
            gene = gene[0:gene.find('\n')]
            if gene in gene2info:
                Cbinfofile.write(gene2info[gene])
                Cbinfofile.write('\n')
            else:
                Cbinfofile.write(gene)
                Cbinfofile.write('\n')
                nomatch += 1

        else:
            Cbinfofile.write(gene)
            Cbinfofile.write('\t')

##    count += 1
##    if count > 2:
##        break    
Cbinfofile.close()

print nomatch
print time.clock() - t_0
