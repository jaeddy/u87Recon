
import itertools, re, time

filename = raw_input('Input HPA file name\n')

NCBIfile = open('Homo_sapiens.gene_info', 'r')

t_0 = time.clock()

ensembl2gene = {}
for line in NCBIfile.readlines():
    cols = line.split('\t')
    geneID = cols[1]
    ensemblID = re.search('ENSG[0-9]+', cols[5])
    if ensemblID:
        ensemblID = ensemblID.group()
        ensembl2gene[ensemblID] = geneID
NCBIfile.close()

HPAfile = open(filename, 'r')
genesfilename = re.sub('.txt', '_genes.txt', filename)
genesfile = open(genesfilename, 'w')

for line in HPAfile.readlines():
    cols = line.split('\t')
    ensemblID = cols[0]
    ensemblID = re.sub('\n', '', ensemblID)
    if ensemblID in ensembl2gene:
        genesfile.write(ensembl2gene[ensemblID])
        genesfile.write('\n')
    else:
        genesfile.write(ensemblID)
        genesfile.write('\n')
print time.clock() - t_0
genesfile.close()
