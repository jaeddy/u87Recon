
import itertools, re, time

##filename = raw_input('Input GEO file name\n')
filename = 'U-87_MG_GEO.txt'

NCBIfile = open('Homo_sapiens.gene_info', 'r')

t_0 = time.clock()

symbol2gene = {}
synonym2gene = {}
for line in NCBIfile.readlines():
    cols = line.split('\t')
    geneID = cols[1]
    symbol = cols[2]
    symbol2gene[symbol] = geneID
    synonyms = cols[4].split('|')
    for synonym in synonyms:
        synonym2gene[synonym] = geneID
NCBIfile.close()

GEOfile = open(filename, 'r')
genesfilename = re.sub('.txt', '_genes.txt', filename)
genesfile = open(genesfilename, 'w')

nomatch = 0
for line in GEOfile.readlines():
    symbol = re.sub('\n', '', line)
    if symbol in symbol2gene:
        #print symbol2gene[symbol]
        genesfile.write(symbol2gene[symbol])
        genesfile.write('\n')
    elif symbol in synonym2gene:
        genesfile.write(synonym2gene[symbol])
        genesfile.write('\n')
    else:
        #print symbol
        genesfile.write(symbol)
        genesfile.write('\n')
        nomatch += 1
genesfile.close()

print nomatch
print time.clock() - t_0
