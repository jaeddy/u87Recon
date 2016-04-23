
import urllib, urlparse, re

##url_l = 'http://www.proteinatlas.org/ENSG00000115977/cell'

##ff = urllib.urlopen(url_l)
ff = open('hpa_source3.htm','r')
cell = 'U-87+MG'
levels = []
Abs = []
for ll in ff.readlines():
    line = re.compile('\s.*cell/.*?<', re.IGNORECASE).finditer(ll)
    for part in line:
        ele = part.group()
        cell_ll = re.search('>\w.*?<',ele).group()
        cell_ll = re.sub('[\<\>]','',cell_ll)
    if cell_ll == cell.replace('+',' '):
        line = re.compile('<td nowrap class="lane".*?</td>', re.IGNORECASE).finditer(ll)
        for part in line:
            ele = part.group()
            icons = re.compile('gif', re.IGNORECASE).finditer(ele)
            level = 0
            for icon in icons:
                level += 1
            levels += [level]
    line = re.search('label">.*?<',ll)
    if line:
        Ab = re.search('\>\w.*?\<',line.group()).group()
        Ab = re.sub('[\<\>]','',Ab)
        if len(Abs) < 3:
            Abs += [Ab]
print levels
print Abs
ff.close()
