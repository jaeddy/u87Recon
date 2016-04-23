######################################################################
#### created by Younhee Ko  | 02.05.2011 | younhee.ko@gmail.com   ####
#### modified by James Eddy | 02.16.2011 | james.a.eddy@gmail.com ####
######################################################################

# parse the human atlas data and get result!

import urllib, urlparse, re, time

# Allow user to select cell line
options = 'Select cell line:\n\
           1 U-87 MG (glioblastoma)\n\
           2 U-138 MG (glioblastoma)\n\
           3 U-251 MG (glioblastoma)\n\n'
option = int(raw_input(options))-1

cells = ['U-87+MG','U-138+MG','U-251+MG']
cell = cells[option]

filename = re.sub('\+','_',cell) + '.txt'
file = open(filename,'w')

url_0 = 'http://www.proteinatlas.org/search/cell_expression:{c};'.format(c=cell)
t_0 = time.clock()

# Collect strongly and moderately expressed proteins for the cell line
strengths = ['Strong','Moderate']
for strength in strengths:
    print 'Collecting {s}ly expressed proteins...\n'.format(s=strength.lower())
    file.write('{s}\n'.format(s=strength))
    url_s = url_0 + strength

    f = urllib.urlopen(url_s)
    for l in f.readlines():
        # Determine total number of pages of proteins
        lines = re.compile('of.*?</p>', re.IGNORECASE).finditer(l)
        for line in lines:
            line = line.group()
            parts = re.compile('[0-9]+\s', re.IGNORECASE).finditer(line)
            for part in parts:
                pages = int(part.group())
        # Determine total number of proteins found
        line = re.search('<p>.*?GENES FOUND</p>', l)
        if line:
            line = line.group()
            genes = re.search('>[0-9].*?\s', line).group()
            genes = int(re.sub('[\>\s]', '', genes))
    f.close()

    # Collect protein data from each page
    add = 0; discard = 0
    for p in range(1,pages+1):
        print 'Page {p} of {total}'.format(p=str(p),total=str(pages))
        url = url_s + '&page=' + str(p)
        f = urllib.urlopen(url)
        for l in f.readlines():
            lines = re.compile('<td ><b><a href="/ENSG.*?tr>', re.IGNORECASE).finditer(l)
            for line in lines:
                line = line.group()
                parts = re.compile('<td.*?</td>', re.IGNORECASE).finditer(line)
                index = 1
                for part in parts:
                    part = part.group()
                    if index == 1:
                        # Collect gene ID
                        pos1 = part.find('"/')
                        pos2 = part.find('">', pos1)
                        ID = part[pos1+2:pos2]
                        # Collect gene name
                        pos1 = part.find('">')
                        pos2 = part.find('<', pos1)
                        name = part[pos1+2:pos2]
                    elif index == 2:
                        # Collect gene synonyms
                        pos1 = part.find('>')
                        pos2 = part.find('</td>',pos1)
                        synonyms = part[pos1+1:pos2]
                        synonyms = synonyms.replace('<br>', ',')
                    elif  index == 3:
                        # Collect gene description
                        pos1 = part.find('>')
                        pos2 = part.find('<', pos1)
                        desc = part[pos1+1:pos2]
                    elif  index == 4:
                        # Collect antibodies
                        Ab_rows = re.compile('\w*?<br>', re.IGNORECASE).finditer(part)
                        Abs = []
                        for Ab_row in Ab_rows:
                            Ab_row = Ab_row.group()
                            Ab = Ab_row.replace('<br>', '')
                            Abs += [Ab]
                    elif index == 5:
                        # Collect validation scores
                        val_rows = re.split('<br>', part)
                        Ab_scores = []
                        for val_row in val_rows:
                            scores = re.compile('/validation.*?.png', re.IGNORECASE).finditer(val_row)
                            Ab_score = 0
                            for score in scores:
                                score = score.group()
                                green = re.search('green', score)
                                if green:
                                    Ab_score += 1
                            Ab_scores += [Ab_score]
                    index += 1
                # Determine the number of valid antibodies
                num_valid = 0
                for i in range(0, len(Ab_scores)):
                    if Ab_scores[i] >= 2:
                        num_valid += 1
                # Determine which antibody or antibodies are measured in the cell line (if multiple are listed)
                if (num_valid > 0) & (num_valid < len(Abs)):
                    # Secondary page reporting measured expression for antibodies in all cell lines
                    url_l = 'http://www.proteinatlas.org/{ID}/cell'.format(ID=ID)
                    f_l = urllib.urlopen(url_l)
                    X_cell = ''; X_counts = []; X_Abs = []
                    for ll in f_l.readlines():
                        # Locate expression levels for the current cell line
                        lines_l = re.compile('\s.*cell/.*?<', re.IGNORECASE).finditer(ll)
                        for line_l in lines_l:
                            line_l = line_l.group()
                            X_cell = re.search('>\w.*?<', line_l).group()
                            X_cell = re.sub('[\<\>]', '', X_cell)
                        # Count cell line replicates in which each antibody is found
                        if X_cell == cell.replace('+', ' '):
                            lines_l = re.compile('<td nowrap class="lane".*?</td>', re.IGNORECASE).finditer(ll)
                            for line_l in lines_l:
                                line_l = line_l.group()
                                X_icons = re.compile('gif', re.IGNORECASE).finditer(line_l)
                                X_count = 0
                                for X_icon in X_icons:
                                    X_count += 1
                                X_counts += [X_count]
                        # Collect names of all antibodies measured in cell lines                                                    
                        line_l = re.search('label">.*?<', ll)
                        if line_l:
                            line_l = line_l.group()
                            X_Ab_label = re.search('\>\w.*?\<', line_l).group()
                            X_Ab = re.sub('[\<\>]', '', X_Ab_label)
                            if X_Ab not in X_Abs:
                                X_Abs += [X_Ab]        
                    f_l.close()
                    # If antibody is not found in both cell line replicates, remove from list
                    i = 0
                    while i < len(X_counts):
                        if X_counts[i] < 2:
                            del(X_counts[i])
                            del(X_Abs[i])
                        else: i += 1
                    i = 0
                    # If antibody is not sufficiently measured in the cell line, remove from original list
                    while i < len(Abs):
                        if Abs[i] not in X_Abs:
                            del(Ab_scores[i])
                            del(Abs[i])
                        else: i += 1
                # If antibody has less than 2 supportive validity scores, remove from list
                num_valid = 0      
                for i in range(0, len(Ab_scores)):
                    if Ab_scores[i] >= 2:
                        num_valid += 1
                # If any valid and measured antibodies remain for the protein, save information
                if num_valid > 0:
                    file.write('{col1}\t{col2}\t{col3}\t{col4}\n'.format(col1=ID, col2=name, col3=synonyms, col4=desc))
                    add += 1
                else: discard += 1
                remain = genes - add - discard
                print '{a} proteins added, {d} discarded, {r} remaining'.format(a=add, d=discard, r=remain)            
        f.close()        
file.close()
print time.clock() - t_0
