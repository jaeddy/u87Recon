function strArray = textArrayRead(fileName,numCol,delim)

% Example inputs:
%  fileName = 'GOPathway.txt';
%  numCol = 306;
%  delim = '\t';

outStr = '';
for i = 1:numCol
    outStr = [outStr,'col',num2str(i),','];
end
outStr(end) = [];

frmtStr = ['%s',repmat(' %s',1,numCol-1)];

eval(['[',outStr,'] = textread(fileName,frmtStr,''emptyvalue'',1,',...
    '''delimiter'',delim);']);

strArray = col1;
for i = 2:numCol
    eval(['strArray = [strArray,col',num2str(i),'];']);
end