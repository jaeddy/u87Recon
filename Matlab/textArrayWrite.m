function textArrayWrite(fileName,strArray,delim)

fid = fopen(fileName,'w');
for i = 1:size(strArray,1)
    for j = 1:size(strArray,2)-1
        fprintf(fid,['%s',delim],strArray{i,j});
    end
    fprintf(fid,'%s\n',strArray{i,end});
end
    
fclose(fid);