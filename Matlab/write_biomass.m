
%%
[a,b] = xlsread('U87 biomass.xlsx','Sheet2','B1:C53');

b = lower(b);
b = regexprep(b,'\-l','\_L');
b = strcat(b,{'[c]'});
b = [cellstr(num2str(a)),b];

%%
fid = fopen('U87_biomass.txt','w');
for i = 1:size(b,1)
    fprintf(fid,'%s\t',num2str(b{i,1}));
    fprintf(fid,'%s\t',b{i,2});
end
    
fclose(fid);