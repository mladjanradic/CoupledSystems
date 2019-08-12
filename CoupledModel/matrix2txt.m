function [] = matrix2txt(header,M,FN)


fid = fopen(FN, 'w') ;


for i=1:length(header)
    fprintf(fid,[header{i} '\t']);
end

fprintf(fid,'\n');

for i=1:size(M,1)
    for j=1:size(M,2)
        fprintf(fid,'%d\t',M(i,j) );
    end
    fprintf(fid,'\n');
end
fclose(fid); 
