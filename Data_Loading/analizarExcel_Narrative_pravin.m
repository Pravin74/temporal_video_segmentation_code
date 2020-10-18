function [GT]=analizarExcel_Narrative_pravin(excelfile, fps)
[~,textA] = xlsread(excelfile);
GT=[];
for i= 2:length(textA)
    GT =[GT ;fps/2 + fps*(3600*str2num(textA{i}(1))+60*str2num(textA{i}(3:4))+str2num(textA{i}(6:7)))];
end
GT = floor(GT);