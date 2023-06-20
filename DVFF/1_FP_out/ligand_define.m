%读取配体文件到矩阵中    ligand format : mol2
%v0.1.0.20220901
%v0.1.1.20220927   函数内打开关闭文件
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v0.1.2.20230616   final
function [ligand] = ligand_define(ligand_path)
fid = fopen(ligand_path,'rt');
flag_ligand = 0;
ligand = [];

while flag_ligand == 0
    tline = fgetl(fid);
    line_split = strsplit(tline);
    if strcmp(line_split{1},'@<TRIPOS>ATOM')
        flag_ligand = 1;
    end
    while flag_ligand == 1
        tline = fgetl(fid);
        if isempty(tline)
            break
        end
        tline = strtrim(tline);
        line_split = strsplit(tline);
        if strcmp(line_split{1},'@<TRIPOS>BOND')
            break
        end
        ligand = [ligand;line_split];
    end
end
fclose(fid);
