%输出蛋白质的菲波那契特征矩阵
%v0.1.1.20221011  修复了一个bug。此前版本无法识别原子类型第一个字符为数字1，2，3的情况。
%v0.1.2.20221011  在上一个版本的基础做了优化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v0.1.3.20230616  final
function [fp_protein]=fp_protein(protein_path)
[prepared_protein,atom_type] = protein_define(protein_path);

point = prepared_protein(:,1:3);
ZA = zeros(size(atom_type,1),1);
fibonacci_out = zeros(size(atom_type,1),1000); 
for i = 1:size(atom_type,1)
    if atom_type{i}(1)=='1'||atom_type{i}(1)=='2'||atom_type{i}(1)=='3'
        atom_type{i}(1)=atom_type{i}(2);
    end                              
    if atom_type{i}(1)=='H'
       ZA(i,1) = 1;    
    elseif atom_type{i}(1)=='C'                                     
       ZA(i,1) = 6;                   
    elseif atom_type{i}(1)=='N'
        ZA(i,1) = 7;
    elseif atom_type{i}(1)=='O'
        ZA(i,1) = 8;                               
    elseif atom_type{i}(1)=='S'
       ZA(i,1) = 16;                          
    else  
        sprintf('You may need to add some other atomic types')
    end
end

protein=[prepared_protein(:,1:3),prepared_protein(:,10),ZA];
Rcutoff = 10;
count = 1;
for i = 1:size(atom_type,1)
    tic
    Sphere_points = fibonacci(point(i,1),point(i,2),point(i,3));
    Sphere_points(:,4) = 0;
    protein_pick = protein;
    for j = size(protein,1):-1:1
        if (sqrt((point(i,1)-protein(j,1))^2+(point(i,2)-protein(j,2))^2+(point(i,3)-protein(j,3))^2)>Rcutoff)
            protein_pick(j,:) = [];
        end
    end
    for j = 1:size(protein_pick,1)
        clear pointj QAj
        pointj(1,1) = protein_pick(j,1);
        pointj(1,2) = protein_pick(j,2);
        pointj(1,3) = protein_pick(j,3);
        QAj = protein_pick(j,4);
        GAPAj = protein_pick(j,5) - QAj;  
        alphaj = 0.5*pi*GAPAj^(4/3);
        distance = dist(Sphere_points,pointj);
        Phi = GAPAj*exp(-alphaj*distance.^2);
        Sphere_points(:,4) = Sphere_points(:,4) + Phi;
    end
    toc
    fibonacci_out_transpose = sort(Sphere_points(:,4),'descend');
    fibonacci_out(count,1:1000) = fibonacci_out_transpose';
    count = count+1;
    
end
fp_protein = [fibonacci_out,point,ZA,(1:size(atom_type,1))'];
