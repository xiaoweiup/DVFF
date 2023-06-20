%v0.1.1.20220909   
%v0.1.2.20230616      final   
function [fp_ligand]=fp_ligand(ligand_path)


ligand_atom = ligand_define(ligand_path);

point = zeros(size(ligand_atom,1),3);                                      %atomic coordinates
ZA = zeros(size(ligand_atom,1),1);                                         %atomic number
charge = zeros(size(ligand_atom,1),1);   
fibonacci_out = zeros(size(ligand_atom,1),1000);                           %fibonacci points

for i = 1:size(ligand_atom,1)
        
    point(i,1) = str2double(ligand_atom{i,3});
    point(i,2) = str2double(ligand_atom{i,4});
    point(i,3) = str2double(ligand_atom{i,5});
    charge(i,1) = str2double(ligand_atom{i,9});
    atom_type = ligand_atom{i,2};        
    if atom_type(1)=='H'
       ZA(i,1) = 1;
    elseif atom_type(1)=='B'
        if size(atom_type,1)>1
            if atom_type(2)=='r'
                ZA(i,1) = 35;
            else
               ZA(i,1) = 5;
            end
        else
            ZA(i,1) = 5;
        end             
    elseif atom_type(1)=='C'                                     
        if size(atom_type,1)>1
            if atom_type(2)=='l'
                ZA(i,1) = 17;
            else
                ZA(i,1) = 6;
            end
        else
            ZA(i,1) = 6;
        end                 
    elseif atom_type(1)=='N'
        ZA(i,1) = 7;
    elseif atom_type(1)=='O'
        ZA(i,1) = 8;      
    elseif atom_type(1)=='F'
        ZA(i,1) = 9;  
    elseif atom_type(1)=='P'
        ZA(i,1) = 15;        
    elseif atom_type(1)=='S'
        ZA(i,1) = 16;                
    elseif atom_type(1)=='I'
        ZA(i,1) = 53;    
    else  
        sprintf('You may need to add some other atomic types')
    end       
end 

count = 1;
for i = 1:size(ligand_atom,1)
        
    Sphere_points = fibonacci(point(i,1),point(i,2),point(i,3));
    Sphere_points(:,4) = 0;
    for j = 1:size(ligand_atom,1)
        clear pointj QAj
        pointj(1,1) = point(j,1);
        pointj(1,2) = point(j,2);
        pointj(1,3) = point(j,3);
        QAj = charge(j,1);
        GAPAj = ZA(j,1) - QAj;  
        alphaj = 0.5*pi*GAPAj^(4/3);
        distance = dist(Sphere_points,pointj);
        Phi = GAPAj*exp(-alphaj*distance.^2);
        Sphere_points(:,4) = Sphere_points(:,4) + Phi;
    end 
    fibonacci_out_transpose = sort(Sphere_points(:,4),'descend');
    fibonacci_out(count,1:1000) = fibonacci_out_transpose';
    count = count+1;
end

fp_ligand = [fibonacci_out,point,ZA,charge,(1:size(ligand_atom,1))'];

