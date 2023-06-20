%GARF step1+step2
%v0.2.0.20221012 将新原子分类嵌入  构建了两个编码器储存输出数据

%v0.2.1.20230616   final
% file_dir='/home/user/Documents/wjm_test/GARF/test/test';
% file_dir='/home/user/Documents/wjm/PDBbind/refined-set3/am1-bcc/r_mono_m_atom';
file_dir='/home/user/Desktop/SCARF/input/general_set/4/4-1';
load('Cluster_H_ligand.mat')
load('Cluster_C_ligand.mat')
load('Cluster_O_ligand.mat')
load('Cluster_N_ligand.mat')
load('Cluster_H_protein.mat')
load('Cluster_C_protein.mat')
load('Cluster_N_protein.mat')
load('Cluster_O_protein.mat')
load('Cluster_S_protein.mat')
RcutoffPL = 10;

ms = 0.1;
max_x = 20;
max_y = 20;
Size_ligand_fibonacci = 46;        
Size_protein_fibonacci = 40; 
Size_PLA_rows = (max_x/ms)*(max_y/ms);
Size_PLA_P_columns = Size_protein_fibonacci*Size_ligand_fibonacci*Size_protein_fibonacci;
PLA_P_distribution = zeros(Size_PLA_rows,Size_PLA_P_columns);
Size_PLA_L_columns = Size_protein_fibonacci*Size_ligand_fibonacci*Size_ligand_fibonacci;
PLA_L_distribution = zeros(Size_PLA_rows,Size_PLA_L_columns);



%%
list=dir(file_dir);
for i=3:size(list,1)  
    tic
    if length(list(i).name)==4
        ligand_name=list(i).name; 
        protein_path=strcat(file_dir,'/',ligand_name(1,1:4),'/',ligand_name(1,1:4),'_protein.pdb');
        ligand_path=strcat(file_dir,'/',ligand_name(1,1:4),'/',ligand_name(1,1:4),'_ligand.mol2');
        
        fp_ligand = fp_ligand(ligand_path);      
        for j=1:size(fp_ligand,1)
            if fp_ligand(j,1004)==1
                [~,idx] = pdist2(Cluster_H_ligand,fp_ligand(j,1:1000),'euclidean','Smallest',1);
                fp_ligand(j,1006) = idx;
            elseif fp_ligand(j,1004)==5               
                fp_ligand(j,1006) = idx+size(Cluster_H_ligand,1);
            elseif fp_ligand(j,1004)==6
                [~,idx] = pdist2(Cluster_C_ligand,fp_ligand(j,1:1000),'euclidean','Smallest',1);
                fp_ligand(j,1006) = 1+idx+size(Cluster_H_ligand,1); 
            elseif fp_ligand(j,1004)==7
                [~,idx] = pdist2(Cluster_N_ligand,fp_ligand(j,1:1000),'euclidean','Smallest',1);
                fp_ligand(j,1006) = 1+idx+size(Cluster_H_ligand,1)+size(Cluster_C_ligand,1);
            elseif fp_ligand(j,1004)==8
                [~,idx] = pdist2(Cluster_O_ligand,fp_ligand(j,1:1000),'euclidean','Smallest',1);
                fp_ligand(j,1006) = 1+idx+size(Cluster_H_ligand,1)+size(Cluster_C_ligand,1)+size(Cluster_N_ligand,1);
            elseif fp_ligand(j,1004)==9
                fp_ligand(j,1006) = 2+size(Cluster_H_ligand,1)+size(Cluster_C_ligand,1)+size(Cluster_N_ligand,1)+size(Cluster_O_ligand,1);    
            elseif fp_ligand(j,1004)==15
                fp_ligand(j,1006) = 3+size(Cluster_H_ligand,1)+size(Cluster_C_ligand,1)+size(Cluster_N_ligand,1)+size(Cluster_O_ligand,1);
            elseif fp_ligand(j,1004)==16
                fp_ligand(j,1006) = 4+size(Cluster_H_ligand,1)+size(Cluster_C_ligand,1)+size(Cluster_N_ligand,1)+size(Cluster_O_ligand,1);
            elseif fp_ligand(j,1004)==17
                fp_ligand(j,1006) = 5+size(Cluster_H_ligand,1)+size(Cluster_C_ligand,1)+size(Cluster_N_ligand,1)+size(Cluster_O_ligand,1);
            elseif fp_ligand(j,1004)==34
                fp_ligand(j,1006) = 6+size(Cluster_H_ligand,1)+size(Cluster_C_ligand,1)+size(Cluster_N_ligand,1)+size(Cluster_O_ligand,1);
            elseif fp_ligand(j,1004)==35
                fp_ligand(j,1006) = 7+size(Cluster_H_ligand,1)+size(Cluster_C_ligand,1)+size(Cluster_N_ligand,1)+size(Cluster_O_ligand,1);
            elseif fp_ligand(j,1004)==53
                fp_ligand(j,1006) = 8+size(Cluster_H_ligand,1)+size(Cluster_C_ligand,1)+size(Cluster_N_ligand,1)+size(Cluster_O_ligand,1);
            end
        end
%         Size_ligand_fibonacci = 8+size(Cluster_H_ligand,1)+size(Cluster_C_ligand,1)+size(Cluster_N_ligand,1)+size(Cluster_O_ligand,1);
        ligand_refine = [fp_ligand(:,1001:1004),(1:size(fp_ligand,1))',fp_ligand(:,1006)];
        
        fp_protein = fp_protein(protein_path);
        for j=1:size(fp_protein,1)
            if fp_protein(j,1004)==1
                [~,idx] = pdist2(Cluster_H_protein,fp_protein(j,1:1000),'euclidean','Smallest',1);
                fp_protein(j,1006) = idx+Size_ligand_fibonacci;
            elseif fp_protein(j,1004)==6
                [~,idx] = pdist2(Cluster_C_protein,fp_protein(j,1:1000),'euclidean','Smallest',1);
                fp_protein(j,1006) = idx+Size_ligand_fibonacci+size(Cluster_H_protein,1);
            elseif fp_protein(j,1004)==7
                [~,idx] = pdist2(Cluster_N_protein,fp_protein(j,1:1000),'euclidean','Smallest',1);
                fp_protein(j,1006) = idx+Size_ligand_fibonacci+size(Cluster_H_protein,1)+size(Cluster_C_protein,1);    
            elseif fp_protein(j,1004)==8
                [~,idx] = pdist2(Cluster_O_protein,fp_protein(j,1:1000),'euclidean','Smallest',1);
                fp_protein(j,1006) = idx+Size_ligand_fibonacci+size(Cluster_H_protein,1)+size(Cluster_C_protein,1)+size(Cluster_N_protein,1);
            elseif fp_protein(j,1004)==16
                [~,idx] = pdist2(Cluster_S_protein,fp_protein(j,1:1000),'euclidean','Smallest',1);
                fp_protein(j,1006) = idx+Size_ligand_fibonacci+size(Cluster_H_protein,1)+size(Cluster_C_protein,1)+size(Cluster_N_protein,1)+size(Cluster_O_protein,1);
            end
        end
        protein_refine = fp_protein(:,1001:1006);
        
        [protein_starterA,protein_starterB] = pocket2find_PL_AFA(protein_refine,ligand_refine,RcutoffPL);
        size_protein_starter=size(protein_starterA,1);
        for j=1:size_protein_starter
            for k=1:size(protein_refine,1)
                if protein_starterA(j,5)~=protein_refine(k,5) && norm(protein_starterA(j,1:3)-protein_refine(k,1:3))<=2
                    protein_adjacent = protein_refine(k,:);
                    distA = sqrt((protein_adjacent(1,1)-protein_starterB(j,1))^2+(protein_adjacent(1,2)-protein_starterB(j,2))^2+(protein_adjacent(1,3)-protein_starterB(j,3))^2);
                    distB = sqrt((protein_starterA(j,1)-protein_starterB(j,1))^2+(protein_starterA(j,2)-protein_starterB(j,2))^2+(protein_starterA(j,3)-protein_starterB(j,3))^2);
                    PLA_rows = floor(distA/ms)*(max_y/ms) + floor(distB/ms) +1 ;
                    PLA_P_columns = (protein_starterA(j,6)-Size_ligand_fibonacci-1)*Size_ligand_fibonacci*Size_protein_fibonacci+(protein_starterB(j,6)-1)*Size_protein_fibonacci+protein_adjacent(1,6)-Size_ligand_fibonacci;
                    PLA_P_distribution(PLA_rows,PLA_P_columns)=PLA_P_distribution(PLA_rows,PLA_P_columns)+1;
                end
                    
            end
            for k=1:size(ligand_refine,1)    
                if protein_starterB(j,5)~=ligand_refine(k,5) && norm(protein_starterB(j,1:3)-ligand_refine(k,1:3))<=2
                    ligand_adjacent = ligand_refine(k,:);
                    distA = sqrt((ligand_adjacent(1,1)-protein_starterA(j,1))^2+(ligand_adjacent(1,2)-protein_starterA(j,2))^2+(ligand_adjacent(1,3)-protein_starterA(j,3))^2);
                    distB = sqrt((protein_starterA(j,1)-protein_starterB(j,1))^2+(protein_starterA(j,2)-protein_starterB(j,2))^2+(protein_starterA(j,3)-protein_starterB(j,3))^2);
                    PLA_rows = floor(distA/ms)*(max_y/ms) + floor(distB/ms) +1 ;
                    PLA_L_columns = (protein_starterA(j,6)-Size_ligand_fibonacci-1)*Size_ligand_fibonacci*Size_ligand_fibonacci+(protein_starterB(j,6)-1)*Size_ligand_fibonacci+ligand_adjacent(1,6);
                    PLA_L_distribution(PLA_rows,PLA_L_columns)=PLA_L_distribution(PLA_rows,PLA_L_columns)+1;
                end
            end
        end
    end
    toc
end


%%
% ligand H1 1         ligand C13 26         protein H1 47           protein N1 73                     
% ligand H2 2         ligand N1 27          protein H2 48           protein N2 74 
% ligand H3 3         ligand N2 28          protein H3 49           protein N3 75 
% ligand H4 4         ligand N3 29          protein H4 50           protein N4 76 
% ligand H5 5         ligand N4 30          protein H5 51           protein N5 77 
% ligand H6 6         ligand N5 31          protein H6 52           protein
% N6 78 -
% ligand H7 7         ligand N6 32          protein H7 53           protein O1 79
% ligand H8 8         ligand N7 33          protein H8 54           protein O2 80
% ligand H9 9         ligand O1 34          protein H9 55           protein O3 81
% ligand H10 10       ligand O2 35          protein H10 56          protein O4 82
% ligand H11 11       ligand O3 36          protein H11 57          protein O5 83
% ligand H12 12       ligand O4 37          protein C1 58           protein S1 84
% ligand B 13         ligand O5 38          protein C2 59           protein S2 85
% ligand C1 14        ligand O6 39          protein C3 60           protein S3 86
% ligand C2 15        ligand F 40           protein C4 61
% ligand C3 16        ligand P 41           protein C5 62
% ligand C4 17        ligand S 42           protein C6 63
% ligand C5 18        ligand Cl 43          protein C7 64
% ligand C6 19        ligand Se 44          protein C8 65
% ligand C7 20        ligand Br 45          protein C9 66
% ligand C8 21        ligand I 46           protein C10 67
% ligand C9 22                              protein C11 68
% ligand C10 23                             protein C12 69
% ligand C11 24                             protein C13 70
% ligand C12 25                             protein C14 71
%                                           protein C15 72
