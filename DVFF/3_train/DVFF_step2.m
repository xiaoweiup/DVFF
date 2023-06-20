%GARF_step3                                                                                                                            %GARF_step2
%v0.1.0.20221024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v0.1.1.20230619
distribution_P_path='/home/user/Documents/wjm/GARF/result/step1/PDBbind/general_set/1-1/PLA_P_distribution.mat';
distribution_L_path='/home/user/Documents/wjm/GARF/result/step1/PDBbind/general_set/1-1/PLA_L_distribution.mat';
ligand_atom_type_total = 46;        
protein_atom_type_total = 40; 
atom_type_total = ligand_atom_type_total+protein_atom_type_total;
Size_PLA_rows = 40000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialization
PLA_P_distA = zeros(200,ligand_atom_type_total*protein_atom_type_total*protein_atom_type_total);
PLA_L_distA = zeros(200,ligand_atom_type_total*protein_atom_type_total*ligand_atom_type_total);
PLA_P_distB = zeros(200,ligand_atom_type_total*protein_atom_type_total*protein_atom_type_total);
PLA_L_distB = zeros(200,ligand_atom_type_total*protein_atom_type_total*ligand_atom_type_total);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load (distribution_P_path)
load (distribution_L_path)

for i=1:200
    PLA_P_distA(i,:)=sum(PLA_P_distribution((i-1)*200+1:i*200,:),1);
    PLA_L_distA(i,:)=sum(PLA_L_distribution((i-1)*200+1:i*200,:),1);
end
for i=1:200
    for j=1:200
        PLA_P_distB(i,:)= PLA_P_distB(i,:)+PLA_P_distribution(200*(j-1)+i,:);
        PLA_L_distB(i,:)= PLA_L_distB(i,:)+PLA_L_distribution(200*(j-1)+i,:);
    end
end

sum_PLA_distB_columns = sum(PLA_P_distB,2) + sum(PLA_L_distB,2);

count_C = zeros(200,atom_type_total);
for i=1:atom_type_total
    for j=1:ligand_atom_type_total*protein_atom_type_total
        if(i<=ligand_atom_type_total)
            count_C(:,i)=count_C(:,i)+PLA_L_distB(:,ligand_atom_type_total*(j-1)+i);
        else
            count_C(:,i)=count_C(:,i)+PLA_P_distB(:,protein_atom_type_total*(j-1)+i-ligand_atom_type_total);
        end
    end
end
P_C = zeros(200,ligand_atom_type_total*protein_atom_type_total*atom_type_total);
for i=1:ligand_atom_type_total*protein_atom_type_total
    P_C(:,(i-1)*atom_type_total+1:i*atom_type_total) = count_C./sum_PLA_distB_columns;   %P(C) 
end
          

count_CS = zeros(200,ligand_atom_type_total*protein_atom_type_total*atom_type_total);
count_SC = zeros(200,ligand_atom_type_total*protein_atom_type_total*atom_type_total);
for i=1:ligand_atom_type_total*protein_atom_type_total
    count_CS(:,atom_type_total*(i-1)+1:atom_type_total*i)=[PLA_L_distB(:,ligand_atom_type_total*(i-1)+1:ligand_atom_type_total*i),PLA_P_distB(:,protein_atom_type_total*(i-1)+1:protein_atom_type_total*i)];
    count_SC(:,atom_type_total*(i-1)+1:atom_type_total*i)=[PLA_L_distA(:,ligand_atom_type_total*(i-1)+1:ligand_atom_type_total*i),PLA_P_distA(:,protein_atom_type_total*(i-1)+1:protein_atom_type_total*i)];
end
for i=1:ligand_atom_type_total*protein_atom_type_total*atom_type_total
    count_CS(:,i)=count_CS(:,i)/sum(count_CS(:,i),1);
    count_SC(:,i)=count_SC(:,i)/sum(count_SC(:,i),1);
end
P_CS = count_CS;                           %P(C|S)
P_SC = count_SC;                           %P(S|C)

P_S = P_SC.*P_C./P_CS;                     %P(S) = P(S|C)*P(C)/P(C|S)  
P_S(isnan(P_S)) = 0;
P_S(isinf(P_S)) = 0;

P_SS = zeros(200,ligand_atom_type_total*protein_atom_type_total);
for i=1:ligand_atom_type_total*protein_atom_type_total
    P_SS(:,i)=sum(P_S(:,(i-1)*atom_type_total+1:i*atom_type_total),2)/atom_type_total;  %result
end
