%protein_dir = 'C:\Users\weijianmin\Desktop\protein_define\184l_protein.pdb';
%读取蛋白质文件到矩阵中    protein format : pdb
%v0.1.1.20230616       final
function [prepared_protein_out,AtomName] = protein_define(protein_path)
%% digitize protein atom types
C=1  ;
Cs=2 ;
CA=3 ;
CB=4 ;
CC=5 ;
CN=6 ;
CR=7 ;
CT=8 ;
CV=9 ;
CW=10;
N=11;
N2=12;
N3=13;
NA=14;
NB=15;
O=16 ;
O2=17;
OH=18;
S=19;
SH=20;
H=21 ;
HO=22;
HC=23;
H1=24;
H2=25;
H3=26;
H4=27;
H5=28;
HA=29;
HS=30;



CAL   = 31	;
AL   = 32	;
MN   = 33  	;
FE   = 34    ;
CO   = 35    ;
NI   = 36    ;
CU   = 37    ;
ZN   = 38    ;
K    = 39 ;
MG   = 40  ;

N_A=0;
D=1;
D2=2;
A=3;
DA=4;

%% digitize residue names
ALA=101;
ARG=102;
ASN=103;
ASP=104;
CYS=105;
GLN=106;
GLU=107;
GLY=108;
HIS=109;
ILE=110;
LEU=111;
LYS=112;
MET=113;
PHE=114;
PRO=115;
SER=116;
THR=117;
TRP=118;
TYR=119;
VAL=120;
%% protein atoms initialization
%read protein atomic information from the pdb file
clear PDBData1 protein_atom1 protein_atom2 protein_atom3 protein_atomname protein_resseq protein_occupancy protein_atom_f protein_atom resName1


fidin=fopen(protein_path);
tline=fgetl(fidin);
k=1;
num=1;
clear AtomName resName chainID protein_atom

while ~feof(fidin)
    
    if length(tline)>6 && ((strcmp(tline(1:4), 'ATOM')==1))
        m=1;
        ct=0;
        while m<length(tline)
            if (strcmp(tline(m), ' ')==1) && ct==0
                ct=1;
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==1
                ct=2;
                m=m+1;
            elseif (strcmp(tline(m), ' ')==1) && ct==2
                ct=3;
                m=m+1;
                
            elseif (strcmp(tline(m), ' ')~=1) && ct==3
                ct=4;
                for mi=1:3
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=5;
                        AtomName{k,1}=tline(m:m+mi-1);
                        m=m+mi;
                        break
                    end
                end
                if ct ==4
                    AtomName{k,1}=tline(m:m+mi);
                    ct=5;
                    m=m+mi+1;
                end
                %cutA(1,2)=m+mi-1;
                
                %m=m+mi;
            elseif (strcmp(tline(m), ' ')==1) && ct==5
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==5
                ct=6;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=7;
                        break
                    end
                end
                RN_test = tline(m:m+mi-1);
                if size(RN_test,2)>3
                    RN_test=RN_test(end-2:end);
                    
                end
                
                resName{k,:}=RN_test;
                
                m=m+mi;
                
            elseif (strcmp(tline(m), ' ')==1) && ct==7
                m=m+3;
                ct=9;
                %%%%%%%%%%%%%%%%%chain ID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %elseif (strcmp(tline(m), ' ')~=1) && ct==7
                %    ct=8;
                %    for mi=1:10
                %        if (strcmp(tline(m+mi), ' ')==1)
                %            ct=9;
                %            break
                %        end
                %    end
                %    chainID{k,1}=tline(m:m+mi-1);
                %    m=m+mi;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif (strcmp(tline(m), ' ')==1) && ct==9
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==9
                ct=10;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=11;
                        break
                    end
                end
                str1=tline(m:m+mi-1);
                resName1(k,1)=str2double(str1(regexp(str1,'\d')));
                
                %resSeq(k,1)=str2num(str1(regexp(str1,'\d')));
                m=m+mi;
                
                
            elseif (strcmp(tline(m), ' ')==1) && ct==11
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==11
                ct=12;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        ct=13;
                        break
                    end
                end
                protein_atom(k,1)=str2double(tline(m:m+mi-1));
                m=m+mi;
                
            elseif (strcmp(tline(m), ' ')==1) && ct==13
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==13
                ct=14;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        ct=15;
                        break
                    end
                end
                protein_atom(k,2)=str2double(tline(m:m+mi-1));
                m=m+mi;
                
            elseif (strcmp(tline(m), ' ')==1) && ct==15
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==15
                ct=16;
                for mi=1:10
                    if length(tline) == (m+mi);
                        protein_atom(k,3)=str2double(tline(m:m+mi));
                        ct=17;
                        break
                        
                    elseif (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        protein_atom(k,3)=str2double(tline(m:m+mi-1));
                        ct=17;
                        break
                    end
                end
                
                m=m+mi;
            elseif ct==17
                break
            else
                m=m+1;
            end
        end
        
        k=k+1;
        tline=fgetl(fidin);
        
    elseif (length(tline)==3) && (strcmp(tline(1:3), 'END')==1)
        break
    else
        tline=fgetl(fidin);
    end
end

if (length(tline)>6)
    if (strcmp(tline(1:3), 'END')~=1) && ((strcmp(tline(1:4), 'ATOM')==1))
    %if (strcmp(tline(1:3), 'END')~=1) && ((strcmp(tline(1:4), 'ATOM')==1)||(strcmp(tline(1:6), 'HETATM')==1))
        m=1;
        ct=0;
        while m<length(tline)
            if (strcmp(tline(m), ' ')==1) && ct==0
                ct=1;
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==1
                ct=2;
                m=m+1;
            elseif (strcmp(tline(m), ' ')==1) && ct==2
                ct=3;
                m=m+1;
                
            elseif (strcmp(tline(m), ' ')~=1) && ct==3
                ct=4;
                for mi=1:3
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=5;
                        AtomName{k,1}=tline(m:m+mi-1);
                        break
                    end
                end
                if ct ==4
                    AtomName{k,1}=tline(m:m+mi);
                    ct=5;
                end
                %cutA(1,2)=m+mi-1;
                
                m=m+mi;
            elseif (strcmp(tline(m), ' ')==1) && ct==5
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==5
                ct=6;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=7;
                        break
                    end
                end
                RN_test = tline(m:m+mi-1);
                if size(RN_test,2)>3
                    RN_test=RN_test(end-2:end);
                    
                end
                
                resName{k,:}=RN_test;
                
                m=m+mi;
                
            elseif (strcmp(tline(m), ' ')==1) && ct==7
                m=m+3;
                ct=9;
                %%%%%%%%%%%%%%%%%chain ID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %elseif (strcmp(tline(m), ' ')~=1) && ct==7
                %    ct=8;
                %    for mi=1:10
                %        if (strcmp(tline(m+mi), ' ')==1)
                %            ct=9;
                %            break
                %        end
                %    end
                %    chainID{k,1}=tline(m:m+mi-1);
                %    m=m+mi;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif (strcmp(tline(m), ' ')==1) && ct==9
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==9
                ct=10;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=11;
                        break
                    end
                end
                str1=tline(m:m+mi-1);
                resName1(k,1)=str2double(str1);
                
                %resSeq(k,1)=str2double(str1(regexp(str1,'\d')));
                m=m+mi;
                
                
            elseif (strcmp(tline(m), ' ')==1) && ct==11
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==11
                ct=12;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        ct=13;
                        break
                    end
                end
                protein_atom(k,1)=str2double(tline(m:m+mi-1));
                m=m+mi;
                
            elseif (strcmp(tline(m), ' ')==1) && ct==13
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==13
                ct=14;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        ct=15;
                        break
                    end
                end
                protein_atom(k,2)=str2double(tline(m:m+mi-1));
                m=m+mi;
                
            elseif (strcmp(tline(m), ' ')==1) && ct==15
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==15
                ct=16;
                for mi=1:10
                    if length(tline) == (m+mi);
                        protein_atom(k,3)=str2double(tline(m:m+mi));
                        ct=17;
                        break
                        
                    elseif (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        protein_atom(k,3)=str2double(tline(m:m+mi-1));
                        ct=17;
                        break
                    end
                end
                %protein_atom(k,3)=str2num(tline(m:m+mi-1));
                m=m+mi;
            elseif ct==17
                break
            else
                m=m+1;
            end
        end
        
        %k=k+1;
        %tline=fgetl(fidin);
    end
end

fclose(fidin);
%%
resSeq = zeros(size(resName,1),1);
for i = 1:size(resName,1)
    if strcmp(resName{i,1}, 'ALA')==1
        resSeq(i,1) = 101;
    elseif strcmp(resName{i,1}, 'ARG')==1
        resSeq(i,1) = 102;
    elseif strcmp(resName{i,1}, 'ASN')==1 || strcmp(resName{i,1}, 'ASX')==1
        resSeq(i,1) = 103;
    elseif strcmp(resName{i,1}, 'ASP')==1
        resSeq(i,1) = 104;
    elseif strcmp(resName{i,1}, 'CYS')==1
        resSeq(i,1) = 105;
    elseif strcmp(resName{i,1}, 'GLN')==1
        resSeq(i,1) = 106;
    elseif strcmp(resName{i,1}, 'GLU')==1 || strcmp(resName{i,1}, 'GLX')==1
        resSeq(i,1) = 107;
    elseif strcmp(resName{i,1}, 'GLY')==1
        resSeq(i,1) = 108;
    elseif strcmp(resName{i,1}, 'HIS')==1
        resSeq(i,1) = 109;
    elseif strcmp(resName{i,1}, 'HID')==1
        resSeq(i,1) = 109.1;
    elseif strcmp(resName{i,1}, 'HIE')==1
        resSeq(i,1) = 109.2;
    elseif strcmp(resName{i,1}, 'HIP')==1
        resSeq(i,1) = 109.3;
    elseif strcmp(resName{i,1}, 'ILE')==1
        resSeq(i,1) = 110;
    elseif strcmp(resName{i,1}, 'LEU')==1
        resSeq(i,1) = 111;
    elseif strcmp(resName{i,1}, 'LYS')==1
        resSeq(i,1) = 112;
    elseif strcmp(resName{i,1}, 'MET')==1
        resSeq(i,1) = 113;
    elseif strcmp(resName{i,1}, 'PHE')==1
        resSeq(i,1) = 114;
    elseif strcmp(resName{i,1}, 'PRO')==1
        resSeq(i,1) = 115;
    elseif strcmp(resName{i,1}, 'SER')==1
        resSeq(i,1) = 116;
    elseif strcmp(resName{i,1}, 'THR')==1
        resSeq(i,1) = 117;
    elseif strcmp(resName{i,1}, 'TRP')==1
        resSeq(i,1) = 118;
    elseif strcmp(resName{i,1}, 'TYR')==1
        resSeq(i,1) = 119;
    elseif strcmp(resName{i,1}, 'VAL')==1
        resSeq(i,1) = 120;
    else
        resSeq(i,1) = 121;
    end
end
%% assign parameters to each protein atom
%clear protein_name protein_atom_radius protein_atom_vdwep protein_acceptor_hb_angle protein_donor_hb_angle protein_atom_charge protein_backbone
protein_name=zeros(size(protein_atom,1),1);
protein_atom_radius=zeros(size(protein_atom,1),1);
protein_atom_vdwep=zeros(size(protein_atom,1),1);
protein_acceptor_hb_angle=zeros(size(protein_atom,1),1);
protein_donor_hb_angle=zeros(size(protein_atom,1));
protein_atom_charge=zeros(size(protein_atom,1),1);
protein_backbone=zeros(size(protein_atom,1),1);

protein_resseq=resSeq;
size_protein_atom=size(protein_atom);
for i=1:size_protein_atom(1,1)
    
    if (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'HOH', 3)==1)% && (length(AtomName{i,1})==1)
        protein_name(i,1)=OH;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.7683;
        protein_atom_vdwep(i)=0.1520;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.6546;
        protein_backbone(i)=0;
        
    elseif (strcmp(AtomName{i,1}, 'HNCA')==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1984;
        protein_backbone(i)=1;
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.7341;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5894;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.3479;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2747;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2637;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4590;
        protein_atom_vdwep(i)=0.0150;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1560;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0007;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0327;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.039;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HG'),2)>0) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0285;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0486;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HD'),2)>0) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0687;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NE', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=N2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.5295;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HE'),2)>0) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.3456;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'CZ', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.8076;
        protein_backbone(i)=0;
        
    elseif (size(strfind(AtomName{i,1}, 'NH1'),3)>0) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=N2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.8627;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'NH2'),3)>0) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=N2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.8627;
        protein_backbone(i)=0;
        
    elseif (size(strfind(AtomName{i,1}, 'HH'),2)>0) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.4478;
        protein_backbone(i)=0;
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0143;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1048;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2041;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0797;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.7130;
        protein_backbone(i)=0;
        
    elseif (size(strfind(AtomName{i,1}, 'OD'),2)>0) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5931;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'ND'),2)>0) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.9191;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HD'),2)>0) && (strncmp(resName{i,1}, 'ASN', 3)==1 || strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.4196;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'OD1', 3)==1 || strncmp(AtomName{i,1}, 'XD1', 3)==1 ) && (strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5931;
        protein_backbone(i)=0;
        resName{i,1} = 'ASN';
        AtomName{i,1} = 'OD1';
    elseif (strncmp(AtomName{i,1}, 'ND2', 3)==1 || strncmp(AtomName{i,1}, 'XD2', 3)==1 ) && (strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.9191;
        protein_backbone(i)=0;
        resName{i,1} = 'ASN';
        AtomName{i,1} = 'ND2';
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5366;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5819;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.5163;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2936;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0381;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.088;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0303;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0122;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.7994;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'OD'),2)>0) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=O2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.8014;
        protein_backbone(i)=0;
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0213;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1124;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1231;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1112;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}(1), 'SG', 1)==1) && (strncmp(resName{i,1}, 'CYS', 3)==1 || strncmp(resName{i,1}, 'CYX', 3)==1)
        %protein_name(i,1)=SH;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=2.000;
        protein_atom_vdwep(i)=0.2500;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_backbone(i)=0;
        sig=0;
        for ip = 1:5
            if i+ip<=size(AtomName,1)
                
                if (size(strfind(AtomName{i+ip,1}, 'HG'),2)==0) && resName1(i+ip,1) == resName1(i,1)
                    protein_atom_charge(i)=-0.1081;
                    protein_name(i,1)=SH;
                    sig=1;
                    break
                end
            elseif i+ip>size(AtomName,1) && sig==0
                sig=0;
                break
            end
            
        end
        
        if sig ==0
            protein_atom_charge(i)=-0.3119;
            protein_name(i,1)=S;
        end
        
    elseif (size(strfind(AtomName{i,1}, 'HG'),2)>0) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=HS;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.1933;
        protein_backbone(i)=0;
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0031;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.085;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0036;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0171;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0645;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HG'),2)>0) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0352;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.6951;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'OE'),2)>0) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.6086;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'NE'),2)>0) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.9407;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HE'),2)>0) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.4251;
        protein_backbone(i)=0;
        
        
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5366;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5819;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.5163;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2936;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0397;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1105;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.056;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0173;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0136;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HG'),2)>0) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0425;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.8054;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'OE'),2)>0) && (strncmp(resName{i,1}, 'GLU', 3)==1 || strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=O2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.8188;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'OE'),2)>0 || size(strfind(AtomName{i,1}, 'XE'),2)>0) && (strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=O2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.8188;
        protein_backbone(i)=0;
        resName{i,1} = 'GLU';
        AtomName{i,1} = 'OE1';
        
        
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'HIS', 3)==1 || strncmp(resName{i,1}, 'HID', 3)==1 || strncmp(resName{i,1}, 'HIE', 3)==1 || strncmp(resName{i,1}, 'HIP', 3)==1)
        Hn1=0;
        Hn2=0;
        if size_protein_atom(1,1) >= i+17
            for ip=1:17
                if ((size(strfind(AtomName{i+ip,1}, 'HD2'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1))
                    Hn1=1;
                end
                if ((size(strfind(AtomName{i+ip,1}, 'HE2'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1))
                    Hn2=1;
                    break
                end
            end
        elseif size_protein_atom(1,1) < i+17
            for ip=1:size_protein_atom(1,1)-i
                if ((size(strfind(AtomName{i+ip,1}, 'HD2'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1))
                    Hn1=1;
                end
                if ((size(strfind(AtomName{i+ip,1}, 'HE2'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1))
                    Hn2=1;
                    break
                end
            end
        end
        ip=0;
        
        
        
        while ip <= size_protein_atom(1,1)-i && resName1(i+ip,1) == resName1(i,1)
            if Hn1 == 1 && Hn2 == 0
                
                
                if (strcmp(AtomName{i+ip,1}, 'N')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=N;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=pi;
                    protein_atom_charge(i+ip)=-0.4157;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'C')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=C;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.5973;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'O')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=O;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.6612;
                    protein_atom_vdwep(i+ip)=0.2100;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.5679;
                    protein_backbone(i+ip)=1;
                elseif (strncmp(AtomName{i+ip,1}, 'OXT', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=O;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.6612;
                    protein_atom_vdwep(i+ip)=0.2100;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.5679;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'H')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=H;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=0.6000;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.2719;
                    protein_backbone(i+ip)=1;
                    
                elseif (strncmp(AtomName{i+ip,1}, 'CA', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=CT;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.1094;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.0188;
                    protein_backbone(i+ip)=1;
                elseif (size(strfind(AtomName{i+ip,1}, 'HA'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=H1;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.3870;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.0881;
                    protein_backbone(i+ip)=1;
                elseif (strncmp(AtomName{i+ip,1}, 'CB', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=CT     ;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.1094;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.0462;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HB'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=HC;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.4870;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.0402;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'ND1', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=NA;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=pi;
                    protein_atom_charge(i+ip)=-0.3811;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HD1'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=H;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=0.6000;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.3649;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'CD2', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=CV;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1292;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HD2'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=H4;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.4090;
                    protein_atom_vdwep(i+ip)=0.0150;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1147;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'CE1', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=CR;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.2057;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HE1'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=H5;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.3590;
                    protein_atom_vdwep(i+ip)=0.0150;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1392;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'NE2', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=NB;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.5727;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'CG', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=CC;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.0266;
                    protein_backbone(i+ip)=0;
                end
                
            elseif Hn1 == 0 && Hn2 == 1
                
                
                if (strcmp(AtomName{i+ip,1}, 'N')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=N;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=pi;
                    protein_atom_charge(i+ip)=-0.4157;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'C')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=C;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.5973;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'O')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=O;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.6612;
                    protein_atom_vdwep(i+ip)=0.2100;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.5679;
                    protein_backbone(i+ip)=1;
                elseif (strncmp(AtomName{i+ip,1}, 'OXT', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=O;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.6612;
                    protein_atom_vdwep(i+ip)=0.2100;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.5679;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'H')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=H;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=0.6000;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.2719;
                    protein_backbone(i+ip)=1;
                    
                    
                elseif (strncmp(AtomName{i+ip,1}, 'CA', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=CT;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.1094;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.0581;
                    protein_backbone(i+ip)=1;
                elseif (size(strfind(AtomName{i+ip,1}, 'HA'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=H1;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.3870;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1360;
                    protein_backbone(i+ip)=1;
                elseif (strncmp(AtomName{i+ip,1}, 'CB', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=CT     ;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.1094;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.0074;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HB'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=HC;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.4870;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.0367;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'ND1', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=NB;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=pi;
                    protein_atom_charge(i+ip)=-0.5432;
                    protein_backbone(i+ip)=0;
                    
                elseif (strncmp(AtomName{i+ip,1}, 'CD2', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=CW;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.2207;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HD2'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=H4;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.4090;
                    protein_atom_vdwep(i+ip)=0.0150;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1862;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'CE1', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=CR;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1635;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HE1'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=H5;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.3590;
                    protein_atom_vdwep(i+ip)=0.0150;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1435;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'NE2', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=NA;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.2795;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HE2'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=H;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.3590;
                    protein_atom_vdwep(i+ip)=0.0150;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.3339;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'CG', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIE', 3)==1)
                    protein_name(i+ip,1)=CC;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1868;
                    protein_backbone(i+ip)=0;
                end
                
                
                
            elseif Hn1 == 1 && Hn2 == 1
                
                
                if (strcmp(AtomName{i+ip,1}, 'N')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=N;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=pi;
                    protein_atom_charge(i+ip)=-0.3479;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'C')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=C;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.7341;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'O')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=O;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.6612;
                    protein_atom_vdwep(i+ip)=0.2100;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.5894;
                    protein_backbone(i+ip)=1;
                elseif (strncmp(AtomName{i+ip,1}, 'OXT', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=O;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.6612;
                    protein_atom_vdwep(i+ip)=0.2100;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.5679;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'H')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=H;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=0.6000;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.2747;
                    protein_backbone(i+ip)=1;
                    
                    
                elseif (strncmp(AtomName{i+ip,1}, 'CA', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=CT;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.1094;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.1354;
                    protein_backbone(i+ip)=1;
                elseif (size(strfind(AtomName{i+ip,1}, 'HA'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=H1;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.3870;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1212;
                    protein_backbone(i+ip)=1;
                elseif (strncmp(AtomName{i+ip,1}, 'CB', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=CT     ;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.1094;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.0414;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HB'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=HC;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.4870;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.0810;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'ND1', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=NA;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=pi;
                    protein_atom_charge(i+ip)=-0.1513;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HD1'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=H;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.4090;
                    protein_atom_vdwep(i+ip)=0.0150;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.3866;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'CD2', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=CW;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.1141;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HD2'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=H4;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.4090;
                    protein_atom_vdwep(i+ip)=0.0150;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.2317;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'CE1', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=CR;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.017;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HE1'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=H;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.3590;
                    protein_atom_vdwep(i+ip)=0.0150;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.2681;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'NE2', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=NA;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.1718;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HE2'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=H;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.3590;
                    protein_atom_vdwep(i+ip)=0.0150;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.3911;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'CG', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HIP', 3)==1)
                    protein_name(i+ip,1)=CC;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.0012;
                    protein_backbone(i+ip)=0;
                end
                
                
            elseif Hn1 == 0 && Hn2 == 0
                
                
                if (strcmp(AtomName{i+ip,1}, 'N')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=N;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=pi;
                    protein_atom_charge(i+ip)=-0.4157;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'C')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=C;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.5973;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'O')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=O;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.6612;
                    protein_atom_vdwep(i+ip)=0.2100;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.5679;
                    protein_backbone(i+ip)=1;
                elseif (strncmp(AtomName{i+ip,1}, 'OXT', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=O;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.6612;
                    protein_atom_vdwep(i+ip)=0.2100;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.5679;
                    protein_backbone(i+ip)=1;
                elseif (strcmp(AtomName{i+ip,1}, 'H')==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=H;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=0.6000;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.2719;
                    protein_backbone(i+ip)=1;
                    
                elseif (strncmp(AtomName{i+ip,1}, 'CA', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=CT;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.1094;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.0188;
                    protein_backbone(i+ip)=1;
                elseif (size(strfind(AtomName{i+ip,1}, 'HA'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=H1;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.3870;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.0881;
                    protein_backbone(i+ip)=1;
                elseif (strncmp(AtomName{i+ip,1}, 'CB', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=CT     ;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.1094;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.0462;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HB'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=HC;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.4870;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.0402;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'ND1', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=NA;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=pi;
                    protein_atom_charge(i+ip)=-0.3811;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HD1'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=H;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=0.6000;
                    protein_atom_vdwep(i+ip)=0.0157;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.3649;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'CD2', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=CV;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1292;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HD2'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=H4;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.4090;
                    protein_atom_vdwep(i+ip)=0.0150;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1147;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'CE1', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=CR;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.2057;
                    protein_backbone(i+ip)=0;
                elseif (size(strfind(AtomName{i+ip,1}, 'HE1'),2)>0) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=H5;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.3590;
                    protein_atom_vdwep(i+ip)=0.0150;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=0.1392;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'NE2', 3)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=NB;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.8240;
                    protein_atom_vdwep(i+ip)=0.1700;
                    protein_acceptor_hb_angle(i+ip)=0.75*pi;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.5727;
                    protein_backbone(i+ip)=0;
                elseif (strncmp(AtomName{i+ip,1}, 'CG', 2)==1) && (strncmp(resName{i+ip,1}, 'HIS', 3)==1 || strncmp(resName{i+ip,1}, 'HID', 3)==1)
                    protein_name(i+ip,1)=CC;
                    protein_name(i+ip,2)=protein_resseq(i+ip);
                    protein_atom_radius(i+ip)=1.9080;
                    protein_atom_vdwep(i+ip)=0.0860;
                    protein_acceptor_hb_angle(i+ip)=0;
                    protein_donor_hb_angle(i+ip)=0;
                    protein_atom_charge(i+ip)=-0.0266;
                    protein_backbone(i+ip)=0;
                end
                
            end
            
            
            
            
            ip=ip+1;
        end
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0597;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0869;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1303;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0187;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG1', 3)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.043;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HG1'),2)>0) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0236;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'CG2', 3)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.3204;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HG2'),2)>0) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0882;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.066;
        protein_backbone(i)=0;
        
    elseif (size(strfind(AtomName{i,1}, 'HD1'),2)>0) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0186;
        protein_backbone(i)=0;
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0518;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0922;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1102;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4590;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0457;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.3531;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HG'),2)>0) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4590;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0361;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.4121;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HD1'),2)>0) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4590;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1000;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.4121;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HD2'),2)>0) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4590;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1000;
        protein_backbone(i)=0;
        
        
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.7341;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5894;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.3479;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2747;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.24;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1426;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0094;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0362;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NZ', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=N3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.875;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.3854;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'HZ', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.34;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0187;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HG'),2)>0) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0103;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0479;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HD'),2)>0) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0621;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0143;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HE'),2)>0) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1135;
        protein_backbone(i)=0;
        
        
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0237;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0880;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0342;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0241;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0018;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HG'),2)>0) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0440;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}(1), 'S', 1)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=S;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=2.0000;
        protein_atom_vdwep(i)=0.2500;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2737;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE', 2)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0536;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HE'),2)>0) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0684;
        protein_backbone(i)=0;
        
        
        
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0024;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0978;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0342;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0295;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0118;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1256;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1256;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HD'),2)>0) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=HA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4590;
        protein_atom_vdwep(i)=0.0150;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1330;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1704;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE2', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1704;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HE'),2)>0) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=HA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4590;
        protein_atom_vdwep(i)=0.0150;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1430;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CZ', 2)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1072;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HZ'),2)>0) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=HA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4590;
        protein_atom_vdwep(i)=0.0150;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1297;
        protein_backbone(i)=0;
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5896;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5748;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.2548;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0266;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0641;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.007;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0253;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0189;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HG'),2)>0) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0213;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0192;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HD'),2)>0) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0391;
        protein_backbone(i)=0;
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0249;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0843;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2117;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0352;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OG', 2)==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=OH;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.7210;
        protein_atom_vdwep(i)=0.2104;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.6546;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'HG', 2)==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=HO;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0;
        protein_atom_vdwep(i)=0;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.4275;
        protein_backbone(i)=0;
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0389;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1007;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.3654;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0043;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'OG1', 3)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=OH;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.7210;
        protein_atom_vdwep(i)=0.2104;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.6761;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'HG1', 2)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=HO;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0;
        protein_atom_vdwep(i)=0;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.4102;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG2', 3)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2438;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HG2'),2)>0) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0642;
        protein_backbone(i)=0;
        
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0275;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1123;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.005;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0339;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=Cs;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1415;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CW;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1638;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HD1'),2)>0) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=H4;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4090;
        protein_atom_vdwep(i)=0.0150;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2062;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CB;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1243;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NE1', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=NA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.3418;
        protein_backbone(i)=0;
    elseif (strcmp(AtomName{i,1}, 'HE1')==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.3412;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CE2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CN;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.138;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE3', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2387;
        protein_backbone(i)=0;
    elseif (strcmp(AtomName{i,1}, 'HE3')==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1700;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CZ2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2601;
        protein_backbone(i)=0;
    elseif (strcmp(AtomName{i,1}, 'HZ2')==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1572;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CZ3', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1972;
        protein_backbone(i)=0;
    elseif (strcmp(AtomName{i,1}, 'HZ3')==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1447;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CH2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1134;
        protein_backbone(i)=0;
    elseif (strcmp(AtomName{i,1}, 'HH2')==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1417;
        protein_backbone(i)=1;
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0014;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0876;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0152;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0295;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0011;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1906;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1906;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HD'),2)>0) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=HA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4590;
        protein_atom_vdwep(i)=0.0150;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1699;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2341;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE2', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2341;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HE'),2)>0) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=HA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4590;
        protein_atom_vdwep(i)=0.0150;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1656;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CZ', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.3226;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OH', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=OH;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.7210;
        protein_atom_vdwep(i)=0.2104;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.5579;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'HH', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=HO;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0;
        protein_atom_vdwep(i)=0;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.3992;
        protein_backbone(i)=0;
        
        
        
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0875;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0876;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2985;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0339;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG1', 3)==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.3192;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG2', 3)==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.3192;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HG'),2)>0) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0791;
        protein_backbone(i)=0;
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'GLY', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'GLY', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'GLY', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'GLY', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'GLY', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0252;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'GLY', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0698;
        protein_backbone(i)=1;
        
        
        
        
        
        
    elseif (strcmp(AtomName{i,1}, 'C')==1) && (strncmp(resName{i,1}, 'ALA', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.0860;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.5973;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'ALA', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5679;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1) && (strncmp(resName{i,1}, 'ALA', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.8240;
        protein_atom_vdwep(i)=0.1700;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.4157;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'H')==1) && (strncmp(resName{i,1}, 'ALA', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2719;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'ALA', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0337;
        protein_backbone(i)=1;
    elseif (size(strfind(AtomName{i,1}, 'HA'),2)>0) && (strncmp(resName{i,1}, 'ALA', 3)==1)
        protein_name(i,1)=H1;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.3870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0823;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'ALA', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1825;
        protein_backbone(i)=0;
    elseif (size(strfind(AtomName{i,1}, 'HB'),2)>0) && (strncmp(resName{i,1}, 'ALA', 3)==1)
        protein_name(i,1)=HC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.4870;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0603;
        protein_backbone(i)=0;
        
        
        
        
    elseif (strncmp(AtomName{i,1}, 'LPD', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_vdwep(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'LPG', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_vdwep(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
        
        
    elseif (strncmp(AtomName{i,1}, 'CE', 2)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.9080;
        protein_atom_vdwep(i)=0.1094;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.045;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OXT', 3)==1)
        protein_name(i,1)=O ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.6612;
        protein_atom_vdwep(i)=0.2100;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.574;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'HOCA', 4)==1)
        protein_name(i,1)=O ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0;
        protein_atom_vdwep(i)=0;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.4422;
        protein_backbone(i)=1;
        
    elseif (strncmp(AtomName{i,1}, 'ZN', 2)==1)
        protein_name(i,1)=ZN;
        protein_name(i,2)=ZN;
        protein_atom_radius(i)=1.2;
        protein_atom_vdwep(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=2;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'H',1)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=0.6000;
        protein_atom_vdwep(i)=0.0157;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1984;
        protein_backbone(i)=1;
        
    end
    
    %protein_resseq(i)=PDBData1.Atom(1,i).resSeq;
    
end

% collect and group the prepared atoms to define the protein molecule 
clear prepared_protein protein_atom_num protein_H prepared_protein_out;
prepared_protein=[protein_atom,protein_name,resName1,protein_backbone,protein_atom_vdwep,protein_atom_radius,protein_atom_charge];
%prepared_protein=[protein_atom,protein_name,resName1,protein_backbone,protein_atom_vdwep,protein_atom_radius,protein_atom_charge,AtomName];

%protein_Namelist = [AtomName resName];

%protein_list_num = protein_list_final_AF_format(protein_Namelist);

prepared_protein_out = prepared_protein;
prepared_protein_out(prepared_protein(:,4) ==0 | prepared_protein(:,5) ==0,:)=[];
AtomName(prepared_protein(:,4) ==0 | prepared_protein(:,5) ==0,:)=[];

size_prepared_protein=size(prepared_protein_out);

protein_atom_num=1:1:size_prepared_protein(1,1);

prepared_protein_out(:,7)=transpose(protein_atom_num);
