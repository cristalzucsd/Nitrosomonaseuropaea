%% Clustergram creation

load('iGC535_HCO3High.mat');
model=iGC535_HCO3High;
load('newSubsystems.mat');
%newSubsystems{1150,1}='Other'
%model=iGC535Plos;

%Low Fructose
model=changeRxnBounds(model,'EX_hco3_e',0,'l');
model=changeRxnBounds(model,'EX_fru_e',-0.014,'l');
model=changeRxnBounds(model,'EX_nh4_e',-0.5,'l');
Fru1=model;

%High Fructose
model=changeRxnBounds(model,'EX_fru_e',-0.746,'l');
Fru2=model;

%Pyruvate
model=changeRxnBounds(model,'EX_fru_e',0,'l');
model=changeRxnBounds(model,'EX_nh4_e',-1000,'l');
model=changeRxnBounds(model,'EX_pyr_e',-0.5,'l');
Pyr=model;

%High HCO3
model=changeRxnBounds(model,'EX_pyr_e',0,'l');
model=changeRxnBounds(model,'EX_hco3_e',-4.3524,'l');
model=changeRxnBounds(model,'EX_nh4_e',-1000,'l');
HCO32=model;

%Low HCO3
model=changeRxnBounds(model,'EX_pyr_e',0,'l');
model=changeRxnBounds(model,'EX_hco3_e',-1.4363,'l');
model=changeRxnBounds(model,'EX_nh4_e',-1000,'l');
HCO31=model;

%HCO3+Methane
model=changeRxnBounds(model,'EX_ch4_e',-1,'l');
model=changeRxnBounds(model,'EX_ch4_e',-1,'u');
ch4=model;

%HCO3+Toluene
model=changeRxnBounds(model,'EX_ch4_e',0,'l');
model=changeRxnBounds(model,'EX_ch4_e',1000,'u');
model=changeRxnBounds(model,'EX_tol_e',-1,'l');
model=changeRxnBounds(model,'EX_tol_e',-1,'u');
toluene=model;

%HCO3+Benzene
model=changeRxnBounds(model,'EX_benz_e',-1,'u');
model=changeRxnBounds(model,'EX_benz_e',-1,'l');
model=changeRxnBounds(model,'EX_tol_e',0,'l');
model=changeRxnBounds(model,'EX_tol_e',1000,'u');
Benzene=model;

%HCO3+Chlorobenzene
model=changeRxnBounds(model,'EX_chlben_e',-1,'u');
model=changeRxnBounds(model,'EX_chlben_e',-1,'l');
model=changeRxnBounds(model,'EX_benz_e',1000,'u');
model=changeRxnBounds(model,'EX_benz_e',0,'l');
Chlorobenzene=model;

%HCO3+Phenol
model=changeRxnBounds(model,'EX_chlben_e',1000,'u');
model=changeRxnBounds(model,'EX_chlben_e',-1000,'l');
model=changeRxnBounds(model,'EX_phenol_e',-1,'u');
model=changeRxnBounds(model,'EX_phenol_e',-1,'l');
Phenol=model;


M1{1,1}=HCO31;
M1{2,1}=HCO32;
M1{3,1}=Fru1;
M1{4,1}=Fru2;
M1{5,1}=Pyr;
M1{6,1}=ch4;
M1{7,1}=toluene;
M1{8,1}=Benzene;
M1{9,1}=Phenol;
M1{10,1}=Chlorobenzene;

%Models benchmarking

FBA=[];
for i=1:numel(M1)
    model=M1{i,1};
    Growth(i,1)=optimizeCbModel(model,'max','one');
    Rate(i,1)=Growth(i,1).f;
    FBA(:,i)=Growth(i,1).x;
end
FBA2=[];
for i=1:length(FBA(1,:))
    for j=1:length(FBA(:,i))
        FBA2(j,i)=FBA(j,i);
    end
end
subsys={};
newFBA=[];
rxnsZ={};
for i=1:length(FBA2(:,1))
    row=FBA2(i,:);
    row(~any(row,2),:)=[];
    if ~isempty(row)
        hitEX=strmatch('Exchange reactions',newSubsystems{i},'exact');
        if isempty(hitEX)
            subsys{end+1,1}=newSubsystems{i};
            rxnsZ{end+1,1}=HCO31.rxns{i};
            Z = zscore(row);
            newFBA(end+1,:)=Z;
        end
    end
end

cgo=clustergram(newFBA,'Symmetric',0,'Colormap',redbluecmap)
set(cgo,'ColumnLabels',{'HCO3','HCO32','Fructose1','Fructose2','Pyruvate','Methane','Toluene','Benzene','Phenol','Chlorobenzene'});
