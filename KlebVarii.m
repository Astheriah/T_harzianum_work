initCobraToolbox
changeCobraSolver('gurobi','all')
cd /mnt/sda1/KlebsiellaNitrogen/
addpath('/mnt/sda1/KlebsiellaNitrogen/')


%Klebsiella variicola diazotroph

% fastaKleb = cleanFasta({'KlebsiellaVariicolaNF.txt'},{'NA'});
% fastaKleb2 =fastaKleb{1,1};
% headerK = {fastaKleb2.Header}.';
% seqK = {fastaKleb2.Sequence}.';
% 
% 
% 
% fastawrite('klebsiellaClean.fasta',headerK,seqK);


folderPath = '/mnt/sda1/KlebsiellaNitrogen/'; 
fileList = dir(folderPath);
fileNames = {fileList.name};
cleanFiles = fileNames(contains(fileNames, '.faa'))';


list2 = {'klebsiellaClean.fasta'};


for i = 1 : length(list2)
    parfor j = 1 : length(cleanFiles)
        blastStruct = getBlast(char(cleanFiles{j,1}), char(cleanFiles{j,1}), list2{i,1}, list2{i,1});
        innerBlast{i,j} = blastStruct;
    end
end



innerBlast2 = innerBlast;
innerBlast3 = innerBlast2;
HitsInfo = {};HitsTables = {};
for k=1:size(innerBlast2,2)
    k
    for z = 1 : size(innerBlast2,1)
        BLAST = innerBlast2{z,k};
        parameters = [1e-5 20 30];
        cont2 = 1;
        SameGenes2={};
        for j=1:length(BLAST(1).evalue)
               if BLAST(1).evalue(j,1)<=parameters(1,1) && BLAST(1).aligLen(j,1)>=parameters(1,2) && BLAST(1).identity(j,1)>=parameters(1,3)
                    SameGenes2{cont2,1}=BLAST(1).fromGenes(j,1);
                    SameGenes2{cont2,2}=BLAST(1).toGenes(j,1);
                    SameGenes2{cont2,3}=BLAST(1).evalue(j,1);
                    SameGenes2{cont2,4}=BLAST(1).identity(j,1);
                    SameGenes2{cont2,5}=BLAST(1).aligLen(j,1);
                    SameGenes2{cont2,6}=BLAST(1).bitscore(j,1);
                    SameGenes2{cont2,7}=BLAST(1).ppos(j,1);
                    cont2=cont2+1;
               end
        
        end
        
        listChar=cellfun(@char,SameGenes2(:,1),'UniformOutput',false);
        [list,A,C]=unique(listChar);
        bestHit = {};
        for i=1:length(list)
            pos=find(strcmp(listChar,list{i,1}));
            identVal=cell2mat(SameGenes2(pos,4));
            proof=SameGenes2(pos,1);
            [getMax,idxMax]=max(identVal);
            bestHit(i,:)=SameGenes2(pos(idxMax),:);
            
        end
    
        innerBlast3{z,k}(1).fromGenes = cellfun(@char,bestHit(:,1),'UniformOutput',false);
        innerBlast3{z,k}(1).toGenes = cellfun(@char,bestHit(:,2),'UniformOutput',false);
        innerBlast3{z,k}(1).evalue = cell2mat(bestHit(:,3));
        innerBlast3{z,k}(1).identity = cell2mat(bestHit(:,4));
        innerBlast3{z,k}(1).aligLen = cell2mat(bestHit(:,5));
        innerBlast3{z,k}(1).bitscore = cell2mat(bestHit(:,6));
        innerBlast3{z,k}(1).ppos = cell2mat(bestHit(:,7));
        HitsInfo{z,k} = bestHit;
        HitsTables{z,k} = cell2table(bestHit,'VariableNames',{'template ID','model ID','evalue','identity','length','bitscore','ppos'});

    end
end

HitsInfo2 = HitsInfo';
lenHits = cellfun(@length,HitsInfo2);
[M,I] = max(lenHits,[],1);
[sortVar,idxSort] = sort(lenHits,'descend');
listOrg = cleanFiles(idxSort);


%link models
folderPath = '/mnt/sda1/KlebsiellaNitrogen/'; 
fileList = dir(folderPath);
fileNames = {fileList.name};
cleanFiles2 = fileNames(contains(fileNames, '.mat'))';
ignoreList= {'workspaceBlastInitial.mat';'workspaceDraftModels.mat';'workspaceCurationKlebsiella.mat';'workspaceCurationKlebsiellaBiologPlates.mat';...
    'workspaceKlebsiellaRxnsAdded.mat';'workspaceModelCuratedMets.mat'};
cleanFiles2 = cleanFiles2(~ismember(cleanFiles2, ignoreList));

for i = 1 : length(cleanFiles2)
    model = load(cleanFiles2{i,1});
    model = getfield(model,erase(cleanFiles2{i,1},'.mat'));
    genes = model.genes;
    fastaData = fastaread(cleanFiles{i,1});
    headerFasta = {fastaData.Header}.';
    seqFasta = {fastaData.Sequence}.';
    headerFasta = strtok(headerFasta);
    [commonGenes,idxGene] = intersect(headerFasta,genes);
    fastawrite(strcat(erase(cleanFiles2{i,1},'.mat'),'.','fasta'),headerFasta(idxGene),seqFasta(idxGene));

end



folderPath = '/mnt/sda1/KlebsiellaNitrogen/'; 
fileList = dir(folderPath);
fileNames = {fileList.name};
cleanFiles3 = fileNames(contains(fileNames, '.fasta'))';
ignoreList= {'klebsiellaClean.fasta'};
cleanFiles3 = cleanFiles3(~ismember(cleanFiles3, ignoreList));


list2 = {'klebsiellaClean.fasta'};

clear innerBlast
for i = 1 : length(list2)
    parfor j = 1 : length(cleanFiles3)
        blastStruct = getBlast(char(cleanFiles3{j,1}), char(cleanFiles3{j,1}), list2{i,1}, list2{i,1});
        innerBlast{i,j} = blastStruct;
    end
end



innerBlast2 = innerBlast;
innerBlast3 = innerBlast2;
HitsInfo = {};HitsTables = {};
for k=1:size(innerBlast2,2)
    k
    for z = 1 : size(innerBlast2,1)
        BLAST = innerBlast2{z,k};
        parameters = [1e-5 20 30];
        cont2 = 1;
        SameGenes2={};
        for j=1:length(BLAST(1).evalue)
               if BLAST(1).evalue(j,1)<=parameters(1,1) && BLAST(1).aligLen(j,1)>=parameters(1,2) && BLAST(1).identity(j,1)>=parameters(1,3)
                    SameGenes2{cont2,1}=BLAST(1).fromGenes(j,1);
                    SameGenes2{cont2,2}=BLAST(1).toGenes(j,1);
                    SameGenes2{cont2,3}=BLAST(1).evalue(j,1);
                    SameGenes2{cont2,4}=BLAST(1).identity(j,1);
                    SameGenes2{cont2,5}=BLAST(1).aligLen(j,1);
                    SameGenes2{cont2,6}=BLAST(1).bitscore(j,1);
                    SameGenes2{cont2,7}=BLAST(1).ppos(j,1);
                    cont2=cont2+1;
               end
        
        end
        
        listChar=cellfun(@char,SameGenes2(:,1),'UniformOutput',false);
        [list,A,C]=unique(listChar);
        bestHit = {};
        for i=1:length(list)
            pos=find(strcmp(listChar,list{i,1}));
            identVal=cell2mat(SameGenes2(pos,4));
            proof=SameGenes2(pos,1);
            [getMax,idxMax]=max(identVal);
            bestHit(i,:)=SameGenes2(pos(idxMax),:);
            
        end
    
        innerBlast3{z,k}(1).fromGenes = cellfun(@char,bestHit(:,1),'UniformOutput',false);
        innerBlast3{z,k}(1).toGenes = cellfun(@char,bestHit(:,2),'UniformOutput',false);
        innerBlast3{z,k}(1).evalue = cell2mat(bestHit(:,3));
        innerBlast3{z,k}(1).identity = cell2mat(bestHit(:,4));
        innerBlast3{z,k}(1).aligLen = cell2mat(bestHit(:,5));
        innerBlast3{z,k}(1).bitscore = cell2mat(bestHit(:,6));
        innerBlast3{z,k}(1).ppos = cell2mat(bestHit(:,7));
        HitsInfo{z,k} = bestHit;
        HitsTables{z,k} = cell2table(bestHit,'VariableNames',{'template ID','model ID','evalue','identity','length','bitscore','ppos'});

    end
end

HitsInfo2 = HitsInfo';
lenHits = cellfun(@length,HitsInfo2);
[M,I] = max(lenHits,[],1);
[sortVar,idxSort] = sort(lenHits,'descend');
listOrg = cleanFiles3(idxSort);


%lets generate template models
folderPath = '/mnt/sda1/KlebsiellaNitrogen/'; 
fileList = dir(folderPath);
fileNames = {fileList.name};
cleanFiles2 = fileNames(contains(fileNames, '.mat'))';
ignoreList= {'workspaceBlastInitial.mat'};
cleanFiles2 = cleanFiles2(~ismember(cleanFiles2, ignoreList));


%fix iML1515
load('iML1515.mat')
tempModel = createMetComp(iML1515,'_');
iML1515 = tempModel;
save('iML1515.mat','iML1515')


%fix iYL1228
load('iYL1228.mat')
tempModel = createMetComp(iYL1228,'_');
iYL1228 = tempModel;
save('iYL1228.mat','iYL1228')


getModelFor=list2{1,1};
preferredOrder={};
strictness=1;
onlyGenesInModels=false();
mapNewGenesToOld=true();
adjParam1 = [1e-10 100 40];

for i = 1 : length(cleanFiles3)
    model = load(cleanFiles2{i,1});
    model = getfield(model,erase(cleanFiles2{i,1},'.mat'));
    model.id=strcat(model.description,'.fasta');
    models= {model};
    blastStructureCode = innerBlast3{1,i};
    draftModel=getModelFromHomology(models,blastStructureCode,getModelFor,preferredOrder,strictness,onlyGenesInModels,adjParam1(1,1),adjParam1(1,2),adjParam1(1,3),mapNewGenesToOld);
    draftModels{i,1} = draftModel;
end

%save('workspaceDraftModels.mat')

metsLength = cellfun(@(x) length(x.mets),draftModels,'UniformOutput',false);
rxnsLength = cellfun(@(x) length(x.rxns),draftModels,'UniformOutput',false);
genesLength = cellfun(@(x) length(x.genes),draftModels,'UniformOutput',false);

statsModels = horzcat(metsLength,rxnsLength,genesLength);

iYL1228.description=strcat(iYL1228.description,'.fasta');
iYL1228 = generateRules(iYL1228);
optimizeCbModel(iYL1228)
modelTemp1 = reconstructionOpt2(iYL1228,list2{1,1},innerBlast3{1,19},adjParam1(1,1),adjParam1(1,2),adjParam1(1,3));
checkGrowth1 = optimizeCbModel(modelTemp1);
[a,b] = exchangeSingleModel(modelTemp1);
[a1,b1] = exchangeSingleModel(iYL1228);
diffRxns = setdiff(iYL1228.rxns,modelTemp1.rxns);
diffRxns(:,2) = iYL1228.grRules(findRxnIDs(iYL1228,diffRxns(:,1)));


blastKleb = HitsInfo2{19,1};

modelTemp2 = iYL1228;
modelTemp2 = changeRxnBounds(modelTemp2,'EX_glc__D_e',-5,'l');
checkGrowth2 = optimizeCbModel(modelTemp2);

for i = 1: length(blastKleb)
    modelTemp2.grRules = strrep(modelTemp2.grRules,blastKleb{i,2},blastKleb{i,1});

end


modelTemp2 =rmfield(modelTemp2,'rules');
modelTemp2 = buildRxnGeneMat(modelTemp2);
modelTemp2 = generateRules(modelTemp2);
modelTemp2 = removeUnusedGenes(modelTemp2);



%let's find the reactions with KPN genes
genesKPN = modelTemp2.genes(find(contains(modelTemp2.genes,'KPN')));

checkGPR = modelTemp2.grRules(find(contains(modelTemp2.grRules,'KPN')));

rxnsDel = {'ALLULPE';'DAPAL';'RBK_Dr';'Rbtt2';'TDPDRE';'TDPDRR';'XYLt2pp'};

modelTemp3 = removeRxns(modelTemp2,rxnsDel);
for i=1:length(modelTemp3.mets)
    checkMet=findRxnsFromMets(modelTemp3,modelTemp3.mets(i));
    if isempty(checkMet)
        listDelete{cont,1}=modelTemp3.mets{i,1};
        cont=cont+1;
        
    end
end

modelTemp3.csense={};

for j=1:length(modelTemp3.mets)
    modelTemp3.csense{j}='E';
end

 modelTemp3.csense=char(modelTemp3.csense);
[a2,b2] = exchangeSingleModel(modelTemp2);

fluxTemp = optimizeCbModel(modelTemp2,'max','one',0);


%let's generate the possible set of reactions to add.
%modelTemp2 will be the initial draft model
draftModels2 = draftModels(1:end-1);
initRxns = modelTemp2.rxns;

for i = 1 : length(draftModels2)
    if i == 1
        [rxnsDif,idxDif] = setdiff(draftModels2{i,1}.rxns,initRxns);
        rxnsDif(:,2) = draftModels2{i,1}.rxnNames(findRxnIDs(draftModels2{i,1},rxnsDif));
        rxnsDif(:,3) = printRxnFormula(draftModels2{i,1},rxnsDif(:,1));
        rxnsDif(:,4) = draftModels2{i,1}.grRules(findRxnIDs(draftModels2{i,1},rxnsDif(:,1)));
        rxnsDif(:,5) = {cleanFiles2{i,1}};
        poolRxns = rxnsDif;
    
    else
        extraRxns = setdiff(draftModels2{i,1}.rxns,poolRxns(:,1));
        [rxnsDif,idxDif] = setdiff(extraRxns(:,1),initRxns);
        rxnsDif(:,2) = draftModels2{i,1}.rxnNames(findRxnIDs(draftModels2{i,1},rxnsDif));
        rxnsDif(:,3) = printRxnFormula(draftModels2{i,1},rxnsDif(:,1));
        rxnsDif(:,4) = draftModels2{i,1}.grRules(findRxnIDs(draftModels2{i,1},rxnsDif(:,1)));
        rxnsDif(:,5) = {cleanFiles2{i,1}};
        poolRxns = vertcat(poolRxns,rxnsDif);
    
    end
    
end

writecell(poolRxns,'poolRxnsCurationKVariicola.csv')

protData = fastaread('klebsiellaClean.fasta');

protHeader = length({protData.Header}.');


%once the curation process is ready, let's add the reaction 

updatedRules = readcell("klebsiella_curationUpdated.xlsx",'Sheet','Sheet1');
mask = cellfun(@ismissing, updatedRules(:,3),'UniformOutput',false); % Use @isnan for only NaN values
mask2 =  cellfun(@length,mask);
mask2 = find(mask2==1);
% Replace missing entries with an empty character vector
updatedRules(mask2,3) = {''}; 

modelTemp3 = modelTemp2;
modelTemp3.grRules = updatedRules(:,3);
modelTemp3 = rmfield(modelTemp3,'rules');
modelTemp3 = rmfield(modelTemp3,'rxnGeneMat');
modelTemp3 = rmfield(modelTemp3,'genes');
modelTemp3 = generateRules(modelTemp3);
modelTemp3 = updateGenes(modelTemp3);

uniqSubs = unique(modelTemp3.subSystems);

rxnsAdd = readcell("rxnsCuratedKlebVariicola2.xlsx",'Sheet','poolRxnsCurationKVariicola');
subsDict = readcell("subsystemsKlebsiella.xlsx",'Sheet','Sheet1');

for i = 1 : length(subsDict)
    getPos = strmatch(subsDict{i,1},modelTemp3.subSystems,'exact');
    modelTemp3.subSystems(getPos) = subsDict(i,2);
end

uniqSubs2 = unique(modelTemp3.subSystems);

redIdx = cell2mat(rxnsAdd(:,6));
redRxns = rxnsAdd(find(redIdx),:);

modelTemp4 = modelTemp3;
fluxTemp = {};
for i = 1 : length(redRxns)
    i
    modelTemp4 = addReaction(modelTemp4,redRxns{i,1},'reactionName',redRxns{i,2},'reactionFormula',redRxns{i,3},'subSystem',redRxns{i,7},'geneRule',redRxns{i,4});
    growthF = optimizeCbModel(modelTemp4);
    growthF.f
    [~,b1] = exchangeSingleModel(modelTemp4);
    fluxTemp(:,i) = table2cell(b1(:,end));
end

[a1,b1] = exchangeSingleModel(modelTemp4);

modelTemp4 = changeRxnBounds(modelTemp4,'HYD4pp',0,'b');
modelTemp4.csense={};
for j=1:length(modelTemp4.mets)
    modelTemp4.csense{j}='E';
end
modelTemp4.csense=char(modelTemp4.csense);

growth2 = optimizeCbModel(modelTemp4,'max','one',0);

metsToAdd = setdiff(modelTemp4.mets, modelTemp3.mets);
%writecell(metsToAdd,'metsCuration.csv')

%add curated metabolites
metsCurated = readcell('metsCuration.xlsx','Sheet','metsCuration');

modelTemp4 = modelTemp3;
for i = 1 : length(metsCurated)
    modelTemp4 = addMetabolite(modelTemp4,metsCurated{i,1},'metName',metsCurated{i,2},'metFormula',metsCurated{i,3},'Charge',metsCurated{i,4});
end

%add reactions again once the metabolites have been added
fluxTemp = {};
for i = 1 : length(redRxns)
    i
    modelTemp4 = addReaction(modelTemp4,redRxns{i,1},'reactionName',redRxns{i,2},'reactionFormula',redRxns{i,3},'subSystem',redRxns{i,7},'geneRule',redRxns{i,4});
    growthF = optimizeCbModel(modelTemp4);
    %growthF.f
    [~,b1] = exchangeSingleModel(modelTemp4);
    fluxTemp(:,i) = table2cell(b1(:,end));
end

[a1,b1] = exchangeSingleModel(modelTemp4);

modelTemp4 = changeRxnBounds(modelTemp4,'HYD4pp',0,'b');
modelTemp4.csense={};
for j=1:length(modelTemp4.mets)
    modelTemp4.csense{j}='E';
end
modelTemp4.csense=char(modelTemp4.csense);

growth2 = optimizeCbModel(modelTemp4,'max','one',0);


%lets include the biolog plates results PM1 and PM2. The experimental growth values have been calculated using a 3 standard deviation threshold. The values are included in the ons to add are included in the files:
%for PM1: "PM1_csv_binary_detailed_new.csv"
%for PM2: "PM2_csv_binary_detailed_new.csv"

PM1 = readcell("PM1_csv_binary_detailed_new.csv");
PM2 = readcell("PM2_csv_binary_detailed_new.csv");

metsKleb = modelTemp4.mets;
metsKleb(:,2)= modelTemp4.metNames;
writecell(metsKleb,'metsKlebsiella.csv')
rxnsKleb = modelTemp4.rxns;
rxnsKleb(:,2) = modelTemp4.rxnNames;
writecell(rxnsKleb,'rxnsKlebsiella.csv')

tablePM1 = readtable("PM1_csv_binary_detailed_new2.csv");
tablePM2 = readtable("PM2_csv_binary_detailed_new2.csv");

tablePM1.mapping(43) = {'2obut_e'};

[modelUpdated1, reportTbl] = processMappedMetabolites(modelTemp4, tablePM1, 'mapping');

writetable(reportTbl, '/mnt/sda1/KlebsiellaNitrogen/metabolite_mapping_exchange_reportPM1.csv');

[modelUpdated2, reportTbl2] = processMappedMetabolites(modelUpdated1, tablePM2, 'mapping');
writetable(reportTbl2, '/mnt/sda1/KlebsiellaNitrogen/metabolite_mapping_exchange_reportPM2.csv');


[summaryPM1, catsPM1] = summarizeMetaboliteReport('metabolite_mapping_exchange_reportPM1.csv', 'PM1');

[summaryPM2, catsPM2] = summarizeMetaboliteReport('metabolite_mapping_exchange_reportPM2.csv', 'PM2');

T = readtable('metabolite_mapping_exchange_reportPM2.csv', 'TextType', 'string');

idx_missing = T.Notes == ...
    "Mapped metabolite not found in extracellular, cytosol, or periplasm";

missingTbl = T(idx_missing, :);

writetable(missingTbl, 'PM2_model_expansion_targets.csv');

[summaryTbl, detailedTbl] = generatePreSimulationStatusSummary( ...
    'metabolite_mapping_exchange_reportPM1.csv', ...
    'metabolite_mapping_exchange_reportPM2.csv', ...
    'PM_preSimulation_status');

%PM2 summary of decisions and priorities

decPM2 =classifyCurationPrioritiesByRow( ...
    'metabolite_mapping_exchange_reportPM2.csv', ...
    'PM2_csv_binary_detailed_new2.csv', ...
    'PM2_curation_decisions_fixed.csv', ...
    true);

sumPM2 = summarizeCurationDecisions(decPM2, 'PM2_curation_decision_summary.csv');

%PM1 summary of decisions and priorities

decPM1 =classifyCurationPrioritiesByRow( ...
    'metabolite_mapping_exchange_reportPM1.csv', ...
    'PM1_csv_binary_detailed_new2.csv', ...
    'PM1_curation_decisions_fixed.csv', ...
    true);
sumPM1 = summarizeCurationDecisions(decPM1, 'PM1_curation_decision_summary.csv');


save('workspaceCurationKlebsiellaBiologPlates.mat')

%After Biolog plates validation, we can proceed to gapfilling and finalizing the model. We can also use the PM1 and PM2 results to further curate the model, especially in terms of the exchange reactions and the associated genes.

%lets run the simple check Biolog plates for PM1 and PM2 conditions with the current modelTemp4. We will use the same growth threshold as for the experimental data to determine if the model predicts growth or not.

%first PM1
modelPM1 = modelTemp4;
[pm1_1,pm1_2] = exchangeSingleModel(modelPM1);
modelPM1 = changeRxnBounds(modelPM1,'EX_glc__D_e',0,'l');

for i = 1 : height(tablePM1)
    modelPM1 = modelTemp4;
    modelPM1 = changeRxnBounds(modelPM1,'EX_glc__D_e',0,'l');
    exchRxn = strcat('EX_',tablePM1.mapping{i});
    exchCheck = findRxnIDs(modelPM1,exchRxn);
    if exchCheck ~=0
        modelPM1 = changeRxnBounds(modelPM1,exchRxn,-5,'l');
        growthF = optimizeCbModel(modelPM1,'max','one',0);
        if growthF.f > 0.01
            tablePM1.predictedGrowth(i) = 1;
        else
            tablePM1.predictedGrowth(i) = 0;
        end
    else
        tablePM1.predictedGrowth(i) = 0;
    end
   
end

writetable(tablePM1,'tablePM1_initial_predictions.csv')


%then PM2
modelPM2 = modelTemp4;
[pm2_1,pm2_2] = exchangeSingleModel(modelPM2);
modelPM2 = changeRxnBounds(modelPM2,'EX_glc__D_e',0,'l');

for i = 1 : height(tablePM2)
    modelPM2 = modelTemp4;
    modelPM2 = changeRxnBounds(modelPM2,'EX_glc__D_e',0,'l');
    exchRxn = strcat('EX_',tablePM2.mapping{i});
    exchCheck = findRxnIDs(modelPM2,exchRxn);
    if exchCheck ~=0
        modelPM2 = changeRxnBounds(modelPM2,exchRxn,-5,'l');
        growthF = optimizeCbModel(modelPM2,'max','one',0);
        if growthF.f > 0.01
            tablePM2.predictedGrowth(i) = 1;
        else
            tablePM2.predictedGrowth(i) = 0;
        end
    else
        tablePM2.predictedGrowth(i) = 0;
    end
   
end

writetable(tablePM2,'tablePM2_initial_predictions.csv')


%PM1 classification 

[Tpm1, summaryPM1] = classifyPredictionOutcomes( ...
    'tablePM1_initial_predictions.csv', 'PM1');


%PM2 classification

[Tpm2, summaryPM2] = classifyPredictionOutcomes( ...
    'tablePM2_initial_predictions.csv', 'PM2');

