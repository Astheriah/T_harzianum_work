%% Para hacer el blast entre .fasta de interes y .faa .mat de organismos

initCobraToolbox();

folderPath = pwd; 
fileList = dir(folderPath);
fileNames = {fileList.name};
%Definir archivo(s).faa 
fileList = dir(fullfile(folderPath, '*_protein_clean_1.faa'));
cleanFiles = {fileList.name}'; % Lista de todos los archivos encontrados

%Definir .fasta de organismo de interes
list2 = {'THM10.fasta'};

% ==================== Blast con .faa de cleanFiles=============================
for i = 1 : length(list2)
    parfor j = 1 : length(cleanFiles)
        blastStruct = getBlast(char(cleanFiles{j,1}), char(cleanFiles{j,1}), list2{i,1}, list2{i,1});
        innerBlast{i,j} = blastStruct;
    end
end

innerBlast2 = innerBlast;
innerBlast3 = innerBlast2;
HitsInfo = {}; HitsTables = {};

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

        listChar = cellfun(@char, SameGenes2(:,1), 'UniformOutput', false);
        [list, ~, ~] = unique(listChar);
        bestHit = {};
        for i=1:length(list)
            pos = find(strcmp(listChar, list{i,1}));
            identVal = cell2mat(SameGenes2(pos,4));
            [~, idxMax] = max(identVal);
            bestHit(i,:) = SameGenes2(pos(idxMax), :);
        end

        innerBlast3{z,k}(1).fromGenes = cellfun(@char, bestHit(:,1), 'UniformOutput', false);
        innerBlast3{z,k}(1).toGenes = cellfun(@char, bestHit(:,2), 'UniformOutput', false);
        innerBlast3{z,k}(1).evalue = cell2mat(bestHit(:,3));
        innerBlast3{z,k}(1).identity = cell2mat(bestHit(:,4));
        innerBlast3{z,k}(1).aligLen = cell2mat(bestHit(:,5));
        innerBlast3{z,k}(1).bitscore = cell2mat(bestHit(:,6));
        innerBlast3{z,k}(1).ppos = cell2mat(bestHit(:,7));
        HitsInfo{z,k} = bestHit;
        HitsTables{z,k} = cell2table(bestHit, 'VariableNames', {'template ID','model ID','evalue','identity','length','bitscore','ppos'});
    end
end

HitsInfo2 = HitsInfo';
lenHits = cellfun(@length, HitsInfo2);
[M, I] = max(lenHits, [], 1);
[sortVar, idxSort] = sort(lenHits, 'descend');
listOrg = cleanFiles(idxSort);

%save(fullfile('workspaces', 'blast_THM10.mat'));

% ==================== Creación de archivos .fasta a partir de modelos metabólicos (.mat) de cleanFiles2 =======
cleanFiles2 = strrep(cleanFiles, '_protein_clean_1.faa', '.mat');

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

% ==================== Segundo blast con los nuevos fasta =============


cleanFiles3 = strrep(cleanFiles2, '.mat', '.fasta');
list2 = {'THM10.fasta'};


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

HitsInfo3 = HitsInfo';
lenHits = cellfun(@length,HitsInfo2);
[M,I] = max(lenHits,[],1);
[sortVar,idxSort] = sort(lenHits,'descend');
listOrg = cleanFiles3(idxSort);

% ======== generacion de modelos =======

cleanFiles2 = strrep(cleanFiles, '_protein_clean_1.faa', '.mat');

archivos_mat = cleanFiles2;   % es una celda

for i = 1:length(archivos_mat)
    nombre = archivos_mat{i};                        
    ruta = fullfile(folderPath, nombre);
    fprintf('Procesando: %s... ', nombre);

    try
        loaded = load(ruta);
        var = erase(nombre, '.mat');

        if isfield(loaded, var)
            modelo = loaded.(var);
        else
            campos = fieldnames(loaded);
            modelo = loaded.(campos{1});
            var = campos{1};
        end

        if isstruct(modelo)
            if ~isfield(modelo, 'metComps') || isempty(modelo.metComps)
                if isfield(modelo, 'mets') && ~isempty(modelo.mets)
                    if contains(modelo.mets{1}, '[')
                        formato = '[';
                    else
                        formato = '_';
                    end
                    modelo = createMetComp(modelo, formato);
                    fprintf('metComps corregido ');
                else
                    fprintf('sin mets ');
                end
            else
                fprintf('metComps OK ');
            end

            if isfield(modelo, 'description')
                orig = modelo.description;
                [~, modelo.description] = fileparts(orig);
                if ~strcmp(orig, modelo.description)
                    fprintf('description: %s→%s ', orig, modelo.description);
                end
            end

            eval([var ' = modelo;']);
            save(ruta, var);
            fprintf('✅\n');
        else
            fprintf('❌ no es estructura\n');
        end
    catch ME
        fprintf('❌ error: %s\n', ME.message);
    end
end
getModelFor=list2{1,1};
preferredOrder={};
strictness=1;
onlyGenesInModels=false();
mapNewGenesToOld=true();
adjParam1 = [1e-10 100 40];

for i = 1 : length(cleanFiles3)
    model = load(cleanFiles2{i,1});
    model = getfield(model,erase(cleanFiles2{i,1},'.mat'));
    [~, baseName] = fileparts(cleanFiles2{i,1});
    model.id=strcat(model.description,'.fasta');
    model = generateGrRules(model);
    models= {model};
    blastStructureCode = innerBlast3{1,i};
    draftModel=getModelFromHomology(models,blastStructureCode,getModelFor,preferredOrder,strictness,onlyGenesInModels,adjParam1(1,1),adjParam1(1,2),adjParam1(1,3),mapNewGenesToOld);
    draftModels{i,1} = draftModel;
end

draftgetModelFromHomology = draftModels{5, 1};


save(fullfile('workspaces', 'workspaceDraftModels_THM10.mat'));

%==========================================================

metsLength = cellfun(@(x) length(x.mets),draftModels,'UniformOutput',false);
rxnsLength = cellfun(@(x) length(x.rxns),draftModels,'UniformOutput',false);
genesLength = cellfun(@(x) length(x.genes),draftModels,'UniformOutput',false);

statsModels = horzcat(metsLength,rxnsLength,genesLength);

%de draftModels selecionar el modelo plantilla (iAP1008)
% Cambiar variable innerBlast3 por la correspondiente (En este caso corresponde a la 5a)
load("iAP1008.mat")

iAP1008.description=strcat(iAP1008.description,'.fasta');
iAP1008 = generateRules(iAP1008);
optimizeCbModel(iAP1008)
modelTemp1 = reconstructionOpt2(iAP1008,list2{1,1},innerBlast3{1,5},adjParam1(1,1),adjParam1(1,2),adjParam1(1,3));
checkGrowth1 = optimizeCbModel(modelTemp1);
[a,b] = exchangeSingleModel(modelTemp1);
[a1,b1] = exchangeSingleModel(iAP1008);
diffRxns = setdiff(iAP1008.rxns,modelTemp1.rxns);
diffRxns(:,2) = iAP1008.grRules(findRxnIDs(iAP1008,diffRxns(:,1)));

%Cambiar Hitsinfo2 por la correspondiente (5a)
blastKleb = HitsInfo2{5,1};

blastKleb2 = cell(size(blastKleb));  
for i = 1:size(blastKleb, 1)
    for j = 1:size(blastKleb, 2)
        valor = blastKleb{i, j};
        while iscell(valor) && isscalar(valor) && ~isempty(valor)
            valor = valor{1};
        end
        blastKleb2{i, j} = valor;
    end
end


load("iAP1008.mat")
modelTemp2 = iAP1008;
modelTemp2 = changeRxnBounds(modelTemp2,'EX_glc__D_e',-5,'l');
checkGrowth2 = optimizeCbModel(modelTemp2);

%======= Remplazo de ===============
% for i = 1: length(blastKleb)
%    modelTemp2.grRules = strrep(modelTemp2.grRules,blastKleb{i,2},blastKleb{i,1});
%end

if istable(blastKleb2)
    oldIDs = cellstr(blastKleb2{:,2});
    newIDs = cellstr(blastKleb2{:,1});
else
    oldIDs = blastKleb2(:,2);
    newIDs = blastKleb2(:,1);
end
% Crear marcadores temporales únicos
n = numel(oldIDs);
placeholders = cell(n,1);
for i = 1:n
    placeholders{i} = sprintf('@@%d@@', i);  % por ejemplo @@1@@
end

for i = 1:n
    pattern = ['(?<![A-Za-z0-9_])' regexptranslate('escape', oldIDs{i}) '(?![A-Za-z0-9_])'];
    modelTemp2.grRules = regexprep(modelTemp2.grRules, pattern, newIDs{i});
end
%==================================================================

modelTemp2 =rmfield(modelTemp2,'rules');
modelTemp2 = buildRxnGeneMat(modelTemp2);
modelTemp2 = generateRules(modelTemp2);
modelTemp2 = removeUnusedGenes(modelTemp2);

genesModelTemp2 = table(modelTemp2.genes, 'VariableNames', {'GeneID'});
genesModelTemp2(1:176, :) = [];
proteinIDs = regexp(genesModelTemp2.GeneID, '(?<=_)\d+', 'match', 'once');
genesModelTemp2.ProteinID = proteinIDs;
writetable(genesModelTemp2, 'genesModelTemp2.csv');


%let's find the reactions with APT genes
genesKPN = modelTemp2.genes(find(contains(modelTemp2.genes,'APT')));

checkGPR = modelTemp2.grRules(find(contains(modelTemp2.grRules,'APT')));

%Poner rxns para eliminar las cuales dentro de checkGPR solo tengas genes
%de iAP1008 es decir que sean genes APT_X

%Para identificar que rxns poner en rxnsDel
idxKPN = find(contains(modelTemp2.grRules, 'APT'));
for i = 1:length(checkGPR)
    idxGlobal = idxKPN(i);
    fprintf('%d. Indice: %d | Reaccion: %s | GPR: %s\n', i, idxGlobal, modelTemp2.rxns{idxGlobal}, checkGPR{i});
end

%CORRER S12_API y obtener ProteinID_EC_info.csv

genesTemp2 = table(modelTemp2.genes, 'VariableNames', {'GeneID'});

genesTemp2Info = readtable('ProteinID_EC_info.csv');

genesTemp2.ProteinID = regexp(genesTemp2.GeneID, '\d+$', 'match', 'once');

if isnumeric(genesTemp2Info.ProteinID)
    genesTemp2Info.ProteinID = string(genesTemp2Info.ProteinID);
end
genesTemp2Info.ProteinID = regexp(genesTemp2Info.ProteinID, '\d+$', 'match', 'once');

colsExtra = setdiff(genesTemp2Info.Properties.VariableNames, 'ProteinID');

for c = 1:numel(colsExtra)
    col = colsExtra{c};
    if isnumeric(genesTemp2Info.(col))
        genesTemp2.(col) = NaN(height(genesTemp2), 1);
    else
        genesTemp2.(col) = repmat({''}, height(genesTemp2), 1);
    end
end

[~, idx] = ismember(genesTemp2.ProteinID, genesTemp2Info.ProteinID);
for c = 1:numel(colsExtra)
    col = colsExtra{c};
    mask = idx > 0;
    genesTemp2.(col)(mask) = genesTemp2Info.(col)(idx(mask));
end


save(fullfile('workspaces', 'workspaceCONTROL.mat'))



















rxnsDel = {'HDAO10x';'GCLDH'};

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

%save(fullfile('workspaces', 'workspaceCONTROL3.mat'))

%fluxTemp = optimizeCbModel(modelTemp2,'max','one',0);


%let's generate the possible set of reactions to add.
%modelTemp2 will be the initial draft model

draftModels2 = draftModels;          
draftModels2(5) = [];

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

writecell(poolRxns,'poolRxnsCuration.csv')

protData = fastaread('THM10.fasta');

protHeader = length({protData.Header}.');


save(fullfile('workspaces', 'workspaceCONTROL4.mat'))

