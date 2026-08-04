%% Para hacer el blast entre .fasta de interes y .faa .mat de organismos
%% plantilla

initCobraToolbox();

folderPath = pwd; 
fileList = dir(folderPath);
fileNames = {fileList.name};
%Definir archivo(s).faa 
cleanFiles = {'iAP1008_protein_clean_1.faa'};
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

% ==================== Creación de archivos .fasta a partir de modelos metabólicos (.mat) de cleanFiles2 =======

cleanFiles2 = {'iAP1008.mat'};

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

cleanFiles3 = {'iAP1008.fasta'};
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

cleanFiles2 = {'iAP1008.mat'};
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

% Guardar todo
save('DraftModelsiAP1008_THM10.mat')

%===================================================

metsLength = cellfun(@(x) length(x.mets),draftModels,'UniformOutput',false);
rxnsLength = cellfun(@(x) length(x.rxns),draftModels,'UniformOutput',false);
genesLength = cellfun(@(x) length(x.genes),draftModels,'UniformOutput',false);

statsModels = horzcat(metsLength,rxnsLength,genesLength);



iAP1008.description=strcat(iAP1008.description,'.fasta');
iAP1008.id = 'iAP1008.fasta'; 
iAP1008 = generateRules(iAP1008);
optimizeCbModel(iAP1008)
modelTemp1 = reconstructionOpt2(iAP1008,list2{1,1},innerBlast3{1,1},adjParam1(1,1),adjParam1(1,2),adjParam1(1,3));
checkGrowth1 = optimizeCbModel(modelTemp1);
[a,b] = exchangeSingleModel(modelTemp1);
[a1,b1] = exchangeSingleModel(iAP1008);
diffRxns = setdiff(iAP1008.rxns,modelTemp1.rxns);
diffRxns(:,2) = iAP1008.grRules(findRxnIDs(iAP1008,diffRxns(:,1)));



blastKleb = HitsInfo2{1,1};

modelTemp2 = iAP1008;
modelTemp2 = changeRxnBounds(modelTemp2,'EX_glc__D_e',-5,'l');
checkGrowth2 = optimizeCbModel(modelTemp2);

for i = 1: length(blastKleb)
    modelTemp2.grRules = strrep(modelTemp2.grRules,blastKleb{i,2},blastKleb{i,1});
end

modelTemp2 =rmfield(modelTemp2,'rules');
modelTemp2 = buildRxnGeneMat(modelTemp2);
modelTemp2 = generateRules(modelTemp2);
modelTemp2 = removeUnusedGenes(modelTemp2);



%let's find the reactions with APT genes
genesKPN = modelTemp2.genes(find(contains(modelTemp2.genes,'APT')));

checkGPR = modelTemp2.grRules(find(contains(modelTemp2.grRules,'APT')));

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