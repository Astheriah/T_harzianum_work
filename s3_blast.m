%% Para poder hacer BLAST utilizando el .fasta (proteinas) del genoma de interes contra todos los archivos .faa de interes (reducidos a solo los que se encuentran en el modelo)

%% RESULTADO: blast_THM10.mat   (Blast respecto a protein_clean_1.faa, TODOS los genes del .faa)
%% RESULTADO: blast_reduced_THM10.mat   (Blast respecto a protein_clan_reduced.faa, solo los genes PRESENTES en el .mat)

folderPath = pwd;

fileList = dir(folderPath);
fileNames = {fileList.name};
cleanFiles = fileNames(contains(fileNames, 'protein_clean_reduced'))';

% Verificar que se encontraron archivos
if isempty(cleanFiles)
    error('No se encontraron archivos con "protein_clean_reduced" en la carpeta actual.');
end

list2 = {'THM10.fasta'};


for i = 1 : length(list2)
    parfor j = 1 : length(cleanFiles)
        blastStruct = getBlast(char(cleanFiles{j,1}), char(cleanFiles{j,1}), list2{i,1}, list2{i,1});
        innerBlast{i,j} = blastStruct;
    end
end

% --- PROCESAMIENTO---
innerBlast2 = innerBlast;
innerBlast3 = innerBlast2;
HitsInfo = {}; HitsTables = {};
for k = 1:size(innerBlast2,2)
    k
    for z = 1 : size(innerBlast2,1)
        BLAST = innerBlast2{z,k};
        parameters = [1e-5 20 30];
        cont2 = 1;
        SameGenes2={};
        for j=1:length(BLAST(1).evalue)
               if BLAST(1).evalue(j,1) <= parameters(1,1) && ...
                  BLAST(1).aligLen(j,1) >= parameters(1,2) && ...
                  BLAST(1).identity(j,1) >= parameters(1,3)
                    SameGenes2{cont2,1} = BLAST(1).fromGenes(j,1);
                    SameGenes2{cont2,2} = BLAST(1).toGenes(j,1);
                    SameGenes2{cont2,3} = BLAST(1).evalue(j,1);
                    SameGenes2{cont2,4} = BLAST(1).identity(j,1);
                    SameGenes2{cont2,5} = BLAST(1).aligLen(j,1);
                    SameGenes2{cont2,6} = BLAST(1).bitscore(j,1);
                    SameGenes2{cont2,7} = BLAST(1).ppos(j,1);
                    cont2 = cont2 + 1;
               end
        end
        
        listChar = cellfun(@char, SameGenes2(:,1), 'UniformOutput', false);
        [list, A, C] = unique(listChar);
        bestHit = {};
        for i = 1:length(list)
            pos = find(strcmp(listChar, list{i,1}));
            identVal = cell2mat(SameGenes2(pos,4));
            proof = SameGenes2(pos,1);
            [getMax, idxMax] = max(identVal);
            bestHit(i,:) = SameGenes2(pos(idxMax),:);
        end
    
        innerBlast3{z,k}(1).fromGenes = cellfun(@char, bestHit(:,1), 'UniformOutput', false);
        innerBlast3{z,k}(1).toGenes   = cellfun(@char, bestHit(:,2), 'UniformOutput', false);
        innerBlast3{z,k}(1).evalue    = cell2mat(bestHit(:,3));
        innerBlast3{z,k}(1).identity  = cell2mat(bestHit(:,4));
        innerBlast3{z,k}(1).aligLen   = cell2mat(bestHit(:,5));
        innerBlast3{z,k}(1).bitscore  = cell2mat(bestHit(:,6));
        innerBlast3{z,k}(1).ppos      = cell2mat(bestHit(:,7));
        HitsInfo{z,k} = bestHit;
        HitsTables{z,k} = cell2table(bestHit, 'VariableNames', {'template ID','model ID','evalue','identity','length','bitscore','ppos'});
    end
end

HitsInfo2 = HitsInfo';
lenHits = cellfun(@length, HitsInfo2);
[M, I] = max(lenHits, [], 1);
[sortVar, idxSort] = sort(lenHits, 'descend');

save('blast_reduced_THM10.mat');

disp('Proceso completado');