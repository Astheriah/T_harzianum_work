%% Para cargar los 2 archivos .mat de s3_blast.m y exportar cada tabla de hits a un archivo CSV
%% Los resultados de _reduced dentro de una carpeta denominada hits_reduced
%% RESULTADO: Archivos CSV con nombre Hits_TH10_vs_MODELO_protein_1.csv

load('blast_THM10.mat');


% Exportar cada tabla de hits a un archivo CSV
for i = 1:size(HitsTables, 1)
    for j = 1:size(HitsTables, 2)
        if ~isempty(HitsTables{i,j}) && istable(HitsTables{i,j})
            % Formato de nombre para el archivo
            refName = strrep(list2{i}, '.fasta', '');
            cleanName = strrep(cleanFiles{j}, '.faa', '');
            cleanName = strrep(cleanName, '_clean', '');
            fileName = sprintf('Hits_%s_vs_%s.csv', refName, cleanName);
            
            % Guardar directamente usando solo el nombre del archivo
            writetable(HitsTables{i,j}, fileName);
            
            fprintf('Exportado: %s (%d hits)\n', fileName, height(HitsTables{i,j}));
        end
    end
end

load('blast_reduced_THM10.mat');

outputDir = 'hits_reduced';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Carpeta creada: %s\n', outputDir);
end

% Exportar en directorio hits_reduced
for i = 1:size(HitsTables, 1)
    for j = 1:size(HitsTables, 2)
        if ~isempty(HitsTables{i,j}) && istable(HitsTables{i,j})
            % Formato de nombre nombre para el archivo
            refName = strrep(list2{i}, '.fasta', '');
            cleanName = strrep(cleanFiles{j}, '.faa', '');
            cleanName = strrep(cleanName, '_clean', '');
            fileName = sprintf('Hits_%s_vs_%s.csv', refName, cleanName);
            
            % Ruta completa incluyendo la carpeta
            filePath = fullfile(outputDir, fileName);
            
            writetable(HitsTables{i,j}, filePath);
            
            fprintf('Exportado: %s (%d hits)\n', fileName, height(HitsTables{i,j}));
        end
    end
end