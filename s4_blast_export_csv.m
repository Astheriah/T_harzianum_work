%% Para cargar los 2 archivos .mat de s3_blast.m y exportar cada tabla de hits a un archivo CSV
%% Los resultados de _reduced dentro de una carpeta denominada hits_reduced
%% RESULTADO: Archivos CSV con nombre Hits_TH10_vs_MODELO_protein_1.csv

% Intenta cargar y exportar el primer archivo
try
    load('blast_THM10.mat');
    fprintf('Cargado blast_THM10.mat correctamente.\n');
    
    % Exportar cada tabla de hits a CSV (en el directorio actual)
    for i = 1:size(HitsTables, 1)
        for j = 1:size(HitsTables, 2)
            if ~isempty(HitsTables{i,j}) && istable(HitsTables{i,j})
                refName = strrep(list2{i}, '.fasta', '');
                cleanName = strrep(cleanFiles{j}, '.faa', '');
                cleanName = strrep(cleanName, '_clean', '');
                fileName = sprintf('Hits_%s_vs_%s.csv', refName, cleanName);
                writetable(HitsTables{i,j}, fileName);
                fprintf('Exportado: %s\n', fileName);
            end
        end
    end
    
catch ME
    fprintf('No se pudo cargar blast_THM10.mat: %s\n', ME.message);
    fprintf('Saltando exportación de este archivo.\n');
end

% Intenta cargar y exportar el segundo archivo
try
    load('blast_reduced_THM10.mat');
    fprintf('Cargado blast_reduced_THM10.mat correctamente.\n');
    
    outputDir = 'hits_reduced';
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
        fprintf('Carpeta creada: %s\n', outputDir);
    end
    
    % Exportar cada tabla de hits a CSV dentro de hits_reduced
    for i = 1:size(HitsTables, 1)
        for j = 1:size(HitsTables, 2)
            if ~isempty(HitsTables{i,j}) && istable(HitsTables{i,j})
                refName = strrep(list2{i}, '.fasta', '');
                cleanName = strrep(cleanFiles{j}, '.faa', '');
                cleanName = strrep(cleanName, '_clean', '');
                fileName = sprintf('Hits_%s_vs_%s.csv', refName, cleanName);
                filePath = fullfile(outputDir, fileName);
                writetable(HitsTables{i,j}, filePath);
                fprintf('Exportado: %s\n', fileName);
            end
        end
    end
    
catch ME
    fprintf('No se pudo cargar blast_reduced_THM10.mat: %s\n', ME.message);
    fprintf('Saltando exportación de este archivo.\n');
end