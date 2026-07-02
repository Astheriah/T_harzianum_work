initCobraToolbox()

% Para obtener todos los genes de los archivos .mat y exportar sus genes en .csv resultando en NOMBREMODELO_genes.csv
% Encontrar todos los archivos .mat en el directorio actual
matFiles = dir('*.mat');

excluir = 'blast_THM10.mat';  

matFiles = matFiles(~strcmp({matFiles.name}, excluir));

% Procesar cada archivo .mat
for fileIdx = 1:length(matFiles)
    try
        matFile = matFiles(fileIdx);
        filename = matFile.name;
        fprintf('\n=== Procesando archivo %d/%d: %s ===\n', fileIdx, length(matFiles), filename);
        
        % Cargar el archivo .mat
        data = load(filename);
        
        % Identificar la variable que contiene el modelo COBRA
        varNames = fieldnames(data);
        model = [];
        for i = 1:length(varNames)
            var = data.(varNames{i});
            if isstruct(var) && (isfield(var, 'genes') || isfield(var, 'rxns') || isfield(var, 'S'))
                model = var;
                modelVarName = varNames{i};
                break;
            end
        end
        
        if isempty(model)
            error('No se encontró una estructura de modelo COBRA en el archivo %s', filename);
        end
        
        % Crear tabla con genes
        if isfield(model, 'geneNames') && ~isempty(model.geneNames)
            geneTable = table(model.genes, model.geneNames, 'VariableNames', {'GeneID', 'GeneName'});
        elseif isfield(model, 'genes')
            geneTable = table(model.genes, 'VariableNames', {'GeneID'});
        else
            error('El modelo no tiene campo genes');
        end
        
        fprintf('Total genes: %d\n', height(geneTable));
        
        % Usar el ID del modelo si está disponible, sino usar nombre del archivo (sin .mat)
        if isfield(model, 'id') && ~isempty(model.id)
            modelName = model.id;
        else
            [~, modelName, ~] = fileparts(filename);
        end
        
        % Limpiar nombre para archivo CSV
        cleanModelName = regexprep(modelName, '[\\/*?:"<>|]', '');
        outputFilename = sprintf('%s_genes.csv', cleanModelName);
        
        % Guardar en la misma carpeta (sin subcarpeta)
        writetable(geneTable, outputFilename);
        
        % Mostrar primeras filas
        if height(geneTable) > 0
            fprintf('Primeros 3 genes:\n');
            disp(geneTable(1:min(3, height(geneTable)), :));
        else
            fprintf('El modelo no contiene genes\n');
        end
        
    catch ME
        fprintf('Error procesando %s: %s\n', filename, ME.message);
    end
end