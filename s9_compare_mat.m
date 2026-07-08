% Comparar NOMBRES de campos de modelos metabólicos y porcentaje de datos
baseFileName = 'iML1515.mat';
baseFilePath = fullfile('compare_model', baseFileName);

% Lista de excepciones: archivos .mat que NO quieres analizar como "otros"
exceptions = {'blast_reduced_THM10.mat'};  

% --- Detectar TODOS los .mat en la carpeta actual ---
allMat = dir('*.mat');
allNames = {allMat.name};

% Filtrar: eliminar el baseFile (solo el nombre, sin ruta) y las excepciones
filesToRemove = [exceptions, {baseFileName}];
idx = ~ismember(allNames, filesToRemove);
otherFiles = allNames(idx);

fprintf('Archivos .mat encontrados: %d\n', length(allMat));
fprintf('Archivo base: %s\n', baseFileName);
fprintf('Excepciones aplicadas: %d\n', length(exceptions));
fprintf('Archivos a comparar (otros): %d\n', length(otherFiles));

% Función para cargar y extraer la estructura del modelo (primera variable)
function model = loadModel(filename)
    data = load(filename);
    fnames = fieldnames(data);
    if isempty(fnames)
        error('El archivo %s no contiene variables.', filename);
    end
    model = data.(fnames{1});  % tomar la primera variable (se asume que es el modelo)
    if ~isstruct(model)
        warning('La variable extraída de %s no es una estructura.', filename);
    end
end

% Función para calcular el porcentaje de elementos "con datos"
function pct = getFieldPercent(val)
    if isempty(val)
        pct = 0;
        return;
    end
    if isnumeric(val)
        if issparse(val)
            nonZeroCount = nnz(val);
            total = numel(val);
            pct = nonZeroCount / total * 100;
        else
            nonZero = (val ~= 0) & ~isnan(val);
            pct = sum(nonZero(:)) / numel(val) * 100;
        end
    elseif iscell(val)
        nonEmpty = cellfun(@(x) ~isempty(x) && ~(ischar(x) && strcmp(x, '')), val);
        pct = sum(nonEmpty(:)) / numel(val) * 100;
    elseif isstruct(val)
        pct = NaN;
    else
        pct = 100;
    end
    pct = double(pct);
    if ~isscalar(pct)
        pct = pct(1);
    end
end

% Cargar modelo base (desde la carpeta compare_model)
baseModel = loadModel(baseFilePath);
baseFields = fieldnames(baseModel);
fprintf('Archivo base: %s\n', baseFileName);
fprintf('Campos del modelo (%d):\n', length(baseFields));
for i = 1:length(baseFields)
    field = baseFields{i};
    pct = getFieldPercent(baseModel.(field));
    if isnan(pct)
        fprintf('  %s: N/A\n', field);
    else
        fprintf('  %s: %.1f%%\n', field, pct);
    end
end
fprintf('\n');

% Inicializar almacenamiento para el resumen final
% Cada fila: {nombreArchivo, listaCamposFaltantes, listaPorcentajes}
missingData = cell(length(otherFiles), 3);

% Recorrer cada archivo a comparar
for i = 1:length(otherFiles)
    otherFile = otherFiles{i};
    fprintf('--- Comparando con: %s ---\n', otherFile);
    
    otherModel = loadModel(otherFile);
    otherFields = fieldnames(otherModel);
    
    commonFields = intersect(baseFields, otherFields);
    onlyInBase = setdiff(baseFields, otherFields);
    onlyInOther = setdiff(otherFields, baseFields);
    
    % ---- Guardar datos de campos faltantes (no presentes en este otro) ----
    if ~isempty(onlyInBase)
        pcts = cell(length(onlyInBase), 1);
        for j = 1:length(onlyInBase)
            field = onlyInBase{j};
            pcts{j} = getFieldPercent(baseModel.(field));
        end
        missingData(i, :) = {otherFile, onlyInBase, pcts};
    else
        missingData(i, :) = {otherFile, {}, {}};
    end
    
    % ---- Campos comunes ----
    fprintf('Campos comunes (%d):\n', length(commonFields));
    if isempty(commonFields)
        fprintf('  (ninguno)\n');
    else
        for j = 1:length(commonFields)
            field = commonFields{j};
            pct = getFieldPercent(baseModel.(field));
            if isnan(pct)
                fprintf('  %s: N/A\n', field);
            else
                fprintf('  %s: %.1f%%\n', field, pct);
            end
        end
    end
    
    % ---- Solo en base (NO presentes en el otro) ----
    fprintf('Campos NO presentes en %s (%d):\n', otherFile, length(onlyInBase));
    if isempty(onlyInBase)
        fprintf('  (ninguno)\n');
    else
        for j = 1:length(onlyInBase)
            field = onlyInBase{j};
            pct = getFieldPercent(baseModel.(field));
            if isnan(pct)
                fprintf('  %s: N/A\n', field);
            else
                fprintf('  %s: %.1f%%\n', field, pct);
            end
        end
    end
    
    % ---- Solo en el otro archivo ----
    fprintf('Solo en %s (%d):\n', otherFile, length(onlyInOther));
    if isempty(onlyInOther)
        fprintf('  (ninguno)\n');
    else
        for j = 1:length(onlyInOther)
            field = onlyInOther{j};
            pct = getFieldPercent(otherModel.(field));
            if isnan(pct)
                fprintf('  %s: N/A\n', field);
            else
                fprintf('  %s: %.1f%%\n', field, pct);
            end
        end
    end
    fprintf('\n');
end

% ---- RESUMEN FINAL: CAMPOS NO PRESENTES EN CADA MODELO ----
fprintf('========== RESUMEN FINAL ==========\n');
fprintf('Campos que están en %s pero NO en cada modelo comparado:\n\n', baseFileName);
for i = 1:length(otherFiles)
    archivo = missingData{i, 1};
    campos = missingData{i, 2};
    pcts = missingData{i, 3};
    if isempty(campos)
        fprintf('Modelo %s: TODOS los campos de %s están presentes.\n', archivo, baseFileName);
    else
        fprintf('Modelo %s (%d campos faltantes):\n', archivo, length(campos));
        for j = 1:length(campos)
            field = campos{j};
            pct = pcts{j};
            if isnan(pct)
                fprintf('  %s: N/A\n', field);
            else
                fprintf('  %s: %.1f%%\n', field, pct);
            end
        end
    end
    fprintf('\n');
end


% Determinar si rxns tienen formato de bigg_id

opts_reac = detectImportOptions('bigg_models_reactions.txt', 'FileType', 'text');
opts_reac.Delimiter = '\t';
T_bigg = readtable('bigg_models_reactions.txt', opts_reac);
T_bigg.Properties.VariableNames = strtrim(T_bigg.Properties.VariableNames)

% --- ANÁLISIS DE REACCIONES (rxns) vs bigg_id ---
% Crear un vector de celdas con todos los bigg_id de la tabla
if istable(T_bigg) && any(strcmp(T_bigg.Properties.VariableNames, 'bigg_id'))
    bigg_ids = T_bigg.bigg_id;
    if iscell(bigg_ids)
        bigg_ids = bigg_ids(:);  % asegurar columna
    else
        % Si es char array, convertir a cell
        bigg_ids = cellstr(bigg_ids);
    end
else
    error('La tabla T_bigg no contiene la columna "bigg_id".');
end

% Función auxiliar para analizar un modelo
function analyzeRxns(model, modelName, bigg_ids)
    if isfield(model, 'rxns')
        rxns = model.rxns;
        if iscell(rxns)
            rxns = rxns(:);
        else
            rxns = cellstr(rxns);
        end
        % Limpiar posibles espacios en blanco
        rxns = strtrim(rxns);
        
        % Comparar con bigg_ids (búsqueda exacta)
        [isFound, idx] = ismember(rxns, bigg_ids);
        numTotal = length(rxns);
        numFound = sum(isFound);
        numMissing = numTotal - numFound;
        pctFound = numFound / numTotal * 100;
        
        fprintf('\n--- ANÁLISIS DE REACCIONES PARA %s ---\n', modelName);
        fprintf('Total de reacciones: %d\n', numTotal);
        fprintf('Reacciones en bigg_id: %d (%.1f%%)\n', numFound, pctFound);
        fprintf('Reacciones NO en bigg_id: %d (%.1f%%)\n', numMissing, 100-pctFound);
        
        if numMissing > 0
            % Listar de algunos que no esten
            missingRxns = rxns(~isFound);
            fprintf('Ejemplos de reacciones no encontradas:\n');
            for k = 1:min(5, length(missingRxns))
                fprintf('  %s\n', missingRxns{k});
            end
        else
            fprintf('Todas las reacciones están en bigg_id.\n');
        end
    else
        fprintf('\n--- ANÁLISIS DE REACCIONES PARA %s ---\n', modelName);
        fprintf('El modelo NO contiene el campo "rxns".\n');
    end
end

% Analizar el modelo base
analyzeRxns(baseModel, baseFileName, bigg_ids);

% Analizar cada otro modelo
for i = 1:length(otherFiles)
    otherModel = loadModel(otherFiles{i});  % reutiliza tu función loadModel
    analyzeRxns(otherModel, otherFiles{i}, bigg_ids);
end