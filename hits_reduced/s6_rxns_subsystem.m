%% Asignar rxns y subsystems a genes a partir de archivos de hits BLAST
% Busca archivos: Hits_THM10_vs_*_protein_reduced.csv
% Para cada uno, carga el modelo .mat desde la carpeta anterior y asigna reacciones y subsistemas.
% Requiere que el modelo tenga rxnGeneMat (se construirá si es posible, pero no usa grRules directamente).

initCobraToolbox()

% 1. Obtener todos los archivos CSV con el patrón específico
archivosCSV = dir('Hits_THM10_vs_*_protein_reduced.csv');
if isempty(archivosCSV)
    error('No se encontraron archivos con el patrón "Hits_THM10_vs_*_protein_reduced.csv"');
end

fprintf('Se encontraron %d archivos para procesar:\n', length(archivosCSV));
for i = 1:length(archivosCSV)
    fprintf('  %d. %s\n', i, archivosCSV(i).name);
end

% Procesar cada archivo
for f = 1:length(archivosCSV)
    archivoCSV = archivosCSV(f).name;
    fprintf('\n==================================================\n');
    fprintf('Procesando archivo: %s\n', archivoCSV);
    
    %% 2. Extraer nombre del modelo desde el nombre del archivo
    tokens = regexp(archivoCSV, 'vs_([^_]+(?:_[^_]+)*?)_protein_reduced\.csv$', 'tokens');
    if isempty(tokens)
        fprintf('  ⚠️ No se pudo extraer el nombre del modelo de %s. Saltando.\n', archivoCSV);
        continue;
    end
    nombreModelo = tokens{1}{1};
    fprintf('  Nombre del modelo detectado: %s\n', nombreModelo);
    
    %% 3. Cargar el modelo .mat desde la carpeta anterior
    archivoMat = fullfile('..', [nombreModelo '.mat']);
    if ~exist(archivoMat, 'file')
        fprintf('  ⚠️ No se encuentra el modelo %s. Saltando.\n', archivoMat);
        continue;
    end
    
    fprintf('  Cargando modelo %s ...\n', archivoMat);
    load(archivoMat);
    
    % Identificar la variable que contiene el modelo COBRA
    vars = whos('-file', archivoMat);
    model = [];
    for i = 1:length(vars)
        tmp = load(archivoMat, vars(i).name);
        if isstruct(tmp.(vars(i).name)) && (isfield(tmp.(vars(i).name), 'genes') || isfield(tmp.(vars(i).name), 'rxns'))
            model = tmp.(vars(i).name);
            fprintf('  Modelo cargado desde variable: %s\n', vars(i).name);
            break;
        end
    end
    
    if isempty(model)
        fprintf('  ⚠️ No se encontró modelo COBRA en %s. Saltando.\n', archivoMat);
        continue;
    end
    
    if ~isfield(model, 'genes') || ~isfield(model, 'rxns')
        fprintf('  ⚠️ El modelo no tiene genes o rxns. Saltando.\n');
        continue;
    end
    
    %% 4. Construir mapa de gen -> reacciones y subsistemas (solo con rxnGeneMat)
    % Verificar o construir rxnGeneMat
    if ~isfield(model, 'rxnGeneMat')
        if isfield(model, 'grRules') || isfield(model, 'rules')
            fprintf('  Construyendo rxnGeneMat desde reglas...\n');
            try
                model = buildRxnGeneMat(model);
                fprintf('  rxnGeneMat construida exitosamente.\n');
            catch ME
                fprintf('  ❌ Error: No se pudo construir rxnGeneMat: %s\n', ME.message);
                fprintf('  El modelo no será procesado (no se permite alternativa con grRules).\n');
                continue;
            end
        else
            fprintf('  ⚠️ El modelo no tiene rxnGeneMat ni reglas para construirla. Saltando.\n');
            continue;
        end
    end
    
    % Determinar campo de subsistemas y convertir a string array
    if isfield(model, 'subSystems')
        subsysRaw = model.subSystems;
    elseif isfield(model, 'rxnSubSystems')
        subsysRaw = model.rxnSubSystems;
    else
        fprintf('  ⚠️ El modelo no tiene campo de subsistemas. Se dejará vacío.\n');
        subsysRaw = cell(length(model.rxns), 1);
    end
    
    % Convertir subsistemas a cell array de strings (uniforme)
    if iscell(subsysRaw)
        subsysList = cellfun(@(x) string(x), subsysRaw, 'UniformOutput', false);
    else
        subsysList = arrayfun(@(x) string(x), subsysRaw, 'UniformOutput', false);
    end
    
    % Inicializar mapas
    geneRxns = containers.Map();   % gen -> lista de reacciones (cell array de strings)
    geneSubsys = containers.Map(); % gen -> lista de subsistemas (string array)
    
    nRxns = length(model.rxns);
    fprintf('  Usando rxnGeneMat para asignar genes a reacciones.\n');
    
    for i = 1:nRxns
        % Subsistema de la reacción i
        subsys = '';
        if i <= length(subsysList) && ~isempty(subsysList{i})
            subsys = subsysList{i};
            % Asegurar que subsys sea string
            if ~isstring(subsys)
                subsys = string(subsys);
            end
        end
        
        % Genes involucrados
        geneIndices = find(model.rxnGeneMat(i, :));
        rxnName = model.rxns{i};
        
        for j = 1:length(geneIndices)
            gen = model.genes{geneIndices(j)};
            
            % --- Añadir reacción ---
            if geneRxns.isKey(gen)
                rxnList = geneRxns(gen);
                if ~any(strcmp(rxnList, rxnName))
                    rxnList{end+1} = rxnName;
                    geneRxns(gen) = rxnList;
                end
            else
                geneRxns(gen) = {rxnName};
            end
            
            % --- Añadir subsistema (usando string arrays para evitar errores de dimensión) ---
            if ~isempty(subsys)
                if geneSubsys.isKey(gen)
                    curr = geneSubsys(gen);
                    % Asegurar que curr sea string array
                    if ~isstring(curr)
                        curr = string(curr);
                    end
                    % Comparar con ismember (robusto)
                    if ~any(ismember(curr, subsys))
                        geneSubsys(gen) = [curr, subsys];
                    end
                else
                    geneSubsys(gen) = subsys;
                end
            else
                if ~geneSubsys.isKey(gen)
                    geneSubsys(gen) = string.empty;
                end
            end
        end
    end
    
    fprintf('  Mapa construido: %d genes con reacciones.\n', geneRxns.Count);
    
    %% 5. Leer el archivo CSV de hits
    fid = fopen(archivoCSV, 'r');
    primeraLinea = fgetl(fid);
    fclose(fid);
    if contains(primeraLinea, '	')
        delimitador = '\t';
    else
        delimitador = ',';
    end
    
    try
        tablaHits = readtable(archivoCSV, 'FileType', 'text', 'Delimiter', delimitador, 'PreserveVariableNames', true);
    catch ME
        fprintf('  ⚠️ Error al leer el CSV: %s. Saltando.\n', ME.message);
        continue;
    end
    
    % Identificar columna 'model ID'
    varNames = tablaHits.Properties.VariableNames;
    colModelID = '';
    for i = 1:length(varNames)
        if contains(lower(varNames{i}), 'model') && contains(lower(varNames{i}), 'id')
            colModelID = varNames{i};
            break;
        end
    end
    if isempty(colModelID) && width(tablaHits) >= 2
        colModelID = varNames{2};
        fprintf('  Usando segunda columna como "model ID": %s\n', colModelID);
    elseif isempty(colModelID)
        fprintf('  ⚠️ No se encontró columna de IDs del modelo. Saltando.\n');
        continue;
    else
        fprintf('  Columna de IDs del modelo: %s\n', colModelID);
    end
    
    % Extraer los IDs (manejo numérico y de cadenas)
    rawIDs = tablaHits.(colModelID);
    if isnumeric(rawIDs)
        modelIDs = strtrim(cellstr(num2str(rawIDs(:))));
    elseif iscell(rawIDs)
        modelIDs = strtrim(cellstr(rawIDs));
    else
        modelIDs = strtrim(cellstr(rawIDs));
    end
    nHits = length(modelIDs);
    fprintf('  Total de hits (genes) en el archivo: %d\n', nHits);
    
    %% 6. Asignar reacciones y subsistemas a cada gen del CSV
    resultados = cell(nHits, 3);
    for i = 1:nHits
        gen = modelIDs{i};
        if geneRxns.isKey(gen)
            reacciones = strjoin(geneRxns(gen), '; ');
            if geneSubsys.isKey(gen)
                subsys_array = geneSubsys(gen);
                % Asegurar que sea string array para strjoin
                if ~isstring(subsys_array)
                    subsys_array = string(subsys_array);
                end
                if isempty(subsys_array)
                    subsistemas = 'Sin subsistema';
                else
                    subsistemas = strjoin(subsys_array, '; ');
                end
            else
                subsistemas = 'Sin subsistema';
            end
            status = 'Encontrado';
        else
            reacciones = 'No encontrado en modelo';
            subsistemas = '';
            status = 'No encontrado';
        end
        resultados{i,1} = gen;
        resultados{i,2} = reacciones;
        resultados{i,3} = subsistemas;
    end
    
    tablaResultados = table(resultados(:,1), resultados(:,2), resultados(:,3), ...
        'VariableNames', {'Gen', 'Reacciones', 'Subsistemas'});
    
    %% 7. Guardar resultados
    [~, nombreBase, ~] = fileparts(archivoCSV);
    archivoSalida = fullfile(pwd, ['rxns_subsystems_' nombreBase '.csv']);
    writetable(tablaResultados, archivoSalida);
    fprintf('  Resultados guardados en: %s\n', archivoSalida);
    
    % Mostrar resumen
    fprintf('\n  --- RESUMEN PARA %s ---\n', nombreModelo);
    numEncontrados = sum(~strcmp(tablaResultados.Reacciones, 'No encontrado en modelo'));
    fprintf('  Genes encontrados en modelo: %d\n', numEncontrados);
    fprintf('  Genes NO encontrados: %d\n', nHits - numEncontrados);
    if nHits > 0
        fprintf('  Porcentaje encontrado: %.1f%%\n', 100 * numEncontrados / nHits);
    end
    
    fprintf('\n  Ejemplo de los primeros 5 resultados:\n');
    disp(tablaResultados(1:min(5, nHits), :));
end

fprintf('\n==================================================\n');
fprintf('Procesamiento completado para %d archivos.\n', length(archivosCSV));