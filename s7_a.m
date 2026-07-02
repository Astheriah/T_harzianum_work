% Script: procesar_hits_con_modelo_COBRA.m
% Calcula cobertura de genes y reacciones, indicando cuáles reacciones
% tienen genes no presentes en el CSV.

clear; clc;
initCobraToolbox()

% ================= CONFIGURACIÓN =================
carpetaCSV = 'hits_reduced';
carpetaMAT = '.';
patronCSV  = 'Hits_*_vs_*_protein_reduced.csv';

% ================= FUNCIONES AUXILIARES =================
function out = convertirACellStr(datos)
    if isempty(datos)
        out = {};
        return;
    end
    if ischar(datos) || isstring(datos)
        out = cellstr(datos);
    elseif iscategorical(datos)
        out = cellstr(datos);
    elseif isnumeric(datos) || islogical(datos)
        out = cellstr(num2str(datos(:)));
    elseif istable(datos)
        nombres = datos.Properties.VariableNames;
        idx = find(contains(nombres, 'gene', 'IgnoreCase', true), 1);
        if isempty(idx)
            idx = 1;
        end
        out = table2cell(datos(:, idx));
    elseif iscell(datos)
        out = datos(:);
    else
        error('Tipo de dato no soportado: %s', class(datos));
    end
    out = strtrim(out);
    out = out(:);
end

function [tabla, delimUsado] = leerCSVconDelimitador(ruta)
    delimiters = {'\t', ',', ';', ' ', 'auto'};
    mejorAncho = 0;
    tablaMejor = [];
    delimUsado = '';
    for d = delimiters
        try
            if strcmp(d{1}, ' ')
                opts = detectImportOptions(ruta, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
            elseif strcmp(d{1}, 'auto')
                opts = detectImportOptions(ruta);
            else
                opts = detectImportOptions(ruta, 'Delimiter', d{1});
            end
            opts.VariableNamingRule = 'preserve';
            opts.VariableTypes = repmat({'string'}, 1, length(opts.VariableNames));
            tbl = readtable(ruta, opts);
            ancho = size(tbl, 2);
            if ancho > mejorAncho
                mejorAncho = ancho;
                tablaMejor = tbl;
                delimUsado = d{1};
                if ancho >= 2
                    break;
                end
            end
        catch
            continue;
        end
    end
    if isempty(tablaMejor)
        error('No se pudo leer el archivo con ningún delimitador probado.');
    end
    tabla = tablaMejor;
end

% ================= PROCESAMIENTO PRINCIPAL =================
archivos = dir(fullfile(carpetaCSV, patronCSV));
if isempty(archivos)
    error('No se encontraron archivos CSV con el patrón "%s" en "%s".', patronCSV, carpetaCSV);
end

numArchivos = length(archivos);
% Columnas: Modelo, IDsCSV, Pres, %Pres, Genes, Faltan, %Falt, TotalRxns, RxnsGPR, SinGPR, NoPresentes, %NoPres, Completas, %Completas
resultados = cell(numArchivos, 14);

for i = 1:numArchivos
    nombreArchivo = archivos(i).name;
    rutaCompleta = fullfile(carpetaCSV, nombreArchivo);
    fprintf('\nProcesando: %s\n', nombreArchivo);
    
    try
        [datos, delim] = leerCSVconDelimitador(rutaCompleta);
        fprintf('  Delimitador detectado: "%s"\n', delim);
    catch ME
        warning('No se pudo leer el archivo %s: %s', nombreArchivo, ME.message);
        continue;
    end
    if size(datos, 2) < 2
        warning('El archivo %s no tiene al menos 2 columnas.', nombreArchivo);
        continue;
    end
    modelIDs = convertirACellStr(datos{:, 2});
    totalIDs = length(modelIDs);
    
    tokens = regexp(nombreArchivo, 'vs_(.*?)_protein', 'tokens', 'once');
    if isempty(tokens)
        warning('No se pudo extraer el nombre del modelo de %s', nombreArchivo);
        continue;
    end
    nombreModelo = tokens{1};
    fprintf('  Nombre del modelo: %s\n', nombreModelo);
    
    archivoMat = fullfile(carpetaMAT, [nombreModelo '.mat']);
    if ~exist(archivoMat, 'file')
        fprintf('  ⚠️ No se encuentra el modelo %s. Saltando.\n', archivoMat);
        continue;
    end
    fprintf('  Cargando modelo %s ...\n', archivoMat);
    load(archivoMat);
    
    vars = whos('-file', archivoMat);
    model = [];
    for j = 1:length(vars)
        tmp = load(archivoMat, vars(j).name);
        if isstruct(tmp.(vars(j).name)) && (isfield(tmp.(vars(j).name), 'genes') || isfield(tmp.(vars(j).name), 'rxns'))
            model = tmp.(vars(j).name);
            fprintf('  Modelo cargado desde variable: %s\n', vars(j).name);
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
    
    % Construir rxnGeneMat si no existe
    if ~isfield(model, 'rxnGeneMat')
        if isfield(model, 'grRules') || isfield(model, 'rules')
            fprintf('  Construyendo rxnGeneMat desde reglas...\n');
            try
                model = buildRxnGeneMat(model);
                fprintf('  rxnGeneMat construida exitosamente.\n');
            catch ME
                fprintf('  ❌ Error: No se pudo construir rxnGeneMat: %s\n', ME.message);
                fprintf('  El modelo no será procesado.\n');
                continue;
            end
        else
            fprintf('  ⚠️ El modelo no tiene rxnGeneMat ni reglas. Saltando.\n');
            continue;
        end
    end
    
    genes = convertirACellStr(model.genes);
    rxns = convertirACellStr(model.rxns);
    rxnGeneMat = model.rxnGeneMat;
    
    % Verificar dimensiones
    if size(rxnGeneMat, 1) ~= length(rxns)
        if size(rxnGeneMat, 2) == length(rxns)
            rxnGeneMat = rxnGeneMat';
            fprintf('  Se transpuso rxnGeneMat.\n');
        else
            fprintf('  ⚠️ Dimensiones inconsistentes. Saltando análisis de reacciones.\n');
            rxnGeneMat = [];
        end
    end
    
    % Presencia de IDs en genes
    presentes = ismember(modelIDs, genes);
    numPresentes = sum(presentes);
    proporcion = numPresentes / totalIDs * 100;
    
    % Genes faltantes en el CSV
    encontradosEnCSV = ismember(genes, modelIDs);
    numGenesNoEnCSV = length(genes) - sum(encontradosEnCSV);
    porcentajeGenesNoEnCSV = numGenesNoEnCSV / length(genes) * 100;
    
    % Análisis de reacciones
    if ~isempty(rxnGeneMat)
        genesFaltantesIdx = ~encontradosEnCSV;
        % Reacciones que tienen al menos un gen no presente en CSV
        reaccionesNoPresentes = any(rxnGeneMat(:, genesFaltantesIdx), 2);
        numReaccionesNoPresentes = nnz(reaccionesNoPresentes);
        totalRxnsConGen = nnz(any(rxnGeneMat, 2));
        totalRxns = length(rxns);
        rxnsSinGPR = totalRxns - totalRxnsConGen;
        reaccionesCompletas = totalRxnsConGen - numReaccionesNoPresentes;
        if totalRxnsConGen > 0
            porcentajeNoPresentes = numReaccionesNoPresentes / totalRxnsConGen * 100;
            porcentajeCompletas = reaccionesCompletas / totalRxnsConGen * 100;
        else
            porcentajeNoPresentes = NaN;
            porcentajeCompletas = NaN;
        end
    else
        numReaccionesNoPresentes = NaN;
        totalRxnsConGen = NaN;
        totalRxns = NaN;
        rxnsSinGPR = NaN;
        reaccionesCompletas = NaN;
        porcentajeNoPresentes = NaN;
        porcentajeCompletas = NaN;
    end
    
    % Almacenar resultados (14 columnas)
    resultados{i,1} = nombreModelo;
    resultados{i,2} = totalIDs;
    resultados{i,3} = numPresentes;
    resultados{i,4} = proporcion;
    resultados{i,5} = length(genes);
    resultados{i,6} = numGenesNoEnCSV;
    resultados{i,7} = porcentajeGenesNoEnCSV;
    resultados{i,8} = totalRxns;
    resultados{i,9} = totalRxnsConGen;
    resultados{i,10} = rxnsSinGPR;
    resultados{i,11} = numReaccionesNoPresentes;
    resultados{i,12} = porcentajeNoPresentes;
    resultados{i,13} = reaccionesCompletas;
    resultados{i,14} = porcentajeCompletas;
    
    % Mostrar en pantalla
    fprintf('  Total de IDs en CSV: %d\n', totalIDs);
    fprintf('  IDs presentes en genes: %d (%.2f%%)\n', numPresentes, proporcion);
    fprintf('  Total de genes en modelo: %d\n', length(genes));
    fprintf('  Genes del modelo NO presentes en CSV: %d (%.2f%%)\n', numGenesNoEnCSV, porcentajeGenesNoEnCSV);
    if ~isempty(rxnGeneMat)
        fprintf('  Total de reacciones en modelo: %d\n', totalRxns);
        fprintf('  Reacciones con GPR: %d\n', totalRxnsConGen);
        fprintf('  Reacciones sin GPR: %d\n', rxnsSinGPR);
        fprintf('  Reacciones con al menos un gen NO presente en CSV: %d (%.2f%% del total con GPR)\n', ...
            numReaccionesNoPresentes, porcentajeNoPresentes);
        fprintf('  Reacciones con todos sus genes presentes en CSV: %d (%.2f%%)\n', ...
            reaccionesCompletas, porcentajeCompletas);
    else
        fprintf('  (No se pudieron analizar reacciones)\n');
    end
end

% ================= RESUMEN FINAL =================
fprintf('\n=== RESUMEN ===\n');
fprintf('%-15s %-10s %-10s %-10s %-10s %-12s %-12s %-12s %-12s %-12s %-15s %-12s %-15s %-12s\n', ...
    'Modelo', 'IDs CSV', 'Pres.', '%Pres.', 'Genes', 'Faltan', '%Falt.', ...
    'TotalRxns', 'RxnsGPR', 'SinGPR', 'NoPresentes', '%NoPres', 'Completas', '%Comp.');
for i = 1:numArchivos
    if ~isempty(resultados{i,1})
        if isnan(resultados{i,8})
            totRx = 'N/A'; gpr = 'N/A'; sinG = 'N/A'; noPres = 'N/A'; pctNo = 'N/A'; comp = 'N/A'; pctComp = 'N/A';
        else
            totRx = num2str(resultados{i,8});
            gpr = num2str(resultados{i,9});
            sinG = num2str(resultados{i,10});
            noPres = num2str(resultados{i,11});
            pctNo = sprintf('%.2f', resultados{i,12});
            comp = num2str(resultados{i,13});
            pctComp = sprintf('%.2f', resultados{i,14});
        end
        fprintf('%-15s %-10d %-10d %-9.2f %-10d %-12d %-11.2f %-12s %-12s %-12s %-15s %-12s %-15s %-12s\n', ...
            resultados{i,1}, resultados{i,2}, resultados{i,3}, resultados{i,4}, ...
            resultados{i,5}, resultados{i,6}, resultados{i,7}, ...
            totRx, gpr, sinG, noPres, pctNo, comp, pctComp);
    end
end


% ================= EXPORTAR RESULTADOS A CSV =================

% Crear tabla vacía con nombres de columna
nombresColumnas = {'Modelo', 'IDsCSV', 'Presentes', 'PorcentajePres', ...
                   'GenesTotales', 'GenesFaltantes', 'PorcentajeFaltantes', ...
                   'TotalRxns', 'RxnsConGPR', 'RxnsSinGPR', ...
                   'RxnsNoPresentes', 'PorcentajeNoPresentes', ...
                   'RxnsCompletas', 'PorcentajeCompletas'};

% Inicializar celdas para la tabla
tablaExport = cell(numArchivos, length(nombresColumnas));

for i = 1:numArchivos
    if ~isempty(resultados{i,1})
        % Copiar datos numéricos y convertir NaN a cadena vacía o 'N/A'
        for j = 1:length(nombresColumnas)
            if j == 1
                % Modelo (string)
                tablaExport{i, j} = resultados{i, j};
            else
                val = resultados{i, j};
                if isnan(val)
                    tablaExport{i, j} = 'N/A';
                elseif ischar(val) || isstring(val)
                    tablaExport{i, j} = char(val);
                else
                    % Números: formatear con dos decimales si es porcentaje
                    if any(j == [4, 7, 12, 14]) % columnas de porcentaje
                        tablaExport{i, j} = sprintf('%.2f', val);
                    else
                        tablaExport{i, j} = num2str(val);
                    end
                end
            end
        end
    end
end

% Convertir a tabla y escribir CSV
T = cell2table(tablaExport, 'VariableNames', nombresColumnas);
writetable(T, 'resultados_1.csv');
