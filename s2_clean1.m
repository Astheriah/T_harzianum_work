
% Para funcionar necesita: iANid1221.mat iANid1221_protein.faa y opcionalmente iANid1221_gtf.gtf (para comparar genes del modelo con genes del GTF)
% donde las terminaciones .mat, _protein.faa _gtf.gtf DEBEN ser estas despues del nombre MISMO nombre del modelo
% Funcionamiento:
% 1) Se exportan los genes del .mat
% 2) Se comparan con los genes del GTF (si existe)
% 3) Se genera un archivo  modelo_protein_clean_reduced.faa con solo las proteínas de los genes en CSV (genes dentro dele modelo)
% 4) Se genera un archivo modelo_protein_clean_1.faa con las proteínas del modelo pero SOLO con las cabeceras de locus_tag

% El funcionamiento genera los archivos si par la generacion del modelo se uso de forma estandar el locus tag como nombre de genes dentro del .mat
% Los que tuvieron error deben corregirse manualmente, ya que no se puede hacer de forma automatica debido a que tienen nombres distintos al locus_tag

% OBJETIVO: generar _protein_clean_1.faa, _protein_clean_reduced.faa con cabeceras de locus_tag (EL MISMO NOMBRE QUE EXISTE DENTRO DE LOS .MAT)

initCobraToolbox()
% Encontrar todos los archivos .mat en el directorio actual
matFiles = dir('*.mat');
excluir = 'blast_THM10.mat';  
matFiles = matFiles(~strcmp({matFiles.name}, excluir));

fprintf('Encontrados %d archivos .mat\n', length(matFiles));

% Estructuras para guardar resúmenes
resumenModelos = {};        % cada fila: {nombreModelo, totalGenes, encontrados, noEncontrados}
modelosSinReducido = {};    % modelos que no generaron el archivo _reduced.faa

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
        
        totalGenes = height(geneTable);
        fprintf('Total genes en modelo: %d\n', totalGenes);
        
        % Usar el ID del modelo si está disponible, sino usar nombre del archivo (sin .mat)
        if isfield(model, 'id') && ~isempty(model.id)
            modelName = model.id;
        else
            [~, modelName, ~] = fileparts(filename);
        end
        
        % Limpiar nombre para archivo CSV
        cleanModelName = regexprep(modelName, '[\\/*?:"<>|]', '');
        outputFilename = sprintf('%s_genes.csv', cleanModelName);
        
        % Guardar CSV con los genes del modelo (en la misma carpeta)
        writetable(geneTable, outputFilename);
        fprintf('✓ Genes exportados a %s\n', outputFilename);
        
        % Mostrar primeras filas del CSV
        if totalGenes > 0
            fprintf('Primeros 3 genes:\n');
            disp(geneTable(1:min(3, totalGenes), :));
        else
            fprintf('El modelo no contiene genes\n');
        end
        
        % --------------------------------------------------------------
        % COMPARACIÓN CON ARCHIVO GTF (mismo nombre base)
        % --------------------------------------------------------------
        [~, baseName, ~] = fileparts(filename);
        gtfFile = [baseName '_gtf.gtf'];
        
        if ~exist(gtfFile, 'file')
            fprintf('⚠️ Archivo GTF no encontrado: %s. No se pudo comparar con GTF.\n', gtfFile);
            resumenModelos(end+1, :) = {modelName, totalGenes, NaN, NaN};
            modelosSinReducido{end+1} = modelName;
        else
            fprintf('\n--- Comparando con GTF: %s ---\n', gtfFile);
            
            % Diccionarios para mapeos
            locusTags = {};
            proteinID_to_locusTag = containers.Map();  % protein_id -> locus_tag
            try
                fid = fopen(gtfFile, 'r');
                if fid == -1
                    error('No se pudo abrir %s', gtfFile);
                end
                while ~feof(fid)
                    line = fgetl(fid);
                    if isempty(line) || line(1) == '#'
                        continue;
                    end
                    fields = strsplit(line, '\t');
                    if length(fields) < 9
                        continue;
                    end
                    feature = fields{3};
                    attributes = fields{9};
                    
                    % Extraer locus_tag y protein_id
                    locus_tag = '';
                    protein_id = '';
                    tokens = regexp(attributes, 'locus_tag "([^"]+)"', 'tokens');
                    if ~isempty(tokens)
                        locus_tag = tokens{1}{1};
                    end
                    tokens = regexp(attributes, 'protein_id "([^"]+)"', 'tokens');
                    if ~isempty(tokens)
                        protein_id = tokens{1}{1};
                    end
                    
                    % Almacenar locus_tag de líneas "gene"
                    if strcmp(feature, 'gene') && ~isempty(locus_tag)
                        locusTags{end+1} = locus_tag;
                    end
                    
                    % Construir mapeo protein_id -> locus_tag (desde líneas que tengan ambos)
                    if ~isempty(protein_id) && ~isempty(locus_tag)
                        proteinID_to_locusTag(protein_id) = locus_tag;
                    end
                end
                fclose(fid);
            catch ME
                fprintf('❌ Error leyendo GTF: %s\n', ME.message);
                continue;
            end
            
            fprintf('Total de genes (locus_tag) en GTF: %d\n', length(locusTags));
            
            % Obtener lista de GeneID del modelo
            geneIDs = geneTable.GeneID;
            if iscell(geneIDs)
                geneIDs = cellstr(geneIDs);
            else
                geneIDs = cellstr(geneIDs);
            end
            geneIDs = strtrim(geneIDs);
            
            locusTagsClean = strtrim(locusTags);
            
            % Comparación: genes del modelo presentes en GTF
            foundInGTF = ismember(geneIDs, locusTagsClean);
            numFound = sum(foundInGTF);
            numNotFound = sum(~foundInGTF);
            
            fprintf('\n--- RESULTADOS DE COMPARACIÓN ---\n');
            fprintf('Genes en CSV: %d\n', length(geneIDs));
            fprintf('Encontrados en GTF (locus_tag): %d\n', numFound);
            fprintf('NO encontrados en GTF: %d\n', numNotFound);
            
            if numNotFound > 0
                notFoundGenes = geneIDs(~foundInGTF);
                fprintf('\nLista de genes NO encontrados en GTF :\n');
                for k = 1:min(2, length(notFoundGenes))
                    fprintf('  - %s\n', notFoundGenes{k});
                end
                if length(notFoundGenes) > 2
                    fprintf('  ... y %d más.\n', length(notFoundGenes)-2);
                end
            end
            
            % --------------------------------------------------------------
            % GENERACIÓN DEL ARCHIVO REDUCIDO (solo proteínas de genes encontrados)
            % --------------------------------------------------------------
            fprintf('\n--- Generando archivo REDUCIDO (solo proteínas de genes en CSV) ---\n');
            faaOriginal = [baseName '_protein.faa'];
            faaReducido = [baseName '_protein_clean_reduced.faa'];
            reducidoGenerado = false;
            
            if ~exist(faaOriginal, 'file')
                fprintf('⚠️ Archivo FASTA original no encontrado: %s. No se puede generar el archivo reducido.\n', faaOriginal);
                modelosSinReducido{end+1} = modelName;
            else
                try
                    [headers, sequences] = readFasta(faaOriginal);
                    totalProteinas = length(headers);
                    
                    % Conjunto de locus_tag que están en los genes encontrados
                    foundLocusTags = geneIDs(foundInGTF);  % son los locus_tag que coincidieron
                    foundSet = containers.Map(foundLocusTags, ones(length(foundLocusTags),1));
                    
                    % Filtrar proteínas: solo aquellas cuyo protein_id mapee a un locus_tag en foundSet
                    kept = 0;
                    eliminated = 0;
                    fidRed = fopen(faaReducido, 'w');
                    if fidRed == -1
                        error('No se pudo crear %s', faaReducido);
                    end
                    
                    for i = 1:length(headers)
                        header = headers{i};
                        protein_id = strtok(header);
                        % Buscar el locus_tag correspondiente a este protein_id
                        if proteinID_to_locusTag.isKey(protein_id)
                            locus = proteinID_to_locusTag(protein_id);
                            if foundSet.isKey(locus)
                                % Escribir esta proteína en el archivo reducido
                                % Usar el locus_tag como cabecera (puede personalizarse)
                                newHeader = locus;
                                fprintf(fidRed, '>%s\n', newHeader);
                                seq = sequences{i};
                                for j = 1:70:length(seq)
                                    fprintf(fidRed, '%s\n', seq(j:min(j+69, end)));
                                end
                                kept = kept + 1;
                            else
                                eliminated = eliminated + 1;
                            end
                        else
                            % No se pudo mapear protein_id -> locus_tag; se elimina
                            eliminated = eliminated + 1;
                        end
                    end
                    fclose(fidRed);
                    reducidoGenerado = true;
                    
                    fprintf('✓ Archivo REDUCIDO guardado como: %s\n', faaReducido);
                    fprintf('Total proteínas procesadas: %d\n', totalProteinas);
                    fprintf('\n--- ESTADÍSTICAS DE FILTRADO PARA %s ---\n', modelName);
                    fprintf('Total de proteínas en FAA original: %d\n', totalProteinas);
                    fprintf('Proteínas MANTENIDAS (en CSV): %d\n', kept);
                    fprintf('Proteínas ELIMINADAS (no en CSV): %d\n', eliminated);
                    fprintf('Porcentaje mantenido: %.2f%%\n', kept/totalProteinas*100);
                    fprintf('Porcentaje eliminado: %.2f%%\n', eliminated/totalProteinas*100);
                    
                catch ME
                    fprintf('❌ Error generando archivo reducido: %s\n', ME.message);
                    modelosSinReducido{end+1} = modelName;
                end
            end
            
            % --------------------------------------------------------------
            % MODIFICACIÓN DEL ARCHIVO _protein.faa SI HAY POCOS NO ENCONTRADOS (1-9)
            % (se mantiene igual que antes)
            % --------------------------------------------------------------
            if numNotFound < 10
                fprintf('\n--- Generando archivo FASTA limpio (solo locus_tag) ---\n');
                faaOriginal = [baseName '_protein.faa'];
                faaNuevo = [baseName '_protein_clean_1.faa'];
                
                if ~exist(faaOriginal, 'file')
                    fprintf('⚠️ Archivo FASTA original no encontrado: %s. No se puede generar el archivo limpio.\n', faaOriginal);
                else
                    try
                        [headers, sequences] = readFasta(faaOriginal);
                        fidOut = fopen(faaNuevo, 'w');
                        if fidOut == -1
                            error('No se pudo crear %s', faaNuevo);
                        end
                        
                        for i = 1:length(headers)
                            header = headers{i};
                            protein_id = strtok(header);
                            if proteinID_to_locusTag.isKey(protein_id)
                                newHeader = proteinID_to_locusTag(protein_id);
                            else
                                newHeader = protein_id;
                                fprintf('  Advertencia: No se encontró locus_tag para %s, se usa ID original\n', protein_id);
                            end
                            fprintf(fidOut, '>%s\n', newHeader);
                            seq = sequences{i};
                            for j = 1:70:length(seq)
                                fprintf(fidOut, '%s\n', seq(j:min(j+69, end)));
                            end
                        end
                        fclose(fidOut);
                        fprintf('✓ Archivo FASTA limpio generado: %s\n', faaNuevo);
                    catch ME
                        fprintf('❌ Error procesando FASTA: %s\n', ME.message);
                    end
                end
            else
                fprintf('\n  ℹ️  Existieron %d genes no encontrados (>=10). No se genera FASTA limpio.\n', numNotFound);
            end
            
            % Almacenar en el resumen global
            resumenModelos(end+1, :) = {modelName, totalGenes, numFound, numNotFound};
        end
        
    catch ME
        fprintf('❌ Error procesando %s: %s\n', filename, ME.message);
        [~, errName, ~] = fileparts(filename);
        resumenModelos(end+1, :) = {errName, NaN, NaN, NaN};
        modelosSinReducido{end+1} = errName;
    end
end

% --------------------------------------------------------------
% RESUMEN GLOBAL DE COMPARACIÓN
% --------------------------------------------------------------
fprintf('\n==================================================\n');
fprintf('           RESUMEN GLOBAL DE COMPARACIÓN\n');
fprintf('==================================================\n');

if ~isempty(resumenModelos)
    resumenTabla = cell2table(resumenModelos, ...
        'VariableNames', {'Modelo', 'TotalGenes', 'Encontrados', 'NoEncontrados'});
    
    disp(resumenTabla);
    
    idxValidos = ~isnan(resumenTabla.NoEncontrados);
    modelosConFaltantes = resumenTabla(idxValidos & resumenTabla.NoEncontrados > 0, :);
    
    if height(modelosConFaltantes) > 0
        fprintf('\n>>> MODELOS CON GENES NO ENCONTRADOS EN GTF <<<\n');
        for i = 1:height(modelosConFaltantes)
            fprintf('  • %s: %d genes no encontrados (de %d totales, %.1f%% éxito)\n', ...
                modelosConFaltantes.Modelo{i}, ...
                modelosConFaltantes.NoEncontrados(i), ...
                modelosConFaltantes.TotalGenes(i), ...
                (1 - modelosConFaltantes.NoEncontrados(i)/modelosConFaltantes.TotalGenes(i))*100);
        end
    else
        fprintf('\n✓ Todos los modelos procesados tienen TODOS sus genes en el GTF correspondiente.\n');
    end
    
    idxSinGTF = isnan(resumenTabla.NoEncontrados);
    if any(idxSinGTF)
        fprintf('\n>>> MODELOS SIN ARCHIVO GTF (no se pudo comparar) <<<\n');
        for i = find(idxSinGTF)'
            fprintf('  • %s\n', resumenTabla.Modelo{i});
        end
    end
else
    fprintf('No se procesó ningún modelo.\n');
end

% --------------------------------------------------------------
% LISTADO DE MODELOS QUE NO GENERARON EL ARCHIVO REDUCIDO
% --------------------------------------------------------------
if ~isempty(modelosSinReducido)
    fprintf('\n==================================================\n');
    fprintf('   MODELOS SIN ARCHIVO _protein_clean_reduced.faa\n');
    fprintf('==================================================\n');
    modelosUnicos = unique(modelosSinReducido);
    for i = 1:length(modelosUnicos)
        fprintf('  • %s\n', modelosUnicos{i});
    end
else
    fprintf('\n✓ Todos los modelos generaron correctamente su archivo _protein_clean_reduced.faa.\n');
end

fprintf('\n=== Procesamiento completado ===\n');

% --------------------------------------------------------------
% FUNCIÓN AUXILIAR PARA LEER FASTA (sin cambios)
% --------------------------------------------------------------
function [headers, sequences] = readFasta(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('No se puede abrir %s', filename);
    end
    headers = {};
    sequences = {};
    currentHeader = '';
    currentSeq = {};
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line)
            continue;
        end
        if line(1) == '>'
            if ~isempty(currentHeader)
                headers{end+1} = currentHeader;
                sequences{end+1} = strjoin(currentSeq, '');
            end
            currentHeader = line(2:end);
            currentSeq = {};
        else
            currentSeq{end+1} = line;
        end
    end
    if ~isempty(currentHeader)
        headers{end+1} = currentHeader;
        sequences{end+1} = strjoin(currentSeq, '');
    end
    fclose(fid);
end