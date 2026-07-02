%Para procesar específicamente los modelos iANig1029.mat y iAP1008.mat y generar archivos FASTA limpios y reducidos basados en los genes de cada modelo

% OBJETIVO: generar _protein_clean_1.faa, _protein_clean_reduced.faa con cabeceras de locus_tag (EL MISMO NOMBRE QUE EXISTE DENTRO DE LOS .MAT)

initCobraToolbox()

%% CONFIGURACIÓN – Lista de modelos a procesar (sin extensión .mat)
modelos = {'iANig1029', 'iAP1008'};   % Añadir aquí los modelos deseados

for idx = 1:length(modelos)
    baseName = modelos{idx};
    modelo_mat = [baseName '.mat'];
    fprintf('\n==================================================\n');
    fprintf('=== PROCESANDO MODELO: %s ===\n', baseName);
    fprintf('==================================================\n');
    
    faaFile = [baseName '_protein.faa'];
    csvFile = [baseName '_genes.csv'];
    
    %% 1. Cargar modelo .mat y exportar genes a CSV
    fprintf('\n--- 1. Cargando modelo y exportando genes ---\n');
    if ~exist(modelo_mat, 'file')
        fprintf('⚠️ Archivo %s no encontrado. Saltando...\n', modelo_mat);
        continue;
    end
    load(modelo_mat);
    vars = whos('-file', modelo_mat);
    model = [];
    for i = 1:length(vars)
        tmp = load(modelo_mat, vars(i).name);
        if isstruct(tmp.(vars(i).name)) && (isfield(tmp.(vars(i).name), 'genes') || isfield(tmp.(vars(i).name), 'rxns'))
            model = tmp.(vars(i).name);
            break;
        end
    end
    if isempty(model)
        fprintf('⚠️ No se encontró modelo COBRA en %s. Saltando...\n', modelo_mat);
        continue;
    end
    
    if isfield(model, 'geneNames') && ~isempty(model.geneNames)
        geneIDs = model.genes;
    elseif isfield(model, 'genes')
        geneIDs = model.genes;
    else
        fprintf('⚠️ El modelo no tiene campo genes. Saltando...\n');
        continue;
    end
    geneIDs = strtrim(cellstr(geneIDs));
    fprintf('Total genes en modelo: %d\n', length(geneIDs));
    
    % Exportar CSV
    geneTable = table(geneIDs, 'VariableNames', {'GeneID'});
    writetable(geneTable, csvFile);
    fprintf('✓ Genes exportados a %s\n', csvFile);
    
    %% 2. Leer el FASTA original (sin transformación)
    if ~exist(faaFile, 'file')
        fprintf('⚠️ Archivo FASTA %s no encontrado. Saltando...\n', faaFile);
        continue;
    end
    fprintf('\n--- 2. Procesando archivo FAA (sin transformación de IDs) ---\n');
    [headers, seqs] = readFasta(faaFile);
    total_proteins = length(headers);
    
    % Archivo de salida completo (todas las proteínas, cabeceras originales)
    output_clean1 = [baseName '_protein_clean_1.faa'];
    fid1 = fopen(output_clean1, 'w');
    if fid1 == -1
        fprintf('⚠️ No se pudo crear %s. Saltando...\n', output_clean1);
        continue;
    end
    
    % Archivo de salida reducido (solo proteínas con genes en el modelo)
    output_reduced = [baseName '_protein_clean_reduced.faa'];
    fidR = fopen(output_reduced, 'w');
    if fidR == -1
        fclose(fid1);
        fprintf('⚠️ No se pudo crear %s. Saltando...\n', output_reduced);
        continue;
    end
    
    % Conjunto de genes del modelo para búsqueda rápida
    model_gene_set = containers.Map(geneIDs, ones(length(geneIDs),1));
    
    kept = 0;      % proteínas que coinciden con algún gen del modelo
    eliminated = 0;
    
    % Guardar IDs de todas las proteínas (para comparación posterior)
    all_ids = cell(total_proteins, 1);
    
    for i = 1:total_proteins
        header = headers{i};
        % Extraer el primer campo (ID) antes del primer espacio
        id_original = strtok(header);
        all_ids{i} = id_original;
        
        % Escribir en el archivo completo (clean_1)
        fprintf(fid1, '>%s\n', header);   % la cabecera original completa
        seq = seqs{i};
        for j = 1:70:length(seq)
            fprintf(fid1, '%s\n', seq(j:min(j+69, end)));
        end
        
        % Escribir en el reducido solo si el ID está en los genes del modelo
        if model_gene_set.isKey(id_original)
            fprintf(fidR, '>%s\n', header);
            for j = 1:70:length(seq)
                fprintf(fidR, '%s\n', seq(j:min(j+69, end)));
            end
            kept = kept + 1;
        else
            eliminated = eliminated + 1;
        end
    end
    fclose(fid1);
    fclose(fidR);
    
    fprintf('✓ Archivo guardado: %s\n', output_clean1);
    fprintf('✓ Archivo REDUCIDO guardado: %s\n', output_reduced);
    fprintf('Total proteínas: %d\n', total_proteins);
    
    % Estadísticas de filtrado
    fprintf('\n--- ESTADÍSTICAS DE FILTRADO PARA %s ---\n', baseName);
    fprintf('Total de proteínas en FAA original: %d\n', total_proteins);
    fprintf('Proteínas MANTENIDAS (en modelo): %d\n', kept);
    fprintf('Proteínas ELIMINADAS (no en modelo): %d\n', eliminated);
    fprintf('Porcentaje mantenido: %.2f%%\n', kept/total_proteins*100);
    fprintf('Porcentaje eliminado: %.2f%%\n', eliminated/total_proteins*100);
    
    %% 3. Verificar qué genes del modelo están presentes en el FASTA original (comparación)
    fprintf('\n--- 3. Comparación con genes del modelo ---\n');
    unique_ids = unique(all_ids);
    fprintf('IDs únicos en FASTA original: %d\n', length(unique_ids));
    
    % Comparar con los geneIDs del modelo
    found_in_fasta = ismember(geneIDs, unique_ids);
    num_found = sum(found_in_fasta);
    num_not_found = sum(~found_in_fasta);
    
    fprintf('Genes del modelo encontrados en FASTA: %d\n', num_found);
    fprintf('Genes del modelo NO encontrados en FASTA: %d\n', num_not_found);
    if num_not_found > 0
        not_found_genes = geneIDs(~found_in_fasta);
        fprintf('\nEjemplos de genes NO encontrados (primeros 10):\n');
        for k = 1:min(10, length(not_found_genes))
            fprintf('  %s\n', not_found_genes{k});
        end
    end
    
    % Ejemplos de genes encontrados
    found_genes = geneIDs(found_in_fasta);
    if ~isempty(found_genes)
        fprintf('\nEjemplos de genes encontrados (primeros 10):\n');
        for k = 1:min(10, length(found_genes))
            fprintf('  %s\n', found_genes{k});
        end
    end
    
    % Estadísticas finales
    fprintf('\n=== ESTADÍSTICAS FINALES ===\n');
    fprintf('Total genes en modelo: %d\n', length(geneIDs));
    fprintf('Presentes en FASTA: %d (%.1f%%)\n', num_found, num_found/length(geneIDs)*100);
    fprintf('Ausentes en FASTA: %d (%.1f%%)\n', num_not_found, num_not_found/length(geneIDs)*100);
    
    fprintf('\n--- Procesamiento de %s completado ---\n', baseName);
end

fprintf('\n=== PROCESAMIENTO DE TODOS LOS MODELOS COMPLETADO ===\n');

%% FUNCIONES AUXILIARES
function [headers, sequences] = readFasta(filename)
    fid = fopen(filename, 'r');
    headers = {}; sequences = {};
    curHead = ''; curSeq = {};
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line), continue; end
        if line(1) == '>'
            if ~isempty(curHead)
                headers{end+1} = curHead;
                sequences{end+1} = strjoin(curSeq, '');
            end
            curHead = line(2:end);
            curSeq = {};
        else
            curSeq{end+1} = line;
        end
    end
    if ~isempty(curHead)
        headers{end+1} = curHead;
        sequences{end+1} = strjoin(curSeq, '');
    end
    fclose(fid);
end