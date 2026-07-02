%Para procesar específicamente el modelo iANid1221.mat y generar un archivo FASTA limpio y reducido basado en los genes del modelo

% OBJETIVO: generar _protein_clean_1.faa, _protein_clean_reduced.faa con cabeceras de locus_tag (EL MISMO NOMBRE QUE EXISTE DENTRO DE LOS .MAT)

initCobraToolbox()

%% CONFIGURACIÓN – Procesar solo el modelo iANid1221
modelo_mat = 'iANid1221.mat';
[~, baseName, ~] = fileparts(modelo_mat);
fprintf('=== PROCESANDO MODELO: %s ===\n', baseName);

faaFile = [baseName '_protein.faa'];
csvFile = [baseName '_genes.csv'];      % se generará desde el modelo

%% 1. Cargar modelo .mat y exportar genes a CSV
fprintf('\n--- 1. Cargando modelo y exportando genes ---\n');
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
    error('No se encontró modelo COBRA en %s', modelo_mat);
end

if isfield(model, 'geneNames') && ~isempty(model.geneNames)
    geneIDs = model.genes;
elseif isfield(model, 'genes')
    geneIDs = model.genes;
else
    error('El modelo no tiene campo genes');
end
geneIDs = strtrim(cellstr(geneIDs));
fprintf('Total genes en modelo: %d\n', length(geneIDs));

% Exportar CSV (opcional)
geneTable = table(geneIDs, 'VariableNames', {'GeneID'});
writetable(geneTable, csvFile);
fprintf('✓ Genes exportados a %s\n', csvFile);

%% 2. Leer el FASTA original y transformar las cabeceras según la regla ANIA_ + 5 dígitos
fprintf('\n--- 2. Procesando archivo FAA y transformando identificadores ---\n');
[headers, seqs] = readFasta(faaFile);
total_proteins = length(headers);

% Función para transformar un identificador (ej. AN0041 -> ANIA_00041)
transform_id = @(id) transform_identifier(id);

% Contadores
transformed = 0;
unchanged = 0;

% Archivo de salida completo (todas las proteínas transformadas)
output_clean1 = [baseName '_protein_clean_1.faa'];
fid1 = fopen(output_clean1, 'w');
if fid1 == -1
    error('No se pudo crear %s', output_clean1);
end

% Archivo de salida reducido (solo proteínas con genes en el modelo)
output_reduced = [baseName '_protein_clean_reduced.faa'];
fidR = fopen(output_reduced, 'w');
if fidR == -1
    error('No se pudo crear %s', output_reduced);
end

% Guardar también los nuevos IDs en una celda para su posterior comparación
new_ids = cell(total_proteins, 1);

% Conjunto de genes del modelo para búsqueda rápida
model_gene_set = containers.Map(geneIDs, ones(length(geneIDs),1));

kept = 0;      % contador para reducido
eliminated = 0;

for i = 1:total_proteins
    header = headers{i};
    id_original = strtok(header);
    new_id = transform_id(id_original);
    new_ids{i} = new_id;
    if ~strcmp(new_id, id_original)
        transformed = transformed + 1;
    else
        unchanged = unchanged + 1;
    end
    % Mantener el resto del header después del ID original
    rest_of_header = header(length(id_original)+1:end);
    
    % Escribir en el archivo completo (clean_1)
    fprintf(fid1, '>%s%s\n', new_id, rest_of_header);
    seq = seqs{i};
    for j = 1:70:length(seq)
        fprintf(fid1, '%s\n', seq(j:min(j+69, end)));
    end
    
    % Escribir en el archivo reducido solo si el nuevo ID está en los genes del modelo
    if model_gene_set.isKey(new_id)
        fprintf(fidR, '>%s%s\n', new_id, rest_of_header);
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

fprintf('Archivo guardado: %s_protein_clean_1.faa\n', baseName);
fprintf('Archivo REDUCIDO guardado: %s_protein_clean_reduced.faa\n', baseName);
fprintf('Total proteínas: %d\n', total_proteins);
fprintf('Identificadores transformados: %d\n', transformed);
fprintf('Identificadores sin cambios: %d\n', unchanged);
fprintf('Tasa de transformación: %.2f%%\n', transformed/total_proteins*100);

% Estadísticas de filtrado
fprintf('\n--- ESTADÍSTICAS DE FILTRADO PARA %s ---\n', baseName);
fprintf('Total de proteínas en FAA original: %d\n', total_proteins);
fprintf('Proteínas MANTENIDAS (en modelo): %d\n', kept);
fprintf('Proteínas ELIMINADAS (no en modelo): %d\n', eliminated);
fprintf('Porcentaje mantenido: %.2f%%\n', kept/total_proteins*100);
fprintf('Porcentaje eliminado: %.2f%%\n', eliminated/total_proteins*100);

%% 3. Verificar qué genes del modelo están presentes en el nuevo FASTA (comparación adicional)
fprintf('\n--- 3. Comparación con genes del modelo ---\n');
% Obtener lista única de IDs transformados (los que aparecen en el FASTA limpio)
unique_new_ids = unique(new_ids);
fprintf('IDs únicos en FASTA limpio: %d\n', length(unique_new_ids));

% Comparar con los geneIDs del modelo
model_geneIDs = geneIDs;
found_in_fasta = ismember(model_geneIDs, unique_new_ids);
num_found = sum(found_in_fasta);
num_not_found = sum(~found_in_fasta);

fprintf('Genes del modelo encontrados en FASTA limpio: %d\n', num_found);
fprintf('Genes del modelo NO encontrados en FASTA limpio: %d\n', num_not_found);
if num_not_found > 0
    not_found_genes = model_geneIDs(~found_in_fasta);
    fprintf('\nEjemplos de genes NO encontrados (primeros 10):\n');
    for k = 1:min(10, length(not_found_genes))
        fprintf('  %s\n', not_found_genes{k});
    end
end

% Ejemplos de genes encontrados
found_genes = model_geneIDs(found_in_fasta);
if ~isempty(found_genes)
    fprintf('\nEjemplos de genes encontrados (primeros 10):\n');
    for k = 1:min(10, length(found_genes))
        fprintf('  %s\n', found_genes{k});
    end
end

% Estadísticas finales
fprintf('\n=== ESTADÍSTICAS FINALES ===\n');
fprintf('Total genes en modelo: %d\n', length(model_geneIDs));
fprintf('Presentes en FASTA limpio: %d (%.1f%%)\n', num_found, num_found/length(model_geneIDs)*100);
fprintf('Ausentes en FASTA limpio: %d (%.1f%%)\n', num_not_found, num_not_found/length(model_geneIDs)*100);

fprintf('\n=== PROCESO COMPLETADO ===\n');

%% FUNCIONES AUXILIARES
function new_id = transform_identifier(id)
    % Transforma identificadores como AN0041 -> ANIA_00041, AN6188 -> ANIA_06188
    pattern = '^([A-Za-z]+)(\d+)$';
    tokens = regexp(id, pattern, 'tokens');
    if isempty(tokens)
        new_id = id;
        return;
    end
    num_str = tokens{1}{2};
    num = str2double(num_str);
    if isnan(num)
        new_id = id;
        return;
    end
    new_id = sprintf('ANIA_%05d', num);
end

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