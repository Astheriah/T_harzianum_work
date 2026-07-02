%Para procesar especificamente el modelo Rt_IFO0880.mat y generar un archivo FASTA limpio y reducido basado en los genes del modelo

% OBJETIVO: generar _protein_clean_1.faa, _protein_clean_reduced.faa con cabeceras de locus_tag (EL MISMO NOMBRE QUE EXISTE DENTRO DE LOS .MAT)

initCobraToolbox()

modelo_mat = 'Rt_IFO0880.mat';
[~, baseName, ~] = fileparts(modelo_mat);
fprintf('=== PROCESANDO MODELO: %s ===\n', baseName);

gtfFile = [baseName '_gtf.gtf'];
faaFile = [baseName '_protein.faa'];
csvFile = [baseName '_genes.csv'];   % Archivo CSV con los GeneID del modelo

%% Leer GTF y construir mapeo protein_id -> gene_id y número
fprintf('\n--- Leyendo GTF y extrayendo mapeo protein_id -> gene_id ---\n');
gtf = readtable(gtfFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'CommentStyle', '#', 'ReadVariableNames', false);
gtf.Properties.VariableNames = {'seqname','source','feature','start','end','score','strand','frame','attributes'};

protein_to_geneid = containers.Map();   % protein_id -> gene_id completo
protein_to_number = containers.Map();   % protein_id -> número extraído

for i = 1:height(gtf)
    if ~ismember(gtf.feature{i}, {'CDS','gene','transcript'}), continue; end
    attrs = gtf.attributes{i};
    prot_id = extract_attr(attrs, 'protein_id');
    gene_id = extract_attr(attrs, 'gene_id');
    if isempty(prot_id) || isempty(gene_id), continue; end
    
    protein_to_geneid(prot_id) = gene_id;
    
    % Extraer el número del gene_id (último bloque numérico)
    num = extract_number(gene_id);
    if ~isempty(num)
        protein_to_number(prot_id) = num;
    end
end
fprintf('Mapeos protein_id -> gene_id: %d\n', protein_to_geneid.Count);
fprintf('Mapeos protein_id -> número: %d\n', protein_to_number.Count);

% Ejemplos
keys = protein_to_geneid.keys;
fprintf('Ejemplos de mapeos:\n');
for k = 1:min(5, length(keys))
    prot = keys{k};
    gene = protein_to_geneid(prot);
    if protein_to_number.isKey(prot)
        num = protein_to_number(prot);
        fprintf('  %s -> %s [número: %s]\n', prot, gene, num);
    else
        fprintf('  %s -> %s [número: ?]\n', prot, gene);
    end
end

%% Generar archivo _protein_clean_1.faa con cabeceras = número (cuando sea posible)
fprintf('\n--- Generando %s_protein_clean_1.faa (cabeceras = número) ---\n', baseName);
[headers, seqs] = readFasta(faaFile);
fid = fopen([baseName '_protein_clean_1.faa'], 'w');
replaced_num = 0;
replaced_gene = 0;
total = length(headers);

for i = 1:total
    prot_id = strtok(headers{i});
    if protein_to_number.isKey(prot_id)
        newHeader = protein_to_number(prot_id);   % usar solo el número
        replaced_num = replaced_num + 1;
    elseif protein_to_geneid.isKey(prot_id)
        newHeader = protein_to_geneid(prot_id);   % fallback: gene_id completo
        replaced_gene = replaced_gene + 1;
    else
        newHeader = prot_id;                      % fallback: protein_id original
    end
    fprintf(fid, '>%s\n', newHeader);
    seq = seqs{i};
    for j = 1:70:length(seq)
        fprintf(fid, '%s\n', seq(j:min(j+69, end)));
    end
end
fclose(fid);
fprintf('Archivo guardado: %s_protein_clean_1.faa\n', baseName);
fprintf('Total proteínas: %d\n', total);
fprintf('  - Cabeceras reemplazadas por número: %d (%.2f%%)\n', replaced_num, replaced_num/total*100);
fprintf('  - Cabeceras reemplazadas por gene_id: %d\n', replaced_gene);
fprintf('  - Cabeceras sin mapeo (protein_id original): %d\n', total - replaced_num - replaced_gene);

%% Leer CSV de genes (Rt_IFO0880_genes.csv) y extraer números de GeneID
fprintf('\n--- Leyendo CSV de genes: %s ---\n', csvFile);
if ~exist(csvFile, 'file')
    error('No se encontró el archivo CSV: %s', csvFile);
end

% Intentar leer CSV (puede tener tabulador o coma)
try
    genesTable = readtable(csvFile, 'FileType', 'text', 'Delimiter', '\t', 'PreserveVariableNames', true);
    fprintf('✓ CSV leído con separador TABULADOR\n');
catch
    try
        genesTable = readtable(csvFile, 'FileType', 'text', 'Delimiter', ',', 'PreserveVariableNames', true);
        fprintf('✓ CSV leído con separador COMA\n');
    catch ME
        error('No se pudo leer el CSV: %s', ME.message);
    end
end

% Determinar columna de GeneID
if any(strcmp(genesTable.Properties.VariableNames, 'GeneID'))
    geneCol = 'GeneID';
elseif any(strcmp(genesTable.Properties.VariableNames, 'gene_id'))
    geneCol = 'gene_id';
else
    geneCol = genesTable.Properties.VariableNames{1};
    fprintf('Usando columna: %s\n', geneCol);
end

% Extraer números de la columna GeneID
rawVals = table2cell(genesTable(:, geneCol));
target_numbers = {};
number_to_original = containers.Map();
for i = 1:length(rawVals)
    val = char(rawVals{i});
    nums = regexp(val, '\d+', 'match');
    if ~isempty(nums)
        num = nums{1};
        target_numbers{end+1} = num;
        number_to_original(num) = val;
    end
end
target_numbers = unique(target_numbers);
fprintf('Números únicos extraídos del CSV (GeneID): %d\n', length(target_numbers));
if ~isempty(target_numbers)
    fprintf('Primeros 10 números: %s\n', strjoin(target_numbers(1:min(10,end))', ', '));
end

%% Leer el archivo _protein_clean_1.faa y extraer los identificadores de cabecera
fprintf('\n--- Leyendo %s_protein_clean_1.faa para extraer cabeceras ---\n', baseName);
cleanFaaFile = [baseName '_protein_clean_1.faa'];
if ~exist(cleanFaaFile, 'file')
    error('No se encontró el archivo FASTA limpio: %s', cleanFaaFile);
end

[cleanHeaders, cleanSeqs] = readFasta(cleanFaaFile);  % ahora necesitamos también las secuencias para el reducido
header_ids = cellfun(@(h) strtrim(h), cleanHeaders, 'UniformOutput', false);
fprintf('Total de entradas en FASTA limpio: %d\n', length(header_ids));

% NUEVO: Generar archivo REDUCIDO solo con las entradas cuyos IDs están en target_numbers
fprintf('\n--- Generando archivo REDUCIDO (%s_protein_clean_reduced.faa) ---\n', baseName);
reducedFile = [baseName '_protein_clean_reduced.faa'];
fidRed = fopen(reducedFile, 'w');
if fidRed == -1
    error('No se pudo crear %s', reducedFile);
end

kept = 0;
discarded = 0;
for i = 1:length(header_ids)
    id = header_ids{i};
    if ismember(id, target_numbers)
        % Escribir esta entrada
        fprintf(fidRed, '>%s\n', id);
        seq = cleanSeqs{i};
        for j = 1:70:length(seq)
            fprintf(fidRed, '%s\n', seq(j:min(j+69, end)));
        end
        kept = kept + 1;
    else
        discarded = discarded + 1;
    end
end
fclose(fidRed);
fprintf('✓ Archivo REDUCIDO guardado como: %s\n', reducedFile);
fprintf('Total proteínas en FASTA limpio: %d\n', length(header_ids));
fprintf('Proteínas MANTENIDAS (en CSV): %d\n', kept);
fprintf('Proteínas ELIMINADAS (no en CSV): %d\n', discarded);
fprintf('Porcentaje mantenido: %.2f%%\n', kept/length(header_ids)*100);
fprintf('Porcentaje eliminado: %.2f%%\n', discarded/length(header_ids)*100);

%% Comparar los números del CSV con los IDs del FASTA limpio (ya no es necesario hacer la comparación por separado, pero la mantenemos por consistencia)
fprintf('\n--- Comparación: números del CSV vs IDs en FASTA limpio ---\n');
found_in_fasta = {};
not_found_in_fasta = {};
for i = 1:length(target_numbers)
    num = target_numbers{i};
    if ismember(num, header_ids)
        found_in_fasta{end+1} = num;
    else
        not_found_in_fasta{end+1} = num;
    end
end

fprintf('Números totales en CSV: %d\n', length(target_numbers));
fprintf('Encontrados en FASTA limpio: %d\n', length(found_in_fasta));
fprintf('NO encontrados en FASTA limpio: %d\n', length(not_found_in_fasta));
if ~isempty(not_found_in_fasta)
    fprintf('Porcentaje de éxito: %.1f%%\n', length(found_in_fasta)/length(target_numbers)*100);
    fprintf('\nPrimeros 10 números NO encontrados:\n');
    for k = 1:min(10, length(not_found_in_fasta))
        fprintf('  %s\n', not_found_in_fasta{k});
    end
else
    fprintf('✓ Todos los números del CSV están presentes en el FASTA limpio.\n');
end

fprintf('\n=== PROCESO COMPLETADO ===\n');

%% FUNCIONES AUXILIARES (sin cambios)
function val = extract_attr(attrsStr, key)
    pattern = [key ' "([^"]+)"'];
    tokens = regexp(attrsStr, pattern, 'tokens');
    if ~isempty(tokens)
        val = tokens{1}{1};
    else
        val = '';
    end
end

function num = extract_number(gene_id)
    tokens = regexp(gene_id, '\d+', 'match');
    if ~isempty(tokens)
        num = tokens{end};
    else
        num = '';
    end
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