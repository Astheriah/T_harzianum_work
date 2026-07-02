%Para procesar específicamente el modelo iPrub22.mat y generar un archivo FASTA limpio y reducido basado en los genes del modelo

% OBJETIVO: generar _protein_clean_1.faa, _protein_clean_reduced.faa con cabeceras de locus_tag (EL MISMO NOMBRE QUE EXISTE DENTRO DE LOS .MAT)

initCobraToolbox()

%% CONFIGURACIÓN – Procesar solo el modelo iPrub22
modelo_mat = 'iPrub22.mat';
[~, baseName, ~] = fileparts(modelo_mat);
fprintf('=== PROCESANDO MODELO: %s ===\n', baseName);

gtfFile = [baseName '_gtf.gtf'];
faaFile = [baseName '_protein.faa'];
csvFile = [baseName '_genes.csv'];      

%% 1. Cargar modelo .mat y obtener los GeneID
fprintf('\n--- 1. Cargando modelo y extrayendo genes ---\n');
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

% Exportar CSV (opcional, con una sola columna GeneID)
geneTable = table(geneIDs, 'VariableNames', {'GeneID'});
writetable(geneTable, csvFile);
fprintf('✓ Genes exportados a %s\n', csvFile);

%% 2. Leer GTF y construir mapeos
fprintf('\n--- 2. Leyendo archivo GTF ---\n');
gtf = readtable(gtfFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'CommentStyle', '#', 'ReadVariableNames', false);
gtf.Properties.VariableNames = {'seqname','source','feature','start','end','score','strand','frame','attributes'};

protein_to_geneid = containers.Map();   % protein_id -> gene_id (sin prefijo)
gene_info_map = containers.Map();       % gene_id -> struct con metadata

for i = 1:height(gtf)
    if ~ismember(gtf.feature{i}, {'CDS','gene','transcript'}), continue; end
    attrs = gtf.attributes{i};
    prot_id = extract_attr(attrs, 'protein_id');
    gene_id = extract_attr(attrs, 'gene_id');
    if isempty(prot_id) || isempty(gene_id), continue; end
    
    protein_to_geneid(prot_id) = gene_id;
    
    gene_name = extract_attr(attrs, 'gene');
    locus_tag = extract_attr(attrs, 'locus_tag');
    product   = extract_attr(attrs, 'product');
    info = struct('gene_name', gene_name, 'locus_tag', locus_tag, 'product', product, 'protein_id', prot_id);
    gene_info_map(gene_id) = info;
end
fprintf('Se encontraron %d mapeos de protein_id a gen\n', protein_to_geneid.Count);
fprintf('Se encontraron %d genes únicos en GTF\n', gene_info_map.Count);

prot_keys = protein_to_geneid.keys;
fprintf('Ejemplos de mapeos encontrados:\n');
for k = 1:min(5, length(prot_keys))
    prot = prot_keys{k};
    gene = protein_to_geneid(prot);
    if gene_info_map.isKey(gene)
        info = gene_info_map(gene);
        fprintf('  %s -> gene_id: %s, locus_tag: %s, product: %s\n', ...
            prot, gene, info.locus_tag, info.product);
    else
        fprintf('  %s -> %s\n', prot, gene);
    end
end

%% 3. Leer el FASTA original y generar estadísticas
fprintf('\n--- 3. Procesando archivo FAA ---\n');
[headers, seqs] = readFasta(faaFile);
total_proteins = length(headers);

replaced = 0;
for i = 1:total_proteins
    prot_id = strtok(headers{i});
    if protein_to_geneid.isKey(prot_id)
        replaced = replaced + 1;
    end
end
fprintf('Total de proteínas en FAA: %d\n', total_proteins);
fprintf('Proteínas con reemplazo: %d\n', replaced);
fprintf('Proteínas sin reemplazo: %d\n', total_proteins - replaced);
fprintf('Tasa de reemplazo: %.2f%%\n', replaced/total_proteins*100);

%% 4. Preparar la lista de genes del modelo (transformados y originales)
% Los geneIDs del modelo pueden tener prefijos como 'gp_', 'gnl_', etc.
% Necesitamos la versión sin prefijo para buscar en el GTF.
transform = @(x) regexprep(x, '^.*?_', '');   % elimina todo hasta el primer '_'
transformed_targets = cellfun(transform, geneIDs, 'UniformOutput', false);
% Lista única de genes a buscar (sin prefijo)
target_names = unique(transformed_targets);
fprintf('Se creó un mapeo de %d genes del modelo (transformados) a GeneID original.\n', length(target_names));
if length(target_names) > 0
    fprintf('Primeros 5 genes transformados: %s\n', strjoin(target_names(1:min(5,end))', ', '));
end

% Mapeo de nombre transformado -> ID original (con prefijo)
transformed_to_original = containers.Map();
for i = 1:length(transformed_targets)
    t = transformed_targets{i};
    o = geneIDs{i};
    if ~transformed_to_original.isKey(t)
        transformed_to_original(t) = o;
    end
end

%% 5. Generar archivo _clean_1.faa con cabeceras = GeneID original (cuando coincide)
fprintf('\n--- 5. Generando archivo con GeneIDs originales -> %s_protein_clean_1.faa ---\n', baseName);
output_clean1 = [baseName '_protein_clean_1.faa'];
fid1 = fopen(output_clean1, 'w');
if fid1 == -1, error('No se pudo crear %s', output_clean1); end

replaced_clean1 = 0;
for i = 1:total_proteins
    prot_id = strtok(headers{i});
    if protein_to_geneid.isKey(prot_id)
        gene_id = protein_to_geneid(prot_id);   % gene_id del GTF (ya sin prefijo)
        % Verificar si este gene_id está en nuestra lista de genes transformados
        if ismember(gene_id, target_names)
            % Obtener el ID original del modelo (con prefijo)
            if transformed_to_original.isKey(gene_id)
                newHeader = transformed_to_original(gene_id);
            else
                newHeader = gene_id;   % fallback (no debería ocurrir)
            end
            replaced_clean1 = replaced_clean1 + 1;
        else
            newHeader = gene_id;   % usar el gene_id del GTF (sin prefijo)
        end
    else
        newHeader = prot_id;   % fallback (protein_id original)
    end
    fprintf(fid1, '>%s\n', newHeader);
    seq = seqs{i};
    for j = 1:70:length(seq)
        fprintf(fid1, '%s\n', seq(j:min(j+69, end)));
    end
end
fclose(fid1);
fprintf('Archivo con GeneIDs guardado como: %s_protein_clean_1.faa\n', baseName);
fprintf('Se reemplazaron %d headers por GeneID original del modelo.\n', replaced_clean1);

%% 5b. Generar archivo REDUCIDO (solo proteínas de genes en CSV)
fprintf('\n--- Generando archivo REDUCIDO (solo proteínas de genes en CSV) ---\n');
output_reduced = [baseName '_protein_clean_reduced.faa'];
fidR = fopen(output_reduced, 'w');
if fidR == -1, error('No se pudo crear %s', output_reduced); end

kept = 0;
eliminated = 0;
for i = 1:total_proteins
    prot_id = strtok(headers{i});
    if protein_to_geneid.isKey(prot_id)
        gene_id = protein_to_geneid(prot_id);
        if ismember(gene_id, target_names)
            % Proteína pertenece a un gen del CSV -> escribir en reducido
            if transformed_to_original.isKey(gene_id)
                newHeader = transformed_to_original(gene_id);   % ID original (con prefijo)
            else
                newHeader = gene_id;
            end
            fprintf(fidR, '>%s\n', newHeader);
            seq = seqs{i};
            for j = 1:70:length(seq)
                fprintf(fidR, '%s\n', seq(j:min(j+69, end)));
            end
            kept = kept + 1;
        else
            eliminated = eliminated + 1;
        end
    else
        eliminated = eliminated + 1;
    end
end
fclose(fidR);
fprintf('✓ Archivo REDUCIDO guardado como: %s\n', output_reduced);
fprintf('Total proteínas procesadas: %d\n', total_proteins);
fprintf('\n--- ESTADÍSTICAS DE FILTRADO PARA %s ---\n', baseName);
fprintf('Total de proteínas en FAA original: %d\n', total_proteins);
fprintf('Proteínas MANTENIDAS (en CSV): %d\n', kept);
fprintf('Proteínas ELIMINADAS (no en CSV): %d\n', eliminated);
fprintf('Porcentaje mantenido: %.2f%%\n', kept/total_proteins*100);
fprintf('Porcentaje eliminado: %.2f%%\n', eliminated/total_proteins*100);

%% 6. Buscar los genes transformados en GTF y mostrar estadísticas
fprintf('\n--- RESULTADOS EN GTF ---\n');
gtf_gene_ids = gene_info_map.keys;   % todos los gene_id del GTF (sin prefijo)
found_in_gtf = {};
not_found_in_gtf = {};
for i = 1:length(target_names)
    name = target_names{i};
    if ismember(name, gtf_gene_ids)
        found_in_gtf{end+1} = name;
    else
        not_found_in_gtf{end+1} = name;
    end
end
fprintf('Genes encontrados en GTF: %d\n', length(found_in_gtf));
fprintf('Genes NO encontrados en GTF: %d\n', length(not_found_in_gtf));

if ~isempty(found_in_gtf)
    fprintf('\nEjemplos de genes encontrados en GTF:\n');
    for k = 1:min(5, length(found_in_gtf))
        name = found_in_gtf{k};
        info = gene_info_map(name);
        fprintf('  %s -> gene_name: ''%s'', locus_tag: ''%s'', product: ''%s'', protein_id: ''%s''\n', ...
            name, info.gene_name, info.locus_tag, info.product, info.protein_id);
    end
end
if ~isempty(not_found_in_gtf)
    fprintf('\nEjemplos de genes NO encontrados en GTF:\n');
    for k = 1:min(5, length(not_found_in_gtf))
        fprintf('  %s\n', not_found_in_gtf{k});
    end
end

%% 7. Búsqueda en FAA (usando protein_id)
fprintf('\n--- RESULTADOS EN FAA ---\n');
gene_to_protein = containers.Map();
gene_keys_list = gene_info_map.keys;
for k = 1:length(gene_keys_list)
    g = gene_keys_list{k};
    info = gene_info_map(g);
    if ~isempty(info.protein_id)
        gene_to_protein(g) = info.protein_id;
    end
end

faa_protein_ids = cellfun(@(h) strtok(h), headers, 'UniformOutput', false);
faa_set = containers.Map(faa_protein_ids, ones(length(faa_protein_ids),1));

found_in_faa = {};
not_found_in_faa = {};
for i = 1:length(found_in_gtf)
    name = found_in_gtf{i};
    if gene_to_protein.isKey(name)
        prot_id = gene_to_protein(name);
        if faa_set.isKey(prot_id)
            found_in_faa{end+1} = name;
        else
            not_found_in_faa{end+1} = name;
        end
    else
        not_found_in_faa{end+1} = name;
    end
end
fprintf('Genes encontrados en FAA: %d\n', length(found_in_faa));
fprintf('Genes NO encontrados en FAA: %d\n', length(not_found_in_faa));

if ~isempty(found_in_faa)
    fprintf('\nEjemplos de genes encontrados en FAA:\n');
    for k = 1:min(5, length(found_in_faa))
        name = found_in_faa{k};
        prot_id = gene_to_protein(name);
        fprintf('  %s -> %s\n', name, prot_id);
    end
end

%% 8. Estadísticas finales
fprintf('\n=== ESTADÍSTICAS FINALES ===\n');
total_names = length(target_names);
fprintf('Total genes en modelo (transformados): %d\n', total_names);
fprintf('Encontrados en GTF: %d (%.1f%%)\n', length(found_in_gtf), length(found_in_gtf)/total_names*100);
fprintf('Encontrados en FAA: %d (%.1f%%)\n', length(found_in_faa), length(found_in_faa)/total_names*100);

fprintf('\n=== PROCESO COMPLETADO ===\n');

%% FUNCIONES AUXILIARES
function val = extract_attr(attrsStr, key)
    pattern = [key ' "([^"]+)"'];
    tokens = regexp(attrsStr, pattern, 'tokens');
    if ~isempty(tokens), val = tokens{1}{1}; else, val = ''; end
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