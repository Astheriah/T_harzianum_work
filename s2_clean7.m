% PROCESAMIENTO PARA iDD1552.mat – VERSIÓN FINAL
% - Ignora genes que empiezan por "AB" (ej. AB1).
% - Para los demás genes (Gglean...), extrae el número, elimina ceros a la izquierda,
%   forma "LOCUS####" y busca ese sufijo al final del locus_tag del GTF.
% - Genera _protein_clean_1.faa y _protein_clean_reduced.faa con cabeceras = GeneID del CSV
%   cuando hay coincidencia; si no, mantiene el protein_id original.

initCobraToolbox()

%% CONFIGURACIÓN
modelo_mat = 'iDD1552.mat';
[~, baseName, ~] = fileparts(modelo_mat);
fprintf('=== PROCESANDO MODELO: %s ===\n', baseName);

gtfFile   = [baseName '_gtf.gtf'];
faaFile   = [baseName '_protein.faa'];
csvFile   = [baseName '_genes.csv'];

%% 1. Cargar modelo y exportar CSV de genes
fprintf('\n--- 1. Cargando modelo y generando CSV de genes ---\n');
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
if ~isfield(model, 'genes')
    error('El modelo no tiene campo genes');
end
geneIDs = strtrim(cellstr(model.genes));
fprintf('Total genes en modelo: %d\n', length(geneIDs));
geneTable = table(geneIDs, 'VariableNames', {'GeneID'});
writetable(geneTable, csvFile);
fprintf('✓ Genes exportados a %s\n', csvFile);

%% 2. Leer CSV y construir sufijos (ignorando genes que empiezan con "AB")
fprintf('\n--- 2. Leyendo CSV y preparando sufijos ---\n');
geneTable = readtable(csvFile, 'Delimiter', ',', 'ReadVariableNames', true);
targetGenes = strtrim(cellstr(geneTable.GeneID));
fprintf('Total genes en el CSV: %d\n', length(targetGenes));

suffix_to_geneid = containers.Map();   % 'LOCUS2178' -> 'Gglean002178'
genesIgnorados = 0;

for i = 1:length(targetGenes)
    geneName = targetGenes{i};

    % ---- IGNORAR SI EMPIEZA POR "AB" ----
    if startsWith(geneName, 'AB')
        genesIgnorados = genesIgnorados + 1;
        % fprintf('   Ignorado (empieza con AB): %s\n', geneName);   % descomentar para ver detalles
        continue;
    end

    numPart = regexp(geneName, '\d+', 'match', 'once');
    if isempty(numPart)
        warning('El gen %s no contiene dígitos, se omite.', geneName);
        continue;
    end
    numValue = str2double(numPart);  % elimina ceros a la izquierda
    locusSuffix = ['LOCUS' num2str(numValue)];
    suffix_to_geneid(locusSuffix) = geneName;
end

fprintf('Se construyeron %d sufijos únicos (ignorando %d genes que empiezan con AB).\n', ...
        suffix_to_geneid.Count, genesIgnorados);

%% 3. Procesar GTF: extraer sufijo del locus_tag y mapear protein_id -> GeneID
fprintf('\n--- 3. Leyendo GTF y vinculando mediante sufijo ---\n');
gtf = readtable(gtfFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'CommentStyle', '#', 'ReadVariableNames', false);
gtf.Properties.VariableNames = {'seqname','source','feature','start','end','score','strand','frame','attributes'};

protein_to_csvGene = containers.Map();   % 'CAI9627292.1' -> 'Gglean002178'
matchedSuffixes = containers.Map('KeyType','char','ValueType','logical');

for i = 1:height(gtf)
    if ~strcmp(gtf.feature{i}, 'CDS'), continue; end
    attrs = gtf.attributes{i};
    locus_val = extract_attr(attrs, 'locus_tag');
    protein_val = extract_attr(attrs, 'protein_id');
    if isempty(locus_val) || isempty(protein_val), continue; end

    % Obtener el sufijo (última parte tras '_')
    parts = strsplit(locus_val, '_');
    suffix = parts{end};   % p.ej. 'LOCUS2178'

    if suffix_to_geneid.isKey(suffix)
        geneName = suffix_to_geneid(suffix);
        % Guardar el mapeo (primera proteína que aparezca para ese sufijo)
        if ~protein_to_csvGene.isKey(protein_val)
            protein_to_csvGene(protein_val) = geneName;
            matchedSuffixes(suffix) = true;
        end
    end
end

allSuffixes = suffix_to_geneid.keys;
unmatchedSuffixes = setdiff(allSuffixes, matchedSuffixes.keys);
fprintf('Sufijos encontrados en GTF: %d\n', matchedSuffixes.Count);
fprintf('Sufijos NO encontrados en GTF: %d\n', length(unmatchedSuffixes));
if ~isempty(unmatchedSuffixes)
    fprintf('Ejemplos de sufijos no encontrados: %s\n', strjoin(unmatchedSuffixes(1:min(5,end)), ', '));
end
fprintf('Total proteínas mapeadas a genes del CSV: %d\n', protein_to_csvGene.Count);

%% 4. Generar _protein_clean_1.faa (cabecera = GeneID si existe mapeo, si no protein_id original)
fprintf('\n--- 4. Generando archivo FASTA limpio ---\n');
[headers, seqs] = readFasta(faaFile);
totalProteins = length(headers);
outputClean1 = [baseName '_protein_clean_1.faa'];
fid1 = fopen(outputClean1, 'w');
if fid1 == -1, error('No se pudo crear %s', outputClean1); end

replaced = 0;
for i = 1:totalProteins
    prot_id = strtok(headers{i});   % primer token (protein_id)
    if protein_to_csvGene.isKey(prot_id)
        newHeader = protein_to_csvGene(prot_id);   % Gglean002178
        replaced = replaced + 1;
    else
        newHeader = prot_id;
    end
    fprintf(fid1, '>%s\n', newHeader);
    seq = seqs{i};
    for j = 1:70:length(seq)
        fprintf(fid1, '%s\n', seq(j:min(j+69, end)));
    end
end
fclose(fid1);
fprintf('Archivo guardado: %s\n', outputClean1);
fprintf('Proteínas totales en FASTA: %d\n', totalProteins);
fprintf('Cabeceras reemplazadas: %d (%.2f%%)\n', replaced, replaced/totalProteins*100);

%% 5. Generar _protein_clean_reduced.faa (solo proteínas de genes en CSV)
fprintf('\n--- 5. Generando archivo REDUCIDO ---\n');
outputReduced = [baseName '_protein_clean_reduced.faa'];
fidR = fopen(outputReduced, 'w');
if fidR == -1, error('No se pudo crear %s', outputReduced); end

kept = 0;
for i = 1:totalProteins
    prot_id = strtok(headers{i});
    if protein_to_csvGene.isKey(prot_id)
        newHeader = protein_to_csvGene(prot_id);
        fprintf(fidR, '>%s\n', newHeader);
        seq = seqs{i};
        for j = 1:70:length(seq)
            fprintf(fidR, '%s\n', seq(j:min(j+69, end)));
        end
        kept = kept + 1;
    end
end
fclose(fidR);
fprintf('Archivo reducido guardado: %s\n', outputReduced);
fprintf('Proteínas mantenidas: %d\n', kept);
fprintf('Proteínas eliminadas: %d\n', totalProteins - kept);
fprintf('Porcentaje mantenido: %.2f%%\n', kept/totalProteins*100);

%% 6. Estadísticas finales
fprintf('\n=== ESTADÍSTICAS FINALES ===\n');
fprintf('Total genes en CSV (todos): %d\n', length(targetGenes));
fprintf('Genes ignorados (empiezan con AB): %d\n', genesIgnorados);
fprintf('Genes considerados para búsqueda: %d\n', length(targetGenes) - genesIgnorados);
fprintf('Sufijos encontrados en GTF: %d\n', matchedSuffixes.Count);
fprintf('Proteínas en FASTA mapeadas: %d\n', kept);
fprintf('=== PROCESO COMPLETADO ===\n');

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