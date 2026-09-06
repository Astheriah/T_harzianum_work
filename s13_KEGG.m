draftModel_genes = readtable('genesModelTemp2.csv');
EC_All = readtable('THM10_EC_info.csv');

EC_All.Inside_Model = ismember(EC_All.ProteinID, draftModel_genes.ProteinID);

% Cargar base de datos BiGG rxns
opts_reac = detectImportOptions('bigg_models_reactions.txt', 'FileType', 'text');
opts_reac.Delimiter = '\t';
T_bigg_reac = readtable('bigg_models_reactions.txt', opts_reac);
T_bigg_reac.Properties.VariableNames = strtrim(T_bigg_reac.Properties.VariableNames);

% Cargar base de datos BiGG mets
opts_met = detectImportOptions('bigg_models_metabolites.txt', 'FileType', 'text');
opts_met.Delimiter = '\t';
T_bigg_met = readtable('bigg_models_metabolites.txt', opts_met);
T_bigg_met.Properties.VariableNames = strtrim(T_bigg_met.Properties.VariableNames);


%Convertir la columna database_links de T_bigg_reac a un arreglo de strings
dbLinks = string(T_bigg_reac.database_links);
dbLinks = dbLinks(~ismissing(dbLinks) & dbLinks ~= "");

for i = 1:height(EC_All)
    kegg = EC_All.KEGG(i);
    % Ignorar valores vacíos o "No hay archivo"
    if ~ismissing(kegg) && kegg ~= "" && kegg ~= "No hay archivo"
        % Verificar si el KEGG aparece como subcadena en alguna de las)
        if any(contains(dbLinks, kegg, 'IgnoreCase', true))
            EC_All.KEGG_in_BiGG(i) = true;
        end
    end
end

fprintf('KEGG presentes en BiGG reactions: %d\n', sum(EC_All.KEGG_in_BiGG));


colNames = T_bigg_reac.Properties.VariableNames;
idxBigg = find(contains(colNames, 'bigg_id', 'IgnoreCase', true), 1);
idxName = find(contains(colNames, 'name', 'IgnoreCase', true), 1);
% Buscar columna de reacción (puede ser 'reaction_string', 'reactionString', 'reaction', etc.)
idxRxn = find(contains(colNames, 'reaction', 'IgnoreCase', true) & ...
              ~contains(colNames, 'database', 'IgnoreCase', true), 1);

if isempty(idxBigg) || isempty(idxName) || isempty(idxRxn)
    error('No se encontraron las columnas necesarias en T_bigg_reac. Revisa los nombres.');
end
biggCol = colNames{idxBigg};
nameCol = colNames{idxName};
rxnCol  = colNames{idxRxn};

EC_All.bigg_id     = strings(height(EC_All), 1);
EC_All.name        = strings(height(EC_All), 1);
EC_All.reaction    = strings(height(EC_All), 1); 

% --- Recorrer solo las filas con KEGG_in_BiGG = true
idxTrue = find(EC_All.KEGG_in_BiGG);
for i = idxTrue'
    kegg = EC_All.KEGG(i);
    % Buscar en T_bigg_reac las filas cuyo database_links contenga este KEGG
    mask = contains(T_bigg_reac.database_links, kegg, 'IgnoreCase', true);
    if any(mask)
        % Extraer los datos
        ids  = T_bigg_reac.(biggCol)(mask);
        names = T_bigg_reac.(nameCol)(mask);
        rxns = T_bigg_reac.(rxnCol)(mask);
        % Concatenar toda la información en una sola cadena
        EC_All.bigg_id(i)     = strjoin(ids, '; ');
        EC_All.name(i)        = strjoin(names, '; ');
        EC_All.reaction(i)    = strjoin(rxns, '; ');
    else
        EC_All.bigg_id(i)     = "";
        EC_All.name(i)        = "";
        EC_All.reaction(i)    = "";
    end
end

idx = ~EC_All.Inside_Model & EC_All.KEGG_in_BiGG;
KEGG_not_in_model = EC_All(idx, :);

fprintf('Filas en KEGG_not_in_model (Inside_Model=false y KEGG_in_BiGG=true): %d\n', height(KEGG_not_in_model));

idx = EC_All.Inside_Model & EC_All.KEGG_in_BiGG;
KEGG_Inside_Model = EC_All(idx, :);

fprintf('Filas en KEGG_Inside_Model (Inside_Model=true y KEGG_in_BiGG=true): %d\n', height(KEGG_Inside_Model));