initCobraToolbox()

%=============================== t_reesei ========================================%
% scerbiomasspseudoreaction = 4006
% treebiomasspseudoreaction = 4007

load("t_reesei.mat")
[x, y] = findMetsFromRxns(t_reesei, t_reesei.rxns(4007))

x_flat = x;
while iscell(x_flat) && any(cellfun(@iscell, x_flat))
    if isscalar(x_flat) && iscell(x_flat{1})
        x_flat = x_flat{1};
    else
        x_flat = cellfun(@(c) c{1}, x_flat, 'UniformOutput', false);
    end
end

y_flat = y;
while iscell(y_flat) && any(cellfun(@iscell, y_flat))
    if isscalar(y_flat) && iscell(y_flat{1})
        y_flat = y_flat{1};
    else
        y_flat = cellfun(@(c) c{1}, y_flat, 'UniformOutput', false);
    end
end
y_num = cell2mat(y_flat);

% --- Convertir a string para comparación segura ---
mets_str = string(t_reesei.mets);          
x_str    = string(x_flat);

% --- Buscar índices en el modelo ---
[esta_en_mets, idx] = ismember(x_str, mets_str);

% Índices válidos (donde el metabolito fue encontrado)
valid_idx = idx(esta_en_mets);

% --- Extraer información adicional del modelo ---
metKEGGID_encontrados = t_reesei.metKEGGID(valid_idx);   
metNames_encontrados  = t_reesei.metNames(valid_idx);    
coef_encontrados      = y_num(esta_en_mets);

% --- Crear tabla con toda la información ---
if ~isempty(valid_idx)
    T1 = table(valid_idx, ...
              x_str(esta_en_mets), ...
              metNames_encontrados, ...
              metKEGGID_encontrados, ...
              coef_encontrados, ...
              'VariableNames', {'Indice', 'ID_Metabolito', 'Nombre', 'KEGG_ID', 'Coeficiente'});
    disp('Guardado en T1');
else
    disp('No se encontraron metabolitos válidos en el modelo.');
end




%=============================== iIsor850 ========================================%
%Biomass_Isor= 1825

load('iIsor850.mat');
[x, y] = findMetsFromRxns(iIsor850, iIsor850.rxns(1825));

x_flat = x;
while iscell(x_flat) && any(cellfun(@iscell, x_flat))
    if isscalar(x_flat) && iscell(x_flat{1})
        x_flat = x_flat{1};                       % Sacamos la celda interna
    else
        x_flat = cellfun(@(c) c{1}, x_flat, 'UniformOutput', false);
    end
end

y_flat = y;
while iscell(y_flat) && any(cellfun(@iscell, y_flat))
    if isscalar(y_flat) && iscell(y_flat{1})
        y_flat = y_flat{1};
    else
        y_flat = cellfun(@(c) c{1}, y_flat, 'UniformOutput', false);
    end
end

y_num = cell2mat(y_flat);
% --- Convertir a string para comparación segura ---
mets_str = string(iIsor850.mets);
x_str    = string(x_flat);
% --- Buscar índices ---
[esta_en_mets, idx] = ismember(x_str, mets_str);

% Índices válidos (donde el metabolito fue encontrado)
valid_idx = idx(esta_en_mets);   % solo los positivos

% Filtrar solo los metabolitos encontrados
valid_idx = idx(esta_en_mets);
metKEGGID_encontrados = iIsor850.metKEGGID(valid_idx);
metNames_encontrados  = iIsor850.metNames(valid_idx);
coef_encontrados      = y_num(esta_en_mets);   % coeficientes correspondientes

% --- Crear tabla con toda la información ---
if ~isempty(valid_idx)
    T2 = table(valid_idx, ...
              x_str(esta_en_mets), ...
              metNames_encontrados, ...
              metKEGGID_encontrados, ...
              coef_encontrados, ...
              'VariableNames', {'Indice', 'ID_Metabolito', 'Nombre', 'KEGG_ID', 'Coeficiente'});
    disp('Guardado en T2');
else
    disp('No se encontraron metabolitos válidos en el modelo.');
end




%=============================== IUL909 ========================================%
% BIOMASS = 1270
% Ex_biomass = 1967
% BIOMASS8020 = 2110
% BIOMASS6040 = 2111
% BIOMASS4060 = 2112
% BIOMASS8020b = 2113
% BIOMASS86040 = 2114
% BIOMASS4060b = 2115
% BIOMASS_glyc = 2148
% BIOMASS_meoh = 2152
% growth = 1184

load('iUL909.mat');
[x, y] = findMetsFromRxns(iUL909, iUL909.rxns(2148));

x_flat = x;
while iscell(x_flat) && any(cellfun(@iscell, x_flat))
    if isscalar(x_flat) && iscell(x_flat{1})
        x_flat = x_flat{1};                       % Sacamos la celda interna
    else
        x_flat = cellfun(@(c) c{1}, x_flat, 'UniformOutput', false);
    end
end

y_flat = y;
while iscell(y_flat) && any(cellfun(@iscell, y_flat))
    if isscalar(y_flat) && iscell(y_flat{1})
        y_flat = y_flat{1};
    else
        y_flat = cellfun(@(c) c{1}, y_flat, 'UniformOutput', false);
    end
end

y_num = cell2mat(y_flat);
% --- Convertir a string para comparación segura ---
mets_str = string(iUL909.mets);
x_str    = string(x_flat);
% --- Buscar índices ---
[esta_en_mets, idx] = ismember(x_str, mets_str);

% Índices válidos (donde el metabolito fue encontrado)
valid_idx = idx(esta_en_mets);   % solo los positivos

% Filtrar solo los metabolitos encontrados
valid_idx = idx(esta_en_mets);
metKEGGID_encontrados = iUL909.metKEGGID(valid_idx);
metNames_encontrados  = iUL909.metNames(valid_idx);
coef_encontrados      = y_num(esta_en_mets);   % coeficientes correspondientes

% --- Crear tabla con toda la información ---
if ~isempty(valid_idx)
    T3 = table(valid_idx, ...
              x_str(esta_en_mets), ...
              metNames_encontrados, ...
              metKEGGID_encontrados, ...
              coef_encontrados, ...
              'VariableNames', {'Indice', 'ID_Metabolito', 'Nombre', 'KEGG_ID', 'Coeficiente'});
    disp('Guardado en T3');
else
    disp('No se encontraron metabolitos válidos en el modelo.');
end

%=============================== iANid1221 ========================================%
% BIOMASS_SC5_notrace = 1381

load('iANid1221.mat');
[x, y] = findMetsFromRxns(iANid1221, iANid1221.rxns(1381));

x_flat = x;
while iscell(x_flat) && any(cellfun(@iscell, x_flat))
    if isscalar(x_flat) && iscell(x_flat{1})
        x_flat = x_flat{1};                       % Sacamos la celda interna
    else
        x_flat = cellfun(@(c) c{1}, x_flat, 'UniformOutput', false);
    end
end

y_flat = y;
while iscell(y_flat) && any(cellfun(@iscell, y_flat))
    if isscalar(y_flat) && iscell(y_flat{1})
        y_flat = y_flat{1};
    else
        y_flat = cellfun(@(c) c{1}, y_flat, 'UniformOutput', false);
    end
end

y_num = cell2mat(y_flat);
% --- Convertir a string para comparación segura ---
mets_str = string(iANid1221.mets);
x_str    = string(x_flat);
% --- Buscar índices ---
[esta_en_mets, idx] = ismember(x_str, mets_str);

% Índices válidos (donde el metabolito fue encontrado)
valid_idx = idx(esta_en_mets);   % solo los positivos

% Filtrar solo los metabolitos encontrados
valid_idx = idx(esta_en_mets);
if isfield(iANid1221, 'metKEGGID')
    metKEGGID_encontrados = iANid1221.metKEGGID(valid_idx);
else
    metKEGGID_encontrados = repmat({''}, size(valid_idx));  % celda de vacíos
end
if isfield(iANid1221, 'metNames')
    metNames_encontrados = iANid1221.metNames(valid_idx);
else
    metNames_encontrados = repmat({''}, size(valid_idx));
end
coef_encontrados      = y_num(esta_en_mets);   % coeficientes correspondientes

% --- Crear tabla con toda la información ---
if ~isempty(valid_idx)
    T4 = table(valid_idx, ...
              x_str(esta_en_mets), ...
              metNames_encontrados, ...
              metKEGGID_encontrados, ...
              coef_encontrados, ...
              'VariableNames', {'Indice', 'ID_Metabolito', 'Nombre', 'KEGG_ID', 'Coeficiente'});
    disp('Guardado en T4');
else
    disp('No se encontraron metabolitos válidos en el modelo.');
end


%=============================== iANig1029 ========================================%
% BIOMASS_SC5_notrace = 1375

load('iANig1029.mat');
[x, y] = findMetsFromRxns(iANig1029, iANig1029.rxns(1375));

x_flat = x;
while iscell(x_flat) && any(cellfun(@iscell, x_flat))
    if isscalar(x_flat) && iscell(x_flat{1})
        x_flat = x_flat{1};                       % Sacamos la celda interna
    else
        x_flat = cellfun(@(c) c{1}, x_flat, 'UniformOutput', false);
    end
end

y_flat = y;
while iscell(y_flat) && any(cellfun(@iscell, y_flat))
    if isscalar(y_flat) && iscell(y_flat{1})
        y_flat = y_flat{1};
    else
        y_flat = cellfun(@(c) c{1}, y_flat, 'UniformOutput', false);
    end
end

y_num = cell2mat(y_flat);
% --- Convertir a string para comparación segura ---
mets_str = string(iANig1029.mets);
x_str    = string(x_flat);
% --- Buscar índices ---
[esta_en_mets, idx] = ismember(x_str, mets_str);

% Índices válidos (donde el metabolito fue encontrado)
valid_idx = idx(esta_en_mets);   % solo los positivos

% Filtrar solo los metabolitos encontrados
valid_idx = idx(esta_en_mets);
if isfield(iANig1029, 'metKEGGID')
    metKEGGID_encontrados = iANid1221.metKEGGID(valid_idx);
else
    metKEGGID_encontrados = repmat({''}, size(valid_idx));  % celda de vacíos
end
if isfield(iANig1029, 'metNames')
    metNames_encontrados = iANid1221.metNames(valid_idx);
else
    metNames_encontrados = repmat({''}, size(valid_idx));
end
coef_encontrados      = y_num(esta_en_mets);   % coeficientes correspondientes

% --- Crear tabla con toda la información ---
if ~isempty(valid_idx)
    T5 = table(valid_idx, ...
              x_str(esta_en_mets), ...
              metNames_encontrados, ...
              metKEGGID_encontrados, ...
              coef_encontrados, ...
              'VariableNames', {'Indice', 'ID_Metabolito', 'Nombre', 'KEGG_ID', 'Coeficiente'});
    disp('Guardado en T5');
else
    disp('No se encontraron metabolitos válidos en el modelo.');
end



%=============================== iYli21 ========================================%
% xBIOMASS = 1846
% biomass_C = 2278

load('iYli21.mat');
[x, y] = findMetsFromRxns(iYli21, iYli21.rxns(2278));

x_flat = x;
while iscell(x_flat) && any(cellfun(@iscell, x_flat))
    if isscalar(x_flat) && iscell(x_flat{1})
        x_flat = x_flat{1};                       % Sacamos la celda interna
    else
        x_flat = cellfun(@(c) c{1}, x_flat, 'UniformOutput', false);
    end
end

y_flat = y;
while iscell(y_flat) && any(cellfun(@iscell, y_flat))
    if isscalar(y_flat) && iscell(y_flat{1})
        y_flat = y_flat{1};
    else
        y_flat = cellfun(@(c) c{1}, y_flat, 'UniformOutput', false);
    end
end

y_num = cell2mat(y_flat);
% --- Convertir a string para comparación segura ---
mets_str = string(iYli21.mets);
x_str    = string(x_flat);
% --- Buscar índices ---
[esta_en_mets, idx] = ismember(x_str, mets_str);

% Índices válidos (donde el metabolito fue encontrado)
valid_idx = idx(esta_en_mets);   % solo los positivos

% Filtrar solo los metabolitos encontrados
valid_idx = idx(esta_en_mets);
if isfield(iYli21, 'metKEGGID')
    metKEGGID_encontrados = iYli21.metKEGGID(valid_idx);
else
    metKEGGID_encontrados = repmat({''}, size(valid_idx));  % celda de vacíos
end
if isfield(iYli21, 'metNames')
    metNames_encontrados = iYli21.metNames(valid_idx);
else
    metNames_encontrados = repmat({''}, size(valid_idx));
end
coef_encontrados      = y_num(esta_en_mets);   % coeficientes correspondientes

% --- Crear tabla con toda la información ---
if ~isempty(valid_idx)
    T7 = table(valid_idx, ...
              x_str(esta_en_mets), ...
              metNames_encontrados, ...
              metKEGGID_encontrados, ...
              coef_encontrados, ...
              'VariableNames', {'Indice', 'ID_Metabolito', 'Nombre', 'KEGG_ID', 'Coeficiente'});
    disp('Guardado en T7');
else
    disp('No se encontraron metabolitos válidos en el modelo.');
end


%=============================== iANig1029 ========================================%
% BIOMASS_SC5_notrace = 1377

load('iAP1008.mat');
[x, y] = findMetsFromRxns(iAP1008, iAP1008.rxns(1377));

x_flat = x;
while iscell(x_flat) && any(cellfun(@iscell, x_flat))
    if isscalar(x_flat) && iscell(x_flat{1})
        x_flat = x_flat{1};                       % Sacamos la celda interna
    else
        x_flat = cellfun(@(c) c{1}, x_flat, 'UniformOutput', false);
    end
end

y_flat = y;
while iscell(y_flat) && any(cellfun(@iscell, y_flat))
    if isscalar(y_flat) && iscell(y_flat{1})
        y_flat = y_flat{1};
    else
        y_flat = cellfun(@(c) c{1}, y_flat, 'UniformOutput', false);
    end
end

y_num = cell2mat(y_flat);
% --- Convertir a string para comparación segura ---
mets_str = string(iAP1008.mets);
x_str    = string(x_flat);
% --- Buscar índices ---
[esta_en_mets, idx] = ismember(x_str, mets_str);

% Índices válidos (donde el metabolito fue encontrado)
valid_idx = idx(esta_en_mets);   % solo los positivos

% Filtrar solo los metabolitos encontrados
valid_idx = idx(esta_en_mets);
if isfield(iAP1008, 'metKEGGID')
    metKEGGID_encontrados = iANid1221.metKEGGID(valid_idx);
else
    metKEGGID_encontrados = repmat({''}, size(valid_idx));  % celda de vacíos
end
if isfield(iAP1008, 'metNames')
    metNames_encontrados = iANid1221.metNames(valid_idx);
else
    metNames_encontrados = repmat({''}, size(valid_idx));
end
coef_encontrados      = y_num(esta_en_mets);   % coeficientes correspondientes

% --- Crear tabla con toda la información ---
if ~isempty(valid_idx)
    T6 = table(valid_idx, ...
              x_str(esta_en_mets), ...
              metNames_encontrados, ...
              metKEGGID_encontrados, ...
              coef_encontrados, ...
              'VariableNames', {'Indice', 'ID_Metabolito', 'Nombre', 'KEGG_ID', 'Coeficiente'});
    disp('Guardado en T6');
else
    disp('No se encontraron metabolitos válidos en el modelo.');
end


%=========================== T. reesei (T1) vs iIlsor850 (T2) ===================
kegg_T1 = string(T1.KEGG_ID);
kegg_T2 = string(T2.KEGG_ID);

% --- Encontrar los KEGG_ID que están en ambas (sin contar vacíos) ---
comunes = intersect(kegg_T1, kegg_T2);
comunes(comunes == "") = [];  

% --- Para cada KEGG común, localizar sus filas en T1 y T2 ---
% Asumimos que un KEGG aparece una sola vez en cada tabla (si no, se expande)
[~, idx_T1] = ismember(comunes, kegg_T1);
[~, idx_T2] = ismember(comunes, kegg_T2);

% --- Extraer las filas correspondientes de cada tabla ---
T1_filas = T1(idx_T1, :);
T2_filas = T2(idx_T2, :);

% --- Renombrar las columnas para evitar conflicto y construir la tabla combinada ---
% Primero renombramos todas las variables de T2 añadiendo '_T2'
vars_T2 = T2.Properties.VariableNames;
for i = 1:numel(vars_T2)
    T2_filas.Properties.VariableNames{vars_T2{i}} = [vars_T2{i} '_T2'];
end

% También renombramos las de T1 añadiendo '_T1' (excepto tal vez KEGG_ID que podría ser la clave)
vars_T1 = T1.Properties.VariableNames;
for i = 1:numel(vars_T1)
    T1_filas.Properties.VariableNames{vars_T1{i}} = [vars_T1{i} '_T1'];
end

% Unir horizontalmente
T1vsT2 = [T1_filas, T2_filas];

% Opcional: mover la columna KEGG_ID común al principio y eliminar duplicada
% (en este caso tenemos KEGG_ID_T1 y KEGG_ID_T2 con el mismo valor)
% Si prefieres una sola columna KEGG_ID, puedes hacer:
T1vsT2.KEGG_ID = comunes;  % agrega columna común
T1vsT2.KEGG_ID_T1 = [];     % eliminar columnas que no sirven
T1vsT2.KEGG_ID_T2 = [];
T1vsT2.Indice_T1 = [];
T1vsT2.Indice_T2 = [];

T1vsT2 = movevars(T1vsT2, 'KEGG_ID', 'Before', 1);

disp('Tabla T1vsT2 creada con las coincidencias de KEGG_ID:');

%=========================== T. reesei (T1) vs iUL909 (T3) ===================
kegg_T1 = string(T1.KEGG_ID);
kegg_T3 = string(T3.KEGG_ID);

% --- Encontrar los KEGG_ID que están en ambas (sin contar vacíos) ---
comunes = intersect(kegg_T1, kegg_T3);
comunes(comunes == "") = [];  

% --- Para cada KEGG común, localizar sus filas en T1 y T3 ---
[~, idx_T1] = ismember(comunes, kegg_T1);
[~, idx_T3] = ismember(comunes, kegg_T3);

% --- Extraer las filas correspondientes de cada tabla ---
T1_filas = T1(idx_T1, :);
T3_filas = T3(idx_T3, :);

% --- Renombrar las columnas para evitar conflicto y construir la tabla combinada ---
vars_T3 = T3.Properties.VariableNames;
for i = 1:numel(vars_T3)
    T3_filas.Properties.VariableNames{vars_T3{i}} = [vars_T3{i} '_T3'];
end

vars_T1 = T1.Properties.VariableNames;
for i = 1:numel(vars_T1)
    T1_filas.Properties.VariableNames{vars_T1{i}} = [vars_T1{i} '_T1'];
end

% Unir horizontalmente
T1vsT3 = [T1_filas, T3_filas];

% Agregar columna común KEGG_ID y eliminar las redundantes
T1vsT3.KEGG_ID = comunes;
T1vsT3.KEGG_ID_T1 = [];
T1vsT3.KEGG_ID_T3 = [];
T1vsT3.Indice_T1 = [];
T1vsT3.Indice_T3 = [];

T1vsT3 = movevars(T1vsT3, 'KEGG_ID', 'Before', 1);

disp('Tabla T1vsT3 creada con las coincidencias de KEGG_ID:');
disp(T1vsT3);

%=========================== T. reesei (T1) vs iANid1221 (T4) ===================
% Coincidencia parcial: Nombre corto de T4 aparece en Nombre de T1

% --- Obtener columnas Nombre como string ---
nombres_T1 = string(T1.Nombre);
nombres_T4 = string(T4.Nombre);

% --- Listas para almacenar los índices de las filas coincidentes ---
iT4_list = [];
iT1_list = [];

% --- Recorrer cada nombre corto de T4 ---
for i = 1:length(nombres_T4)
    nombre_corto = nombres_T4(i);
    if nombre_corto == ""
        continue;  % ignorar vacíos
    end
    % Buscar todas las filas de T1 cuyo Nombre contenga el nombre corto
    idx_T1 = find(contains(nombres_T1, nombre_corto, 'IgnoreCase', true));
    if ~isempty(idx_T1)
        n = length(idx_T1);
        iT4_list = [iT4_list; repmat(i, n, 1)];
        iT1_list = [iT1_list; idx_T1(:)];
    end
end

% --- Extraer las filas correspondientes de cada tabla ---
T1_filas = T1(iT1_list, :);
T4_filas = T4(iT4_list, :);

% --- Renombrar columnas con sufijos _T1 y _T4 ---
vars_T1 = T1.Properties.VariableNames;
for j = 1:numel(vars_T1)
    T1_filas.Properties.VariableNames{vars_T1{j}} = [vars_T1{j} '_T1'];
end

vars_T4 = T4.Properties.VariableNames;
for j = 1:numel(vars_T4)
    T4_filas.Properties.VariableNames{vars_T4{j}} = [vars_T4{j} '_T4'];
end

% --- Unir horizontalmente ---
T1vsT4 = [T1_filas, T4_filas];

% --- Agregar columna común con el nombre corto que provocó el match (Nombre_Corto_T4) ---
T1vsT4.Nombre_Corto_T4 = nombres_T4(iT4_list);

% --- Eliminar columnas Indice_T1 e Indice_T4 (si existen) ---
if any(strcmp(T1vsT4.Properties.VariableNames, 'Indice_T1'))
    T1vsT4.Indice_T1 = [];
end
if any(strcmp(T1vsT4.Properties.VariableNames, 'Indice_T4'))
    T1vsT4.Indice_T4 = [];
end

% --- Mover Nombre_Corto_T4 al principio para claridad ---
T1vsT4 = movevars(T1vsT4, 'Nombre_Corto_T4', 'Before', 1);

disp('Tabla T1vsT4 creada con coincidencias parciales de Nombre:');
disp(T1vsT4);

%=========================== T. reesei (T1) vs iANig1029 (T5) ===================
% Coincidencia parcial: Nombre corto de T5 aparece en Nombre de T1

% --- Obtener columnas Nombre como string ---
nombres_T1 = string(T1.Nombre);
nombres_T5 = string(T5.Nombre);

% --- Listas para almacenar los índices de las filas coincidentes ---
iT5_list = [];
iT1_list = [];

% --- Recorrer cada nombre corto de T5 ---
for i = 1:length(nombres_T5)
    nombre_corto = nombres_T5(i);
    if nombre_corto == ""
        continue;  % ignorar vacíos
    end
    % Buscar todas las filas de T1 cuyo Nombre contenga el nombre corto
    idx_T1 = find(contains(nombres_T1, nombre_corto, 'IgnoreCase', true));
    if ~isempty(idx_T1)
        n = length(idx_T1);
        iT5_list = [iT5_list; repmat(i, n, 1)];
        iT1_list = [iT1_list; idx_T1(:)];
    end
end

% --- Extraer las filas correspondientes de cada tabla ---
T1_filas = T1(iT1_list, :);
T5_filas = T5(iT5_list, :);

% --- Renombrar columnas con sufijos _T1 y _T5 ---
vars_T1 = T1.Properties.VariableNames;
for j = 1:numel(vars_T1)
    T1_filas.Properties.VariableNames{vars_T1{j}} = [vars_T1{j} '_T1'];
end

vars_T5 = T5.Properties.VariableNames;
for j = 1:numel(vars_T5)
    T5_filas.Properties.VariableNames{vars_T5{j}} = [vars_T5{j} '_T5'];
end

% --- Unir horizontalmente ---
T1vsT5 = [T1_filas, T5_filas];

% --- Agregar columna común con el nombre corto que provocó el match (Nombre_Corto_T5) ---
T1vsT5.Nombre_Corto_T5 = nombres_T5(iT5_list);

% --- Eliminar columnas Indice_T1 e Indice_T5 (si existen) ---
if any(strcmp(T1vsT5.Properties.VariableNames, 'Indice_T1'))
    T1vsT5.Indice_T1 = [];
end
if any(strcmp(T1vsT5.Properties.VariableNames, 'Indice_T5'))
    T1vsT5.Indice_T5 = [];
end

% --- Mover Nombre_Corto_T5 al principio para claridad ---
T1vsT5 = movevars(T1vsT5, 'Nombre_Corto_T5', 'Before', 1);

disp('Tabla T1vsT5 creada con coincidencias parciales de Nombre:');
disp(T1vsT5);

%=========================== T. reesei (T1) vs iAP1008 (T6) ===================
% Coincidencia parcial: Nombre corto de T6 aparece en Nombre de T1

% --- Obtener columnas Nombre como string ---
nombres_T1 = string(T1.Nombre);
nombres_T6 = string(T6.Nombre);

% --- Listas para almacenar los índices de las filas coincidentes ---
iT6_list = [];
iT1_list = [];

% --- Recorrer cada nombre corto de T6 ---
for i = 1:length(nombres_T6)
    nombre_corto = nombres_T6(i);
    if nombre_corto == ""
        continue;  % ignorar vacíos
    end
    % Buscar todas las filas de T1 cuyo Nombre contenga el nombre corto
    idx_T1 = find(contains(nombres_T1, nombre_corto, 'IgnoreCase', true));
    if ~isempty(idx_T1)
        n = length(idx_T1);
        iT6_list = [iT6_list; repmat(i, n, 1)];
        iT1_list = [iT1_list; idx_T1(:)];
    end
end

% --- Extraer las filas correspondientes de cada tabla ---
T1_filas = T1(iT1_list, :);
T6_filas = T6(iT6_list, :);

% --- Renombrar columnas con sufijos _T1 y _T6 ---
vars_T1 = T1.Properties.VariableNames;
for j = 1:numel(vars_T1)
    T1_filas.Properties.VariableNames{vars_T1{j}} = [vars_T1{j} '_T1'];
end

vars_T6 = T6.Properties.VariableNames;
for j = 1:numel(vars_T6)
    T6_filas.Properties.VariableNames{vars_T6{j}} = [vars_T6{j} '_T6'];
end

% --- Unir horizontalmente ---
T1vsT6 = [T1_filas, T6_filas];

% --- Agregar columna común con el nombre corto que provocó el match (Nombre_Corto_T6) ---
T1vsT6.Nombre_Corto_T6 = nombres_T6(iT6_list);

% --- Eliminar columnas Indice_T1 e Indice_T6 (si existen) ---
if any(strcmp(T1vsT6.Properties.VariableNames, 'Indice_T1'))
    T1vsT6.Indice_T1 = [];
end
if any(strcmp(T1vsT6.Properties.VariableNames, 'Indice_T6'))
    T1vsT6.Indice_T6 = [];
end

% --- Mover Nombre_Corto_T6 al principio para claridad ---
T1vsT6 = movevars(T1vsT6, 'Nombre_Corto_T6', 'Before', 1);

disp('Tabla T1vsT6 creada con coincidencias parciales de Nombre:');
disp(T1vsT6);



%=========================== T. reesei (T1) vs iYli21 (T7) ===================
% Coincidencia parcial: Nombre corto extraído de T7 (antes del primer '_')
% aparece en Nombre de T1

% --- Obtener columnas Nombre como string ---
nombres_T1 = string(T1.Nombre);
nombres_T7_original = string(T7.Nombre);

% --- Preprocesar nombres de T7: tomar solo la parte antes del primer '_' ---
nombres_T7_corto = strings(size(nombres_T7_original));
for k = 1:length(nombres_T7_original)
    str = nombres_T7_original(k);
    idx_us = strfind(str, '_');
    if ~isempty(idx_us)
        nombres_T7_corto(k) = extractBefore(str, idx_us(1));
    else
        nombres_T7_corto(k) = str;   % sin guion bajo, se queda igual
    end
end

% --- Listas para almacenar los índices de las filas coincidentes ---
iT7_list = [];
iT1_list = [];

% --- Recorrer cada nombre corto de T7 ---
for i = 1:length(nombres_T7_corto)
    nombre_corto = nombres_T7_corto(i);
    if nombre_corto == ""
        continue;  % ignorar vacíos
    end
    % Buscar todas las filas de T1 cuyo Nombre contenga el nombre corto
    idx_T1 = find(contains(nombres_T1, nombre_corto, 'IgnoreCase', true));
    if ~isempty(idx_T1)
        n = length(idx_T1);
        iT7_list = [iT7_list; repmat(i, n, 1)];
        iT1_list = [iT1_list; idx_T1(:)];
    end
end

% --- Extraer las filas correspondientes de cada tabla ---
T1_filas = T1(iT1_list, :);
T7_filas = T7(iT7_list, :);

% --- Renombrar columnas con sufijos _T1 y _T7 ---
vars_T1 = T1.Properties.VariableNames;
for j = 1:numel(vars_T1)
    T1_filas.Properties.VariableNames{vars_T1{j}} = [vars_T1{j} '_T1'];
end

vars_T7 = T7.Properties.VariableNames;
for j = 1:numel(vars_T7)
    T7_filas.Properties.VariableNames{vars_T7{j}} = [vars_T7{j} '_T7'];
end

% --- Unir horizontalmente ---
T1vsT7 = [T1_filas, T7_filas];

% --- Agregar columna con el nombre corto extraído que provocó el match ---
T1vsT7.Nombre_Corto_T7 = nombres_T7_corto(iT7_list);

% --- Eliminar columnas Indice_T1 e Indice_T7 (si existen) ---
if any(strcmp(T1vsT7.Properties.VariableNames, 'Indice_T1'))
    T1vsT7.Indice_T1 = [];
end
if any(strcmp(T1vsT7.Properties.VariableNames, 'Indice_T7'))
    T1vsT7.Indice_T7 = [];
end

% --- Mover Nombre_Corto_T7 al principio ---
T1vsT7 = movevars(T1vsT7, 'Nombre_Corto_T7', 'Before', 1);

disp('Tabla T1vsT7 creada con coincidencias parciales (nombre corto extraído)');




% ========================================================
% Filtrar tablas T1vsT4, T1vsT5, T1vsT6:
% Para cada Nombre_Corto_Tx, conservar solo la fila con
% la menor diferencia absoluta (Dif_Abs_Coef).
% Luego recalcular las métricas de similitud.
% ========================================================
if ~ismember('Dif_Abs_Coef', T1vsT4.Properties.VariableNames)
    T1vsT4.Dif_Abs_Coef = abs(T1vsT4.Coeficiente_T1 - T1vsT4.Coeficiente_T4);
end
if ~ismember('Dif_Abs_Coef', T1vsT5.Properties.VariableNames)
    T1vsT5.Dif_Abs_Coef = abs(T1vsT5.Coeficiente_T1 - T1vsT5.Coeficiente_T5);
end
if ~ismember('Dif_Abs_Coef', T1vsT6.Properties.VariableNames)
    T1vsT6.Dif_Abs_Coef = abs(T1vsT6.Coeficiente_T1 - T1vsT6.Coeficiente_T6);
end

% Tablas originales y nombres de columna clave
tablas = {T1vsT4, T1vsT5, T1vsT6};
sufijos = {'T4', 'T5', 'T6'};
col_nombre = @(s) ['Nombre_Corto_' s];   % ej: 'Nombre_Corto_T4'
col_dif    = 'Dif_Abs_Coef';

% Inicializar celdas para las nuevas tablas
T_filt = cell(1,3);
for idx = 1:3
    T = tablas{idx};
    colID  = col_nombre(sufijos{idx});
    
    % --- Agrupar por el identificador corto y seleccionar la menor diferencia ---
    ids = T.(colID);
    [ids_unicos, ~, grupos] = unique(ids, 'stable');
    filas_sel = false(height(T), 1);
    for g = 1:length(ids_unicos)
        filas_grupo = find(grupos == g);
        difs = T.(col_dif)(filas_grupo);
        [~, pos_min] = min(difs);
        filas_sel(filas_grupo(pos_min)) = true;
    end
    T_filt{idx} = T(filas_sel, :);
end

% Asignar a variables con nombre v2
T1vsT4_v2 = T_filt{1};
T1vsT5_v2 = T_filt{2};
T1vsT6_v2 = T_filt{3};

% --------------------------------------------------------
% Recalcular métricas para las tablas filtradas
% --------------------------------------------------------


tablas_v2 = {T1vsT4_v2, T1vsT5_v2, T1vsT6_v2};
modelos   = {'iANid1221 (T4) v2', 'iANig1029 (T5) v2', 'iAP1008 (T6) v2'};
resultados_v2 = cell(1,3);
tolerancia = 1e-12;

for k = 1:3
    T = tablas_v2{k};
    col_T1 = 'Coeficiente_T1';
    col_Tk = ['Coeficiente_' sufijos{k}];
    
    if ~all(ismember({col_T1, col_Tk}, T.Properties.VariableNames))
        fprintf('La tabla de %s no tiene las columnas necesarias.\n', modelos{k});
        continue;
    end
    
    coef_T1 = T.(col_T1);
    coef_Tk = T.(col_Tk);
    coef_T1 = coef_T1(:);
    coef_Tk = coef_Tk(:);
    
    validos = ~isnan(coef_T1) & ~isnan(coef_Tk);
    if sum(validos) < 3
        warning('Pocos datos válidos en %s', modelos{k});
        continue;
    end
    x = coef_T1(validos);
    y = coef_Tk(validos);
    
    [r_pearson,  p_pearson]  = corr(x, y, 'Type', 'Pearson');
    [r_spearman, p_spearman] = corr(x, y, 'Type', 'Spearman');
    
    diferencias = abs(x - y);
    MAE  = mean(diferencias);
    RMSE = sqrt(mean((x - y).^2));
    R2   = r_pearson^2;
    
    iguales = diferencias < tolerancia;
    num_iguales = sum(iguales);
    total = length(x);
    porcentaje_igual = 100 * num_iguales / total;
    
    fprintf('\n==================================\n');
    fprintf('Comparación con %s (tras filtrar duplicados)\n', modelos{k});
    fprintf('==================================\n');
    fprintf('Filas originales: %d → filas tras filtro: %d\n', height(tablas{k}), total);
    fprintf('Coeficientes exactamente iguales: %d de %d (%.1f%%)\n', ...
            num_iguales, total, porcentaje_igual);
    fprintf('Correlación de Pearson : %.4f (p = %.4g)\n', r_pearson, p_pearson);
    fprintf('Correlación de Spearman: %.4f (p = %.4g)\n', r_spearman, p_spearman);
    fprintf('R² (lineal)            : %.4f\n', R2);
    fprintf('MAE (error abs. medio) : %.4f\n', MAE);
    fprintf('RMSE                   : %.4f\n', RMSE);
    
    resultados_v2{k} = table({modelos{k}}, r_pearson, p_pearson, r_spearman, p_spearman, ...
                             R2, MAE, RMSE, num_iguales, total, porcentaje_igual, ...
                             'VariableNames', {'Modelo', 'Pearson_r', 'Pearson_p', ...
                             'Spearman_r', 'Spearman_p', 'R2', 'MAE', 'RMSE', ...
                             'Coef_Iguales', 'Total_Comparados', 'Porcentaje_Iguales'});
end

% Mostrar resumen final
if ~isempty(resultados_v2)
    res_validos_v2 = resultados_v2(~cellfun(@isempty, resultados_v2));
    Resumen_Coeficientes_v2 = vertcat(res_validos_v2{:});
    disp('Resumen de similitud de coeficientes (TABLAS FILTRADAS):');
    disp(Resumen_Coeficientes_v2);
end