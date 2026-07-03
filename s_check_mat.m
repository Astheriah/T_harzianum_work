% =====================================================
% VERIFICAR Y LUEGO CORREGIR MODELOS
% =====================================================

folderPath = pwd;
archivos_mat = dir(fullfile(folderPath, '*.mat'));
ignoreList = {'blast_reduced_THM10.mat', 'blast_THM10.mat'};

modelos_a_corregir = {};

% PRIMERO: IDENTIFICAR MODELOS A CORREGIR
fprintf('=== IDENTIFICANDO MODELOS A CORREGIR ===\n\n');

for i = 1:length(archivos_mat)
    nombre_archivo = archivos_mat(i).name;
    ruta_completa = fullfile(folderPath, nombre_archivo);
    
    if ismember(nombre_archivo, ignoreList)
        continue;
    end
    
    try
        loaded = load(ruta_completa);
        nombre_variable = erase(nombre_archivo, '.mat');
        
        if isfield(loaded, nombre_variable)
            model = loaded.(nombre_variable);
        else
            nombres_campos = fieldnames(loaded);
            model = loaded.(nombres_campos{1});
        end
        
        if ~isfield(model, 'metComps') || isempty(model.metComps)
            modelos_a_corregir{end+1} = nombre_archivo;
            fprintf('❌ %s - Falta metComps\n', nombre_archivo);
        else
            fprintf('✅ %s - OK\n', nombre_archivo);
        end
        
    catch ME
        fprintf('⚠️  %s - Error: %s\n', nombre_archivo, ME.message);
    end
end

% SEGUNDO: PREGUNTAR ANTES DE CORREGIR
if ~isempty(modelos_a_corregir)
    fprintf('\n📋 Modelos que necesitan corrección (%d):\n', length(modelos_a_corregir));
    for i = 1:length(modelos_a_corregir)
        fprintf('  %d. %s\n', i, modelos_a_corregir{i});
    end
    
    respuesta = input('\n¿Corregir estos modelos? (s/n): ', 's');
    
    if lower(respuesta) == 's'
        fprintf('\n=== APLICANDO CORRECCIONES ===\n');
        
        for i = 1:length(modelos_a_corregir)
            nombre_archivo = modelos_a_corregir{i};
            ruta_completa = fullfile(folderPath, nombre_archivo);
            
            fprintf('🛠️  Corrigiendo: %s... ', nombre_archivo);
            
            try
                % Cargar el archivo
                loaded = load(ruta_completa);
                
                % Determinar el nombre de la variable del modelo
                if isfield(loaded, erase(nombre_archivo, '.mat'))
                    nombre_variable = erase(nombre_archivo, '.mat');
                else
                    nombres_campos = fieldnames(loaded);
                    nombre_variable = nombres_campos{1};
                end
                
                % Obtener el modelo
                model = loaded.(nombre_variable);
                
                % Detectar formato de compartimento
                if ~isempty(model.mets)
                    primer_met = model.mets{1};
                    if contains(primer_met, '[')
                        compartment_char = '[';
                    else
                        compartment_char = '_';
                    end
                    
                    % Aplicar corrección
                    model_corregido = createMetComp(model, compartment_char);
                    
                    % SOBRESCRIBIR el modelo original manteniendo el mismo nombre de variable
                    eval([nombre_variable ' = model_corregido;']);
                    
                    % Guardar usando el nombre de variable original
                    save(ruta_completa, nombre_variable);
                    
                    fprintf('✅ (formato: %s)\n', compartment_char);
                else
                    fprintf('⚠️  no tiene mets\n');
                end
                
            catch ME
                fprintf('❌ %s\n', ME.message);
            end
        end
        
        fprintf('\n=== CORRECCIÓN COMPLETADA ===\n');
    else
        fprintf('\n=== CORRECCIÓN CANCELADA ===\n');
    end
else
    fprintf('\n✅ Todos los modelos ya tienen metComps\n');
end

% =====================================================
% VERIFICAR Y LUEGO CORREGIR MODELOS
% =====================================================

folderPath = pwd;
archivos_mat = dir(fullfile(folderPath, '*.mat'));
ignoreList = {'blast_reduced_THM10.mat', 'blast_THM10.mat'};

modelos_a_corregir = {};

% PRIMERO: IDENTIFICAR MODELOS A CORREGIR
fprintf('=== IDENTIFICANDO MODELOS A CORREGIR ===\n\n');

for i = 1:length(archivos_mat)
    nombre_archivo = archivos_mat(i).name;
    ruta_completa = fullfile(folderPath, nombre_archivo);
    
    if ismember(nombre_archivo, ignoreList)
        continue;
    end
    
    try
        loaded = load(ruta_completa);
        nombre_variable = erase(nombre_archivo, '.mat');
        
        % Determinar qué variable contiene el modelo
        if isfield(loaded, nombre_variable)
            modelo = loaded.(nombre_variable);
        else
            nombres_campos = fieldnames(loaded);
            if ~isempty(nombres_campos)
                modelo = loaded.(nombres_campos{1});
            else
                fprintf('⚠️  %s - Archivo vacío\n', nombre_archivo);
                continue;
            end
        end
        
        % Verificar si es una estructura
        if isstruct(modelo)
            if ~isfield(modelo, 'metComps') || isempty(modelo.metComps)
                modelos_a_corregir{end+1} = nombre_archivo;
                fprintf('❌ %s - Falta metComps (estructura)\n', nombre_archivo);
            else
                fprintf('✅ %s - OK (estructura)\n', nombre_archivo);
            end
        else
            fprintf('⚠️  %s - No es estructura (es %s)\n', nombre_archivo, class(modelo));
        end
        
    catch ME
        fprintf('⚠️  %s - Error: %s\n', nombre_archivo, ME.message);
    end
end

% SEGUNDO: PREGUNTAR ANTES DE CORREGIR
if ~isempty(modelos_a_corregir)
    fprintf('\n📋 Modelos que necesitan corrección (%d):\n', length(modelos_a_corregir));
    for i = 1:length(modelos_a_corregir)
        fprintf('  %d. %s\n', i, modelos_a_corregir{i});
    end
    
    respuesta = input('\n¿Corregir estos modelos? (s/n): ', 's');
    
    if lower(respuesta) == 's'
        fprintf('\n=== APLICANDO CORRECCIONES ===\n');
        
        for i = 1:length(modelos_a_corregir)
            nombre_archivo = modelos_a_corregir{i};
            ruta_completa = fullfile(folderPath, nombre_archivo);
            
            fprintf('🛠️  Corrigiendo: %s... ', nombre_archivo);
            
            try
                % Cargar el archivo
                loaded = load(ruta_completa);
                
                % Determinar el nombre de la variable
                if isfield(loaded, erase(nombre_archivo, '.mat'))
                    nombre_variable = erase(nombre_archivo, '.mat');
                else
                    nombres_campos = fieldnames(loaded);
                    nombre_variable = nombres_campos{1};
                end
                
                % Obtener el modelo
                modelo = loaded.(nombre_variable);
                
                % VERIFICAR TIPO DE DATO
                if ~isstruct(modelo)
                    fprintf('❌ No es estructura (es %s)\n', class(modelo));
                    continue;
                end
                
                % Verificar que tenga mets
                if ~isfield(modelo, 'mets') || isempty(modelo.mets)
                    fprintf('❌ No tiene campo mets\n');
                    continue;
                end
                
                % Detectar formato de compartimento
                primer_met = modelo.mets{1};
                if contains(primer_met, '[')
                    compartment_char = '[';
                    fprintf('(formato []) ');
                elseif contains(primer_met, '_')
                    compartment_char = '_';
                    fprintf('(formato _) ');
                else
                    compartment_char = '_';
                    fprintf('(formato default _) ');
                end
                
                % APLICAR CORRECCIÓN
                modelo_corregido = createMetComp(modelo, compartment_char);
                
                % SOBRESCRIBIR la variable original
                eval([nombre_variable ' = modelo_corregido;']);
                
                % Guardar
                save(ruta_completa, nombre_variable);
                
                fprintf('✅ Corregido\n');
                
            catch ME
                fprintf('❌ %s\n', ME.message);
            end
        end
        
        fprintf('\n=== CORRECCIÓN COMPLETADA ===\n');
    else
        fprintf('\n=== CORRECCIÓN CANCELADA ===\n');
    end
else
    fprintf('\n✅ Todos los modelos ya tienen metComps\n');
end

%% ===== CORREGIR METCOMPS Y DESCRIPTION (y ID) =====
folderPath = pwd;
archivos = dir(fullfile(folderPath, '*.mat'));
ignoreList = {'blast_reduced_THM10.mat', 'blast_THM10.mat'};
archivos = archivos(~ismember({archivos.name}, ignoreList));

for i = 1:length(archivos)
    nombre = archivos(i).name;
    ruta = fullfile(folderPath, nombre);
    baseName = erase(nombre, '.mat');   % nombre del modelo (sin extensión)
    fprintf('Procesando: %s... ', nombre);
    
    try
        loaded = load(ruta);
        % Identificar la variable que contiene el modelo (estructura con genes/rxns)
        campos = fieldnames(loaded);
        idx = find(structfun(@isstruct, loaded), 1);
        if isempty(idx)
            error('No se encontró ninguna estructura en %s', nombre);
        end
        var = campos{idx};
        modelo = loaded.(var);
        
        if ~isstruct(modelo)
            error('La variable %s no es una estructura', var);
        end
        
        % ---------- Corregir metComps ----------
        if ~isfield(modelo, 'metComps') || isempty(modelo.metComps)
            if isfield(modelo, 'mets') && ~isempty(modelo.mets)
                if contains(modelo.mets{1}, '[')
                    formato = '[';
                else
                    formato = '_';
                end
                modelo = createMetComp(modelo, formato);
                fprintf('metComps corregido ');
            else
                fprintf('sin mets ');
            end
        else
            fprintf('metComps OK ');
        end
        
        % ---------- Forzar description al nombre base ----------
        original_desc = '';
        if isfield(modelo, 'description')
            original_desc = modelo.description;
        end
        modelo.description = baseName;
        
        % ---------- Forzar id al nombre base ----------
        original_id = '';
        if isfield(modelo, 'id')
            original_id = modelo.id;
        end
        modelo.id = baseName;   % siempre se asigna (o sobreescribe)
        % Si también existe el campo 'modelID', actualizarlo por coherencia
        if isfield(modelo, 'modelID')
            modelo.modelID = baseName;
        end
        
        % Mostrar cambios en description e id
        fprintf('description: %s → %s, id: %s → %s ', ...
            original_desc, baseName, original_id, baseName);
        
        % Guardar el modelo actualizado
        eval([var ' = modelo;']);
        save(ruta, var);
        fprintf('✅\n');
        fprintf('  %s: description "%s" → "%s", id "%s" → "%s"\n', ...
            nombre, original_desc, baseName, original_id, baseName);
        
    catch ME
        fprintf('❌ error: %s\n', ME.message);
    end
end