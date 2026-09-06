% =====================================================
% CORREGIR AUTOMÁTICAMENTE MODELOS (metComps, description, id)
% =====================================================
% Procesa TODOS los .mat de la carpeta (excepto los ignorados):
%   - Si falta metComps, lo corrige con createMetComp.
%   - Siempre asigna description, id y modelID = nombre del archivo.
%   - Guarda el modelo con el nombre del archivo como variable interna.
% =====================================================

folderPath = pwd;
archivos_mat = dir(fullfile(folderPath, '*.mat'));
ignoreList = {'blast_reduced_THM10.mat', 'blast_THM10.mat'};
archivos_mat = archivos_mat(~ismember({archivos_mat.name}, ignoreList));

fprintf('=== PROCESANDO MODELOS ===\n\n');

for i = 1:length(archivos_mat)
    nombre_archivo = archivos_mat(i).name;
    ruta_completa = fullfile(folderPath, nombre_archivo);
    baseName = erase(nombre_archivo, '.mat');
    
    fprintf('🔄 Procesando: %s... ', nombre_archivo);
    
    try
        % ---------- CARGAR ----------
        loaded = load(ruta_completa);
        
        % Obtener el modelo (por nombre o primera estructura)
        if isfield(loaded, baseName)
            modelo = loaded.(baseName);
        else
            campos = fieldnames(loaded);
            idx = find(structfun(@isstruct, loaded), 1);
            if isempty(idx)
                error('No se encontró estructura');
            end
            modelo = loaded.(campos{idx});
        end
        
        if ~isstruct(modelo)
            error('La variable no es una estructura');
        end
        
        % ---------- CORREGIR metComps SI ES NECESARIO ----------
        if ~isfield(modelo, 'mets') || isempty(modelo.mets)
            error('El modelo no tiene campo mets');
        end
        
        if ~isfield(modelo, 'metComps') || isempty(modelo.metComps)
            % Detectar formato de compartimento
            primer_met = modelo.mets{1};
            if contains(primer_met, '[')
                compartimento = '[';
            elseif contains(primer_met, '_')
                compartimento = '_';
            else
                compartimento = '_';
            end
            
            modelo = createMetComp(modelo, compartimento);
            fprintf('metComps corregido (formato %s) ', compartimento);
        else
            fprintf('metComps OK ');
        end
        
        % ---------- ACTUALIZAR description, id, modelID (SIEMPRE) ----------
        modelo.description = baseName;
        modelo.id = baseName;
        if isfield(modelo, 'modelID')
            modelo.modelID = baseName;
        end
        
        % ---------- GUARDAR CON EL NOMBRE DEL ARCHIVO COMO VARIABLE ----------
        eval([baseName ' = modelo;']);
        save(ruta_completa, baseName);
        
        fprintf('✅ guardado como variable "%s"\n', baseName);
        
    catch ME
        fprintf('❌ Error: %s\n', ME.message);
    end
end

fprintf('\n=== PROCESAMIENTO COMPLETADO ===\n');