initCobraToolbox()

% Procesar todos los archivos .xml y sbml en el directorio actual para transformar en archivos .mat y poder procesarlos mas rapido
archivos = [dir('*.xml'); dir('*.sbml')];

for i = 1:length(archivos)
    try
        filename = archivos(i).name;
        
        % Cargar modelo
        modelo_cargado = readCbModel(filename);
        
        % Obtener nombre del modelo
        if isfield(modelo_cargado, 'id') && ~isempty(modelo_cargado.id)
            modelName = modelo_cargado.id;
        else
            [~, modelName, ~] = fileparts(filename);
        end
        
        % Limpiar nombre
        cleanModelName = regexprep(modelName, '[\\/*?:"<>|]', '');
        outputFilename = sprintf('%s.mat', cleanModelName);
        outputPath = outputFilename;  % Guarda en el directorio actual
        
        eval([cleanModelName ' = modelo_cargado;']);
        
        % Guardar con el nombre original del modelo
        save(outputPath, cleanModelName);
        
        fprintf('✓ %s -> %s.mat (variable: %s)\n', filename, cleanModelName, cleanModelName);
        
    catch ME
        fprintf('✗ Error en %s: %s\n', filename, ME.message);
    end
end