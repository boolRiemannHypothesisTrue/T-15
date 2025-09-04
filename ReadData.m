function processed = ReadData(filename, divider1, divider2)
% READDATA Обработка многоканальных данных из MAT-файла
%   Поддерживает файлы с:
%   - одной матрицей (столбцы = каналы)
%   - несколькими векторными переменными (каждая = канал)

    % Универсальная загрузка данных
    data_struct = load(filename);
    vars = fieldnames(data_struct);
    
    % Автоматическое определение формата данных
    if numel(vars) == 1
        % Случай 1: Одна матрица со всеми каналами
        raw_data = data_struct.(vars{1});
    else
        % Случай 2: Множество векторных переменных (каналов)
        % Сортировка для сохранения порядка каналов
        [~, idx] = sort_nat(vars); % Требуется функция sort_nat
        raw_data = cell2mat(cellfun(@(v) data_struct.(v), vars(idx), 'UniformOutput', false));
    end
    
    %% Параметры обработки (без изменений) %%
    [num_points, num_channels] = size(raw_data);
    baseline_samples = min(500, num_points); % Защита от коротких сигналов
    post_corr_samples = 2000;        
    decimation_factor = 200;      
    adc_coeff = 0.00031245;
    
    %% Обработка каналов %%
    new_length = ceil(num_points / decimation_factor);
    total_gain = (1/divider1) * divider2;
    
    % Первичная коррекция нуля
    for ch = 1:num_channels
        last_samples = raw_data(end-baseline_samples+1:end, ch);
        raw_data(:, ch) = raw_data(:, ch) - median(last_samples);
    end
    
    % Децимация и масштабирование
    processed = zeros(new_length, num_channels);
    for ch = 1:num_channels
        smoothed = smoothdata(raw_data(:,ch), 'gaussian', 20);
        scaled = smoothed * adc_coeff * total_gain;
        processed(:,ch) = decimate(scaled, decimation_factor);
    end
    
    % Финальная коррекция нуля
    for ch = 1:num_channels
        last_processed = processed(end-post_corr_samples+1:end, ch);
        processed(:,ch) = processed(:,ch) - median(last_processed);
    end
end

% Локальная функция для натуральной сортировки (если нет в пути)
function [sorted, idx] = sort_nat(cellarray)
    % Натуральная сортировка строк с числами
    [~, ~, ~, matches] = regexp(cellarray, '(\d+)|(\D+)');
    splits = cellfun(@(x) [x{:}], matches, 'UniformOutput', false);
    splits = cellfun(@(x) cellfun(@(y) [y,0], x, 'UniformOutput', false), splits, 'UniformOutput', false);
    [~, idx] = sortrows(splits);
    sorted = cellarray(idx);
end