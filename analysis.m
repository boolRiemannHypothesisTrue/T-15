%% 
addpath(genpath('C:\Users\MSI\Desktop\Курчатовский институт\Т-10 Анализ рядов')) % вспомогательные функции и все такое
%%

data = ReadData("T15MD_4452.mat",0.5,200/2.7);
s1 = data(:,1);
dt = 1e-3;
t = (1:length(s1)) * dt;
t = t';
plot(t,s1)

%%

idx = find(t > 0.26 & t < 1.7);

s1_cleaned = s1(idx);
t_cleaned = t(idx);

%% MODWT DETRENDING 

t = t_cleaned;
y = s1_cleaned(:);

wname = 'db8';
level = floor(log2(length(y))) - 1;
w = modwt(y, wname, level);

% Анализ энергии по уровням
energy = sum(w.^2, 2);
energy_ratio = energy / sum(energy);

% figure;
% bar(1:level+1, energy_ratio,'k');
% xlabel('Уровень MODWT', 'FontSize', 14);
% ylabel('Доля энергии', 'FontSize', 14);
% title('Энергия по уровням', 'FontSize', 14);
% set(gca, 'FontSize', 16, 'LineWidth', 2)
% set(gcf, 'Color', 'white')

% Адаптивный выбор тренд-уровней (берём где > 2 % энергии) . Можно
% поиграться с цифрой(больше/меньше)
trend_levels = find(energy_ratio > 0.05);

% Формируем тренд
w_trend = zeros(size(w));
w_trend(trend_levels, :) = w(trend_levels, :);

trend = imodwt(w_trend, wname);
trend = trend(:);
y_detrended = y - trend;

% Визуализация
figure;
hold on
plot(t, y, 'b');
plot(t, trend, 'g','LineWidth',2);
plot(t, y_detrended, 'k');
grid on;

hh = legend('$Signal$', '$Trend$', '$Detrended~signal$');
tt = title('$MODWT~Detrending$');
set(hh, 'Interpreter','latex')
set(gca, 'FontSize', 16, 'LineWidth', 2)
set(gcf, 'Color', 'white')


labels = [tt,hh,...
          xlabel('$Time,\ \mu s$', 'FontSize', 14), ...
          ylabel('$Ion\ saturation\ current,\ A$', 'FontSize', 14)];
set(labels, 'Interpreter', 'latex');
clear t wname level energy energy_ratio y w w_trend trend trend_levels


%%  Histogram

fp = s1_cleaned;
figure
bins = 100; % increase for better skewness/kurtosis visualization

histogram(fp, bins); % вместо histfit используем histogram

ylabel('$Number\ of\ observations$', Interpreter='latex')
xlabel('$Floating\ potential$', Interpreter='latex')
legend('Floating potential')
title("$ Floating\ potential\ histogram $", Interpreter="latex")
set(gca, 'FontSize', 16, 'LineWidth', 2) 
set(gcf, 'Color', 'white')


%% Skewness/ Kurtosis

skew = skewness(fp); % коэффициент ассиметрии
kurt = kurtosis(fp); % коэффициент эксцесса 

%% Fourier spectrum


warning off
s = fp';
Ts = 1e-6;
Fs = 1 / Ts;
n = length(s);

% Расчёт спектра
y = fft(s);
fshift = (-n/2:n/2-1)*(Fs/n);
yshift = fftshift(y);
z = abs(yshift);

% Положительные частоты
positive_idx = fshift > 0;
f_pos = fshift(positive_idx);

%  f^{-1}
x_ref = linspace(1e5, 1e7, 200);   % частотный диапазон
C = 1e8;                           % масштаб
y_ref = C ./ x_ref;                % f^{-1}

% График
figure
loglog(fshift, z,'k')
hold on
loglog(x_ref, y_ref, 'k', 'LineWidth', 1.5)
labels = [title('$Fourier\ spectrum\ of\ signal$', 'FontSize', 14), ...
          xlabel('$f, Hz$', 'FontSize', 14), ...
          ylabel('$s(f), rel.units$', 'FontSize', 14)];
          

% Подпись над линией f^{-1}
xtxt = x_ref(round(end/6));           % середина линии по X
ytxt = y_ref(round(end/2)) *15;     % чуть выше линии по Y
text(xtxt, ytxt, '$s(f) \sim f^{-1}$', ...
     'FontSize', 14, 'Color', 'k', 'Interpreter', 'latex')
set(labels, 'Interpreter', 'latex');
set(gca, 'FontSize', 16, 'LineWidth', 2);
set(gcf, 'Color', 'white');

savefig(fullfile("Fourier_spectrum_decline"))

%% Mask for wavelet transform 


% % Выбираем для анализа интервал по времени 
% 
% time_min = 300; % dt не более 5
% time_max = 305;
% % Создаем маску для индексов
% mask = (time >= time_min) & (time <= time_max);
% % Извлекаем данные
% time_for_wavelet = time(mask);
% fp_for_wavelet= fp(mask);




% %%%%%%%%%%%%%%%%%%%%%%%%%% UNCOMMENT FOR ABS. LVL WAVELET TRANSFORMS
% %% Wavelet MHAT. ABS.LVL
% %  Scale 30
% % 
% figure
% 
% tic
% subplot(2,1,1)
%  cwt(fp_for_wavelet,1:0.01:30, 'mexh','abslvl',[100 1000]);
% title("MHAT, scale 30")
% toc
% 
% %  Scale 100
% subplot(2,1,2)
% tic
% 
%  cwt(fp_for_wavelet,1:0.1:100, 'mexh','abslvl',[100 1000]);
% title("MHAT, scale 100")
% toc
% 
% 
% %% Wavelet Morlet ABS.LVL
% figure
% tic
% subplot(2,1,1)
%  cwt(fp_for_wavelet,1:0.01:30, 'morl','abslvl',[100 1000]);
% 
% toc
% % Wavelet Morlet. Scale 100
% 
% tic
% subplot(2,1,2)
%  cwt(fp_for_wavelet,1:0.1:100, 'morl','abslvl',[100 1000]);
%  toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% END UNCOMMENT

%% Wavelet transform -> Skeleton

figure 
fp_for_wavelet = s1_cleaned;
t = t_cleaned; % временная ось

% Вейвлет и шкала
wname = 'morl'; % например Морле
scales2 = 1:0.01:100;
cwt_coeffs2 = cwt(fp_for_wavelet, scales2, wname);

% Берем модуль
cfs_abs = abs(cwt_coeffs2);

% Матрица скелета
skeleton = zeros(size(cfs_abs));

for k = 1:size(cfs_abs,2)   % по времени
    % Максимумы
    [~,locsMax] = findpeaks(cfs_abs(:,k));
    skeleton(locsMax,k) = 1; 
    
    % Минимумы
    [~,locsMin] = findpeaks(-cfs_abs(:,k));
    skeleton(locsMin,k) = -1; 
end

% Отображение скелетограммы
imagesc(t, scales2, skeleton);
axis xy;
colormap([0 0 1; 1 1 1; 1 0 0]); % синий=минимумы, белый=фон, красный=максимумы
set(gca, 'YScale', 'log', 'YTick', 2.^(0:10));
ylabel('Scale, μs');
xlabel('Time, μs ');
title(['Skeleton of Wavelet (' wname ') coefficients']);

set(gca, 'FontSize', 16, 'LineWidth', 2)
set(gcf, 'Color', 'white')

%savefig(fullfile(results_path,['Wavelet ' num2str(wname)]))

%% 
dim = 3;
[~,lag] = phaseSpaceReconstruction(s1_cleaned,[],dim);
Np = 300;
correlationDimension(s1_cleaned,lag,dim,"NumPoints",Np);
%%
ts = s1_cleaned;
[~,elag,edim] = phaseSpaceReconstruction(ts)

%%
