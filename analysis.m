%% 
clear all 
clc
addpath(genpath('C:\Users\MSI\Documents\GitHub\T-15')) % вспомогательные функции и все такое
addpath(genpath("C:\Users\MSI\Documents\GitHub\T-10")) % для мультифрактального спектра, тоже дополнительные функции
cd('C:\Users\MSI\Documents\GitHub\T-15')
%%
    load("loaded_data.mat")

if ~exist("s1","var")

    data = ReadData("T15MD_4452.mat",0.5,200/2.7);
    s1 = data(:,1);
    dt = 1e-3;
    t = (1:length(s1)) * dt;
    t = t';

    plot(t,s1,'k')

    
    tt = title('$Signal$');
    
    set(gca, 'FontSize', 16, 'LineWidth', 2)
    set(gcf, 'Color', 'white')
    
    
    labels = [tt,...
              xlabel('$Time,\  ms$', 'FontSize', 14), ...
              ylabel('$Ion\ saturation\ current,\ A$', 'FontSize', 14)];
    set(labels, 'Interpreter', 'latex');


end

if ~exist("loaded_data.mat","file")
    
    save loaded_data

end

% интервал где есть плазма, определяется визуально

idx = find(t > 0.26 & t < 1.7);

s1_cleaned = s1(idx);
t_cleaned = t(idx);

%% MODWT DETRENDING 

t = t_cleaned;
y = s1_cleaned(:);

wname = 'db6';
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
clear  wname level energy energy_ratio y w w_trend trend trend_levels


%% PDF

fp = y_detrended;

% Данные (нормированные)
x = (fp - mean(fp)) / std(fp);  

skew = skewness(x); % коэффициент ассиметрии
kurt = kurtosis(x); % коэффициент эксцесса 

% Оценка плотности с помощью ядра
[f_pdf, xi] = ksdensity(x,NumPoints=length(x)/2);


figure; hold on;


plot(xi, f_pdf, 'o', 'MarkerSize', 8, 'Color', 'k');

% Гауссиан
x_gauss = linspace(min(xi), max(xi), 500);
y_gauss = (1/sqrt(2*pi)) * exp(-0.5 * x_gauss.^2);
plot(x_gauss, y_gauss, 'k', 'LineWidth', 2);

xlabel('$Normalized\ signal\ values $', 'Interpreter','latex');
ylabel('$Probability\ density$', 'Interpreter','latex');




titlestr = ['$PDF\ of\ signal\ vs\ Gaussian(skew = ' num2str(skew,2) ', kurt = ' num2str(kurt,2) ')$'];
title(titlestr, 'Interpreter','latex');


set(gca,'FontSize',16,'LineWidth',2);
set(gcf,'Color','white');
legend({'$Experimental\ PDF$ ','$Standard\ Gaussian$'}, 'Interpreter','latex');
grid on
box on




%% Fourier Spectrum
% 
% Fs = 1/dt;                % частота дискретизации
% N = length(s1);           % число точек
% FP_fft = fft(s1);         % прямое преобразование Фурье
% FP_mag = abs(FP_fft)/N;   % амплитуда (нормируем)
% f = (0:N-1)*(Fs/N);       % частоты от 0 до Fs
% 
% % --- основной спектр ---
% figure;
% subplot(2,1,1)
% loglog(f(1:floor(N/2)), 2*FP_mag(1:floor(N/2)), 'k', 'LineWidth', 1.5); 
% set(gca,'XScale','log','YScale','log','FontSize',16,'LineWidth',2); 
% set(gcf,'Color','white');
% labels = [title('$Fourier\ spectrum\ of\ signal\ (log-log\ scale)$', 'FontSize', 14), ...
%           xlabel('$f, Hz$', 'FontSize', 14), ...
%           ylabel('$s(f), rel.units$', 'FontSize', 14)];
% set(labels, 'Interpreter', 'latex');
% %grid on;
% 
% % --- данные для половины спектра ---
% f_half = f(1:floor(N/2));
% s_half = 2*FP_mag(1:floor(N/2));
% 
% lf = log10(f_half);
% ls = log10(s_half);
% 
% % --- фигура с кусочными фитами ---
% 
% % --- интервалы для фита (задаёт пользователь) ---
% intervals = [1e-2 1e0;  
%              1e0 1e2;  
%              1e2 1e3];  
% 
% subplot(2,1,2)
% hold on;
% colors = lines(size(intervals,1));
% legtxt = cell(size(intervals,1),1);
% 
% for k = 1:size(intervals,1)
%     fmin = intervals(k,1);
%     fmax = intervals(k,2);
% 
%     idx = (f_half >= fmin) & (f_half <= fmax);
% 
%     % Линейная аппроксимация в log-log
%     p = polyfit(lf(idx), ls(idx), 1);
%     alpha = p(1); % наклон
% 
%     % Вычисляем аппроксимацию
%     f_fit = logspace(log10(fmin), log10(fmax), 200);
%     s_fit = 10.^(polyval(p, log10(f_fit)));
% 
%     % Рисуем log-log линию
%     loglog(f_fit, s_fit, 'Color', colors(k,:), 'LineWidth', 2);
% 
%     % Подпись для легенды
%     legtxt{k} = sprintf('$%.1f\\ \\leq f(Hz) \\leq %.1f$,  $\\alpha=%.2f$', fmin, fmax, abs(alpha));
% end
% 
% set(gca,'XScale','log','YScale','log','FontSize',16,'LineWidth',2); 
% set(gcf,'Color','white');
% labels = [title('$Piecewise\ linear\ fit\ of\ spectrum$', 'FontSize', 14), ...
%           xlabel('$f, Hz$', 'FontSize', 14), ...
%           ylabel('$s(f), rel.units$', 'FontSize', 14)];
% set(labels, 'Interpreter', 'latex');
% 
% legend(legtxt, 'Interpreter','latex','FontSize',12,'Location','best');
% %grid on;
% 
% % в конце кода
% subplot(2,1,1);
% ax1 = gca;
% subplot(2,1,2);
% box on;
% ax2 = gca;
% 
% linkaxes([ax1 ax2], 'x');  % синхронизируем ось X
% linkaxes([ax1 ax2], 'y');  % синхронизируем ось X


%% 
dim = 3;
[~,lag] = phaseSpaceReconstruction(s1_cleaned,[],dim);
Np = 300;
correlationDimension(s1_cleaned,lag,dim,"NumPoints",Np);
%%
[~,elag,edim] = phaseSpaceReconstruction(fp);

%%
phaseSpaceReconstruction(fp,elag,edim)

%%
Np = 300;
correlationDimension(s1_cleaned,lag,dim,"NumPoints",Np);