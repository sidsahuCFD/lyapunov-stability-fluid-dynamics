clear all
close all
clc

CLV2x = load("CLV2x.txt");
CLV2y = load("CLV2y.txt");
N = size(CLV2x,1)*2;
CLV2 = zeros(N,18000);
CLV2(1:2:N,:) = CLV2x;
CLV2(2:2:N,:) = CLV2y;
clear CLV2x CLV2y


% Y = load("../CLV_collection.txt");
% Y = reshape(U,[168506 size(U,1)/168506]);
writematrix(mean(CLV2,2),'CLV2Mean.txt');
CLV2 = CLV2(:,7000:14000);
CLV2 = CLV2 - mean(CLV2,2);
% writematrix(CLV2,'CLV2.txt');
Y = transpose(CLV2);
N = size(Y,2)
%% If post-processing

clear all
close all
clc

spod_energy = load("results_CLV2/nfft256_novlp128_nblks61/spod_energy.mat")
L = spod_energy.L;
f = spod_energy.f;


%% Memory-efficient SPOD version that stores and reloads FFT blocks from hard drive.
%   In this example, we use a function handle to provide individual
%   snapshots to SPOD(_). Additionally, we ask SPOD(_) to save the FFT
%   blocks on hard drive instead of keeping all data in memory. These
%   features enable the SPOD of very large data sets but require more
%   computing time and additional hard drive space. We reduce the
%   additional storage requirenment by saving only a few modes at selected
%   frequencies.

opts.savefft    = true;             % save FFT blocks insteasd of keeping them in memory
opts.deletefft  = false;            % keep FFT blocks in this example for demonstration purposes
opts.savedir    = 'results_CLV2';        % save results to 'results' folder in the current directory
opts.savefreqs  = 1:1:120          % save modes frequencies of indices [10 15 20]
opts.nt         = 2000;             % use 2000 snapshots using XFUN             
% opts.mean       = p_mean;           % provide a long-time mean
opts.nsave      = 6;                % save the 5 most energetic modes

%   weights 
W = load("W.txt");

%   Use function handle to GETJET(_) to provide snapshots to SPOD(_). The
%   function file getjet.m can be found in the 'utils' folder. You can use 
%   getjet.m as a template to interface your own data with SPOD(_). A 
%   default (Hamming) window of length 256 with 128 snaphots overlap is
%   used in this example.
dt = 0.3;
[L,P,f] = spod(Y,256,W,128,dt,opts);

%% Plot the SPOD spectrum and some modes as before.
%   Note that P is a function handle that loads the corresponding modes from
%   hard drive since we are in FFT saving mode (OPTS.savefft is true), and
%   that the spectrum is restricted to the 5 frequencies specified through
%   OPTS.savefreqs. 

set(groot, ...
    'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex');

colors = {[0, 0.4470, 0.7410], ...  % Blue (SPOD mode 1)
          [0.8500, 0.3250, 0.0980], ...  % Orange (SPOD mode 2)
          [0.9290, 0.6940, 0.1250], ...  % Yellow (SPOD mode 3)
          [0.4940, 0.1840, 0.5560], ...  % Purple (SPOD mode 4)
          [0.4660, 0.6740, 0.1880], ...  % Green (SPOD mode 5)
          [0.3010, 0.7450, 0.9330]};  % Cyan/Light Blue (SPOD mode 6)

% define the two bands and their colours
clr1  = [0.00, 0.30, 0.00];   % dark green
clr2  = [0.30, 0.00, 0.30];   % dark purple
clrSub = [1 0.4 0.4];   % darker purple for sub‐harmonic

band1 = [0.057, 0.068];
band2 = [0.155, 0.174];
bandSub = [0.073 0.089];

ff = figure 
for i = 1:6
    semilogy(f, L(:,i), '-d', 'LineWidth', 1.5, 'MarkerFaceColor', colors{i}, 'Color', colors{i});
    hold on;
end

xlabel('frequency'), ylabel('SPOD mode energy')
xlim([0 0.5])
ylim([1e-6 6e-4])
grid on
legend('SPOD mode 1','SPOD mode 2','SPOD mode 3','SPOD mode 4','SPOD mode 5','SPOD mode 6','FontSize',18,'Location','eastoutside')
xlabel('$St$','Interpreter','Latex','Fontsize',18);
ylabel("SPOD mode energy",'Interpreter','Latex','Fontsize',18)
yl = ylim;
% p1 = patch([band1, fliplr(band1)], [yl(1) yl(1) yl(2) yl(2)], clr1, ...
%            'FaceAlpha',0.3,'EdgeColor','none','DisplayName','Lower frequency band');
% p2 = patch([band2, fliplr(band2)], [yl(1) yl(1) yl(2) yl(2)], clr2, ...
%            'FaceAlpha',0.3,'EdgeColor','none','DisplayName','Upper frequency band');
p3 = patch([bandSub,fliplr(bandSub)],[yl(1) yl(1) yl(2) yl(2)], clrSub, ...
           'FaceAlpha',0.3, 'EdgeColor','none', 'DisplayName','Subharmonic band');
% uistack([p2 p3  p1], 'bottom')



        P   = L(:,1);
        col = [0, 0.4470, 0.7410];

    m      = f >= bandSub(1) & f <= bandSub(2);
    [A1,i1] = max(P(m)); fvals = f(m); St1 = fvals(i1);
    plot(St1, A1, 'd', ...
         'MarkerEdgeColor',[0, 0.4470, 0.7410], 'MarkerFaceColor',col, ...
         'HandleVisibility','off','MarkerSize',12)
    text(St1, A1, sprintf(' $St$=%.3f',St1), ...
         'VerticalAlignment','bottom','FontSize',12)

    m      = f >= band2(1) & f <= band2(2);
    [A2,i2] = max(P(m)); fvals = f(m); St2 = fvals(i2);
    % plot(St2, A2, 's', ...
    %      'MarkerEdgeColor',col, 'MarkerFaceColor',col, ...
    %      'HandleVisibility','off')
    % text(St2, A2, sprintf(' $St$=%.3f',St2), ...
    %      'VerticalAlignment','bottom','FontSize',12)

 
grid on
% set(gca,'FontSize',14,'FontName','Courier')
legend('Location','northeast')
set(gcf,'Position',[10 10 1200 900])
fs = 24;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', fs)

% peaks at 0.078
%%

L(:,1)
[~,indx] = max(L(:,1));
writematrix(L(:,1:6),'L.txt');
writematrix(f,'f.txt');

%% Generating an animation of the modes

%   Note how all wavepackets travel at approximately the same phase
%   speed c_ph. The reason is that their streamwise wavenumber k_x changes 
%   with frequency such that c_ph = omega/k_x is approximately constant.
nt      = 50;
T       = 1/f(2);              % period of the 10th frequency
time    = linspace(0,T,nt);     % animate over one period

file = matfile(['results_CLV2/nfft256_novlp128_nblks61/spod_f' num2str(7,'%.4i')])
for ti = 1:nt
            f13m1ani(ti,:) = real(squeeze(file.Psi(:,1)*exp(2i*pi*f(2)*time(ti))));
            f13m2ani(ti,:) = real(squeeze(file.Psi(:,2)*exp(2i*pi*f(2)*time(ti))));
end

file = matfile(['results/nfft256_novlp128_nblks61/spod_f' num2str(12,'%.4i')])
for ti = 1:nt
            f29m1ani(ti,:) = real(squeeze(file.Psi(:,1)*exp(2i*pi*f(2)*time(ti))));
            f29m2ani(ti,:) = real(squeeze(file.Psi(:,2)*exp(2i*pi*f(2)*time(ti))));
end

writematrix(transpose(f13m1ani),'f13m1ani.txt','Delimiter','tab');
writematrix(transpose(f13m2ani),'f13m2ani.txt','Delimiter','tab');
writematrix(transpose(f29m1ani),'f29m1ani.txt','Delimiter','tab');
writematrix(transpose(f29m2ani),'f29m2ani.txt','Delimiter','tab');

%% Static eigenvectors in the frequency space

file = matfile(['results_CLV2/nfft256_novlp128_nblks61/spod_f' num2str(7,'%.4i')])

set(groot, ...
    'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex');

fCLV2SPODmode10078 = real(squeeze(file.Psi(:,1)));
CLV2SPOD1y = fCLV2SPODmode10078(1:2:end);
% CLV2SPOD1y = sqrt(fCLV2SPODmode10078(1:2:end).^2 + fCLV2SPODmode10078(2:2:end).^2);

x = load('../CLV_visualisation/X.txt');
y = load('../CLV_visualisation/Y.txt');

fig = figure(4);
sgtitle('$St = 0.078$','FontSize',20)

subplot(2,2,1)

x_length = [-2.5 11.5]; y_length = [-4 4];
pixels = 200;
% Define a regular grid covering the range of your data:
xq = linspace(-4, max(x), pixels);  % Change '200' as needed
yq = linspace(-8, 8, pixels);
[Xgrid, Ygrid] = meshgrid(xq, yq);
   
F = scatteredInterpolant(x, y, double(CLV2SPOD1y) , 'linear');
Z = F(Xgrid, Ygrid); % Evaluate the interpolant on the grid:
contourf(xq,yq,Z,50,'Edgecolor','flat'); hold on;
pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
plot(pgon,'FaceColor',[0 0 0],'FaceAlpha',1,'EdgeColor','none'); hold on;
pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
plot(pgon,'FaceColor',[0 0 0],'FaceAlpha',1,'EdgeColor','none');
shading interp; 
axis equal; 
view(2); 
colorbar;
colorbar('eastoutside','TickLabelInterpreter','latex');
title('$u$ $\Phi_{1r}^{0.078}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);

shading interp; axis equal; view(2); c = colorbar('eastoutside');
% set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%0.3f')));
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
% colormap(flipud(gray(256)));  % This inverts the grayscale
addpath('../../../../MATLAB_codes/BrewerMap/')
n = 256;                     % number of colors
cmap = brewermap(n,'RdBu');  % get RdBu
cmap = flipud(cmap); 
m = max(abs(Z(:)));
caxis([-m m])
colormap(cmap)
colorbar('eastoutside');
colorbar('eastoutside','TickLabelInterpreter','latex');



subplot(2,2,2)

CLV2SPOD1y = fCLV2SPODmode10078(2:2:end);
F = scatteredInterpolant(x, y, double(CLV2SPOD1y) , 'linear');
Z = F(Xgrid, Ygrid); % Evaluate the interpolant on the grid:
contourf(xq,yq,Z,50,'Edgecolor','flat'); hold on;
pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
plot(pgon,'FaceColor',[0 0 0],'FaceAlpha',1,'EdgeColor','none'); hold on;
pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
plot(pgon,'FaceColor',[0 0 0],'FaceAlpha',1,'EdgeColor','none');
shading interp; 
axis equal; 
view(2); 
colorbar;
title('$v$ $\Phi_{1r}^{0.078}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
shading interp; axis equal; view(2); c = colorbar('eastoutside');
% set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%0.3f')));
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
% colormap(flipud(gray(256)));  % This inverts the grayscale
addpath('../../../../MATLAB_codes/BrewerMap/')
n = 256;                     % number of colors
cmap = brewermap(n,'RdBu');  % get RdBu
cmap = flipud(cmap); 
m = max(abs(Z(:)));
caxis([-m m])
colormap(cmap)
colorbar('eastoutside');
colorbar('eastoutside','TickLabelInterpreter','latex');


fCLV2SPODmode10078 = imag(squeeze(file.Psi(:,1)));
CLV2SPOD1y = fCLV2SPODmode10078(1:2:end);

subplot(2,2,3)

x_length = [-2.5 11.5]; y_length = [-4 4];
% Define a regular grid covering the range of your data:
xq = linspace(-4, max(x), pixels);  % Change '200' as needed
yq = linspace(-8, 8, pixels);
[Xgrid, Ygrid] = meshgrid(xq, yq);
   
F = scatteredInterpolant(x, y, double(CLV2SPOD1y) , 'linear');
Z = F(Xgrid, Ygrid); % Evaluate the interpolant on the grid:
contourf(xq,yq,Z,50,'Edgecolor','flat'); hold on;
pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
plot(pgon,'FaceColor',[0 0 0],'FaceAlpha',1,'EdgeColor','none'); hold on;
pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
plot(pgon,'FaceColor',[0 0 0],'FaceAlpha',1,'EdgeColor','none');
shading interp; 
axis equal; 
view(2); 
colorbar;
colorbar('eastoutside','TickLabelInterpreter','latex');
title('$u$ $\Phi_{1i}^{0.078}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);

shading interp; axis equal; view(2); c = colorbar('eastoutside');
% set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%0.3f')));
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
% colormap(flipud(gray(256)));  % This inverts the grayscale
addpath('../../../../MATLAB_codes/BrewerMap/')
n = 256;                     % number of colors
cmap = brewermap(n,'RdBu');  % get RdBu
cmap = flipud(cmap); 
m = max(abs(Z(:)));
caxis([-m m])
colormap(cmap)
colorbar('eastoutside');
colorbar('eastoutside','TickLabelInterpreter','latex');



subplot(2,2,4)

CLV2SPOD1y = fCLV2SPODmode10078(2:2:end);
F = scatteredInterpolant(x, y, double(CLV2SPOD1y) , 'linear');
Z = F(Xgrid, Ygrid); % Evaluate the interpolant on the grid:
contourf(xq,yq,Z,50,'Edgecolor','flat'); hold on;
pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
plot(pgon,'FaceColor',[0 0 0],'FaceAlpha',1,'EdgeColor','none'); hold on;
pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
plot(pgon,'FaceColor',[0 0 0],'FaceAlpha',1,'EdgeColor','none');
shading interp; 
axis equal; 
view(2); 
colorbar;
title('$v$ $\Phi_{1i}^{0.078}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
shading interp; axis equal; view(2); c = colorbar('eastoutside');
% set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%0.3f')));
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
% colormap(flipud(gray(256)));  % This inverts the grayscale
addpath('../../../../MATLAB_codes/BrewerMap/')
n = 256;                     % number of colors
cmap = brewermap(n,'RdBu');  % get RdBu
cmap = flipud(cmap); 
m = max(abs(Z(:)));
caxis([-m m])
colormap(cmap)
colorbar('eastoutside');
colorbar('eastoutside','TickLabelInterpreter','latex');

fs = 16;
set(gcf,'Position',[10 10 1400 600])
set(findall(gcf, '-property', 'FontSize'), 'FontSize', fs)

print(gcf, 'SPODcontours_CLV2_XY', '-dpng', '-r500');    % dpi PNG
exportgraphics(gcf, 'SPODcontours_CLV2_XY.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none');
%%

xq = linspace(-4, max(x), pixels);  % Change '200' as needed
yq = 0;
[Xgrid, Ygrid] = meshgrid(xq, yq);
   
F = scatteredInterpolant(x, y, double(CLV2SPOD1y) , 'linear');
Z = F(Xgrid, Ygrid); % Evaluate the interpolant on the grid:

plot(xq,abs(Z))