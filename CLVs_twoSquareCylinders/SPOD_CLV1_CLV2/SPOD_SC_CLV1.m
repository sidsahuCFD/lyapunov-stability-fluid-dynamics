clear all
close all
clc

CLV1x = load("CLV1x.txt");
CLV1y = load("CLV1y.txt");
% CLV1 = load("CLV1.txt");
N = size(CLV1x,1)*2;
CLV1 = zeros(N,18000);
CLV1(1:2:N,:) = CLV1x;
CLV1(2:2:N,:) = CLV1y;
clear CLV1x CLV1y


% Y = load("../CLV_collection.txt");
% Y = reshape(U,[168506 size(U,1)/168506]);
writematrix(mean(CLV1,2),'CLV1Mean.txt');
CLV1 = CLV1(:,7000:14000);
CLV1 = CLV1 - mean(CLV1,2);
% writematrix(CLV2,'CLV1.txt');
Y = transpose(CLV1);
N = size(Y,2)

%% If post-processing

clear all; close all; clc

% spod_energy = load("results_CLV1/nfft512_novlp256_nblks30/spod_energy.mat")
spod_energy = load("results_CLV1/nfft256_novlp128_nblks61/spod_energy.mat")
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
opts.savedir    = 'results';        % save results to 'results' folder in the current directory
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
[L,P,f,~,A] = spod(Y,256,W,128,dt,opts);

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
p1 = patch([band1, fliplr(band1)], [yl(1) yl(1) yl(2) yl(2)], clr1, ...
           'FaceAlpha',0.3,'EdgeColor','none','DisplayName','Lower frequency band');
p2 = patch([band2, fliplr(band2)], [yl(1) yl(1) yl(2) yl(2)], clr2, ...
           'FaceAlpha',0.3,'EdgeColor','none','DisplayName','Upper frequency band');
% p3 = patch([bandSub,fliplr(bandSub)],[yl(1) yl(1) yl(2) yl(2)], clrSub, ...
%            'FaceAlpha',0.3, 'EdgeColor','none', 'DisplayName','Subharmonic band');
uistack([p2  p1], 'bottom')



        P   = L(:,1);
        col = [0, 0.4470, 0.7410];

    m      = f >= 0.05 & f <= band1(2);
    [A1,i1] = max(P(m)); fvals = f(m); St1 = fvals(i1);
    plot(St1, A1, 'd', ...
         'MarkerEdgeColor',[0, 0.4470, 0.7410], 'MarkerFaceColor',col, ...
         'HandleVisibility','off','MarkerSize',12)
    text(St1, A1, sprintf(' $St$=%.3f',St1), ...
         'VerticalAlignment','bottom','FontSize',12)

    m      = f >= band2(1) & f <= band2(2);
    [A2,i2] = max(P(m)); fvals = f(m); St2 = fvals(i2);
    plot(St2, A2, 'd', ...
         'MarkerEdgeColor',col, 'MarkerFaceColor',col, ...
         'HandleVisibility','off','MarkerSize',12)
    text(St2, A2, sprintf(' $St$=%.3f',St2), ...
         'VerticalAlignment','bottom','FontSize',12)

 
grid on
% set(gca,'FontSize',14,'FontName','Courier')
legend('Location','northeast')
set(gcf,'Position',[10 10 1200 900])
fs = 24;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', fs)


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

file = matfile(['results_CLV1/nfft256_novlp128_nblks61/spod_f' num2str(6,'%.4i')])
for ti = 1:nt
            % f13m1ani = real(file.Psi(:,1));
            % f13m2ani = imag(file.Psi(:,1));
            f13m1ani(ti,:) = real(squeeze(file.Psi(:,1)*exp(2i*pi*f(6)*time(ti))));
            f13m2ani(ti,:) = real(squeeze(file.Psi(:,2)*exp(2i*pi*f(6)*time(ti))));
end

file = matfile(['results_CLV1/nfft256_novlp128_nblks61/spod_f' num2str(13,'%.4i')])
% for ti = 1:nt
            f29m1ani = real(file.Psi(:,1));
            f29m2ani = imag(file.Psi(:,1));
% end

writematrix(transpose(f13m1ani),'f13m1ani.txt','Delimiter','tab');
writematrix(transpose(f13m2ani),'f13m2ani.txt','Delimiter','tab');
writematrix(transpose(f29m1ani),'f29m1ani.txt','Delimiter','tab');
writematrix(transpose(f29m2ani),'f29m2ani.txt','Delimiter','tab');

%% Generating an animation of the modes to locate missing crescent structures

%   Note how all wavepackets travel at approximately the same phase
%   speed c_ph. The reason is that their streamwise wavenumber k_x changes 
%   with frequency such that c_ph = omega/k_x is approximately constant.
nt      = 50;
T       = 1/f(2);              % period of the 10th frequency
time    = linspace(0,T,nt);     % animate over one period

file = matfile(['results_CLV1/nfft256_novlp128_nblks61/spod_f' num2str(11,'%.4i')])
            f13m1ani = real(file.Psi(:,1));

file = matfile(['results_CLV1/nfft256_novlp128_nblks61/spod_f' num2str(12,'%.4i')])
            f13m2ani = real(file.Psi(:,1));

file = matfile(['results_CLV1/nfft256_novlp128_nblks61/spod_f' num2str(13,'%.4i')])
            f29m1ani = real(file.Psi(:,1));

file = matfile(['results_CLV1/nfft256_novlp128_nblks61/spod_f' num2str(14,'%.4i')])
            f29m2ani = real(file.Psi(:,1));

writematrix(transpose(f13m1ani),'f13m1ani.txt','Delimiter','tab');
writematrix(transpose(f13m2ani),'f13m2ani.txt','Delimiter','tab');
writematrix(transpose(f29m1ani),'f29m1ani.txt','Delimiter','tab');
writematrix(transpose(f29m2ani),'f29m2ani.txt','Delimiter','tab');

%% Animation

load("colorRdBu.mat")
x = load('X.txt');
y = load('Y.txt');

x_length = [-2.5 11.5]; y_length = [-4 4];
pixels = 200;
% Define a regular grid covering the range of your data:
xq = linspace(-4, max(x), pixels);  % Change '200' as needed
yq = linspace(-8, 8, pixels);
[Xgrid, Ygrid] = meshgrid(xq, yq);

for t_i= 1:50
   
F = scatteredInterpolant(x, y, double(f13m1ani(t_i,1:2:end)'), 'linear');
Z = F(Xgrid, Ygrid);

%--- plot filled contours
figure(1);
contourf(xq, yq, Z, 50, 'EdgeColor', 'none'); 
hold on;

%--- overlay black rectangles
plot( polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]), 'FaceColor','k','FaceAlpha',1,'EdgeColor','none' );
plot( polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]), 'FaceColor','k','FaceAlpha',1,'EdgeColor','none' );

%--- axis, view, labels & ticks
axis equal off;
view(2);
xlim([xq(1) xq(end)]);
ylim([yq(1) yq(end)]);
xticks(0:4:40);
yticks(-30:4:30);
xlabel('x');
ylabel('y');
title('Instantaneous');

%--- smooth 512-step RdBu colormap
m = max(abs(Z(:)));
caxis([-m m]);
map512 = interp1(1:11, RdBu11, linspace(1,11,511), 'linear');
colormap(map512);
colorbar('eastoutside');
pause(0.1)
%--- save frame
saveas(gcf, sprintf('animation_CLV1x_SPODmode1_0013/a1%03d.png', t_i));

end
%% Static eigenvectors in the frequency space (Real part)

pixels = 200;
file = matfile(['results_CLV1/nfft256_novlp128_nblks61/spod_f' num2str(5,'%.4i')])

set(groot, ...
    'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex');

fCLV1SPODmodeLF = real(squeeze(file.Psi(:,1)));
% fCLV1SPODmodeLF = imag(squeeze(file.Psi(:,1)));
CLV1SPOD1yLF = fCLV1SPODmodeLF(1:2:end);
% CLV2SPOD1y = sqrt(fCLV2SPODmode10078(1:2:end).^2 + fCLV2SPODmode10078(2:2:end).^2);

x = load('../CLV_visualisation/X.txt');
y = load('../CLV_visualisation/Y.txt');

fig = figure(4);

ax(1) = subplot(2,2,1);

x_length = [-2.5 11.5]; y_length = [-4 4];
% Define a regular grid covering the range of your data:
xq = linspace(-4, max(x), pixels);  % Change '200' as needed
yq = linspace(-8, 8, pixels);
[Xgrid, Ygrid] = meshgrid(xq, yq);
F = scatteredInterpolant(x, y, double(CLV1SPOD1yLF) , 'linear');
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
title('$u$ $\Phi_{1r}^{0.052}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
shading interp; axis equal; view(2); c = colorbar('eastoutside');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
addpath('/media/sidsahu/Elements/PhD_files/MATLAB_codes/BrewerMap/')
n = 256;                     % number of colors
cmap = brewermap(n,'RdBu');  % get RdBu
cmap = flipud(cmap); 
m = max(abs(Z(:)));
caxis([-m m])
colormap(cmap)
colorbar('eastoutside');
colorbar('eastoutside','TickLabelInterpreter','latex');


ax(2) = subplot(2,2,2)

CLV1SPOD1yLF = fCLV1SPODmodeLF(2:2:end);
F = scatteredInterpolant(x, y, double(CLV1SPOD1yLF) , 'linear');
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
title('$v$ $\Phi_{1r}^{0.052}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
shading interp; axis equal; view(2); c = colorbar('eastoutside');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
addpath('../../../../MATLAB_codes/BrewerMap/')
n = 256;                     % number of colors
cmap = brewermap(n,'RdBu');  % get RdBu
cmap = flipud(cmap); 
m = max(abs(Z(:)));
caxis([-m m])
colormap(cmap)
colorbar('eastoutside');



file = matfile(['results_CLV1/nfft256_novlp128_nblks61/spod_f' num2str(13,'%.4i')])
fCLV1SPODmodeHF = real(squeeze(file.Psi(:,1)));
% fCLV1SPODmodeHF = imag(squeeze(file.Psi(:,1)));
CLV1SPOD1yHF = fCLV1SPODmodeHF(1:2:end);
colorbar('eastoutside','TickLabelInterpreter','latex');

ax(3) = subplot(2,2,3)

CLV1SPOD1yHF = fCLV1SPODmodeHF(1:2:end);
F = scatteredInterpolant(x, y, double(CLV1SPOD1yHF) , 'linear');
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
title('$u$ $\Phi_{1r}^{0.156}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
shading interp; axis equal; view(2); c = colorbar('eastoutside');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
addpath('../../../../MATLAB_codes/BrewerMap/')
n = 256;                     % number of colors
cmap = brewermap(n,'RdBu');  % get RdBu
cmap = flipud(cmap); 
m = max(abs(Z(:)));
caxis([-m m])
colormap(cmap)
colorbar('eastoutside');
colorbar('eastoutside','TickLabelInterpreter','latex');


ax(4) = subplot(2,2,4)

CLV1SPOD1yHF = fCLV1SPODmodeHF(2:2:end);
F = scatteredInterpolant(x, y, double(CLV1SPOD1yHF) , 'linear');
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
title('$v$ $\Phi_{1r}^{0.156}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
shading interp; axis equal; view(2); c = colorbar('eastoutside');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
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

% Titles
St_row = [St1, St2];  

for row = 1:2
  % indices of the two panels in this row
  idx = (row-1)*2 + [1 2];

  % get positions ( [x y w h] in normalized figure units )
  p1 = get(ax(idx(1)),'Position');
  p2 = get(ax(idx(2)),'Position');

  % left edge is x of left panel; 
  % width is from that left edge to the right edge of the right panel
  left = p1(1);
  width = (p2(1)+p2(3)) - 0.05; 

  % vertical placement: 
  %   top row → just above p1(2)+p1(4)
  %   bottom row → just below p1(2)
  if row==1
    y = p1(2) + p1(4) + 0.01;    % tweak the “+0.02” as you like
    vAlign = 'bottom';
  else
    y = p1(2) - 0.05;            % a little below the bottom subplots
    vAlign = 'top';
  end

  % now draw the centered textbox
  annotation('textbox', [left, y, width, 0.03], ...
             'String', sprintf('$St$ = %.3f', St_row(row)), ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment',  vAlign, ...
             'LineStyle', 'none', ...
             'FontSize', fs+4,'Interpreter','latex');
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize', fs)
print(gcf, 'SPODcontours_CLV1_XY', '-dpng', '-r500');    % dpi PNG
exportgraphics(gcf, 'SPODcontours_CLV1_XY.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none');

%% Imaginary part 

file = matfile(['results_CLV1/nfft256_novlp128_nblks61/spod_f' num2str(5,'%.4i')])

set(groot, ...
    'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex');

fCLV1SPODmodeLF = imag(squeeze(file.Psi(:,1)));
% fCLV1SPODmodeLF = imag(squeeze(file.Psi(:,1)));
CLV1SPOD1yLF = fCLV1SPODmodeLF(1:2:end);
% CLV2SPOD1y = sqrt(fCLV2SPODmode10078(1:2:end).^2 + fCLV2SPODmode10078(2:2:end).^2);

x = load('../CLV_visualisation/X.txt');
y = load('../CLV_visualisation/Y.txt');

fig = figure(4);

ax(1) = subplot(2,2,1)

x_length = [-2.5 11.5]; y_length = [-4 4];
% Define a regular grid covering the range of your data:
xq = linspace(-4, max(x), pixels);  % Change '200' as needed
yq = linspace(-8, 8, pixels);
[Xgrid, Ygrid] = meshgrid(xq, yq);
F = scatteredInterpolant(x, y, double(CLV1SPOD1yLF) , 'linear');
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
title('$u$ $\Phi_{1i}^{0.052}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
shading interp; axis equal; view(2); c = colorbar('eastoutside');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
addpath('/media/sidsahu/Elements/PhD_files/MATLAB_codes')
n = 256;                     % number of colors
cmap = brewermap(n,'RdBu');  % get RdBu
cmap = flipud(cmap); 
m = max(abs(Z(:)));
caxis([-m m])
colormap(cmap)
colorbar('eastoutside');
colorbar('eastoutside','TickLabelInterpreter','latex');


ax(2) = subplot(2,2,2)

CLV1SPOD1yLF = fCLV1SPODmodeLF(2:2:end);
F = scatteredInterpolant(x, y, double(CLV1SPOD1yLF) , 'linear');
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
title('$v$ $\Phi_{1i}^{0.052}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
shading interp; axis equal; view(2); c = colorbar('eastoutside');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
addpath('../../../../MATLAB_codes/BrewerMap/')
n = 256;                     % number of colors
cmap = brewermap(n,'RdBu');  % get RdBu
cmap = flipud(cmap); 
m = max(abs(Z(:)));
caxis([-m m])
colormap(cmap)
colorbar('eastoutside');



file = matfile(['results_CLV1/nfft256_novlp128_nblks61/spod_f' num2str(13,'%.4i')])
fCLV1SPODmodeHF = imag(squeeze(file.Psi(:,1)));
% fCLV1SPODmodeHF = imag(squeeze(file.Psi(:,1)));
CLV1SPOD1yHF = fCLV1SPODmodeHF(1:2:end);
colorbar('eastoutside','TickLabelInterpreter','latex');

ax(3) = subplot(2,2,3)

CLV1SPOD1yHF = fCLV1SPODmodeHF(1:2:end);
F = scatteredInterpolant(x, y, double(CLV1SPOD1yHF) , 'linear');
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
title('$u$ $\Phi_{1i}^{0.156}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
shading interp; axis equal; view(2); c = colorbar('eastoutside');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
addpath('../../../../MATLAB_codes/BrewerMap/')
n = 256;                     % number of colors
cmap = brewermap(n,'RdBu');  % get RdBu
cmap = flipud(cmap); 
m = max(abs(Z(:)));
caxis([-m m])
colormap(cmap)
colorbar('eastoutside');
colorbar('eastoutside','TickLabelInterpreter','latex');


ax(4) = subplot(2,2,4)

CLV1SPOD1yHF = fCLV1SPODmodeHF(2:2:end);
F = scatteredInterpolant(x, y, double(CLV1SPOD1yHF) , 'linear');
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
title('$v$ $\Phi_{1i}^{0.156}$');
xlabel('$x$');
ylabel('$y$');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
shading interp; axis equal; view(2); c = colorbar('eastoutside');
xlim([xq(1) xq(end)]); ylim([yq(1) yq(end)]);
xticks(0:4:40); yticks(-30:4:30);
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

% Titles
St_row = [St1, St2];  

for row = 1:2
  % indices of the two panels in this row
  idx = (row-1)*2 + [1 2];

  % get positions ( [x y w h] in normalized figure units )
  p1 = get(ax(idx(1)),'Position');
  p2 = get(ax(idx(2)),'Position');

  % left edge is x of left panel; 
  % width is from that left edge to the right edge of the right panel
  left = p1(1);
  width = (p2(1)+p2(3)) - 0.05; 

  % vertical placement: 
  %   top row → just above p1(2)+p1(4)
  %   bottom row → just below p1(2)
  if row==1
    y = p1(2) + p1(4) + 0.01;    % tweak the “+0.02” as you like
    vAlign = 'bottom';
  else
    y = p1(2) - 0.05;            % a little below the bottom subplots
    vAlign = 'top';
  end

  % now draw the centered textbox
  annotation('textbox', [left, y, width, 0.03], ...
             'String', sprintf('$St$ = %.3f', St_row(row)), ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment',  vAlign, ...
             'LineStyle', 'none', ...
             'FontSize', fs+4,'Interpreter','latex');
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize', fs)
print(gcf, 'SPODcontours_CLV1_XY_imag', '-dpng', '-r500');    % dpi PNG
exportgraphics(gcf, 'SPODcontours_CLV1_XY_imag.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none');