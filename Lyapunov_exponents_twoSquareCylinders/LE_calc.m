clear all
close all
clc

stept = 0.3;
load("R_full.txt")

k=1;
T_final = 5400;

for i=1:(size(R_full,1)/6)
    k = (i-1)*6;
    LE(i,1) = i*stept;
    if(i==1)
        LE(i,2) = log(R_full(k+1,1))/LE(i,1);
        LE(i,3) = log(R_full(k+2,2))/LE(i,1);
        LE(i,4) = log(R_full(k+3,3))/LE(i,1);
        LE(i,5) = log(R_full(k+4,4))/LE(i,1);
        LE(i,6) = log(R_full(k+5,5))/LE(i,1);
        LE(i,7) = log(R_full(k+6,6))/LE(i,1);
    else
        LE(i,2) = LE(i-1,2)*(i-1)/i + log(R_full(k+1,1))/LE(i,1);
        LE(i,3) = LE(i-1,3)*(i-1)/i + log(R_full(k+2,2))/LE(i,1);
        LE(i,4) = LE(i-1,4)*(i-1)/i + log(R_full(k+3,3))/LE(i,1);
        LE(i,5) = LE(i-1,5)*(i-1)/i + log(R_full(k+4,4))/LE(i,1);
        LE(i,6) = LE(i-1,6)*(i-1)/i + log(R_full(k+5,5))/LE(i,1);
        LE(i,7) = LE(i-1,7)*(i-1)/i + log(R_full(k+6,6))/LE(i,1);
    end
end

LE(end,1)

%% 

fgh = figure(1);

zm = 0.01;

yBox = [zm, zm, 0.15, 0.15, 0];
xBox = [0, T_final, T_final, 0, 0];
retained = patch(xBox, yBox, 'red', 'FaceColor', 'red', 'EdgeColor','none','FaceAlpha', 0.1);

yBox = [-zm, -zm, -0.12, -0.12, 0];
xBox = [0, T_final, T_final, 0, 0];
truncated = patch(xBox, yBox, 'green', 'FaceColor','green', 'EdgeColor','none','FaceAlpha', 0.2);

yBox = [-zm, -zm, zm, zm, 0];
xBox = [0, T_final, T_final, 0, 0];
neutral = patch(xBox, yBox, 'yellow', 'FaceColor', 'yellow', 'EdgeColor','none','FaceAlpha', 0.5);
hold on;

legend([retained, truncated, neutral],'Unstable manifold','Stable manifold','Neutral manifold','Fontsize',14,'interpreter','latex'); 

gap = 50;
for i = 1:6
    plot(LE(1:gap:end,1), LE(1:gap:end,i+1),'-', 'DisplayName', "LE " + i, 'LineWidth', 1.5)
    % plot(LE(1:100:end,1), LE(1:100:end,i+1),'-')
    hold on
end
set(gcf,'Renderer','opengl', ...
        'GraphicsSmoothing','on');
xlabel('Time','Interpreter','latex')
ylabel('Lyapunov exponent','Interpreter','latex')
title("Lyapunov spectrum, $\mathbf{t_{step}}$ = " + stept, 'Interpreter','latex')
legend('show', 'Location', 'southeast')
grid on
grid minor

% Add dashed vertical lines at x = 2100 and x = 4500
% Add dashed vertical lines without appearing in the legend
% xline(2100, '--', 'LineWidth', 2.5, 'HandleVisibility', 'off');
% xline(4500, '--', 'LineWidth', 2.5, 'HandleVisibility', 'off');

ylim([-0.08 0.15])
% writematrix(LE, 'LE_1.txt', 'Delimiter', 'tab');
xlim([0 5400])
legend('Location', 'eastoutside', 'FontSize', 18)
set(gca, 'Fontsize', 18); 
set(gca, 'FontSize', 18, 'FontName', 'Courier')
fgh.Position = [680 458 1300 500];

% Save the figure
% saveas(gcf, 'LEspectrumvsTime.png')
print(gcf, 'LEvstime1.png', '-dpng', '-r500');

%%
figure(3)

yBox = [zm, zm, 0.12, 0.12, 0];
xBox = [1, 6, 6, 1, 1];
retained = patch(xBox, yBox, 'red', 'FaceColor', 'red', 'EdgeColor','none','FaceAlpha', 0.1);

yBox = [-zm, -zm, -0.12, -0.12, 0];
xBox = [1, 6, 6, 1, 1];
truncated = patch(xBox, yBox, 'green', 'FaceColor','green', 'EdgeColor','none','FaceAlpha', 0.2);

yBox = [-zm, -zm, zm, zm, 0];
xBox = [1, 6, 6, 1, 1];
neutral = patch(xBox, yBox, 'yellow', 'FaceColor', 'yellow', 'EdgeColor','none','FaceAlpha', 0.5);
hold on;


plot1 = plot(1:6,sort(LE(end,2:7),2,'descend'),'-','LineWidth',1.5,'Marker','diamond','MarkerSize',10,'MarkerFaceColor','blue','DisplayName','$\mathbf{t_{step}}=0.3$');
hold on
plot2 = plot(1:6, [0.1036    0.0214   -0.0033   -0.0361   -0.0552   -0.0677],'-','LineWidth',1.5,'Marker','square','MarkerSize',10,'MarkerFaceColor','none','DisplayName','$\mathbf{t_{step}}=10$');
xlabel('Index','Interpreter','latex')
ylabel('Lyapunov exponent','Interpreter','latex')
title('Lyapunov spectrum ','Interpreter','latex')
grid on

% Combine all legend entries
legend([plot1, plot2, retained, truncated, neutral], ...
       {'$\mathbf{t_{step}}=0.3$', '$\mathbf{t_{step}}=10$', ...
        'Unstable Manifold', 'Stable Manifold', 'Neutral Manifold'}, ...
        'FontSize',14,'Interpreter','latex');
    
set(gca,'FontSize',14); 
xlim([1 6])
ylim([-0.12 0.12])

