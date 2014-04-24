%% Plot Doyle-Fuller-Newman Model Results
%   Created May 23, 2012 by Scott Moura
close all;

figure(1)
clf

subplot(3,1,1)
plot(t,I,'LineWidth',2)
ylabel('Current')

subplot(3,1,2)
plot(t,Volt,'LineWidth',2)
ylabel('Voltage')

subplot(3,1,3)
plot(t,SOC,'LineWidth',2)
ylabel('SOC')
xlabel('Time [sec]')

pause;

%% Animation

for k = 1:NT

    figure(2)
    set(gcf,'Position',[103 3 1151 673]);
    
    subplot(3,3,1)
    cla
    plot(1:Nn,c_ss_n(:,k)/p.c_s_n_max,1:Nn,c_avg_n(:,k)/p.c_s_n_max,'LineWidth',2);
%     plot(1:Nn,phi_s_n(:,k),'LineWidth',2);
    xlim([1 Nn])
    ylim([min(min(c_ss_n/p.c_s_n_max)) max(max(c_ss_n/p.c_s_n_max))])
%     ylim([min(min(phi_s_n)) max(max(phi_s_n))])
    ylabel('$$c_{ss,n}(x,t)$$','interpreter','latex')
%     ylabel('$$\phi_{s,n}$$','interpreter','latex')
    
    subplot(3,3,3)
    cla 
    plot(1:Np,c_ss_p(:,k)/p.c_s_p_max,1:Np,c_avg_p(:,k)/p.c_s_p_max,'LineWidth',2);
%     plot(1:Np,phi_s_p(:,k),'LineWidth',2)
    xlim([1 Np])
    ylim([min(min(c_ss_p/p.c_s_p_max)) max(max(c_ss_p/p.c_s_p_max))])
%     ylim([min(min(phi_s_p)) max(max(phi_s_p))])
    ylabel('$$c_{ss,p}(x,t)$$','interpreter','latex')
%     ylabel('$$\phi_{s,p}$$','interpreter','latex')
    
    jall = [jn(:,k); zeros(p.Nxs-1,1); jp(:,k)];
    etaall = [eta_n(:,k); zeros(p.Nxs-1,1); eta_p(:,k)];
    ieall = [i_en(:,k); I(k)*ones(p.Nxs-1,1); i_ep(:,k)];

    subplot(3,3,[4 5 6])
    cla
    plot(1:Nx,jall,'LineWidth',2);
    xlim([1 Nx])
%     ylim([min(min(i_en)), max(max(i_en))])
    ylabel('$$i_{e,n}(x,t)$$','interpreter','latex')

%     subplot(3,3,6)
%     cla
%     plot(1:Np,jp(:,k),'LineWidth',2);
%     xlim([1 Np])
% %     ylim([min(min(i_ep)), max(max(i_ep))])
%     ylabel('$$i_{e,p}(x,t)$$','interpreter','latex')

    subplot(3,3,[7 8 9])
    cla
    plot(1:(p.Nx+1), c_ex(:,k),'LineWidth',2);
%     plot(1:(p.Nx-3), phi_e(:,k),'LineWidth',2)
%     plot(1:Nx,etaall, 'LineWidth',2)
    xlim([1 p.Nx+1])
    ylabel('$$c_e(x,t)$$','interpreter','latex')
%     ylabel('$$\phi_e$$','interpreter','latex')
    ylim([min(min(c_ex)), max(max(c_ex))])
%     ylim([min(min(phi_e)), max(max(phi_e))])

    pause(0.1);

end

%% Animation for Seminars
Icrate = I/35;

[c_sr_n, c_sr_p] = sim_cs(p,t,jn,jp,csn0,csp0);
th_sr_n = c_sr_n / p.c_s_n_max;
th_sr_p = c_sr_p / p.c_s_p_max;

%%
figure(3)
clf;

% vidObj = VideoWriter('img/dfn_5Cdispulse.avi');
% vidObj.FrameRate = 10;
% vidObj.Quality = 100;
% open(vidObj);

for k = 3:NT
    
    figure(3)
    set(gcf,'Position',[160 -16 850 700],'PaperPositionMode','auto')
    
    % Anode concentration
    vn = linspace(0,0.9,101);
    subplot(3,3,1)
    contourf(fliplr(th_sr_n(:,:,k))',vn);
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    ylabel('Space, r','FontSize',18)
    title('$$c_s^-(x,r,t)$$','FontSize',18,'Interpreter','latex')
    set(gca,'FontSize',16)
    set(gca,'Position',[0.08 0.63 0.35 0.33])
    
    % Cathode concentration
    vp = linspace(0.5,1,101);
    subplot(3,3,3)
    contourf(fliplr(th_sr_p(:,:,k))',vp);
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    ylabel('Space, r','FontSize',18)
    title('$$c_s^+(x,r,t)$$','FontSize',18,'Interpreter','latex')
    set(gca,'FontSize',16)
    set(gca,'Position',[0.63 0.63 0.35 0.33])
    
    % Electrolyte concentration
    subplot(3,3,[4 5 6])
    plot(linspace(0,1,Nx+4),c_ex(:,k)/1e3,'LineWidth',2);
    xlim([0 1])
    ylim([0 1.5])
    xlabh2 = xlabel('Space, x/L','FontSize',18);
    set(xlabh2,'Position',[0.5 -0.15 1]);
    ylabel('$$c_e(x,t)$$','FontSize',18,'Interpreter','latex')
    set(gca,'FontSize',16)
    set(gca,'Position',[0.08 0.34 0.9 0.25])
    
    % Current
    subplot(3,3,[7 8 9])
    plot(t,Icrate,'LineWidth',2); hold on;
%     plot(t(k), Icrate(k), 'ro', 'MarkerSize',16, 'MarkerFaceColor','r'); hold off;
    plot(t,Volt,'Color',[0 0.5 0],'LineWidth',2);
    plot([t(k) t(k)], [0 6],'k--'); hold off;
    xlim([0 t(end)]);
    xlabh = xlabel('Time, t','FontSize',18);
    set(xlabh, 'Position',[120, -0.5, 1])
    ylabel('Current (B) | Voltage (G)','FontSize',18)
    set(gca,'FontSize',16)
    set(gca,'Position',[0.08 0.08 0.9 0.18])
    
%     F = getframe(gcf);
%     writeVideo(vidObj,F);
    
end

% close(vidObj);