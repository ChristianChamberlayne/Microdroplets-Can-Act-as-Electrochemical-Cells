

function getRadialPicture(RadialSolutions)
    ThermalVoltage25 = 0.025679600000000;

    %convert everything um and kIons units
    N_A = 6.0221409 *10^23; 
    convertC = N_A / 1000 * 10^-15;                                              % converts M to kIOns / um^3
    convertV = 10^12;  
   
    
    %remove data
    
    K_w = RadialSolutions{1}{1}(2);
    
    
    
for i = 3
    %set up results
    xint = cell2mat(RadialSolutions{i}(3));
    Sxint = cell2mat(RadialSolutions{i}(2));
    RatioPos = RadialSolutions{i}{1}(4);
    RatioNeg = RadialSolutions{i}{1}(5);
    DropletSurfaceAreaNM2 = 4 * pi * xint(end)^2 * 10^6;  
    

    %find ion distributions
    H = Sxint(1,:);
    OH = K_w ./ Sxint(1,:);
    P = RatioPos * H;
    N = RatioNeg * OH;
    
    linefade = linspace(1,.2,5);
    
    
    
    %set up ion concentration figure
     
    fig1 = figure(1);
   
    hold off
    
    plot(xint, H / convertC, 'DisplayName', 'H^+ (M)','LineWidth',2)
    hold on
    plot(xint, OH / convertC, 'DisplayName', 'OH^- (M)','LineWidth',2)
    plot(xint, P / convertC, 'DisplayName', '+ ion (M)','LineWidth',2)
    plot(xint, N / convertC, 'DisplayName', '- ion (M)','LineWidth',2)
    set(gca,'fontsize', 18)
    ylabel('Concentration (M)');
    xlabel('Radial Distance (\mum)');
    set(gca,'Yscale','log');
    ylim([10^-10 10^-2])
    title(sprintf('changing surface charge')); 
    legend('show','Location', 'northwest','fontsize', 12)
    movieFrames(6-i) = getframe;
    
    fig2 = figure(2);
    hold on
        
    line(4*i - 3) = plot(xint, H / convertC, 'DisplayName', 'H^+ (M)','LineWidth',1);
    line(4*i - 3).Color = [1 0 0 linefade(i)];
    line(4*i - 2) = plot(xint, OH / convertC, 'DisplayName', 'OH^- (M)','LineWidth',1);
    line(4*i - 2).Color = [0 0 1 linefade(i)];
    line(4*i - 1) =plot(xint, P / convertC, 'DisplayName', '+ ion (M)','LineWidth',1);
    line(4*i - 1).Color = [0 1 0 linefade(i)];
    line(4*i) = plot(xint, N / convertC, 'DisplayName', '- ion (M)','LineWidth',1);
    line(4*i).Color = [1 0 1 linefade(i)];
    set(gca,'fontsize', 18)
    ylabel('Concentration (M)');
    xlabel('Radial Distance (\mum)');
    set(gca,'Yscale','log');
    ylim([10^-10 10^-2])
    title(sprintf('changing surface charge')); 
    legend('show','Location', 'northwest','fontsize', 12)
    legend({'H^+ (M)';'OH^- (M)';'+ ion (M)';'- ion (M)'})
    
    [TotalNS(i), TotalP(i), TotalS(i), U_elec(i), TotalOH(i), TotalH(i)] = DropletValues(RadialSolutions{i}{:});
    
end

%{
figure(3)
hold on
    plot(TotalS*1000/DropletSurfaceAreaNM2, TotalH*1000, 'DisplayName', 'H^+','LineWidth',2)
    %plot(TotalS, TotalOH-TotalOH(3), 'DisplayName', 'OH^- ','LineWidth',2)
    %plot(TotalS, TotalP-TotalP(3), 'DisplayName', '+ ion ','LineWidth',2)
    %plot(TotalS, (TotalNS - TotalS)-(TotalNS(3) - TotalS(3)), 'DisplayName', '- ion','LineWidth',2)
    %plot(TotalS, TotalNS-TotalNS(3), 'DisplayName', '- ion + surface ions (M)','LineWidth',2)
    set(gca,'fontsize', 18)
    ylabel('Ion Count in Microdroplet');
    xlabel('Radial Distance (\mum)');
    title(sprintf('Redox Reaction Effect on H^+ Conc.')); 
    legend('show','Location', 'northwest','fontsize', 12)
    hold off
%}

 %find center-to-wall potential
            E = ThermalVoltage25*Sxint(2,:)./Sxint(1,:);                   %in V/um
            V = -cumtrapz(xint,E);                                         %in V
figure(4)
hold on
    plot(xint, E * 10^6,'LineWidth',2)
    set(gca,'fontsize', 18)
    xlabel('Radial Distance (\mum)');
    ylabel('E (V/m)');
    set(gca,'Yscale','log');
    %ylim([10^-10 10^-2])
    title(sprintf('Electric Field Strength')); 
    
figure(5)
hold on
    plot(xint, V * 1000,'LineWidth',2)
    set(gca,'fontsize', 18)
    xlabel('Radial Distance (\mum)');
    ylabel('\phi (mV)');
    %set(gca,'Yscale','log');
    %ylim([10^-10 10^-2])
    title(sprintf('Electric Potential')); 


figure(1)

%{
vw = VideoWriter('radial_video.avi');  %taking a guess that you intend to modify the filename each time you write a video
vw.FrameRate = 1;
open(vw);
writeVideo(vw, movieFrames);
close(vw)
%}
end