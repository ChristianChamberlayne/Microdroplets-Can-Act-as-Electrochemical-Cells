%Code used to make figures 2, 4 and 8

%start timer
tic()

%use 1*10^-7 for figures 2A and 4A and figure 8
%use 1*10^-6 for figures 2B and 4B
TargetSalt_Molar = 1*10^-7;

%use .01 for figures 2 and 4
%use 1 for figure 8
TargetS_ionspernm2 = 0.01;

R = 2.5;                                                   % radius of microdroplet

% all figures in paper use SurfactantIsSalt = true
SurfactantIsSalt = true;
SurfactantIsAcid = false;

manyGraphs = false;
verbose = true;
guessPoint = [80, 10]; %in kions 
PointsInDerivative = 11;   %use odd numbers so a central value exists

%use .005 for figures 2 and 4
%use .000001 for figure 8
DerivativeStepRange = .005; 

%constants at 25C
epsilon = 78.4;                                                             %dielectric for water at 25C 
epsilon0 = 8.854*10^-12;                                                    %in F/m   
ThermalVoltage25 = 25.6796 *10^-3;                                          %thermal voltage in V at 25C; 
F = 96485.33;                                                               %in C/mol
N_A = 6.0221409 *10^23;                                                     %ions per mole
K = F / (epsilon*epsilon0*ThermalVoltage25);                                %solves for constant K from other constants
K_w = 10^-14;                                                               %K_W in M^2

%calc volume (L)
DropletVolumeL = 4/3 * pi * R.^3 *10^-15;  
DropletSurfaceAreaNM2 = 4 * pi * R.^2 * 10^6;  

%convert everything um and kIons units
convertC = N_A / 1000 * 10^-15;                                              % converts M to kIOns / um^3
convertV = 10^12;                                                             % converts V (kg?m2?s?3?A?1) to new units (kg?um2?s?3?A?1)  
convert_kIons2M = 1000 / N_A / DropletVolumeL;
convert_kIons2IonsPerNM2 = 1000 / DropletSurfaceAreaNM2;

%converting constants used
K_w = K_w * convertC * convertC;                                               
K = K * 10^6 * 1000 / N_A;                                                  % convert from m/mole to um / kIons 
ThermalVoltage25 = ThermalVoltage25 * convertV;
F = F * 1000 / N_A;



if SurfactantIsSalt && SurfactantIsAcid
        disp('please pick only one: is the surfactant an acid or its salt?')
        exit()
        
elseif SurfactantIsSalt
       TargetP_kions =  (TargetSalt_Molar .* DropletVolumeL .* N_A + TargetS_ionspernm2 .* DropletSurfaceAreaNM2)./1000;
       TargetNS_kions = (TargetSalt_Molar .* DropletVolumeL .* N_A + TargetS_ionspernm2 .* DropletSurfaceAreaNM2)./1000;
       TargetS_kions =  (TargetS_ionspernm2 * DropletSurfaceAreaNM2)/1000;
       
elseif SurfactantIsAcid
       TargetP_kions =  (TargetSalt_Molar .* DropletVolumeL .* N_A)./1000;
       TargetNS_kions = (TargetSalt_Molar .* DropletVolumeL .* N_A + TargetS_ionspernm2 .* DropletSurfaceAreaNM2)./1000;
       TargetS_kions =  (TargetS_ionspernm2 * DropletSurfaceAreaNM2)/1000;    
       
else
        disp('please pick one: is the surfactant an acid or its salt?')
        exit()
end

[U,storedPaths,RadialSolutions]  = Solve_ion_to_ion_V4(TargetP_kions, TargetNS_kions, TargetS_kions, R, guessPoint,PointsInDerivative,DerivativeStepRange, manyGraphs,verbose);

figure()
plot(R * 2, squeeze(U(:,1,:))*1000,'o-','linewidth',1.5,'markersize',8)
legend('Total Voltage', 'Electrostatics Voltage', 'Water Recombination Voltge')
ylabel('Voltage mV')
xlabel('Droplet Diameter \mum') 


if SurfactantIsAcid
textbox = sprintf(' Surface Charge Density:      %#.2g ions per nm^2 \n Average Salt Concentration: %#.2g M \n Acidic Surfactant', TargetS_ionspernm2,TargetSalt_Molar) 
end
if SurfactantIsSalt
textbox = sprintf(' Surface Charge Density:      %#.2g ions per nm^2 \n Average Salt Concentration: %#.2g M \n Surfactant Salt', TargetS_ionspernm2, TargetSalt_Molar) 
end
   
annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', textbox)

%begin graphs
%get average plot
middle = ceil(PointsInDerivative/2); 
if rem(PointsInDerivative,2) == 0; disp('please use an odd value for PointsInDerivative, rounding down to determine "middle" point'); end

xint = cell2mat(RadialSolutions{middle}(3));
Sxint = cell2mat(RadialSolutions{middle}(2));
RatioPos = RadialSolutions{middle}{1}(4);
RatioNeg = RadialSolutions{middle}{1}(5);

%find ion distributions
middleH = Sxint(1,:);
middleOH = K_w ./ Sxint(1,:);
middleP = RatioPos * middleH;
middleN = RatioNeg * middleOH;

%used for graph change in ions vs S later
for i = 1:PointsInDerivative
[TotalNS(i), TotalP(i), TotalS(i), U_elec(i), TotalOH(i), TotalH(i)] = DropletValues(RadialSolutions{i}{:});
 TotalN(i) = TotalNS(i) - TotalS(i);
end


i = ceil(PointsInDerivative/2);
    %set up results
    xint = cell2mat(RadialSolutions{i}(3));
    Sxint = cell2mat(RadialSolutions{i}(2));
    RatioPos = RadialSolutions{i}{1}(4);
    RatioNeg = RadialSolutions{i}{1}(5);
    
    %find ion distributions
    H = Sxint(1,:);
    OH = K_w ./ Sxint(1,:);
    P = RatioPos * H;
    N = RatioNeg * OH;
    
    
    % graph change in ion conc radially
    linefade = linspace(1,.2,PointsInDerivative);
     
    fig2 = figure(2);
    
    plot(xint, OH / convertC, 'DisplayName', 'OH^- (M)','LineWidth',1);
    hold on
    plot(xint, H / convertC, 'DisplayName', 'H^+ (M)','LineWidth',1);
    plot(xint, P / convertC, 'DisplayName', '+ ion (M)','LineWidth',1);
    plot(xint, N / convertC, 'DisplayName', '- ion (M)','LineWidth',1);
    set(gca,'fontsize', 18)
    ylabel('Concentration (M)');
    xlabel('Radial Distance (\mum)');
    set(gca,'Yscale','log');
    ylim([10^-10 10^-2])
    %title(sprintf('changing surface charge')); 
    legend('show','Location', 'northwest','fontsize', 12)
    

%graph change in ions vs S 

figure(4) 
plot(TotalS*1000/DropletSurfaceAreaNM2,(TotalOH - TotalOH(middle))*1000,'DisplayName', 'OH^- ions','LineWidth',2,'Marker','o')
hold on
plot(TotalS*1000/DropletSurfaceAreaNM2,(TotalH - TotalH(middle))*1000,'DisplayName', 'H^+ ions','LineWidth',2,'Marker','x')
plot(TotalS*1000/DropletSurfaceAreaNM2,(TotalP - TotalP(middle))*1000,'DisplayName', 'Positive strong electrolyte','LineWidth',2,'Marker','o')
plot(TotalS*1000/DropletSurfaceAreaNM2,(TotalN - TotalN(middle))*1000,'DisplayName', 'Negative strong electrolyte','LineWidth',2,'Marker','o')
plot(TotalS*1000/DropletSurfaceAreaNM2,(TotalS - TotalS(middle))*1000,'DisplayName', 'Surface ions','LineWidth',2,'Marker','o')
set(gca,'fontsize', 18)
legend()
ylabel('Change in ion count');
xlabel('Surface charge (ions/nm^2)');
title(sprintf('Changing surface charge'));

%graph H ions and U 
figure(5) 
set(gca,'fontsize', 18)
ylabel('Electrostatic energy (fJ)');
%plot in 10^-15 J i.e. fJ, was originally in pJ
plot(TotalS*1000/DropletSurfaceAreaNM2, U_elec * 10^3,'DisplayName', 'Electrostatic energy','LineWidth',2,'Marker','o')

hold on 

yyaxis right
plot(TotalS*1000/DropletSurfaceAreaNM2,TotalH*1000,'DisplayName', 'H^+ ions','LineWidth',2,'Marker','o')
ylabel('H^+ ion count');

xlabel('Surface charge (ions/nm^2)');
title(sprintf('Changing surface charge'));
legend('show','Location', 'northwest','fontsize', 12)
%record time of run
toc()