%start clock
tic()

starting_surface_ions_per_nm2 = .01;
end_surface_ions_per_nm2 = .00001;
PointsAlongGraph = 31;

StartingN_Molar = 10^-7;
%this code assumes P ions = N ions + surfaceions. This does not have to be true

R = 2.5;

manyGraphs = false;
verbose = true;

%settings for solutions
guessPoint = [50 15];
PointsInDerivative = 21;

%constants at 25C

epsilon = 78.4;                                                             %dielectric for water at 25C 
epsilon0 = 8.854*10^-12;                                                    %in F/m   
ThermalVoltage25 = 25.6796 *10^-3;                                          %thermal voltage in V at 25C; 
F = 96485.33;                                                               %in C/mol
N_A = 6.0221409 *10^23;                                                     %ions per mole
K = F / (epsilon*epsilon0*ThermalVoltage25);                                %solves for constant K from other constants
K_w = 10^-14;                                                               %K_W in M^2

%calc volume (L and um^3)
DropletVolumeL = 4/3 * pi * R^3 *10^-15;
DropletVolumeUm3 = 4/3 * pi * R^3; 
DropletSurfaceAreaNM2 = 4 * pi * R^2 * 10^6;  

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


%_________________________________________________________________________
%solve for solution constants P S and NS ions and convert to kions units
TargetS_ions_per_nm2 = logspace(log10(starting_surface_ions_per_nm2),log10(end_surface_ions_per_nm2),PointsAlongGraph);
%TargetS_ions_per_nm2 = [linspace(starting_surface_ions_per_nm2,end_surface_ions_per_nm2,PointsAlongGraph) .001 * [7 6 5 4 3 2]];
%PointsAlongGraph = length(TargetS_ions_per_nm2);
TargetS_kions = TargetS_ions_per_nm2 .* DropletSurfaceAreaNM2 ./ 1000;

StartingN_kions = StartingN_Molar / convert_kIons2M;
StartingS_kions = TargetS_kions(1);
StartingP_kions = StartingS_kions + StartingN_kions;

%Note: other target P definitions can be used, this one assumes an even number of strong electrolytes in precursor solution
TargetP_kions = TargetS_kions + StartingN_kions; % in kions
TargetDiffPS_kions = StartingP_kions - StartingS_kions % in kions
TargetN_kions = StartingN_kions;
%_______________________________________________________________________
%parameters to pass to functions
DropParam = [R, K, K_w, F];
params = [ThermalVoltage25, K_w, F];

U = cell(1,PointsAlongGraph);
storePath = cell(1,PointsAlongGraph);
RadialSolutions = cell(1,PointsAlongGraph);


parfor k = 1:PointsAlongGraph
       [U{k},storePath{k},RadialSolutions{k}] = Solve_ion_to_nuetral(TargetN_kions, TargetDiffPS_kions, TargetS_kions(k), R, guessPoint,PointsInDerivative,manyGraphs,verbose);
end

for k = 1:PointsAlongGraph
    
    if ~any(isnan(U{k}(:))) && ~isempty(U{k})
        UgradV(k)      = U{k}(1,1);
        UgradV_elec(k) = U{k}(2,1);
        UgradV_H2O(k)  = U{k}(3,1);

        U_r2(:,k) = [U{k}(1,2);U{k}(2,2);U{k}(3,2);];
    else
        UgradV(k)      = NaN;
        UgradV_elec(k) = NaN;
        UgradV_H2O(k)  = NaN;

        U_r2(:,k) = [NaN;NaN;NaN];
    end
end

figure(106)
    plot(TargetS_ions_per_nm2,UgradV *1000  ,'-o','DisplayName', 'microdroplet voltage'    ,'LineWidth',2)
    hold on
    plot(TargetS_ions_per_nm2,UgradV_elec * 1000,'-o','DisplayName', 'electrostatics'          ,'LineWidth',2,'color','green')
    plot(TargetS_ions_per_nm2,UgradV_H2O * 1000 ,'-o','DisplayName', 'water recombination'     ,'LineWidth',2, 'color', 'red')
      
    ylabel('Microdroplet reaction driving force (mV)')
    xlabel('Surface-bound ion density (ions nm^-^2)')
    title('Reaction driving force along reaction path')
    hold off
    
    


