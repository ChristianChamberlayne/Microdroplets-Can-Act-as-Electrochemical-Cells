
% HAS ITS OWN CONSTANTS STORED. BE CAREFUL CHANGING UNITS

function [U,storePath,RadialSolutions] = Solve1D_V3(TargetP_kions, TargetNS_kions, TargetS_kions, R, guessPoint,manyGraphs,verbose)
%_________________________________________________________________________
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

%settings for solutions
precisionP =  10^-3; %axis X
precisionNS =  10^-3; %axis Y
maxloops = 50;
maxSolveTries = 5;
maxStepRestrictionsBeforeReset = 5;
PointsInDerivative = 5;   % later, make sure this matches number of points in TargetS
DerivativeStepRange = .001;
showSolutionGraph = manyGraphs;  % this is for trouble shooting midfit function, turn off to speed up solve. 

TargetP = TargetP_kions; % in kions
TargetNS = TargetNS_kions; % in kions
TargetS = TargetS_kions * linspace(1 - DerivativeStepRange,1 + DerivativeStepRange,PointsInDerivative); 

%Check inputs are good
if TargetNS_kions < TargetS_kions
    disp('Impossible targets. \n Surface ions can not exceed the sum of surface ions + negative ions. \aborting solve of this point')
     U = [NaN, NaN;...
          NaN, NaN;...
          NaN, NaN];
     storePath = NaN;
     RadialSolutions = NaN;
    return
end

if TargetNS <= TargetS(end)
    disp('Derivative step size stepped to an impossible target. \n Squishing nearby derivative step sizes to avoid target surface ions exceeding the sum of surface ions + negative ions. \n continuing solve of this point...')
    maxDerivativeStep =  (TargetNS_kions - TargetS_kions) * .99;
    TargetS(ceil(PointsInDerivative/2):end) = linspace(TargetS(ceil(PointsInDerivative/2)),TargetS_kions + maxDerivativeStep,ceil(PointsInDerivative/2)); 
end

%parameters to pass to functions
DropParam = [R, K, K_w, F];
params = [ThermalVoltage25, K_w, F];


%declare loop output variables
MolesNegativeIons = zeros(1,maxloops);       %M_n + M_s
MolesPositiveIons = zeros(1,maxloops);       %M_p
MolesSurfaceIons = zeros(1,maxloops);        %M_s
StepsToSolve = ones(1,PointsInDerivative) * PointsInDerivative; %this gets reset if it uses less than max loops

%Set Path of first 2 points (2 points are needed to extrapolate new points on path)
PathRatioPos = [.99*guessPoint(1),guessPoint(1)];
rGuess = 48/50 * R;

%Calc first 2 points for first loop
    if verbose,fprintf('Starting Run. Loops Tracking Below:\n j = 1: i:'),end 
for i = 1:2
    if verbose, fprintf('%i ',i), end
        for k = 1:maxSolveTries 
        %calculate the function for RatioNeg dependence on RatioPos that would result with given targets.
        %this equation enforces that the net charge on the microdroplet is 0. Modify for charged droplets 
        RatioNeg = -(TargetNS - TargetS(1)) / (TargetS(1) - TargetP * (1/PathRatioPos(i) + 1) + TargetNS - TargetS(1));

        %calculate the function for innerConc. uses the average H value of the droplet
        %this equation enforces that the net charge on the microdroplet is 0 to get a valid inner conc value. Modify for charged droplets.
        innerConc = (TargetS_kions + sqrt(TargetS_kions^2 + 4 * (1 + PathRatioPos(i)) * K_w * (1 + RatioNeg)))/ (2*(1+PathRatioPos(i))*DropletVolumeUm3);        


        [Sxint, xint, sucessfulSolve,rGuess] = midfit_PTarget_V3(TargetP_kions,PathRatioPos(i),RatioNeg,DropParam,innerConc,rGuess,showSolutionGraph);
        %saves r value as guess for next solve. speeds things up a bit. 
        %check solution OK
        if sucessfulSolve 
            break
        elseif k ~= maxSolveTries
            PathRatioPos(i) = PathRatioPos(i) / 2;            
        else
                %deal with this later
                sprintf('starting points give bad solution, redo starting points')
                U = [NaN, NaN;...
                     NaN, NaN;...
                     NaN, NaN];
                 storePath = NaN;
                 RadialSolutions = NaN;
                return
        end
     end
%get droplet values
[MolesNegativeIons(i), MolesPositiveIons(i), MolesSurfaceIons(i), ~, ~, ~]...
 = DropletValues([params, PathRatioPos(i), RatioNeg],Sxint, xint, sucessfulSolve);
end

%loops over each derivative point
for j = 1:PointsInDerivative

if (j ~= 1)  %skip getting points from previous loop solution if first loop 
    if verbose,fprintf('\n j = %i: i:',j), end 
    
    %Set last point of previous loop + 2 new scaled points as new starting set of points for next loop. 
    PathRatioPos = [PathRatioPos(StepsToSolve(j-1))*.99,PathRatioPos(StepsToSolve(j-1))];
    
    %solve this first set of new points for next loop starting point
    for i = 1:2
        if verbose,fprintf('%i ',i),end
        %calculate the function for RatioNeg dependence on RatioPos that would result with given targets.
        %this equation enforces that the net charge on the microdroplet is 0. Modify for charged droplets 
        RatioNeg = -(TargetNS - TargetS(j)) / (TargetS(j) - TargetP * (1/PathRatioPos(i) + 1) + TargetNS - TargetS(j));
        
        %calculate the function for innerConc. uses the average H value of the droplet
        %this equation enforces that the net charge on the microdroplet is 0 to get a valid inner conc value. Modify for charged droplets.
        innerConc = (TargetS_kions + sqrt(TargetS_kions^2 + 4 * (1 + PathRatioPos(i)) * K_w * (1 + RatioNeg)))/ (2*(1+PathRatioPos(i))*DropletVolumeUm3);        
        
        [Sxint, xint, sucessfulSolve,rGuess] = midfit_PTarget_V3(TargetP_kions,PathRatioPos(i),RatioNeg,DropParam,innerConc,rGuess,showSolutionGraph);
        %check solution OK
        if sucessfulSolve == false
            
            %deal with this later
            sprintf('Initial points give bad solution for loop j = %i, suggest a smaller step size for EdgeH',j)
            U = [NaN, NaN;...
                 NaN, NaN;...
                 NaN, NaN];
            storePath = NaN;
            RadialSolutions = NaN;
            return
        end
        %get droplet values
        [MolesNegativeIons(i), MolesPositiveIons(i), MolesSurfaceIons(i), ~, ~, ~]...
         = DropletValues([params, PathRatioPos(i), RatioNeg],Sxint, xint, sucessfulSolve);
    end
end


for i = 3:maxloops
    
    if verbose,fprintf('%i ',i),end
    
           
    %linear extrapolation from last two point's values to next guess
    pointExtrapolation = (TargetNS - MolesNegativeIons(i-1))*(PathRatioPos(i-2) - PathRatioPos(i-1))/(MolesNegativeIons(i-2) - MolesNegativeIons(i-1)) + PathRatioPos(i-1);

    %make sure guess isn't extrapolating beyond previous max guesses. 
    %if pointExtrapolation > RatioPGuessMax
    %   pointExtrapolation = (RatioPGuessMax + PathRatioPos(i)) / 2; %if so, just go half way toward that edge 
    
    % is step more than doubling
    if pointExtrapolation > PathRatioPos(i-1) * 2   
        pointExtrapolation = PathRatioPos(i-1) * 2; % if so: just double
    
    %make sure guess isn't extrapolating beyond previous min guesses.
    %if pointExtrapolation < RatioPGuessMin
    %   pointExtrapolation = (RatioPGuessMin + PathRatioPos(i)) / 2; %if so, just go half way toward that edge 
        
    % is step more than half value
    elseif pointExtrapolation < PathRatioPos(i-1) / 2 
        pointExtrapolation = PathRatioPos(i-1) / 2; % if so: just half
    end
    
    for k = 1:maxSolveTries   %break out of this loop if solve is sucessfull. 
    
        %move point to next location
        PathRatioPos(i) = pointExtrapolation;
        
        %calculate the function for RatioNeg dependence on RatioPos that would result with given targets.
        %this equation enforces that the net charge on the microdroplet is 0. Modify for charged droplets 
        RatioNeg = -(TargetNS - TargetS(j)) / (TargetS(j) - TargetP * (1/PathRatioPos(i) + 1) + TargetNS - TargetS(j));

        %calculate the function for innerConc. uses the average H value of the droplet
        %this equation enforces that the net charge on the microdroplet is 0 to get a valid inner conc value. Modify for charged droplets.
        innerConc = (TargetS_kions + sqrt(TargetS_kions^2 + 4 * (1 + PathRatioPos(i)) * K_w * (1 + RatioNeg)))/ (2*(1+PathRatioPos(i))*DropletVolumeUm3);        
 
        %calculate point
        [Sxint, xint, sucessfulSolve,rGuess] = midfit_PTarget_V3(TargetP_kions,PathRatioPos(i),RatioNeg,DropParam,innerConc,rGuess,showSolutionGraph);

        %check solution OK
        if sucessfulSolve == true
                %get droplet values
                %note final set of variables gets overwritten each loop and is only used for the final
                 %resulting point
                [MolesNegativeIons(i), MolesPositiveIons(i), MolesSurfaceIons(i), ~, ~, ~]...
                 = DropletValues([params, PathRatioPos(i), RatioNeg],Sxint, xint, sucessfulSolve);
                                          
                if isnan(MolesNegativeIons(i)) || isnan(MolesPositiveIons(i)) || isnan(MolesSurfaceIons(i))
                    %bad solution, take smaller step (currently 1/10th of previous step)
                    if verbose,fprintf('(k = %i) ',k),end
                    pointExtrapolation = .1 * (pointExtrapolation(1) - PathRatioPos(i-1)) + PathRatioPos(i-1);
                else
                    break
                end
                
        else %ie bad solution, take smaller step (currently 1/10th of previous step)
                if verbose,fprintf('(k = %i) ',k),end
                pointExtrapolation = .1 * (pointExtrapolation(1) - PathRatioPos(i-1)) + PathRatioPos(i-1);
                rGuess = 48/50 * R; %reset rGuess 
        end
        
        if(k == maxSolveTries) 
            fprintf('WARNING: k hit max loops, will end solve of this point. Location j,i: (%i,%i) \n in function for TargetP_kions %g, TargetNS_kions %g, TargetS_kions %g, R %g',...
                    j,i,TargetP_kions, TargetNS_kions, TargetS_kions, R) 
            U = [NaN, NaN;...
                 NaN, NaN;...
                 NaN, NaN];
            storePath = NaN;
            RadialSolutions = NaN;
            return
        end 
        
    end
    

    %check if close enough to correct values or last loop
    if abs(MolesNegativeIons(i) - TargetNS) < precisionNS 
%&& abs(MolesPositiveIons(i) - TargetP) < precisionP) || i == maxloops
        
        if i == maxloops, disp('hit max steps allowed on path, increase step limit'),end
            
        %store i
        StepsToSolve(j) = i;

        %get droplet values
           [FinalMolesNegativeIons(j), FinalMolesPositiveIons(j), FinalMolesSurfaceIons(j), FinalU_elec(j), FinalTotalOH(j), FinalTotalH(j)]...
           = DropletValues([params, PathRatioPos(i), RatioNeg],Sxint, xint, sucessfulSolve);

         %back into J
           FinalU_elec(j) = FinalU_elec(j) * 10^-12; 

        %Enthalpy of nuetralization for strong acid base is -57.62 kJ/mol at 25C
           FinalU_recOH(j) = FinalTotalOH(j) * 1000 / N_A * 57.62 * 1000; %in J
           FinalU_recH(j)  = FinalTotalH(j) * 1000 / N_A * 57.62 * 1000; %in J
           FinalU(j) = FinalU_elec(j) + FinalU_recOH(j);    
       break
    end
    if manyGraphs
        figure(101)
        plot(1:i,PathRatioPos)
        title(sprintf('solve path for j = %i \n targets P NS and S are: %.2e %.2e %.2e', j, TargetP, TargetNS, TargetS(j)))
        drawnow()

        figure(102)
        plot(1:i,[MolesNegativeIons(1:i); MolesPositiveIons(1:i); MolesSurfaceIons(1:i)])
        title(sprintf('solve path for j = %i', j))
        title(sprintf('solve path for j = %i \n targets P NS and S are: %.2e %.2e %.2e', j, TargetP, TargetNS, TargetS(j)))
        drawnow()
        
        figure(103)
        plot(1:i,MolesNegativeIons(1:i)/TargetNS)
        title(sprintf('solve path for j = %i', j))
        title(sprintf('solve path for j = %i \n targets NS is: %.2e \n ratio NS/targetNS ', j, TargetNS))
        drawnow()
        
        figure(104)
        plot(PathRatioPos,MolesNegativeIons(1:i)/TargetNS,'-o')
        title(sprintf('solved curve for j = %i \n targets NS is: %.2e \n ratio NS/targetNS ', j, TargetNS))
        drawnow()
    end 
end
    %stuff from each j loops that get returned
    storePath{j} = [PathRatioPos; TargetP_kions*ones(1,i); MolesNegativeIons(1:i); MolesPositiveIons(1:i); MolesSurfaceIons(1:i)];
    RadialSolutions{j} = {[params, PathRatioPos(i), RatioNeg],Sxint, xint, sucessfulSolve};   
end

if verbose,fprintf('\n'),end

slope = polyfit(FinalMolesSurfaceIons,FinalU,1);
CorrelationCoef = corrcoef(FinalMolesSurfaceIons,FinalU);
Rsquared_U = CorrelationCoef(1,2)*CorrelationCoef(1,2);
UgradV = slope(1) .* N_A./1000 ./ (F / 1000 * N_A);

slope = polyfit(FinalMolesSurfaceIons,FinalU_elec,1);
CorrelationCoef = corrcoef(FinalMolesSurfaceIons,FinalU_elec);
Rsquared_U_elec = CorrelationCoef(1,2)*CorrelationCoef(1,2);
UgradV_elec = slope(1) .* N_A./1000 ./ (F / 1000 * N_A);

slope = polyfit(FinalMolesSurfaceIons,FinalU_recOH,1);
CorrelationCoef = corrcoef(FinalMolesSurfaceIons,FinalU_recOH);
Rsquared_U_recOH = CorrelationCoef(1,2)*CorrelationCoef(1,2);
UgradV_H2O = slope(1) .* N_A./1000 ./ (F / 1000 * N_A);

U = [UgradV,      Rsquared_U;...
     UgradV_elec, Rsquared_U_elec;...
     UgradV_H2O,  Rsquared_U_recOH];

if manyGraphs
    figure()
    plot(FinalMolesSurfaceIons,FinalU)
    figure()
    plot(FinalMolesSurfaceIons,FinalU_elec)
    figure()
    plot(FinalMolesSurfaceIons,FinalU_recOH)
end

%DropletINFO
ions_P_molar = FinalMolesPositiveIons(1) * convert_kIons2M;
ions_N_molar = (FinalMolesNegativeIons(1) - mean(FinalMolesSurfaceIons)) * convert_kIons2M;
ions_surface_Charge_Density = mean(FinalMolesSurfaceIons) * convert_kIons2IonsPerNM2;

fprintf('Droplet Radious: %i \x03BCm \n', R)
fprintf('Ion Molarity: P = %.3e M \n              N = %.3e M \n', ions_P_molar, ions_N_molar)
fprintf('Surface Charge Density: %.3g ions per nm^2 \n \n', ions_surface_Charge_Density)

fprintf('Voltage:            %#.3g mV (R^2 = %f) \n', UgradV * 1000, Rsquared_U)
fprintf('Electrostatics      %#.3g mV (R^2 = %f) \n', UgradV_elec * 1000, Rsquared_U_elec)
fprintf('Water Recombination %#.3g mV (R^2 = %f) \n', UgradV_H2O * 1000, Rsquared_U_recOH)


%get time of run
%toc()
end    

