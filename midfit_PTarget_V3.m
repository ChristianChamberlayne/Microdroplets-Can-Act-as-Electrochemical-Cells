%V3: uses slope instead of innerConc

function [Sxint, xint, sucessfulSolve,r] = midfit_PTarget_V3(TargetP,RatioPos,RatioNeg,DropParam,innerConc, rGuess, showSolutionGraph)

arguments
    
    TargetP    (1,1) double;
    RatioPos    (1,1) double;
    RatioNeg    (1,1) double;
    DropParam   (1,4) double;
    innerConc   (1,1) double;
    rGuess      (1,1) double;
    showSolutionGraph (1,1) logical;
   
end

warning('off','MATLAB:ode45:IntegrationTolNotMet') %this will get triggered a bunch for solutions that shoot to inf.

R = DropParam(1); %radious 
K = DropParam(2);
K_w = DropParam(3);
%F = DropParam(4);
options = odeset('RelTol',1e-10,'AbsTol',1e-12,'Refine',10);

slope = 1;
slopeInSolution = false;
tolerance = 10^-6;
maxbadtries = 3;

loops = 50;
rmax = R;
rmin = 0;

if isnan(rGuess), rGuess = 48/50*R; end

path_r = zeros(length(loops));
path_r(1) = rGuess;
path_r(2) = (R - rGuess) * .1 + rGuess;
path_P = zeros(length(loops));

if showSolutionGraph, figure(1001), end


for i = 1:2
    clear InHconc InXint OutHconc OutXint Pconc xint PShell TotalP
   %take first 2 points 
   for k = 1:maxbadtries
    %check slope exists on curve
    if ~slopeInSolution
        [InHconc, InXint, sucessfulSolve] = midfitInner_PTarget_slope(R, slope, RatioPos, RatioNeg, K, K_w);

        if length(InHconc) == 1 %did it fail to solve
            if showSolutionGraph, disp('Midfit_shooting: midfit inner failed to solve values at start, resetting another attempt with higher slope'),end
            slope = slope * 10;
            slopeInSolution = false;
            rmax = R;
            rmin = 0;
            path_r(1) = 1 * R/2;
            path_r(2) = 3 * R/4;
                            
            if k == maxbadtries
                if showSolutionGraph, disp('Midfit_shooting: midfit inner failed to solve intial values after several attempts'),end
                Sxint =  NaN;
                xint =  NaN;
                sucessfulSolve = false;
                r = NaN;
                return
            end
            continue
        end
        
        PShell = 4 * pi * InXint .* InXint .* InHconc(1,:) * RatioPos; 
        TotalP = trapz(InXint,PShell);

        if TotalP >= TargetP
            %requested slope does not exist in solution, reduce slope and resolve
            slope = slope * .05; % does not go down by 10 to avoid an up and down by 10x loop in k.
            slopeInSolution = false;
            continue
        else
            slopeInSolution = true;
        end
    end
        
   clear InHconc InXint OutHconc OutXint Pconc xint PShell TotalP
      
   [InHconc, InXint, sucessfulSolve] = midfitInner_PTarget_slope(path_r(i), slope, RatioPos, RatioNeg, K, K_w); 
   
   %did it fail to solve?   
   if length(InHconc) == 1 || ~sucessfulSolve
            if showSolutionGraph, disp('Midfit_shooting: midfit inner failed to solve values at start, resetting another attempt with higher slope'),end
            slope = slope * 10;
            slopeInSolution = false;
            rmax = R;
            rmin = 0;
            path_r(1) = 1 * R/2;
            path_r(2) = 3 * R/4;
                            
            if k == maxbadtries
                if showSolutionGraph, disp('Midfit_shooting: midfit inner failed to solve intial values after several attempts'),end
                Sxint =  NaN;
                xint =  NaN;
                sucessfulSolve = false;
                r = NaN;
                return
            end  
   else
       break
   end
   end
   
   r_logc = log(InHconc(1,end));
   r_logslope = InHconc(2,end)/InHconc(1,end);
   
   
   sol = ode45(@(x,y)bvp4ode(x,y,RatioPos, RatioNeg,K,K_w),[path_r(i) R],[r_logc r_logslope],options);
   
%did outer part reach the end?
    if sol.x(end) == R
                
%combine inner and outer solutions and find P total
        OutHconc(1,:) = exp(sol.y(1,:));
        OutXint = sol.x;
        Pconc = [InHconc(1,:) OutHconc(1,2:end)]*RatioPos;
        xint =  [InXint(1,:) OutXint(1,2:end)];
        PShell = 4 * pi * xint .* xint .* Pconc; 
        TotalP = trapz(xint,PShell);
        path_P(i) = TotalP(end);

        if abs(path_P(i) - TargetP)/TargetP < tolerance
            OutSxint(1,:) = exp(sol.y(1,:));
            OutSxint(2,:) = OutSxint(1,:) .* sol.y(2,:);
            OutXint = sol.x;
            Sxint =  [InHconc(1,:) OutSxint(1,2:end); InHconc(2,:) OutSxint(2,2:end)];
            xint =  [InXint(1,:) OutXint(1,2:end)];
            r = path_r(i);
            return
                       
        elseif path_P(i) > TargetP
            rmin = path_r(i);
        else 
            rmax = path_r(i);
        end
        
%did not solve out fully
    else
            path_P(i) = NaN;
        if sol.y(1,end) > 0
            rmin = path_r(i);
        else 
            rmax = path_r(i);
        end
    end
end

if isnan(path_P(2)) && isnan(path_P(1)) %both points don't exist 
    if sol.y(1,end) > 0
            path_r(3) = (rmax + path_r(2)) / 2;
    else 
            path_r(3) = (rmin + path_r(2)) / 2;
    end
    
elseif isnan(path_P(1)) %just second point exists
    if  path_P(2) > TargetP
            path_r(3) = (rmax + path_r(2)) / 2;
    else 
            path_r(3) = (rmin + path_r(2)) / 2;
    end
    
elseif isnan(path_P(2)) %just first point exists
    if  path_P(1) > TargetP
            path_r(3) = (rmax + path_r(1)) / 2;
    else 
            path_r(3) = (rmin + path_r(1)) / 2;
    end
    
else % both points exist, do extrapolation.
    
    %linear extrapolation from last two point's values to next guess
    path_r(3) = (TargetP - path_P(2))*(path_r(1) - path_r(2))/(path_P(1) - path_P(2)) + path_r(2);

    %max sure guess isn't extrapolating to outside the droplet or beyond previous guesses. 
    if path_r(3) > rmax 
        %if it does, just go half way toward that edge 
        path_r(3) = (rmax + path_r(2)) / 2;
    elseif path_r(3) < rmin 
        %if it does, just go half way toward that edge 
        path_r(i+1) = (rmin + path_r(2)) / 2;
    end
end


for i = 3:loops
   clear InHconc InXint OutHconc OutXint Pconc xint PShell TotalP
   for k = 1:maxbadtries
    [InHconc, InXint, sucessfulSolve] = midfitInner_PTarget_slope(path_r(i), slope,RatioPos, RatioNeg, K, K_w); 
   
   %did it fail to solve?   
   if length(InHconc) == 1 || ~sucessfulSolve
            if showSolutionGraph, disp('Midfit_shooting: midfit inner failed to solve values, resetting another attempt'),end
            path_r(i) = (rmax - rmin) * k/(maxbadtries+1) + rmin;
            
            if k == maxbadtries
                if showSolutionGraph, disp('Midfit_shooting: midfit inner failed to solve values after many attempts'),end
                Sxint =  NaN;
                xint =  NaN;
                sucessfulSolve = false;
                r = NaN;
                return;
            end        
   else
       break
   end       
   end
   
   r_logc = log(InHconc(1,end));
   r_logslope = InHconc(2,end)/InHconc(1,end);
  
   sol = ode45(@(x,y)bvp4ode(x,y,RatioPos, RatioNeg,K,K_w),[path_r(i) R],[r_logc r_logslope],options);
   
%did outer part reach the end?
    if sol.x(end) == R 
        
%combine inner and outer solutions and find P total
        OutHconc(1,:) = exp(sol.y(1,:));
        OutXint = sol.x;
        Pconc = [InHconc(1,:) OutHconc(1,2:end)]*RatioPos;
        xint =  [InXint(1,:) OutXint(1,2:end)];
        PShell = 4 * pi * xint .* xint .* Pconc ; 
        TotalP = trapz(xint,PShell);
        path_P(i) = TotalP(end);
                
%Determine if solution is within tolerance
        if abs(path_P(i) - TargetP)/TargetP < tolerance
            OutSxint(1,:) = exp(sol.y(1,:));
            OutSxint(2,:) = OutSxint(1,:) .* sol.y(2,:);
            OutXint = sol.x;
            Sxint =  [InHconc(1,:) OutSxint(1,2:end); InHconc(2,:) OutSxint(2,2:end)];
            xint =  [InXint(1,:) OutXint(1,2:end)];
            r = path_r(i);
            %fprintf('(%i) ',length(OutXint))
            return

%are there 2 points to extrapolate from? If not, step halfway toward an edge   
        elseif isnan(path_P(i-1)) 
            if path_P(i) > TargetP
                rmin = path_r(i);
                path_r(i+1) = (rmax + path_r(i)) / 2;
            else 
                rmax = path_r(i);
                path_r(i+1) = (rmin + path_r(i)) / 2;
            end

% is above target?
        elseif path_P(i) > TargetP
            rmin = path_r(i);
            %linear extrapolation from last two point's values
            path_r(i+1) = (TargetP - path_P(i))*(path_r(i-1) - path_r(i))/(path_P(i-1) - path_P(i)) + path_r(i);
            %max sure guess isn't extrapolating to outside the droplet or beyond previous guesses. 
            if  path_r(i+1) > rmax 
                %if it does, just go half way toward that edge 
                path_r(i+1) = (rmax + path_r(i)) / 2;
            end
            
% else is below target
        else 
            rmax = path_r(i);
            %linear extrapolation from last two point's values
            path_r(i+1) = (TargetP - path_P(i))*(path_r(i-1) - path_r(i))/(path_P(i-1) - path_P(i)) + path_r(i);
            %max sure guess isn't extropolating to negative radious or beyond previous guesses. 
            if  path_r(i+1) < rmin 
                %if it does, just go half way toward that edge 
                path_r(i+1) = (rmin + path_r(i)) / 2;
            end
        end
        
        
%did not solve out fully
    else
        path_P(i) = NaN;
        if sol.y(1,end) > 0
            rmin = path_r(i);
            path_r(i+1) = (rmax + path_r(i)) / 2;
        else 
            rmax = path_r(i);
            path_r(i+1) = (rmin + path_r(i)) / 2;
        end
    end    
    
    
    if showSolutionGraph
        title('targetP search: H conc value')
        GraphYOut = exp(sol.y(1,:));
        GraphXOut = sol.x;
        GraphY =  [InHconc(1,:) GraphYOut(2:end)];
        GraphX =  [InXint GraphXOut(2:end)];
        figure(1001)
        hold off
        semilogy(GraphX,GraphY)
        hold on
        semilogy(GraphXOut,GraphYOut)

        figure(1002)
        title('targetP search: P value')
        semilogy(1:i, path_P,'-o')
        %drawnow()
        
        figure(1003)
        title('targetP search: r value ')
        
        trail = 5;
        if i < trail, trail = i-1; end
        semilogy(path_r(end-trail:i),'-o')
        %drawnow()
        
        
        %assumes kion um unit set, only effects this graph if not corrected
        convertC = 6.022140900000000e+05;
               
        GraphYOut = exp(sol.y(1,:));
        GraphXOut = sol.x;
        GraphY =  [InHconc(1,:) GraphYOut(2:end)];
        GraphX =  [InXint GraphXOut(2:end)];
        figure(1004)
        hold off
        semilogy(GraphX,GraphY * RatioPos / convertC)
        hold on
        semilogy(GraphXOut,GraphYOut * RatioPos / convertC)
        title('targetP search: P conc value')
        drawnow()
    end
    
end

%did it not solve in max loops trys? 
Sxint =  NaN;
xint =  NaN;
sucessfulSolve = false;
r = NaN;
if showSolutionGraph, disp('Midfit_shooting: midfit inner failed to solve within max tries'),end
return

end


function dydx = bvp4ode(x,y,RatioPos, RatioNeg,K,K_w) % equation being solved 

%this is the Log(c) equation 
%No seperate singular term!!!! This is the full thing
dydx = [y(2) 
        K * ((1 + RatioPos) * exp(y(1)) - (1 + RatioNeg) * K_w / (exp(y(1)))) - 2 * y(2)./x];
       
end