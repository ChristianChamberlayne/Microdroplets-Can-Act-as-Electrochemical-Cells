%uses x^5 guess instead of line to speed up the solve



function [Sxint, xint, sucessfulSolve] = midfitInner_PTarget_slope(R, slope,RatioPos, RatioNeg, K, K_w)

arguments
    R        (1,1) double; 
    slope    (1,1) double;
    RatioPos (1,1) double;
    RatioNeg (1,1) double;
    K        (1,1) double;
    K_w      (1,1) double;
end


logEdgeC = log(slope);

%makes intial guess for differencial solver
%this guess is x^5 shaped connecting the BC
xmesh = linspace(0,R,1000);
guessSolve = bvpinit(xmesh, @(x)guess(x,slope,R));

%makes S matrix for differencial solver
S = [0 0; 0 -2];
options = bvpset('SingularTerm',S,'NMax',50000,'RelTol',1e-4);
try
%solves equation
sol = bvp4c(@(x,y)node(x,y, RatioPos, RatioNeg, K, K_w), @(x,y)BC(x,y,slope), guessSolve, options);
catch
   Sxint = NaN;
   xint = NaN;
   sucessfulSolve = false;
   return
end
%undo Log equation results

Sxint(1,:) = exp(sol.y(1,:));
Sxint(2,:) = Sxint(1,:) .* sol.y(2,:);
xint = sol.x;
sucessfulSolve = true;

%check for solution errors that come up if the input numbers give too small a doublelayer  
    for point = Sxint(1,:)
        if point < 0
            %disp('Error: Inner function: Negative concentration, parameters out of solvable area, will display questionable results')
            sucessfulSolve = false;
            break
        end
    end

    
    for point = Sxint(2,:)./Sxint(1,:)
        if point < 0
            %disp('Error: Inner function: Negative E, parameters out of solvable area, will display questionable results')
            sucessfulSolve = false;
            break
        end
    end
end

function dydx = node(x,y,RatioPos, RatioNeg,K,K_w) % equation being solved 

%this is the Log(c) equation
dydx = [y(2) 
        K * ((1 + RatioPos) * exp(y(1)) - (1 + RatioNeg) * K_w / (exp(y(1)))) ];
       
end
%-------------------------------------------
function res = BC(ya,yb,slope) % boundary conditions

res = [ya(2)
       yb(2) - slope];
end
%-------------------------------------------
function initY = guess(x,slope,R) % initial guess for solution
%this is x^5 shaped
initY = [slope * (x/R)^5 * R / 5
         slope * (x/R)^4];
end