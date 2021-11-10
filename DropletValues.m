% returns unit dependtent on params units. Typically for this code params will be in units of kions and um.

function [MolesNegativeIons, MolesPositiveIons, MolesSurface, U_elec, TotalOH, TotalH] = DropletValues(params,Sxint, xint, sucessfulSolve)
        if sucessfulSolve
            ThermalVoltage25 = params(1);
            K_w = params(2);
            F = params(3);
            RatioPos = params(4);
            RatioNeg = params(5);

            %find center-to-wall potential
            E = ThermalVoltage25*Sxint(2,:)./Sxint(1,:);
            V = -cumtrapz(xint,E);

            %find ion distributions
            H = Sxint(1,:);
            OH = K_w ./ Sxint(1,:);
            P = RatioPos * H;
            N = RatioNeg * OH;

            %find total ions of each
            HShell = 4 * pi * xint .* xint .* H ; 
            TotalH = trapz(xint,HShell);
            OHShell = 4 * pi * xint .* xint .* OH ; 
            TotalOH = trapz(xint,OHShell);
            TotalP = RatioPos * TotalH;
            TotalN = RatioNeg * TotalOH;

            %find rho
            rho = F * (H + P - OH - N);

            %find surface charges
            rhoShell = 4 * pi * xint .* xint .* rho ; 
            surfaceCharges = -1 * trapz(xint,rhoShell);
            TotalSurface = surfaceCharges / F * -1; 

            %find Energy of System inside droplet 
            U_inShell = 4 * pi * xint .* xint .* rho .* V;
            U_in = 1/2 * trapz(xint,U_inShell);

            %find energy of outer shell
            U_out = 1/2 * surfaceCharges * V(end);

            %find total internal energy
            U_elec = U_in + U_out;
                     
            MolesNegativeIons = TotalN + TotalSurface;
            MolesPositiveIons = TotalP;
            MolesSurface = TotalSurface;
            
        else 
            MolesNegativeIons = NaN;
            MolesPositiveIons = NaN;
            MolesSurface = NaN;
            U_elec = NaN;
            TotalOH = NaN;
            TotalH = NaN;
        end
end