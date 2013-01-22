classdef CStrat < handle
    %CSTRAT Numerical solver of Stratonovich equations
    
    properties
        % Scalars
        nx 
        T
        Dxxm1
        sqrtDxx
        
        % Vectors
        dX
        Xmin
        Xmax
        lenX
        
        Xestmin
        Xestmax
        jXestmin
        jXestmax
        jXestmin_old
        jXestmax_old        
        lenXest
        
        Xextrmin
        Xextrmax
        jXextrmin
        jXextrmax
        lenXextr        
        
        LastIndexShift
        
        % Matrix
        Pest
        Pextr
        ExpInExtr
        
        % Cell Arrays
        X
        Xest
        Xextr 
    end
    
    methods
        function SF = CStrat(dX, Xmin, Xmax, P0)
            SF.nx = length(dX);
            SF.dX = dX;
            SF.Xmin = Xmin;
            SF.Xmax = Xmax;
            SF.Pest = P0;
            SF.Pextr = P0*0;
        end
        
        function SF = setAprioriParameters(SF, T, Dxx)
            SF.T = T;
            SF.Dxxm1 = 1/Dxx;
            SF.sqrtDxx = sqrt(Dxx);
            SF.resizeX();
        end
        
        function resizeX(SF)
            for j = 1:SF.nx 
                SF.X{j} = SF.Xmin(j):SF.dX(j):SF.Xmax(j);
                SF.lenX(j) = length(SF.X{j});
            end
            SF.resizeXest(SF.Xmin, SF.Xmax);
            
            % exp(-0.5 * (x(end) - x(end))^2 / D )
            SF.ExpInExtr = nan(SF.lenX(SF.nx), SF.lenX(SF.nx));
            for i = 1:SF.lenX(SF.nx)
                for j = 1:SF.lenX(SF.nx)
                    if j>i
                        SF.ExpInExtr(i, j) = exp( -0.5*(SF.X{SF.nx}(i) -  SF.X{SF.nx}(j)).^2*SF.Dxxm1 );
                    elseif i == j
                        SF.ExpInExtr(i, j)  = 1;
                    else
                        SF.ExpInExtr(i, j) = SF.ExpInExtr(j, i);
                    end
                end
            end

            % Indexes for X_k
            mnoj = 1;
            for n = 2:SF.nx
                mnoj = mnoj * SF.lenX(n-1);
            end
            SF.LastIndexShift = (0:(SF.lenX(SF.nx)-1))*mnoj;
            
        end
        
        function resizeXest(SF, Xmin, Xmax)
            
            SF.jXestmin_old = SF.jXestmin;
            SF.jXestmax_old = SF.jXestmax;
            
            SF.jXestmin = floor(Xmin./SF.dX);
            SF.Xestmin = Xmin.*SF.dX;
            
            SF.jXestmax = ceil(Xmax./SF.dX);
            SF.Xestmax = Xmax.*SF.dX;
            
            for n = 1:SF.nx
                SF.Xestmin = SF.jXestmin*SF.dX(n);
                SF.Xestmax = SF.jXestmax*SF.dX(n);
                SF.Xest{n} = SF.Xestmin(n):SF.dX(n):SF.Xestmax(n);
                SF.lenXest(n) = SF.jXestmax(n) - SF.jXestmin(n) + 1;
            end
        end
        
        function shiftPest(SF)
            s=sprintf('Pest_new = zeros(SF.lenXest(1)');
            for n = 2:SF.nx
                s = [s ', ' double2str(SF.lenXest(n))]; 
            end
            s = [s ');'];
            eval(s);
        end
        
        function resizeXextr(SF)
            
            for n = (1:SF.nx-1)
                SF.Xextrmax(n) = SF.Xestmax(n-1) + SF.Xestmax(n-1)*SF.T;
            end
            SF.Xextrmax(SF.nx) = SF.Xestmax(SF.nx) + 5*SF.sqrtDxx;
            SF.Xextrmin(SF.nx) = SF.Xestmin(SF.nx) - 5*SF.sqrtDxx;
            
            SF.jXestmin = floor(SF.Xextrmin./SF.dX);
            SF.Xestmin = SF.Xextrmin.*SF.dX;
            
            SF.jXestmax = ceil(SF.Xextrmax./SF.dX);
            SF.Xestmax = SF.Xextrmax.*SF.dX;
            
            for n = 1:SF.nx
                SF.Xextrmin = SF.jXextrmin*SF.dX(n);
                SF.Xextrmax = SF.jXextrmax*SF.dX(n);
                SF.Xextr{n} = SF.Xextrmin(n):SF.dX(n):SF.Xextrmax(n);
                SF.lenXextr(n) = SF.jXextrmax(n) - SF.jXextrmin(n) + 1;
            end
        end        

        
        function Extrapolate(SF)
            
            SF.Pextr = SF.Pextr*0;
            Xkm1 = nan(SF.nx, 1);
            i = nan(1, SF.nx-1);
            for n = 1:(SF.nx-1)
                i(n) = 1;
            end
            linind = 1; % linear index for X_k(1:end-1)
            not_end_of_int = 1;
            while not_end_of_int

                inde_i = linind + SF.LastIndexShift;
                
                tin = tic;
                for jnx = 1:SF.lenX(SF.nx) % index for X_{k-1}(end)
                    
                    Xkm1(SF.nx) = SF.X{SF.nx}(jnx);
                    for n = (SF.nx-1):-1:1
                        Xkm1(n) = SF.X{n}(i(n)) - Xkm1(n+1)*SF.T; % X_{k-1}(n) in pest = X_{k}(n) - X_{k-1}(n-1)T => .... cycle
                    end
                    
                    if SF.IsInPestArea(Xkm1)
                        pest = SF.getPDValue(SF.Pest, Xkm1);
                        
                        if pest > 0
                            SF.Pextr(inde_i) = SF.Pextr(inde_i) + pest*SF.ExpInExtr(jnx, 1:SF.lenX(SF.nx));
                        end
                    end
                    
                end
                tout = toc(tin);
               
                % Next X_{k} 
                i(1) = i(1) + 1;
                for n = 1:(SF.nx-1)
                    if i(n) <= SF.lenX(n)
                        break;
                    else
                        if n < (SF.nx-1)
                            i(n) = 1;
                            i(n + 1) = i(n + 1) + 1;
                        else
                            i(n) = i(n) - 1;
                            not_end_of_int = 0; % all done
                        end
                    end
                end
                linind = SF.getLinIndexLastButOne(i);
            end
           SF.Pextr = SF.normPD(SF.Pextr);
        end
        
        function Observe(SF, lnL)
            ln_pest = log(SF.separatePD(SF.Pextr, 120)) + lnL;
            pest = SF.separatePD_ln(ln_pest, 120);
            SF.Pest = SF.normPD(pest);
        end
        
        function Ps = separatePD(SF, P, dyn_dB)
            maxP = max(P);
            for i = 2:SF.nx
                maxP = max(maxP);
            end
            minP = maxP*exp(-dyn_dB);
            Ps = (P > minP).*P  + (P <= minP)*minP;
            Ps = Ps / minP; % min -> 1
        end
        
        function LnPs = separatePD_ln(SF, LnP, dyn_dB)
            MaxLnP = max(LnP);
            for i = 2:SF.nx
                MaxLnP = max(MaxLnP);
            end
            MinLnPs = MaxLnP - dyn_dB;
            LnPs = (LnP > MinLnPs ).*LnP  + (LnP <= MinLnPs )*MinLnPs ;
            LnPs = LnPs - MinLnPs;
            LnPs = exp(LnPs).*(LnPs > MinLnPs); % *0, if  arg = argmax - dB_est
        end

        function Pnorm = normPD(SF, P)
            inte = sum(P);
            dXX = SF.dX(1);
            for i = 2:SF.nx
                inte = sum(inte);
                dXX = dXX * SF.dX(i);
            end
            Pnorm = P / inte * dXX;
        end

        % P(Xp) - ?
        function p = getPDValue(SF, P, Xp)
            nxp = SF.lenX;
            Grad = zeros(1, SF.nx);
            stepX = zeros(SF.nx, 1);
            for n = 1:SF.nx
                difx = Xp(n) - SF.X{n};
                [absstep nxp(n)] = min(abs(difx));
                if (difx(nxp(n)) < 0) && (nxp(n) > 1)
                    nxp(n) = nxp(n) - 1;
                end
            end
%             for n = 1:SF.nx
%                 step_ind = zeros(1, SF.nx);
%                 step_ind(n) = 1;
%                 if nxp(n) < SF.lenX(n)
%                     Grad(n) = (P(SF.getLinIndex(nxp + step_ind)) - P(SF.getLinIndex(nxp))) / SF.dX(n);
%                     stepX(n) = Xp(n) - SF.X{n}(nxp(n));
%                 else
%                     Grad(n) = 0;
%                     stepX(n) = 0;
%                 end
%             end
            p = P(SF.getLinIndex(nxp)) + Grad*stepX;
        end
        
        % a(i1, i2, i3) -> a(i1 + (i2-1)*() + ...)
        function i = getLinIndex(SF, np)
            i = np(1);
            mnoj = 1;
            for n = 2:SF.nx
                mnoj = mnoj * SF.lenX(n-1);
                i = i + (np(n)-1)*mnoj;
            end
        end
        
        function i = getLinIndexLastButOne(SF, np)
            i = np(1);
            mnoj = 1;
            for n = 2:(SF.nx-1)
                mnoj = mnoj * SF.lenX(n-1);
                i = i + (np(n)-1)*mnoj;
            end
        end        
        
        function isIn = IsInPestArea(SF, Xp)
            isIn = 1;
            for n = 1:SF.nx
                if (Xp(n) < SF.Xestmin(n)) || (Xp(n) > SF.Xestmax(n))
                    isIn = 0;
                    return;
                end
            end
        end
        
%         function setPestArea(SF)
%             SF.Pest > 0
%         end
        
        function Resa = getResults(SF)
            Resa.X = SF.X;
            Resa.Pest = SF.Pest;
            Resa.Pextr = SF.Pextr;
        end
    end
    
end

