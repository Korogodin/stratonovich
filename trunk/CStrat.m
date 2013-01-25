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
        
        Xestmin
        Xestmax
        lenXest
        
        Xextrmin
        Xextrmax
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
        Pest_j
    end
    
    methods
        function SF = CStrat(T, Dxx, Xmin, Xmax, dX, P0)
            SF.nx = length(dX);
            SF.dX = dX;
            SF.Pest = P0;
            SF.T = T;
            SF.Dxxm1 = 1/Dxx;
            SF.sqrtDxx = sqrt(Dxx);
            SF.setXest(Xmin, Xmax);
            SF.resizeXextr();            
        end
       
        function setXest(SF, Xmin, Xmax)
            jXestmin = floor(Xmin./SF.dX);
            SF.Xestmin = Xmin.*SF.dX;
            
            jXestmax = ceil(Xmax./SF.dX);
            SF.Xestmax = Xmax.*SF.dX;
            
            SF.Xestmin = jXestmin.*SF.dX;
            SF.Xestmax = jXestmax.*SF.dX;
            for n = 1:SF.nx
                SF.Xest{n} = SF.Xestmin(n):SF.dX(n):SF.Xestmax(n);
            end
            SF.lenXest = jXestmax - jXestmin + 1;
            
            SF.Xextrmax = 0*SF.Xestmax;
            SF.Xextrmin = 0*SF.Xestmin;
        end
        
        function lnPest_cuted = cutXextr2Xest(SF, lnPest, imin, imax)
            
            SF.Xestmin = [SF.Xextr{1}(imin(1)); SF.Xextr{2}(imin(2))];
            SF.Xestmax = [SF.Xextr{1}(imax(1)); SF.Xextr{2}(imax(2))];
            
            SF.Xestmin = round(SF.Xestmin./SF.dX).*SF.dX;
            SF.Xestmax = round(SF.Xestmax./SF.dX).*SF.dX;
            
            SF.lenXest = imax - imin + 1;
            
            for n = 1:SF.nx
                SF.Xest{n} = SF.Xestmin(n):SF.dX(n):SF.Xestmax(n);
            end            
            
            s=sprintf('lnPest_cuted = zeros(SF.lenXest(1)');
            for n = 2:SF.nx
                s = [s ', ' num2str(SF.lenXest(n))]; 
            end
            s = [s ');'];
            eval(s);
            
            if SF.nx == 2
                for i = 1:SF.lenXest(1)
                    for j = 1:SF.lenXest(2)
                        lnPest_cuted(i, j) = lnPest_cuted(i, j) + lnPest( imin(1) - 1 + i , imin(2) - 1 + j );
                    end
                end
            end
        end
        
        
        function resizeXextr(SF)
            
            for n = (1:SF.nx-1)
                SF.Xextrmax(n) = SF.Xestmax(n) + SF.Xestmax(n+1)*SF.T;
                SF.Xextrmin(n) = SF.Xestmin(n) + SF.Xestmin(n+1)*SF.T;
            end
            SF.Xextrmax(SF.nx) = SF.Xestmax(SF.nx) + 5*SF.sqrtDxx;
            SF.Xextrmin(SF.nx) = SF.Xestmin(SF.nx) - 5*SF.sqrtDxx;
            
            jXextrmax = ceil(SF.Xextrmax./SF.dX);
            SF.Xextrmax = jXextrmax.*SF.dX;
            jXextrmin = floor(SF.Xextrmin./SF.dX);
            SF.Xextrmin = jXextrmin.*SF.dX;
            
            for n = 1:SF.nx
                SF.Xextr{n} = SF.Xextrmin(n):SF.dX(n):SF.Xextrmax(n);
            end
            SF.lenXextr = jXextrmax - jXextrmin + 1;
            
            s=sprintf('SF.Pextr = zeros(SF.lenXextr(1)');
            for n = 2:SF.nx
                s = [s ', ' num2str(SF.lenXextr(n))]; 
            end
            s = [s ');'];
            eval(s);            
            
            SF.ExpInExtr = nan(SF.lenXextr(SF.nx), SF.lenXest(SF.nx));
            for i = 1:SF.lenXextr(SF.nx)
                for j = 1:SF.lenXest(SF.nx)
                        SF.ExpInExtr(j, i) = exp( -0.5*(SF.Xextr{SF.nx}(i) -  SF.Xest{SF.nx}(j)).^2*SF.Dxxm1 );
                end
            end

            % Indexes for X_k
            mnoj = 1;
            for n = 2:SF.nx
                mnoj = mnoj * SF.lenXextr(n-1);
            end
            SF.LastIndexShift = (0:(SF.lenXextr(SF.nx)-1))*mnoj;            
        end        

        
        function Extrapolate(SF)
            
            SF.resizeXextr();
            Xkm1 = nan(SF.nx, 1);
            i = nan(1, SF.nx-1);
            for n = 1:(SF.nx-1)
                i(n) = 1;
            end
            linind = 1; % linear index for X_k(1:end-1)
            not_end_of_int = 1;
            while not_end_of_int

                inde_i = linind + SF.LastIndexShift;
                
                for jnx = 1:SF.lenXest(SF.nx) % index for X_{k-1}(end)
                    Xkm1(SF.nx) = SF.Xest{SF.nx}(jnx);
                    for n = (SF.nx-1):-1:1
                        Xkm1(n) = SF.Xextr{n}(i(n)) - Xkm1(n+1)*SF.T; % X_{k-1}(n) in pest = X_{k}(n) - X_{k-1}(n-1)T => .... cycle
                    end
                    if SF.IsInPestArea(Xkm1)
                        pest = SF.getPestValue(Xkm1);
                        if pest > 0
                            SF.Pextr(inde_i) = SF.Pextr(inde_i) + pest*SF.ExpInExtr(jnx, 1:SF.lenXextr(SF.nx));
                        end
                    end
                end
              
                % Next X_{k} 
                i(1) = i(1) + 1;
                for n = 1:(SF.nx-1)
                    if i(n) <= SF.lenXextr(n)
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
        
        function nmax = getnArgMaxPest(SF)
            if SF.nx == 2
                [a b] = max(SF.Pest);
                [nul, n2] = max(a);
                n1 = b(n2);
                nmax = [n1, n2];
            else
                disp('Need some function there... getArgMaxPest');
                error; 
            end
        end
        
        function calcPest_j(SF)
            for j = 1:SF.nx
                SF.Pest_j{j} = zeros(SF.lenXest(j), 1);
            end
            if SF.nx == 2
                for j1 = 1:SF.lenXest(1)
                    SF.Pest_j{1}(j1) = sum(SF.Pest(j1, :));
                end
                SF.Pest_j{1} = SF.Pest_j{1} / sum(SF.Pest_j{1}) / SF.dX(1);
                for j2 = 1:SF.lenXest(2)
                    SF.Pest_j{2}(j2) = sum(SF.Pest(:, j2));
                end
                SF.Pest_j{2} = SF.Pest_j{2} / sum(SF.Pest_j{2}) / SF.dX(2);
            else
                disp('Need some function there... calcPest_j');
                error; 
            end
        end
        
        function mX = getMeanPest(SF)
            mX = zeros(SF.nx, 1);
            for n = 1:SF.nx
                mX(n) = SF.Xest{n}*SF.Pest_j{n};
            end
            mX = mX.*SF.dX;
        end        
        
        function Observe(SF, lnL)
            ln_pest = log(SF.separatePD(SF.Pextr, 120)) + lnL;
            pest = SF.separatePD_ln(ln_pest, 10);
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
        
        function pest = separatePD_ln(SF, LnP, dyn_dB)
            MaxLnP = max(LnP);
            for i = 2:SF.nx
                MaxLnP = max(MaxLnP);
            end
            MinLnPs = MaxLnP - dyn_dB;
            LnPs = (LnP > MinLnPs ).*LnP  + (LnP <= MinLnPs )*MinLnPs ;
            thereisPD = (LnP >= MinLnPs);
            LnPs = LnPs - MinLnPs;
            LnPs = SF.setPestArea(LnPs, thereisPD);
%             pest = exp(LnPs).*thereisPD; % *0, if  arg = argmax - dB_est
            pest = exp(LnPs); % *0, if  arg = argmax - dB_est
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
        function p = getPestValue(SF, Xp)
            nxp = SF.lenXest;
            Grad = zeros(1, SF.nx);
            stepX = zeros(SF.nx, 1);
%             tin = tic;
            for n = 1:SF.nx
                difx = Xp(n) - SF.Xest{n};
                [absstep nxp(n)] = min(abs(difx));
                if (difx(nxp(n)) < 0) && (nxp(n) > 1)
                    nxp(n) = nxp(n) - 1;
                end
            end
%             tout1 = toc(tin);
            if SF.nx == 2
                if nxp(1) < SF.lenXest(1)
                    Grad(1) = (SF.Pest(nxp(1) + 1, nxp(2)) - SF.Pest(nxp(1), nxp(2))) / SF.dX(1);
                    stepX(1) = Xp(1) - SF.Xest{1}(nxp(1));
                else
                    Grad(1) = 0;
                    stepX(1) = 0;
                end
                if nxp(2) < SF.lenXest(2)
                    Grad(2) = (SF.Pest(nxp(1), nxp(2) + 1) - SF.Pest(nxp(1), nxp(2))) / SF.dX(2);
                    stepX(2) = Xp(2) - SF.Xest{2}(nxp(2));
                else
                    Grad(2) = 0;
                    stepX(2) = 0;
                end
            else
                for n = 1:SF.nx
                    step_ind = zeros(SF.nx, 1);
                    step_ind(n) = 1;
                    if nxp(n) < SF.lenXest(n)
                        Grad(n) = (SF.Pest(SF.getLinIndexEst(nxp + step_ind)) - SF.Pest(SF.getLinIndexEst(nxp))) / SF.dX(n);
                        stepX(n) = Xp(n) - SF.Xest{n}(nxp(n));
                    else
                        Grad(n) = 0;
                        stepX(n) = 0;
                    end
                end
            end
%             tout2 = toc(tin);
            if SF.nx == 2
                p = SF.Pest(nxp(1), nxp(2)) + Grad*stepX;
            else
                p = SF.Pest(SF.getLinIndexEst(nxp)) + Grad*stepX;
            end
%             tout3 = toc(tin);
        end
        
        % a(i1, i2, i3) -> a(i1 + (i2-1)*() + ...)
        function i = getLinIndexEst(SF, np)
            i = np(1);
            mnoj = 1;
            for n = 2:SF.nx
                mnoj = mnoj * SF.lenXest(n-1);
                i = i + (np(n)-1)*mnoj;
            end
        end
        function i = getLinIndexExtr(SF, np)
            i = np(1);
            mnoj = 1;
            for n = 2:SF.nx
                mnoj = mnoj * SF.lenXextr(n-1);
                i = i + (np(n)-1)*mnoj;
            end
        end
        
        
        function i = getLinIndexLastButOne(SF, np)
            i = np(1);
            mnoj = 1;
            for n = 2:(SF.nx-1)
                mnoj = mnoj * SF.lenXextr(n-1);
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
        
        function cutedLnPs = setPestArea(SF, LnPs, thereisPD)
               
            imin = ones(SF.nx, 1);
            imax = SF.lenXextr;
            if SF.nx == 2
                for n = 1:SF.lenXextr(1)
                    if sum(thereisPD(n, :))
                        imin(1) = n;
                        break;
                    end
                end
                
                for n = 1:SF.lenXextr(2)
                    if sum(thereisPD(:, n))
                        imin(2) = n;
                        break;
                    end
                end
                
                for n = SF.lenXextr(1):-1:1
                    if sum(thereisPD(n, :))
                        imax(1) = n;
                        break;
                    end
                end
                
                for n = SF.lenXextr(2):-1:1
                    if sum(thereisPD(:, n))
                        imax(2) = n;
                        break;
                    end
                end
                
            elseif SF.nx == 3
                
            end
            cutedLnPs = cutXextr2Xest(SF, LnPs, imin, imax);
        end
    end
    
end

