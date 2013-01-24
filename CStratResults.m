classdef CStratResults < handle
    %CSTRATRESULTS Results accumulator for CStrat
    
    properties
        % Scalar
        dBtoNp = 8.685889638; % = 1 Np / 1 dB
        k
        nx
        T
        LogMode
        
        % Vector
        nArgMax
        ArgMax
        
        % Matrix
        Pest
        Pextr
        lnL
        P0
        
        % Cells
        Xest
        Xextr
        X0
    end
    
    methods
        function SFRes = CStratResults(K, nx, T)
            SFRes.Pest = nan(1,1);
            SFRes.Pextr = nan(1,1);
            SFRes.Xest = cell(1, nx);
            SFRes.Xextr = cell(1, nx);
            SFRes.k = 0;
            SFRes.nx = nx;
            SFRes.T = T;
            for n = 1:SFRes.nx
                SFRes.nArgMax{n} = nan(1, K);
                SFRes.ArgMax{n} = nan(1, K);
            end

        end
        
        function takeResults(SFRes, SF)
            SFRes.k = SFRes.k + 1;
            ind_argmax = SF.getnArgMaxPest();
            for n = 1:SFRes.nx
                SFRes.Xest{n} = SF.Xest{n};
                SFRes.Xextr{n} = SF.Xextr{n};
                SFRes.nArgMax{n}(SFRes.k) = ind_argmax(n);
                SFRes.ArgMax{n}(SFRes.k) = SF.Xest{n}(ind_argmax(n));
            end
            SFRes.Pest = SF.Pest;
            SFRes.Pextr = SF.Pextr;
        end
        
        function setLnL(SFRes, lnL)
            SFRes.lnL = lnL;
        end
        
        function setP0(SFRes, Xmin, Xmax, dX, P0)
            for n = 1:SFRes.nx
                SFRes.X0{n} = Xmin(n):dX(n):Xmax(n);
            end
            SFRes.P0 = P0;
        end
        
        function setLogMode(SF, LogMode)
            SF.LogMode = LogMode;
        end
    end
    
end

