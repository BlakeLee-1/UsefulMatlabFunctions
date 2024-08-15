function output = RtoT_cal_InputRT(R,T, option, order,interpopt,modelopt)

% This function processes files which contain datacells


%% Feed this function:
% --------------------------------------------------
% fileroot : where to look for files
% filenames : a cell of filenames to use in the fit
% tempCellName : name of the datacell entry i.e. (BathTemp etc)corresponding to T
% resCellName : name of the datacell entry corresponding to R
% option : what to do with the data (see below)
% order : order of the log-space polynomial fit
% interpopt : Allows for interpolating input dataset so that points are
% distributed evenly and certain regions don't get overfit.

%% Load and Sort all the relevant data


%% fit model (cheby)
switch modelopt
    case 'Polynomial'
        spc = 2; inc = (10^(-spc));
        Xi = round(min(log(R)),spc);
        Xf = round(max(log(R)),spc);
        X = [Xi+inc:inc:Xf-inc];
        Y = Interp1NonUnique(log(R),log(T),X);
        
        p = polyfit(X,Y,order);
        f = fittype('a*exp(polyval(p,log(x)))');
        pfit = cfit(f,1,p);
        output=pfit;
        c=pfit;
    case 'Chebyshev'
        spc = 2; inc = (10^(-spc));
        ZL = round(min(log10(R)),spc)-inc;
        ZU = round(max(log10(R)),spc)+inc;
        letter = {'a','b','c','d','e','f','g','h','k','l','m','n','o','p','q','r','s','t','u','v','w','z'};
        
        Xstr =  ['(2*log10(x)-' num2str(ZL + ZU) ')/',...
            '(' num2str(ZU - ZL) ')'];
        Tstr = [];
        for j = 0:order
            Tstr = [Tstr '+' letter{j+1} '*' 'chebyshevT(' num2str(j) ',' Xstr ')'];
        end
        ft = fittype(Tstr);
        [xData, yData] = prepareCurveData( R , T );
        if strcmp(interpopt,'Yes')
            interpRs=linspace(min(R),max(R),3000);
            interpT=Interp1NonUnique(R,T,interpRs);
            [xData, yData] = prepareCurveData( interpRs , interpT );
        end
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        c = fit( xData, yData, ft, opts );
        % chebfit = cfit(ft,coeffvalues(c));
        output=c;
    case 'Chebyshev Feed LogRes'
        spc = 2; inc = (10^(-spc));
        ZL = round(min(log10(R)),spc)-inc;
        ZU = round(max(log10(R)),spc)+inc;
        letter = {'a','b','c','d','e','f','g','h','k','l','m','n','o','p','q','r','s','t','u','v','w','z'};
        
        Xstr =  ['(2*x-' num2str(ZL + ZU) ')/',...
            '(' num2str(ZU - ZL) ')'];
        Tstr = [];
        for j = 0:order
            Tstr = [Tstr '+' letter{j+1} '*' 'chebyshevT(' num2str(j) ',' Xstr ')'];
        end
        ft = fittype(Tstr);
        [xData, yData] = prepareCurveData( log10(R) , T );
        if strcmp(interpopt,'Yes')
            interpRs=linspace(min(R),max(R),3000);
            interpT=Interp1NonUnique(R,T,interpRs);
            [xData, yData] = prepareCurveData( interpRs , interpT );
        end
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        c = fit( xData, yData, ft, opts );
        % chebfit = cfit(ft,coeffvalues(c));
        output=c;
        
    case 'Smoothing Spline'
        [xData, yData] = prepareCurveData( R , T );
        ft = fittype( 'smoothingspline' );
        opts = fitoptions( 'Method', 'SmoothingSpline' );
        opts.SmoothingParam = 0.837805372983563;
        if strcmp(interpopt,'Yes')
            interpRs=linspace(min(R),max(R),3000);
            interpT=Interp1NonUnique(R,T,interpRs);
            [xData, yData] = prepareCurveData( interpRs , interpT );
        end
        % Fit model to data.
        c = fit( xData, yData, ft, opts );
        output = c;
    otherwise 
        display(sprintf('Please select a valid option:\n -Polynomial-\n -Chebyshev- \n -Chebyshev Feed LogRes-\n -Smoothing Spline-')); 
end




%% Processing
switch option
    
    case 'plotRvT'
        
        h1 = figure;
        hold on;
        subplot(2,1,1); hold on;
        plot(R,T,'gsq');
        %plot(R,pfit(R),'r--');
        plot(R,c(R),'r--');
        ylabel('T_{Samp}'); xlabel('R_{CX}');
        legend({'Raw Data','ChebyFit'});
        subplot(2,1,2); hold on;
        try
            plot(T,T-c(R)); title('Residuals');
        catch
            plot(T,T-c(R)'); title('Residuals');
        end
        xlabel('T [K]');
        output = c;
        
    case 'fitRvT'
        
        h1 = figure;
        
        subplot(2,1,1); hold on;
        plot(log(R),log(T),'sq','Color','b'); hold on
        title('log(R_{CX}) v log(T_{Samp})')
        ylabel('log(T_{Samp})'); xlabel('log(R_{CXK})');
        plot(X,Y,'g.')
        plot(X,polyval(p,X),'r--');
        legend({'Raw Data','Interpolation','Polyfit'})
        
        subplot(2,1,2)
        plot(log(R),(T-pfit(R)')./T,'k+')
        title('Fit Model Residuals')
        xlabel('log(R_{Samp})'); ylabel('\Delta T/T');
        
        output = c;
        
    case 'getfit'
        
        output = c;
        
    case 'plotSens'
        h1 = figure;
        hold on
        plot(T,abs(1./differentiate(chfit,R)),'b');
        ylabel('|dR/dT| [\Omega/K]'); xlabel('T [K]');
        
    otherwise
        
        display('Invalid option: See (RtoT_cal.m) for ref.')
        
end

end