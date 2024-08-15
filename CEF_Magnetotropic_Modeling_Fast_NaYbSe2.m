function output = CEF_Magnetotropic_Modeling_Fast_NaYbSe2(option, params)


%% Fixed angle/Temp ZF resonance frequencies
amp = 9.95;
T_RTM = [4,10,20,40,80,160]; 
w0 = [49.5241, 49.5276, 49.5238, 49.5215, 49.5149, 49.4941, 49.4440, 49.3838];

%% 2nd Deriv coefficients
C1 = [1/12, 	-2/3, 	0,      2/3, 	-1/12];
C2 = [-1/12 	4/3 	-5/2 	4/3 	-1/12];

%% Spin 7/2 matrices, {|7/2,7/2>.,...,|7/2,-7/2>| basis

ID = [1,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0;
    0,0,0,1,0,0,0,0;
    0,0,0,0,1,0,0,0;
    0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,1];

Jz = [ 7/2,   0,   0,   0,   0,   0,   0,   0;
    0, 5/2,   0,   0,   0,   0,   0,   0;
    0,   0, 3/2,   0,   0,   0,   0,   0;
    0,   0,   0, 1/2,   0,   0,   0,   0;
    0,   0,   0,   0,-1/2,   0,   0,   0;
    0,   0,   0,   0,   0,-3/2,   0,   0;
    0,   0,   0,   0,   0,   0,-5/2,   0;
    0,   0,   0,   0,   0,   0,   0,-7/2];

Jp =  [       0,sqrt( 7),       0,       0,       0,       0,       0,       0;
    0,       0,sqrt(12),       0,       0,       0,       0,       0;
    0,       0,       0,sqrt(15),       0,       0,       0,       0;
    0,       0,       0,       0,       4,       0,       0,       0;
    0,       0,       0,       0,       0,sqrt(15),       0,       0;
    0,       0,       0,       0,       0,       0,sqrt(12),       0;
    0,       0,       0,       0,       0,       0,       0,sqrt( 7);
    0,       0,       0,       0,       0,       0,       0,       0];

Jm =  [       0,       0,       0,       0,       0,       0,       0,       0;
    sqrt( 7),       0,       0,       0,       0,       0,       0,       0;
    0,sqrt(12),       0,       0,       0,       0,       0,       0;
    0,       0,sqrt(15),       0,       0,       0,       0,       0;
    0,       0,       0,       4,       0,       0,       0,       0;
    0,       0,       0,       0,sqrt(15),       0,       0,       0;
    0,       0,       0,       0,       0,sqrt(12),       0,       0;
    0,       0,       0,       0,       0,       0,sqrt( 7),       0];

Jx = (Jp+Jm)/2;

Jy = (Jp-Jm)/2i;

%% Steven's Operators
J=7/2; X = J*(J+1); A = Jp*Jp*Jp + Jm*Jm*Jm;

O20 = 3*Jz*Jz - X*ID;
O40 = 35*power(Jz,4) - (30*X - 25)*Jz*Jz + (3*X*X - 6*X)*ID;
O60 = 231*power(Jz,6) - (315*X-735)*power(Jz,4) + (105*X*X - 525*X +294)*power(Jz,2) - (5*X*X*X + 40*X*X -60*X)*ID;
O43 = (1/4)*( (A)*Jz + Jz*(A) );
O63 = (1/4)*( A*(11*power(Jz,3) - (3*X + 59)*Jz ) + (11*power(Jz,3) -(3*X + 59)*Jz)*A );
O66 = (1/2)*(Jp*Jp*Jp*Jp*Jp*Jp + Jm*Jm*Jm*Jm*Jm*Jm);

fac = [1,10,1,1e3,1e2,1e2];

HCEFf = @(Bmn) Bmn(1)*O20/fac(1) + Bmn(2)*O40/fac(2) + Bmn(3)*O43/fac(3)...
    + Bmn(4)*O60/fac(4) + Bmn(5)*O63/fac(5) + Bmn(6)*O66/fac(6);

muB = 5.78828e-2; % [meV/T];
kB  = 8.617e-2  ; % [meV/K];
gJ = 8/7; % L=3, S=1/2, J=7/2 g-lande factor;
A1 = 6*kB/(muB*gJ);
C0 = 2.0416;

E10 = 15.8;
E20 = 24.3;
E30 = 30.5;

%% Curve Functions

    function out = Mab_vH_alt(HCEF,E0,Jxx,Jzz,t,H,dH)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mx = zeros(s);
        mzT0 = zeros(s);
        mzT1 = zeros(s);
        
        for h = [H(1)-2*dH, H(1)-dH, H, H(end)+dH, H(end)+2*dH]
            
            HpZ = HCEF - gJ*muB*h*Jx;
            F1(n) = -(1./B)*log(trace(expm(-B*(HpZ - E0*ID)) ));
            n = n+1;
        end
        
        M1 = -(C1(1)*F1(1:end-4) + C1(2)*F1(2:end-3) + C1(3)*F1(3:end-2) + C1(4)*F1(4:end-1) + C1(5)*F1(5:end))/dH;
        
        %X1 = H + Jxx*A1*mx;
        out = M1;%Interp1NonUnique([0 X1], [0 gJ*muB.*mx], H);
    end

    function out = Mc_vH_alt(HCEF,E0,Jxx,Jzz,t,H,dH)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mx = zeros(s);
        mzT0 = zeros(s);
        mzT1 = zeros(s);
        
        for h = [H(1)-2*dH, H(1)-dH, H, H(end)+dH, H(end)+2*dH]
            
            HpZ = HCEF - gJ*muB*h*Jz;
            F1(n) = -(1./B)*log(trace(expm(-B*(HpZ - E0*ID)) ));
            n = n+1;
        end
        
        M1 = -(C1(1)*F1(1:end-4) + C1(2)*F1(2:end-3) + C1(3)*F1(3:end-2) + C1(4)*F1(4:end-1) + C1(5)*F1(5:end))/dH;
        
        %X1 = H + Jxx*A1*mx;
        out = M1;%Interp1NonUnique([0 X1], [0 gJ*muB.*mx], H);
    end

    function out = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mx = zeros(s);
        mzT0 = zeros(s);
        mzT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jx;
            mx(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            n = n+1;
        end
        
        X1 = H + Jxx*A1*mx;
        out = Interp1NonUnique([0 X1], [0 gJ*muB.*mx], H);
    end

    function out = Mc_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mz = zeros(s);
        mxT0 = zeros(s);
        mxT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jz;
            mz(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            n = n+1;
        end
        
        X1 = H + Jzz*A1*mz;
        out = Interp1NonUnique([0 X1], [0 gJ*muB.*mz], H);
    end

    function out = kab_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mx = zeros(s);
        mzT0 = zeros(s);
        mzT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jx;
            mx(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            mzT0(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            HpZ = HpZ - gJ*muB*hs*Jz;
            mzT1(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            n = n+1;
        end
        
        X1 = H + Jxx*A1*mx; X2 = power(X1,2);
        chiT = (mzT1-mzT0)./(hs + Jzz*A1*(mzT1 - mzT0));
        
        out = Interp1NonUnique([0 X1], [0 gJ*muB*(X1.*mx - X2.*chiT)], H);
    end

    function out = kc_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mz = zeros(s);
        mxT0 = zeros(s);
        mxT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jz;
            mz(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            mxT0(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            HpZ = HpZ - gJ*muB*hs*Jx;
            mxT1(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            n = n+1;
        end
        
        X1 = H + Jzz*A1*mz; X2 = power(X1,2);
        chiT = (mxT1-mxT0)./(hs + Jxx*A1*(mxT1 - mxT0));
        
        out = Interp1NonUnique([0 X1], [0 gJ*muB*(X1.*mz - X2.*chiT)], H);
    end

    function out = XTc_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mx = zeros(s);
        mzT0 = zeros(s);
        mzT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jx;
            mx(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            mzT0(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            HpZ = HpZ - gJ*muB*hs*Jz;
            mzT1(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            n = n+1;
        end
        
        X1 = H + Jxx*A1*mx; 
        chiT = (mzT1-mzT0)./(hs + Jzz*A1*(mzT1 - mzT0));
        
        out = Interp1NonUnique([0 X1], gJ*muB*[chiT(1) chiT], H);
    end

    function out = XTab_vH(HCEF,E0,Jxx,Jzz,t,H)
        
        B = 1./(kB*t); hs = .01;
        n = 1; s = size(H);
        mz = zeros(s);
        mxT0 = zeros(s);
        mxT1 = zeros(s);
        
        for h = H
            
            HpZ = HCEF - gJ*muB*h*Jz;
            mz(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            mxT0(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            HpZ = HpZ - gJ*muB*hs*Jx;
            mxT1(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            
            n = n+1;
        end
        
        X1 = H + Jzz*A1*mz;
        chiT = (mxT1-mxT0)./(hs + Jxx*A1*(mxT1 - mxT0));
        
        out = Interp1NonUnique([0 X1], gJ*muB*[chiT(1) chiT], H);
    end

    function out = Xc_inv_vT(HCEF,E0,Jaa,T)
        
        hs = .01;
        HpZ = HCEF - gJ*muB*hs*Jz;
        
        n = 1; m = zeros(size(T));
        for t1 = T
            B = 1./(kB*t1);
            m(n) = trace( Jz*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            n = n+1;
        end
        
        out = C0*((gJ*muB/kB)*hs./m + 6*Jaa);
        
    end

    function out = Xab_inv_vT(HCEF,E0,Jaa,T)
        
        hs = .01;
        HpZ = HCEF - gJ*muB*hs*Jx;
        
        n = 1; m = zeros(size(T));
        for t1 = T
            B = 1./(kB*t1);
            m(n) = trace( Jx*expm(-B*(HpZ - E0*ID)) )./trace( expm(-B*(HpZ - E0*ID)));
            n = n+1;
        end
        
        out = C0*((gJ*muB/kB)*hs./m + 6*Jaa);
        
    end

    M_params = []; H_single = [];

%% Minimization Function

    function out = CummErr(p0)
        
        
        HCEF1 = HCEFf(p0(1:6));
        Jxx1 = p0(7); Jzz1 = p0(8); chi0 = p0(9);

        [P1,D1] = eig(HCEF1 + Jz*1e-10);

        En = real([D1(1,1) D1(2,2) D1(3,3) D1(4,4) D1(5,5) D1(6,6) D1(7,7) D1(8,8)]);
        Ev1 = sort(En);
        E1 = Ev1(1);

        for n = 1:4
            ind1 = find(En==Ev1(2*(n-1) + 1));
            ind2 = find(En==Ev1(2*n));

            ev1{n} = P1(:,ind1(1));
            ev2{n} = P1(:,ind2(end));
        end
        
        for n3 = 2:8
             out = out + sum(power(kd_int{1}{n3} - kab_vH(HCEF1,E1,Jxx1,Jzz1,T_RTM(n3),Hd_int{1}{n3}), 2));
             out = out + sum(power(kd_int{2}{n3} -  kc_vH(HCEF1,E1,Jxx1,Jzz1,T_RTM(n3),Hd_int{2}{n3}), 2));
         end 


        out = out + 10*(sum(power(Xa1 - XaDatInt,2)) + sum(power(Xc1 - XcDatInt,2))) ;
%

    end

 
%% Load data
temp = load('C:\Users\LeeLabLaptop\Documents\ML_Dec21_RTM\NaYbSe2_RTM.mat');

for nn = 1:6
    Hd_int{1}{nn} = temp.Dat{nn}.H;
    kd_int{1}{nn} = temp.Dat{nn}.Df;
    
    Hd_int{2}{nn} = temp.Dat{nn+6}.H;
    kd_int{2}{nn} = temp.Dat{nn+6}.Df;
end


%% Analysis

switch option
    
    case 'DC Susc'
        
        Jxx = params(7);
        Jzz = params(8);
        
        HCEF = HCEFf(params(1:6));
        [P, D] = eig(HCEF);
        
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on;
        T = 1:300;
        %plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'csq')
        %plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'msq')
        plot(T, Xab_inv_vT(HCEF,E0,Jxx,T),'b')
        plot(T,  Xc_inv_vT(HCEF,E0,Jzz,T),'r')
        
    case 'DC Susc model diff'
        
        Jxx = params(7);
        Jzz = params(8);
        
        HCEF = HCEFf(params(1:6));
        [P, D] = eig(HCEF);
        
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on;
        T = [4,6,12,20,30,50,70,100:50:250];
        plot(T, XaDatInt - Xab_inv_vT(HCEF,E0,Jxx,T), 'bsq-')
        plot(T, XcDatInt -  Xc_inv_vT(HCEF,E0,Jzz,T), 'rsq-')
    
    case 'principle k'
        
%         Jxx = params(7);
%         Jzz = params(8);
%         
%         HCEF = HCEFf(params(1:6));
%         [P, D] = eig(HCEF);
%         Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
%         E0 = Ev(1);
%         
        figure; hold on;
        
        
        
        H = 0:60; nn = 1; colz = varycolor(6);
        for t = T_RTM
            
            subplot(1,2,2); hold on;
            plot(Hd_int{2}{nn}, kd_int{2}{nn},'sq-','color',colz(nn,:));
            %plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H))
            
            subplot(1,2,1); hold on;
            plot(Hd_int{1}{nn}, kd_int{1}{nn},'sq-','color',colz(nn,:));
            %plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H))
            
            nn = nn + 1;
            
        end
        
        subplot(1,2,1); grid on; box on; xlabel('H [T]'); ylabel('k_a [meV]'); xlim([0,60]); %ylim([0,2.5]);
        subplot(1,2,2); grid on; box on; xlabel('H [T]'); ylabel('k_c [meV]'); xlim([0,60]); %([0,2.5]);
        
    case 'Interp df out'
        
        % pi/2 
        dat1 = load('C:\Users\LeeLabLaptop\Documents\ML_Dec21_RTM\NaYbSe2_RTM_AB.mat');
        dat2 = load('C:\Users\LeeLabLaptop\Documents\ML_Dec21_RTM\NaYbSe2_RTM_C.mat');
        
        startind = 380;
        endind = 870; %(startind:endind)
        
        for n2 = 1:6
            Xd{1}{n2} = dat1.NaYbSe2_RTM_AB(n2).Field(startind:endind);
            Yd{1}{n2} = dat1.NaYbSe2_RTM_AB(n2).frequency(startind:endind);
            
            Xd{2}{n2} = dat2.NaYbSe2_RTM_C(n2).Field(startind:endind);
            Yd{2}{n2} = dat2.NaYbSe2_RTM_C(n2).frequency(startind:endind);
        end
        
        
     
        
        T_RTM = [4,10,20,40,80,160]; o = 0;
       
        
        xi = 2:2:60;
        n4 = 1;
        for n1 = [1,2]

            for n2 = [1:6]
             
                %display(sum(isnan(Xd{n1}{n2})));
                x = Xd{n1}{n2};
                y = Yd{n1}{n2};
                
                ind = find(x <= 2);
                freq0{n1}(n2) = mean(y(ind));
                
                [y ind] = rmoutliers(y,'movmedian',20,'ThresholdFactor',.5);
                x = x(not(ind));
                y = (y - freq0{n1}(n2))/(1e3*amp);
                ind = find(x == max(x));
                x = [x; 60]; y = [y; y(ind(1))];
                
                yi = Interp1NonUnique(x,smooth(y,10),xi);
               
                output{n4}.H  = [0 xi];
                output{n4}.Df = [0 yi];
                n4 = n4+1;
            end
        end
             
    case 'H-Dept. k'
        
        % pi/2 
        dat1 = load('C:\Users\LeeLabLaptop\Documents\ML_Dec21_RTM\NaYbSe2_RTM_AB.mat');
        dat2 = load('C:\Users\LeeLabLaptop\Documents\ML_Dec21_RTM\NaYbSe2_RTM_C.mat');
        
        startind = 380;
        endind = 870; %(startind:endind)
        
        for n2 = 1:6
            Xd{1}{n2} = dat1.NaYbSe2_RTM_AB(n2).Field(startind:endind);
            Yd{1}{n2} = dat1.NaYbSe2_RTM_AB(n2).frequency(startind:endind);
            
            Xd{2}{n2} = dat2.NaYbSe2_RTM_C(n2).Field(startind:endind);
            Yd{2}{n2} = dat2.NaYbSe2_RTM_C(n2).frequency(startind:endind);
        end
        
        
        
        figure; hold on;
        
        T_RTM = [4,10,20,40,80,160]; o = 0;
        
        colz = flipud(varycolor(6));
        
        xi = 2:2:60;
        
        for n1 = [1,2]
            subplot(1,2,n1); hold on;
            for n2 = [1:6]
             
                %display(sum(isnan(Xd{n1}{n2})));
                x = Xd{n1}{n2};
                y = Yd{n1}{n2};
                
                ind = find(x <= 2);
                freq0{n1}(n2) = mean(y(ind));
                
                [y ind] = rmoutliers(y,'movmedian',20,'ThresholdFactor',.5);
                x = x(not(ind));
                y = (y - freq0{n1}(n2))/(1e3*amp);
                ind = find(x == max(x));
                x = [x; 60]; y = [y; y(ind(1))];
                
                yi = Interp1NonUnique(x,smooth(y,10),xi);
                
                %plot( x, (y-freq0{n1}(n2))/(1e3*amp/w0(n2)) + o*(n2-1),'o','color',colz(n2-1,:)
                plot( x, y, 'sq','color', [.75,.75,.75]);
                
                plot( [0 xi], [0 yi],'sq-','color',colz(n2,:),...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5)%,...
                    %'markersize',);
            end
        end
        
%         H = 1:60;
%         
%         Jxx = params(7);
%         Jzz = params(8);
%         
%         HCEF = HCEFf(params(1:6));
%         [P, D] = eig(HCEF);
%         Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
%         E0 = Ev(1);
%         
%         n2 = 1;
%         for t = [4,6,12,20, 30,50,70]
%             subplot(1,2,2); hold on;
%             Y = [0 kc_vH(HCEF,E0,Jxx,Jzz,t,H) + o*n2];
%             plot([0 H],Y,'k--','displayname' ,[num2str(t) ' K model'],...
%                 'linewidth',2)
%             subplot(1,2,1); hold on;
%             Y = [0 kab_vH(HCEF,E0,Jxx,Jzz,t,H) + o*n2];
%             plot([0 H],Y,'k--','displayname' ,[num2str(t) ' K  model'],...
%                 'linewidth',2)
%             n2 = n2+1;
%         end
        
        subplot(1,2,2); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,2,1); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
    case 'H-Dept. k LF limit'
        
        % pi/2
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        T_RTM = [1.5,4,6,12,20, 30,50,70]; o = .75;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        for n1 = [1,2]
            subplot(1,2,n1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                plot( x, (y-freq0{n1}(n2))/(1e3*amp/w0(n2)) + o*(n2-1),'o','color',colz(n2-1,:),...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
            end
        end
        
        H = 1:30;
        
        Jxx = params(7);
        Jzz = params(8);
        
        HCEF = HCEFf(params(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        n2 = 1;
        for t = [4,6,12,20,30,50,70]
            subplot(1,2,2); hold on;
            Xa = power(gJ*muB,2)*(C0/kB)./XaDatInt(n2);
            Xc = power(gJ*muB,2)*(C0/kB)./XcDatInt(n2);
            Y = [0 -power(H,2)*(Xa-Xc) ] + o*n2;
            plot([0 H],Y,'k--','displayname' ,[num2str(t) ' K model'],...
                'linewidth',2)
            subplot(1,2,1); hold on;
            Y = [0 power(H,2)*(Xa-Xc) ] + o*n2;
            plot([0 H],Y,'k--','displayname' ,[num2str(t) ' K  model'],...
                'linewidth',2)
            n2 = n2+1;
        end
        
        subplot(1,2,2); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,2,1); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
    
    case 'MvH data'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {[.7,.7,.7],'r','b','g','c','m','y','k'};
        
        subplot(1,2,1); hold on;
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            Y = gJ*Mset0(ind);
            %plot(H_single,Y, '-', 'color', colz{i+1})
            plot( H_single(2:end-1), diff(diff(Y)) , '-', 'color', colz{i+1})
        end
        
         subplot(1,2,2); hold on;
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            Y = gJ*Mset0(ind);
            %plot(H_single,Y, '-', 'color', colz{i+1-7})
            plot( H_single(2:end-1), diff(diff(Y)) , '-', 'color', colz{i-7+1})
        end  
        
        gxx0 = 3.5237; gzz0 = 1.3336;
        
        for t = [4,6,12,20,30,50,70]
            subplot(1,2,1);
            plot([1,1]*.66/(gxx0*muB/(2*kB*t)), [-3,1]*1e-5 , 'k--') 
            
            subplot(1,2,2);
            plot([1,1]*.66/(gzz0*muB/(2*kB*t)), [-3,1]*1e-5 , 'k--') 
        end
         
    case 'MvT model'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {'r','b','g','c','m','y','k'};
%         
        
        T = [.1:.1:5 6:150]; i = 1;
        for H = [1,5]
            
            
            n1 = 1;
            for t = T
                Y1(n1) = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
                Y2(n1) =  Mc_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
                n1 = n1+1;
            end
            subplot(1,2,1); hold on;  set(gca,'fontsize',20);
            plot(T,1./(Y1/H),'color',colz{i},'displayname', ['H = ' num2str(H) ' T']);
            
            subplot(1,2,2); hold on;  set(gca,'fontsize',20);
            plot(T,1./(Y2/H),'color',colz{i},'displayname', ['H = ' num2str(H) ' T']);
            
            i = i + 1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('m_{ab} [\mu_B]'); xlabel('T [K]'); legend('location','northeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('m_c [\mu_B]'); xlabel('T [K]'); legend('location','northeast')    
           
    case 'nmag'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {'r','b','g','c','m','y','k'};
         m0 = 4;
        
        T = 1:150; i = 1;
        for H = [1,3,6,9,12,15,18]
            
            
            n1 = 1;
            for t = T
                Y1(n1) = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
                Y2(n1) =  Mc_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
                n1 = n1+1;
            end
            subplot(1,2,1); hold on;  set(gca,'fontsize',20);
            plot(T,1-Y1/m0,'color',colz{i},'displayname', ['H = ' num2str(H) ' T']);
            
            subplot(1,2,2); hold on;  set(gca,'fontsize',20);
            plot(T,1-Y2/m0,'color',colz{i},'displayname', ['H = ' num2str(H) ' T']);
            
            i = i + 1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('n_{mag} = 1- m_{ab}/m_s'); xlabel('T [K]'); legend('location','northeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('n_{mag} = 1 - m_{c}/m_s'); xlabel('T [K]'); legend('location','northeast')    
        
    case 'nmag vs H'
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {[.7,.7,.7],'r','b','g','c','m','y','k'};
        
      ms = 4;
%         
        H = 0:.1:60; i = 1;
        for t = [1.5,4,6,12,20,30,50,70]
            
            subplot(1,2,2); hold on;
            
            Y = Mc_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
            plot(H,1-Y/ms,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            subplot(1,2,1); hold on;
            
            Y = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
            plot(H,1-Y/ms,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            i = i +1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('n_{mag} = 1-m_{ab}/m_s'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('n_{mag} = 1-m_{c}/m_s'); xlabel('\mu_0H [T]'); legend('location','southeast')  
               
    case 'MvH model alt'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {[.7,.7,.7],'r','b','g','c','m','y','k'};
        
        dH = 0.01;
        H = 0:dH:60; i = 1;
        for t = [1.5,4,6,12,20,30,50,70]
            
            
            subplot(1,2,1); hold on;
            
            Y = Mab_vH_alt(HCEF,E0,Jxx,Jzz,t,H,dH)/muB;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            
            subplot(1,2,2); hold on;
            
            Y = Mc_vH_alt(HCEF,E0,Jxx,Jzz,t,H,dH)/muB;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            
            i = i +1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('m_{ab} [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('m_c [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
    case 'MvH model'
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {[.7,.7,.7],'r','b','g','c','m','y','k'};
        
        subplot(1,2,1); hold on;
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            plot( [0 H_single], gJ*[0; Mset0(ind)], 'sq', 'color', colz{i+1})
            
            %plot([0 H_single], power(gJ,2)*power(power(C0*(muB/kB),-1)*XaDatInt(i),-1)*[0 H_single],'k--','displayname',[num2str(T_RTM(i+1))] )
        end
        
         subplot(1,2,2); hold on;
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            plot( [0 H_single], gJ*[0; Mset0(ind)], 'sq', 'color', colz{i-7+1})
            %plot([0 H_single], power(gJ,2)*power(power(C0*(muB/kB),-1)*XcDatInt(i-7),-1)*[0 H_single],'k--','displayname',[num2str(T_RTM(i+1-7))] )
        end
%         
        
        H = 0:.1:60; i = 1;
        for t = [2,4,6,8,14,20,30,50];
            
            subplot(1,2,2); hold on;
            
            Y = Mc_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
            
            output.Mc{i} = Y;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            subplot(1,2,1); hold on;
            
            Y = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)/muB;
            
            output.Ma{i} = Y;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            i = i +1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('m_{ab} [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('m_c [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
    case 'fmin MvH'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        
        f = @(x) CummErr_MvH([x, Jxx, Jzz, 0]);
        
        p0 = params(1:6);
        pf = fminsearchbnd(f,[p0],[-1,-1,-1,-1,-1,-1], [1,1,1,1,1,1]);
        %pf = fminsearchbnd(f,[p0, Jxx, Jzz, 0],[-1,-1,-1,-1,-1,-1,Jxx,0,-10], [1,1,1,1,1,1,Jxx,10,10]);
        Jxx = pf(7);  Jzz = pf(8);
        display(pf); X0 = pf(9);
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        figure; hold on; colz = {[.7,.7,.7],'r','b','g','c','m','y','k'};
        
        subplot(1,2,1); hold on;
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            plot( [0 H_single], gJ*[0; Mset0(ind)], 'sq-', 'color', colz{i+1})
        end
        
         subplot(1,2,2); hold on;
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            plot( [0 H_single], gJ*[0; Mset0(ind)], 'sq-', 'color', colz{i-7+1})
        end
%         
        
        H = 0:.1:7; i = 1;
        for t = [1.5,4,6,12,20,30,50,70]
            
            subplot(1,2,2); hold on;
            
            Y = Mc_vH(HCEF,E0,Jxx,Jzz,t,H)/muB - X0*H;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            subplot(1,2,1); hold on;
            
            Y = Mab_vH(HCEF,E0,Jxx,Jzz,t,H)/muB - X0*H ;
            plot(H,Y,'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K'])
            %d2Y = (C2(1)*Y(1:end-4) + C2(2)*Y(2:end-3) + C2(3)*Y(3:end-2) + C2(4)*Y(4:end-1) + C2(5)*Y(5:end));
            %plot(H(2:end-1),diff(diff(smooth(Y,3))),'-','linewidth',2, 'color', colz{i},'displayname', [num2str(t) ' K']);
            
            i = i +1;
        end
        subplot(1,2,1); grid on; box on;
        ylabel('m_{ab} [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')
        
        subplot(1,2,2); grid on; box on;
        ylabel('m_c [\mu_B]'); xlabel('\mu_0H [T]'); legend('location','southeast')    
        
    case 'Trans Succ Analytical'
        
        Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        
        Dat2 = CEF_Spectrum_D3d('Perturbative expression',p0);

        ax0 = Dat2.ax0{1}; ax1 = Dat2.ax0{2}; ax2 = Dat2.ax0{3}; ax3 = Dat2.ax0{4};
        az0 = Dat2.az0{1}; az1 = Dat2.az0{2}; az2 = Dat2.az0{3}; az3 = Dat2.az0{4};
        
        Bx0 = 2*kB*ax0/power(muB*gJ,2);
        Bx1 = 2*kB*ax1/power(muB*gJ,2);
        Bx2 = 2*kB*ax2/power(muB*gJ,2);
        Bx3 = 2*kB*ax3/power(muB*gJ,2);
        
        Bz0 = 2*kB*az0/power(muB*gJ,2);
        Bz1 = 2*kB*az1/power(muB*gJ,2);
        Bz2 = 2*kB*az2/power(muB*gJ,2);
        Bz3 = 2*kB*az3/power(muB*gJ,2);
        
        gxx0 = Dat2.gxx{1}; gxx1 = Dat2.gxx{2}; gxx2 = Dat2.gxx{3}; gxx3 = Dat2.gxx{4};
        gzz0 = Dat2.gzz{1}; gzz1 = Dat2.gzz{2}; gzz2 = Dat2.gzz{3}; gzz3 = Dat2.gzz{4};
        
        Ax0 = power( gxx0/(2*gJ),2);
        Ax1 = power( gxx1/(2*gJ),2);
        Ax2 = power( gxx2/(2*gJ),2);
        Ax3 = power( gxx3/(2*gJ),2);
        
        Az0 = power( gzz0/(2*gJ),2);
        Az1 = power( gzz1/(2*gJ),2);
        Az2 = power( gzz2/(2*gJ),2);
        Az3 = power( gzz3/(2*gJ),2);
        
        D10 = Dat2.D{1}; D20 = Dat2.D{2}; D30 = Dat2.D{3};
        
        T = 1:300;
        
        E10 = exp(-Dat2.D{1}./(kB.*T));
        E20 = exp(-Dat2.D{2}./(kB.*T));
        E30 = exp(-Dat2.D{3}./(kB.*T));
        
        Xinv1 = 2.0416*(T.*(1 + E10 + E20 + E30)./(...
            Ax0 - Bx0*T + (Ax1 - Bx1*T).*E10 + (Ax2 - Bx2*T).*E20 + (Ax3 - Bx3*T).*E30) + 6*Jxx);
        
        Xinv2 = 2.0416*(T.*(1 + E10 + E20 + E30)./(...
            Az0 - Bz0*T + (Az1 - Bz1*T).*E10 + (Az2 - Bz2*T).*E20 + (Az3 - Bz3*T).*E30) + 6*Jzz);
        
        %figure; hold on;
        %plot(T, Xinv1, 'r')
        %plot(T, Xinv2, 'b')
        
        
        D10 = Dat2.D{1}./kB;
        D20 = Dat2.D{2}./kB;
        D30 = Dat2.D{3}./kB;
        
        XT_inv1 = @(t,h) 2.0416*(t.*(  exp(-az0*power(h,2)./(kB*t)).*cosh(muB*h*gzz0./(2*kB*t))...
            + exp(-az1*power(h,2)./(kB*t)).*cosh(muB*h*gzz1./(2*kB*t)).*exp(-D10./t)...
            + exp(-az2*power(h,2)./(kB*t)).*cosh(muB*h*gzz2./(2*kB*t)).*exp(-D20./t)...
            + exp(-az3*power(h,2)./(kB*t)).*cosh(muB*h*gzz3./(2*kB*t)).*exp(-D30./t))./(...
            (Ax0.*sinh(muB*h*gzz0./(2*kB*t))./(muB*h*gzz0./(2*kB*t)) - Bx0*t.*cosh(muB*h*gzz0./(2*kB*t)) ).*exp(-az0*power(h,2)./(kB*t))...
            + (Ax1.*sinh(muB*h*gzz1./(2*kB*t))./(muB*h*gzz1./(2*kB*t)) - Bx1*t.*cosh(muB*h*gzz1./(2*kB*t)) ).*exp(-az1*power(h,2)./(kB*t)).*exp(-D10./t)...
            + (Ax2.*sinh(muB*h*gzz2./(2*kB*t))./(muB*h*gzz2./(2*kB*t)) - Bx2*t.*cosh(muB*h*gzz2./(2*kB*t)) ).*exp(-az2*power(h,2)./(kB*t)).*exp(-D20./t)...
            + (Ax3.*sinh(muB*h*gzz3./(2*kB*t))./(muB*h*gzz3./(2*kB*t)) - Bx3*t.*cosh(muB*h*gzz3./(2*kB*t)) ).*exp(-az3*power(h,2)./(kB*t)).*exp(-D30./t) ) + TCWx);
        
        XT_inv2 = @(t,h) 2.0416*(t.*(  exp(-ax0*power(h,2)./(kB*t)).*cosh(muB*h*gxx0./(2*kB*t))...
            + exp(-ax1*power(h,2)./(kB*t)).*cosh(muB*h*gxx1./(2*kB*t)).*exp(-D10./t)...
            + exp(-ax2*power(h,2)./(kB*t)).*cosh(muB*h*gxx2./(2*kB*t)).*exp(-D20./t)...
            + exp(-ax3*power(h,2)./(kB*t)).*cosh(muB*h*gxx3./(2*kB*t)).*exp(-D30./t))./(...
            (Az0.*sinh(muB*h*gxx0./(2*kB*t))./(muB*h*gxx0./(2*kB*t)) - Bz0*t.*cosh(muB*h*gxx0./(2*kB*t)) ).*exp(-ax0*power(h,2)./(kB*t))...
            + (Az1.*sinh(muB*h*gxx1./(2*kB*t))./(muB*h*gxx1./(2*kB*t)) - Bz1*t.*cosh(muB*h*gxx1./(2*kB*t)) ).*exp(-ax1*power(h,2)./(kB*t)).*exp(-D10./t)...
            + (Az2.*sinh(muB*h*gxx2./(2*kB*t))./(muB*h*gxx2./(2*kB*t)) - Bz2*t.*cosh(muB*h*gxx2./(2*kB*t)) ).*exp(-ax2*power(h,2)./(kB*t)).*exp(-D20./t)...
            + (Az3.*sinh(muB*h*gxx3./(2*kB*t))./(muB*h*gxx3./(2*kB*t)) - Bz3*t.*cosh(muB*h*gxx3./(2*kB*t)) ).*exp(-ax3*power(h,2)./(kB*t)).*exp(-D30./t) ) + TCWz);
        
        
        
        
        h = 0:.5:60;
        
        figure; hold on;
        for t = [4,6,12,20,30,50,70]
            subplot(1,2,1); hold on;
            plot(h, (.0104/.1002)./XT_inv1(t,h))
            
            subplot(1,2,2); hold on;
            plot(h, (.0104/.1002)./XT_inv2(t,h))
        end
        
    case 'Transervse Susc'
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        p0 = params(1:6);
        pf = p0;
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            %Ma{i} = muB*gJ*[Mset0(ind)];
            %Xa(i) = muB*gJ*f_mag(M_params(i,:),.1)./.1;
        end
        
        %subplot(1,2,2); hold on;
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            %Mc{i-7} = muB*gJ*[Mset0(ind)];
            %Xc(i-7) = muB*gJ*f_mag(M_params(i,:),.1)./.1;

        end
        
        
        figure; hold on
        H = 0:60; i = 1;
        for t = [4]
            
            subplot(1,2,1); hold on;
            plot([H],[power(H,2).*XTc_vH(HCEF,E0,Jxx,Jzz,t,H)],'linewidth',2)
            
            subplot(1,2,2); hold on;
            plot([H],[power(H,2).*XTab_vH(HCEF,E0,Jxx,Jzz,t,H)],'linewidth',2)
            
            i = i+1;
        end
          
    case 'model'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        o = 0;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        
        for i = 1:7
            ind = (i-1)*L1+1:i*L1;
            Xa(i) = muB*gJ*f_mag(M_params(i,:),.1)./.1;
        end
        for i = 8:14
            ind = (i-1)*L1+1:i*L1;
            Xc(i-7) = muB*gJ*f_mag(M_params(i,:),.1)./.1;
        end
        X_1 = 3:10;
        X_2 = power(X_1,2);
        
        for n1 = [1,2]
            subplot(1,3,n1+1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                %ind = find(x == min(x));
                %freq0{n1}(n2) = y(ind);
                
                ind = find(x <= 12);
                k_susc = power(-1,n1-1)*X_2.*(Xa(n2-1) - Xc(n2-1));    
                k_LFint = Interp1NonUnique(x(ind),y(ind),X_1); 
                g = @(x) sum(power( (k_LFint - x*1e3)./(1e3*amp/x) - k_susc,2));
                
                freq0{n1}(n2) = fminsearch(g,w0(n2));
                
                k_dat = (y-freq0{n1}(n2)*1e3)/(1e3*amp/w0(n2));
                plot( x, k_dat + o*(n2-1),'o','color',[.5,.5,.5],...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
                
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
                %display(sum(isnan(kd_int{n1}{n2} )))
                
                %plot(Hd_int{n1}{n2}, kd_int{n1}{n2}, 'sq', 'color', colz(n2-1,:))
            end
        end
        
        subplot(1,3,3); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,2); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,1); hold on
        grid on; box on; set(gca,'fontsize',16);
        xlabel('T [K]')
        ylabel('\chi^{-1}')
        plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'bsq')
        plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'rsq')
        
        
        Jxx = params(7); Jzz = params(8);
        %Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        display([Jxx,Jzz])
        f = @(x) CummErr([x, Jxx, Jzz, 0]);
        
        p0 = params(1:6);
        pf = p0;%fminsearchbnd(f,p0,[-1,-1,-1,-1,-1,-1], [1,1,1,1,1,1]);
        
        %display(pf)
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        
        H = 0:60;
        for t = [4,6,12,20,30,50,70]
            
            subplot(1,3,3); hold on;
            plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
            subplot(1,3,2); hold on;
            plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
        end
        
        T = 1:300;
        subplot(1,3,1)
        plot(T, Xab_inv_vT(HCEF,E0,Jxx,T),'b')
        plot(T,  Xc_inv_vT(HCEF,E0,Jzz,T),'r')
    
    case 'fmin susc only'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        f= figure; hold on;
        f.Position = [100 100 2000 800];
        
        o = 0;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        for n1 = [1,2]
            subplot(1,3,n1+1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                plot( x, k_dat + o*(n2-1),'o','color',[.5,.5,.5],...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
                
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
                %display(sum(isnan(kd_int{n1}{n2} )))
                
                %plot(Hd_int{n1}{n2}, kd_int{n1}{n2}, 'sq', 'color', colz(n2-1,:))
            end
        end
        
        subplot(1,3,3); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,2); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,1); hold on
        grid on; box on; set(gca,'fontsize',16);
        xlabel('T [K]')
        ylabel('\chi^{-1}')
        plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'bsq')
        plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'rsq')
        
        
        
        Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        f = @(x) CummErrSuscOnly([x, Jxx, Jzz, 0]);
        
        p0 = params(1:6);
        pf = fminsearchbnd(f,p0, 2*[-1,-1,-1,-1,-1,-1], 2*[1,1,1,1,1,1]);
        
        display(pf)
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        H = 0:60;
        for t = [4,6,12,20,30,50,70]
            
            subplot(1,3,3); hold on;
            plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
            subplot(1,3,2); hold on;
            plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
        end
        
        T = 1:300;
        subplot(1,3,1)
        plot(T, Xab_inv_vT(HCEF,E0,Jxx,T),'b')
        plot(T,  Xc_inv_vT(HCEF,E0,Jzz,T),'r')    
        
    case 'fmin k'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        o = 0;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        for n1 = [1,2]
            subplot(1,3,n1+1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                plot( x, k_dat + o*(n2-1),'o','color',[.5,.5,.5],...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
                
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
                %display(sum(isnan(kd_int{n1}{n2} )))
                
                %plot(Hd_int{n1}{n2}, kd_int{n1}{n2}, 'sq', 'color', colz(n2-1,:))
            end
        end
        
        subplot(1,3,3); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,2); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,1); hold on
        grid on; box on; set(gca,'fontsize',16);
        xlabel('T [K]')
        ylabel('\chi^{-1}')
        plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'bsq')
        plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'rsq')
        
        
        
        Jxx = params(7);
        Jzz = params(8);
        %LTparams(1)/6; Jzz = LTparams(4)/6;
        
        f = @(x) CummErr([x, Jxx, Jzz, 0]);
        
        p0 = params(1:6);
        pf = fminsearchbnd(f,p0,[-1,-1,-1,-1,-1,-1], [1,1,1,1,1,1]);
        
        display(pf)
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        
        H = 0:60;
        for t = [4,6,12,20,30,50,70]
            
            subplot(1,3,3); hold on;
            plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
            subplot(1,3,2); hold on;
            plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
        end
        
        T = 1:300;
        subplot(1,3,1)
        plot(T, Xab_inv_vT(HCEF,E0,Jxx,T),'b')
        plot(T,  Xc_inv_vT(HCEF,E0,Jzz,T),'r')
        
    case 'fmin k no plot'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        for n1 = [1,2]
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );

            end
        end

        Jxx = LTparams(1)/6; Jzz = LTparams(4)/6;
        
        f = @(x) CummErr([x, Jxx, Jzz, 0]);
        
        p0 = params(1:6);
        pf = fminsearchbnd(f,p0,2*[-1,-1,-1,-1,-1,-1], 2*[1,1,1,1,1,1]);
        
       % display(pf)
        
        output.p = pf;
        output.Res = f(pf);
        
    case 'get Res'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        for n1 = [1,2]
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
      
            end
        end
        
        f = @(x) CummErr([x 0]);
        
        p0 = params(1:8);
        output = f(p0);
          
    case 'fmin k float Jaa'
        
        Hset1 = [2:3:29];
        Hset2 = [4:3:58];
        
        Hd_int{1} = {Hset2, Hset2, Hset1, Hset1, Hset1, Hset2, Hset1, Hset2};
        Hd_int{2} = {Hset2, Hset2, Hset1, Hset2, Hset1, Hset2, Hset1, Hset2};
        
        dat{1}{1} = load('p048_031020_TDO004.dat');
        dat{1}{2} = load('p043_031020_TDO003.dat');
        dat{1}{3} = load('p050_031020_TDO002.dat');
        dat{1}{4} = load('p052_031020_TDO001.dat');
        dat{1}{5} = load('p058_031020_TDO001.dat');
        dat{1}{6} = load('p067_031020_TDO001.dat');
        dat{1}{7} = load('p079_031020_TDO001.dat');
        dat{1}{8} = load('p092_031020_TDO001.dat');
        
        % 0
        dat{2}{1} = load('p059_031120_TDO001.dat');
        dat{2}{2} = load('p042_030920_TDO002.dat');
        dat{2}{3} = load('p057_031120_TDO001.dat');
        dat{2}{4} = load('p055_031120_TDO001.dat');
        dat{2}{5} = load('p050_031120_TDO001.dat');
        dat{2}{6} = load('p044_031120_TDO001.dat');
        dat{2}{7} = load('p037_031120_TDO001.dat');
        dat{2}{8} = load('p097_031020_TDO001.dat');
        
        startind = 251;
        endind = 2749-450;
        
        figure; hold on;
        
        o = 0;
        
        colz = [0,0,1;...
            0.3,0.75,0.93;...
            0.47,0.67,0.19;...
            0,1,0;...
            0.87,0.49,0;...
            0.75,0,0.75;...
            1,0,0;...
            ];
        
        for n1 = [1,2]
            subplot(1,3,n1+1); hold on;
            for n2 = [2:8]
                x = dat{n1}{n2}(startind:end-endind,1);
                y = dat{n1}{n2}(startind:end-endind,2);
                
                ind = find(x == min(x));
                freq0{n1}(n2) = y(ind);
                
                k_dat = (y-freq0{n1}(n2))/(1e3*amp/w0(n2));
                plot( x, k_dat + o*(n2-1),'o','color',[.5,.5,.5],...
                    'displayname', [num2str(T_RTM(n2)) ' K Data'],...
                    'markerfacecolor',[1,1,1],'linewidth',1.5,...
                    'markersize',8);
                
                kd_int{n1}{n2} = Interp1NonUnique(x, k_dat, Hd_int{n1}{n2} );
                %display(sum(isnan(kd_int{n1}{n2} )))
                
                %plot(Hd_int{n1}{n2}, kd_int{n1}{n2}, 'sq', 'color', colz(n2-1,:))
            end
        end
        
        subplot(1,3,3); title('\theta = 0, H || c')
        %legend('location','southwest');
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,2); title('\theta = \pi/2, H || ab')
        grid on; box on; set(gca,'fontsize',16);
        xlabel('\mu_0H [T]')
        ylabel('\Deltaf [Hz]')
        
        subplot(1,3,1); hold on
        grid on; box on; set(gca,'fontsize',16);
        xlabel('T [K]')
        ylabel('\chi^{-1}')
        plot([4,6,12,20,30,50,70,100:50:250], XaDatInt, 'bsq')
        plot([4,6,12,20,30,50,70,100:50:250], XcDatInt, 'rsq')
        
        
        f = @(x) CummErr([x 0]);
        
        p0 = params(1:8);
        pf = fminsearchbnd(f,p0,[-1,-1,-1,-1,-1,-1,0,0], [1,1,1,1,1,1,10,10]);
        
        display(pf)
        
        HCEF = HCEFf(pf(1:6));
        [P, D] = eig(HCEF);
        Ev = sort(real([D(1,1) D(2,2) D(3,3) D(4,4) D(5,5) D(6,6) D(7,7) D(8,8)]));
        E0 = Ev(1);
        
        Jxx = pf(7); Jzz = pf(8);
        
        H = 0:60;
        for t = [4,6,12,20,30,50,70]
            
            subplot(1,3,3); hold on;
            plot(H,kc_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
            subplot(1,3,2); hold on;
            plot(H,kab_vH(HCEF,E0,Jxx,Jzz,t,H),'linewidth',2)
            
        end
        
        T = 1:300;
        subplot(1,3,1)
        plot(T, Xab_inv_vT(HCEF,E0,Jxx,T),'b')
        plot(T,  Xc_inv_vT(HCEF,E0,Jzz,T),'r')
        
        output = pf;
        
           
end


end