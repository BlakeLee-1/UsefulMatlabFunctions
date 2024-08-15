BathCell = load('C:\Users\LeeLabLaptop\Documents\Minhyea Lee Research\Maglab_March2019\CXfits_March2019\SDCXcell.mat');
BathFit = BathCell.CXcell{2};

deg = 7;
Dat = load('CX_R.csv'); x = Dat(:,1); y = Dat(:,2);
HotFitOG = RtoT_cal_inputRT(x,y,'getfit',deg,'Yes','Polynomial',[]);

Dat = load('C:\Users\LeeLabLaptop\Documents\Minhyea Lee Research\Maglab_Feb20\CXcell_QUT_v1');
ColdFitOG = Dat.CXcell_QUT{1};

% DAQ.CXcell = {P Cal, M Cal, N Cal, SD Bath Cal, L Bath Cal}
DAQ.CXcell = {HotFitOG, ColdFitOG, BathFit};

load('HtrFit_2017_9_1_15_59_51.mat')


KYbSe2_PPMS_Nov21_driverv0( 'Plot_KvT_all', DAQ.CXcell, HtrFit, 1,1)

AnalyzeKvB_KYbSe2_Nov21('Plot KxxvB',DAQ.CXcell, HtrFit,1)