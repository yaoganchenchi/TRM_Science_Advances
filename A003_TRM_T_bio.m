function A003_TRM_T_bio()
% MATLAB R2019b
% Dated: July-19-2020
% Author: Chi Chen @ Boston University
% Email: chenchi@bu.edu

case_res = 'f05_g17';
this_season = 'yearly';
path_out = sprintf('/projectnb/amazondr/data12/cliveg/chenchi/LST_LAI/CLM5_TRM_output/%s/%s/',case_res,this_season);
if ~(exist(path_out,'dir') ==7)
    mkdir(path_out);
end

lc_mask = loadMatData('/projectnb/amazondr/data12/cliveg/chenchi/MCD12C1/output_c5/MCD12C1.CMG050.C5.IGBP.MODE.LC.2001.2012.Allveg.mat');
lai_trend = loadMatData('/projectnb/amazondr/data12/cliveg/chenchi/LST_LAI/TrendAnalysis/MCD_LAI_C6_0_50CMG_monthly_2000_2014_AnnualAve/Trend_2000_2014/trendMatMK/MCD_LAI_C6_0_50CMG_monthly_2000_2014_AnnualAve_trend_MK.mat') ;
lai_sig = loadMatData('/projectnb/amazondr/data12/cliveg/chenchi/LST_LAI/TrendAnalysis/MCD_LAI_C6_0_50CMG_monthly_2000_2014_AnnualAve/Trend_2000_2014/trendMatMK/MCD_LAI_C6_0_50CMG_monthly_2000_2014_AnnualAve_pValue_MK.mat') ;
lai_trend(lai_sig>0.1 | lai_sig<0)=nan;

% get constants
constants.sigma_sbc =  5.67*10^-8; % stefan-boltzmann constant
constants.Lv = 2.501*10^6; % latent heat of vaporization/condensation
constants.Cp=1005; % specific heat capacity, Cp=Cv+R
constants.Rv=461.5;% gas constant of water vapour
constants.R_gas = 287.058; %dry air gas constant
constants.Pstd = 101325; % standard pressure
constants.Avogadro = 6.02214*10^26; % avogadro's number
constants.k_bol = 1.38065*10^-23; % boltzmann constant
constants.MWda = 28.966; % molecular weight of dry air
constants.syear=2000;
constants.eyear=2014;
constants.path_output = path_out;
constants.case_res = case_res;

nyear = constants.eyear - constants.syear +1;

LST2PARA = loadMatData(sprintf('%sLST2PARA_median_sensitivity_%04d_%04d_%s.mat',path_out,constants.syear,constants.eyear,case_res));
PARA2LAI = loadMatData(sprintf('%sPARA2LAI_median_sensitivity_%04d_%04d_%s.mat',path_out,constants.syear,constants.eyear,case_res));

temp = fieldnames(LST2PARA);
temp2 = LST2PARA.(temp{1});
[clm_row,clm_col] = size(temp2);
res_lat = 180/clm_row;
res_lon = 360/clm_col;
R_clm =  makerefmat(-180+res_lon/2, 90-res_lat/2, res_lon, -res_lat);

% individual sensitivity
dTs_dALBEDO_temp=LST2PARA.dTs_dALBEDO_median;
dTs_dALBEDO_temp(1,1)=-999;
[dTs_dALBEDO] = resizem(dTs_dALBEDO_temp,[360 720],R_clm);
dTs_dALBEDO(lc_mask~=1) = nan;
output.dTs_dALBEDO=dTs_dALBEDO;

dTs_dRA_temp=LST2PARA.dTs_dRA_median;
dTs_dRA_temp(1,1)=-999;
[dTs_dRA] = resizem(dTs_dRA_temp,[360 720],R_clm);
dTs_dRA(lc_mask~=1) = nan;
output.dTs_dRA=dTs_dRA;

dTs_dRS_temp=LST2PARA.dTs_dRS_median;
dTs_dRS_temp(1,1)=-999;
[dTs_dRS] = resizem(dTs_dRS_temp,[360 720],R_clm);
dTs_dRS(lc_mask~=1) = nan;
output.dTs_dRS=dTs_dRS;

dALBEDO_dLAI_temp=PARA2LAI.dALBEDO_dLAI_median;
dALBEDO_dLAI_temp(1,1)=-999;
[dALBEDO_dLAI] = resizem(dALBEDO_dLAI_temp,[360 720],R_clm);
dALBEDO_dLAI(lc_mask~=1) = nan;
output.dALBEDO_dLAI=dALBEDO_dLAI;


dRA_dLAI_temp=PARA2LAI.dRA_dLAI_median;
dRA_dLAI_temp(1,1)=-999;
[dRA_dLAI] = resizem(dRA_dLAI_temp,[360 720],R_clm);
dRA_dLAI(lc_mask~=1) = nan;
output.dRA_dLAI=dRA_dLAI;

dRS_dLAI_temp=PARA2LAI.dRS_dLAI_median;
dRS_dLAI_temp(1,1)=-999;
[dRS_dLAI] = resizem(dRS_dLAI_temp,[360 720],R_clm);
dRS_dLAI(lc_mask~=1) = nan;
output.dRS_dLAI=dRS_dLAI;


%% sen/delta Ts due to biophyscial parameter becasue of changes in LAI
% total bio sensitivity
sen_Ts_bio_lai_temp = LST2PARA.dTs_dALBEDO_median.*PARA2LAI.dALBEDO_dLAI_median + LST2PARA.dTs_dRA_median.*PARA2LAI.dRA_dLAI_median + LST2PARA.dTs_dRS_median.*PARA2LAI.dRS_dLAI_median + LST2PARA.dTs_dG_median.*PARA2LAI.dG_dLAI_median + LST2PARA.dTs_dEMIS_median.*PARA2LAI.dEMIS_dLAI_median;
sen_Ts_bio_lai_temp(1,1)=-999;
[sen_Ts_bio_lai_with_G_EMIS] = resizem(sen_Ts_bio_lai_temp,[360 720],R_clm);
sen_Ts_bio_lai_with_G_EMIS(lc_mask~=1) = nan;
% total bio LAI-LST change
delta_Ts_bio_lai_with_G_EMIS = sen_Ts_bio_lai_with_G_EMIS.*lai_trend*nyear;

% bio albedo sensitivity
sen_Ts_bio_ALBEDO_lai_temp = LST2PARA.dTs_dALBEDO_median.*PARA2LAI.dALBEDO_dLAI_median;
sen_Ts_bio_ALBEDO_lai_temp(1,1)=-999;
[sen_Ts_bio_ALBEDO_lai] = resizem(sen_Ts_bio_ALBEDO_lai_temp,[360 720],R_clm);
sen_Ts_bio_ALBEDO_lai(lc_mask~=1) = nan;
% bio ALBEDO LAI-LST change
delta_Ts_bio_ALBEDO_lai = sen_Ts_bio_ALBEDO_lai.*lai_trend*nyear;

% bio RA sensitivity
sen_Ts_bio_RA_lai_temp = LST2PARA.dTs_dRA_median.*PARA2LAI.dRA_dLAI_median;
sen_Ts_bio_RA_lai_temp(1,1)=-999;
[sen_Ts_bio_RA_lai] = resizem(sen_Ts_bio_RA_lai_temp,[360 720],R_clm);
sen_Ts_bio_RA_lai(lc_mask~=1) = nan;
% bio RA LAI-LST change
delta_Ts_bio_RA_lai = sen_Ts_bio_RA_lai.*lai_trend*nyear;

% bio RS sensitivity
sen_Ts_bio_RS_lai_temp = LST2PARA.dTs_dRS_median.*PARA2LAI.dRS_dLAI_median;
sen_Ts_bio_RS_lai_temp(1,1)=-999;
[sen_Ts_bio_RS_lai] = resizem(sen_Ts_bio_RS_lai_temp,[360 720],R_clm);
sen_Ts_bio_RS_lai(lc_mask~=1) = nan;
% bio RS LAI-LST change
delta_Ts_bio_RS_lai = sen_Ts_bio_RS_lai.*lai_trend*nyear;

% bio GH sensitivity
sen_Ts_bio_G_lai_temp = LST2PARA.dTs_dG_median.*PARA2LAI.dG_dLAI_median;
sen_Ts_bio_G_lai_temp(1,1)=-999;
[sen_Ts_bio_G_lai] = resizem(sen_Ts_bio_G_lai_temp,[360 720],R_clm);
sen_Ts_bio_G_lai(lc_mask~=1) = nan;
% bio G LAI-LST change
delta_Ts_bio_G_lai = sen_Ts_bio_G_lai.*lai_trend*nyear;

% bio EMIS sensitivity
sen_Ts_bio_EMIS_lai_temp = LST2PARA.dTs_dEMIS_median.*PARA2LAI.dEMIS_dLAI_median;
sen_Ts_bio_EMIS_lai_temp(1,1)=-999;
[sen_Ts_bio_EMIS_lai] = resizem(sen_Ts_bio_EMIS_lai_temp,[360 720],R_clm);
sen_Ts_bio_EMIS_lai(lc_mask~=1) = nan;
% bio EMIS LAI-LST change
delta_Ts_bio_EMIS_lai = sen_Ts_bio_EMIS_lai.*lai_trend*nyear;


output.sen_Ts_bio_lai_with_G_EMIS=sen_Ts_bio_lai_with_G_EMIS;
output.delta_Ts_bio_lai_with_G_EMIS=delta_Ts_bio_lai_with_G_EMIS;

output.sen_Ts_bio_ALBEDO_lai=sen_Ts_bio_ALBEDO_lai;
output.delta_Ts_bio_ALBEDO_lai=delta_Ts_bio_ALBEDO_lai;
output.sen_Ts_bio_RA_lai=sen_Ts_bio_RA_lai;
output.delta_Ts_bio_RA_lai=delta_Ts_bio_RA_lai;
output.sen_Ts_bio_RS_lai=sen_Ts_bio_RS_lai;
output.delta_Ts_bio_RS_lai=delta_Ts_bio_RS_lai;
output.sen_Ts_bio_G_lai=sen_Ts_bio_G_lai;
output.delta_Ts_bio_G_lai=delta_Ts_bio_G_lai;
output.sen_Ts_bio_EMIS_lai=sen_Ts_bio_EMIS_lai;
output.delta_Ts_bio_EMIS_lai=delta_Ts_bio_EMIS_lai;

savename = sprintf('%sTRM_Ts_bio_lai_results_%04d_%04d.mat',path_out,constants.syear,constants.eyear);
save(savename,'output','-v7.3');

end
