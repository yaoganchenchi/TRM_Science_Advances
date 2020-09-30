function A001_dTs_to_dBio()
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
constants.this_season=this_season;

% get dim from sample data
if strcmp (constants.case_res,'f05_g17')
    this_case_name = 'I2000Clm50Sp.trend.f05_g17';
end
path_in2 = sprintf('/projectnb/amazondr/data12/cliveg/chenchi/LST_LAI/CLM5_TRM_input/%s/%s/',this_case_name,this_season);
LWin = loadMatData(sprintf('%sFLDS/%s.FLDS.%04d.mat',path_in2,this_case_name,2000));

% start calculation
nyear = constants.eyear - constants.syear +1;
icount=0;
[nrow,ncol]=size(LWin);
dTs_dALBEDO_stack = nan(nrow,ncol,nyear);
dTs_dRA_stack = nan(nrow,ncol,nyear);
dTs_dRS_stack = nan(nrow,ncol,nyear);
dTs_dEMIS_stack = nan(nrow,ncol,nyear);
dTs_dSWin_stack = nan(nrow,ncol,nyear);
dTs_dLWin_stack = nan(nrow,ncol,nyear);
dTs_dQA_stack = nan(nrow,ncol,nyear);
dTs_dTA_stack = nan(nrow,ncol,nyear);
dTs_dG_stack = nan(nrow,ncol,nyear);
for iyear = constants.syear:constants.eyear
    icount=icount+1;
    % get LST sensitivities
    [dTs_dALBEDO,dTs_dRA,dTs_dRS,dTs_dEMIS,dTs_dSWin,dTs_dLWin,dTs_dQA,dTs_dTA,dTs_dG]=get_TRM_LST_sensitivities(iyear,constants);
    dTs_dALBEDO_stack(:,:,icount)=dTs_dALBEDO;
    dTs_dRA_stack(:,:,icount)=dTs_dRA;
    dTs_dRS_stack(:,:,icount)=dTs_dRS;
    dTs_dEMIS_stack(:,:,icount)=dTs_dEMIS;
    dTs_dSWin_stack(:,:,icount)=dTs_dSWin;
    dTs_dLWin_stack(:,:,icount)=dTs_dLWin;
    dTs_dQA_stack(:,:,icount)=dTs_dQA;
    dTs_dTA_stack(:,:,icount)=dTs_dTA;
    dTs_dG_stack(:,:,icount)=dTs_dG;
end

%% median sensitivity 2000-2014
output.dTs_dALBEDO_median = nanmedian(dTs_dALBEDO_stack,3);
output.dTs_dRA_median = nanmedian(dTs_dRA_stack,3);
output.dTs_dRS_median = nanmedian(dTs_dRS_stack,3);
output.dTs_dEMIS_median = nanmedian(dTs_dEMIS_stack,3);
output.dTs_dSWin_median = nanmedian(dTs_dSWin_stack,3);
output.dTs_dLWin_median = nanmedian(dTs_dLWin_stack,3);
output.dTs_dQA_median = nanmedian(dTs_dQA_stack,3);
output.dTs_dTA_median = nanmedian(dTs_dTA_stack,3);
output.dTs_dG_median = nanmedian(dTs_dG_stack,3);
savename = sprintf('%sLST2PARA_median_sensitivity_%04d_%04d_%s.mat',path_out,constants.syear,constants.eyear,case_res);
save(savename,'output','-v7.3');
clear output
end

%% get LST sensitivies to other parameters
function [dTs_dALBEDO,dTs_dRA,dTs_dRS,dTs_dEMIS,dTs_dSWin,dTs_dLWin,dTs_dQA,dTs_dTA,dTs_dG]=get_TRM_LST_sensitivities(iyear,constants)
addpath('/usr3/graduate/chenchi/matlab_tool/')
addpath('/projectnb/amazondr/data12/cliveg/chenchi/LST_LAI/Code_CLM5_analysis/')
if strcmp (constants.case_res,'f05_g17')
    this_case_name = 'I2000Clm50Sp.control.climatology.f05_g17';
end
path_in2 = sprintf('/projectnb/amazondr/data12/cliveg/chenchi/LST_LAI/CLM5_TRM_input/%s/%s/',this_case_name,constants.this_season);


% load data
LWin = loadMatData(sprintf('%sFLDS/%s.FLDS.%04d.mat',path_in2,this_case_name,iyear));
ALBEDO = loadMatData(sprintf('%sALBEDO/%s.ALBEDO.%04d.mat',path_in2,this_case_name,iyear));
RA = loadMatData(sprintf('%sRA1/%s.RA1.%04d.mat',path_in2,this_case_name,iyear));
RS = loadMatData(sprintf('%sRS1/%s.RS1.%04d.mat',path_in2,this_case_name,iyear));

SWin = loadMatData(sprintf('%sFSDS/%s.FSDS.%04d.mat',path_in2,this_case_name,iyear));
EMIS = loadMatData(sprintf('%sEMIS/%s.EMIS.%04d.mat',path_in2,this_case_name,iyear));
RHOA = loadMatData(sprintf('%sRHOA/%s.RHOA.%04d.mat',path_in2,this_case_name,iyear));
PS = loadMatData(sprintf('%sPBOT/%s.PBOT.%04d.mat',path_in2,this_case_name,iyear));
TA = loadMatData(sprintf('%sTBOT/%s.TBOT.%04d.mat',path_in2,this_case_name,iyear));
QA = loadMatData(sprintf('%sQBOT/%s.QBOT.%04d.mat',path_in2,this_case_name,iyear));
GHLAND = loadMatData(sprintf('%sGH/%s.GH.%04d.mat',path_in2,this_case_name,iyear));

% mask out the data
mask_rs = isnan(RS);
LWin(mask_rs)=nan;
ALBEDO(mask_rs)=nan;
RA(mask_rs)=nan;
RS(mask_rs)=nan;
SWin(mask_rs)=nan;
EMIS(mask_rs)=nan;
RHOA(mask_rs)=nan;
PS(mask_rs)=nan;
TA(mask_rs)=nan;
QA(mask_rs)=nan;
GHLAND(mask_rs)=nan;


% constants, convert from constants inputs structure
[nrow,ncol]=size(LWin);
fnames = fieldnames(constants);
for ii = 1:length(fnames)
    this_field_name = fnames{ii};
    eval([this_field_name,'=constants.(this_field_name);']);   
end

% substitutions
lambda0 = 1./(4*sigma_sbc.*EMIS.*(TA.^3));
Rn_star = SWin.*(1-ALBEDO)+EMIS.*LWin-sigma_sbc*EMIS.*TA.^4;
delta=get_delta(TA); % e sat at air temperature
R0 = RHOA.*Cp.*lambda0;
gamma = Cp.* PS ./ (0.622*Lv);
ftrm = R0./RA.*(1+delta./gamma.*(RA./(RA+RS)));

% intermedia derivatives
dftrm_dRA = -R0./(RA.^2).*(1+delta./gamma.* (RA./(RA+RS)).^2);
[e_sat] = get_e_sat(TA);
Q_sat_Ta = 0.622.*e_sat./(PS-0.378.*e_sat);

dftrm_dRS = -delta./gamma.*R0./(RA+RS).^2;
dRn_star_dTA = -4*sigma_sbc.*EMIS.*(TA.^3);
dQ_sat_Ta_dTA = 0.622.*PS ./ ((PS-0.378.*e_sat).^2) .*delta; % check it
dR0_dTA = -RHOA.*Cp.*(3/4)./(sigma_sbc.*EMIS.*(TA.^4));
ddelta_dTA =get_delta_dTa(TA);
dftrm_dTA = dR0_dTA./RA .* (1+delta./gamma.*RA./(RA+RS)) + R0./RA .* ddelta_dTA./gamma.*RA./(RA+RS);
dlambda0_dTA = -3/4 ./(sigma_sbc.*EMIS.*(TA.^4));

% final derivatives - sensitivity
del_q = Q_sat_Ta-QA; % a temporary variable

dTs_dALBEDO = -lambda0.*SWin./(1+ftrm);
dTs_dRA = lambda0.*RHOA.*Lv.*(del_q)./(RA+RS).^2./(1+ftrm)  - dftrm_dRA.*lambda0.*(Rn_star-GHLAND-RHOA.*Lv.*(del_q)./(RA+RS))./(1+ftrm).^2;
dTs_dRS = lambda0.*RHOA.*Lv.*(del_q)./(RA+RS).^2./(1+ftrm)  - dftrm_dRS.*lambda0.*(Rn_star-GHLAND-RHOA.*Lv.*(del_q)./(RA+RS))./(1+ftrm).^2;
dTs_dEMIS = lambda0.*( LWin-sigma_sbc.*(TA.^4) - ((Rn_star-GHLAND) - RHOA.*Lv./(RA+RS).*del_q  )./(EMIS.*(1+ftrm)))./(1+ftrm); % corrected on Aug 11
dTs_dSWin = lambda0.*(1-ALBEDO) ./(1+ftrm);
dTs_dLWin = lambda0.*EMIS ./(1+ftrm);
dTs_dQA = RHOA.*Lv./(RA+RS)./(1+ftrm).*lambda0;
dTs_dTA = lambda0.*(dRn_star_dTA-RHOA.*Lv./(RA+RS).*dQ_sat_Ta_dTA)./(1+ftrm) + dlambda0_dTA.*(Rn_star-GHLAND-RHOA.*Lv./(RA+RS).*(del_q))./(1+ftrm) - dftrm_dTA.*lambda0 .*(Rn_star-GHLAND-RHOA.*Lv./(RA+RS).*(del_q))./(1+ftrm).^2 + 1;
dTs_dG = -lambda0./(1+ftrm);

end


%% get saturation vapor pressure
function [es_sat]=get_e_sat(temp)
temp_gr_0 = temp; temp_gr_0(temp<273.15)=nan; % get >=0 
es_sat_gr_0 = 0.61078*exp(17.27*(temp_gr_0-273.15)./(temp_gr_0-273.15+237.3))*1000;

temp_ls_0 = temp; temp_ls_0(temp>=273.15)=nan; % get <0
es_sat_ls_0 = 0.61078*exp(21.875*(temp_ls_0-273.15)./(temp_ls_0-273.15+265.5))*1000;
es_sat = es_sat_gr_0;
es_sat(isnan(es_sat_gr_0))=es_sat_ls_0(isnan(es_sat_gr_0));
end

%% get delta based on es_sat
function [delta]=get_delta(TA)
% syms TA
% esat1 = 0.61078*exp(17.27*(TA-273.15)./(TA-273.15+237.3))*1000;
% esat2 = 0.61078*exp(21.875*(TA-273.15)./(TA-273.15+265.5))*1000;
% diff_T1 = diff(esat1,TA);
% diff_T2 = diff(esat2,TA);

this_TA = TA;this_TA(TA<273.15)=nan; % T>=0 C
delta1 = (30539.*exp(((1727.*this_TA)./100 - 9434601./2000)./(this_TA - 717/20)).*(1727./(100.*(this_TA - 717./20)) - ((1727.*this_TA)./100 - 9434601./2000)./(this_TA - 717./20).^2))./50;

this_TA = TA;this_TA(TA>=273.15)=nan; % T<0 C
delta2 = (30539.*exp(((175.*this_TA)./8 - 191205./32)./(this_TA - 153./20)).*(175./(8*(this_TA - 153./20)) - ((175.*this_TA)./8 - 191205./32)./(this_TA - 153/20).^2))./50;
delta = delta1;
delta(isnan(delta1))=delta2(isnan(delta1));
end


%% get Ddelta/DTa based on es_sat
function [ddelta_dTA]=get_delta_dTa(TA)
% syms TA
% delta1 = (30539*exp(((1727*TA)/100 - 9434601/2000)/(TA - 717/20))*(1727/(100*(TA - 717/20)) - ((1727*TA)/100 - 9434601/2000)/(TA - 717/20)^2))/50;
% delta2 = (30539*exp(((175*TA)/8 - 191205/32)/(TA - 153/20))*(175/(8*(TA - 153/20)) - ((175*TA)/8 - 191205/32)/(TA - 153/20)^2))/50;
% diff_T1 = diff(delta1,TA);
% diff_T2 = diff(delta2,TA);

ddelta_dTA1 = (30539.*exp(((1727.*TA)./100 - 5186726751463539./1099511627776)./(TA - 717./20)).*(1727./(100.*(TA - 717./20)) - ((1727.*TA)./100 - 5186726751463539./1099511627776)./(TA - 717./20).^2).*(1727./(100.*TA - 3585) - ((1727.*TA)./100 - 5186726751463539./1099511627776)./(TA - 717./20).^2))./50 - (30539.*exp(((1727.*TA)./100 - 5186726751463539./1099511627776)./(TA - 717./20)).*(1727./(100.*(TA - 717./20).^2) + 172700./(100.*TA - 3585).^2 - (2.*((1727.*TA)./100 - 5186726751463539./1099511627776))./(TA - 717./20).^3))./50;
ddelta_dTA1(TA<273.15)=nan; % T>=0 C

ddelta_dTA2 = (30539.*exp(((175.*TA)./8 - 191205./32)./(TA - 153./20)).*(175./(8.*(TA - 153./20)) - ((175.*TA)./8 - 191205./32)./(TA - 153./20).^2).*(175./(8.*TA - 306./5) - ((175.*TA)./8 - 191205./32)./(TA - 153./20).^2))./50 - (30539.*exp(((175.*TA)./8 - 191205./32)./(TA - 153./20)).*(175./(8.*(TA - 153./20).^2) + 1400./(8.*TA - 306./5).^2 - (2.*((175.*TA)./8 - 191205./32))./(TA - 153./20).^3))./50;
ddelta_dTA2(TA>=273.15)=nan; % T<0 C

ddelta_dTA = ddelta_dTA1;
ddelta_dTA(isnan(ddelta_dTA1))=ddelta_dTA2(isnan(ddelta_dTA1));
end
