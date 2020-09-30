function A002_dBio_to_dLAI()
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
    this_case_name = 'I2000Clm50Sp.control.climatology.f05_g17';
end
path_in2 = sprintf('/projectnb/amazondr/data12/cliveg/chenchi/LST_LAI/CLM5_TRM_input/%s/%s/',this_case_name,this_season);
LWin = loadMatData(sprintf('%sFLDS/%s.FLDS.%04d.mat',path_in2,this_case_name,2000));

% start calculation
nyear = constants.eyear - constants.syear +1;
icount=0;
[nrow,ncol]=size(LWin);
dALBEDO_dLAI_stack = nan(nrow,ncol,nyear);
dRA_dLAI_stack = nan(nrow,ncol,nyear);
dRS_dLAI_stack = nan(nrow,ncol,nyear);
dEMIS_dLAI_stack = nan(nrow,ncol,nyear);
dG_dLAI_stack = nan(nrow,ncol,nyear);

for iyear = constants.syear:constants.eyear
    icount=icount+1;
    % get LST sensitivities
    [dALBEDO_dLAI,dRA_dLAI,dRS_dLAI,dEMIS_dLAI,dG_dLAI]=get_Para_LAI_sensitivities(iyear,constants);
    dALBEDO_dLAI_stack(:,:,icount)=dALBEDO_dLAI;
    dRA_dLAI_stack(:,:,icount)=dRA_dLAI;
    dRS_dLAI_stack(:,:,icount)=dRS_dLAI;
    dEMIS_dLAI_stack(:,:,icount)=dEMIS_dLAI;
    dG_dLAI_stack(:,:,icount)=dG_dLAI;
    check=1;
end

%% median sensitivity 2000-2017
output.dALBEDO_dLAI_median = nanmedian(dALBEDO_dLAI_stack,3);
output.dRA_dLAI_median = nanmedian(dRA_dLAI_stack,3);
output.dRS_dLAI_median = nanmedian(dRS_dLAI_stack,3);
output.dEMIS_dLAI_median = nanmedian(dEMIS_dLAI_stack,3);
output.dG_dLAI_median = nanmedian(dG_dLAI_stack,3);
savename = sprintf('%sPARA2LAI_median_sensitivity_%04d_%04d_%s.mat',path_out,constants.syear,constants.eyear,case_res);
save(savename,'output','-v7.3');
clear output

end


%% get other parameter sensitivty to LAI
function [dALBEDO_dLAI_out,dRA_dLAI_out,dRS_dLAI_out,dEMIS_dLAI_out,dG_dLAI_out]=get_Para_LAI_sensitivities(iyear,constants)
addpath('/usr3/graduate/chenchi/matlab_tool/')
addpath('/projectnb/amazondr/data12/cliveg/chenchi/LST_LAI/Code_CLM5_analysis/')

% constants, convert from constants inputs structure

fnames = fieldnames(constants);
for ii = 1:length(fnames)
    this_field_name = fnames{ii};
    eval([this_field_name,'=constants.(this_field_name);']);   
end

if strcmp (constants.case_res,'f05_g17')
    control_case = 'I2000Clm50Sp.control.climatology.f05_g17';
    sensitivies_case = {'I2000Clm50Sp.sensitivity098.f05_g17','I2000Clm50Sp.sensitivity102.f05_g17'};
end

conrol_path_in = sprintf('/projectnb/amazondr/data12/cliveg/chenchi/LST_LAI/CLM5_TRM_input/%s/%s/',control_case,constants.this_season);
for ii = 1:length(sensitivies_case)
    sen_path_in{ii} = sprintf('/projectnb/amazondr/data12/cliveg/chenchi/LST_LAI/CLM5_TRM_input/%s/%s/',sensitivies_case{ii},constants.this_season);
end

% load control data
ALBEDO = loadMatData(sprintf('%sALBEDO/%s.ALBEDO.%04d.mat',conrol_path_in,control_case,iyear));
RA = loadMatData(sprintf('%sRA1/%s.RA1.%04d.mat',conrol_path_in,control_case,iyear));
RS = loadMatData(sprintf('%sRS1/%s.RS1.%04d.mat',conrol_path_in,control_case,iyear));
TLAI = loadMatData(sprintf('%sTLAI/%s.TLAI.%04d.mat',conrol_path_in,control_case,iyear));
EMIS = loadMatData(sprintf('%sEMIS/%s.EMIS.%04d.mat',conrol_path_in,control_case,iyear));
GH = loadMatData(sprintf('%sGH/%s.GH.%04d.mat',conrol_path_in,control_case,iyear));

[nrow,ncol]=size(LWin);
dALBEDO_dLAI = nan(nrow,ncol,length(sensitivies_case));
dRA_dLAI = nan(nrow,ncol,length(sensitivies_case));
dRS_dLAI = nan(nrow,ncol,length(sensitivies_case));
dEMIS_dLAI = nan(nrow,ncol,length(sensitivies_case));
dG_dLAI = nan(nrow,ncol,length(sensitivies_case));
icount=0;
for ii=1:length(sensitivies_case)
    icount=icount+1;
    % load sensitivites data
    ALBEDO_sen = loadMatData(sprintf('%sALBEDO/%s.ALBEDO.%04d.mat',sen_path_in{ii},sensitivies_case{ii},iyear));
    RA_sen = loadMatData(sprintf('%sRA1/%s.RA1.%04d.mat',sen_path_in{ii},sensitivies_case{ii},iyear));
    RS_sen = loadMatData(sprintf('%sRS1/%s.RS1.%04d.mat',sen_path_in{ii},sensitivies_case{ii},iyear));
    TLAI_sen = loadMatData(sprintf('%sTLAI/%s.TLAI.%04d.mat',sen_path_in{ii},sensitivies_case{ii},iyear));
    EMIS_sen = loadMatData(sprintf('%sEMIS/%s.EMIS.%04d.mat',sen_path_in{ii},sensitivies_case{ii},iyear));
    GH_sen = loadMatData(sprintf('%sGH/%s.GH.%04d.mat',sen_path_in{ii},sensitivies_case{ii},iyear));
    
    dALBEDO_dLAI(:,:,icount) = (ALBEDO_sen-ALBEDO)./(TLAI_sen-TLAI);
    dRA_dLAI(:,:,icount) = (RA_sen-RA)./(TLAI_sen-TLAI);
    dRS_dLAI(:,:,icount) = (RS_sen-RS)./(TLAI_sen-TLAI);
    dEMIS_dLAI(:,:,icount) = (EMIS_sen-EMIS)./(TLAI_sen-TLAI);
    dG_dLAI(:,:,icount) = (GH_sen-GH)./(TLAI_sen-TLAI);  
end
dALBEDO_dLAI_out = nanmean(dALBEDO_dLAI,3);
dRA_dLAI_out = nanmean(dRA_dLAI,3);
dRS_dLAI_out = nanmean(dRS_dLAI,3);
dEMIS_dLAI_out = nanmean(dEMIS_dLAI,3);
dG_dLAI_out = nanmean(dG_dLAI,3);
end

