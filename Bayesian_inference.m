%========================================================================
%Script was used in "Non-selective response inhibition in Go/NoGo task:
%Bayesian analysis of fMRI data" paper

%Bayesian parameter inference (Friston & Penny, 2003; Penny & Ridgway, 2013)
%The "ROPE-only" decision rule (Kruschke, 2018)

%Masharipov Ruslan, october, 2019
%Institute of Human Brain of RAS, St. Petersburg, Russia
%Neuroimaging lab
%masharipov@ihb.spb.ru
%========================================================================

% Before running the script use SPM12 to:
% 1) Create one-sample t-test at the second level for contrast e.g. (Cond A - Cond B)
% 2) Estimate model using method: Classical
% 3) Estimate model using method: Bayesian 2nd-level
% 4) Load SPM.mat
%========================================================================


%set path
path = SPM.swd;
cd(path)

%read Posterior Beta
XYZ  = SPM.xVol.XYZ;
iXYZ = cumprod([1,SPM.xVol.DIM(1:2)'])*XYZ - sum(cumprod(SPM.xVol.DIM(1:2)'));
cB    = spm_data_read(SPM(1).VCbeta,'xyz',XYZ); 

%compute Posterior Variance
%choose contrast c = +1 or c = -1
c=1;
VcB   = c'*SPM.PPM.Cby*c;
for j = 1:length(SPM.PPM.l)
    l   = spm_data_read(SPM.VHp(j),'xyz',XYZ);
    VcB = VcB + (c'*SPM.PPM.dC{j}*c)*(l - SPM.PPM.l(j));
end
%post_SD = sqrt(VcB);                

%prior SD
prior_SD = full(sqrt(c'*SPM.PPM.Cb*c));

%choose effect size threshold
%default SPM12 threshold:
%one standard deviation of the prior variance of the contrast (Friston & Penny, 2003)  
ES = prior_SD;
type = 'PSC';

%compute PPM
PPM_pos_eff = 1 - spm_Ncdf(ES,cB,VcB);
    name_1 = {strrep(['_positive_effect_[' num2str(ES) '_' type ']'],'.',',')};
    descr_1={['ES threshold = ' num2str(ES) ' ' type]};
    
PPM_neg_eff = spm_Ncdf(-ES,cB,VcB);
    name_2 = {strrep(['_negative_effect_[' num2str(ES) '_' type ']'],'.',',')};
    descr_2={['ES threshold = ' num2str(ES) ' ' type]};

PPM_null = 1 - PPM_pos_eff - PPM_neg_eff;
    name_3 = {strrep(['_Null_effect_[' num2str(ES) '_' type ']'],'.',',')};
    descr_3={['ES threshold = ' num2str(ES) ' ' type]};

%compute Log Posterior Odds
Log_Post_Odds_pos_eff = log( (PPM_pos_eff+eps) ./ (1 - PPM_pos_eff+eps) ); 
Log_Post_Odds_neg_eff = log( (PPM_neg_eff+eps) ./ (1 - PPM_neg_eff+eps) );                  
Log_Post_Odds_null = log( (PPM_null+eps) ./ (1 - PPM_null+eps) );

%name, description, values
info=struct ('name', [name_1 name_2 name_3],...
             'description', [descr_1 descr_2 descr_3]);

PPM_all=[PPM_pos_eff; PPM_neg_eff; PPM_null];
Log_Post_Odds_all=[Log_Post_Odds_pos_eff; Log_Post_Odds_neg_eff; Log_Post_Odds_null];
    
%hdr
hdr = spm_vol([path '\Cbeta_0001.nii']);

%mask
mask = spm_read_vols(hdr);
mask(~isnan(mask)) = 0;

%save PPM
for j=1:3
        hdr.fname = [path '\PPM' info(j).name '.nii'];
        hdr.descrip = [info(j).description];    
        hdr.private.descrip = [info(j).description];
        tmp           = mask;
        tmp(iXYZ)     = PPM_all(j,:);
        spm_write_vol(hdr,tmp);
        clear tmp 
end
    
%save Log Posterior Odds
for j=1:3
        hdr.fname = [path '\Log_Post_Odds' info(j).name '.nii'];
        hdr.descrip = [info(j).description];    
        hdr.private.descrip = [info(j).description];
        tmp           = mask;
        tmp(iXYZ)     = Log_Post_Odds_all(j,:);
        spm_write_vol(hdr,tmp);
        clear tmp         
end
    
clear
