clear

TEEEEEEEEEEEEEEST

%load SPM.mat
[spmmatfile] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
load(spmmatfile);
%set path
path = SPM.swd;
cd(path)

%read Posterior Beta
XYZ  = SPM.xVol.XYZ;
iXYZ = cumprod([1,SPM.xVol.DIM(1:2)'])*XYZ - sum(cumprod(SPM.xVol.DIM(1:2)'));
cB    = spm_data_read(SPM(1).VCbeta,'xyz',XYZ); 

%compute Posterior SD (Taylor approximation)
c=1;
VcB   = c'*SPM.PPM.Cby*c;
l   = spm_data_read(SPM.VHp(1),'xyz',XYZ);
VcB = VcB + (c'*SPM.PPM.dC{1}*c)*(l - SPM.PPM.l(1));
post_SD = sqrt(VcB);                

%Prior SD and Pooled SD
prior_SD = full(sqrt(SPM.PPM.Cb));
pooled_SD = sqrt((post_SD.^2+prior_SD.^2)./2);

%Plot histograms of 1) Posterior SD, 2) Pooled SD, 3) Posterior Beta
figure('Position', [10 500 900 600])
subplot(3,1,1)
hist(post_SD,100)
hold on;
line([prior_SD, prior_SD], ylim, 'LineWidth', 2, 'Color', 'r');
txt = ['\leftarrow Prior SD =' num2str(prior_SD)];
text(prior_SD,500,txt,'FontSize',12,'Color','r')
%xlim([0, 0.1]);        
title('Histogram of Posterior SD', 'FontSize', 12);
subplot(3,1,2)
hist(pooled_SD,100)
hold on;
line([prior_SD, prior_SD], ylim, 'LineWidth', 2, 'Color', 'r');
txt = ['\leftarrow Prior SD =' num2str(prior_SD)];
text(prior_SD,1500,txt,'FontSize',12,'Color','r')
%xlim([0, 0.1]);        
title('Histogram of Pooled SD', 'FontSize', 12);
subplot(3,1,3)
hist(cB,100)
hold on;
line([prior_SD, prior_SD], ylim, 'LineWidth', 2, 'Color', 'r');
txt = ['\leftarrow Prior SD =' num2str(prior_SD)];
text(prior_SD,1500,txt,'FontSize',12,'Color','r')
title('Histogram of Posterior Beta', 'FontSize', 12);

%Effect Size Threshold
figure('Position', [100 80 660 325])
spm_input(['Prior SD =' num2str(prior_SD)],1,'d');
ES_type = spm_input('Effect Size Threshold Type:','+1','Prior SD|Post SD|Pooled SD|PSC',[0,1,2,3],4);
ES_value = spm_input('Effect Size Value',3);
spm_input('Wait...',4,'d');

if ES_type == 0
    ES = ES_value*prior_SD;
    type = 'Prior_SD';
elseif ES_type == 1
    ES = ES_value.*post_SD;
    type = 'Post_SD';
elseif ES_type == 2
    ES = ES_value.*pooled_SD;
    type = 'Pooled_SD';
else
    ES = ES_value;
    type = 'PSC';
end

%compute PPM
PPM_pos_eff = 1 - spm_Ncdf(ES,cB,VcB);
    name_1 = {strrep(['_positive_effect_[' num2str(ES_value) '_' type ']'],'.',',')}; descr_1={['ES threshold = ' num2str(ES_value) ' ' type]};
    
PPM_neg_eff = spm_Ncdf(-ES,cB,VcB);
    name_2 = {strrep(['_negative_effect_[' num2str(ES_value) '_' type ']'],'.',',')}; descr_2={['ES threshold = ' num2str(ES_value) ' ' type]};

PPM_null = 1 - PPM_pos_eff - PPM_neg_eff;
    name_3 = {strrep(['_Null_effect_[' num2str(ES_value) '_' type ']'],'.',',')}; descr_3={['ES threshold = ' num2str(ES_value) ' ' type]};

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
    
%save Posterior SD image
        hdr.fname = [path '\Posterior_SD.nii'];
        hdr.descrip = ['Posterior_SD'];    
        hdr.private.descrip = ['Posterior_SD'];
        tmp           = mask;
        tmp(iXYZ)     = post_SD;
        spm_write_vol(hdr,tmp);
        clear tmp
%save Pooled SD image
        hdr.fname = [path '\Pooled_SD.nii'];
        hdr.descrip = ['Pooled_SD'];    
        hdr.private.descrip = ['Pooled_SD'];
        tmp           = mask;
        tmp(iXYZ)     = pooled_SD;
        spm_write_vol(hdr,tmp);
        clear tmp
%save Cohen's d = post_beta/post_SD
        hdr.fname = [path '\Cohens_d_post_SD.nii'];
        hdr.descrip = ['Cohens d = post_beta/post_SD'];    
        hdr.private.descrip = ['Cohens d = post_beta/post_SD'];
        tmp           = mask;
        tmp(iXYZ)     = cB./post_SD;
        spm_write_vol(hdr,tmp);
        clear tmp
%save Cohen's d = post_beta/pooled_SD
        hdr.fname = [path '\Cohens_d_pooled_SD.nii'];
        hdr.descrip = ['Cohens d = post_beta/pooled_SD'];    
        hdr.private.descrip = ['Cohens d = post_beta/pooled_SD'];
        tmp           = mask;
        tmp(iXYZ)     = cB./pooled_SD;
        spm_write_vol(hdr,tmp);
        clear tmp

%HDI 95% limits     
HDImax = spm_invNcdf(0.975,cB,VcB);
HDImin = spm_invNcdf(0.025,cB,VcB);
%ROPE
ROPE_max = ES;
ROPE_min = -ES;

if ES_type == 0 | ES_type == 3; %For Prior_SD and PSC
    %1) Entire HDI falls within the ROPE (Accept Null)
    for i=1:length(cB)
        if HDImin(1,i)>ROPE_min && HDImax(1,i)<ROPE_max 
            HDI_in_ROPE(1,i) = cB(1,i);
        else
            HDI_in_ROPE(1,i) = NaN;
        end
    end
    
    %2) Entire HDI falls outside the ROPE (Positive Effect)
    for i=1:length(cB)
        if HDImin(1,i)>ROPE_max 
            HDI_outside_ROPE_pos_eff(1,i) = cB(1,i);
        else
            HDI_outside_ROPE_pos_eff(1,i) = NaN;
        end
    end    
    
    %3) Entire HDI falls outside the ROPE (Negative Effect)
    for i=1:length(cB)
        if HDImax(1,i)<ROPE_min 
            HDI_outside_ROPE_neg_eff(1,i) = cB(1,i);
        else
            HDI_outside_ROPE_neg_eff(1,i) = NaN;
        end
    end 
    
    %3) HDI partially overlaps the ROPE
    for i=1:length(cB)
        if ROPE_max<HDImax(1,i) &&  ROPE_max>HDImin(1,i) && ROPE_min<HDImin(1,i) 
            HDI_overlaps_ROPE(1,i) = cB(1,i);
        elseif ROPE_min<HDImax(1,i) &&  ROPE_min>HDImin(1,i) && ROPE_max>HDImax(1,i)
            HDI_overlaps_ROPE(1,i) = cB(1,i);
        else
            HDI_overlaps_ROPE(1,i) = 0;
        end
    end
    
    %4) HDI wider than the ROPE
    for i=1:length(cB)
        if HDImax(1,i)>ROPE_max &&  ROPE_min>HDImin(1,i) 
            HDI_wider_ROPE(1,i) = cB(1,i);
        else
            HDI_wider_ROPE(1,i) = 0;
        end
    end        
else %For Post_SD and Pooled_SD
   %1) Entire HDI falls within the ROPE (Accept Null)
    for i=1:length(cB)
        if HDImin(1,i)>ROPE_min(1,i) && HDImax(1,i)<ROPE_max(1,i) 
            HDI_in_ROPE(1,i) = cB(1,i);
        else
            HDI_in_ROPE(1,i) = NaN;
        end
    end
    
%2) Entire HDI falls outside the ROPE (Positive Effect)
    for i=1:length(cB)
        if HDImin(1,i)>ROPE_max(1,i) 
            HDI_outside_ROPE_pos_eff(1,i) = cB(1,i);
        else
            HDI_outside_ROPE_pos_eff(1,i) = NaN;
        end
    end    
    
 %3) Entire HDI falls outside the ROPE (Negative Effect)
    for i=1:length(cB)
        if HDImax(1,i)<ROPE_min(1,i) 
            HDI_outside_ROPE_neg_eff(1,i) = cB(1,i);
        else
            HDI_outside_ROPE_neg_eff(1,i) = NaN;
        end
    end 
    
 %3) HDI partially overlaps the ROPE
    for i=1:length(cB)
        if ROPE_max(1,i)<HDImax(1,i) &&  ROPE_max(1,i)>HDImin(1,i) && ROPE_min(1,i)<HDImin(1,i) 
            HDI_overlaps_ROPE(1,i) = cB(1,i);
        elseif ROPE_min(1,i)<HDImax(1,i) &&  ROPE_min(1,i)>HDImin(1,i) && ROPE_max(1,i)>HDImax(1,i)
            HDI_overlaps_ROPE(1,i) = cB(1,i);
        else
            HDI_overlaps_ROPE(1,i) = 0;
        end
    end
    
  %4) HDI wider than the ROPE
    for i=1:length(cB)
        if HDImax(1,i)>ROPE_max(1,i) &&  ROPE_min(1,i)>HDImin(1,i) 
            HDI_wider_ROPE(1,i) = cB(1,i);
        else
            HDI_wider_ROPE(1,i) = 0;
        end
    end  
end
    
%5) HDI wider or partially overlaps the ROPE
HDI_overlaps_or_wider_ROPE = HDI_overlaps_ROPE + HDI_wider_ROPE;

%change zeros to NaN for cases 3), 4), 5)
HDI_overlaps_ROPE(HDI_overlaps_ROPE == 0) = NaN;
HDI_wider_ROPE(HDI_wider_ROPE == 0) = NaN;
HDI_overlaps_or_wider_ROPE(HDI_overlaps_or_wider_ROPE == 0) = NaN;
    
%name, description, values
nm_1={['HDI_in_ROPE_[' num2str(ES_value) '_' type ']']}; dscr_1={['Posterior Beta. ROPE: [-' num2str(ES_value) ':+' num2str(ES_value) ']' type]}; 
nm_2={['HDI_outside_ROPE_Positive_Effect_[' num2str(ES_value) '_' type ']']}; dscr_2={['Posterior Beta. ROPE: [-' num2str(ES_value) ':+' num2str(ES_value) ']' type]}; 
nm_3={['HDI_outside_ROPE_Negative_Effect_[' num2str(ES_value) '_' type ']']}; dscr_3={['Posterior Beta. ROPE: [-' num2str(ES_value) ':+' num2str(ES_value) ']' type]}; 
nm_4={['HDI_overlaps_ROPE_[' num2str(ES_value) '_' type ']']}; dscr_4={['Posterior Beta. ROPE: [-' num2str(ES_value) ':+' num2str(ES_value) ']' type]};
nm_5={['HDI_wider_ROPE_[' num2str(ES_value) '_' type ']']}; dscr_5={['Posterior Beta. ROPE: [-' num2str(ES_value) ':+' num2str(ES_value) ']' type]}; 
nm_6={['HDI_overlaps_or_wider_ROPE_[' num2str(ES_value) '_' type ']']}; dscr_6={['Posterior Beta. ROPE: [-' num2str(ES_value) ':+' num2str(ES_value) ']' type]}; 

inf=struct ('name', [nm_1; nm_2; nm_3; nm_4; nm_5; nm_6],...
        'description', [dscr_1; dscr_2; dscr_3; dscr_4; dscr_5; dscr_6]);
    
All=[HDI_in_ROPE; HDI_outside_ROPE_pos_eff; HDI_outside_ROPE_neg_eff; HDI_overlaps_ROPE; HDI_wider_ROPE; HDI_overlaps_or_wider_ROPE];

%save Posterior Beta using HDI+ROPE decision rule
status = exist('HDI_ROPE');
if status == 0
    mkdir 'HDI_ROPE'; 
end

for k=1:6
        hdr.fname = [path '\HDI_ROPE\' inf(k).name '.nii'];
        hdr.descrip = [inf(k).description];    
        hdr.private.descrip = [inf(k).description];
        tmp           = mask;
        tmp(iXYZ)     = All(k,:);
        spm_write_vol(hdr,tmp);
        clear tmp 
end

%done
spm_input('Done',5,'d');

%Plot Histograms of 1)PPM positive effect, 2)PPM negative effect, 3)PPM null effect
figure('Position', [950 50 900 1050])
subplot(3,1,1)
hist(PPM_pos_eff,100)
hold on;
line([0.95, 0.95], ylim, 'LineWidth', 2, 'Color', 'r');
txt = ['PPM =' num2str(0.95) '\rightarrow'];
text(0.81,800,txt,'FontSize',12,'Color','r')
xlim([0, 1]);
ylim([0, 1000]);
title('Histogram of PPM "Positive effect"', 'FontSize', 12);
subplot(3,1,2)
hist(PPM_neg_eff,100)
hold on;
line([0.95, 0.95], ylim, 'LineWidth', 2, 'Color', 'r');
txt = ['PPM =' num2str(0.95) '\rightarrow'];
text(0.81,800,txt,'FontSize',12,'Color','r')
xlim([0, 1]);
ylim([0, 1000]);
title('Histogram of PPM "Negative effect"', 'FontSize', 12);
subplot(3,1,3)
hist(PPM_null,100)
hold on;
line([0.95, 0.95], ylim, 'LineWidth', 2, 'Color', 'r');
txt = ['PPM =' num2str(0.95) '\rightarrow'];
text(0.81,800,txt,'FontSize',12,'Color','r')
xlim([0, 1]);
ylim([0, 1000]);
title('Histogram of PPM "Null effect"', 'FontSize', 12);

%count voxels using HDI+ROPE decision rule
whole_brain = length(cB);
null = nnz(~isnan(HDI_in_ROPE));
pos_eff = nnz(~isnan(HDI_outside_ROPE_pos_eff));
neg_eff = nnz(~isnan(HDI_outside_ROPE_neg_eff));
undecided = whole_brain - null - pos_eff - neg_eff;
%coint voxels using PPM(Posterior Beta>ES)>95% rule
null_2 = sum(PPM_null>0.95);
pos_eff_2 = sum(PPM_pos_eff>0.95);
neg_eff_2 = sum(PPM_neg_eff>0.95);
undecided_2 = whole_brain - null_2 - pos_eff_2 - neg_eff_2;

%Pie Chart
figure('Position', [200 80 800 325])
subplot(1,2,1)
x = [null, pos_eff, neg_eff, undecided];
labels = {'Null','PosEff','NegEff','Undec'};
color = ([0 1 0;      %// green for null
          1 0 0;      %// red for positive
          0 0 1;      %// blue for negative
          1 1 1]);  %// white for undecided
color = color((x~=0),:);
labels = labels((x~=0));
x=x(x~=0);
p = pie(x);
legend(labels,'Location','southoutside','Orientation','horizontal')
title('HDI+ROPE decision rule')
subplot(1,2,2)
x = [null_2, pos_eff_2, neg_eff_2, undecided_2];
labels = {'Null','PosEff','NegEff','Undec'};
labels = labels((x~=0));
x=x(x~=0);
p = pie(x);
legend(labels,'Location','southoutside','Orientation','horizontal')
colormap(color);
title('PPM(Posterior Beta>ES)>95% rule')
set(findobj(gcf,'type','text'),'fontsize',14);

%save Pie Chart
saveas(gcf,['[' num2str(ES_type) ']_' num2str(ES_value) type '.png']);


%FIND BEST EFFECT SIZE THRESHOLD
%ES and Number of Pos/Neg/Null/Undecided voxels dependency
%PPM(Posterior Beta>ES)>95% rule
figure('Position', [100 80 660 325])
spm_input(['Plot dependency between'],1,'d');
spm_input(['ES threshold'],'+1','d');
spm_input(['and'],'+1','d');
spm_input(['Number of Pos/Neg/Null/Undecided voxels?'],'+1','d');
plot_dep = spm_input('','+1','Yes|No',[1,2],2);
if plot_dep == 1
    %compute PPM
    ES1 = spm_input('ES Threshold [min:step:max]','+1');
    for i=1:length(ES1)
        PPM_pos_eff1(i,:) = 1 - spm_Ncdf(ES1(i),cB,VcB);
        PPM_neg_eff1(i,:) = spm_Ncdf(-ES1(i),cB,VcB);
        PPM_null1(i,:) = 1 - PPM_pos_eff1(i,:) - PPM_neg_eff1(i,:);

        null_3(i) = sum(PPM_null1(i,:)>0.95);
        pos_eff_3(i) = sum(PPM_pos_eff1(i,:)>0.95);
        neg_eff_3(i) = sum(PPM_neg_eff1(i,:)>0.95);
        undecided_3(i) = whole_brain - null_3(i) - pos_eff_3(i) - neg_eff_3(i);
        %percent of voxels
        null_3(i) = 100*null_3(i)/whole_brain;
        pos_eff_3(i) = 100*pos_eff_3(i)/whole_brain;
        neg_eff_3(i) = 100*neg_eff_3(i)/whole_brain;
        undecided_3(i) = 100*undecided_3(i)/whole_brain;
    end    
    figure('Position', [10 500 900 600])
    plot(ES1,undecided_3,ES1,null_3,ES1,pos_eff_3,ES1,neg_eff_3,'LineWidth',2)
    title('PPM(Posterior Beta>ES)>95% rule')
    xlabel('ES threshold, BOLD signal %') 
    ylabel('% of voxels') 
    legend('Undecided','Null','Positive Effect','Negative Effect')
    ax = gca;
    ax.FontSize = 13;
    ay = gca;
    ay.FontSize = 13;
    %save
    saveas(gca,['ESthreshold_dependency.png']);
end
%clear

