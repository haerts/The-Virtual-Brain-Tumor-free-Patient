%% TVBii PSE

% Settings
scale = 68;
tm='thr0';
dm='Log';

%% PSE 
main='/home/hannelore/Documents/ANALYSES/TVB_post/subjects_VS/';
cd(main)
sublist=dir('PAT*');
n=length(sublist);
results_fold=('/home/hannelore/Documents/ANALYSES/TVB_post/results_TVBii_VS');
if ~exist(results_fold)
    mkdir(results_fold)
end

for sub=1:n

    %% Load empirical data
            
    subID = sublist(sub).name;
    disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
            
    in_fold = fullfile('/home/hannelore/Documents/ANALYSES/BTC_prepro/subjects/preop', subID);
    out_fold = (['/home/hannelore/Documents/ANALYSES/TVB_post/subjects_VS/' subID]);

    % FC t1
    cd(fullfile(in_fold, 'fmri'));
    load([subID '_fMRI_new.mat']);
    FC_emp = FC_cc_DK68;
    FC_emp = weight_conversion(FC_emp, 'autofix'); %put diagonal to zero
    FC_emp_z=atanh(FC_emp);
    clear PAT* FC_cc* FC_mi ROI*

    % SC t1
    cd(out_fold);
    SCtmp=load(['TVBiiInput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'SC');
    SC=cell2mat(struct2cell(SCtmp));
    clear SCtmp
    
    
    %% Compare empirical structure vs. function
            
    upperIdx = ones([scale scale]);
    upperIdx = triu(upperIdx,1);
    len=length(find(upperIdx));

    FC_emp_zv = reshape(FC_emp_z(~~upperIdx), [len 1]);
    SCv = reshape(SC(~~upperIdx), [len 1]);

    SC_FC_corP = corr(SCv, FC_emp_zv);


    %% Load simulated BOLD 

    params_path=([out_fold '/params']);
    sim_path=([out_fold '/output_' tm '_dist' dm]);
    paramslist=dir([params_path '/param_set_*']);
    m=length(paramslist);
    
    TVBii_PSE = nan(m,9);
    TVBii_Ji = nan(scale,m);

    for index=1:m
        % Load parameters
        params = dlmread([params_path '/param_set_' num2str(index)]);

        % Load simulated TS
        s=dir([sim_path '/BOLD_param_set_' num2str(index) '.txt']);

        if exist([sim_path '/BOLD_param_set_' num2str(index) '.txt']) == 0
            disp('No file found: putting nan in output')
            SC_FC_corP(2)=nan;
            emp_sim_cor(2)=nan;
            tmp_Ji_max = nan(68,1);
            J = mean(tmp_Ji_max);

        elseif s.bytes ~= 0
            tmp_Ji_max = dlmread([sim_path '/BOLD_param_set_' num2str(index) '.txt']);
            tmp_Ji_max = tmp_Ji_max(1:scale,2);
            J = median(tmp_Ji_max);

            TS_sim = dlmread([sim_path '/BOLD_param_set_' num2str(index) '.txt']);
            TS_sim = TS_sim(end-179:end,:);
            FC_sim = corr(TS_sim);
            FC_sim = weight_conversion(FC_sim, 'autofix');
            FC_sim_z = atanh(FC_sim);
            FC_sim_zv = reshape(FC_sim_z(~~upperIdx), [len 1]);

            % Correlate simulated and empirical BOLD
            emp_sim_corP = corr(FC_emp_zv, FC_sim_zv);
            emp_sim_corS = corr(FC_emp_zv, FC_sim_zv, 'Type', 'Spearman');


        elseif s.bytes == 0
            disp('Empty file found: putting nan in output')
            emp_sim_cor=nan;
            tmp_Ji_max = nan(68,1);
            J = mean(tmp_Ji_max);

        end;
          
        % Save output
        TVBii_PSE(index, 1) = params(1); %number of nodes
        TVBii_PSE(index, 2) = params(2); %G
        TVBii_PSE(index, 3) = params(3); %J_NMDA
        TVBii_PSE(index, 4) = params(4); %W+
        TVBii_PSE(index, 5) = J; %median Ji 
        TVBii_PSE(index, 6) = params(9); %TR
        TVBii_PSE(index, 7) = SC_FC_corP; %SC_FCemp correlation 
        TVBii_PSE(index, 8) = emp_sim_corP; %FCemp_FCsim cor (Pearson)
        TVBii_PSE(index, 9) = emp_sim_corS; %FCemp_FCsim cor (Spearman)

        TVBii_Ji(:,index)=tmp_Ji_max;

        clear FC_sim* TS_sim J Ji emp_sim_cor* params s
    end
             
    % Determine optimal G
    [maxcor,maxind]=max(TVBii_PSE(:,8));
    maxG=TVBii_PSE(maxind,2);
    
    % Visually check optimal G value
    cd(results_fold)
    plot(TVBii_PSE(:,2), TVBii_PSE(:,8)); 
    hold on
    plot(TVBii_PSE(:,2), TVBii_PSE(:,9), 'Color', 'm'); 
    hold on
    plot(maxG, maxcor, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10)
    ylabel('Correlation FCemp-FCsim'); xlabel('G'); title(subID);
    legend('Pearson', 'Spearman', 'Pearson max G', 'Location', 'southwest');
    print(['G_FCemp-FCsim_' subID '_VS_pre_' tm], '-dpng')
    close all
    
    % Save output 
    cd(out_fold)
    save(sprintf('TVBiiOutput_%s_scale%d_dist%s',tm,scale,dm), 'FC_emp', 'FC_emp_z', 'FC_emp_zv', 'SC', 'SCv','TVBii_PSE', 'TVBii_Ji');
    clear ans TVBii_PSE TVBii_Ji params* sim_path SC_FC_cor* FC* SC* max* 
    clear in_fold index len out_fold upperIdx m tmp*
            
end 


%% Comparison between TVBii_pre and TVBii_VS results (thrA) -- skip
% Comparison of PSE using SC obtained from MSMT-CSD and SS3T-CSD, both
% using absolute threshold

for sub=1:n
    subID = sublist(sub).name;
    disp(['Processing ' subID]);   
    results_fold_SCorig=(['/home/hannelore/Documents/ANALYSES/TVB_pre/subjects/' subID '/TVBii_' subID '_G-Ji/TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat']);
    results_fold_SCVS=(['/home/hannelore/Documents/ANALYSES/TVB_post/subjects_VS/' subID '/TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat']);

    load(results_fold_SCorig, 'TVBii_PSE'); TVBii_PSE_orig=TVBii_PSE; clear TVBii_PSE
    load(results_fold_SCVS, 'TVBii_PSE'); TVBii_PSE_VS=TVBii_PSE; clear TVBii_PSE
    
    % Visually check optimal G value for both SC construction methods
    cd(results_fold)
    plot(TVBii_PSE_orig(:,2), TVBii_PSE_orig(:,9)); %9=emp_sim_cor_Pearson 
    hold on
    plot(TVBii_PSE_VS(:,2), TVBii_PSE_VS(:,8), 'Color', 'm'); %8=emp_sim_cor_Pearson 
    ylabel('Correlation FCemp-FCsim'); xlabel('G'); title(subID);
    legend('SC ~ MSMT-CSD', 'SC ~ SS3T-CSD', 'Location', 'southeast');
    print(['G_FCemp-FCsim_' subID '_pre_MSMT-SS3T_' tm], '-dpng')
    close all
end



  
%% Save Gmax, SC-FC, FCemp-FCsim, and J 

tm='thrA';
results_fold='/home/hannelore/Documents/ANALYSES/TVB_post/results_TVBii_VS';
TVBii_results=nan(n,5);

for sub=1:n
    subID = sublist(sub).name;
    disp(['Processing ' subID]);      
    out_fold = (['/home/hannelore/Documents/ANALYSES/TVB_post/subjects_VS/' subID]);
    cd(out_fold);

    load([out_fold '/TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm], 'TVBii_PSE');
    results_sorted = sortrows(TVBii_PSE, -8);

    TVBii_results(sub, 1)=sub; %subject number (alphabetical order)
    TVBii_results(sub, 2)=results_sorted(1,7); %SCt1-FCempt1
    TVBii_results(sub, 3)=results_sorted(1,8); %FCempt1-FCsimt1 max corP
    TVBii_results(sub, 4)=results_sorted(1,2); %G
    TVBii_results(sub, 5)=results_sorted(1,5); %median J  
    
    clear TVBii_PSE results_sorted
end 
    
% Save
cd(results_fold)
save(['results_PSE_' tm], 'TVBii_results');

