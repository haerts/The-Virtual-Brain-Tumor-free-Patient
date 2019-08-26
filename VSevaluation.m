%% Virtual surgery evaluation
% Evaluate performance of TVBii prediction of post-surgical brain dynamics
% after virtual surgery

% Settings
scale = 68;
tm='thr0';
dm='Log';

main='/home/hannelore/Documents/ANALYSES/TVB_post/subjects_VS/';
cd(main)
sublist=dir('PAT*');
n=length(sublist);
results_fold=('/home/hannelore/Documents/ANALYSES/TVB_post/results_TVBii_VS');
cd(results_fold)
load(['results_PSE_' tm]) 
    %col1=subID nr
    %col2=correlation between SC (pre-op, reconstructed with SS3T-CSD)
    %       and FCemp (pre-op)
    %col3=correlation between FCsim (pre-op, based on SC SS3T-CSD) and
    %       FCemp (pre-op)
    %col4=optimal G pre-op
    %col5=optimal mean J pre-op

for sub=1:7
    
    subID = sublist(sub).name;
    subIDpost = ([subID(1:5) 'T2']);
    disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
    
    in_fold = fullfile('/home/hannelore/Documents/ANALYSES/BTC_prepro/subjects/postop', subIDpost);
    out_fold = (['/home/hannelore/Documents/ANALYSES/TVB_post/subjects_VS/' subID]);

    
    %% Evaluation of parameters
    %--> was optimal J used or optimized again for this SC?
    
    sim_path=([out_fold '/output']);
    tmp_Ji_max = dlmread([sim_path '/BOLD_param_set_OptimalPre_' tm '.txt']);
    tmp_Ji_max = tmp_Ji_max(1:scale,2);
    J = mean(tmp_Ji_max);
    

    %% Load data

    % FCemp post-surgery
    cd(fullfile(in_fold, 'fmri'));
    load([subIDpost '_fMRI_new.mat']);
    FC_emp = FC_cc_DK68;
    FC_emp = weight_conversion(FC_emp, 'autofix'); %put diagonal to zero
    FC_emp_z=atanh(FC_emp);
    clear PAT* FC_cc* FC_mi ROI*

    %SC pre-surgery, after virtual surgery
    cd(out_fold);
    SCtmp=load(['TVBiiInput_' tm '_scale' num2str(scale) '_dist' dm '_VS.mat'], 'SC');
    SC=cell2mat(struct2cell(SCtmp));
    clear SCtmp
    
    %FCsim after virtual surgery
    TS_sim = dlmread([sim_path '/BOLD_param_set_OptimalPre_' tm '.txt']);
    TS_sim = TS_sim(end-179:end,:);
    FC_sim = corr(TS_sim);
    FC_sim = weight_conversion(FC_sim, 'autofix');
    FC_sim_z = atanh(FC_sim);
    
    
    
    %% Compare empirical vs. predicted FC
    
    upperIdx = ones([scale scale]);
    upperIdx = triu(upperIdx,1);
    len=length(find(upperIdx));
    
    FC_sim_zv = reshape(FC_sim_z(~~upperIdx), [len 1]);
    FC_emp_zv = reshape(FC_emp_z(~~upperIdx), [len 1]);
    SCv = reshape(SC(~~upperIdx), [len 1]);
            
    % Correlate simulated and empirical BOLD
    FCemp_FCsim_corP = corr(FC_emp_zv, FC_sim_zv);
    
    % Correlate structure (pre, after VS) with FCemp post-op
    SC_FC_corP = corr(SCv, FC_emp_zv);

    
    %% Save output
    TVBii_results(sub, 6) = J; 
        %mean J 
    TVBii_results(sub, 7) = SC_FC_corP; 
        %correlation between SC (pre, after virtual surgery) and FCemp 
        %       post-surgery 
    TVBii_results(sub, 8) = FCemp_FCsim_corP; 
        %correlation between FCsim (pre, after virtual surgery) and FCemp 
        %       post-surgery 
    

    clear FC* SC* TS_sim J subID* tmp* 
end

cd(results_fold)
save(['results_PSE_' tm], 'TVBii_results')



      
