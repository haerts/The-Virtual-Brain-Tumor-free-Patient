%% TVBii PSE

% Settings
scale = 68;
tm='thrA';
dm='Log';

%% Loop over all subjects
main='/home/hannelore/Documents/ANALYSES/TVB_post/subjects/';
sublist=dir(main); sublist=sublist(3:end-1,:);
n=length(sublist);
results_fold='/home/hannelore/Documents/ANALYSES/TVB_post/results_TVBii';

for sub=1:n

    %% Load empirical data
            
    subID = sublist(sub).name;
    subID_t1 = [subID(1:5), 'T1'];
    disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
            
    in_fold = fullfile('/home/hannelore/Documents/ANALYSES/BTC_prepro/subjects/postop', subID);
    in_fold_t1 = fullfile('/home/hannelore/Documents/ANALYSES/BTC_prepro/subjects/preop', subID_t1);
    out_fold = (['/home/hannelore/Documents/ANALYSES/TVB_post/subjects/' subID '/TVBii_' subID '_G-Ji']);

    % FC t2
    cd(fullfile(in_fold, 'fmri'));
    load([subID '_fMRI_new.mat']);
    FC_emp = FC_cc_DK68;
    FC_emp = weight_conversion(FC_emp, 'autofix'); %put diagonal to zero
    FC_emp_z=atanh(FC_emp);
    clear CON* PAT* FC_cc* FC_mi ROI*

    % SC t2
    cd(out_fold);
    SCtmp=load(['TVBiiInput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], ['SC' tm 'n']);
    SC=cell2mat(struct2cell(SCtmp));
    clear SCtmp
    
    % SC t1
    cd([in_fold_t1 '/dwi'])
    SCtmp=load(['SC_sift1_30M_scale' num2str(scale) '_thrAn_dist' dm '.mat'], ['SC' tm 'n']);
    SC_t1=cell2mat(struct2cell(SCtmp));
    clear SCtmp


    %% Compare empirical structure vs. function
            
    upperIdx = ones([scale scale]);
    upperIdx = triu(upperIdx,1);
    len=length(find(upperIdx));

    FC_emp_zv = reshape(FC_emp_z(~~upperIdx), [len 1]);
    SCv = reshape(SC(~~upperIdx), [len 1]);
    SCv_t1 = reshape(SC_t1(~~upperIdx), [len 1]);

    SC_FC_corP = corr(SCv, FC_emp_zv);
    SCt1_FC_corP = corr(SCv_t1, FC_emp_zv);


    %% Load simulated BOLD 

    params_path=([out_fold '/params']);
    sim_path=([out_fold '/output_' tm '_dist' dm]);
    paramslist=dir([params_path '/param_set_*']);
    m=length(paramslist);
    
    TVBii_PSE = nan(m,10);
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
        TVBii_PSE(index, 5) = J; %median Ji (TVB pre this was mean Ji)
        TVBii_PSE(index, 6) = params(9); %TR
        TVBii_PSE(index, 7) = SC_FC_corP; %SCt2_FCempt2 correlation 
        TVBii_PSE(index, 8) = SCt1_FC_corP; %SCt1_FCempt2 correlation 
        TVBii_PSE(index, 9) = emp_sim_corP; %FCempt2_FCsimt2 cor (Pearson)
        TVBii_PSE(index, 10) = emp_sim_corS; %FCempt2_FCsimt2 cor (Spearman)

        TVBii_Ji(:,index)=tmp_Ji_max;

        clear FC_sim* TS_sim J Ji emp_sim_cor* params s
    end
             
    % Determine optimal G
    [maxcor,maxind]=max(TVBii_PSE(:,9));
    maxG=TVBii_PSE(maxind,2);
    
    % Visually check optimal G value
    cd(results_fold)
    plot(TVBii_PSE(:,2), TVBii_PSE(:,9)); 
    hold on
    plot(TVBii_PSE(:,2), TVBii_PSE(:,10))%, 'Color', 'm'); 
    hold on
    plot(maxG, maxcor, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10)
    ylabel('Correlation FCemp-FCsim'); xlabel('G'); title(subID);
    legend('Pearson', 'Spearman', 'Pearson max G', 'Location', 'southwest');
    print(['G_FCemp-FCsim_' subID], '-dpng')
    close all
    
    % Save output 
    cd(out_fold)
    save(sprintf('TVBiiOutput_%s_scale%d_dist%s',tm,scale,dm), 'FC_emp', 'FC_emp_z', 'FC_emp_zv', 'SC', 'SCv', 'SC_t1', 'SCv_t1', 'TVBii_PSE', 'TVBii_Ji');
    clear TVBii_PSE TVBii_Ji params* sim_path SC_FC_cor* FC* SC* 
    clear in_fold index len out_fold upperIdx m
            
end 

    
%% Save Gmax, SC-FC, FCemp-FCsim, and J 

results_fold='/home/hannelore/Documents/ANALYSES/TVB_post/results_TVBii';
TVBii_results=nan(n,6);

for sub=1:n
    subID = sublist(sub).name;
    disp(['Processing ' subID]);      
    out_fold = (['/home/hannelore/Documents/ANALYSES/TVB_post/subjects/' subID '/TVBii_' subID '_G-Ji']);
    cd(out_fold);

    load([out_fold '/TVBiiOutput_thrA_scale68_distLog'], 'TVBii_PSE');
    results_sorted = sortrows(TVBii_PSE, -9);

    TVBii_results(sub, 1)=sub; %subject number (alphabetical order)
    TVBii_results(sub, 2)=results_sorted(1,7); %SCt2-FCempt2
    TVBii_results(sub, 3)=results_sorted(1,8); %SCt1-FCempt2
    TVBii_results(sub, 4)=results_sorted(1,9); %FCempt2-FCsimt2 max corP
    TVBii_results(sub, 5)=results_sorted(1,2); %G
    TVBii_results(sub, 6)=results_sorted(1,5); %median J  
    
    clear TVBii_PSE results_sorted
end 
    
% Save
cd(results_fold)
save('results_G-J', 'TVBii_results');


%% Plot some results

% Model fit: SCt1-FCempt2, SCt2-FCempt2, FCsimt2-FCempt2
plot(TVBii_results(:,1), TVBii_results(:,4), 'v', 'Color', 'c', 'MarkerFaceColor', 'c')
axis([0 n+1 0 0.65]);
hold on
plot(TVBii_results(:,1), TVBii_results(:,2), '^', 'Color', 'm', 'MarkerFaceColor', 'm'); 
hold on
plot(TVBii_results(:,1), TVBii_results(:,3), '^', 'Color', 'b', 'MarkerFaceColor', 'b')
legend({'FCsim t2 - FCemp t2','SC t2 - FCemp t2', 'SC t1 - FCemp t2'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
legend('boxoff')
xlabel('subjects'); ylabel('Pearson correlation coefficient'); 
%print('ModelFit_SCt1-SCt2-FCsimt2', '-dpng')

% Optimal G
plot(TVBii_results(:,5), '*')
xlabel('Subjects'), ylabel('G'); axis([0 n+1 1 2.6])
%print('G', '-dpng');

% Optimal J
plot(TVBii_results(:,6), '*')
xlabel('Subjects'), ylabel('median J'); axis([0 n+1 1.15 1.65])
%print('Jmd', '-dpng');

% G-J
plot(TVBii_results(:,5), TVBii_results(:,6), '*');
xlabel('G'), ylabel('J'); 
%print('G_Jmd', '-dpng');



%% Check Ji distribution 

for index=1:n
    % Load data
    subID=sublist(index).name;
    disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
    load([main subID '/TVBii_' subID '_G-Ji/TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'TVBii_Ji', 'TVBii_PSE')

    % Identify max
    [maxval, maxind]=max(TVBii_PSE(:,9));
    Ji_max=TVBii_Ji(:,maxind);
    hist(Ji_max)
    save([main subID '/TVBii_' subID '_G-Ji/TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'Ji_max', 'maxind', '-append')
               
    clear Ji_max TVBii* maxval maxind
end
%--> distribution is skewed, more values around one, some larger


%% Median Ji in tumor vs. non-tumor regions
% In controls, copy whole-brain Ji to JiT & JiNT

% Tumor nodes
tumornodes_fold='/home/hannelore/Documents/ANALYSES/TVB_post/tumorrois_scale68.csv';
tumornodes=csvread(tumornodes_fold);
clear tumornodes_fold

for sub=1:n    
    % First get Ji values corresponding to max FCsim-FCemp
    subID = sublist(sub).name;
    disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
    load([main subID '/TVBii_' subID '_G-Ji/TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'Ji_max')
        
    % Now identify tumor and nontumor nodes, and get median Ji in these
    % regions
    sub_tumornodes = find(tumornodes(:,sub));
    Ji_tumornodes = Ji_max(sub_tumornodes,:);
    J_tumornodes_md(sub,1) = median(Ji_tumornodes);
    
    sub_nontumornodes = find(tumornodes(:,sub)==0);
    Ji_nontumornodes = Ji_max(sub_nontumornodes,:);
    J_nontumornodes_md(sub,1) = median(Ji_nontumornodes); 
    
    clear sub_* Ji*
end


%% Median Ji after regressing out SC strength -- wrong! regress out region size instead!

for sub=2:n
   subID = sublist(sub).name;
   disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
   
   % Load Ji data and SC strength values for optimal G value
   load([main subID '/TVBii_' subID '_G-Ji/TVBiiOutput_thrA_scale68_distLog.mat'], 'maxind', 'Ji_max');
   results_fold=([main subID '/TVBii_' subID '_G-Ji/output_thrA_distLog']);
   J = dlmread([results_fold '/BOLD_param_set_' num2str(maxind) '.txt']); 
   tmp_Ji_max = J(1:scale,2);
   tmp_Ji_max_doublecheck = Ji_max - tmp_Ji_max;
   disp(['    Difference between Ji max and Ji check is ' num2str(mean(tmp_Ji_max_doublecheck))])
   clear tmp*
 
   Ji_strength = J(1:scale,1);
   B=regress(Ji_max, Ji_strength);
   res(:,sub)=Ji_max-B*Ji_strength;   
   
   % Median JiBrain res
   J_res(sub)=median(res(:,sub));
   
   % Median JiTumor res
   sub_tumornodes = find(tumornodes(:,sub));
   JiT_res(sub) = median(res(sub_tumornodes,sub));
   
   % Median JiNonTumor res
   sub_nontumornodes = find(tumornodes(:,sub)==0);
   JiNT_res(sub) = median(res(sub_nontumornodes,sub));
   
   clear J B Ji_* maxind sub_*
end
    
J_res=J_res';
JiT_res=JiT_res';
JiNT_res=JiNT_res';
      

%% Median Ji after regressing out region size!

load('/home/hannelore/Documents/ANALYSES/TVB_post/ROIsizes.mat', 'RoiSize')

for sub=2:n
   subID = sublist(sub).name;
   disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
   
   % Load Ji data for optimal G value and region size
   load([main subID '/TVBii_' subID '_G-Ji/TVBiiOutput_thrA_scale68_distLog.mat'], 'maxind', 'Ji_max');
   results_fold=([main subID '/TVBii_' subID '_G-Ji/output_thrA_distLog']);
   J = dlmread([results_fold '/BOLD_param_set_' num2str(maxind) '.txt']); 
   tmp_Ji_max = J(1:scale,2);
   tmp_Ji_max_doublecheck = Ji_max - tmp_Ji_max;
   disp(['    Difference between Ji max and Ji check is ' num2str(mean(tmp_Ji_max_doublecheck))])
   clear tmp*
 
   sub_roisize = RoiSize(:,sub);
   B=regress(Ji_max, sub_roisize);
   res(:,sub)=Ji_max-B*sub_roisize;   
   
   % Median JiBrain res
   J_res(sub)=median(res(:,sub));
   
   % Median JiTumor res
   sub_tumornodes = find(tumornodes(:,sub));
   JiT_res(sub) = median(res(sub_tumornodes,sub));
   
   % Median JiNonTumor res
   sub_nontumornodes = find(tumornodes(:,sub)==0);
   JiNT_res(sub) = median(res(sub_nontumornodes,sub));
   
   clear J B Ji_* maxind sub_*
end
    
J_res=J_res';
JiT_res=JiT_res';
JiNT_res=JiNT_res';
