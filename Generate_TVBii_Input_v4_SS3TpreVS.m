
%% ************************* INPUT DATA FOR TVB **************************%
% Version 4: 
% --threshold because otherwise degree distribution isn't decaying
% --normalize because weights cannot be >1 for TVB
%
% SS3TpreVS: pre-operative SC reconstructed using SS3T-CSD, and tracts
% intersecting resection mask removed (ie virtual surgery performed)


%% Prep

% Folders and subjects
prepro_fold='/home/hannelore/Documents/ANALYSES/BTC_prepro/subjects/preop';
TVB_fold='/home/hannelore/Documents/ANALYSES/TVB_post/subjects_VS';
SCchecks_fold='/home/hannelore/Documents/ANALYSES/BTC_prepro/SC_SanityChecks/scale68_SS3T_pre';
cd(TVB_fold);
sublist=dir('PAT*');
m=length(sublist);

% Parameters
dm = 'Log';
scale = 68;
SCversion = 'SCcount_SS3T_5ttnorm_sift1_exclRM';
tm='thrA';


%% Now loop over every subject to adjust TVB input:

sparsity=nan(m,1);
ncomps=nan(m,1);

for index=1:m
    %% Go to subject dwi folder
    subname=sublist(index).name;
    fprintf('Processing %s, dist %s \n\n',subname,dm);
    path2results=fullfile(TVB_fold, subname);
    if ~exist([path2results,'/input'])
        mkdir (path2results, 'input');
    end 
    if ~exist([path2results,'/output'])
        mkdir (path2results, 'output');
    end  
    path2input=fullfile(prepro_fold, subname, 'dwi');
    
    
    %% Load and adapt TVB input
    
    if exist([path2input '/' SCversion '.csv']) == 2
        % Load SC weights
        if exist([path2input '/' SCversion '.mat']) == 0
            SC=dlmread([path2input '/' SCversion '.csv']);
            SC=SC+SC';
            save([path2input '/' SCversion, 'SC']);
        else 
            load([path2input '/' SCversion '.mat']);
        end
                      
        % SC matrix size 
        SCsize = size(SC,1);
           
        % Absolute threshold
        if strcmp(tm,'thrA')==1
            SC=threshold_absolute(SC, 20);
        end
        sparsity(index,1)=sum(find(SC(:))>0) / (SCsize*SCsize);
        % Number of components
        ncomps(index,1)=length(unique(get_components(SC)));
            
        % Histogram of weights
        cd(SCchecks_fold)
        upperIdx = ones([SCsize SCsize]);
        upperIdx = triu(upperIdx,1);
        len=length(find(upperIdx));
        SCv=reshape(SC(~~upperIdx), [len 1]);
        SCvl=log(SCv(SCv>0,:));
        hist(SCvl,20);
        print(['SCweights_' subname '_' tm '_VS'], '-dpng')
        close all
        
        % Degree dist
        degree = degrees_und(SC);
        hist(degree)
        print(['DegreeDist_' subname '_' tm '_VS'], '-dpng')
        close all             
        
        % Normalize weights by dividing by 75K, hence connection  
        % weights represents proportion of connections connecting any 
        % two ROIs + all weights <1
        SC=SC ./ 75000;
            
        % Distances
        SC_dist=-log(SC);
        SC_dist(isinf(SC_dist)) = 0;
        

        %% Save TVB input files
       
        cd(path2results)
        sc_cap_file = [path2results '/input/' subname '_scale' num2str(scale) '_' tm '_dist' dm 'VS_SC_strengths.txt'];
        sc_dist_file = [path2results '/input/' subname '_scale' num2str(scale) '_' tm '_dist' dm 'VS_SC_distances.txt'];
        sc_id_file = [path2results '/input/' subname '_scale' num2str(scale) '_' tm '_dist' dm 'VS_SC_regionids.txt'];
        save( sprintf('TVBiiInput_%s_scale%d_dist%s_VS',tm,scale,dm), 'SC','SC_dist')

        % Write number of nodes into SC files
        dlmwrite(sc_cap_file,SCsize);
        dlmwrite(sc_dist_file,SCsize);
        dlmwrite(sc_id_file,SCsize);

        % Write maximum distance into dist file
        maxdist=max(SC_dist(:));
        dlmwrite(sc_dist_file,maxdist,'delimiter',' ','-append');

        % Write the rest into respective files
        for i = 1:SCsize,
            inpregs=find(SC(i,:)>0);
            inpcaps=SC(i,inpregs);
            inpdists=SC_dist(i,inpregs);

            inpregs=inpregs-1; 
            cap_line=[(i-1) length(inpregs)];
            dist_line=[(i-1) length(inpregs)];
            inp_line=[(i-1) length(inpregs)];
            dlmwrite(sc_cap_file,cap_line,'delimiter',' ','-append');
            dlmwrite(sc_dist_file,dist_line,'delimiter',' ','-append');
            dlmwrite(sc_id_file,inp_line,'delimiter',' ','-append');

            dlmwrite(sc_cap_file,inpcaps,'delimiter',' ','-append','precision','%.8f');
            dlmwrite(sc_dist_file,inpdists,'delimiter',' ','-append','precision','%.8f');
            dlmwrite(sc_id_file,inpregs,'delimiter',' ','-append');
        end
        
    clear SC SC_dist SCsize SCthr* cap_line dd degree* dist_line
    clear i index inp* maxdist sc_* subname sumSC thr upperIdx
     
        
    % In case no SC matrix is available, skip    
    elseif exist([path2input '/SCcount_SS3T_5ttnorm_sift1.mat']) == 0
        disp('No SC matrix found!');
        SCtot(index,:)=nan;
        SCmax(index,:)=nan;
        weights_sparsity_sift1(index,1) = nan;
    end 
end




%% Plot some group stats
cd(SCchecks_fold)

% Number of components
plot(ncomps, '*', 'Color', 'b'); axis([0 m+1 0 3]);
xlabel('subjects'); ylabel('Number of components in SC matrix');
print(['NComps_' tm '_VS'], '-dpng');

% Sparsity of matrices
plot(sparsity, '*', 'Color', 'b'); axis([0 m+1 0.4 0.6]);
xlabel('subjects'); ylabel('Density SC matrix');
print(['Density_' tm '_VS'], '-dpng');
        