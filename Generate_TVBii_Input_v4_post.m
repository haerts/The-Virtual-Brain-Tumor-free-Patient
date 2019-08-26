
%% ************************* INPUT DATA FOR TVB **************************%
% Version 4: thresholded + normalized 
% --threshold SC matrix in 2 ways because otherwise degree distribution
% isn't decaying 
% --normalize because weights cannot be >1 for TVB


%% Prep

% Folders and subjects
TVB_fold='/home/hannelore/Documents/ANALYSES/TVB_post/subjects_post';
cd(TVB_fold);
sublist=dir('*'); sublist=sublist(3:end);
m=length(sublist);

% Parameters
dm = 'Log';
scale = 68;

SCchecks_fold=(['/home/hannelore/Documents/ANALYSES/BTC_prepro/SC_SanityChecks/scale' num2str(scale)]);


%% Now loop over every subject to adjust TVB input:

sparsity=nan(m,1);
sparsity_thrA=nan(m,1);
SCtot=nan(m,1);
SCmax=nan(m,1);
ncomps=nan(m,1);
ncomps_thrA=nan(m,1);

for index=1:m
   
    %% Go to subject dwi folder
    subname=sublist(index).name;
    fprintf('Processing %s, dist %s \n\n',subname,dm);
    path2results=fullfile(TVB_fold, subname, ['TVBii_' subname '_G-Ji']);
    if ~exist(path2results,'file')
        mkdir (path2results);
        mkdir (path2results, 'input');
        mkdir (path2results, 'output');
    end  
    path2input=(['/home/hannelore/Documents/ANALYSES/BTC_prepro/subjects/postop/' subname '/dwi']);

    
    %% Load and adapt TVB input
    
    if exist([path2input '/SCcount_sift1_30M_' num2str(scale)]) == 2
        % Load SC weights
        load([path2input '/SCcount_sift1_30M_' num2str(scale)]);
            
        % SC matrix size 
        SCsize = size(SC,1);
            
        % SC basic stats
        SCtot(index,:)=sum(SC(:)); 
        SCmax(index,:)=max(SC(:));
        sparsity(index,1) = sum(find(SC(:))>0) / (SCsize*SCsize);
            
        % Number of components: should be 1 (all connected), check if not!
        ncomps(index,1)=length(unique(get_components(SC)));
            
        % Histogram of weights: should be nicely distributed
        cd(SCchecks_fold);
        upperIdx = ones([SCsize SCsize]);
        upperIdx = triu(upperIdx,1);
        len=length(find(upperIdx));
        SCv=reshape(SC(~~upperIdx), [len 1]);
        SCvl=log(SCv(SCv>0,:));
        hist(SCvl,20);
        %print(['SCweights_' subname '_sift1_30M'], '-dpng');
        
        % Degree dist: should be decaying, normal is also okay
        degree = degrees_und(SC);
        hist(degree)
        %print(['DegreeDist_' subname '_sift1_30M'], '-dpng')
            %still no powerlaw...
                
            
        % Absolute threshold: adjust (as high as possible, 
		% without making network disconnected)
        SCthrA=threshold_absolute(SC, 5);
        sparsity_thrA(index,1)=sum(find(SCthrA(:))>0) / (SCsize*SCsize);
        ncomps_thrA(index,1)=length(unique(get_components(SCthrA))); 
        SCthrAv=reshape(SCthrA(~~upperIdx), [len 1]);
        SCthrAvl=log(SCthrAv(SCthrAv>0,:));
        hist(SCthrAvl,20); 
        %print(['SCweights_' subname '_sift1_30M_thrA'], '-dpng');
        degreethrA = degrees_und(SCthrA);
        hist(degreethrA)
        %print(['DegreeDist_' subname '_sift1_30M_thrA'], '-dpng')
            
        
        % Normalize weights by dividing by SC max <<across>> subjects
		% Goal: all weights between 0 and 1 (round SCmax)
        SC=SC ./ 75000;
            
        % Distances
        SC_dist=-log(SC);
        SC_dist(isinf(SC_dist)) = 0;
        
        
        %% Save TVB input files
            
        cd(path2results)
        sc_cap_file = [path2results '/input/' subname '_scale' num2str(scale) '_thr0_dist' dm '_SC_strengths.txt'];
        sc_dist_file = [path2results '/input/' subname '_scale' num2str(scale) '_thr0_dist' dm '_SC_distances.txt'];
        sc_id_file = [path2results '/input/' subname '_scale' num2str(scale) '_thr0_dist' dm '_SC_regionids.txt'];
        save( sprintf('TVBiiInput_thrA_scale%d_dist%s',scale,dm), 'SC','SC_dist')

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
        clear len SCv SCvl

        
    % In case no SC matrix is available, skip    
    elseif exist([path2input '/SCcount_sift1_30M_' num2str(scale) '_noACT.mat']) == 0
        disp('No SC matrix found!');
        SCtot(index,:)=nan;
        SCmax(index,:)=nan;
        weights_sparsity_sift1(index,1) = nan;
    end 
end


% Plot some group stats
cd(SCchecks_fold)

% Number of components
plot(ncomps, '*', 'Color', 'b'); axis([0 m+1 0 2]);
hold on
plot(ncomps_thrA, '*', 'Color', 'c');
legend('Unthresholded', 'Absolute thr (SC=5)');
xlabel('subjects'); ylabel('Number of components in SC matrix');
print('NComps_sift1_30M', '-dpng');

% Sparsity of matrices
plot(sparsity, '*','Color', 'b'); axis([0 m+1 0.3 0.8]);
hold on
plot(sparsity_thrA, '*', 'Color', 'c');
legend('Unthresholded', 'Absolute thr (SC=5)');
xlabel('subjects'); ylabel('Density SC matrix');
print('Density_sift1_30M', '-dpng');


% SCtot & max
bar(SCtot); print('SCtot_sift1_30M', '-dpng')
bar(SCmax); xlabel('subjects'), ylabel('SCmax'); 
print('SCmax_sift1_30M', '-dpng');

        
