%% CreateVariantFiles.m
%This script makes the variants from spatial correlation maps, excluding
%regions with low signal. 
%Written by Brian Kraus. Edited by Diana Perez.

clear all

%% Paths
%change paths
workbenchdir = '/Applications/workbench/bin_macosx64/';
leftsurf = '/Users/briankraus/Desktop/Quest_Files/Atlases/32k_ConteAtlas_v2_distribute/Conte69.L.midthickness.32k_fs_LR.surf.gii';
rightsurf = '/Users/briankraus/Desktop/Quest_Files/Atlases/32k_ConteAtlas_v2_distribute/Conte69.R.midthickness.32k_fs_LR.surf.gii';
dirpath = '/Users/briankraus/Desktop/';
resttxtname = 'MSC_rest_alltask_samedata_all_varmaps.txt';
tasktxtname = 'MSC_alltask_samedata_all_varmaps.txt';
SNRpath = '/Users/briankraus/Desktop/MSC_SNR_Maps/CIFTI_Files/MSC_Template/Individual_Session_Map/';
outfilepath = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Matched_Task_Data/TaskCat_Split_Half_Data/';

%%
threshold = 2.5;  %% Thresholds used to calculate variants (lowest % or correlation values)
SNRexclusion = 1;  %% Toggles whether to exclude variants based on SNR
MatchedMaps = 1;  %% Toggles whether to use the most task/rest data from the subject with the least amount of data
MatchedCatTasks = 1; %% Toggles whether to use files that are combined and matched to individual tasks
SplitHalf = 1;  %% Toggles whether to use split-half data

%%
    % reads file paths, sub numbers split-halves from txt files    
    [task_files, subjects1, tasksplithalf] = textread([dirpath tasktxtname],'%s%s%s');
    [rest_files, subjects2, restsplithalf] = textread([dirpath resttxtname],'%s%s%s');

% main for-loop, creates variant maps
    for x = 1:length(rest_files)

        subject = subjects2{x};
        
        %reads cifti files from the txt files
        cifti_rest = ft_read_cifti_mod(rest_files{x});
        cifti_task = ft_read_cifti_mod(task_files{x}); 

        %%
        % if you want to exclude low signal regions, this will exclude by SNR
        if SNRexclusion == 1
            SNRmap = ft_read_cifti_mod([SNRpath subject '/' subject '__SNRMap__AllDataConcatenated.dscalar.nii']);
            
            SNRmap.data = SNRmap.data(1:59412,:);
            SNRexclude = find(SNRmap.data < 750);

            cifti_rest.data(SNRexclude,1) = NaN;
            cifti_task.data(SNRexclude,1) = NaN;
        end
        
        %% Makes variant maps
        %this part makes the variant mask by finding the vertices that are
        %below the threshold, making those vertices = 1 and all other
        %vertices = 0
        cifti_task_threshold = find(cifti_task.data < prctile(cifti_task.data,threshold));
        cifti_rest_threshold = find(cifti_rest.data < prctile(cifti_rest.data,threshold));

        cifti_rest_thresh_dat = zeros(size(cifti_rest.data));
        cifti_rest_thresh_dat(cifti_rest_threshold,1) = 1;

        cifti_task_thresh_dat = zeros(size(cifti_task.data));
        cifti_task_thresh_dat(cifti_task_threshold,1) = 1;
        
        %%
        % this does exactly what the six lines above do?
        cifti_rest_final_dat = zeros(size(cifti_rest.data));
        cifti_task_final_dat = zeros(size(cifti_task.data));

        for w = 1:length(cifti_rest.data)

            if cifti_rest_thresh_dat(w) == 1

                cifti_rest_final_dat(w) = 1;

            end

            if cifti_task_thresh_dat(w) == 1

                cifti_task_final_dat(w) = 1;

            end
        end
        %%
        %This creates the output file names for rest and task 
        outfilerest = strrep(rest_files{x}, 'vs_120_allsubs_corr_cortex_corr', ['ThresholdedVariantMap_SNRExclude_' num2str(threshold)]);
        outfiletask = strrep(task_files{x}, 'cortex_vs_120_allsubs_corr_cortex_corr', ['ThresholdedVariantMap_SNRExclude_' num2str(threshold)]);
        
        strrest = ['SNRExclude_' restsplithalf{x}];
        strtask = ['SNRExclude_' tasksplithalf{x}];
        outfilewbtask = [outfilepath subject '/' subject '_matcheddata_Variant_Size_' strtask '_' num2str(threshold) '.dtseries.nii'];
        outfilewbrest = [outfilepath subject '/' subject '_matcheddata_REST_Variant_Size_' strrest '_' num2str(threshold) '.dtseries.nii'];

        %this creates and writes the file in cifti format
        cifti_rest.data = cifti_rest_final_dat;
        cifti_task.data = cifti_task_final_dat;

        ft_write_cifti_mod(outfilerest, cifti_rest)
        ft_write_cifti_mod(outfiletask, cifti_task)

        % This numbers that variants, clustering adjacent vertices that =1
        % together and assigning each cluster a number
        system([workbenchdir 'wb_command -cifti-find-clusters ' outfilerest ' 0 0 0 0 COLUMN ' outfilewbrest ' -left-surface ' leftsurf ' -right-surface ' rightsurf])
        system([workbenchdir 'wb_command -cifti-find-clusters ' outfiletask ' 0 0 0 0 COLUMN ' outfilewbtask ' -left-surface ' leftsurf ' -right-surface ' rightsurf])

    end

