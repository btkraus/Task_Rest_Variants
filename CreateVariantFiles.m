%% CreateVariantFiles.m
%This script makes the variants from spatial correlation maps, excluding
%regions with low signal. 
%Written by Brian Kraus. Edited by Diana Perez.

clear all

%% Paths
%change paths
workbenchdir = '/Applications/workbench/bin_macosx64/';
leftsurf = '/Users/diana/Desktop/Dependencies/32k_ConteAtlas_v2_distribute/Conte69.L.midthickness.32k_fs_LR.surf.gii';
rightsurf = '/Users/diana/Desktop/Dependencies/32k_ConteAtlas_v2_distribute/Conte69.R.midthickness.32k_fs_LR.surf.gii';
dirpath = '/Users/diana/Desktop/Reliability_figure/';
resttxtname = 'MSC_spCorr_Odd_rest.txt';
tasktxtname = 'MSC_spCorr_Odd_task.txt';
SNRpath = '/Users/diana/Desktop/Reliability_figure/SNR_Maps/';
outfilepath = '/Users/diana/Desktop/Reliability_figure/';

%%
threshold = 5;  %% Thresholds used to calculate variants (lowest % or correlation values)
SNRexclusion = 1;  %% Toggles whether to exclude variants based on SNR, 1 = exclude, 0 = don't exclude
ExcludeBySize = 1;

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
            SNRmap = ft_read_cifti_mod([SNRpath '/' subject '_SNRMap_REST_MSCTemplate_AllSessions.4dfp.img_LR_surf_subcort_333_32k_fsLR.dscalar.nii']);
            
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
        outfilewbtask = [outfilepath '/' subject '_matcheddata_Variant_Size_' strtask '_' num2str(threshold) '.dtseries.nii'];
        outfilewbrest = [outfilepath '/' subject '_matcheddata_REST_Variant_Size_' strrest '_' num2str(threshold) '.dtseries.nii'];

        %this creates and writes the file in cifti format
        cifti_rest.data = cifti_rest_final_dat;
        cifti_task.data = cifti_task_final_dat;
        
        ft_write_cifti_mod(outfilerest, cifti_rest)
        ft_write_cifti_mod(outfiletask, cifti_task)
        
        % This numbers that variants, clustering adjacent vertices that =1
        % together and assigning each cluster a number
        system([workbenchdir 'wb_command -cifti-find-clusters ' outfilerest ' 0 0 0 0 COLUMN ' outfilewbrest ' -left-surface ' leftsurf ' -right-surface ' rightsurf])
        system([workbenchdir 'wb_command -cifti-find-clusters ' outfiletask ' 0 0 0 0 COLUMN ' outfilewbtask ' -left-surface ' leftsurf ' -right-surface ' rightsurf])

        % writes file with numbered clusters in preparation for size
        % exclusion
        cifti_rest = ft_read_cifti_mod(outfilewbrest);
        cifti_task = ft_read_cifti_mod(outfilewbtask);

        if ExcludeBySize == 1             
            % exclusion criteria is set to 15 vertices (any variant less
            % than 15 vertices big will be excluded
            [cifti_rest.data, cifti_task.data] = ExcludeVariantSize(cifti_rest.data, cifti_task.data, subject, threshold, 50);
            
        end 
        
        % writes size-excluded variant masks
        ft_write_cifti_mod(outfilerest, cifti_rest)
        ft_write_cifti_mod(outfiletask, cifti_task)
    end

