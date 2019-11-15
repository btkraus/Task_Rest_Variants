%% things to address
%Reliability maps and reliability exclusions


clear all


%%
%change paths
workbenchdir = '/Applications/workbench/bin_macosx64/';
leftsurf = '/Users/briankraus/Desktop/Quest_Files/Atlases/32k_ConteAtlas_v2_distribute/Conte69.L.midthickness.32k_fs_LR.surf.gii';
rightsurf = '/Users/briankraus/Desktop/Quest_Files/Atlases/32k_ConteAtlas_v2_distribute/Conte69.R.midthickness.32k_fs_LR.surf.gii';

%%
CreateVariantFiles = 0;  %% Toggles whether to create files with variant size
threshold = 2.5;  %% Thresholds used to calculate variants (lowest % or correlation values)
ExcludeVariantSize = 1;  %% Toggles whether to exclude variants based on size
SNRexclusion = 1;  %% Toggles whether to exclude variants based on SNR
ReliabilityExclusion = 0;  %% Toggles whether to exclude variants based on low reliability
ExclusionCriteria = 15;  %% Toggles what minimum variant size should be in vertices
FullMaps = 0;  %% Toggles whether to use all the data from task/rest for each subject
MatchedMaps = 1;  %% Toggles whether to use the most task/rest data from the subject with the least amount of data
ConsecMaps = 0;  %% Toggles whether to use maps that are consecutively sampled instead of randomly sampled
IndividTasks = 0; %% Toggles whether to use files that are split into individual tasks
MatchedCatTasks = 1; %% Toggles whether to use files that are combined and matched to individual tasks
SplitHalf = 1;  %% Toggles whether to use split-half data
AbsoluteThresholds = 0;  %% Toggles whether to use absolute or percent thresholds to categorize variants

%%
% reads file paths, sub numbers split-halves from txt files; change paths
if CreateVariantFiles == 1
        
    [task_files, subjects1, tasksplithalf] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_samedata_all_varmaps.txt','%s%s%s');

    [rest_files, subjects2, restsplithalf] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_alltask_all_varmaps.txt','%s%s%s');


% main for-loop, 
    for x = 1:length(rest_files)

        subject = subjects2{x};

        cifti_rest = ft_read_cifti_mod(rest_files{x});
        cifti_task = ft_read_cifti_mod(task_files{x}); 

        if SNRexclusion == 1

            SNRmap = ft_read_cifti_mod(['/Users/briankraus/Desktop/MSC_SNR_Maps/CIFTI_Files/MSC_Template/Individual_Session_Map/' subject '/' subject '__SNRMap__AllDataConcatenated.dscalar.nii']);
            %this is commented out in my version of the script
            %Reliabilitymap = ft_read_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Reliability_Data/WholeBrain_Reliability_Correlations/' subject '_SplitHalf_Reliability_Map_vs_' subject '_All_REST_EvenSessions_cortex_cortex_corr.dtseries.nii']);

            SNRmap.data = SNRmap.data(1:59412,:);
            SNRexclude = find(SNRmap.data < 750);
            %Reliabilityexclude = find(Reliabilitymap.data < .8);

            cifti_rest.data(SNRexclude,1) = NaN;
            cifti_task.data(SNRexclude,1) = NaN;

        end


        cifti_task_threshold = find(cifti_task.data < prctile(cifti_task.data,threshold));
        cifti_rest_threshold = find(cifti_rest.data < prctile(cifti_rest.data,threshold));

        cifti_rest_thresh_dat = zeros(size(cifti_rest.data));
        cifti_rest_thresh_dat(cifti_rest_threshold,1) = 1;

        cifti_task_thresh_dat = zeros(size(cifti_task.data));
        cifti_task_thresh_dat(cifti_task_threshold,1) = 1;

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


        outfilerest = strrep(rest_files{x}, 'vs_120_allsubs_corr_cortex_corr', ['ThresholdedVariantMap_SNRExclude_' num2str(threshold)]);
        outfiletask = strrep(task_files{x}, 'cortex_vs_120_allsubs_corr_cortex_corr', ['ThresholdedVariantMap_SNRExclude_' num2str(threshold)]);
        strrest = ['SNRExclude_' restsplithalf{x}];
        strtask = ['SNRExclude_' tasksplithalf{x}];

        outfilewbtask = ['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Matched_Task_Data/TaskCat_Split_Half_Data/' subject '/' subject '_matcheddata_Variant_Size_' strtask '_' num2str(threshold) '.dtseries.nii'];
        outfilewbrest = ['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Matched_Rest_Data/TaskCat_Split_Half_Data/' subject '/' subject '_matcheddata_REST_Variant_Size_' strrest '_' num2str(threshold) '.dtseries.nii'];

        cifti_rest.data = cifti_rest_final_dat;
        cifti_task.data = cifti_task_final_dat;

        ft_write_cifti_mod(outfilerest, cifti_rest)
        ft_write_cifti_mod(outfiletask, cifti_task)

        system([workbenchdir 'wb_command -cifti-find-clusters ' outfilerest ' 0 0 0 0 COLUMN ' outfilewbrest ' -left-surface ' leftsurf ' -right-surface ' rightsurf])

        system([workbenchdir 'wb_command -cifti-find-clusters ' outfiletask ' 0 0 0 0 COLUMN ' outfilewbtask ' -left-surface ' leftsurf ' -right-surface ' rightsurf])

    end
    
end

%%
if ExcludeVariantSize == 1
    %change paths        
    [task_files, subjects1, thresholds1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_matched_variants_absolute_thresholds_SNRexclude.txt','%s%s%s');

    [rest_files, subjects2, thresholds2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_matched_variants_absolute_thresholds_SNRexclude.txt','%s%s%s');
        
    
    for x = 1:length(rest_files)

        subject = subjects2{x};
        
        cifti_rest = ft_read_cifti_mod(rest_files{x});
        cifti_task = ft_read_cifti_mod(task_files{x});
        
        rest_sizes = [];
        task_sizes = [];
        
        % same values as in cifti_rest.data but w/o repetitions, of
        % vertices?
        vars_rest = unique(cifti_rest.data);
        vars_task = unique(cifti_task.data);
        
        
        %this is counting vertices for what? variants?
        for q = 1:length(vars_rest)
            
            vertcount = 0;
            
            for r = 1:length(cifti_rest.data)
                
                if cifti_rest.data(r) == vars_rest(q)
                    
                    vertcount = vertcount + 1;
                    
                end
            end
            
            rest_sizes = [rest_sizes; vars_rest(q) vertcount];
            
        end
        
        for q = 1:length(vars_task)
            
            vertcount = 0;
            
            for r = 1:length(cifti_task.data)
                
                if cifti_task.data(r) == vars_task(q)
                    
                    vertcount = vertcount + 1;
                    
                end
            end
            
            task_sizes = [task_sizes; vars_task(q) vertcount];
            
        end
        
        %removing vertices belonging to variants that are not at least 15
        %vertices big?
        for i = 2:size(rest_sizes,1)
            
            if rest_sizes(i,2) < ExclusionCriteria
                
                removeverts = find(cifti_rest.data == rest_sizes(i,1));
                
                cifti_rest.data(removeverts,1) = 0;
                
            else
                
                setverts = find(cifti_rest.data == rest_sizes(i,1));
                
                cifti_rest.data(setverts,1) = 1;
                
            end
        end
        
        for i = 2:size(task_sizes,1)
            
            if task_sizes(i,2) < ExclusionCriteria
                
                removeverts = find(cifti_task.data == task_sizes(i,1));
                
                cifti_task.data(removeverts,1) = 0;
                
            else
                
                setverts = find(cifti_task.data == task_sizes(i,1));
                
                cifti_task.data(setverts,1) = 1;
                
            end
        end
        
        
        %sets file name    
        outfilerest = strrep(rest_files{x}, 'SNRExclude', ['SNRExclude_SizeExclude_' num2str(threshold)]);
        outfiletask = strrep(task_files{x}, 'SNRExclude', ['SNRExclude_SizeExclude_' num2str(threshold)]);

        %writes file in cifti format
        ft_write_cifti_mod(outfilerest, cifti_rest)
        ft_write_cifti_mod(outfiletask, cifti_task)
        
        
    end
    
end

