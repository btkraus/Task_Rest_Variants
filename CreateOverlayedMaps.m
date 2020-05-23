%% Create maps showing overlap or magnitiude differences between states

% This script is for creating maps showing topographic differences in
% overlap or magnitude across states. For creating magnitude maps, a
% variant map (correlation of individual to group) is necessary. For
% applying a size exclusion or creating a map of spatial overlap, a parcel
% map (map with a unique identifier the vertices of each variant) is
% necessary. These maps are read from space delimited text files in the format:
% 'pathtofile subjectID task'. The input data can be excluded by size or by SNR
% mask. Separate input files are necessary for each state.



clear all

%% Set variables for script

SeparateThresholds = 0;  %% Toggles whether to create separate thresholds for task/rest (set to 1)
ShowMagnitudes = 1;  %% Toggles whether to save overlay maps as magnitude differences in variants (set to 1)
SNRExclude = 1;  %% Toggles whether to remove low SNR regions from consideration (set to 1)
SizeExclude = 1;  %% Toggles whether to remove variants below a certain size from consideration (set to 1)
minsize = 50;     %% Minimum size for exclusion of variants (vertices)
threshold = 5;  %% Threshold of variants to create maps from

outputdir = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/';  %% Path for saving overlayed maps
templatepath = '/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Rest_Data/Full_Rest/Full_Data/MSC01/MSC01_REST_AllSessions_vs_120_allsubs_corr_cortex_corr.dtseries.nii';  %% Path where templace cifti is located
SNRMaskPath = '/Users/briankraus/Desktop/MSC_SNR_Maps/CIFTI_Files/MSC_Template/Individual_Session_Map/';  %% Path to SNR masks




if SNRExclude == 1 && SizeExclude == 1 && ShowMagnitudes == 1       %% Text files containing paths to cifti files
    
    % Variant maps
    
    [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_matched_consec_all_varmaps.txt','%s%s%s');
    
    [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_matched_consec_all_varmaps.txt','%s%s%s');
    
    % Parcel maps
    
    [task_masks, sub1, t1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(threshold) '_consec_matched_variants_SNRexclude.txt'],'%s%s%s');
    
    [rest_masks, sub2, t2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(threshold) '_consec_matched_variants_SNRexclude.txt'],'%s%s%s');
    
elseif SNRExclude == 1 && SizeExclude == 1
    
    % Parcel maps
    
    [task_files, subjects1, tasks1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(threshold) '_consec_matched_variants_SNRexclude.txt'],'%s%s%s');
    
    [rest_files, subjects2, tasks2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(threshold) '_consec_matched_variants_SNRexclude.txt'],'%s%s%s');
    
elseif SNRExclude == 1 && ShowMagnitudes == 1
    
    % Variant maps
    
    [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_matched_consec_all_varmaps.txt','%s%s%s');
    
    [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_matched_consec_all_varmaps.txt','%s%s%s');
    
else
    
    % Parcel maps
    
    [task_files, subjects1, tasks1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(threshold) '_consec_matched_variants_SNRexclude.txt'],'%s%s%s');
    
    [rest_files, subjects2, tasks2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(threshold) '_consec_matched_variants_SNRexclude.txt'],'%s%s%s');
    
end


%% Run script


subs = subjects1;   % Get subject IDs

for s = 1:numel(subs)   % Loop through subjects
    
    subject = subs{s};
    task_file = task_files{s};
    rest_file = rest_files{s};
    
    cifti_task = ft_read_cifti_mod(task_file);
    cifti_rest = ft_read_cifti_mod(rest_file);
    
    cifti_rest_thresh_dat = cifti_rest.data;
    cifti_task_thresh_dat = cifti_task.data;
    
    if SizeExclude == 1 && ShowMagnitudes == 1    % Loads parcel files for size exclusion if necessary
        
        task_mask_file = task_masks{s};
        rest_mask_file = rest_masks{s};
        
        task_mask = ft_read_cifti_mod(task_mask_file);
        rest_mask = ft_read_cifti_mod(rest_mask_file);
        
        cifti_task_mask = task_mask.data;
        cifti_rest_mask = rest_mask.data;
        
    end
    
    if SNRExclude == 1      % Performs SNR Exclusion
        
        SNRmask = ft_read_cifti_mod([SNRMaskPath subject '/' subject '__SNRMap__AllDataConcatenated.dscalar.nii']);
        
        SNRmask.data = SNRmask.data(1:59412,:);   % Select cortical SNR only
        
        LowSNR = find(SNRmask.data < 750);      % Find vertices with SNR less than 750
        
    end
    
    if SizeExclude == 1 && ShowMagnitudes == 1      %% Perform size exclusion based on parcel files (magnitude maps only)
        
        allvars_rest = unique(cifti_rest_mask);
        allvars_rest(1) = [];
        allvars_task = unique(cifti_task_mask);
        allvars_task(1) = [];
        
        removevars_rest = [];
        removevars_task = [];
        
        for i = 1:length(allvars_rest)
            
            if length(find(cifti_rest_mask == allvars_rest(i))) < minsize
                
                removevars_rest = [removevars_rest allvars_rest(i)];
                
            end
        end
        
        for h = 1:length(cifti_rest_mask)
            
            if ismember(cifti_rest_mask(h),removevars_rest)
                
                cifti_rest_mask(h) = 0;
                
            elseif cifti_rest_mask(h) > 0
                
                cifti_rest_mask(h) = 1;
                
            end
        end
        
        for i = 1:length(allvars_task)
            
            if length(find(cifti_task_mask == allvars_task(i))) < minsize
                
                removevars_task = [removevars_task allvars_task(i)];
                
            end
        end
        
        for h = 1:length(cifti_task_mask)
            
            if ismember(cifti_task_mask(h),removevars_task)
                
                cifti_task_mask(h) = 0;
                
            elseif cifti_task_mask(h) > 0
                
                cifti_task_mask(h) = 1;
                
            end
        end
        
        nonvariants = find(cifti_rest_mask ~= 1 & cifti_task_mask ~= 1);
        
    else
        
        nonvariants = [];
        
    end
    
    if SeparateThresholds == 1 || ShowMagnitudes == 1   %% Zero out non-variants in magnitude map
        
        if SizeExclude == 1 && ShowMagnitudes == 1
            
            cifti_rest_thresh_dat(nonvariants,:) = 0;
            cifti_task_thresh_dat(nonvariants,:) = 0;
            
        elseif SizeExclude == 1      %% Perform size exclusion based on parcel files (overlap maps only)
            
            allvars_rest = unique(cifti_rest_thresh_dat);
            allvars_rest(1) = [];
            allvars_task = unique(cifti_task_thresh_dat);
            allvars_task(1) = [];
            
            removevars_rest = [];
            removevars_task = [];
            
            for i = 1:length(allvars_rest)
                
                if length(find(cifti_rest_thresh_dat == allvars_rest(i))) < minsize
                    
                    removevars_rest = [removevars_rest allvars_rest(i)];
                    
                end
            end
            
            for h = 1:length(cifti_rest_thresh_dat)
                
                if ismember(cifti_rest_thresh_dat(h),removevars_rest)
                    
                    cifti_rest_thresh_dat(h) = 0;
                    
                elseif cifti_rest_thresh_dat(h) > 0
                    
                    cifti_rest_thresh_dat(h) = 1;
                    
                end
            end
            
            for i = 1:length(allvars_task)
                
                if length(find(cifti_task_thresh_dat == allvars_task(i))) < minsize
                    
                    removevars_task = [removevars_task allvars_task(i)];
                    
                end
            end
            
            for h = 1:length(cifti_task_thresh_dat)
                
                if ismember(cifti_task_thresh_dat(h),removevars_task)
                    
                    cifti_task_thresh_dat(h) = 0;
                    
                elseif cifti_task_thresh_dat(h) > 0
                    
                    cifti_task_thresh_dat(h) = 1;
                    
                end
            end
            
            overlaymapcomp = mean([cifti_rest_thresh_dat'; cifti_task_thresh_dat'],1)';
            
            overlap = find(overlaymapcomp == 1);
            
            restunique = setxor(find(cifti_rest_thresh_dat == 1), overlap);
            
            taskunique = setxor(find(cifti_task_thresh_dat == 1), overlap);
            
        elseif SizeExclude == 1
            
            cifti_rest_thresh_dat(nonvariantsrest,:) = 0;
            cifti_task_thresh_dat(nonvariantstask,:) = 0;
            
            overlaymapcomp = mean([cifti_rest_thresh_dat'; cifti_task_thresh_dat'],1)';
            
            overlap = find(overlaymapcomp == 1);
            
            restunique = setxor(find(cifti_rest_thresh_dat == 1), overlap);
            
            taskunique = setxor(find(cifti_task_thresh_dat == 1), overlap);
            
        else
            
            overlaymapcomp = mean([cifti_rest_thresh_dat'; cifti_task_thresh_dat'],1)';
            
            overlap = find(overlaymapcomp == 1);
            
            restunique = setxor(find(cifti_rest_thresh_dat == 1), overlap);
            
            taskunique = setxor(find(cifti_task_thresh_dat == 1), overlap);
            
        end
        
        if SeparateThresholds == 1      %% Set values for each threshold
            
            combinedoutput = zeros(size(cifti_rest.data));
            combinedoutput(overlap,:) = 1;
            combinedoutput(restunique,:) = .7;
            combinedoutput(taskunique,:) = .5;
            
        elseif ShowMagnitudes == 1      %% Perform substraction for magnitude maps
            
            magnitude = abs((cifti_rest_thresh_dat - cifti_task_thresh_dat));
            %magnitude(nonvariants,:) = NaN;
            
        end
    end
    
    template = ft_read_cifti_mod(templatepath);     %% Read in template file for saving cifti
    template.data = [];
    
    
    if SNRExclude == 1 && SizeExclude == 1 && ShowMagnitudes == 1       %% Save map with filename according to settings used
        
        template.data = magnitude;
        
        ft_write_cifti_mod([outputdir subject 'Task_vs_Rest_Magnitude_Consec_SNRAndSizeExcluded_Overlay_Map_' num2str(threshold) '_Percent'],template)
        
        clear template
        
    elseif SNRExclude == 1 && ShowMagnitudes == 1
        
        template.data = magnitude;
        
        ft_write_cifti_mod([outputdir subject 'Task_vs_Rest_Magnitude_Consec_SNRExcluded_Overlay_Map_' num2str(threshold) '_Percent'],template)
        
        clear template
        
    elseif ShowMagnitudes == 1
        
        template.data = magnitude;
        
        ft_write_cifti_mod([outputdir subject '_Task_vs_Rest_Magnitude_Overlay_Map'],template)
        
        clear template
        
    elseif SeparateThresholds == 1 && SNRExclude == 1 && SizeExclude == 1
        
        template.data = combinedoutput;
        
        ft_write_cifti_mod([outputdir subject 'Task_vs_Rest_Threshold_' num2str(threshold) '_Consec_SNRAndSizeExcluded_Overlay_Map'],template)
        
        clear template
        
    elseif SeparateThresholds == 1 && SNRExclude == 1
        
        template.data = combinedoutput;
        
        ft_write_cifti_mod([outputdir subject '_Task_vs_Rest_Thresholded_SNRExcluded_Overlay_Map'],template)
        
        clear template
        
    elseif SeparateThresholds == 1
        
        template.data = combinedoutput;
        
        ft_write_cifti_mod([outputdir subject '_Task_vs_Rest_Thresholded_Overlay_Map'],template)
        
        clear template
        
    else
        
        overlaymap = [cifti_rest_thresh_dat'; cifti_task_thresh_dat'];
        template.data = mean(overlaymap,1)';
        
        ft_write_cifti_mod([outputdir subject '_Task_vs_Rest_Overlay_Map'],template)
        clear template
        
    end
end

    
    
       