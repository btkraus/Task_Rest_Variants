
clear all

AcrossSubjects = 0;  %% Toggles whether to overlay maps across subjects
WithinSubjects = 1;  %% Toggles whether to overlap maps within subjects
MSCreference = 0;  %% Toggles whether to use MSC referenced data or 120 subject referenced data
SeparateMaps = 0; %% Toggles whether to create separate maps for task/rest
SeparateThresholds = 1;  %% Toggles whether to create separate thresholds for task/rest
ShowMagnitudes = 0;  %% Toggles whether to save overlay maps as magnitude differences in variants
GroupSNRMasks = 0;  %% Toggles whether to exclude areas of low SNR via group SNR masks
IndividSNRMasks = 1;  %% Toggles whether to exclude areas of low SNR via individual SNR masks
IndividtSNRMasks = 0;  %% Toggles whether to exclude areas of low tSNR via individual SNR masks
ReliabilityMasks = 0;  %% Toggles whether to exclude areas of low reliability via reliability maps
SNRExclude = 1;  %% Toggles whether to use data that has already had low SNR removed
ReliableExclude = 0;  %% Toggles whether to use data that has already had low reliability removed
SizeExclude = 1;  %% Toggles whether to use data that has already had small variants removed
MatchedMaps = 1;  %% Toggles whether to use all data within a subject or matched amounts of data

AbsoluteThresholds = 0;  %% Toggles whether to use absolute value or percent thresholds for data
thresholds = [2.5 10];  %% Thresholds to create overlap maps from
usedeciles = 0;     %% Toggles whether to calculate deciles
matrixmaps = 0;     %% Toggles whether to create maps for comparison with matching matrix

% deciles = [0 10;      %% Percentage Deciles
%            10 20;
%            20 30;
%            30 40;
%            40 50;
%            50 60;
%            60 70;
%            70 80;
%            80 90;
%            90 100];
       
% deciles = [-1 .1;       %% Absolute Deciles
%            .1 .2;
%            .2 .3;
%            .3 .4;
%            .4 .5;
%            .5 .6;
%            .6 .7;
%            .7 1];
       
       
for g = 1:numel(thresholds)
%for g = 1:size(deciles,1)
    
    threshold = thresholds(g);
    
    if usedeciles == 1
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_even_matched.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_even_matched.txt','%s%s%s');
        
    elseif matrixmaps == 1 && SizeExclude == 1
        
        [task_files, subjects1, tasks1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(thresholds(g)) '_matched_variants_SNRexclude.txt'],'%s%s%s');

        [rest_files, subjects2, tasks2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(thresholds(g)) '_matched_variants_SNRexclude.txt'],'%s%s%s');
    
        [task_masks, sub1, t1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(thresholds(g)) '_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');

        [rest_masks, sub2, t2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(thresholds(g)) '_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');
        
    elseif matrixmaps == 1
        
        [task_files, subjects1, tasks1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(thresholds(g)) '_matched_variants_SNRexclude.txt'],'%s%s%s');

        [rest_files, subjects2, tasks2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(thresholds(g)) '_matched_variants_SNRexclude.txt'],'%s%s%s');
        
    elseif SNRExclude == 1 && SizeExclude == 1 && ShowMagnitudes == 1 && MatchedMaps == 1 && numel(thresholds) > 1
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_matched_rest.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_matched_task.txt','%s%s%s');
    
        [task_masks, sub1, t1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(thresholds(g)) '_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');

        [rest_masks, sub2, t2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(thresholds(g)) '_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');
        
    elseif SNRExclude == 1 && ShowMagnitudes == 1 && MatchedMaps == 1 && numel(thresholds) > 1
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_matched_rest.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_matched_task.txt','%s%s%s');
        
    elseif SNRExclude == 1 && ReliableExclude == 1 && SizeExclude == 1
    
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps_variants_bothexclude_sizeexclude.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps_variants_bothexclude_sizeexclude.txt','%s%s%s');
        
    elseif SNRExclude == 1 && SizeExclude == 1 && numel(thresholds) > 1 && AbsoluteThresholds == 1
        
        [task_files, subjects1, threshold1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/Absolute_Thresholds/MSC_alltask_varmaps_' num2str(threshold) '.txt'],'%s%s%s');

        [rest_files, subjects2, threshold2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/Absolute_Thresholds/MSC_rest_varmaps_' num2str(threshold) '.txt'],'%s%s%s');
    
        [task_masks, sub1, t1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/Absolute_Thresholds/MSC_alltask_varmaps_variants_SNRexclude_sizeexclude_' num2str(threshold) '.txt'],'%s%s%s');

        [rest_masks, sub2, t2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/Absolute_Thresholds/MSC_rest_varmaps_variants_SNRexclude_sizeexclude_' num2str(threshold) '.txt'],'%s%s%s');
        
    elseif SNRExclude == 1 && SizeExclude == 1 && numel(thresholds) > 1
        
        [task_files, subjects1, threshold1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/MSC_alltask_varmaps_' num2str(threshold) '.txt'],'%s%s%s');

        [rest_files, subjects2, threshold2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/MSC_rest_varmaps_' num2str(threshold) '.txt'],'%s%s%s');
    
        [task_masks, sub1, t1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/MSC_alltask_varmaps_variants_SNRexclude_sizeexclude_' num2str(threshold) '.txt'],'%s%s%s');

        [rest_masks, sub2, t2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/MSC_rest_varmaps_variants_SNRexclude_sizeexclude_' num2str(threshold) '.txt'],'%s%s%s');
        
    elseif SNRExclude == 1 && numel(thresholds) > 1 && AbsoluteThresholds == 1
        
        [task_files, subjects1, threshold1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/Absolute_Thresholds/MSC_alltask_varmaps_' num2str(threshold) '.txt'],'%s%s%s');

        [rest_files, subjects2, threshold2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/Absolute_Thresholds/MSC_rest_varmaps_' num2str(threshold) '.txt'],'%s%s%s');
        
    elseif SNRExclude == 1 && numel(thresholds) > 1
        
        [task_files, subjects1, threshold1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/MSC_alltask_varmaps_' num2str(threshold) '.txt'],'%s%s%s');

        [rest_files, subjects2, threshold2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/Thresholded_Maps/MSC_rest_varmaps_' num2str(threshold) '.txt'],'%s%s%s');
        
    elseif SNRExclude == 1 && SizeExclude == 1 && ShowMagnitudes == 1 && MatchedMaps == 1
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_matched.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_matched.txt','%s%s%s');
    
        [task_masks, sub1, t1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_matched_variants_SNRexclude_sizeexclude.txt','%s%s%s');

        [rest_masks, sub2, t2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_matched_variants_SNRexclude_sizeexclude.txt','%s%s%s');
    
    elseif SNRExclude == 1 && SizeExclude == 1 && ShowMagnitudes == 1
    
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps.txt','%s%s%s');
    
        [task_masks, sub1, t1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps_variants_SNRexclude_sizeexclude.txt','%s%s%s');

        [rest_masks, sub2, t2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps_variants_SNRexclude_sizeexclude.txt','%s%s%s');
        
    elseif SNRExclude == 1 && SizeExclude == 1 && MatchedMaps == 1
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_matched_variants_SNRexclude_sizeexclude.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_matched_variants_SNRexclude_sizeexclude.txt','%s%s%s');
    
    elseif SNRExclude == 1 && SizeExclude == 1
    
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps_variants_SNRexclude_sizeexclude.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps_variants_SNRexclude_sizeexclude.txt','%s%s%s');
    
    elseif ReliableExclude == 1 && SizeExclude == 1 && ShowMagnitudes == 1
    
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps.txt','%s%s%s');
    
        [task_masks, sub1, t1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps_variants_Reliableexclude_sizeexclude.txt','%s%s%s');

        [rest_masks, sub2, t2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps_variants_Reliableexclude_sizeexclude.txt','%s%s%s');
    
    elseif ReliableExclude == 1 && SizeExclude == 1
    
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps_variants_Reliableexclude_sizeexclude.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps_variants_Reliableexclude_sizeexclude.txt','%s%s%s');
    
    elseif SizeExclude == 1 && ShowMagnitudes == 1
    
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps.txt','%s%s%s');
    
        [task_masks, sub1, t1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_variants_SizeOnly_sizeexclude.txt','%s%s%s');

        [rest_masks, sub2, t2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_variants_SizeOnly_sizeexclude.txt','%s%s%s');
    
    elseif SNRExclude == 1 && ReliableExclude == 1
    
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps_BothSNRexclude.txt','%s%s%s');
        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps_BothSNRexclude.txt','%s%s%s');
        
    elseif SNRExclude == 1 && ShowMagnitudes == 1
        
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps.txt','%s%s%s');
        
    elseif SNRExclude == 1
    
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_SNRexclude.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_SNRexclude.txt','%s%s%s');
    
    elseif ReliableExclude == 1
    
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps_Reliableexclude.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps_Reliableexclude.txt','%s%s%s');
    
    elseif SizeExclude == 1
    
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps_variants_SizeOnly_sizeexclude.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps_variants_SizeOnly_sizeexclude.txt','%s%s%s');

    elseif MSCreference == 1
    
        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps_MSCRef.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps_MSCRef.txt','%s%s%s');
    
    else

        [task_files, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_alltask_varmaps.txt','%s%s%s');

        [rest_files, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/MSC_rest_varmaps.txt','%s%s%s');

    end


    subs = unique(subjects1);

    if AcrossSubjects == 1
    
        AllSubsTask = [];
        AllSubsRest = [];
    
    end
    
    for s = 1:numel(subs)
    
    
        task_filenums = find(strcmp(subjects1, subs{s}))';    %% Gets file/task indices for current subject
        rest_filenums = find(strcmp(subjects2, subs{s}))';    %% Gets file/task indices for current subject
        
        for t = task_filenums
        
            subject = subs{s};
            task_file = task_files{t};
            rest_file = rest_files{t};
            %task = tasks1{t};
            %task2 = tasks2{t};
        
            cifti_task = ft_read_cifti_mod(task_file);
            cifti_rest = ft_read_cifti_mod(rest_file);
        
            if (ShowMagnitudes == 1 && (SizeExclude == 1)) || (numel(thresholds) > 1 && SizeExclude == 1)
            
                task_mask_file = task_masks{t};
                rest_mask_file = rest_masks{t};
            
                cifti_task_mask = ft_read_cifti_mod(task_mask_file);
                cifti_rest_mask = ft_read_cifti_mod(rest_mask_file);
            
            end
        
            if GroupSNRMasks == 1
    
                SNRmask = ft_read_cifti_mod('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/bottomBrainMask.dtseries.nii'); 
            
                LowSNR = find(SNRmask.data == 1);
            
            elseif IndividSNRMasks == 1
            
                SNRmask = ft_read_cifti_mod(['/Users/briankraus/Desktop/MSC_SNR_Maps/CIFTI_Files/MSC_Template/Individual_Session_Map/' subject '/' subject '__SNRMap__AllDataConcatenated.dscalar.nii']); 
            
                SNRmask.data = SNRmask.data(1:59412,:);
                
                LowSNR = find(SNRmask.data < 750);
            
            elseif IndividtSNRMasks == 1
            
                SNRmask = ft_read_cifti_mod(['/Users/briankraus/Desktop/MSC_SNR_Maps/CIFTI_Files/MSC_Template/Individual_Session_Map/' subject '/' subject '__tSNRMap__AllDataConcatenated.dscalar.nii']); 
            
                LowSNR = find(SNRmask.data < 25);
            
            elseif ReliabilityMasks == 1
            
                Reliabilitymask = ft_read_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Reliability_Data/WholeBrain_Reliability_Correlations/' subject '_SplitHalf_Reliability_Map_vs_' subject '_All_REST_EvenSessions_cortex_cortex_corr.dtseries.nii']);
            
                LowSNR = find(Reliabilitymask.data <= .8);
            
            end
        
        
            if (GroupSNRMasks + IndividSNRMasks + IndividtSNRMasks + ReliabilityMasks) == 1
                
                %if usedeciles == 1
                    
                    cifti_rest.data(LowSNR,:) = NaN;
                    cifti_task.data(LowSNR,:) = NaN;
                    
                %else
        
                    %cifti_task_thresh_dat(LowSNR,:) = NaN;
                    %cifti_rest_thresh_dat(LowSNR,:) = NaN;
                    
                %end
            
            end
        
        %end       % End of task loop if only one file per subject
        
        if usedeciles == 1
            
            if AbsoluteThresholds == 1
                
                cifti_rest_thresh = deciles(g,:);
                cifti_task_thresh = deciles(g,:);
                
            else
            
                cifti_rest_thresh = prctile(cifti_rest.data,deciles(g,:));
                cifti_task_thresh = prctile(cifti_task.data,deciles(g,:));
            
            end
            
            if g == 1
                
                cifti_rest_remove = find(cifti_rest.data >= cifti_rest_thresh(2));
                cifti_task_remove = find(cifti_task.data >= cifti_task_thresh(2));
                 
            elseif g > 1 && g < size(deciles,1)
            
                cifti_rest_remove = find(cifti_rest.data < cifti_rest_thresh(1) | cifti_rest.data >= cifti_rest_thresh(2));
                cifti_task_remove = find(cifti_task.data < cifti_task_thresh(1) | cifti_task.data >= cifti_task_thresh(2));

            elseif g == size(deciles,1)
                
                cifti_rest_remove = find(cifti_rest.data < cifti_rest_thresh(1));
                cifti_task_remove = find(cifti_task.data < cifti_task_thresh(1));
                
            end
            
            cifti_rest_thresh_dat = cifti_rest.data;
            cifti_task_thresh_dat = cifti_task.data;
            
            
    
        elseif ((SizeExclude == 1) && ShowMagnitudes == 1) || (numel(thresholds) > 1 && SizeExclude == 1)
        
            cifti_rest_mask = cifti_rest_mask.data;
            cifti_task_mask = cifti_task_mask.data;
        
            cifti_rest_thresh_dat = cifti_rest.data;
            cifti_task_thresh_dat = cifti_task.data;
    
        elseif SizeExclude == 1 || SNRExclude == 1 || ReliableExclude == 1
        
            cifti_rest_thresh_dat = cifti_rest.data;
            cifti_task_thresh_dat = cifti_task.data;
            
        end
    
        if WithinSubjects == 1
            
            if SizeExclude == 1 && (numel(thresholds) > 1 || matrixmaps == 1) && ShowMagnitudes == 0
                    
            	nonvariantsrest = find(cifti_rest_mask ~= 1); 
                nonvariantstask = find(cifti_task_mask ~= 1);
                
            elseif ShowMagnitudes == 1
                
                nonvariants = find(cifti_rest_mask ~= 1 & cifti_task_mask ~= 1);
                
            else
                
                nonvariants = [];
                    
            end
        
            if (SeparateThresholds == 1 || ShowMagnitudes == 1)
                
                if usedeciles == 1
                
                    cifti_rest_thresh_dat(cifti_rest_remove,:) = NaN;
                    cifti_task_thresh_dat(cifti_task_remove,:) = NaN;
                
                    overlaymapcomp = mean([cifti_rest_thresh_dat'; cifti_task_thresh_dat'],1)';
                    
                    overlap = find(isnan(overlaymapcomp) == 0 & overlaymapcomp ~= 0);
                    
                    restunique = setxor(find(isnan(cifti_rest_thresh_dat) == 0 & cifti_rest_thresh_dat ~= 0), overlap);
            
                    taskunique = setxor(find(isnan(cifti_task_thresh_dat) == 0 & cifti_task_thresh_dat ~= 0), overlap);
                    
                
                elseif SizeExclude == 1 && ShowMagnitudes == 1
                    
                    cifti_rest_thresh_dat(nonvariants,:) = 0;
                    cifti_task_thresh_dat(nonvariants,:) = 0;
                    
                elseif SizeExclude == 1
                
                    cifti_rest_thresh_dat(nonvariantsrest,:) = 0;
                    cifti_task_thresh_dat(nonvariantstask,:) = 0;
                    
                    if matrixmaps == 1
                        
                        replacevalsrest = find(cifti_rest_thresh_dat > 0);
                        replacevalstask = find(cifti_task_thresh_dat > 0);
                        
                        cifti_rest_thresh_dat(replacevalsrest) = 1;
                        cifti_task_thresh_dat(replacevalstask) = 1;
                        
                    end
            
                    overlaymapcomp = mean([cifti_rest_thresh_dat'; cifti_task_thresh_dat'],1)';
            
                    overlap = find(overlaymapcomp == 1);
            
                    restunique = setxor(find(cifti_rest_thresh_dat == 1), overlap);
            
                    taskunique = setxor(find(cifti_task_thresh_dat == 1), overlap);
                    
                else
                    
                    if matrixmaps == 1
                        
                        replacevalsrest = find(cifti_rest_thresh_dat > 0);
                        replacevalstask = find(cifti_task_thresh_dat > 0);
                        
                        cifti_rest_thresh_dat(replacevalsrest) = 1;
                        cifti_task_thresh_dat(replacevalstask) = 1;
                        
                    end
            
                    overlaymapcomp = mean([cifti_rest_thresh_dat'; cifti_task_thresh_dat'],1)';
            
                    overlap = find(overlaymapcomp == 1);
            
                    restunique = setxor(find(cifti_rest_thresh_dat == 1), overlap);
            
                    taskunique = setxor(find(cifti_task_thresh_dat == 1), overlap);

                end
                    
                if SeparateThresholds == 1
            
                    combinedoutput = zeros(size(cifti_rest.data));
                    combinedoutput(overlap,:) = 1;
                    combinedoutput(restunique,:) = .7;
                    combinedoutput(taskunique,:) = .5;
                
                elseif ShowMagnitudes == 1

                    magnitude = abs((cifti_rest_thresh_dat - cifti_task_thresh_dat));
                    %magnitude(nonvariants,:) = NaN;
                
                
                else
                
                    magnitude = abs((cifti_rest.data - cifti_task.data));
                    %magnitude(nonvariants,:) = NaN;
                
                end
            
            end
        
            template = ft_read_cifti_mod('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Rest_Data/Full_Rest/Full_Data/MSC01/MSC01_REST_AllSessions_vs_120_allsubs_corr_cortex_corr.dtseries.nii');
            template.data = [];
        
        
            if SeparateMaps == 1 && MSCreference == 1 && GroupSNRMasks == 1
            
                template.data = cifti_task_thresh_dat;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_Only_GroupSNRMask_Overlay_MSCRef_Map'],template)
            
                template.data = cifti_rest_thresh_dat;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Rest_Only_SNRMask_Overlay_MSCRef_Map'],template)
            
                clear template
        
            elseif SeparateMaps == 1 && MSCreference == 1
            
                template.data = cifti_task_thresh_dat;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_Only_Overlay_MSCRef_Map'],template)
            
                template.data = cifti_rest_thresh_dat;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Rest_Only_Overlay_MSCRef_Map'],template)
            
                clear template
            
            elseif SeparateMaps == 1 && GroupSNRMasks == 1
            
                template.data = cifti_task_thresh_dat;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_Only_GroupSNRMask_Overlay_Map'],template)
            
                template.data = cifti_rest_thresh_dat;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Rest_Only_SNRMask_Overlay_Map'],template)
            
                clear template    
             
            elseif SeparateMaps == 1
            
                template.data = cifti_task_thresh_dat;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_Only_Overlay_Map'],template)
            
                template.data = cifti_rest_thresh_dat;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Rest_Only_Overlay_Map'],template)
            
                clear template
                
            elseif SNRExclude == 1 && SizeExclude == 1 && ShowMagnitudes == 1 && MatchedMaps == 1 && numel(thresholds) > 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject 'MatchedData_Task_vs_Rest_Magnitude_SNRAndSizeExcluded_Overlay_Map_'  num2str(threshold) '_Percent'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && MSCreference == 1 && GroupSNRMasks == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_GroupSNRMask_Overlay_MSCRef_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && MSCreference == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_Overlay_MSCRef_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && SNRExclude == 1 && ReliableExclude == 1 && SizeExclude == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_BothSNRAndSizeExcluded_Overlay_Map'],template)
            
                clear template
                
            elseif SNRExclude == 1 && SizeExclude == 1 && ShowMagnitudes == 1 && MatchedMaps == 1
                
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject 'MatchedData_Task_vs_Rest_Magnitude_SNRAndSizeExcluded_Overlay_Map_' num2str(threshold) '_Percent'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && SNRExclude == 1 && SizeExclude == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_SNRAndSizeExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && ReliableExclude == 1 && SizeExclude == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_ReliabilityAndSizeExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && SizeExclude == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_SizeExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && SNRExclude == 1 && ReliableExclude == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_BothSNRExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && SNRExclude == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_SNRExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && ReliableExclude == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_ReliabilityExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && GroupSNRMasks == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_GroupSNRMask_Overlay_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && ReliabilityMasks == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_ReliabilityMask_Overlay_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && IndividSNRMasks == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_IndividSNRMask_Overlay_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1 && IndividtSNRMasks == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_IndividtSNRMask_Overlay_Map'],template)
            
                clear template
            
            elseif ShowMagnitudes == 1
            
                template.data = magnitude;
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Magnitude_Overlay_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && MSCreference == 1 && GroupSNRMasks == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_GroupSNRMask_Overlay_MSCRef_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && MSCreference == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_Overlay_MSCRef_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && SNRExclude == 1 && ReliableExclude == 1 && SizeExclude == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_BothSNRAndSizeExcluded_Overlay_Map'],template)
            
                clear template
                
            elseif SeparateThresholds == 1 && SNRExclude == 1 && SizeExclude == 1 && matrixmaps == 1
                
                template.data = combinedoutput;
                
             	ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_BothSNRAndSizeExcluded_Overlay_Map_' num2str(thresholds(g)) '_pct_MatrixMatch'],template)
            
             	clear template
                
            elseif SeparateThresholds == 1 && SNRExclude == 1 && matrixmaps == 1
                
                template.data = combinedoutput;
                
             	ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_SNRExcluded_Overlay_Map_' num2str(thresholds(g)) '_pct_MatrixMatch'],template)
            
             	clear template
                
            elseif SeparateThresholds == 1 && SNRExclude == 1 && usedeciles == 1 && AbsoluteThresholds == 1
                
                template.data = combinedoutput;
                
             	ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/Decile_Maps/' subject '_Task_vs_Rest_Thresholded_SNRExcluded_EvenSessions_Overlay_Map_' num2str(deciles(g,2)) '_Decile_Abs'],template)
            
             	clear template
                
            elseif SeparateThresholds == 1 && SNRExclude == 1 && usedeciles == 1
                
                template.data = combinedoutput;
                
             	ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/Decile_Maps/' subject '_Task_vs_Rest_Thresholded_SNRExcluded_EvenSessions_Overlay_Map_' num2str(deciles(g,2)) '_Decile'],template)
            
             	clear template
                
            elseif SeparateThresholds == 1 && SNRExclude == 1 && SizeExclude == 1 && numel(thresholds) > 1 && AbsoluteThresholds == 1
                
                if mod(t,2) == 1
                
                    template.data = combinedoutput;
            
                    ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/Thresholded_Maps/Absolute_Thresholds/SNR_Size_Excluded/' subject '_Task_vs_Rest_Thresholded_SNRAndSizeExcluded_EvenSessions_Overlay_Map_' num2str(threshold) '_Abs'],template)
            
                    clear template
                
                else
                    
                    template.data = combinedoutput;
            
                    ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/Thresholded_Maps/Absolute_Thresholds/SNR_Size_Excluded/' subject '_Task_vs_Rest_Thresholded_SNRAndSizeExcluded_OddSessions_Overlay_Map_' num2str(threshold) '_Abs'],template)
            
                    clear template
                    
                end
                
            elseif SeparateThresholds == 1 && SNRExclude == 1 && SizeExclude == 1 && numel(thresholds) > 1
                
                if mod(t,2) == 1
                
                    template.data = combinedoutput;
            
                    ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/Thresholded_Maps/' subject '_Task_vs_Rest_Thresholded_SNRAndSizeExcluded_EvenSessions_Overlay_Map_' num2str(threshold) '_Percent'],template)
            
                    clear template
                
                else
                    
                    template.data = combinedoutput;
            
                    ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/Thresholded_Maps/' subject '_Task_vs_Rest_Thresholded_SNRAndSizeExcluded_OddSessions_Overlay_Map'],template)
            
                    clear template
                    
                end
                
            elseif SeparateThresholds == 1 && SNRExclude == 1 && SizeExclude == 1 && MatchedMaps == 1
                
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject 'MatchedData_Task_vs_Rest_Thresholded_SNRAndSizeExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && SNRExclude == 1 && SizeExclude == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_SNRAndSizeExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && ReliableExclude == 1 && SizeExclude == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_ReliabilityAndSizeExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && SizeExclude == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_SizeExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && SNRExclude == 1 && ReliableExclude == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_BothSNRExcluded_Overlay_Map'],template)
            
                clear template
                
            elseif SeparateThresholds == 1 && SNRExclude == 1 && numel(thresholds) > 1 && AbsoluteThresholds == 1
                
                if mod(t,2) == 1
                
                    template.data = combinedoutput;
            
                    ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/Thresholded_Maps/Absolute_Thresholds/SNR_Excluded/' subject '_Task_vs_Rest_Thresholded_SNRExcluded_EvenSessions_Overlay_Map_' num2str(threshold) '_Abs'],template)
            
                    clear template
                    
                else
                    
                    template.data = combinedoutput;
            
                    ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/Thresholded_Maps/Absolute_Thresholds/SNR_Excluded/' subject '_Task_vs_Rest_Thresholded_SNRExcluded_OddSessions_Overlay_Map_' num2str(threshold) '_Abs'],template)
            
                    clear template
                    
                end
                
            elseif SeparateThresholds == 1 && SNRExclude == 1 && numel(thresholds) > 1
                
                if mod(t,2) == 1
                
                    template.data = combinedoutput;
            
                    ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/Thresholded_Maps/' subject '_Task_vs_Rest_Thresholded_SNRExcluded_EvenSessions_Overlay_Map_' num2str(threshold) '_Percent'],template)
            
                    clear template
                    
                else
                    
                    template.data = combinedoutput;
            
                    ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/Thresholded_Maps/' subject '_Task_vs_Rest_Thresholded_SNRExcluded_OddSessions_Overlay_Map_' num2str(threshold) '_Percent'],template)
            
                    clear template
                    
                end
            
            elseif SeparateThresholds == 1 && SNRExclude == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_SNRExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && ReliableExclude == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_ReliabilityExcluded_Overlay_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && GroupSNRMasks == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_GroupSNRMask_Overlay_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && ReliabilityMasks == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_ReliabilityMask_Overlay_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && IndividSNRMasks == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_IndividSNRMask_Overlay_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1 && IndividtSNRMasks == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_IndividtSNRMask_Overlay_Map'],template)
            
                clear template
            
            elseif SeparateThresholds == 1
            
                template.data = combinedoutput;
            
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Thresholded_Overlay_Map'],template)
            
                clear template
            
            elseif MSCreference == 1 && GroupSNRMasks == 1
            
                overlaymap = [cifti_rest_thresh_dat'; cifti_task_thresh_dat'];
                template.data = mean(overlaymap,1)';
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_GroupSNRMask_Overlay_MSCRef_Map'],template)
                clear template
            
            elseif MSCreference == 1
            
                overlaymap = [cifti_rest_thresh_dat'; cifti_task_thresh_dat'];
                template.data = mean(overlaymap,1)';
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Overlay_MSCRef_Map'],template)
                clear template
            
            elseif GroupSNRMasks == 1
            
                overlaymap = [cifti_rest_thresh_dat'; cifti_task_thresh_dat'];
                template.data = mean(overlaymap,1)';
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_GroupSNRMask_Overlay_Map'],template)
                clear template
            
            elseif ReliabilityMasks == 1
            
                overlaymap = [cifti_rest_thresh_dat'; cifti_task_thresh_dat'];
                template.data = mean(overlaymap,1)';
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_ReliabilityMask_Overlay_Map'],template)
                clear template
            
            else
            
                overlaymap = [cifti_rest_thresh_dat'; cifti_task_thresh_dat'];
                template.data = mean(overlaymap,1)';
        
                ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Within_Subjects/' subject '_Task_vs_Rest_Overlay_Map'],template)
                clear template
            
            end
        
        else
        
            AllSubsTask = [AllSubsTask; cifti_task_thresh_dat'];
            AllSubsRest = [AllSubsRest; cifti_rest_thresh_dat'];
        
        end
        end             %% End of task loop for multiple files per subject

    if AcrossSubjects == 1
    
        template = ft_read_cifti_mod('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Rest_Data/Full_Rest/Full_Data/MSC01/MSC01_REST_AllSessions_vs_120_allsubs_corr_cortex_corr.dtseries.nii');
        template.data = [];
    
        AllSubsTask = mean(AllSubsTask,1);
        AllSubsRest = mean(AllSubsRest,1);
        overlaymap = [AllSubsTask; AllSubsRest];
    
        if AcrossSubjects == 1 && MSCreference == 1 && SNRMasks == 1
        
            template.data = mean(AllSubsTask,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Task_SNRMask_Overlay_MSCRef_Map'],template)
    
            template.data = mean(AllSubsRest,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Rest_SNRMask_Overlay_MSCRef_Map'],template)
    
            template.data = mean(overlaymap,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Task_vs_Rest_SNRMask_Overlay_MSCRef_Map'],template)
    
            clear template
    
        elseif AcrossSubjects == 1 && MSCreference == 1
        
            template.data = mean(AllSubsTask,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Task_Overlay_MSCRef_Map'],template)
    
            template.data = mean(AllSubsRest,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Rest_Overlay_MSCRef_Map'],template)
    
            template.data = mean(overlaymap,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Task_vs_Rest_Overlay_MSCRef_Map'],template)
    
            clear template
        
        elseif SNRMasks == 1
        
            template.data = mean(AllSubsTask,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Task_SNRMask_Overlay_Map'],template)
    
            template.data = mean(AllSubsRest,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Rest_SNRMask_Overlay_Map'],template)
    
            template.data = mean(overlaymap,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Task_vs_Rest_SNRMask_Overlay_Map'],template)
    
            clear template
        
        else
    
            template.data = mean(AllSubsTask,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Task_Overlay_Map'],template)
    
            template.data = mean(AllSubsRest,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Rest_Overlay_Map'],template)
    
            template.data = mean(overlaymap,1)';
    
            ft_write_cifti_mod(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Overlayed_Maps/Between_Subjects/AllSubjects_Combined_Task_vs_Rest_Overlay_Map'],template)
    
            clear template
            
        end
            
        end
    end
end

    
    
       