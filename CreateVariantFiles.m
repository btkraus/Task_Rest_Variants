%% CreateVariantFiles.m
%This script makes the variants from spatial correlation maps, excluding
%regions with low signal. 
%Dependencies: SNR maps, wb_command, 32k_ConteAtlas_v2_distribute,cifti-matlab-master, some fieldtrip
%functions (ft_read_cifti_mod, ft_write_cifti_mod).
%Written by Brian Kraus. Edited by Diana Perez.

clear all

%% Variables
ExcludeBySize = 1;
SplitHalf = 1;
ExclusionCriteria = 15;  %minimum number of vertices required for variant to not be excluded
threshold = 2.5;  %% Thresholds used to calculate variants (lowest % or correlation values)
SNRexclusion = 1;  %% Toggles whether to exclude variants based on SNR, 1 = exclude, 0 = don't exclude
%% Paths
workbenchdir = '/Applications/workbench/bin_macosx64/';
leftsurf = '/Users/dianaperez/Box/Dependencies/32k_ConteAtlas_v2_distribute/Conte69.L.midthickness.32k_fs_LR.surf.gii';
rightsurf = '/Users/dianaperez/Box/Dependencies/32k_ConteAtlas_v2_distribute/Conte69.R.midthickness.32k_fs_LR.surf.gii';
dirpath = '/Users/dianaperez/Documents/GitHub/Task_Rest_Variants/';
SNRpath = '/Users/dianaperez/Box/Latest_Analysis_Replication/SNR_maps/';
outfilepath = '/Users/dianaperez/Box/Latest_Analysis_Replication/';
resttxtname = 'MSC_rest_spCorrMaps.txt';
tasktxtname = 'MSC_task_spCorrMaps.txt';
addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master/'));

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
        SNRmap = ft_read_cifti_mod([SNRpath subject '_SNRMap_REST_MSCTemplate_AllSessions.4dfp.img_LR_surf_subcort_333_32k_fsLR.dscalar.nii']);

        SNRmap.data = SNRmap.data(1:59412,:);
        SNRexclude = find(SNRmap.data < 750);

        cifti_rest.data(SNRexclude,1) = NaN;
        cifti_task.data(SNRexclude,1) = NaN;
    end

    %% Makes variant mask 
    
    %finds vertices below the threshold percentile 
    cifti_task_threshold = find(cifti_task.data < prctile(cifti_task.data,threshold));
    cifti_rest_threshold = find(cifti_rest.data < prctile(cifti_rest.data,threshold));

    % sets all values to zero except for those that are below the threshold
    cifti_rest_thresh_dat = zeros(size(cifti_rest.data));
    cifti_rest_thresh_dat(cifti_rest_threshold,1) = 1;

    cifti_task_thresh_dat = zeros(size(cifti_task.data));
    cifti_task_thresh_dat(cifti_task_threshold,1) = 1;

    cifti_rest_final_dat = zeros(size(cifti_rest.data));
    cifti_task_final_dat = zeros(size(cifti_task.data));

    % sets value = 1 if vertex is below the threshold
    for w = 1:length(cifti_rest.data)
        if cifti_rest_thresh_dat(w) == 1
            cifti_rest_final_dat(w) = 1;
        end
        if cifti_task_thresh_dat(w) == 1
            cifti_task_final_dat(w) = 1;
        end
    end
    
    % creates strings to be included in file name for outputs
    if SNRexclusion == 1
        if SplitHalf == 1
            strrest = ['SNRExclude_' restsplithalf{x}];
            strtask = ['SNRExclude_' tasksplithalf{x}];
        else
            strrest = ['SNRExclude_'];
            strtask = ['SNRExclude_'];
        end
    else
        if SplitHalf == 1
            strrest = [restsplithalf{x}];
            strtask = [tasksplithalf{x}];
        else
            strrest = ['_'];
            strtask = ['_'];
        end
    end
    
    % creates subfolders for thresholded variant maps
    subFolderRest = [outfilepath 'Thresholded_Variant_Maps/Rest'];
    subFolderTask = [outfilepath 'Thresholded_Variant_Maps/Task'];
    if ~isfolder(subFolderRest)
        mkdir(subFolderRest)
    end
    if ~isfolder(subFolderTask)
        mkdir(subFolderTask)
    end
    
    %This creates the output file names for rest and task 
    outfilerest = [subFolderRest '/' subject '_Rest_ThresholdedVariantMap_' strrest '_' num2str(threshold) '.dtseries.nii'];
    outfiletask = [subFolderTask '/' subject '_Task_ThresholdedVariantMap_' strtask '_' num2str(threshold) '.dtseries.nii'];

    %this creates and writes the file in cifti format
    cifti_rest.data = cifti_rest_final_dat;
    cifti_task.data = cifti_task_final_dat;
    ft_write_cifti_mod(outfilerest, cifti_rest)
    ft_write_cifti_mod(outfiletask, cifti_task)
    
    % creates subfolder for size excluded variant maps
    subFolderRest = [outfilepath 'SizeExcluded_Variant_Maps/Rest'];
    subFolderTask = [outfilepath 'SizeExcluded_Variant_Maps/Task'];
    if ~isfolder(subFolderRest)
        mkdir(subFolderRest)
    end
    if ~isfolder(subFolderTask)
        mkdir(subFolderTask)
    end
    
    %sets names for output files
    outfilewbtask = [subFolderTask '/' subject '_Task_Variant_SizeExcluded_' strtask '_' num2str(threshold) '.dtseries.nii'];
    outfilewbrest = [subFolderRest '/' subject '_Rest_Variant_SizeExcluded_' strrest '_' num2str(threshold) '.dtseries.nii'];

    % This numbers that variants, clustering adjacent vertices that =1
    % together and assigning each cluster a number
    system([workbenchdir 'wb_command -cifti-find-clusters ' outfilerest ' 0 0 0 0 COLUMN ' outfilewbrest ' -left-surface ' leftsurf ' -right-surface ' rightsurf])
    system([workbenchdir 'wb_command -cifti-find-clusters ' outfiletask ' 0 0 0 0 COLUMN ' outfilewbtask ' -left-surface ' leftsurf ' -right-surface ' rightsurf])

    cifti_rest = ft_read_cifti_mod(outfilewbrest);
    cifti_task = ft_read_cifti_mod(outfilewbtask);

    if ExcludeBySize == 1 
        [cifti_rest.data, cifti_task.data] = ExcludeVariantSize(cifti_rest.data, cifti_task.data, subject, threshold, ExclusionCriteria);
    end
    
    %% writes output file with numbered variants
    ft_write_cifti_mod(outfilewbrest, cifti_rest)
    ft_write_cifti_mod(outfilewbtask, cifti_task)
  
    %% Creates txt file outputs for next script
    if SplitHalf == 1
        if strcmp(restsplithalf(x),'Even')
            resttxtline = [outfilewbrest ' ' subject ' Even'];            
            if ~exist('evenresttxtfile')
                evenresttxtfile = resttxtline;
            else
                evenresttxtfile = [evenresttxtfile; resttxtline];
            end
        elseif strcmp(restsplithalf(x), 'Odd_')
            resttxtline = [outfilewbrest ' ' subject ' Odd_'];
            if ~exist('oddresttxtfile')
                oddresttxtfile = resttxtline;
            else
                oddresttxtfile = [oddresttxtfile; resttxtline];
            end
        end
        if strcmp(tasksplithalf(x),'Even')
            tasktxtline = [outfilewbtask ' ' subject ' Even'];
            if ~exist('eventasktxtfile')
                eventasktxtfile = tasktxtline;
            else
                eventasktxtfile = [eventasktxtfile; tasktxtline];
            end
        elseif strcmp(tasksplithalf(x), 'Odd_')
            tasktxtline = [outfilewbtask ' ' subject ' Odd_'];
            if ~exist('oddtasktxtfile')
                oddtasktxtfile = tasktxtline;
            else
                oddtasktxtfile = [oddtasktxtfile; tasktxtline];
            end
        end
    else
        resttxtline = [outfilewbrest ' ' subject ' ' threshold];
        tasktxtline = [outfilewbtask ' ' subject ' ' threshold];
        if ~exist('restouttxtfile')
            resttxtfile = resttxtline;
        else
            resttxtfile = [resttxtfile; resttxtline];
        end
        if ~exist('taskouttxtfile')
            tasktxtfile = tasktxtline;
        else
            tasktxtfile = [tasktxtfile; tasktxtline];
        end
    end
end
%% writes txt files
if SplitHalf == 1
    writematrix(evenresttxtfile,sprintf('%sMSC_rest_varMaps_Even.txt',outfilepath));
    writematrix(eventasktxtfile,sprintf('%sMSC_task_varMaps_Even.txt',outfilepath));
    writematrix(oddresttxtfile,sprintf('%sMSC_rest_varMaps_Odd.txt',outfilepath));
    writematrix(oddtasktxtfile,sprintf('%sMSC_task_varMaps_Odd.txt',outfilepath));
else
    writematrix(resttxtfile,sprintf('%sMSC_rest_varMaps.txt',outfilepath));
    writematrix(tasktxtfile,sprintf('%sMSC_task_varMaps.txt',outfilepath));
end

function [cifti_rest_data cifti_task_data] = ExcludeVariantSize(cifti_rest_data, cifti_task_data, subject, threshold, exclusion_criteria)

    %This function excludes variants that are less than the minimum number
    %of adjacent vertices specified by exclusion_criteria.
    % WRITE DESCRITION OF INPUTS AND OUTPUTS WITH SOME DETAILS OF WHAT FORMAT OF INPUTS IS SUPPOSED TO BE
    %. Written by Brian Kraus. Edited by Diana Perez.     

        rest_sizes = [];
        task_sizes = [];
        
        % same values as in cifti_rest.data but w/o repetitions, of
        % vertices?
        vars_rest = unique(cifti_rest_data);
        vars_task = unique(cifti_task_data);
        
        
        %counts vertices to determine if variant meet exclusion criteria
        for q = 1:length(vars_rest)
            vertcount = 0;
            for r = 1:length(cifti_rest_data)
                if cifti_rest_data(r) == vars_rest(q)
                    vertcount = vertcount + 1;
                end
            end
            rest_sizes = [rest_sizes; vars_rest(q) vertcount];
        end
        
        for q = 1:length(vars_task)
            vertcount = 0;
            for r = 1:length(cifti_task_data)
                if cifti_task_data(r) == vars_task(q)                    
                    vertcount = vertcount + 1;                    
                end
            end            
            task_sizes = [task_sizes; vars_task(q) vertcount];            
        end
        
        %removing vertices belonging to variants that are not at least 15
        %vertices big
        for i = 2:size(rest_sizes,1)            
            if rest_sizes(i,2) < exclusion_criteria                
                removeverts = find(cifti_rest_data == rest_sizes(i,1));                
                cifti_rest_data(removeverts,1) = 0;                
            else                
                setverts = find(cifti_rest_data == rest_sizes(i,1));                
                cifti_rest_data(setverts,1) = 1;                
            end
        end
        
        for i = 2:size(task_sizes,1)            
            if task_sizes(i,2) < exclusion_criteria                
                removeverts = find(cifti_task_data == task_sizes(i,1));                
                cifti_task_data(removeverts,1) = 0;                
            else                
                setverts = find(cifti_task_data == task_sizes(i,1));                
                cifti_task_data(setverts,1) = 1;                
            end
        end
    end
