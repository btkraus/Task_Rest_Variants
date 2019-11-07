
%%% Calculate percent overlap for each subject's variants across states

clear all

SplitHalf = 0;
SizeExclude = 1;
threshold = 2.5;
minsize = 50;
usemasks = 0;

VarMatchFinal = [];

load('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Variant_Network_Matching/Matched_Variant_Struct')

if SizeExclude == 1 && SplitHalf == 1
    
	[task_files_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_even_matched_variants_SizeOnly.txt','%s%s%s');

	[rest_files_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_even_matched_variants_SizeOnly.txt','%s%s%s');
    
	[task_files_odd, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_odd_matched_variants_SizeOnly.txt','%s%s%s');

	[rest_files_odd, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_odd_matched_variants_SizeOnly.txt','%s%s%s');
    
	[task_masks_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_even_matched_variants_SizeOnly_sizeexclude.txt','%s%s%s');

	[rest_masks_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_even_matched_variants_SizeOnly_sizeexclude.txt','%s%s%s');
    
	[task_masks_odd, sub1, t1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_odd_matched_variants_SizeOnly_sizeexclude.txt','%s%s%s');

	[rest_masks_odd, sub2, t2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_odd_matched_variants_SizeOnly_sizeexclude.txt','%s%s%s');
    
elseif SizeExclude == 1
    
	[task_files, subjects1, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(threshold) '_matched_variants_SNRexclude.txt'],'%s%s%s');

	[rest_files, subjects2, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(threshold) '_matched_variants_SNRexclude.txt'],'%s%s%s');
    
	[task_masks, sub1, t1] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(threshold) '_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');

	[rest_masks, sub2, t2] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(threshold) '_matched_variants_SNRexclude_sizeexclude.txt'],'%s%s%s');

elseif SplitHalf == 1
    
 	[task_files_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_even_matched_variants_SNRexclude.txt','%s%s%s');

   	[rest_files_even, ~, ~] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_even_matched_variants_SNRexclude.txt','%s%s%s');
    
  	[task_files_odd, subjects1, tasks1] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_splithalf_odd_matched_variants_SNRexclude.txt','%s%s%s');

  	[rest_files_odd, subjects2, tasks2] = textread('/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_splithalf_odd_matched_variants_SNRexclude.txt','%s%s%s');
    
else
    
	[task_files, subjects1, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_alltask_varmaps_' num2str(threshold) '_variants_SNRexclude_matched.txt'],'%s%s%s');

	[rest_files, subjects2, ~] = textread(['/Users/briankraus/Desktop/Correct_Variant_Maps/MSC_Data/Text_Lists/MSC_rest_varmaps_' num2str(threshold) '_variants_SNRexclude_matched.txt'],'%s%s%s');
    
end

if SplitHalf == 1
    
	nfiles = length(rest_files_even);
    
else
    
 	nfiles = length(rest_files);
    
end

for x = 1:nfiles
    
    subject = subjects2{x};
    
  	if SplitHalf == 1
        
    	cifti_rest_even = ft_read_cifti_mod(rest_files_even{x});
      	cifti_task_even = ft_read_cifti_mod(task_files_even{x});
       	cifti_rest_odd = ft_read_cifti_mod(rest_files_odd{x});
       	cifti_task_odd = ft_read_cifti_mod(task_files_odd{x});
        
    else
        
      	cifti_rest = ft_read_cifti_mod(rest_files{x});
      	cifti_task = ft_read_cifti_mod(task_files{x});
        
    end

   	SNRmask = ft_read_cifti_mod(['/Users/briankraus/Desktop/MSC_SNR_Maps/CIFTI_Files/MSC_Template/Individual_Session_Map/' subject '/' subject '__SNRMap__AllDataConcatenated.dscalar.nii']); 
    
    SNRmask.data = SNRmask.data(1:59412,:);
	LowSNR = find(SNRmask.data < 750);
    
    if SplitHalf == 1
            
        cifti_rest_even.data(LowSNR,:) = 0;
        cifti_task_even.data(LowSNR,:) = 0;
        cifti_rest_odd.data(LowSNR,:) = 0;
        cifti_task_odd.data(LowSNR,:) = 0;
        
    else

        cifti_rest.data(LowSNR,:) = 0;
        cifti_task.data(LowSNR,:) = 0;
        
    end
    
  	if SizeExclude == 1         %% Apply exclusion masks for maps that exclude size
        
        if usemasks == 1
        
            if SplitHalf == 1
            
                cifti_rest_mask_even = ft_read_cifti_mod(rest_masks_even{x});
                cifti_task_mask_even = ft_read_cifti_mod(task_masks_even{x});
                cifti_rest_mask_odd = ft_read_cifti_mod(rest_masks_odd{x});
                cifti_task_mask_odd = ft_read_cifti_mod(task_masks_odd{x});
        
                for d = 1:length(cifti_rest_mask_even.data)
            
                    if cifti_rest_mask_even.data(d) == 0
                
                        cifti_rest_even.data(d) = 0;
                
                    end
                
                    if cifti_rest_mask_odd.data(d) == 0
                    
                        cifti_rest_odd.data(d) = 0;
                    
                    end  
                end
        
                for e = 1:length(cifti_task_mask_even.data)
            
                    if cifti_task_mask_even.data(e) == 0
                
                        cifti_task_even.data(e) = 0;
                
                    end
                
                    if cifti_task_mask_odd.data(e) == 0
                    
                        cifti_task_odd.data(e) = 0;
                    
                    end
                end
            
            else
        
                cifti_rest_mask = ft_read_cifti_mod(rest_masks{x});
                cifti_task_mask = ft_read_cifti_mod(task_masks{x});
        
                for d = 1:length(cifti_rest_mask.data)
            
                    if cifti_rest_mask.data(d) == 0
                
                        cifti_rest.data(d) = 0;
                
                    end
                end
        
                for e = 1:length(cifti_task_mask.data)
            
                    if cifti_task_mask.data(e) == 0
                
                        cifti_task.data(e) = 0;
                
                    end
                end
            end
            
        else
            
            if SplitHalf == 1
                
                allvars_rest_even = unique(cifti_rest_even.data);
                allvars_rest_even(1) = [];
                allvars_rest_odd = unique(cifti_rest_odd.data);
                allvars_rest_odd(1) = [];
                allvars_task_even = unique(cifti_task_even.data);
                allvars_task_even(1) = [];
                allvars_task_odd = unique(cifti_task_odd.data);
                allvars_task_odd(1) = [];
                
                removevars_rest_even = [];
                removevars_rest_odd = [];
                removevars_task_even = [];
                removevars_task_odd = [];
                
                for g = 1:length(allvars_rest_even)
                    
                    if length(find(cifti_rest_even.data == allvars_rest_even(g))) < minsize
                        
                        removevars_rest_even = [removevars_rest_even allvars_rest_even(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_rest_even.data)
                    
                    if ismember(cifti_rest_even.data(h),removevars_rest_even)
                        
                        cifti_rest_even.data(h) = 0;
                        
                    end
                end
                
                for g = 1:length(allvars_rest_odd)
                    
                    if length(find(cifti_rest_odd.data == allvars_rest_odd(g))) < minsize
                        
                        removevars_rest_odd = [removevars_rest_odd allvars_rest_odd(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_rest_odd.data)
                    
                    if ismember(cifti_rest_odd.data(h),removevars_rest_odd)
                        
                        cifti_rest_odd.data(h) = 0;
                        
                    end
                end
                
                for g = 1:length(allvars_task_even)
                    
                    if length(find(cifti_task_even.data == allvars_task_even(g))) < minsize
                        
                        removevars_task_even = [removevars_task_even allvars_task_even(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_task_even.data)
                    
                    if ismember(cifti_task_even.data(h),removevars_task_even)
                        
                        cifti_task_even.data(h) = 0;
                        
                    end
                end
                
                for g = 1:length(allvars_task_odd)
                    
                    if length(find(cifti_task_odd.data == allvars_task_odd(g))) < minsize
                        
                        removevars_task_odd = [removevars_task_odd allvars_task_odd(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_task_odd.data)
                    
                    if ismember(cifti_task_odd.data(h),removevars_task_odd)
                        
                        cifti_task_odd.data(h) = 0;
                        
                    end
                end
                
            else
                
            	allvars_rest = unique(cifti_rest.data);
                allvars_rest(1) = [];
                allvars_task = unique(cifti_task.data);
                allvars_task(1) = [];
                
                removevars_rest = [];
                removevars_task = [];
                
                for g = 1:length(allvars_rest)
                    
                    if length(find(cifti_rest.data == allvars_rest(g))) < minsize
                        
                        removevars_rest = [removevars_rest allvars_rest(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_rest.data)
                    
                    if ismember(cifti_rest.data(h),removevars_rest)
                        
                        cifti_rest.data(h) = 0;
                        
                    end
                end
                
                for g = 1:length(allvars_task)
                    
                    if length(find(cifti_task.data == allvars_task(g))) < minsize
                        
                        removevars_task = [removevars_task allvars_task(g)];
                        
                    end
                end
                
                for h = 1:length(cifti_task.data)
                    
                    if ismember(cifti_task.data(h),removevars_task)
                        
                        cifti_task.data(h) = 0;
                        
                    end
                end
                
            end
        end
    end

	if SplitHalf == 1       %% Get number of variants
        
     	vars_rest_even = unique(cifti_rest_even.data);
    	vars_rest_even(1) = [];
      	vars_rest_odd = unique(cifti_rest_odd.data);
        vars_rest_odd(1) = [];
    	vars_task_even = unique(cifti_task_even.data);
    	vars_task_even(1) = [];
     	vars_task_odd = unique(cifti_task_odd.data);
     	vars_task_odd(1) = [];
        
    else
        
      	vars_rest = unique(cifti_rest.data);
       	vars_rest(1) = [];
      	vars_task = unique(cifti_task.data);
       	vars_task(1) = [];
        
    end

    VarMatchSub = 0;
    
    for a = 1:length(vars_rest)
        
        for b = 1:size(Variant_Match_Struct,2)
            
            if strcmp(Variant_Match_Struct(b).State, 'Rest') && Variant_Match_Struct(b).Matched == 1
                
                if vars_rest(a) == Variant_Match_Struct(b).Rest_Var_Num && strcmp(subject, Variant_Match_Struct(b).Subject)
                
                    RestNet = Variant_Match_Struct(b).Network;
                    TaskVars = Variant_Match_Struct(b).Task_Var_Num;
                    
                    for c = 1:length(TaskVars)
                        
                        for d = 1:size(Variant_Match_Struct,2)
                        
                            if Variant_Match_Struct(d).Matched == 1 && strcmp(Variant_Match_Struct(d).State, 'Task') && strcmp(subject, Variant_Match_Struct(d).Subject) && TaskVars(c) == Variant_Match_Struct(d).Task_Var_Num
                                
                                if strcmp(Variant_Match_Struct(d).Network, RestNet)
                                    
                                    VarMatchSub = VarMatchSub + 1;
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    VarMatchFinal = [VarMatchFinal VarMatchSub/round((length(vars_rest)+length(vars_task))/2)];
    
end
        
                


    
    
    
    