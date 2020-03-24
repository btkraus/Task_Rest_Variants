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
        

        
