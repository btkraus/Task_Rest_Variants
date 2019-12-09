function [cifti_rest.data, cifti_task.data] = ExcludeVariantSize(cifti_rest.data, cifti_task.data, subject, threshold)

    %This function excludes variants that are less than 15 adjacent vertices
    %big. WRITE DESCRITION OF INPUTS AND OUTPUTS WITH SOME DETAILS OF WHAT FORMAT OF INPUTS IS SUPPOSED TO BE
    %. Written by Brian Kraus. Edited by Diana Perez.     
    
    
        %minimum number of vertices required for variant; MAKE THIS ONE OF THE INPUTS
        ExclusionCriteria = 15;
        
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
        %vertices big
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
 end
        

        
