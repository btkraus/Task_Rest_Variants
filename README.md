# Task_Rest_Variants
DEPENDENCIES: 
MSC Data: https://openneuro.org/datasets/ds000224/versions/1.0.0
cifti-matlab-master: https://github.com/Washington-University/cifti-matlab
32k_ConteAtlas_v2: https://github.com/gablab/brain_plots/tree/master/32k_ConteAtlas_v2
Connectome workbench: https://www.humanconnectome.org/software/get-connectome-workbench

1. Create SNR Maps?
2.    Create variant maps using MakeTaskDconns.m and MakeTaskDconnsRestLocal.m for task and rest data, respectively. 
    a.    This will need to be ran in a high processing cluster/high power computer 
    b. The following are the conditions for the analysis presented in the paper:
        iv.    CortexOnly = 1;
        v.    SplitHalf = 1;
        vi.    MatchData = 1;
        vii.    RandSample = 1;
        viii.    MatchAcrossSubs = 1;
        ix.    MatchAcrossTasks = 0;
        x.    ConcatenateTasks = 1;
        xi.    MakeDconn = 0; (Dconns are large connectivity matrices files, only make this value 1 if you want them saved on your local computer)
        xii.    MakeVariantMap = 1; (This is the main output of this script and the input for the next script)
        xiii.    ConcatenateSplitHalf = 0;
        xiv.    SaveTimeseries = 1;
    c.   Change the paths in the script to directories in your local computer before running the script 
    d.    The output of this script will be 18 files, (an even and an odd file for each subject) with names like this: MSC01_allTaskCat_matcheddata_Final_REST_EvenSessions_Rand_vs_120_allsubs_corr_cortex_corr.dtseries.nii
3.    Make .txt file containing the path for the variant map, subject number, and whether it’s even or odd (rest or task?) 
    a.    Example:
    
    /Users/dianaperez/Desktop/Brian_Analysis/Run3/VariantFiles/Rest/Even/MSC01_matcheddata_REST_Variant_Size_SNRExclude_Even_2.5.dtseries.nii MSC01 Even    /Users/dianaperez/Desktop/Brian_Analysis/Run3/VariantFiles/Rest/Odd/MSC01_matcheddata_REST_Variant_Size_SNRExclude_Odd_2.5.dtseries.nii MSC01 Odd
    ...
    This for all subjects (except MSC08), separately for rest and task data 


4.    Create Thresholded Variant Maps: Use the script VariantSizeExclusion.m to create variant maps containing only the variants with the lowest 2.5% correlation coefficients based on the variant maps created in step 2 
    a.    This script uses the the .txt files created in step 3 containing paths to the variant maps and  the SNR maps created in step 1
    b.    The following are the conditions for this script:
        i.     CreateVariantFiles = 1;
        ii.    thresholds = [2.5];  
        iii.    ExcludeVariantSize = 0;
        iv.    SNRexclusion = 1;  
        v.    ReliabilityExclusion = 0; 
        xvi.    ExclusionCriteria = 15; 
        vii.    FullMaps = 0;    
        viii.    MatchedMaps = 1; 
        ix.    ConsecMaps = 0;  
        x.    IndividTasks = 0;  
        xi.    MatchedCatTasks = 1;
        xii.    SplitHalf = 1;  
        xiii.    AbsoluteThresholds = 0;
    c.    The output files will be in the directory that you specified in lines 178-179 (lines will change)
5.    Make .txt files containing the path for the thresholded variant files, subject number, and whether it’s even or odd as in step 3
6.    Make Size Excluded Variant Files: Use the same script (VariantSizeExclusion.m) to exclude variants that do not have at least 15 vertices
    a.    This script will use the .txt files made in step 5  the .txt files containing paths to the thresholded variant maps and the SNR maps 
    b.    Change CreateVariantFiles to 0 ad ExcludeVariantSize to 1 so that it looks like this:
        i.    CreateVariantFiles = 0;
        ii.    thresholds = [2.5];  
        iii.    ExcludeVariantSize = 1;
        iv.    SNRexclusion = 1;  
        v.    ReliabilityExclusion = 0; 
        vi.    ExclusionCriteria = 15; 
        vii.    FullMaps = 0;    
        viii.    MatchedMaps = 1; 
        ix.    ConsecMaps = 0;  
        x.    IndividTasks = 0;  
        xi.    MatchedCatTasks = 1;
        xii.    SplitHalf = 1;  
        xiii.    AbsoluteThresholds = 0;
    c.    The script will output the files in the same directory as the thresholded variant files
7.    Make .txt files with the paths to the size excluded variant files 
8.    Use the script QuantifyVariantOverlap.m to calculate the overlap between rest and task variants per subject and create a plot
    a.    This script needs 8 .txt files
        1.    Variant maps rest even
        2.    Variant maps task even
        3.    Variant maps rest odd
        4.    Variant maps task odd
        5.    Size excluded variant files rest even
        6.    Size excluded variant files rest odd
        7.    Size excluded variant files task even
        8.    Size excluded variant files task odd
    b.    The plot will be save in the output directory specified in line 32 (line might change)


