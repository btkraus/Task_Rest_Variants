parpool('local', 20)     %% Name of cluster profile for batch job


clear all


SNR = 0;  %% Toggles whether to calculate SNR with signal measure only
tSNR = 1;  %% Toggles whether to calculate tSNR incorporating noise
tmasks = 0;  %% Toggles whether to apply tmask to each session for SNR maps
ConcatenateSubs = 0;  %% Toggles whether to concatenate subs for a group mask
MSCtemplate = 1;  %% Toggles whether to use MSC template or generic template
SessionTemplate = 1;  %% Toggles whether to use a separate MSC template for each session

outdir = '/projects/b1081/Brian_MSC/dconn_task_files/SNR_Maps/';
%subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
subs = {'MSC01'};


disp(sprintf('Job Submitted: %s', datestr(now)));


disp(sprintf('Job Started: %s', datestr(now)));

if SessionTemplate == 1
    
    outdir= strcat(outdir, 'Session_Maps/');
    
end


catData = [];

for i=1:numel(subs)
    
    
    % Create template path for resting data
    MSCciftidir = ['/projects/b1081/MSC/TaskFC/FCProc_' subs{i} '_mem_pass2/cifti_timeseries_normalwall_native_freesurf'];
    
    
    % Get correct rest vcids
    restdir = dir(['/projects/b1081/MSC/MSCdata_v1/' subs{i} '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall_native_freesurf']);
    vcidlist = restdir(~cellfun(@isempty,strfind({restdir.name},'dtseries')));
    vcidlist = vcidlist(~cellfun(@isempty,strfind({vcidlist.name},'vc')));
            
    disp(sprintf('%i resting sessions found for subject %s: %s', size(vcidlist,1), subs{i}, datestr(now)));
    
    if tmasks == 1
    
        % Load rest tmasks
        load(['/projects/b1081/MSC/MSCdata_v1/' subs{i} '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/QC.mat'])
        
    end
    
    % Load and concatentate data
    for k=1:length(vcidlist)
        
        vcidname = strsplit(vcidlist(k).name, '_');
    	vcid = vcidname{1};
        
        if strcmp(subs{i}, 'MSC02') && k == 9
            
            vcid = strcat(vcidname{1}, '_A');
            
        end
            
    	disp(sprintf('Loading resting state timeseries 4dfp file for session %s for subject %s task %s: %s', vcid, subs{i}, datestr(now)));
        
        data = ['/projects/b1081/MSC/MSCdata_v1/' subs{i} '/Functionals/' vcid '/bold1/' vcid '_b1_faln_dbnd_xr3d_uwrp_atl.4dfp.img'];
        [inputdata] = read_4dfpimg(data);
        [v, f, etype] = fcimage_attributes(data); %v= voxel size, f=number of frames or fourth dimension, etype = endian type (way data is stored)
        
        if tmasks == 1
    
            resttmask = QC(k).tmask;
                    
            disp(sprintf('tmask for rest file has %i good sample points, %s', sum(resttmask), datestr(now)));
            
            inputdata = inputdata(:,logical(resttmask));
            
        end
        
        if MSCtemplate == 1 && SessionTemplate == 1
            
            disp('Calculating SNR')
            
        	if SNR == 1
            
            	MeanSNR = mean(inputdata,2);
    
          	elseif tSNR == 1

            	signal = mean(inputdata,2);
        
                noise = std(inputdata,0,2);
    
                MeanSNR = signal./noise;
                   
            end
                
          	if SNR == 1 && tmasks == 1
        
             	outname = ([subs{i} '_' 'SNRMap_tmasks_REST_MSCTemplate_' vcid '.4dfp.img']);
        
          	elseif SNR == 1
        
               	outname = ([subs{i} '_' 'SNRMap_REST_MSCTemplate_' vcid '.4dfp.img']);
        
           	elseif SNR == 1 && tmasks == 1
        
              	outname = ([subs{i} '_' 'tSNRMap_tmasks_REST_MSCTemplate_' vcid '.4dfp.img']);
        
          	elseif tSNR == 1
        
               	outname = ([subs{i} '_' 'tSNRMap_REST_MSCTemplate_' vcid '.4dfp.img']);
        
            end
                
           	out_data = MeanSNR;
           	fout = [outdir outname];
    
           	disp('Writing .4dfp file')
    
           	write_4dfpimg(out_data,fout,etype);
          	write_4dfpifh(fout,1,etype); %note that the 1 denotes this is only 1 volume large; etype should be the same as when the data was loaded
    
          	disp('Mapping volume to surface')
                
          	map_vol_to_surface_MSCspecific(fout,subs{i},vcid)
            
            clear out_data
    
            niftiout = strrep(fout, '.4dfp.img', '.nii');
    
            disp('Creating NIFTI from .4dfp')
    
            system(['nifti_4dfp -n ' fout ' ' niftiout]);
                
                
        else
            
            disp(sprintf('Adding %i sample points to catData, %s', size(inputdata,2), datestr(now)));
        
        
            catData = [catData inputdata];
        
        end
        
    end
%end        %% End subject loop for concatenated data

    if SessionTemplate == 0
        
        disp('Calculating SNR')

        if SNR == 1
            
            MeanSNR = mean(catData,2);
    
        elseif tSNR == 1

            signal = mean(catData,2);
        
            noise = std(catData,0,2);
    
            MeanSNR = signal./noise;
                   
        end
    
        if SNR == 1 && tmasks == 1 && ConcatenateSubs == 1
        
            outname = ('AllSubs_SNRMap_tmasks_REST_AllSessions.4dfp.img');
    
        elseif SNR == 1 && ConcatenateSubs == 1
        
            outname = ('AllSubs_SNRMap_REST_AllSessions.4dfp.img');
        
        elseif tSNR == 1 && tmasks == 1 && ConcatenateSubs == 1
        
            outname = ('AllSubs_tSNRMap_tmasks_REST_AllSessions.4dfp.img');
        
        elseif tSNR == 1 && ConcatenateSubs == 1
        
            outname = ('AllSubs_tSNRMap_REST_AllSessions.4dfp.img');
        
        elseif SNR == 1 && tmasks == 1 && MSCtemplate == 1
            
            outname = ([subs{i} '_' 'SNRMap_tmasks_REST_MSCTemplate_AllSessions.4dfp.img']);
            
        elseif SNR == 1 && tmasks == 1
        
            outname = ([subs{i} '_' 'SNRMap_tmasks_REST_AllSessions.4dfp.img']);
        
        elseif SNR == 1 && MSCtemplate == 1
        
            outname = ([subs{i} '_' 'SNRMap_REST_MSCTemplate_AllSessions.4dfp.img']);
            
        elseif SNR == 1
        
            outname = ([subs{i} '_' 'SNRMap_REST_AllSessions.4dfp.img']);
            
        elseif tSNR == 1 && tmasks == 1 && MSCtemplate == 1
        
            outname = ([subs{i} '_' 'tSNRMap_tmasks_REST_MSCTemplate_AllSessions.4dfp.img']);
        
        elseif tSNR == 1 && tmasks == 1
        
            outname = ([subs{i} '_' 'tSNRMap_tmasks_REST_AllSessions.4dfp.img']);
            
        elseif tSNR == 1 && MSCtemplate == 1
        
            outname = ([subs{i} '_' 'tSNRMap_REST_MSCTemplate_AllSessions.4dfp.img']);    
        
        elseif tSNR == 1
        
            outname = ([subs{i} '_' 'tSNRMap_REST_AllSessions.4dfp.img']);
        
        end
    
    
        out_data = MeanSNR;
        fout = [outdir outname];
    
        disp('Writing .4dfp file')
    
        write_4dfpimg(out_data,fout,etype);
        write_4dfpifh(fout,1,etype); %note that the 1 denotes this is only 1 volume large; etype should be the same as when the data was loaded
    
        disp('Mapping volume to surface')
    
    
        if MSCtemplate == 1     %% Use MSC-specific template
        
            sessionvox = strsplit(vcidlist(1).name, '_');
        
            map_vol_to_surface_MSCspecific(fout,subs{i},char(sessionvox(1)))
        
        
        elseif MSCtemplate == 0     %% Use generic templated
    
            map_vol_to_surface(fout,'both','ribbon-constrained','711-2B'); %%% this function maps volume data to the the surface using a group estimate (it's not as precise as the individualized method we usually use for the MSC, but is good for a quick look at the data) - it can be found in scripts/WorkbenchScripts/map_vol_to_surface.m;
            %both = both hemispheres or just one
            %ribbon-constrained = how the interpolation is done (within the gray matter ribbon, or other options, including no interpolation/averaging
            %711-2B = group template space that the data is in; WashU data is usually in 711-2B, but other datasets are often in MN
        
        end
    
        clear out_data
    
        niftiout = strrep(fout, '.4dfp.img', '.nii');
    
        disp('Creating NIFTI from .4dfp')
    
        system(['nifti_4dfp -n ' fout ' ' niftiout]);
    
    end


    catData = [];
    
end        %% End subject loop for individual data
    
    
    
    