The included files and scripts were used in: Kraus et al., 2021.
Network variants are similar between task and rest states. NIMG. This code
will allow users to match sampling of rest data between tasks (if desired), create 
network variants, find the spatial overlap of variants, assign them to network templates, 
and other subsidiary analyses included in the paper.

The first scripts you will run are in the create variant map folder.

The scripts for creating a variant map use a structure to match data, so if this is desired
one can be created for your data using CreateTmaskStruct_MSC.m. The next step is to create
the variant maps. The following scripts can be used to create variant maps (and other outputs)
with certain data matching constraints:

-Make_Rest_Dconns_MSC.m
-Make_Task_Dconns_MSC.m
-Make_Rest_Dconns_Times_MSC.m
-Make_Task_Dconns_Times_MSC.m

More information on the specific functions of these scripts is available in their documentation.
Once a variant map has been generated, the next step is to threshold it using the 
threshold_variant_maps_wrapper.m function (and threshold_variant_maps.m). Once the threshold is
applied, the code in the analyze variant scripts folder can be used to perform any analysis 
reported in the paper. For convenience, the scripts are listed in the order they should be run 
to recreate all of the figures in the paper from after threshold_variant_maps_wrapper.m is run.

-Table S3
1.CountVertices.m

-Figures 2, S1, S2, S3.
1.TaskReliabilityAnalysis.m or RestReliabilityAnalysis.m

-Figures 3, S4, S9.
1.Make_State_Overlap_Maps.m
2.BreakUpOverlapMaps.m (S9 only)
3.CalcTaskActivationsByState.m (S9 only)

-Figures 4, S5, S6, S7.
1.DiceCorrelations_CombinedTasks.m

-Figure S8
1.Rotate_Parcels_Variant_Magnitude.m

-Figure 5, S10, S11, S13
1.AssignNetworksByVertex.m
2.ComputeDiceCorrelationsNetworkVerticesConsensus.m

-Figure 6
1.PlotSpatialOverlapBins.m

-Figure 7, Figure S12
1.DiceCorrelations_IndividTasks.m

Additionally, supporting .mat or .dtseries.nii (CIFTI) files are included 
in with this code. Please cite our paper when using any of these scripts  
(Kraus et al., 2021. Network variants are similar between task and rest states.
NIMG). For any questions, problems, or bugs found, please email 
btkraus@u.northwestern.edu. Thank you!