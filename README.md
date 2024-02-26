# TRAP-Rabies-Analysis
Code and data related to Figures 2 and S2 in Ryan et al., 2024

D1-cre, A2a-cre, and FosTRAP-creER mice were first injected with a "helper" virus into the striatum (green fluorophore - "starter" cells). A modified, G-deleted form of the rabies virus was then injected into the same target in the striatum to label presynaptic (red fluorophore) neurons across the brain. 

Rabies_Analysis.m performs analysis of data collected using NeuroInfo software to detect pre-synaptic, rabies-labeled neurons, map them to the allen brain atlas and normalize to the number of co-infected (green/red) cells in the striatum for each mouse. Saved into a struct (Rabies.mat).

Rabies_Summary.mat is a struct that contains summary normalized data for each animal ("data") and for each cohort ("group_stats") used to generate figures 2 and S2.

Raw and processed data can be found here:
