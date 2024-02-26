%Rabies analysis code used to generate data in Figures 2, S2 in Ryan et al.,
%2024.

%Outputs from NeuroInfo (excel sheets) are first read into matlab and
%compiled into a struct named "Rabies". The following code normalizes this
%data using several normalization methods. The normalization used in the
%2024 paper is the "Colocalization" normalization which takes the
%rabies-labeled red cells in each brain region and divides it by the total
%number of co-infected (red/green) cells in the striatum.

%Created by Michael B. Ryan 01/01/20


%% Normalization to Green Starter Cells, Colocalized Green/Red Striatal Cells, All Red Cells Brainwide, all Red Cells Brainwide - Red Striatum, and all Red Striatal Cells
groups = fieldnames(Rabies);
for i = 1:length(groups)
    animals = fieldnames(Rabies.(groups{i}));
    for a = 1:length(animals)
        animal = Rabies.(groups{i}).(animals{a});
        
        %Test for overexpression of starter cells outside the striatum
        Str_starter_prop = (animal.data{34,4})/(animal.data{3,4}); %divide all green "starter" cells in striatum but total brain
        Total_Str_starter = (animal.data{34,4}); %all green cells in striatum
        Coloc_Str_Starter = (animal.data{34,6}); %all colocalized (green/red) cells in striatum
        
        Rabies.(groups{i}).(animals{a}).info.Str_Starter_Prop = [Str_starter_prop];
        Rabies.(groups{i}).(animals{a}).info.Total_Str_Starter = [Total_Str_starter];
        Rabies.(groups{i}).(animals{a}).info.Coloc_Str_Starter = [Coloc_Str_Starter];
        
        %If striatal starter expression is over 85% perform analysis,
        %otherwise skip this animal
        if Str_starter_prop > 0.85
            %Normalize to the number of green starter cells
            Str_starter = (animal.data{34,4});
            Starter_Normalized_data = animal.data(:,1:2);
            Starter_Normalized_data{1,3} = 'Proportion of Red Cells';
            %Normalize to colocalized red and green striatal cells
            Coloc_starter = (animal.data{34,6});
            Coloc_Normalized_data = animal.data(:,1:2);
            Coloc_Normalized_data{1,3} = 'Proportion of Red Cells';
            %Normalize to all red brainwide cells
            All_red_cells = (animal.data{3,3})-(animal.data{34,3})+(animal.data{34,5});
            All_Red_Normalized_data = animal.data(:,1:2);
            All_Red_Normalized_data{1,3} = 'Proportion of Red Cells';
            %Normalize to all red brainwide - red striatum cells
            NonStr_red_cells = (animal.data{3,3})-(animal.data{34,3});
            NonStr_Red_Normalized_data = animal.data(:,1:2);
            NonStr_Red_Normalized_data{1,3} = 'Proportion of Red Cells';
            %Normalize to all red brainwide - colocalized striatum cells
            All_red_minus_coloc_cells = (animal.data{3,3})-(animal.data{34,3})+((animal.data{34,5})-(animal.data{34,6}));
            All_Red_Minus_Coloc_Normalized_data = animal.data(:,1:2);
            All_Red_Minus_Coloc_Normalized_data{1,3} = 'Proportion of Red Cells';
            %Normalize to all red striatal cells
            Str_red_cells = (animal.data{34,5});
            Str_Red_Normalized_data = animal.data(:,1:2);
            Str_Red_Normalized_data{1,3} = 'Proportion of Red Cells';
            %Normalize to all red cortical cells
            Ctx_red_cells = (animal.data{5,3});
            Ctx_Red_Normalized_data = animal.data(:,1:2);
            Ctx_Red_Normalized_data{1,3} = 'Proportion of Red Cells';
            
            for b = [2:33,35:length(animal.data)] %use cortical detection settings for whole brain
                Starter_Normalized_data{b,3} = [animal.data{b,3}]'/Str_starter;
                Coloc_Normalized_data{b,3} = [animal.data{b,3}]'/Coloc_starter;
                All_Red_Normalized_data{b,3} = [animal.data{b,3}]'/All_red_cells;
                NonStr_Red_Normalized_data{b,3} = [animal.data{b,3}]'/NonStr_red_cells;
                Str_Red_Normalized_data{b,3} = [animal.data{b,3}]'/Str_red_cells;
                All_Red_Minus_Coloc_Normalized_data{b,3} = [animal.data{b,3}]'/All_red_minus_coloc_cells;
                Ctx_Red_Normalized_data{b,3} = [animal.data{b,3}]'/Ctx_red_cells;           
            end
            
            for b = [34] %use separate detection settings for the striatal cells
                Starter_Normalized_data{b,3} = [animal.data{b,5}]'/Str_starter;
                Coloc_Normalized_data{b,3} = [animal.data{b,5}]'/Coloc_starter;
                All_Red_Normalized_data{b,3} = [animal.data{b,5}]'/All_red_cells;
                NonStr_Red_Normalized_data{b,3} = [animal.data{b,5}]'/NonStr_red_cells;
                Str_Red_Normalized_data{b,3} = [animal.data{b,5}]'/Str_red_cells;
                All_Red_Minus_Coloc_Normalized_data{b,3} = [animal.data{b,5}]'/All_red_minus_coloc_cells;
                Ctx_Red_Normalized_data{b,3} = [animal.data{b,5}]'/Ctx_red_cells;
            end
            
            Rabies.(groups{i}).(animals{a}).starter_normalized_data = Starter_Normalized_data;
            Rabies.(groups{i}).(animals{a}).coloc_normalized_data = Coloc_Normalized_data;
            Rabies.(groups{i}).(animals{a}).all_red_normalized_data = All_Red_Normalized_data;
            Rabies.(groups{i}).(animals{a}).nonstr_red_normalized_data = NonStr_Red_Normalized_data;
            Rabies.(groups{i}).(animals{a}).str_red_normalized_data = Str_Red_Normalized_data;
            Rabies.(groups{i}).(animals{a}).all_red_minus_coloc_normalized_data = All_Red_Minus_Coloc_Normalized_data;
            Rabies.(groups{i}).(animals{a}).ctx_red_normalized_data = Ctx_Red_Normalized_data;
            Rabies.(groups{i}).(animals{a}).info.NonStr_Red_Coloc_Normalized = NonStr_red_cells/Coloc_starter;
        else
        end
    end
end
%% Summary Structs: compile summary data following colocalized normalization for each mouse (quantification used for Figure 2, S2)
groups = fieldnames(Rabies);
brain_regions = {[],[];'Brain Region','Hemisphere';'root','left';'root','right';'Cerebral cortex','left';'Cerebral cortex','right';'Somatomotor areas','left';'Somatomotor areas','right';'Primary motor area','left';'Primary motor area','right';'Secondary motor area','left';'Secondary motor area','right';'Somatosensory areas','left';'Somatosensory areas','right';'Primary somatosensory area','left';'Primary somatosensory area','right';'Primary somatosensory area/ nose','right';'Primary somatosensory area/ barrel field','right';'Primary somatosensory area/ lower limb','right';'Primary somatosensory area/ mouth','right';'Primary somatosensory area/ upper limb','right';'Primary somatosensory area/ trunk','right';'Supplemental somatosensory area','right';'Gustatory areas','right';'Visceral area','right';'Anterior cingulate area','left';'Anterior cingulate area','right';'Anterior cingulate area/ dorsal part','right';'Anterior cingulate area/ ventral part','right';'Orbital area','left';'Orbital area','right';'Agranular insular area','left';'Agranular insular area','right';'Olfactory areas','right';'Striatum','right';'Central amygdalar nucleus','right';'Globus pallidus/ external segment','right';'Pallidum/ ventral region','right';'Thalamus','right';'Ventral group of the dorsal thalamus','right';'Ventral medial nucleus of the thalamus','right';'Thalamus/ polymodal association cortex related','right';'Intralaminar nuclei of the dorsal thalamus','right';'Central medial nucleus of the thalamus','right';'Paracentral nucleus','right';'Parafascicular nucleus','right';'Substantia nigra/ reticular part','right';'Midbrain raphe nuclei','right';'Dorsal nucleus raphe','right';'Primary auditory area','right';'Primary visual area','right';'Prelimbic area','right';'Infralimbic area','right';'Posterior parietal association areas','right';'Piriform area','right';'Hippocampal formation','right';'Lateral amygdalar nucleus','left';'Lateral amygdalar nucleus','right';'Basolateral amygdalar nucleus','left';'Basolateral amygdalar nucleus','right';'Striatum dorsal region','left';'Striatum dorsal region','right';'Striatum ventral region','left';'Striatum ventral region','right';'Nucleus accumbens','left';'Nucleus accumbens','right';'Anterior amygdalar area','left';'Anterior amygdalar area','right';'Central amygdalar nucleus','left';'Central amygdalar nucleus','right';'Medial amygdalar nucleus','left';'Medial amygdalar nucleus','right';'Pallidum','left';'Pallidum','right';'Pallidum/ dorsal region','left';'Pallidum/ dorsal region','right';'Globus pallidus/ internal segment','left';'Globus pallidus/ internal segment','right';'Bed nuclei of the stria terminalis','left';'Bed nuclei of the stria terminalis','right';'Medial habenula','left';'Medial habenula','right';'Lateral habenula','left';'Lateral habenula','right';'Hypothalamus','left';'Hypothalamus','right';'Subthalamic nucleus','left';'Subthalamic nucleus','right';'Zona incerta','left';'Zona incerta','right';'Ventral tegmental area','left';'Ventral tegmental area','right';'Substantia nigra/ compact part','left';'Substantia nigra/ compact part','right';'Cerebellum','left';'Cerebellum','right';'Cerebellar nuclei','left';'Cerebellar nuclei','right'};

Coloc_Normalized_Summary = brain_regions;

for i = 1:length(groups)
    animals = fieldnames(Rabies.(groups{i}));
    for a = 1:length(animals)
        try
            animal = Rabies.(groups{i}).(animals{a});
            summary_column = size(Coloc_Normalized_Summary,2)+1;
            Coloc_Normalized_Summary{1,summary_column} = groups{i};
            Coloc_Normalized_Summary{2,summary_column} = animals{a};

            for j = 2:(length(Coloc_Normalized_Summary)-1)
                Coloc_Normalized_Summary{j+1,summary_column} = [animal.coloc_normalized_data{j,3}];
            end
        catch
            sprintf('Animal did not meet inclusion criteria')
        end
    end
end

%save to Summary struct
Rabies_Summary.Coloc_Normalized.data = Coloc_Normalized_Summary;

%% Group Stats
D1_Healthy_columns = find(strcmp('D1_Healthy',[Coloc_Normalized_Summary(1,1:size(Coloc_Normalized_Summary,2))]));
D1_6OHDA_columns = find(strcmp('D1_6OHDA',[Coloc_Normalized_Summary(1,1:size(Coloc_Normalized_Summary,2))]));
D1_LID_columns = find(strcmp('D1_LID',[Coloc_Normalized_Summary(1,1:size(Coloc_Normalized_Summary,2))]));
A2a_Healthy_columns = find(strcmp('A2a_Healthy',[Coloc_Normalized_Summary(1,1:size(Coloc_Normalized_Summary,2))]));
A2a_6OHDA_columns = find(strcmp('A2a_6OHDA',[Coloc_Normalized_Summary(1,1:size(Coloc_Normalized_Summary,2))]));
A2a_LID_columns = find(strcmp('A2a_LID',[Coloc_Normalized_Summary(1,1:size(Coloc_Normalized_Summary,2))]));
TRAP_LID_columns = find(strcmp('TRAP_LID',[Coloc_Normalized_Summary(1,1:size(Coloc_Normalized_Summary,2))]));

brain_regions = {'Brain Region','Hemisphere';'root','left';'root','right';'Cerebral cortex','left';'Cerebral cortex','right';'Somatomotor areas','left';'Somatomotor areas','right';'Primary motor area','left';'Primary motor area','right';'Secondary motor area','left';'Secondary motor area','right';'Somatosensory areas','left';'Somatosensory areas','right';'Primary somatosensory area','left';'Primary somatosensory area','right';'Primary somatosensory area/ nose','right';'Primary somatosensory area/ barrel field','right';'Primary somatosensory area/ lower limb','right';'Primary somatosensory area/ mouth','right';'Primary somatosensory area/ upper limb','right';'Primary somatosensory area/ trunk','right';'Supplemental somatosensory area','right';'Gustatory areas','right';'Visceral area','right';'Anterior cingulate area','left';'Anterior cingulate area','right';'Anterior cingulate area/ dorsal part','right';'Anterior cingulate area/ ventral part','right';'Orbital area','left';'Orbital area','right';'Agranular insular area','left';'Agranular insular area','right';'Olfactory areas','right';'Striatum','right';'Central amygdalar nucleus','right';'Globus pallidus/ external segment','right';'Pallidum/ ventral region','right';'Thalamus','right';'Ventral group of the dorsal thalamus','right';'Ventral medial nucleus of the thalamus','right';'Thalamus/ polymodal association cortex related','right';'Intralaminar nuclei of the dorsal thalamus','right';'Central medial nucleus of the thalamus','right';'Paracentral nucleus','right';'Parafascicular nucleus','right';'Substantia nigra/ reticular part','right';'Midbrain raphe nuclei','right';'Dorsal nucleus raphe','right';'Primary auditory area','right';'Primary visual area','right';'Prelimbic area','right';'Infralimbic area','right';'Posterior parietal association areas','right';'Piriform area','right';'Hippocampal formation','right';'Lateral amygdalar nucleus','left';'Lateral amygdalar nucleus','right';'Basolateral amygdalar nucleus','left';'Basolateral amygdalar nucleus','right';'Striatum dorsal region','left';'Striatum dorsal region','right';'Striatum ventral region','left';'Striatum ventral region','right';'Nucleus accumbens','left';'Nucleus accumbens','right';'Anterior amygdalar area','left';'Anterior amygdalar area','right';'Central amygdalar nucleus','left';'Central amygdalar nucleus','right';'Medial amygdalar nucleus','left';'Medial amygdalar nucleus','right';'Pallidum','left';'Pallidum','right';'Pallidum/ dorsal region','left';'Pallidum/ dorsal region','right';'Globus pallidus/ internal segment','left';'Globus pallidus/ internal segment','right';'Bed nuclei of the stria terminalis','left';'Bed nuclei of the stria terminalis','right';'Medial habenula','left';'Medial habenula','right';'Lateral habenula','left';'Lateral habenula','right';'Hypothalamus','left';'Hypothalamus','right';'Subthalamic nucleus','left';'Subthalamic nucleus','right';'Zona incerta','left';'Zona incerta','right';'Ventral tegmental area','left';'Ventral tegmental area','right';'Substantia nigra/ compact part','left';'Substantia nigra/ compact part','right';'Cerebellum','left';'Cerebellum','right';'Cerebellar nuclei','left';'Cerebellar nuclei','right'};

Group_summary = brain_regions;
try
    for i = 3:length(Rabies_Summary.Coloc_Normalized.data)
        D1_Healthy(i-2,1) = mean([Rabies_Summary.Coloc_Normalized.data{i,D1_Healthy_columns}]);
        D1_6OHDA(i-2,1) = mean([Rabies_Summary.Coloc_Normalized.data{i,D1_6OHDA_columns}]);
        D1_LID(i-2,1) = mean([Rabies_Summary.Coloc_Normalized.data{i,D1_LID_columns}]);
        A2a_Healthy(i-2,1) = mean([Rabies_Summary.Coloc_Normalized.data{i,A2a_Healthy_columns}]);
        A2a_6OHDA(i-2,1) = mean([Rabies_Summary.Coloc_Normalized.data{i,A2a_6OHDA_columns}]);
        A2a_LID(i-2,1)= mean([Rabies_Summary.Coloc_Normalized.data{i,A2a_LID_columns}]);
        TRAP_LID(i-2,1) = mean([Rabies_Summary.Coloc_Normalized.data{i,TRAP_LID_columns}]);
        
        D1_Healthy(i-2,2) = std([Rabies_Summary.Coloc_Normalized.data{i,D1_Healthy_columns}])/sqrt(length([Rabies_Summary.Coloc_Normalized.data{i,D1_Healthy_columns}]));
        D1_6OHDA(i-2,2) = std([Rabies_Summary.Coloc_Normalized.data{i,D1_6OHDA_columns}])/sqrt(length([Rabies_Summary.Coloc_Normalized.data{i,D1_6OHDA_columns}]));
        D1_LID(i-2,2) = std([Rabies_Summary.Coloc_Normalized.data{i,D1_LID_columns}])/sqrt(length([Rabies_Summary.Coloc_Normalized.data{i,D1_LID_columns}]));
        A2a_Healthy(i-2,2) = std([Rabies_Summary.Coloc_Normalized.data{i,A2a_Healthy_columns}])/sqrt(length([Rabies_Summary.Coloc_Normalized.data{i,A2a_Healthy_columns}]));
        A2a_6OHDA(i-2,2) = std([Rabies_Summary.Coloc_Normalized.data{i,A2a_6OHDA_columns}])/sqrt(length([Rabies_Summary.Coloc_Normalized.data{i,A2a_6OHDA_columns}]));
        A2a_LID(i-2,2)= std([Rabies_Summary.Coloc_Normalized.data{i,A2a_LID_columns}])/sqrt(length([Rabies_Summary.Coloc_Normalized.data{i,A2a_LID_columns}]))
        TRAP_LID(i-2,2) = std([Rabies_Summary.Coloc_Normalized.data{i,TRAP_LID_columns}])/sqrt(length([Rabies_Summary.Coloc_Normalized.data{TRAP_LID_columns}]))
    end
catch
    sprintf('Animal did not meet inclusion criteria')
end

Group_summary(1,3:16) = {'D1_Healthy_Avg','D1_Healthy_SEM','D1_6OHDA_Avg','D1_6OHDA_SEM','D1_LID_Avg','D1_LID_SEM','A2a_Healthy_Avg','A2a_Healthy_SEM','A2a_6OHDA_Avg','A2a_6OHDA_SEM','A2a_LID_Avg','A2a_LID_SEM','TRAP_LID_Avg','TRAP_LID_SEM'}
Group_summary(2:end,3:16) = num2cell([D1_Healthy,D1_6OHDA,D1_LID,A2a_Healthy,A2a_6OHDA,A2a_LID,TRAP_LID]);

%remove contralateral hemisphere from Summary struct
contralateral_indices = find(strcmp('left',[Group_summary(:,2)]))
Group_summary(contralateral_indices,:) = []

%save to Summary struct
Rabies_Summary.Coloc_Normalized.group_stats = Group_summary
