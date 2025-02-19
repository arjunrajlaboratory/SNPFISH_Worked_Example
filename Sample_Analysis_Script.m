%% Worked Example Sample Analysis Script

% This sample contains a sample SNP FISH gene with two alleles.
% The gene in this case is Dusp6 in a mouse embryonic fibroblast.
% The guide probe is in the gfp wavelength, and the two alleles are
% in the Cy3 and Cy5 wavelenghts. This sample was selected because you 
% can observe transcriptional bursting from both alleles, one colocalizing
% with one allele and the other with the other. 


%% Read YAML FILE
expData = ReadYaml('readme.yaml');



%% Sample snpMap Definition
snpMap.channels = {'gfp', 'tmr', 'cy'};
snpMap.names = {expData.channels.gfp.probe, ...
                expData.channels.tmr.probe, ...
                expData.channels.cy.probe,};

%% Gaussian Fitting

dataAdder = improc2.processing.DataAdder();
unprocessedFittedData = improc2.nodeProcs.TwoStageSpotFitProcessedData();

dataAdder.addDataToObject(unprocessedFittedData, 'gfp', 'gfp:Fitted')
dataAdder.addDataToObject(unprocessedFittedData, 'tmr', 'tmr:Fitted')
dataAdder.addDataToObject(unprocessedFittedData, 'cy', 'cy:Fitted')

dataAdder.repeatForAllObjectsAndQuit();

improc2.processing.updateAll()


%% Perform SNP Colocalization

dataAdder = improc2.processing.DataAdder();
dataAdder.addDataToObject(improc2.nodeProcs.SNPColocalizer(snpMap, ...
    'zAllow', 3), snpMap.channels, 'snpColoc2'); % IF YOU ARE RUNNING THE ANALYSIS YOURSELF AND WANT TO ADD A NEW NODE, CHANGE THE LABEL FROM 'snpColoc2' to 'YOURNAMEHERE' or whatever you like. Note: the node with label 'snpColoc' is an old version of the SNPColocalizer node analysis
dataAdder.repeatForAllObjectsAndQuit();

improc2.processing.updateAll();


%% Create MOLTEN Data tables for Export

tools = improc2.launchImageObjectTools();


cellCounts = [];
tools.iterator.goToFirstObject
cellID = 1;

guideAll = []; 
cyAll = [];
tmrAll = [];

while(tools.iterator.continueIteration)
    results = tools.objectHandle.getData('snpColoc2');
    
    guide = struct2table(results.data.gfp);
    guide.cellID = ones(height(guide),1) * cellID;
    cy = struct2table(results.data.cy);
    cy.cellID = ones(height(cy),1) * cellID;
    tmr = struct2table(results.data.tmr);
    tmr.cellID = ones(height(tmr),1) * cellID;
    
    if isempty(guideAll)
    guideAll = guide;
    else
    guideAll = vertcat(guideAll, guide);
    end
    
    if isempty(cyAll)
    cyAll = cy;
    else
    cyAll = vertcat(cyAll, cy);
    end
    
    if isempty(tmrAll)
        tmrAll = tmr; 
    else
    tmrAll = vertcat(tmrAll, tmr);
    end
    
    tools.iterator.goToNextObject;
    cellID = cellID + 1;
end


%% Export to CSV

writetable(guideAll, 'sample_SNPFISH_genetable.csv')

%% Visualize Complete SNPFISH Object data structure

tools.objectHandle.view


%% View Sample Overlay


fighandle = figure;

img = tools.objectHandle.getData([snpMap.channels{1}, ':Spots']).getImage;
imgmax = max(img,[],3);  %MAX MERGE
imgmax = imadjust(imgmax, stretchlim(imgmax, [0 0.995]));


figure(fighandle)
imshow(imgmax,[]);  %Plot image
hold on
plot(guide.position(:,1), guide.position(:,2), 'wo'); 
plot(guide.position(guide.labels == snpMap.names{2},1), guide.position(guide.labels == snpMap.names{2},2),'co','markersize',6); %Plot the co_localized spots
plot(guide.position(guide.labels ==  snpMap.names{3},1), guide.position(guide.labels == snpMap.names{3},2),'Color', [1, 0.65, 0.2], ...
    'LineStyle', 'none', 'Marker','o','markersize',6); %Plot the co_localized spots
plot(guide.position(guide.labels == '3-color',1), guide.position(guide.labels == '3-color',2),'mo','markersize',6);
legend({snpMap.names{1}, snpMap.names{2},snpMap.names{3},'3-color'})
hold off