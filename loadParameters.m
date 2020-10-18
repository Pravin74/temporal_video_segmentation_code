%%%%%
%
% File for loading the main parameters for R-Clustering
%version_concepts
%%%%%%%%%%%%%%%%%%%%%%

% clear all, close all

%% Load paths

addpath('Data_Loading','Features_Extraction','Adwin','Concepts','Evaluation','Features_Preprocessing','GraphCuts','PCA','Tests','SpectralClust','Utils')

%% Data loading
% directorio_im = 'D:/LIFELOG_DATASETS'; % SHARED PC
% directorio_im = '/media/HDD_2TB/R-Clustering/Demo/test_data';
% directorio_im = '/media/HDD_2TB/DATASETS/Peleg_data/';
directorio_im = '/media/HDD_3TB/DATASETS/EDUB-Seg';
% directorio_im = ''; % put your own datasets location

% camera = {'SenseCam', 'SenseCam', 'SenseCam', 'SenseCam', 'SenseCam'};
% camera = {'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'SenseCam', 'SenseCam', 'SenseCam', 'SenseCam', 'SenseCam'};
% folders={'Day1','Day2','Day3','Day4','Day6'};
% folders={'Petia1', 'Petia2', 'Estefania1', 'Estefania2', 'Mariella', 'Day1','Day2','Day3','Day4','Day6'};
% folders={'Petia1', 'Petia2', 'Estefania1', 'Estefania2', 'Mariella'};
% formats={'.JPG','.JPG','.JPG','.JPG','.JPG'};
% formats={'.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.JPG','.JPG','.JPG','.JPG','.JPG'};

%%% New Narrative datasets
% camera = {'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative'};
% folders = {'Maya1', 'Marc1'};
% formats = {'.jpg', '.jpg', '.jpg', '.jpg', '.jpg'};
% folders = {'Maya1', 'Maya2', 'Maya3', 'Marc1', 'Estefania3'};
%folders = {'Maya1', 'Maya2', 'Maya3', 'Marc1', 'Estefania3', 'Petia1', 'Petia2', 'Estefania1', 'Estefania2', 'Mariella', 'Day1','Day2','Day3','Day4','Day6'};
%camera = {'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'SenseCam', 'SenseCam', 'SenseCam', 'SenseCam', 'SenseCam'};
% folders = {'Maya1', 'Maya2', 'Maya3', 'Marc1', 'Petia1', 'Petia2', 'Estefania1', 'Estefania2', 'Estefania3', 'Mariella', 'Day1','Day2','Day3','Day4','Day6'};
%formats={'.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.JPG','.JPG','.JPG','.JPG','.JPG'};

%%% EDUB-Seg parts 2 and 3
camera = {'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative'};
formats = {'.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.jpg'};
folders = {'Pedro4', 'Estefania5', 'Marc2', 'Marc3', 'Marc4', 'MarcC1', 'Pedro1', 'Pedro2', 'Pedro3', 'Estefania4'};

% folders = {'Marc4'};

% camera = {'Narrative', 'SenseCam', 'SenseCam', 'SenseCam', 'SenseCam', 'SenseCam'};
% folders = {'Chetan_MergedFolders', 'Yair_MergedFolders'};
% folders = {'Chetan_MergedFolders', 'Yair_MergedFoldersAll'};
% folders = {'FramesVideo_Huji_Chetan_4_Dinner_Part1', 'FramesVideo_Huji_Yair_Sitting_Eating1'};
% folders = {'FramesVideo_Huji_Yair_Sitting_Eating1'};
% formats = {'.jpg','.jpg'};% '.JPG','.JPG','.JPG','.JPG','.JPG'};
% camera = {'Narrative', 'Narrative'};

% camera = {'Narrative'};
% folders = {'Maya1'};
% formats = {'.jpg'};

%directorio_results = 'D:/R-Clustering_Results'; % SHARED PC
%directorio_results = '/media/lifelogging/HDD_2TB/R-Clustering_Data/R-ClusteringResultsMOPCNN';
%directorio_results = '/media/lifelogging/HDD_2TB/R-Clustering_Data/R-Clustering_Results_IBPRIA_new';
%directorio_results = '/media/lifelogging/HDD_2TB/R-Clustering_Data/R-Clustering_Results_concepts_v2';
%directorio_results = '/media/lifelogging/HDD_2TB/R-Clustering_Data/R-Clustering_Results_concepts_v2_allLSDA';
%directorio_results = '/media/lifelogging/HDD_2TB/R-Clustering_Data/R-Clustering_Results_concepts_v3';
%directorio_results = '/media/HDD_2TB/R-Clustering_Data/R-Clustering_Results_short_Peleg';
%directorio_results = '/media/HDD_2TB/R-Clustering_Data/R-Clustering_Results_long_Peleg';
directorio_results = '/media/HDD_2TB/SR-Clustering-master_Data/EDUB-Seg_ALL_LSDA';
%directorio_results = '/media/lifelogging/HDD_2TB/R-Clustering_Data/R-Clustering_Results_concepts_v3_smoothed_50';
%directorio_results = '/Volumes/SHARED HD/R-Clustering Results'; % MARC PC
%  directorio_results = '../Results/Spectral_GC'; % EST PC
% directorio_results = ''; % put your own results location


%% R-Clustering parameters
clus_type = 'Both1'; % Clustering type used before the GraphCuts. 
                        % It can take the following values:
                        %           'Clustering' : Clustering + GC
                        %           'Both1' : Clustering + Adwin + GC
                        %           'Spectral' : Spectral + GC
                        %           'Both2' : Spectral + Adwin + GC
                        
paramsPCA.minVarPCA=0.95;
paramsPCA.standarizePCA=false;
paramsfeatures.type = 'CNNconcepts'; %'CNNconcepts';%MOPCNN'; %CNN ....    

% In case of using concept-based features
version_concepts = 2;

%% Clustering parameters
methods_indx={'ward','centroid','complete','weighted','single','median','average'};
% methods_indx={'centroid'};
cut_indx=(0.2:0.2:1.2);
% cut_indx = [0.45];
paramsPCA.usePCA_Clustering = true;

%% Spectral Clustering 
%paramsSpec.NN = false;
%paramsSpec.Sig = true;
%paramsSpec.Eps = false;

NN=5; 
Sig=1; 
Eps=1; 

sim_matrix={'Sigma','NN','Epsilon'};

k_valuesSp=6:4:36;

paramsPCA.usePCA_Spect = true;

%% Adwin parameters
pnorm = 2;
confidence = 0.1;
paramsPCA.usePCA_Adwin = true;

%% GraphCuts parametersData_Loading
paramsPCA.usePCA_GC = false;
window_len = 11;
% window_len = 50;

W_unary = 0.1;      % 0 <= W_unary <= 1 for evalType == 1
W_pairwise = 0.5;   % 0 <= W_pairwise <= 1 for evalType == 1

nUnaryDivisions = 11; % number of equally spaces W_unary values for evalType == 2
nPairwiseDivisions = 11; % number of equally spaced W_pairwise values for evalType == 2

evalType = 2; % 1 = single test, 2 = iterative W increase
doEvaluation = true; % plot precision/recall and f-measure when performing single test

%% Evaluation parameters
tol=5; % tolerance for the final evaluation

% tol = 10; %%%% REMOVE!!!

minImCl=0; % (deprecated)
plotFigResults = false;



% Additional parameters
min_length_merge = 0 % (default 0) minimum length of segments, if smaller, will be merged to most similar adjacent segment.
