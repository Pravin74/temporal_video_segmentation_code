%%%%%
%
% File for loading the main parameters for R-Clustering
%
%%%%%%%%%%%%%%%%%%%%%%

% clear all, close all

%% Load paths

addpath('/media/lifelogging/HDD_2TB/R-Clustering/Data_Loading','/media/lifelogging/HDD_2TB/R-Clustering/Features_Extraction','/media/lifelogging/HDD_2TB/R-Clustering/Adwin','/media/lifelogging/HDD_2TB/R-Clustering/Concepts','/media/lifelogging/HDD_2TB/R-Clustering/Evaluation','/media/lifelogging/HDD_2TB/R-Clustering/Features_Preprocessing','/media/lifelogging/HDD_2TB/R-Clustering/GraphCuts','/media/lifelogging/HDD_2TB/R-Clustering/PCA','/media/lifelogging/HDD_2TB/R-Clustering/Tests','/media/lifelogging/HDD_2TB/R-Clustering/SpectralClust','/media/lifelogging/HDD_2TB/R-Clustering/Utils')

%% Data loading
% directorio_im = 'D:/LIFELOG_DATASETS'; % SHARED PC
directorio_im = '/media/lifelogging/HDD_2TB/LIFELOG_DATASETS'; % SHARED PC
% directorio_im = '/Volumes/SHARED HD/Video Summarization Project Data Sets/R-Clustering'; % MARC PC
% directorio_im='/Users/estefaniatalaveramartinez/Desktop/LifeLogging/IbPRIA/Sets'; % EST PC
% directorio_im = ''; % put your own datasets location

camera = {'Narrative', 'Narrative', 'Narrative', 'Narrative', 'Narrative', 'SenseCam', 'SenseCam', 'SenseCam', 'SenseCam', 'SenseCam'};
folders={'Estefania1', 'Estefania2', 'Petia1', 'Petia2', 'Mariella', 'Day1','Day2','Day3','Day4','Day6'};
% folders={'Petia2'};
formats={'.jpg', '.jpg', '.jpg', '.jpg', '.jpg', '.JPG','.JPG','.JPG','.JPG','.JPG'};

%directorio_results = 'D:/R-Clustering_Results'; % SHARED PC
%directorio_results = '/media/lifelogging/HDD_2TB/R-ClusteringResultsMOPCNN';
directorio_results = '/media/lifelogging/HDD_2TB/R-Clustering_Results_concepts_v3';
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
paramsfeatures.type = 'CNNconcepts';%MOPCNN'; %CNN ....                       
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

%% GraphCuts parameters
paramsPCA.usePCA_GC = false;
window_len = 11;

W_unary = 0.1;      % 0 <= W_unary <= 1 for evalType == 1
W_pairwise = 0.5;   % 0 <= W_pairwise <= 1 for evalType == 1

nUnaryDivisions = 5; % number of equally spaces W_unary values for evalType == 2
nPairwiseDivisions = 5; % number of equally spaced W_pairwise values for evalType == 2

evalType = 2; % 1 = single test, 2 = iterative W increase
doEvaluation = true; % plot precision/recall and f-measure when performing single test

%% Evaluation parameters
tol=5; % tolerance for the final evaluation
minImCl=0; % (deprecated)
plotFigResults = false;
