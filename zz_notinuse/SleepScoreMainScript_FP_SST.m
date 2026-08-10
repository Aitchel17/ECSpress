%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-TurnerScoringResults
%________________________________________________________________________________________________________________________
%
% Purpose:
%________________________________________________________________________________________________________________________
% addpath(genpath('C:\Users\mfh5734\OneDrive - The Pennsylvania State University\Documents\Research_Codes\FiberPhotometry\Data-Analysis-master'))
% function SleepScoreMainScript_FP_GRABNE()
% Clear workspace/Load in file names for various analysis
clear; clc; close all
disp('Loading necessary file names...'); disp(' ')
baselineType = 'manualSelection';
startingDirectory = cd;
% Create training data set for each animal
% cd to the animal's bilateral imaging folder to load the baseline structure
baselineDirectory = [startingDirectory];% '\Bilateral Imaging\'];
cd(baselineDirectory)       
% load the baseline structure
baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}';
baselinesFileID = char(baselinesFile);
load(baselinesFileID)
cd(startingDirectory)
% cd to the animal's training set folder
trainingDirectory = [startingDirectory]; % '\Training Data\'];
cd(trainingDirectory)
%% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% measure length of the matrix based on the length of the data
load(procDataFileIDs(1,:))
trialDuration = ProcData.notes.trialDuration_sec;
BinSize = 5; % 5 seconds bins
NBins = trialDuration/BinSize; % duration divided by 5 seconds
%%
% add sleep parameters (each behavior we care about during sleep)
AddSleepParameters_FP_GRABNE_SingleFiber(procDataFileIDs,RestingBaselines,baselineType,NBins)
%AddSleepParameters_FP_GRABNE_SingleFiber_modifiedIndexing(procDataFileIDs,RestingBaselines,baselineType,NBins)
% create a table of values for sleep scoring model
CreateModelDataSet_FP_GRABNE_SingleFiber(procDataFileIDs,NBins)
%%
% create manual decisions for each 5 second bin
CreateTrainingDataSet_FP_GRABNE(procDataFileIDs,NBins)
% combine the existing training set decisions with any sleep parameter changes
UpdateTrainingDataSets_FP(procDataFileIDs)
cd(startingDirectory)
%% create microarousal labels for NREM Sleep Only
% CreateMicroArousalDataSet_FP_GRABNE(procDataFileIDs)
%% Change Arousal to Microarousal if the arousal is less than 20s
% AddMicroArousalsParameters_FP_GRABNE_SingleFiber(procDataFileIDs,RestingBaselines,baselineType,trialDuration)
%% 
% MicroArousalTimeMax = 1;
% MicroArousalTimeMin = 1;
% MicroArousal_DataSeparation(startingDirectory,trainingDirectory,baselineDirectory,MicroArousalTimeMax,MicroArousalTimeMin);
%%
% Train Models - cycle through each data set and update any necessary parameters
%  [animalID] = TrainSleepModels_FP_GRABNE;%TrainSleepModels_FP_Shak_opt;
[animalID] = GetAnimalID();
%% Sleep score an animal's data set and create a SleepData.mat structure for classification
modelNames = {'Manual'};%{'Forest','SVM','Manual','Ensemble'};
SleepData = [];
% cd to the animal's bilateral imaging folder
cd(baselineDirectory)
% load the baseline structure
baselinesFileStruct = dir('*_RestingBaselines.mat');
baselinesFile = {baselinesFileStruct.name}; 
baselinesFileID = char(baselinesFile);
load(baselinesFileID)
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% add sleep parameters (each behavior we care about during sleep)
% AddSleepParameters_FP_GRABNE(procDataFileIDs,RestingBaselines,baselineType,NBins)
% create a table of values for sleep scoring model
% CreateModelDataSet_FP_GRABNE(procDataFileIDs,NBins)
% character list of all ModelData files
modelDataFileStruct = dir('*_ModelData.mat');
modelDataFiles = {modelDataFileStruct.name}';
modelDataFileIDs = char(modelDataFiles);
for c = 1:length(modelNames)
    modelName = modelNames{1,c};
    [ScoringResults] = PredictBehaviorEvents_FP_GRABNE(animalID,startingDirectory,baselineDirectory,modelDataFileIDs,modelName);
    ApplySleepLogical_FP(startingDirectory,trainingDirectory,baselineDirectory,modelName,ScoringResults)
    NREMsleepTime = 30;   % seconds
    REMsleepTime = 60;   % seconds
    [SleepData] = CreateSleepData_FP_GRABNE_Auto(startingDirectory,trainingDirectory,baselineDirectory,NREMsleepTime,REMsleepTime,modelName,SleepData);
end 
cd(baselineDirectory)
save([animalID '_SleepData.mat'],'SleepData')
cd(startingDirectory)

disp('Sleep Scoring analysis complete'); disp(' ')
