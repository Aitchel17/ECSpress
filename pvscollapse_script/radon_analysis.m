% load file
savepath = 'E:\OneDrive - The Pennsylvania State University\2023_ECSpress\01_primary_analysis\matlab_analysis';
file1 = ECSanalysis(savepath);
file1.analog = file1.loadanalog;
file1.stackch1 = file1.loadstack("ch1");
file1.stackch2 = file1.loadstack("ch2");


