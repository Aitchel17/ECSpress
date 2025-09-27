roilist = primary_datastruct.roilist;
roilist = roilist.addroi(preprocessed_ch1,'extraparenchyma','polygon');
roilist = roilist.modifyroi(preprocessed_ch2,'extraparenchyma');
roilist = roilist.modifyroi(preprocessed_ch1,'extraparenchyma');
%%
save(fullfile(savepath,"roilist.mat"),'roilist')
