function session_path = dirs_sessionpath(volume_serial, directory)
%DIRS_SESSIONPATH  Today's path of a session the dirtable recorded by volume serial and Directory.
%   Caller: dirs_hashstate
%   Directory holds the letter of the day it was scanned; only its part below \tmp\ is trusted.
%
% IN   volume_serial  N x 1 str   the VolumeSerial column
%      directory      N x 1 str   the Directory column
% OUT  session_path   N x 1 str   <letter today>\tmp\<part below tmp>, "" where that volume is not
%                                 mounted or the serial is blank

volume = dirs_volumeletter();
volume_serial = string(volume_serial);
directory = string(directory);
session_path = strings(size(directory));
for k = 1:numel(directory)
    letter = volume.letter(volume.serial == volume_serial(k));
    if isempty(letter) || strlength(letter) == 0
        continue
    end
    relative = extractAfter(directory(k), '\tmp\');
    session_path(k) = fullfile(letter, "tmp", relative);
end
end
