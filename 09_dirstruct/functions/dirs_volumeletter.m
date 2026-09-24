function volume = dirs_volumeletter()
%DIRS_VOLUMELETTER  Every fixed disk mounted today, with the role its beacon file declares.
%   Caller: dirs_central, dirs_mapstruct, dirs_sessionpath
%   A data volume carries <letter>\dirstruct_beacon.txt holding one word: work, archive
%   or mirror. Several disks may share a role; the beacon sits at the root, which
%   dirs_mirror never copies. No letter, serial or label is in code.
%
% OUT  volume   table   letter | serial | label | role, one row per fixed disk;
%                       role is "" for a disk without a beacon

beacon_name = "dirstruct_beacon.txt";
row = struct('letter', {}, 'serial', {}, 'label', {}, 'role', {});

% asked through .NET, no shell
NET.addAssembly('System.Management');
query = 'SELECT DeviceID, VolumeSerialNumber, VolumeName FROM Win32_LogicalDisk WHERE DriveType = 3';
searcher = System.Management.ManagementObjectSearcher(query);
disks = searcher.Get().GetEnumerator();
while disks.MoveNext()
    disk = disks.Current;
    letter = string(char(disk.Item('DeviceID')));
    serial = string(char(disk.Item('VolumeSerialNumber')));
    label = string(char(disk.Item('VolumeName')));
    role = "";
    beacon_path = fullfile(letter + "\", beacon_name);
    if isfile(beacon_path)
        role = strtrim(string(fileread(beacon_path)));
    end
    row(end + 1) = struct('letter', letter, 'serial', serial, 'label', label, 'role', role); %#ok<AGROW>
end
volume = struct2table(row);
end
