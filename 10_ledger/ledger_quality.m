function record = ledger_quality(sessiondir, grade, note)
%LEDGER_QUALITY  One session's grade, read or written, keyed the way every table is.
%   Called from review_session with the figures still open: look, then grade. The
%   numbers a table already carries are NOT what this is for -- bv_range says the
%   diameter span and cond_sd says the scatter, and re-typing those as a grade only
%   makes a second copy that can disagree. What goes in the note is what no column
%   can see: where the kymograph drifts, which boundary the FWHM loses, whether the
%   vessel sits at an angle.
%
%   THE NOTE IS WRITTEN IN ENGLISH. The register is read alongside code and its
%   lines end up in figure captions, and MATLAB writes csv through the system
%   codepage on Windows unless told otherwise.
%
%   S  goes in a representative figure. Clean from start to end, boundaries steady
%   A  used in the analysis as it is. Local defects that do not move the measurement
%   B  used, never representative. The defect is visible and the note says where
%   C  dropped once more data arrives. Kept now because the sample is worth having
%
%   IN   sessiondir  1x1 str      the session folder, as pasted into review_session
%        grade       1x1 str      "S" | "A" | "B" | "C". Omit both to read
%        note        1x1 str      why, in English. Required with a grade -- a bare
%                                 letter sends the next reader back to the figure,
%                                 and that second look IS the whole of the grade
%   OUT  record      0x0 | 1-row table   what the register holds for this session
%
%   Example
%     ledger_quality(sessiondir)
%     ledger_quality(sessiondir, "C", "kymograph drifts after 20 min, lower BV edge lost")
    arguments
        sessiondir (1,1) string
        grade      (1,1) string {mustBeMember(grade, ["", "S", "A", "B", "C"])} = ""
        note       (1,1) string = ""
    end

    dirs = dirs_central();
    dataset = string(getenv('ECSPRESS_DATASET'));
    if dataset == ""
        dataset = "merged_igkl_igkltdt";
    end
    register_path = fullfile(dirs.secondary_root, dataset, "session_quality.csv");
    register = read_register(register_path);

    identity = session_identity(dirs.secondary_root, dataset, sessiondir);
    on_session = register.key == identity.key;

    if grade == ""
        record = register(on_session, :);
        report(record, identity, register);
        return
    end
    if note == ""
        error('ledger_quality:noNote', ...
            'grade %s needs a note. what does the figure show that the columns cannot', grade);
    end

    row = identity;
    row.grade = grade;
    row.note = note;
    row.graded_on = datetime('now', 'Format', 'yyyy-MM-dd HH:mm');
    if any(on_session)
        fprintf('%s was %s, now %s\n', identity.key, register.grade(on_session), grade);
        register(on_session, :) = row;
    else
        register = [register; row];
    end
    write_register(register, register_path);
    record = row;
end

function identity = session_identity(secondary_root, dataset, sessiondir)
%SESSION_IDENTITY  The four key columns for this folder, looked up rather than parsed.
%   The dirtable already split every path into MouseID / Date / SessionType /
%   SessionID once. Parsing the folder name a second time here would be a second
%   rule for the same thing, and the two would drift.
    sheet_path = fullfile(secondary_root, dataset, dataset + "_dirtable.xlsx");
    if ~isfile(sheet_path)
        error('ledger_quality:noDirtable', ...
            'no %s. run centralize_primary first', sheet_path);
    end
    dirtable = readtable(sheet_path, 'Sheet', 'reference', 'TextType', 'string');
    wanted = strip(string(sessiondir), 'right', filesep);
    on_row = strcmpi(strip(string(dirtable.Directory), 'right', filesep), wanted);
    if ~any(on_row)
        error('ledger_quality:noSession', ...
            '%s is not a row of %s', sessiondir, sheet_path);
    end
    identity = dirtable(find(on_row, 1), ["MouseID", "Date", "SessionType", "SessionID"]);
    identity.Directory = wanted;
    identity.key = util_sessionkey(identity);
    identity = identity(:, ["key", "MouseID", "Date", "SessionType", "SessionID", "Directory"]);
end

function register = read_register(register_path)
%READ_REGISTER  The register as it stands, or an empty one with the right columns.
    if isfile(register_path)
        % every column read back as TEXT. SessionID is "005" and readtable would
        % otherwise infer a double from it and hand back 5, which no longer matches
        % the key the rest of the pipeline joins on
        opts = detectImportOptions(register_path, 'TextType', 'string');
        opts = setvartype(opts, opts.VariableNames, 'string');
        register = readtable(register_path, opts);
        register.graded_on = datetime(register.graded_on, 'Format', 'yyyy-MM-dd HH:mm');
        return
    end
    register = table(strings(0,1), strings(0,1), strings(0,1), strings(0,1), ...
        strings(0,1), strings(0,1), strings(0,1), strings(0,1), ...
        datetime.empty(0,1), 'VariableNames', ["key", "MouseID", "Date", ...
        "SessionType", "SessionID", "Directory", "grade", "note", "graded_on"]);
end

function write_register(register, register_path)
%WRITE_REGISTER  Backed up first. This is hand-scored work, like sleep_score.
    if isfile(register_path)
        stamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmm'));
        backup_path = replace(register_path, ".csv", "_backup_" + stamp + ".csv");
        copyfile(register_path, backup_path);
    end
    writetable(register, register_path, 'Encoding', 'UTF-8');
    fprintf('%d sessions graded, %s\n', height(register), register_path);
end

function report(record, identity, register)
%REPORT  What is held for this session, and the shape of the register around it.
    if isempty(record)
        fprintf('%s  not graded yet\n', identity.key);
    else
        fprintf('%s  %s  %s  (%s)\n', record.key, record.grade, record.note, ...
            string(record.graded_on));
    end
    if isempty(register)
        return
    end
    counted = strings(1, 0);
    for letter = ["S", "A", "B", "C"]
        counted(end+1) = letter + " " + sum(register.grade == letter); %#ok<AGROW>
    end
    fprintf('   register holds %d : %s\n', height(register), strjoin(counted, " | "));
end
