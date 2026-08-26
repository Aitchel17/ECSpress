function [row, failure] = centralize_scansession(i, source_path, product, old)
%CENTRALIZE_SCANSESSION  One session's row for one product, or 0x0 if it has none.
%   IN   i            1x1 double      which row of live_table this is
%        source_path  1x1 str         the .mat under the session
%        product      1x1 str         which stripper its contents take
%        old          0|1 row table   what the last run wrote for this session
%   OUT  row          0x0 | 1x1 struct  the five fields are declared below
%        failure      1x1 str         "" unless the load threw
    row = struct('live_i', {}, 'data', {}, 'source_bytes', {}, ...
        'source_modified', {}, 'reused', {});
    failure = "";
    if ~isfile(source_path)
        return
    end
    listing = dir(source_path);

    % the source has not moved since the last run read it, so its row still stands
    unmoved = ~isempty(old) && old.source_bytes == listing.bytes ...
        && old.source_modified == listing.datenum;
    if ~unmoved
        try
            loaded = load(source_path);
            field_list = fieldnames(loaded);
            if isscalar(field_list)
                content = loaded.(field_list{1});
            else
                content = loaded;
            end
            if isnumeric(content) || isempty(content)
                error('centralize_scansession:classMissing', ...
                    'loaded as %s %s, not an object -- its class is not on the path: %s', ...
                    mat2str(size(content)), class(content), source_path);
            end
            data = centralize_strip(product, content);
            row(1).live_i = i;
            row(1).data = data;
            row(1).source_bytes = listing.bytes;
            row(1).source_modified = listing.datenum;
            row(1).reused = false;
            return
        catch err
            % a source that will not load is no reason to drop what was already
            % extracted from it, so the old row stands and the failure is named
            failure = string(err.message);
        end
    end
    if isempty(old)
        return
    end
    row(1).live_i = i;
    row(1).data = old.data{1};
    row(1).source_bytes = old.source_bytes;
    row(1).source_modified = old.source_modified;
    row(1).reused = true;
end
