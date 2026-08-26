function [score, primary_analog, stamp] = centralize_scoring(score_table, key_score, ...
    analog_table, key_analog, session_key, pax_stamp)
%CENTRALIZE_SCORING  Which scoring splits this session, and when its sources moved.
%   A sleep score wins. Almost every recording carries an analysis_analog.mat,
%   sleep ones included, and theirs has an empty air table -- running
%   awake_integration over it asks a table with no columns for Duration. So an
%   analog counts as a scoring only when its air table has rows.
%
%   The stamp covers the sources the row is actually built from, and only those: a
%   sleep recording that also carries an analog analysis does not read it, so a
%   change there must not force a recompute.
%   IN   score_table     table       centralized_sleep_score
%        key_score       Nx1 str     util_sessionkey over it
%        analog_table    table       centralized_analysis_analog
%        key_analog      Nx1 str     util_sessionkey over it
%        session_key     1x1 str     the session being asked about
%        pax_stamp       1x1 double  source_modified of its paxfwhm
%   OUT  score           1x1 struct | []
%        primary_analog  1x1 struct | []
%        stamp           1x1 double | []   [] = nothing to split this session by
    score = [];
    primary_analog = [];
    stamp = [];

    hit_score = find(key_score == session_key, 1);
    if ~isempty(hit_score)
        score = score_table.data{hit_score};
        stamp = pax_stamp + score_table.source_modified(hit_score);
        return
    end

    hit_analog = find(key_analog == session_key, 1);
    if isempty(hit_analog)
        return
    end
    candidate = analog_table.data{hit_analog};
    if ~isfield(candidate, 'airtable') || isempty(candidate.airtable)
        return
    end
    primary_analog = candidate;
    stamp = pax_stamp + analog_table.source_modified(hit_analog);
end
