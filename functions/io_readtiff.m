function [stack, info] = io_readtiff(tiff_path)
%IO_READTIFF  .tiff to stack, plus what its tags say
%   Reads back what 00_mdfExtractor writes, pages in file order. A page loop and not
%   tiffreadVolume, which is far slower on a long stack; see CLAUDE_LOG.md
%   Caller: piv_zstack
%
% IN   tiff_path         1 x n char       the TIFF
% OUT  stack             H x W x N        every page, in the class the file stores
%      info.n_page       1 x 1 double     pages in the file
%      info.resolution   1 x 3 double     [x y z] step per pixel and per page; x, y as the TIFF RATIONAL stores them, z to six significant digits
%      info.unit         1 x 3 string     the unit of each, as the file says
%      info.description  1 x 1 struct     the ImageJ block, one field per key
    page_info = imfinfo(tiff_path);
    first_tags = page_info(1);
    info.n_page = numel(page_info);
    info.description = parse_description(first_tags.ImageDescription);
    % mdf.tiftag stores um as pixels per cm, io_postsavetiff pixels per unit with no ResolutionUnit
    switch first_tags.ResolutionUnit
        case 'Centimeter'
            xy_step = 10000 ./ [first_tags.XResolution, first_tags.YResolution];
        case 'None'
            xy_step = 1 ./ [first_tags.XResolution, first_tags.YResolution];
        otherwise
            error('io_readtiff:resolutionUnit', '%s stores its resolution per %s', tiff_path, first_tags.ResolutionUnit);
    end
    % a z page axis carries spacing, a time axis finterval; both writers put zunit beside either
    if isfield(info.description, 'spacing')
        z_step = info.description.spacing;
    else
        z_step = info.description.finterval;
    end
    info.resolution = [xy_step, z_step];
    % yunit entered the writer on 2026-09-02; before that, and in io_postsavetiff, y shares x's unit
    if isfield(info.description, 'yunit')
        y_unit = info.description.yunit;
    else
        y_unit = info.description.unit;
    end
    info.unit = string({info.description.unit, y_unit, info.description.zunit});

    tiff_obj = Tiff(tiff_path, 'r');
    first_page = tiff_obj.read();
    [H, W] = size(first_page);
    stack = zeros(H, W, info.n_page, 'like', first_page);
    stack(:, :, 1) = first_page;

    report_every = max(1, floor(info.n_page / 10));
    for page_idx = 2:info.n_page
        tiff_obj.nextDirectory();
        stack(:, :, page_idx) = tiff_obj.read();
        if mod(page_idx, report_every) == 0
            fprintf('%.0f%% read\n', page_idx * 100 / info.n_page);
        end
    end
    tiff_obj.close();
    fprintf('read %d pages %d x %d %s from %s\n', info.n_page, size(stack, 1), size(stack, 2), class(stack), tiff_path);
end

function description = parse_description(text)
%PARSE_DESCRIPTION  ImageJ's key=value lines as a struct; a value that reads as a number is one
%
% IN   text         1 x n char     the ImageDescription tag
% OUT  description  1 x 1 struct   one field per key
    description = struct();
    line_list = splitlines(string(text));
    for description_line = line_list'
        if ~contains(description_line, '=')
            continue
        end
        key = char(extractBefore(description_line, '='));
        value_text = extractAfter(description_line, '=');
        value = str2double(value_text);
        if isnan(value)
            description.(key) = char(value_text);
        else
            description.(key) = value;
        end
    end
end
