function ui_tile_kymograph_selection(fig,fnames_array)
allAxes = flipud(findall(fig, 'type', 'axes'));
numPlots = numel(fnames_array);

axes_fields = strings(numPlots*2, 1);
for i = 1:numPlots
    axes_fields(2*i-1) = fnames_array{i}+"_bv";
    axes_fields(2*i)   = fnames_array{i}+"_csf";
end

for idx = 1:numel(allAxes)
    ax = allAxes(idx);
    set(ax, 'UserData', axes_fields(idx), ...
            'ButtonDownFcn', @axesClickCallback, ...
            'PickableParts', 'all', 'HitTest', 'on');

    children = ax.Children;
    set(children, 'HitTest', 'off', 'PickableParts', 'none');
end

% 클릭된 항목 저장할 빈 리스트 준비
setappdata(fig, 'clicked_list', strings(0));

% 실시간 표시 UI 생성
listfig = uifigure('Name','Clicked Axes List', ...
                   'Position',[50,50,250,300]);
listbox = uilistbox(listfig,'Position',[10,10,230,280]);

% 클릭 이벤트 콜백
    function axesClickCallback(src,~)
        clicked_fieldname = src.UserData;

        % 뒤에 _bv, _csf 제거
        clean_name = regexprep(clicked_fieldname,'(_bv|_csf)$','');

        % 클릭한 항목 리스트에 추가 or 삭제
        clicked_list = getappdata(fig,'clicked_list');

        if any(clicked_list==clean_name)
            % 이미 리스트에 있으면 제거
            clicked_list(clicked_list==clean_name) = [];
            disp("Removed: " + clean_name);
        else
            % 없으면 추가
            clicked_list(end+1) = clean_name;
            disp("Added: " + clean_name);
        end

        % 리스트 업데이트
        setappdata(fig, 'clicked_list', clicked_list);

        % UI 창에도 리스트 업데이트
        listbox.Items = clicked_list;
        drawnow;
    end
end


