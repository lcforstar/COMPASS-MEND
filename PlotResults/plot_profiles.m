function fig=plot_profiles(outputDir,acceptedStep,visibility)
%PLOT_PROFILES View five native model-state profiles from a completed output.
% fig=plot_profiles(outputDir) displays the last written accepted row.
% fig=plot_profiles(outputDir,acceptedStep,'off') renders without a window.
% Variables retain the units and node/layer geometry stored in NetCDF.
% This viewer does not apply the observation-matching or20cm dry-mass operator
% and does not calculate observation statistics. A row index is labelled as
% an accepted KT; inherited NetCDF time-coordinate conventions are not changed.
% No model or site file is read or changed. PNG/PDF go to outputDir/figures.
% Editor: Wenjie Liu; Date: 2026-10-31.
    if nargin<1||isempty(outputDir)
        root=fileparts(fileparts(mfilename('fullpath')));
        outputDir=fullfile(root,'data','output','changwu');
    end
    if nargin<2,acceptedStep=[];end
    if nargin<3,visibility='on';end
    groups={'soil_water','soil_nitrogen','soil_nitrogen','soil_organic','soil_organic'};
    variables={'soil_moisture','nh4_total','no3_total','soc_observation_profile','doc_carbon_profile'};
    labels={'Water','NH4-N','NO3-N','SOC','DOC (free D)'};
    values=cell(5,1);depth=cell(5,1);units=cell(5,1);steps=zeros(5,1);
    for k=1:5
        file=fullfile(outputDir,'nc',[groups{k} '.nc']);
        assert(isfile(file),'plotProfiles:missingOutput','Missing native output: %s',file);
        info=ncinfo(file,variables{k});names={info.Dimensions.Name};
        t=find(strcmp(names,'time'));
        assert(isscalar(t));
        steps(k)=info.Dimensions(t).Length;
        if isempty(acceptedStep),whichStep=steps(k);else,whichStep=acceptedStep;end
        assert(isscalar(whichStep)&&isfinite(whichStep)&&whichStep==fix(whichStep) ...
            &&whichStep>=1&&whichStep<=steps(k),'plotProfiles:step','Requested accepted row is unavailable.');
        first=ones(1,numel(names));count=[info.Dimensions.Length];
        first(t)=whichStep;count(t)=1;
        values{k}=double(squeeze(ncread(file,variables{k},first,count)));
        if any(strcmp(names,'node')),coordinate='depth_node';
        elseif any(strcmp(names,'layer')),coordinate='depth_layer';
        else,error('plotProfiles:geometry','No native node/layer dimension.');end
        depth{k}=double(ncread(file,coordinate));
        assert(isvector(values{k})&&numel(values{k})==numel(depth{k}) ...
            &&all(isfinite(values{k}))&&all(diff(depth{k})>0));
        units{k}=char(ncreadatt(file,variables{k},'units'));
    end
    assert(all(steps==steps(1)),'plotProfiles:unequalEndpoints','Groups have different written endpoints.');
    if isempty(acceptedStep),acceptedStep=steps(1);end
    fig=figure('Color','w','Visible',visibility,'Position',[100 100 1450 650]);
    layout=tiledlayout(fig,1,5,'TileSpacing','compact','Padding','compact');
    for k=1:5
        ax=nexttile(layout);plot(ax,values{k},depth{k},'-','Color',[0 .447 .698],'LineWidth',1.5);
        set(ax,'YDir','reverse');grid(ax,'on');title(ax,labels{k});
        xlabel(ax,units{k},'Interpreter','none');
        if k==1,ylabel(ax,'Depth (cm)');end
    end
    title(layout,sprintf('Changwu winter wheat | accepted KT %d | native output units',acceptedStep));
    folder=fullfile(outputDir,'figures');if ~isfolder(folder),mkdir(folder);end
    stem=fullfile(folder,sprintf('native_profiles_KT%05d',acceptedStep));
    exportgraphics(fig,[stem '.png'],'Resolution',200);
    exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
end
