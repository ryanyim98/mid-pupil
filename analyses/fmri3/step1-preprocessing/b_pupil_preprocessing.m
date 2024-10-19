clear all; clc;


cd('~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/derivatives');
addpath(genpath('~/Desktop/VRMID-analysis/mid-pupil/analyses/functions/functions-fmri'));
physio_data = 'all_subjs_pupil_preMATLAB.csv';
T = readtable(physio_data);
T.ps = T.ps/1000;

% T.xp = str2double(T.xp);
% T.yp = str2double(T.yp);
% T.time = str2double(T.time);
%%
subplot(1,3,1);
histogram(T.ps);

subplot(1,3,2);
histogram(T.xp);

subplot(1,3,3);
histogram(T.yp);
%%
hist3([T.xp T.yp],[101 101]);
% imagesc(values.')
% colorbar

%% convert file structure
clear alldata;
subjects_list = readtable('../subjects_list/subjects-pupil.txt','ReadVariableNames',false);
subjects = [];

for s = 1:length(subjects_list.Var1)
    subjects = [subjects;char(subjects_list.Var1(s)),'b1'];
    subjects = [subjects;char(subjects_list.Var1(s)),'b2'];
    subjects = [subjects;char(subjects_list.Var1(s)),'b3'];
    subjects = [subjects;char(subjects_list.Var1(s)),'b4'];
end

nsubjects = size(subjects,1)

%
sr_old = 1000;

for s = 1:length(subjects)
    sub = subjects(s,:);
    disp(["processing participant " + sub]);

    %downsample the physio data!
    temp = T(string(cell2mat(T.subject)) == sub,:);

    alldata.sampling_rate = sr_old;

    alldata.subjectdata(s).ID = sub;

    pupil_position = temp(:,["xp","yp"]);
    pupil_size = temp(:,["ps","blink"]);
    beh = temp(:,["trial","blink","sacc","sacc_ampl","fix","event"]);

    alldata.subjectdata(s).beh = beh;
    alldata.subjectdata(s).Physio.pupil_position = pupil_position;
    alldata.subjectdata(s).Physio.pupil_size = pupil_size;

    times = temp.time;
    times_conv = (times - min(times))/sr_old;
    alldata.subjectdata(s).nTimePoints = length(times);
    alldata.subjectdata(s).total_time = (times(length(times)) - times(1))/sr_old;
    alldata.subjectdata(s).beh.times = times_conv;

end

clc;

%% plot subj's blinking to get a sense of the amount of gap expansion needed

my_time = [1,10];
my_time = (sr_old*my_time(1)):(sr_old*my_time(2));
for sub = 1:3
    ps = alldata.subjectdata(sub).Physio.pupil_size.ps(my_time);
    blink = alldata.subjectdata(sub).Physio.pupil_size.blink(my_time);
    blink(blink==0) = NaN;
    subplot(3,1,sub);
    plot(my_time,ps);
    % plot(my_time,blink+5,'r','LineWidth',2);
    ylim([min(ps),max(ps)]);grid on; hold off;
end

%50 fwd, 50 back
%% filter setting

filter = 'low'; %'low'

sr=1000;

sampling_duration = seconds(1/sr);

subjects = cell2mat({alldata.subjectdata.ID}');

% clear val;

% visualization of filters
filter_order = 4;

if filter == "low"
    filter_range = 4; %Hz
elseif filter == "bandpass"
    filter_range = [0.5 4];
end

% bandpass
[b,a] = butter(filter_order,filter_range/(sr/2), filter);

% Plot the frequency response
% freqz(b, a, [], sr);

%% Preprocessing steps
% (1) gap expansion
for p = 1:nsubjects
    disp(["...participant number "+ p + "..."])
    %left
    alldata.subjectdata(p).Physio.LeftPDil.data.gapExpand = alldata.subjectdata(p).Physio.pupil_size.ps;
    [alldata.subjectdata(p).Physio.LeftPDil.data.gapExpand, ...
        alldata.subjectdata(p).Physio.LeftPDil.valid_id] = expandGaps(alldata.subjectdata(p).Physio.LeftPDil.data.gapExpand,...
            sr);
end


% (2) pupil dilation speed filtering
for p = 1:nsubjects
    disp(["...participant number "+ p + "..."])
    %left
    RAW = alldata.subjectdata(p).Physio.LeftPDil.data.gapExpand;
    RAW = pupilSpeedfilter(RAW,sr);
    alldata.subjectdata(p).Physio.LeftPDil.data.speedFilter = RAW;
end

clear RAW

% remove 'island data'
clc;
for p = 1:nsubjects
    %left
    RAW = alldata.subjectdata(p).Physio.LeftPDil.data.speedFilter;
    %plot(RAW(1600:1700)); hold on;
    alldata.subjectdata(p).Physio.LeftPDil.valid_id = pupilIslandRemove(RAW,sr,alldata.subjectdata(p).Physio.LeftPDil.valid_id);
    RAW(alldata.subjectdata(p).Physio.LeftPDil.valid_id == 0) = NaN;

    %plot(RAW(1600:1700));
    alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove = RAW;
end

clear RAW

%% get null file id before interpolation
subject_to_keep = [];
subject_to_drop = [];
subject_to_drop_id = [];
for p = 1:nsubjects
    if sum(~isnan(alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove)) & p ~= 120
        subject_to_keep = [subject_to_keep,p];
    else
        subject_to_drop = [subject_to_drop;alldata.subjectdata(p).ID];
        subject_to_drop_id = [subject_to_drop_id,p];
    end
end

%% continue preprocessing
% (3) interpolation
% data has to be fully interpolated to be filtered.
% fully interpolate while noting where should be later excluded in long_blink_masked
% note the start and end; they are not supposed to be interpolated;
nsubjects = length(subject_to_keep);

for s = 1:nsubjects
    p = subject_to_keep(s);
    disp(["...participant number "+ p + "..."])
    %left
    RAW = alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove;

    [alldata.subjectdata(p).Physio.LeftPDil.data.interpolation, ...
        alldata.subjectdata(p).Physio.LeftPDil.data.missing_data] = interpolate_data(RAW,sr);
    alldata.subjectdata(p).Physio.LeftPDil.valid_id = long_blink_mask(alldata.subjectdata(p).Physio.LeftPDil.valid_id,...
        alldata.subjectdata(p).Physio.LeftPDil.data.missing_data);

    %after removal of 'long blinks', calculate valid sample proportion
    alldata.subjectdata(p).Physio.ValidSample = mean(alldata.subjectdata(p).Physio.LeftPDil.valid_id);
end

for s = 1:length(subject_to_drop_id)
    p = subject_to_drop_id(s);
    alldata.subjectdata(p).Physio.ValidSample = 0;
end

clear RAW;

% (4) filtering (has to be full continuous signal to be filtered; NaN
% doesn't work
for s = 1:nsubjects
    p = subject_to_keep(s);
    RAW = alldata.subjectdata(p).Physio.LeftPDil.data.interpolation;
    alldata.subjectdata(p).Physio.LeftPDil.data.filter = filtfilt(b,a,double(RAW));
end

% (5) exclude data based on long blink mask
for s = 1:nsubjects
    p = subject_to_keep(s);

    RAW= alldata.subjectdata(p).Physio.LeftPDil.data.filter;
    %remove interpolation at start and end of recording
    na_ind=alldata.subjectdata(p).Physio.LeftPDil.data.missing_data;
    if na_ind.start(1) == 1
        RAW(na_ind.start(1):na_ind.end(1)) = NaN;
    end
    last_na_seg_id = length(na_ind.end);
    if na_ind.end(length(na_ind.end)) == length(RAW)
        RAW(na_ind.start(last_na_seg_id):na_ind.end(last_na_seg_id)) = NaN;
    end
    alldata.subjectdata(p).Physio.LeftPDil.data.out = RAW;
end

%% get list of "bad" participant
val = [];

for i = 1:nsubjects
    val(i,1) = alldata.subjectdata(i).Physio.ValidSample;
end

alldata.bad_participants = subjects(val(:,1) < 0.5,:)

subject_to_drop = unique([subject_to_drop; alldata.bad_participants],'rows');

writematrix(subject_to_drop,"../subjects_list/subject_to_drop_pupil.txt")

%% examine abnormal values
for p = 1:nsubjects
    if sum(~isnan(alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove)) && ... %js b3
            ~ismember(string(alldata.subjectdata(p).ID),subject_to_drop)
        makefigure(30,10);
        plot(alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove,'lineWidth',2); hold on;
        plot(alldata.subjectdata(p).Physio.LeftPDil.data.out + 0.1,'lineWidth',2);
        ylim([min(alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove) max(alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove)]);
        title(subjects(p,:));
        xlim([1 20*sr]);
        print(['~/Desktop/VRMID-analysis/mid-pupil/figures/fmri3/examine_data/',subjects(p,:),'.tiff'],'-dtiff','-r100');
    end
end

close all;
%% some people's data needs further cleaning; there was some extremely small values

% for p = [1,29,30,31,32,44,86,89,90,91,92,118,119,121]
%     alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove(alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove <= 5) = NaN;
% %     if p == 1
% %         subplot(2,3,1);
% %     elseif p < 33
% %         subplot(2,3,p-27);
% %     else
% %         subplot(2,3,6);
% %     end
% %     plot(alldata.subjectdata(p).Physio.LeftPDil.data.speedFilter);
% end

%% plot correlation between max and min pupil size (should be highly correlated)
pupil_range = [];
for s = 1:nsubjects
    p = subject_to_keep(s);
    pupil_range(p,:) = [min(alldata.subjectdata(p).Physio.LeftPDil.data.out) max(alldata.subjectdata(p).Physio.LeftPDil.data.out)];
end

plot(pupil_range(:,1),pupil_range(:,2),'ko');
ylim([5 13]);
xlim([5 13]);

%% examplary figure
sr = sr_old;
cm = brewermap(6,'Set2'); %colour map
subj = 124;

RAWl = alldata.subjectdata(subj).Physio.pupil_size.ps;
blks = alldata.subjectdata(subj).Physio.pupil_size.blink;
blks(blks==0) = NaN;
blks = blks+min(RAWl);
vis_duration = [1 21]; % seconds

ind = (vis_duration(1)*sr_old):(1+vis_duration(2)*sr_old);
RAW_visl = RAWl(ind);
blks = blks(ind);
ind = (vis_duration(1)*sr):(1+vis_duration(2)*sr);
my_ylim = [min(RAW_visl) max(RAW_visl)];

RAW_vis2l = alldata.subjectdata(subj).Physio.LeftPDil.data.gapExpand(ind);

RAW_vis3l = alldata.subjectdata(subj).Physio.LeftPDil.data.speedFilter(ind);

RAW_vis4l = alldata.subjectdata(subj).Physio.LeftPDil.data.interpolation(ind);

RAW_vis5l = alldata.subjectdata(subj).Physio.LeftPDil.data.filter(ind);

RAW_vis6l = alldata.subjectdata(subj).Physio.LeftPDil.data.out(ind);

makefigure(18,18);

subplot(2,2,1);
plot(ind/sr,RAW_visl,'lineWidth',3); hold on;
plot(ind/sr,RAW_vis2l+0.1,'lineWidth',3); hold on;
plot(ind/sr,blks,'r','lineWidth',3); hold on;
legend({'raw','gapExpand','blinks'},'Location','southeast');
ax = gca;
colororder(ax,cm([1 2],:));
title("left eye");
ylim([my_ylim(1),my_ylim(2)+0.1]);
xlim(vis_duration);


subplot(2,2,2);
plot(ind/sr,RAW_vis2l,'lineWidth',3); hold on;
plot(ind/sr,RAW_vis3l+0.1,'lineWidth',3); hold on;
legend({'gapExpand','speedfilter'},'Location','southeast');
ax = gca;
colororder(ax,cm([2 3],:));
ylim([my_ylim(1),my_ylim(2)+0.1]);
xlim(vis_duration);


subplot(2,2,3);
% plot(ind/120,RAW_vis3,'k', 'lineWidth',3); hold on;
plot(ind/sr,RAW_vis4l,'lineWidth',3); hold on;
plot(ind/sr,RAW_vis5l+0.1,'lineWidth',3); hold on;
legend({'interpolation','filtering'},'Location','southeast');
ax = gca;
colororder(ax,cm([4 5],:));
ylim([my_ylim(1),my_ylim(2)+0.1]);
xlim(vis_duration);

subplot(2,2,4);
% plot(ind/120,RAW_vis3,'k', 'lineWidth',3); hold on;
plot(ind/sr,RAW_vis6l,'k', 'lineWidth',1); hold on;
legend({'used data left','used data right'},'Location','southeast');
ylim([my_ylim(1),my_ylim(2)+0.1]);
xlim(vis_duration);

print(['~/Desktop/VRMID-analysis/mid-pupil/figures/fmri3/example_pupil_sub',num2str(p),'.tiff'],'-dtiff','-r300');

%% (6) resample to 100 Hz
sr_old = 1000;
sr_new = 200;

for s = 1:nsubjects
    p = subject_to_keep(s);
    alldata.subjectdata(p).beh = downsample(alldata.subjectdata(p).beh,sr_old/sr_new);
    alldata.subjectdata(p).Physio.pupil_position = downsample(alldata.subjectdata(p).Physio.pupil_position,sr_old/sr_new);
    alldata.subjectdata(p).Physio.pupil_size = downsample(alldata.subjectdata(p).Physio.pupil_size,sr_old/sr_new);
    alldata.subjectdata(p).Physio.LeftPDil.valid_id = downsample(alldata.subjectdata(p).Physio.LeftPDil.valid_id,sr_old/sr_new);
    alldata.subjectdata(p).Physio.LeftPDil.valid_id = downsample(alldata.subjectdata(p).Physio.LeftPDil.valid_id,sr_old/sr_new);
    alldata.subjectdata(p).Physio.LeftPDil.data.out = downsample(alldata.subjectdata(p).Physio.LeftPDil.data.out,sr_old/sr_new);
end

alldata.sampling_rate = sr_new;

%% (7) reshape to R shape and export data
clc;
TdataOut = [];

for s = 1:nsubjects
    p = subject_to_keep(s);
    disp(["processing participant " + p]);
    if ~isempty(alldata.subjectdata(p).beh)

            temp_TdataOut = [table(repmat(alldata.subjectdata(p).ID,[height(alldata.subjectdata(p).beh),1]),'VariableNames',{'subject'}), alldata.subjectdata(p).beh, ...
                array2table([alldata.subjectdata(p).Physio.LeftPDil.data.out, alldata.subjectdata(p).Physio.pupil_size.blink, ...
                alldata.subjectdata(p).Physio.pupil_position.xp, alldata.subjectdata(p).Physio.pupil_position.yp],...
                'VariableNames',{'pupilDiameter','blinks','pupil_x','pupil_y'})];

        TdataOut = [TdataOut; temp_TdataOut];


    end
end

% this takes a while
writetable(TdataOut, './pupillometry.csv');
clc;
