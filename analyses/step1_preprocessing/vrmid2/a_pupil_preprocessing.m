clear all; clc;

cd '~/Desktop/VRMID-analysis/mid-pupil/data/vrmid2/raw'
physio_data = 'physio_120hz.csv';
addpath(genpath('~/Desktop/VRMID-analysis/mid-pupil/analyses/step1_preprocessing/functions/functions-vrmid'));
T = readtable(physio_data);
%%
fmt = "hh:mm:ss.SSSSS";
sr = 120; %sampling rate in Hz
sampling_duration = seconds(1/sr);

% make structure
subjects = unique(T.Subject);
nsubjects = length(subjects);

check_columns = T(T.LeftPDil == 1 | T.RightPDil == 1 ,:);

check_columns2 = T(T.LeftPDil > 1 & T.LeftPDil <2 | T.RightPDil > 1 & T.RightPDil <2 ,:);

check_columns3 = T(T.LeftPDil > 1 & T.LeftPDil <= 0 | T.RightPDil > 1 & T.RightPDil <= 0 ,:);
%% some rightPDil which should be -1 was written as 1. this also affect AvgPDil
T.RightPDil(T.RightPDil == 1) = -1;
T.LeftPDil(T.LeftPDil == 1) = -1;
%
subplot(2,2,1);
ndhist(T.RightPPos_x(T.RightPPos_x ~= -1 & T.RightPPos_y ~= -1),...
    T.RightPPos_y(T.RightPPos_x ~= -1 & T.RightPPos_y ~= -1));

subplot(2,2,2);
ndhist(T.LeftPPos_x(T.LeftPPos_x ~= -1 & T.LeftPPos_y ~= -1),...
    T.LeftPPos_y(T.LeftPPos_x ~= -1 & T.LeftPPos_y ~= -1));

subplot(2,2,3:4);
histogram(T.RightPPos_x(T.RightPPos_x ~= -1 & T.RightPPos_y ~= -1)); hold on;
histogram(T.LeftPPos_x(T.LeftPPos_x ~= -1 & T.LeftPPos_y ~= -1));
legend('Right','Left');
%%
subplot(1,2,1);
histogram(T.LeftPDil(T.LeftPDil ~= -1));
subplot(1,2,2);
histogram(T.RightPDil(T.RightPDil ~= -1));

%%
subplot(1,2,1);
histogram(T.LeftPPos_x(T.LeftPPos_x ~= -1));
subplot(1,2,2);
histogram(T.RightPPos_x(T.RightPPos_x ~= -1));
%% convert file structure
clear alldata;

for s = 1:nsubjects
    sub = subjects(s);
    temp = T(string(cell2mat(T.Subject)) == string(sub{1}),:);
    alldata.subjectdata(s).ID = sub{1};
    times = temp.Time;
    pupil_position = temp(:,["LeftPPos_x","LeftPPos_y","LeftPPos_c",...
        "RightPPos_x","RightPPos_y","RightPPos_c"]);
    pupil_size = temp(:,["LeftOpen","LeftOpen_c","LeftPDil",...
        "LeftPDil_c","RightOpen","RightOpen_c","RightPDil","RightPDil_c"]);
    % go through every seconds
    [times_conv, pupil_position_conv, pupil_size_conv] = convert_time(times,...
        fmt,sr,pupil_position,pupil_size);
    alldata.sampling_rate = sr;
    alldata.subjectdata(s).TimePoints_orig = length(times);
    alldata.subjectdata(s).Times = times_conv;
    alldata.subjectdata(s).total_time = seconds(length(times_conv.time)/sr);
    alldata.subjectdata(s).Physio.pupil_position = pupil_position_conv;
    alldata.subjectdata(s).Physio.pupil_size = pupil_size_conv;
    alldata.subjectdata(s).Physio.ValidSample.L = sum(pupil_size_conv.LeftPDil ~= -1)/length(times);
    alldata.subjectdata(s).Physio.ValidSample.R = sum(pupil_size_conv.RightPDil ~= -1)/length(times);
%     [M I] = max(pupil_size_conv.LeftPDil_c) %weird...
end

%% plot subj's blinking to get a sense of the amount of gap expansion needed

my_time = [1,10];
my_time = (sr*my_time(1)):(sr*my_time(2));
for sub = 1:3
    ps = alldata.subjectdata(sub).Physio.pupil_size.LeftPDil(my_time);
    blink = alldata.subjectdata(sub).Physio.pupil_size.LeftOpen(my_time);
    blink(blink==0) = NaN;
    subplot(1,3,sub);
    plot(my_time,ps);
    % plot(my_time,blink+5,'r','LineWidth',2);
    ylim([min(ps),max(ps)]);grid on; hold off;
end


% get invalid participants
val = [];

for i = 1:nsubjects
    val(i,1) = alldata.subjectdata(i).Physio.ValidSample.L;
    val(i,2) = alldata.subjectdata(i).Physio.ValidSample.R;
end
alldata.bad_participants = subjects(val(:,1) < 0.5,:);
clear val;


%% set up filter
filter = 'low'; %'low'

sr = 120; %sampling rate in Hz
sampling_duration = seconds(1/sr);

% visualization of filters
filter_order = 3;

if filter == "low"
    filter_range = 4; %Hz
elseif filter == "bandpass"
    filter_range = [0.5 4];
end


% bandpass
[b,a] = butter(filter_order,filter_range/(sr/2), filter);

% Plot the frequency response
freqz(b, a, [], sr);
%% -1 into nan and gap expansion
for p = 1:nsubjects
    disp(["...participant number "+ p + "..."])
    %left
    alldata.subjectdata(p).Physio.LeftPDil.data.gapExpand = alldata.subjectdata(p).Physio.pupil_size.LeftPDil;
    alldata.subjectdata(p).Physio.LeftPDil.data.gapExpand(alldata.subjectdata(p).Physio.LeftPDil.data.gapExpand == -1) = NaN;
    [alldata.subjectdata(p).Physio.LeftPDil.data.gapExpand, ...
        alldata.subjectdata(p).Physio.LeftPDil.valid_id] = expandGaps(alldata.subjectdata(p).Physio.LeftPDil.data.gapExpand,...
            sr);
    %right
    alldata.subjectdata(p).Physio.RightPDil.data.gapExpand = alldata.subjectdata(p).Physio.pupil_size.RightPDil;
    alldata.subjectdata(p).Physio.RightPDil.data.gapExpand(alldata.subjectdata(p).Physio.RightPDil.data.gapExpand == -1) = NaN;
    [alldata.subjectdata(p).Physio.RightPDil.data.gapExpand, ...
        alldata.subjectdata(p).Physio.RightPDil.valid_id] = expandGaps(alldata.subjectdata(p).Physio.RightPDil.data.gapExpand,...
            sr);
end

%% pupil dilation speed filtering
clc;
for p = 1:nsubjects
    disp(["...participant number "+ p + "..."])
    %left
    RAW = alldata.subjectdata(p).Physio.LeftPDil.data.gapExpand;
    alldata.subjectdata(p).Physio.LeftPDil.valid_id = pupilSpeedfilter(RAW,sr,alldata.subjectdata(p).Physio.LeftPDil.valid_id);
    RAW(alldata.subjectdata(p).Physio.LeftPDil.valid_id == 0) = NaN;
    alldata.subjectdata(p).Physio.LeftPDil.data.speedFilter = RAW;


    %right
    RAW = alldata.subjectdata(p).Physio.RightPDil.data.gapExpand;
    alldata.subjectdata(p).Physio.RightPDil.valid_id = pupilSpeedfilter(RAW,sr,alldata.subjectdata(p).Physio.RightPDil.valid_id);
    RAW(alldata.subjectdata(p).Physio.RightPDil.valid_id == 0) = NaN;
    alldata.subjectdata(p).Physio.RightPDil.data.speedFilter = RAW;
end

clear RAW

%% remove 'island data'
clc;

for p = 1:nsubjects
    %left
    RAW = alldata.subjectdata(p).Physio.LeftPDil.data.speedFilter;
    %plot(RAW(1600:1700)); hold on;
    alldata.subjectdata(p).Physio.LeftPDil.valid_id = pupilIslandRemove(RAW,sr,alldata.subjectdata(p).Physio.LeftPDil.valid_id);
    RAW(alldata.subjectdata(p).Physio.LeftPDil.valid_id == 0) = NaN;

    %plot(RAW(1600:1700));
    alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove = RAW;


    %right
    RAW = alldata.subjectdata(p).Physio.RightPDil.data.speedFilter;
    alldata.subjectdata(p).Physio.RightPDil.valid_id = pupilIslandRemove(RAW,sr,alldata.subjectdata(p).Physio.RightPDil.valid_id);
    RAW(alldata.subjectdata(p).Physio.RightPDil.valid_id == 0) = NaN;
    alldata.subjectdata(p).Physio.RightPDil.data.islandRemove = RAW;
end

clear RAW

%% get null file id before interpolating
subject_to_keep = [];
subject_to_drop = [];
subject_to_drop_id = [];
for p = 1:nsubjects
    if prod(isnan(alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove)) || prod(isnan(alldata.subjectdata(p).Physio.RightPDil.data.islandRemove))
        subject_to_drop = [subject_to_drop;alldata.subjectdata(p).ID];
        subject_to_drop_id = [subject_to_drop_id,p];
    else
        subject_to_keep = [subject_to_keep,p];
    end
end

subjects = cell2mat({alldata.subjectdata.ID}');
subject_to_drop = unique([subject_to_drop; alldata.bad_participants],'rows');
writematrix(cell2mat(subject_to_drop),"../subjects_list/subject_to_drop_pupil.txt")

%% interpolation
% data has to be fully interpolated to be filtered.
% fully interpolate while noting where should be later excluded in long_blink_masked
% note the start and end; they are not supposed to be interpolated;

for p = 1:nsubjects
    disp(["...participant number "+ p + "..."])
    %left
    RAW = alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove;

    [alldata.subjectdata(p).Physio.LeftPDil.data.interpolation, ...
        alldata.subjectdata(p).Physio.LeftPDil.data.missing_data] = interpolate_data(RAW,sr);
    alldata.subjectdata(p).Physio.LeftPDil.valid_id = long_blink_mask(alldata.subjectdata(p).Physio.LeftPDil.valid_id,...
        alldata.subjectdata(p).Physio.LeftPDil.data.missing_data);

    %right
    RAW = alldata.subjectdata(p).Physio.RightPDil.data.islandRemove;

    [alldata.subjectdata(p).Physio.RightPDil.data.interpolation, ...
        alldata.subjectdata(p).Physio.RightPDil.data.missing_data] = interpolate_data(RAW,sr);
    alldata.subjectdata(p).Physio.RightPDil.valid_id = long_blink_mask(alldata.subjectdata(p).Physio.RightPDil.valid_id,...
        alldata.subjectdata(p).Physio.RightPDil.data.missing_data);

end

clear RAW;
% filtering
for p = 1:nsubjects
    RAW = alldata.subjectdata(p).Physio.LeftPDil.data.interpolation;
    alldata.subjectdata(p).Physio.LeftPDil.data.filter = filtfilt(b,a,double(RAW));

    RAW = alldata.subjectdata(p).Physio.RightPDil.data.interpolation;
    alldata.subjectdata(p).Physio.RightPDil.data.filter = filtfilt(b,a,double(RAW));
end

%% exclude data based on mask
for p = 1:nsubjects
    %left
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
%      alldata.subjectdata(p).Physio.LeftPDil.data.out(...
%          alldata.subjectdata(p).Physio.LeftPDil.valid_id == 0) = NaN;

    %right
    RAW= alldata.subjectdata(p).Physio.RightPDil.data.filter;
    %remove interpolation at start and end of recording
    na_ind=alldata.subjectdata(p).Physio.RightPDil.data.missing_data;
    if na_ind.start(1) == 1
        RAW(na_ind.start(1):na_ind.end(1)) = NaN;
    end
    last_na_seg_id = length(na_ind.end);
    if na_ind.end(length(na_ind.end)) == length(RAW)
        RAW(na_ind.start(last_na_seg_id):na_ind.end(last_na_seg_id)) = NaN;
    end
    alldata.subjectdata(p).Physio.RightPDil.data.out = RAW;
%      alldata.subjectdata(p).Physio.RightPDil.data.out(...
%          alldata.subjectdata(p).Physio.RightPDil.valid_id == 0) = NaN;
end

% for p = 1:nsubjects
%     p
%     [max(alldata.subjectdata(p).Physio.LeftPDil.data.out), max(alldata.subjectdata(p).Physio.RightPDil.data.out)]
% end

%% examplary figure
sr = 120;
cm = brewermap(8,'Set1'); %colour map

p = 15; %subject no
RAWl = alldata.subjectdata(p).Physio.pupil_size.LeftPDil;
RAWl(RAWl == -1) = NaN;

RAWr = alldata.subjectdata(p).Physio.pupil_size.RightPDil;
RAWr(RAWr == -1) = NaN;

vis_duration = [1 20]; % seconds
ind = (vis_duration(1)*sr):(1+vis_duration(2)*sr);

RAW_visl = RAWl(ind);
RAW_visr = RAWr(ind);

RAW_vis2l = alldata.subjectdata(p).Physio.LeftPDil.data.gapExpand(ind)+0.01;
RAW_vis2r = alldata.subjectdata(p).Physio.RightPDil.data.gapExpand(ind)+0.01;

RAW_vis3l = alldata.subjectdata(p).Physio.LeftPDil.data.speedFilter(ind);
RAW_vis3r = alldata.subjectdata(p).Physio.RightPDil.data.speedFilter(ind);

RAW_vis4l = alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove(ind)+0.01;
RAW_vis4r = alldata.subjectdata(p).Physio.RightPDil.data.islandRemove(ind)+0.01;

RAW_vis5l = alldata.subjectdata(p).Physio.LeftPDil.data.interpolation(ind);
RAW_vis5r = alldata.subjectdata(p).Physio.RightPDil.data.interpolation(ind);

RAW_vis6l = alldata.subjectdata(p).Physio.LeftPDil.data.filter(ind);
RAW_vis6r = alldata.subjectdata(p).Physio.RightPDil.data.filter(ind);

RAW_vis7l = alldata.subjectdata(p).Physio.LeftPDil.data.out(ind)+0.01;
RAW_vis7r = alldata.subjectdata(p).Physio.RightPDil.data.out(ind)+0.01;

makefigure(18,18);

subplot(5,2,1);
plot(ind/120,RAW_visl,'lineWidth',2); hold on;
plot(ind/120,RAW_vis2l+0.01,'lineWidth',2); hold on;
legend({'raw','gapExpand'},'Location','southeast');
ax = gca;
colororder(ax,cm([1 2],:));
title("left eye");
ylim([min(RAW_visl),max(RAW_visl)]);
xlim(vis_duration);

subplot(5,2,2);
plot(ind/120,RAW_visr,'lineWidth',2); hold on;
plot(ind/120,RAW_vis2r+0.01,'lineWidth',2); hold on;
legend({'raw','gapExpand'},'Location','southeast');
title("right eye");
ax = gca;
colororder(ax,cm([1 2],:));
ylim([min(RAW_visr),max(RAW_visr)]);
xlim(vis_duration);

subplot(5,2,3);
plot(ind/120,RAW_vis2l,'lineWidth',2); hold on;
plot(ind/120,RAW_vis3l+0.01,'lineWidth',2); hold on;
legend({'gapExpand','speedfilter'},'Location','southeast');
ax = gca;
colororder(ax,cm([2 3],:));
ylim([min(RAW_visl),max(RAW_visl)]);
xlim(vis_duration);

subplot(5,2,4);
plot(ind/120,RAW_vis2r,'lineWidth',2); hold on;
plot(ind/120,RAW_vis3r+0.01, 'lineWidth',2); hold on;
legend({'gapExpand','speedfilter'},'Location','southeast');
ax = gca;
colororder(ax,cm([2 3],:));
ylim([min(RAW_visr),max(RAW_visr)]);
xlim(vis_duration);

subplot(5,2,5);
% plot(ind/120,RAW_vis3,'k', 'lineWidth',3); hold on;
plot(ind/120,RAW_vis3l,'lineWidth',2); hold on;
plot(ind/120,RAW_vis4l+0.01,'lineWidth',2); hold on;
legend({'speedfilter','islandremove'},'Location','southeast');
ax = gca;
colororder(ax,cm([4 5],:));
ylim([min(RAW_visl),max(RAW_visl)]);
xlim(vis_duration);

subplot(5,2,6);
% plot(ind/120,RAW_vis3,'k', 'lineWidth',3); hold on;
plot(ind/120,RAW_vis3r,'lineWidth',2); hold on;
plot(ind/120,RAW_vis4r+0.01,'lineWidth',2); hold on;
legend({'speedfilter','islandremove'},'Location','southeast');
ax = gca;
colororder(ax,cm([4 5],:));
ylim([min(RAW_visr),max(RAW_visr)]);
xlim(vis_duration);

subplot(5,2,7);
% plot(ind/120,RAW_vis3,'k', 'lineWidth',3); hold on;
plot(ind/120,RAW_vis4l,'lineWidth',2); hold on;
plot(ind/120,RAW_vis5l+0.01,'lineWidth',2); hold on;
legend({'interpolation','filtering'},'Location','southeast');
ax = gca;
colororder(ax,cm([6 7],:));
ylim([min(RAW_visl),max(RAW_visl)]);
xlim(vis_duration);

subplot(5,2,8);
% plot(ind/120,RAW_vis3,'k', 'lineWidth',3); hold on;
plot(ind/120,RAW_vis4r,'lineWidth',2); hold on;
plot(ind/120,RAW_vis5r+0.01,'lineWidth',2); hold on;
legend({'interpolation','filtering'},'Location','southeast');
ax = gca;
colororder(ax,cm([6 7],:));
ylim([min(RAW_visr),max(RAW_visr)]);
xlim(vis_duration);

subplot(5,2,9);
% plot(ind/120,RAW_vis3,'k', 'lineWidth',3); hold on;
plot(ind/120,RAW_vis7l,'k', 'lineWidth',1); hold on;
plot(ind/120,RAW_vis7r+0.01,'k', 'lineWidth',2); hold on;
% legend({'used data left','used data right'},'Location','NorthEastOutside');
ylim([min([RAW_visl;RAW_visr]),max([RAW_visl;RAW_visr])]);
xlim(vis_duration);

print(['~/Desktop/VRMID-analysis/mid-pupil/figures/vrmid2/example_pupil_sub',num2str(p),'.tiff'],'-dtiff','-r300');
%% get avg pupil size for where both eyes were available
for p = 1:nsubjects
    disp(["...participant number "+ p + "..."])

    alldata.subjectdata(p).Physio.AvgPDil.data.out = (alldata.subjectdata(p).Physio.LeftPDil.data.out + ...
        alldata.subjectdata(p).Physio.RightPDil.data.out)/2;

    alldata.subjectdata(p).Physio.AvgPDil.data.out(isnan(alldata.subjectdata(p).Physio.LeftPDil.data.out) & ...
        isnan(alldata.subjectdata(p).Physio.RightPDil.data.out)) = NaN;

end

%% examine abnormal values
for p = 1:nsubjects
    if sum(~isnan(alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove)) %js b3
        makefigure(40,10);
        plot(alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove,'lineWidth',2); hold on;
        plot(alldata.subjectdata(p).Physio.LeftPDil.data.out + 0.1,'lineWidth',2);
        ylim([min(alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove) max(alldata.subjectdata(p).Physio.LeftPDil.data.islandRemove)]);
        xlim([1 20*sr]);
        title(subjects(p,:));
        print(['~/Desktop/VRMID-analysis/mid-pupil/figures/vrmid2/examine_data/',subjects(p,:),'.tiff'],'-dtiff','-r100');
    end
end

close all;

%% combine behavioral data
clc; %clear all;

% read in behavioral data
beh_data = 'per_second_data.csv';
T = readtable(beh_data);

% get only useful info
clear dataOut;

% %one subjects' id was written wrong
% T.Subject(string(T.Subject) == 'mr120612') = {'mr230612'};
%
clear dataOut;
dataOut.bad_participants = alldata.bad_participants;
p=1;
for nsub = 1:nsubjects
    if alldata.subjectdata(nsub).ID ~= string(cell2mat(dataOut.bad_participants))
        dataOut.subject(p).ID = alldata.subjectdata(nsub).ID;
        dataOut.subject(p).original_time_points = alldata.subjectdata(nsub).TimePoints_orig;
        dataOut.subject(p).total_time_in_sec = alldata.subjectdata(nsub).total_time(1);
        dataOut.subject(p).pupil_L = alldata.subjectdata(nsub).Physio.LeftPDil.data.out;
        dataOut.subject(p).pupil_R = alldata.subjectdata(nsub).Physio.RightPDil.data.out;
        dataOut.subject(p).pupil_Avg = alldata.subjectdata(nsub).Physio.AvgPDil.data.out;
        dataOut.subject(p).Times_ms = alldata.subjectdata(nsub).Times.timess;
        dataOut.subject(p).Times_s = alldata.subjectdata(nsub).Times.time;
        for i = 1:length(dataOut.subject(p).Times_s)
            if dataOut.subject(p).Times_s(i) < duration('07:00:00')
                dataOut.subject(p).Times_s(i) = dataOut.subject(p).Times_s(i) + duration('12:00:00');
            end
        end
        p=p+1;
    end
end

%
nsubjects = size(dataOut.subject(:),1);
% this takes a long time. ~15 min
for p = 1:nsubjects
    disp(["processing participant " + p]);
    clear beh_temp;
    beh_temp = T(string(T.Subject) == dataOut.subject(p).ID,:);
    beh_temp = sortrows(beh_temp,5); %sort by Time_str

    unique_times =unique(beh_temp.Time_str);

    for t = 1:length(unique_times) %needs to be Time not Time_str because of the
                                    % 1pm 13pm issue
        current_sec = unique_times(t);

        beh_this_sec = beh_temp(beh_temp.Time_str == current_sec,:);

        ind = find(dataOut.subject(p).Times_s == current_sec);
        dataOut.subject(p).behavior(ind,:) = repmat(beh_this_sec,[length(ind) 1]);
    end

    dataOut.subject(p).behavior = sortrows(dataOut.subject(p).behavior,5); %sort by Time_str
    physio_ind_retain = ismember(dataOut.subject(p).Times_s,beh_temp.Time_str);

    dataOut.subject(p).pupil_L =dataOut.subject(p).pupil_L(physio_ind_retain, :);
    dataOut.subject(p).pupil_R =dataOut.subject(p).pupil_R(physio_ind_retain, :);
    dataOut.subject(p).pupil_Avg =dataOut.subject(p).pupil_Avg(physio_ind_retain, :);
    dataOut.subject(p).Times_ms =dataOut.subject(p).Times_ms(physio_ind_retain, :);
    dataOut.subject(p).Times_s =dataOut.subject(p).Times_s(physio_ind_retain, :);

    if dataOut.subject(p).ID == 'mm230613' | dataOut.subject(p).ID == 'sw230608'
        beh_ind_retain = ismember(dataOut.subject(p).behavior.Time_str,unique(dataOut.subject(p).Times_s));
        dataOut.subject(p).behavior = dataOut.subject(p).behavior(beh_ind_retain,:);
    end
end

% reshape to R shape
clc;
TdataOut = [];
for p = 1:nsubjects
    disp(["processing participant " + p]);
    temp_TdataOut = [dataOut.subject(p).behavior, ...
        array2table([dataOut.subject(p).pupil_L, dataOut.subject(p).pupil_R,...
        dataOut.subject(p).pupil_Avg],...
        'VariableNames',{'pupil_L','pupil_R','pupil_Avg'})];
    TdataOut = [TdataOut; temp_TdataOut];
end

plot(temp_TdataOut.pupil_L(1:10000),temp_TdataOut.pupil_R(1:10000),'k.');

%
writetable(TdataOut, '../derivatives/pupillometry.csv');
clear TdataOut;
