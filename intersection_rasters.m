function intersection_rasters(stimulus_coding, behavioral_readout, plot_marginals, use_hidden_feature, signal_intensity)
% This function runs the simulations used in figures 3 and S1 of the paper
% "Cracking the Neural Code for Sensory Perception by Combining Statistics,
% Intervention and Behavior", by S. Panzeri, C. D. Harvey, E. Piasini, P.
% E. Latham and T. Fellin.
%
% -------------------------------------------------------------------------
% Examples:
%
% To reproduce the figures in the paper, use the following settings.
%
% Fig 3A: intersection_rasters('gaussian', 'difference')
% Fig 3B: intersection_rasters('gaussian', 'r1')
% Fig 3C: intersection_rasters('gaussian', 'sum')
% Fig S1A3: intersection_rasters('elongated_gaussian', 'difference', false, true, 0.01) 
% Fig S1A4: intersection_rasters('elongated_gaussian', 'difference', false, true, 0.06) 
% Fig S1A5: intersection_rasters('elongated_gaussian', 'difference', false, true, 0.2) 
% Fig S1B: intersection_rasters('gaussian', 'r1', true)
% -------------------------------------------------------------------------
% Usage:
%
% stimulus_coding: this parameter controls the general mechanism by which
% the stimulus is encoded in the neural response (r1,r2). Possible values:
% "gaussian", "gaussian_r1", "elongated_gaussian", "uniform_r1r2",
% "uniform_r1", "uniform_r1_nonlinear", "uniform_r1r2_nonlinear". See below
% in the code for a description of each of these.
%
% behavioral_readout: this parameter controls the orientation of the
% decision boundary. Possible values: 'r1' (vertical), 'r2' (horizontal),
% 'sum' (main diagonal, choice uses r1+r2), 'difference (other diagonal,
% choice uses r1-r2).
%
% 'plot_marginals': boolean flag controlling whether to add marginals to
% the scatter plot. If this is true, information values as in Fig S1D are
% also computed.
%
% 'use_hidden_feature': boolean flag controlling whether to bias the choice
% with a hidden neural feature (r3) with similar stimulus tuning to r1+r2.
%
% 'signal_intensity': 'rho' parameter controlling stimulus discriminability
% as defined in the SI.

if nargin < 3
    plot_marginals = false;
end
if nargin < 4
    use_hidden_feature = false;
end
if nargin < 5
    signal_intensity = NaN;
end

npoints=1000000; % number of points generated for smooth marginal estimation
fraction_npoints_joint = 0.0001; % thinning fraction for clearer joint distribution visualisation
thinned_npoints = npoints * fraction_npoints_joint;

% define covariance matrix for gaussian stimulus coding
switch stimulus_coding
    case {'gaussian', 'gaussian_r1'}
        gaussian_covariance = [0.02, -0.005; -0.005, 0.02];
    case 'elongated_gaussian'
        s_plus = 0.18; % corresponds to sigma_+ in the paper
        s_minus = 0.07; % corresponds to sigma_- in the paper
        gaussian_covariance = [[s_plus^2+s_minus^2, s_plus^2-s_minus^2]; [s_plus^2-s_minus^2, s_plus^2+s_minus^2]]./4;
    otherwise
        gaussian_covariance = [];
end

if use_hidden_feature
    s_3 = 0.1; % standard deviation of hidden neural feature r3 (sigma_3 in the paper)
else
    s_3 = NaN;
end
hidden_feature_variance = s_3^2;

% color settings
color_stim_1 = 'g';
color_stim_2 = 'b';
color_stim_boundary = 'k';
ylorbr = brewermap(9,'Reds');
color_choice_1 = ylorbr(4,:);
color_choice_2 = ylorbr(7,:);
color_choice_boundary = 'r';

% other plotting parameters
ms = 10; % size of dots
linecircle = 2 ; % line thickness for open circles
linewidth_marginal = 3; % line tickness of the marginal probabilities
linewidth_plot = 4; % line tickness of the sensory and decision boundaries
density_x_range = 0:0.001:1; % values for which to estimate the marginal KDE
density_bandwidth = 0.005; % marginal KDE bandwidth
density_kernel = 'box'; % type of kernel to use for marginal density estimation
if use_hidden_feature
    ms = 7;
    linecircle = 1 ;
    linewidth_plot = 2;
end


% generate and plot data for joint pdf
[~, ~, ~, ~,...
    xdata1_left, xdata1_right, xdata2_left, xdata2_right,...
    ydata1_left, ydata1_right, ydata2_left, ydata2_right] = generate_data(thinned_npoints, stimulus_coding, behavioral_readout, signal_intensity, gaussian_covariance, hidden_feature_variance);


% plot with green and blue for s=1,2, filled circles for correct and open
% circles for wrong trials
if plot_marginals
    figure(1)
    clf
    set(gcf, 'Position', [0,0,550,500]);
    subplot('Position', [0.3, 0.10, 0.6, 0.6]);
else
    figure(2)
    clf
    set(gcf, 'Position', [0, 0, 422, 370]);
end
box on;
hold on;
axis([0 1 0 1]);
h1 = plot(xdata1_left, ydata1_left,'o');
set(h1,'MarkerEdgeColor',color_stim_1,'MarkerFaceColor',color_stim_1,'markersize',ms)
h4 = plot(xdata2_right, ydata2_right,'o');
set(h4,'MarkerEdgeColor',color_stim_2,'MarkerFaceColor',color_stim_2,'markersize',ms)
h2 = plot(xdata1_right, ydata1_right,'o');
set(h2,'MarkerEdgeColor',color_stim_1,'MarkerFaceColor','none','markersize',ms,'LineWidth',linecircle)
h3 = plot(xdata2_left, ydata2_left,'o');
set(h3,'MarkerEdgeColor',color_stim_2, 'MarkerFaceColor','none','markersize',ms,'LineWidth',linecircle)
% plot stimulus boundary line
switch stimulus_coding
    case {'gaussian', 'elongated_gaussian', 'uniform_r1r2'}
        plot((0:0.1:1),1-(0:0.1:1),'LineStyle','--','Color', color_stim_boundary,'linewidth',linewidth_plot);
    case {'gaussian_r1', 'uniform_r1'}
        plot([0.5, 0.5], [0, 1], 'LineStyle','--','Color', color_stim_boundary, 'linewidth',linewidth_plot);
end
% plot choice boundary line
switch behavioral_readout
    case 'r1'
        plot([0.5, 0.5], [0, 1], 'Color', color_choice_boundary,'LineStyle', '--', 'LineWidth', linewidth_plot);
    case 'sum'
        plot((0:0.1:1),1-(0:0.1:1) + 0.005,'Color', color_choice_boundary, 'LineStyle', '--', 'linewidth', linewidth_plot);
    case 'r2'
        plot([0, 1], [0.5, 0.5], 'Color', color_choice_boundary,'LineStyle', '--', 'LineWidth', linewidth_plot);
    case 'difference'
        plot((0:0.1:1),(0:0.1:1) + 0.005,'Color', color_choice_boundary, 'LineStyle', '--', 'linewidth', linewidth_plot);%SPX
end
set(gca,'XTick',[0 1]);
set(gca,'YTick',[0 1]);
set(gca,'fontsize',16,...
    'XTick', [], 'XTickLabel', {},...
    'YTick', [], 'YTickLabel', {});

label_color_1 = 'k';
label_color_2 = 'k';
if plot_marginals
    y_label_pos_1 = [1.05, 0.10];
    y_label_pos_2 = [1.05, 0.90];
    y_label_pos_3 = [1.10, 0.50];
else
    y_label_pos_1 = [-0.05, 0.10];
    y_label_pos_2 = [-0.05, 0.90];
    y_label_pos_3 = [-0.08, 0.50];
end

% Place the new labels
manual_ylabels(1) = text(y_label_pos_1(1), y_label_pos_1(2), 'Low', 'FontSize', 14,'color', label_color_1);
manual_ylabels(2) = text(y_label_pos_2(1), y_label_pos_2(2), 'High', 'FontSize', 14,'color', label_color_1);
manual_ylabels(3) = text(y_label_pos_3(1), y_label_pos_3(2), 'r_2', 'FontSize', 16,'color', label_color_1);
set(manual_ylabels,'Rotation',90,'HorizontalAlignment','center')
manual_xlabels(1) = text(0.10, -0.05, 'Low', 'FontSize', 14,'color', label_color_2);
manual_xlabels(2) = text(0.90, -0.05, 'High', 'FontSize', 14,'color', label_color_2);
manual_xlabels(3) = text(0.5, -0.1, 'r_1', 'FontSize', 16,'color', label_color_2);


set(manual_xlabels,'HorizontalAlignment','center')


if plot_marginals
    % generate and plot marginal probabilities
    [xdata1, xdata2, ydata1, ydata2,...
        xdata1_left, xdata1_right, xdata2_left, xdata2_right,...
        ydata1_left, ydata1_right, ydata2_left, ydata2_right,...
        faithful_idx_1, left_idx_1,...
        faithful_idx_2, left_idx_2] = generate_data(npoints, stimulus_coding, behavioral_readout, signal_intensity, gaussian_covariance, hidden_feature_variance);
    
    % compute intersection on real data and null ditribution (shuffled
    % data)
    [I_stim, I_choice, II, fII] = intersection_statistics(...
        faithful_idx_1, left_idx_1,...
        faithful_idx_2, left_idx_2, false);
    [I_stim_sh, I_choice_sh, II_sh, fII_sh] = intersection_statistics(...
        faithful_idx_1, left_idx_1,...
        faithful_idx_2, left_idx_2, true);
    fprintf('Intersection statistics: I_stim=%f, I_choice=%f, fII=%f, II=%f\n', I_stim, I_choice, fII, II);
    fprintf('Chance levels: I_stim=%f, I_choice=%f, fII=%f, II=%f\n', I_stim_sh, I_choice_sh, fII_sh, II_sh);
    
    
    % x axis marginal
    subplot('Position', [0.3, 0.76, 0.6, 0.08]); %SPX
    %stimulus panel
    hold on;
    [s1x_y, s1x_x] = ksdensity(xdata1, density_x_range, 'Kernel', density_kernel, 'Bandwidth', density_bandwidth);
    plot(s1x_x,s1x_y, 'Color', color_stim_1,'linewidth',linewidth_marginal, 'Clipping', 'off');
    [s2x_y, s2x_x] = ksdensity(xdata2, density_x_range, 'Kernel', density_kernel, 'Bandwidth', density_bandwidth);
    plot(s2x_x,s2x_y, 'Color', color_stim_2,'linewidth',linewidth_marginal, 'Clipping', 'off');
    axis([0 1 0 max(s1x_y)*1.01])
    set(gca, 'XTick', [], 'XTickLabel', '', 'YTick', [], 'YTickLabel', '');
    text(0.1, label_position(s1x_x, s1x_y, 0.1), 's=1', 'Color', color_stim_1, 'FontSize', 14, 'HorizontalAlignment', 'center');
    text(0.9, label_position(s2x_x, s2x_y, 0.9), 's=2', 'Color', color_stim_2, 'FontSize', 14, 'HorizontalAlignment', 'center');
    
    %choice panel
    subplot('Position', [0.3, 0.86, 0.6, 0.08]); %SPX
    hold on;
    [cleftx_y, cleftx_x] = ksdensity([xdata1_left xdata2_left], density_x_range, 'Kernel', density_kernel, 'Bandwidth', density_bandwidth);
    plot(cleftx_x,cleftx_y,'LineStyle', '-', 'Color', color_choice_1,'linewidth',linewidth_marginal, 'Clipping', 'off');
    [crightx_y, crightx_x] = ksdensity([xdata1_right xdata2_right], density_x_range, 'Kernel', density_kernel, 'Bandwidth', density_bandwidth);
    plot(crightx_x,crightx_y,'LineStyle', '-', 'Color', color_choice_2,'linewidth',linewidth_marginal, 'Clipping', 'off');
    axis([0 1 0 max(cleftx_y)*1.01])
    set(gca, 'XTick', [], 'XTickLabel', '', 'YTick', [], 'YTickLabel', '');
    text(0.1, label_position(cleftx_x, cleftx_y, 0.1), 'c=1', 'Color', color_choice_1, 'FontSize', 14, 'HorizontalAlignment', 'center');
    text(0.9, label_position(crightx_x, crightx_y, 0.9), 'c=2', 'Color', color_choice_2, 'FontSize', 14, 'HorizontalAlignment', 'center');
    
    % y axis marginal
    %stimulus panel
    subplot('Position', [0.17, 0.10, 0.08, 0.6])  %SPX
    hold on;
    [s1y_y, s1y_x] = ksdensity(ydata1, density_x_range, 'Kernel', density_kernel, 'Bandwidth', density_bandwidth);
    plot(s1y_y,s1y_x,'Color', color_stim_1,'linewidth',linewidth_marginal, 'Clipping', 'off');
    [s2y_y, s2y_x] = ksdensity(ydata2, density_x_range, 'Kernel', density_kernel, 'Bandwidth', density_bandwidth);
    plot(s2y_y,s2y_x,'Color', color_stim_2,'linewidth',linewidth_marginal, 'Clipping', 'off');
    axis([0 max(s1y_y)*1.01 0 1])
    set(gca, 'XTick', [], 'XTickLabel', '', 'YTick', [], 'YTickLabel', '', 'YAxisLocation', 'right', 'Xdir', 'reverse');
    text(label_position(s1y_x, s1y_y, 0.2), 0.2, 's=1', 'Color', color_stim_1, 'FontSize', 14, 'HorizontalAlignment', 'center','Rotation', 90);
    text(label_position(s2y_x, s2y_y, 0.8), 0.8, 's=2', 'Color', color_stim_2, 'FontSize', 14, 'HorizontalAlignment', 'center','Rotation', 90);
    
    
    %choice panel
    subplot('Position', [0.07, 0.10, 0.08, 0.6]) %SPX
    hold on;
    [clefty_y, clefty_x] = ksdensity([ydata1_left ydata2_left], density_x_range, 'Kernel', density_kernel, 'Bandwidth', density_bandwidth);
    plot(clefty_y,clefty_x,'LineStyle', '-', 'Color', color_choice_1,'linewidth',linewidth_marginal);
    [crighty_y, crighty_x] = ksdensity([ydata1_right ydata2_right], density_x_range, 'Kernel', density_kernel, 'Bandwidth', density_bandwidth);
    plot(crighty_y,crighty_x,'LineStyle', '-', 'Color', color_choice_2, 'linewidth',linewidth_marginal);
    axis([0 max(clefty_y)*1.01 0 1])
    set(gca, 'XTick', [], 'XTickLabel', '', 'YTick', [], 'YTickLabel', '', 'YAxisLocation', 'right', 'Xdir', 'reverse');
    text(label_position(clefty_x, clefty_y, 0.2), 0.2, 'c=1', 'Color', color_choice_1, 'FontSize', 14, 'HorizontalAlignment', 'center','Rotation', 90);
    text(label_position(crighty_x, crighty_y, 0.8), 0.8, 'c=2', 'Color', color_choice_2, 'FontSize', 14, 'HorizontalAlignment', 'center','Rotation', 90);

    figure(3)
    clf
    set(gcf, 'Position', [1129, 709, 250, 130]);
    bar([I_stim, II, II_sh], 'k');
    set(gca,...
        'YLim', [I_stim/2 - 0.02, 1],...
        'YTick', [0, 0.5, 1],...
        'XTickLabel', {'I(S)', 'II', 'II_chance'});
end

if use_hidden_feature && strcmp(stimulus_coding, 'elongated_gaussian')
    % plot neurometric and psychometric functions in a separate figure
    % (analytical form)
    figure(4)
    set(gcf, 'Position', [955, 891, 280, 100])
    clf
    s = 0:0.001:0.8;
    
    neurometric_plus = (1/2) * (1+erf(sqrt(2)*s/s_plus));
    % the following line would compute the neurometric function of either
    % r1 or r2 taken individually (they are identical)
    % neurometric_r1_or_r2 = (1/2) * (1+erf(s/sqrt((s_plus^2+s_minus^2)/2))); 
    psychometric = (1/2) * (1+erf(s/sqrt(2*(s_minus^2+s_3^2))));

    semilogx(s,neurometric_plus, 'linewidth', 2, 'color', 'k')
    hold on
    semilogx(s,psychometric, 'linestyle', '--', 'linewidth', 2, 'color', 'r')
    ylim([0.5, 1])
    ax = gca;
    set(ax, 'YTick', [0.5, 0.75, 1])
    ylim([0.48, 1.02])
    xlabel('Stimulus magnitude/coherence')
    ylabel('Fraction correct')
end

end

function y_pos = label_position(x_data, y_data, x_pos)
    % utility function for axis label positioning
    line_pos = y_data(x_data==x_pos);
    vertical = strcmp(get(gca, 'Xdir'), 'reverse');
    if vertical
        distance_from_line = 0.75;
    else
        distance_from_line = 0.5;
    end
    if line_pos < distance_from_line
        y_pos = line_pos + distance_from_line;
    else
        y_pos = line_pos - distance_from_line;
    end
end

function [xdata1, xdata2, ydata1, ydata2,...
    xdata1_left, xdata1_right, xdata2_left, xdata2_right,...
    ydata1_left, ydata1_right, ydata2_left, ydata2_right,...
    faithful_idx_1, left_idx_1, faithful_idx_2, left_idx_2] = ...
    generate_data(npoints, stimulus_coding, behavioral_readout,signal_intensity, gaussian_covariance, hidden_feature_variance)
    % Generate the data needed for the scatter plots. Keep track of the
    % stimulus, response and faithfulness of each trial.


    switch stimulus_coding
        case 'uniform_r1r2'
            % this is a raster of uncorrelated neuron activity with
            % information in both rate (r1) and time (r2).
            %
            % the probability satisfies the fact that the probablity to stimulus 1 is
            % one minus the probability to stimulus 2 for each point. In this way
            % P(r1r2) = p(r1)P(r2)
            % For stimulus 1 the joint response has probability alpha to be in the
            % lower half and 1-alpha to be in the upper half
            %for stimulus 2 the joint response has probability 1-alpha to be in the
            % lower half and alpha to be in the upper half
            
            alpha = 0.9; % probability of responses being in the correct triangle
            %(lower triangle for stimulus 1 and upper triangle for stimulus 2)
            
            
            NF = round(npoints*alpha);  % number of Faithful trials
            NM = round(npoints*(1-alpha));  % number of Misleading trials
            
            %generate Faithful trials to both stimuli
            xdata= rand(1,2*NF);
            ydata= rand(1,2*NF);
            sumdata=xdata+ydata;
            F1x=xdata(sumdata <=1);
            F1y=ydata(sumdata <=1);
            F2x=xdata(sumdata >1);
            F2y=ydata(sumdata >1);
            
            clear xdata ydata;
            
            %generate Misleading trials to both stimuli
            xdata= rand(1,2*NM);
            ydata= rand(1,2*NM);
            sumdata=xdata+ydata;
            M1x=xdata(sumdata >1);
            M1y=ydata(sumdata >1);
            M2x=xdata(sumdata <=1);
            M2y=ydata(sumdata <=1);
            
            xdata1 = [F1x M1x];
            ydata1 = [F1y M1y];
            xdata2 = [F2x M2x];
            ydata2 = [F2y M2y];
            
            
        case 'gaussian_r1'
            % this is a raster of correlated neuron activity (gaussian
            % distributions in r1 and r2) with information only in r1
            if isnan(signal_intensity)
                signal_intensity = 0.25;
            end
            
            %parameteters for the simulation
            mu1 = [0.5 - signal_intensity, 0.5];
            sigma1=gaussian_covariance;
            data1 = mvnrnd(mu1,sigma1,npoints);
            
            
            mu2 = [0.5 + signal_intensity, 0.5];
            sigma2=gaussian_covariance;
            data2 = mvnrnd(mu2,sigma2,npoints);
            
            
            xdata1= data1(:,1)';
            ydata1= data1(:,2)';
            F1x=xdata1(xdata1 <=0.5);
            F1y=ydata1(xdata1 <=0.5);
            M1x=xdata1(xdata1 >0.5);
            M1y=ydata1(xdata1 >0.5);
            
            xdata2= data2(:,1)';
            ydata2= data2(:,2)';
            F2x=xdata2(xdata2 >0.5);
            F2y=ydata2(xdata2 >0.5);
            M2x=xdata2(xdata2 <=0.5);
            M2y=ydata2(xdata2 <=0.5);
            
            % re-order xdata and ydata to make sure faithful trials are
            % always listed before misleading trials
            xdata1 = [F1x M1x];
            ydata1 = [F1y M1y];
            xdata2 = [F2x M2x];
            ydata2 = [F2y M2y];
            
            
        case 'uniform_r1_nonlinear'
            % this is a raster of uncorrelated neuron activity with
            % information in rate (r1) and with rate having behavior
            % information but encoded in a nonlinear way (s=2 is encoded in
            % a limited range of rates; rates higher or lower than that
            % encode s=1).
            %
            % the probability satisfies the fact that the probablity to stimulus 1 is
            % one minus the probability to stimulus 2 for each point. In this way
            % P(r1r2) = p(r1)P(r2)
            % For stimulus 1 the joint response has probability alpha to be outside the
            % central stripe and 1-alpha to be outside the central stripe
            
            alpha = 0.9; % probability of responses being in the correct part
            
            
            NF = round(npoints*alpha);  % number of Faithful trials
            NM= round(npoints*(1-alpha));  % number of Misleading trials
            
            %generate Faithful trials to both stimuli
            xdata= rand(1,2*NF);
            ydata= rand(1,2*NF);
            %rhodata=sqrt(((xdata-0.5).^2)+((ydata-0.5).^2));
            F1x=xdata((xdata <=0.25)|(xdata >0.75));
            F1y=ydata((xdata <=0.25)|(xdata >0.75));
            F2x=xdata((xdata >0.25)&(xdata <=0.75));
            F2y=ydata((xdata >0.25)&(xdata <=0.75));
            
            clear xdata ydata;
            
            %generate Misleading trials to both stimuli
            xdata= rand(1,2*NM);
            ydata= rand(1,2*NM);
            M2x=xdata((xdata <=0.25)|(xdata >0.75));
            M2y=ydata((xdata <=0.25)|(xdata >0.75));
            M1x=xdata((xdata >0.25)&(xdata <=0.75));
            M1y=ydata((xdata >0.25)&(xdata <=0.75));
            
            xdata1 = [F1x M1x];
            ydata1 = [F1y M1y];
            xdata2 = [F2x M2x];
            ydata2 = [F2y M2y];
            
            
        case 'uniform_r1r2_nonlinear'
            % this is a raster of uncorrelated neuron activity with
            % information in both rate (r1) and time (r2) and with both
            % codes having behavior information but encoded in a nonlinear
            % way
            %
            % the probability satisfies the fact that the probablity to stimulus 1 is
            % one minus the probability to stimulus 2 for each point. In this way
            % P(r1r2) = p(r1)P(r2)
            % For stimulus 1 the joint response has probability alpha to be inside the
            % circle and 1-alpha to be inside the circle
            
            alpha = 0.9; % probability of responses being in the correct part
            
            
            NF = round(npoints*alpha);  % number of Faithful trials
            NM= round(npoints*(1-alpha));  % number of Misleading trials
            
            %generate Faithful trials to both stimuli
            xdata= rand(1,2*NF);
            ydata= rand(1,2*NF);
            rhodata=sqrt(((xdata-0.5).^2)+((ydata-0.5).^2));
            F1x=xdata(rhodata <=0.25);
            F1y=ydata(rhodata <=0.25);
            F2x=xdata(rhodata >0.25);
            F2y=ydata(rhodata >0.25);
            
            clear xdata ydata;
            
            %generate Misleading trials to both stimuli
            xdata= rand(1,2*NM);
            ydata= rand(1,2*NM);
            rhodata=sqrt(((xdata-0.5).^2)+((ydata-0.5).^2));
            M1x=xdata(rhodata >0.25);
            M1y=ydata(rhodata >0.25);
            M2x=xdata(rhodata <=0.25);
            M2y=ydata(rhodata <=0.25);
            
            xdata1 = [F1x M1x];
            ydata1 = [F1y M1y];
            xdata2 = [F2x M2x];
            ydata2 = [F2y M2y];
            
        case 'uniform_r1'
            % this is a raster of uncorrelated neuron activity with information only in  rate
            %
            % the probability satisfies the fact that the probablity to stimulus 1 is
            % one minus the probaility to stimulus 2 for each point. In this way
            % P(r1r2) = p(r1)P(r2)
            % For stimulus 1 the joint response has probability alpha to be
            % above the diagonal and 1-alpha to be below the diagonal
            
            alpha = 0.9; % probability of responses being in the correct part
            
            
            NF = round(npoints*alpha);  % number of Faithful trials
            NM= round(npoints*(1-alpha));  % number of Misleading trials
            
            %generate Faithful trials to both stimuli
            xdata= rand(1,2*NF);
            ydata= rand(1,2*NF);
            %rhodata=sqrt(((xdata-0.5).^2)+((ydata-0.5).^2));
            F1x=xdata(xdata <=0.5);
            F1y=ydata(xdata <=0.5);
            F2x=xdata(xdata >0.5);
            F2y=ydata(xdata >0.5);
            
            clear xdata ydata;
            
            %generate Misleading trials to both stimuli
            xdata= rand(1,2*NM);
            ydata= rand(1,2*NM);
            M2x=xdata(xdata <=0.5);
            M2y=ydata(xdata <=0.5);
            M1x=xdata(xdata >0.5);
            M1y=ydata(xdata >0.5);
            
            xdata1 = [F1x M1x];
            ydata1 = [F1y M1y];
            xdata2 = [F2x M2x];
            ydata2 = [F2y M2y];
            
        case 'gaussian'
            % this is a raster of stimulus being encoded by gaussian
            % distributions in rate and time, with noise correlations given
            % by gaussian_covariance
            
            if isnan(signal_intensity)
                signal_intensity = 0.1;
            end
            %parameteters for the simulation
            mu1 = [0.5, 0.5] - signal_intensity;
            mu2 = [0.5, 0.5] + signal_intensity;
            data1 = mvnrnd(mu1,gaussian_covariance,npoints);
            data2 = mvnrnd(mu2,gaussian_covariance,npoints);
            
            xdata1= data1(:,1)';
            ydata1= data1(:,2)';
            sumdata1=xdata1+ydata1;
            F1x=xdata1(sumdata1 <=1);
            F1y=ydata1(sumdata1 <=1);
            M1x=xdata1(sumdata1 >1);
            M1y=ydata1(sumdata1 >1);
            
            xdata2= data2(:,1)';
            ydata2= data2(:,2)';
            sumdata2=xdata2+ydata2;
            F2x=xdata2(sumdata2 >1);
            F2y=ydata2(sumdata2 >1);
            M2x=xdata2(sumdata2 <=1);
            M2y=ydata2(sumdata2 <=1);
            
            % re-order xdata and ydata to make sure faithful trials are
            % always listed before misleading trials
            xdata1 = [F1x M1x];
            ydata1 = [F1y M1y];
            xdata2 = [F2x M2x];
            ydata2 = [F2y M2y];
            
        case 'elongated_gaussian'
            % this is a raster of stimulus being encoded by gaussian
            % distributions in r1 and r2, with strong noise correlations
            % between r1 and r2 ("elongated" gaussians) but no noise
            % correlations between (r1+r2) and (r1-r2) (so the principal
            % components of the stimulus-conditional distributions are
            % diagonal with respect to the r1,r2 axes).
            
            if isnan(signal_intensity)
                signal_intensity = 0.1;
            end
            
            r_mean_1 = [0.5, 0.5] - signal_intensity;
            r_mean_2 = [0.5, 0.5] + signal_intensity;
            
            r_1 = mvnrnd(r_mean_1, gaussian_covariance, npoints);% response to stimulus 1
            r_2 = mvnrnd(r_mean_2, gaussian_covariance, npoints);% response to stimulus 2
            
            f_1 = r_1(:,1) + r_1(:,2) <= 1; % faithfulness
            f_2 = r_2(:,1) + r_2(:,2) >= 1;
            
            F1x = r_1(cat(2, f_1, false(size(f_1))))';
            F1y = r_1(cat(2, false(size(f_1)), f_1))';
            M1x = r_1(cat(2, ~f_1, false(size(f_1))))';
            M1y = r_1(cat(2, false(size(f_1)), ~f_1))';
            F2x = r_2(cat(2, f_2, false(size(f_2))))';
            F2y = r_2(cat(2, false(size(f_2)), f_2))';
            M2x = r_2(cat(2, ~f_2, false(size(f_2))))';
            M2y = r_2(cat(2, false(size(f_2)), ~f_2))';

            xdata1 = [F1x M1x];
            ydata1 = [F1y M1y];
            xdata2 = [F2x M2x];
            ydata2 = [F2y M2y];
            
    end

    switch behavioral_readout
        case 'r1'
            %choice decision based only on x (left if x data less that 0.5)
            decision_variable_1 = xdata1 - 0.5;
            decision_variable_2 = xdata2 - 0.5;
            
        case 'sum'
            %choice decision based on the diagonal
            sumdata1= xdata1+ydata1;
            sumdata2= xdata2+ydata2;
            
            decision_variable_1 = sumdata1 - 1;
            decision_variable_2 = sumdata2 - 1;
        case 'r2'
            %choice decision based only on y (left if y data less than 0.5)
            decision_variable_1 = ydata1 - 0.5;
            decision_variable_2 = ydata2 - 0.5;
        case 'difference'
            %choice decision based on the orthogonal diagonal
            sumdata1= xdata1-ydata1;
            sumdata2= xdata2-ydata2;
            
            decision_variable_1 = sumdata1;
            decision_variable_2 = sumdata2;
    end
    
    if ~isnan(hidden_feature_variance)
        b_1 = normrnd(-signal_intensity, sqrt(hidden_feature_variance), size(decision_variable_1)); % "hidden" neural feature contributing to choice (r3)
        b_2 = normrnd(signal_intensity, sqrt(hidden_feature_variance), size(decision_variable_2));
        
        decision_variable_1 = decision_variable_1 + b_1;
        decision_variable_2 = decision_variable_2 + b_2;
    end
        
    left_idx_1 = decision_variable_1 < 0;
    left_idx_2 = decision_variable_2 < 0;
    
    xdata1_left= xdata1(left_idx_1);
    xdata1_right= xdata1(~left_idx_1);
    ydata1_left= ydata1(left_idx_1);
    ydata1_right= ydata1(~left_idx_1);
    xdata2_left= xdata2(left_idx_2);
    xdata2_right= xdata2(~left_idx_2);
    ydata2_left= ydata2(left_idx_2);
    ydata2_right= ydata2(~left_idx_2);
    
    n_faithful_1 = numel(F1x);
    n_faithful_2 = numel(F2x);
    faithful_idx_1 = ((1:numel(xdata1)) <= n_faithful_1);
    faithful_idx_2 = ((1:numel(xdata2)) <= n_faithful_2);
    

end
