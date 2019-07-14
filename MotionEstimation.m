%------------------------------------------------------
% RBG-based motionanalysis 
%                   Oline A. Ranum S2018
%------------------------------------------------------

%--------------------------------------------
% Parameters to be set by user
%--------------------------------------------

Vid_Names = {'vid1.mp4','vid2.mp4','vid3.mp4'};                 
                                                                 % Regular_values:
threshold_R        = 0.9;                                        % 0.9 for bottom pendulum  [red channel]   
threshold_G        = 0.85;                                       % 0.85 for upper pendulum  [green channel]
threshold_B        = 0.8;                                        % 0.8 for reference point  [blue channel]


%------------------------------------------------------
% Debugg and efficiency regulators
% On = true; Off = false
%------------------------------------------------------

warning('off');

global print_flag; print_flag = true;
global show_imgs; show_imgs = false;                            % Show path of pendulum (regulates efficiency of motion-estimation)
global debug_flag; debug_flag = false;                          % Allows manual handeling of noice
global PANIC_MODE; PANIC_MODE = false;

weights_R = [4 1 1];                                             % Colorweights for relative RGB-values
weights_G = [1 4 1];                                             % Can be adjusted to improve noice reduction in a specific channel

COM_distance_threshold = 200;                                    % Maximum traveling distance (in pixels) allowed between two frames.
color_distance_threshold = 50;                                   % Maximum color distance allowed between two frames.


%------------------------------------------------------
% Run Video Editor - Neccesary firs round
%------------------------------------------------------

FLAG = false;
while ~FLAG                                                     % Flag testing for valid userpoint
   str = input('Run Videoeditor Y/N:   ','s');
   if ismember(str, {'Y','y','n','N'})
      FLAG = true;
   end
end

if (str == 'Y' || str == 'y')
     [nStart, nStop, xpos, ypos] = VideoEditor(Vid_Names);
end

%------------------------------------------------------
% Creating reference point from blue channel
%------------------------------------------------------
A = read(VideoReader(Vid_Names{1}),nStart(1));
[ref_x, ref_y] = reference(A(ypos, xpos, :), threshold_B);

%------------------------------------------------------
% Estimation and plotting of motion
%------------------------------------------------------

for i=1:length(Vid_Names)
    %--------------------------------------------------
    % Analysing vid #i 
    %--------------------------------------------------
    vid  = VideoReader([Vid_Names{i}]);
    fpr  = vid.frameRate;  
    indx = nStart:nStop;
    
    [PR_x, PR_y, t_P] = motion_analysis('red', vid, indx, xpos, ypos, threshold_R, fpr, weights_R, COM_distance_threshold, color_distance_threshold);
    [PG_x, PG_y, t_G] = motion_analysis('green', vid, indx, xpos, ypos, threshold_G, fpr, weights_G, COM_distance_threshold, color_distance_threshold);

    [theta_R, theta_G] = vinkelanalyse(PR_x, PR_y, ref_x, ref_y, PG_x, PG_y);
    [PR_x, PR_y, t_R,PG_x, PG_y, t_G, theta_R, theta_G] = adapt(PR_x, PR_y, t_G,PR_x, PG_y, t_G, theta_G, theta_R);
    %------------------------------------------------------
    % Plotting
    %------------------------------------------------------
    legendInfo{i} = ['X = ' Vid_Names{i}];
    colormap winter;
    
    
    figure(11); hold on; 
    p1 = patch(PR_y,PR_x,t_R(1:length(PR_x)),'EdgeColor','interp','Marker','x','MarkerFaceColor','flat','FaceColor', 'none');
   % title('O2'); ylabel('X [pxs]'); xlabel('Y [pxs]');
    title('Motion of chaotic pendulum','fontsize',14);xlabel('Y [px]','fontsize',14); ylabel('X [px]','fontsize',14)
    c = colorbar('fontsize', 14);
    c.Label.String = 'Time [s]';
    ax = gca;
    ax.FontSize = 12; 
    
    figure(12); hold on;
    plot(t_R,theta_R,'-x')
    title('Angular motion of lower joint','fontsize',14);xlabel('Time [s]', 'fontsize',14); ylabel('\Theta [rad]','fontsize',14)
    omega_u = diff(theta_G)./diff(t_G);
    omega_d = diff(theta_R)./diff(t_G);
    ax = gca;
    ax.FontSize = 12; 
    
    figure(13); hold on;
    plot(t_G,theta_G,'-x');
    title('Angular motion of upper joint','fontsize',14);xlabel('Time [s]','fontsize',14); ylabel('\Theta [rad]','fontsize',14);
    ax = gca;
    ax.FontSize = 12; 
end
figure(12); legend(legendInfo);
figure(13); legend(legendInfo);
close(1)

%------------------------------------------------------
% Functions
%------------------------------------------------------
function col_plot(col_array)
    R=col_array(:,1); R_diff=abs(R(2:end)-R(1:end-1));
    G=col_array(:,2); G_diff=abs(G(2:end)-G(1:end-1));
    B=col_array(:,3); B_diff=abs(B(2:end)-B(1:end-1));
    figure(); plot(R, 'r');
    hold on;  plot(G, 'g');
    hold on;  plot(B, 'b');
    
    disp('HIGEHST  and  AVERAGE  one-frame differences:');
    disp(['Red:     ', num2str(max(R_diff)), '   ', num2str(mean(R_diff))]);
    disp(['Green:   ', num2str(max(G_diff)), '   ', num2str(mean(G_diff))]);
    disp(['Blue:    ', num2str(max(B_diff)), '   ', num2str(mean(B_diff))]);
end

function color_sample = get_colors(img, color_name)
    AreYouDoneYet = false;
    while ~AreYouDoneYet
        fig = figure(116);
        imshow(uint8(img));
        title(['Please click on the ', color_name, ' color. Press BACKSPACE to undo last selection. Press ENTER when done.']);
        [X, Y] = getpts(fig);
        if size(X,1) == 1
            AreYouDoneYet = true;
        else
            disp('Too many or few input points. Please Try again');
        end
    end
    close(116);
    x=round(X(1)); y=round(Y(1));
    color_sample=squeeze(img(y,x,:));
end

function [ref_x, ref_y] = reference(A, threshold)
    B = get_colors(A, 'blue');
    B = double(B);
    A2 = double(A);
    R_img = zeros(size(A2,1), size(A2,2));
    R_img(:,:) = sqrt((A2(:,:,1) - B(1)).^2 + (A2(:,:,2) - B(2)).^2 + (A2(:,:,3) - B(3)).^2 );
    R_img(:,:) = 1 - R_img/max(R_img(:));   
    R_bin = R_img > threshold;
    R_bin = bwareafilt(R_bin, 1);
    figure(1);
    title('Reference point in blue channel');
    imshow(R_bin);
    c = regionprops(R_bin,'centroid'); 
    ref_x = c.Centroid(1,2);
    ref_y = c.Centroid(1,1); 
    close(1)
end 

function [bpos_x, bpos_y, time] = motion_analysis(color_name, vid, indx, xpos, ypos, threshold, fpr, RGB_factors, COM_distance_threshold, color_distance_threshold)
        global debug_flag; global print_flag; global show_imgs; global PANIC_MODE;
        fargekanal = double(get_colors(read(vid,indx(1)-1), color_name));
        nBilder = indx(end)-indx(1)+1; 
        bpos_x = zeros(nBilder,1); 
        bpos_y = zeros(nBilder,1); 
        time = zeros(nBilder,1);
        col_array = zeros(size(indx,2)+1,3);
        col_array(1,:) = fargekanal.';
        for i=1:nBilder
            bilde = double(read(vid,indx(1)+i-1)); 
            bilde = bilde(ypos,xpos,:);
            R_img = zeros(size(bilde,1), size(bilde,2));
            R_img(:,:) = RGB_factors(1)*abs(bilde(:,:,1) - fargekanal(1)) + RGB_factors(2)*abs(bilde(:,:,2)-fargekanal(2)) + RGB_factors(3)*abs(bilde(:,:,3)-fargekanal(3));
            R_img(:,:) = 1 - R_img/(sum(RGB_factors)*255);
            R_bin = R_img > threshold;
            %R_bin = bwareaopen(R_bin, 500);
            R_bin = bwareafilt(R_bin, 1);

            c = regionprops(R_bin,'centroid');
            if size(c,1) == 0
                disp('No object in frame');
                bpos_x(i) = NaN;
                bpos_y(i) = NaN;
            elseif size(c,1) > 1
                if print_flag disp('Multiple objects in frame, using object 1'); end
                bpos_x(i) = c(1).Centroid(1,2);
                bpos_y(i) = c(1).Centroid(1,1);
            else
                bpos_x(i) = c.Centroid(1,2);
                bpos_y(i) = c.Centroid(1,1);
            end
            new_fargekanal = squeeze(bilde(round(bpos_x(i)), round(bpos_y(i)), :));
            color_diff = sum(abs(new_fargekanal(:) - fargekanal(:)));
            if i==1; distance=0; else; distance = sqrt((bpos_x(i) - bpos_x(i-1))^2 + (bpos_y(i) - bpos_y(i-1))^2);end

            if color_diff > color_distance_threshold || distance > COM_distance_threshold
                if print_flag || debug_flag
                    if color_diff>color_distance_threshold; disp(['New color of ', mat2str(new_fargekanal(:).') ,' has color difference of ', num2str(color_diff), ', exceeding limit.']); end
                    if distance>COM_distance_threshold; disp(['Distance at ', num2str(distance), ',which is above threshold']);end
                end
                if debug_flag
                    bilde2 = double(read(vid,indx(1)+i-1));
                    bilde2 = bilde2(ypos,xpos,:);
                    figure();imagesc(uint8(bilde2));title(["frame ", i]);
                    hold on; plot(bpos_y(i), bpos_x(i), 'r*')
                    figure();imshow(R_bin);title(["frame ", i]);
                    hold on; plot(bpos_y(i), bpos_x(i), 'r*')
                    pause();
                    close();close();
                end
                if PANIC_MODE
                    disp("Plz help. Resample color.")
                    fargekanal(:) = double(get_colors(bilde, color_name));
                end
            else
                fargekanal(:) = new_fargekanal(:);
                
            end
            
            if show_imgs
                figure(1);
                imshow(R_bin);
                hold on; plot(bpos_y(i), bpos_x(i), 'r*');
                title('Analysis of motion');
                drawnow; 
            end
            if print_flag; fprintf(1,'Frame: %3i    Mass Center(x,y): [%3i, %3i]    New Color: (%3i %3i %3i)    Color Dist: %3i    COM Dist: %3i\n',i, round(bpos_x(i)), round(bpos_y(i)), fargekanal(1), fargekanal(2), fargekanal(3), color_diff, round(distance)); end

            col_array(i+1,:) = fargekanal(:);
            
            time(i) = i/fpr;
        end
    bpos_x(end) = NaN; 
    bpos_y(end) = NaN;
    
    if debug_flag; col_plot(col_array); pause(); end
end

function [theta1, theta2] = vinkelanalyse(P_x, P_y, ref_x, ref_y, P2_x, P2_y)
    theta1 = zeros(length(P_x),1);
    theta2 = zeros(length(P2_x),1);
    
    for l=1:length(theta1)
        theta1(l) = asin((P_x(l)-ref_x)/sqrt((P_x(l)-ref_x)^2 + (P_y(l)-ref_y)^2));
        theta2(l) = asin((P2_x(l)-ref_x)/sqrt((P2_x(l)-ref_x)^2 + (P2_y(l)-ref_y)^2));
    end
    
    for k=2:length(theta1)
        if abs(theta1(k-1) - theta1(k)) > 150
            theta1(k) = - theta1(k);
        end
        if abs(theta2(k-1) - theta2(k)) > 150
            theta2(k) = - theta2(k);
        end
    end
end 

function [P_x, P_y, t_P,P2_x, P2_y, t_P2, theta1, theta2] = adapt(P_x, P_y, t_P,P2_x, P2_y, t_P2, theta1, theta2)
    traceback = 20;    
    t2 = theta2(1);
    if t2 < 0 
        for i=1:length(P_x)
                if theta2(i) > 0  
                    P_x = P_x(i-traceback:end); 
                    P_y = P_y(i-traceback:end);
                    t_P = t_P(1: length(P_y));
                    P2_x = P2_x(i-traceback:end); 
                    P2_y = P2_y(i-traceback:end);
                    t_P2 = t_P2(1:length(P_y));
                    theta1 = theta1(i-traceback:end);
                    theta2 = theta2(i-traceback:end);
                    break 
                end
                
        end
    else
               for i=1:length(P_x)
                if theta2(i) < 0  
                    P_x = P_x(i-traceback:end); 
                    P_y = P_y(i-traceback:end);
                    t_P = t_P(1: length(P_y));
                    P2_x = P2_x(i-traceback:end); 
                    P2_y = P2_y(i-traceback:end);
                    t_P2 = t_P2(1:length(P_y));
                    theta1 = theta1(i-traceback:end);
                    theta2 = theta2(i-traceback:end);
                    break 
                end
               end      
        end
    P_x(end) = NaN; P_y(end) = NaN; P2_x(end) = NaN; P2_y(end) = NaN; t_P(end) = NaN; t_P2(end) = NaN; theta1(end) = NaN; theta2(end) = NaN;
end 