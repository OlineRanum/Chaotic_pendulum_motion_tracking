%--------------------------------------------
% Videoeditor for cropping and cutting of frames
%                       Oline A. Ranum S2018
%--------------------------------------------


function [nStart, nStop, xpos, ypos] = VideoEditor(Vid_Names)

    nStart = [1]; nStop  = [1];
    warning('off');
    n = 10;                                         % Number of steps skipped in videocutter
 
    %------------------------------------------------------
    % Cut video
    %------------------------------------------------------
     for j=1:length(Vid_Names)
        vid =  VideoReader([Vid_Names{j}]);
        NumberOfFrames = vid.NumberOfFrames;
     figure(1);
     for t = drange(1:n:NumberOfFrames)
            imshow(imrotate(read(vid,t),270));
            title(['Frame:  ' num2str(t)]);
            drawnow;
     end
     close all;   

        FLAG = false;
        while ~FLAG
           str2 = input(sprintf(['\nCut video length of video ', num2str(j), ' Y/N:   ']),'s');
           if ismember(str2, {'Y','y','n','N'})
              FLAG = true;
           end
        end
        while (str2 == 'Y' || str2 == 'y') 
            FLAG = false;
            while ~FLAG
               start = input('Choose startindex:                 ');
               stop  = input('Choose stopindex:                  ');
               if isempty(start) start=2; end
               if isempty(stop)  stop = NumberOfFrames+1; end
               
               if (ismember(start, [1:NumberOfFrames]) && ismember(stop, [1:NumberOfFrames]) && start < stop)
                            FLAG = true;      
               else
                   fprintf('\nIndex out of range \n')
               end
            end
            
            figure(1);
            for m = drange(start:n:stop)
                imshow(imrotate(read(vid,m),270));
                title(['Framerate:  ' num2str(m)]);
                drawnow;
            end
            
            close all;
            nStart(j) = start; nStop(j)  = stop;
           
           FLAG = false; 
           while ~FLAG
                str2 = input('\nRecut video Y/N:                   ','s');
                if ismember(str2, {'Y','y','n','N'})
                    FLAG = true;
                end
            end
            
            
        end 
     end 
     
    %------------------------------------------------------
    % Crop video
    %------------------------------------------------------
    img = read(VideoReader([Vid_Names{1}]),nStart(1));
    figure(1);
    image(imrotate(img(:,:,:),90));
    title('Crop video frame'); xlabel('[Pixel]'); ylabel('[Pixel]');
    yLimits = get(gca,'YLim')-0.0005*1e3;  
    yTicks = yLimits(2)-get(gca,'YTick');  
    set(gca,'YTickLabel',num2str(yTicks.'));


    FLAG = false;    
     while ~FLAG
           str = input('\nCrop video frame Y/N:              ', 's');
           if ismember(str, {'Y','y','n','N'})
              FLAG = true;
           end
        end
    
    
    
    while (str == 'Y' || str == 'y')
            FLAG = false;
            while ~FLAG
               image(imrotate(img(:,:,:),90));
               title('Crop video frame'); xlabel('[Pixel]'); ylabel('[Pixel]');
               yLimits = get(gca,'YLim')-0.0005*1e3;  
               yTicks = yLimits(2)-get(gca,'YTick');  
               set(gca,'YTickLabel',num2str(yTicks.'));
               ypos = input('xrange (ex. 100:500) =             ');
               xpos = input('yrange (ex. 200:900) =             ');
               
               if isempty(ypos) 
                   1:vid.Width+1; end
               if isempty(xpos)
                   xpos = 1:vid.Height+1; end
               if (ypos(1) == 1 || ypos(1) == 0)
                   ypos(1) = 2; end 
               if (xpos(1) == 1 || xpos(1) == 0)
                   xpos(1) = 2;end
               if (isvector(ypos) == 0)
                   ypos = input('xrange must be of format x1:x2) =             '); end
               if (isvector(ypos) == 0)
                   xpos = input('yrange must be of format y1:y2) =             '); end
               W = vid.Width; H = vid.Height;
               if (xpos(1) > 1 && ypos(1) > 1 && xpos(end) <= W && ypos(end) <= H && xpos(1)< xpos(end) && ypos(1) < ypos(end))
                            FLAG = true;      
               else
                   fprintf('\nIndexError \n \n')
               end
            end
            figure(1);
            image(imrotate(img(ypos,xpos,:),90));
            title('Crop video frame'); xlabel('[Pixel]'); ylabel('[Pixel]');
            yLimits = get(gca,'YLim')-0.0005*1e3;  
            yTicks = yLimits(2)-get(gca,'YTick');  
            set(gca,'YTickLabel',num2str(yTicks.'));
           FLAG = false; 
           while ~FLAG
                str = input('Recrop video   Y/N:                ', 's');
                if ismember(str, {'Y','y','n','N'})
                    FLAG = true;
                end
           end
        close all;
    end 
     
    end