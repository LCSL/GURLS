classdef ConsoleProgressBar < handle
    %ConsoleProgressBar Console progress bar for long-running operations
    %
    % Description:
    %   This class creates a console progress bar (status bar) 
    %   for long-running operations.
    %
    %   It is possible to measure the elapsed and remained time of progress.
    %
    %
    % Usage:
    %   cpb = ConsoleProgressBar();
    %
    %
    % Example:
    %   % Create Instance
    %   cpb = ConsoleProgressBar();
    %
    %   % Set progress bar parameters
    %   cpb.setLeftMargin(1);   % progress bar left margin
    %   cpb.setTopMargin(1);    % rows margin
    %   
    %   cpb.setLength(40);      % progress bar length: [.....]
    %   cpb.setMinimum(0);      % minimum value of progress range [min max]
    %   cpb.setMaximum(100);    % maximum value of progress range [min max]
    %   
    %   cpb.start();
    %   
    %   for k = 0:100
    %       userText = sprintf('Progress: [%d/%d]', k, 100);
    %       
    %       cpb.setValue(k);       	% update progress value
    %       cpb.setText(userText)   % update user text
    %       
    %       pause(0.025)
    %   end
    %   
    %   cpb.stop();
    %
    %
    % See also WAITBAR
    %
    
    % ---------------------------------------------------------------------
    %   Version   :  0.3
    %   Author    :  Evgeny Prilepin aka iroln <esp.home@gmail.com>
    %   Created   :  03.02.11
    %   Updated   :  05.02.11
    %
    %   Copyright :  (C) 2011 Evgeny Prilepin
    % ---------------------------------------------------------------------
    
    
    properties (GetAccess = public, SetAccess = private)
        % GetAccess public properties
        
        leftMargin = 1                  % Left margin in characters
        topMargin = 0                   % Top margin in rows
        progressLength = 50             % Length of progress bar in characters
        minimum = 0                     % Minimum progress value
        maximum = 100                   % Maximum progress value
        value = 0                       % Current progress value
        
        text = ''                       % User text
        isTextVisible = true            % User text visible flag
        textPosition = 'right'          % User text position 'left' or 'right'
        
        isPercentVisible = true         % Percent text visible flag
        percentPosition = 'left'        % Percent text position 'left' or 'right'
        
        isElapsedTimeVisible = false    % Elapsed time text visible flag
        elapsedTimePosition = 'right'   % Elapsed time text position
        
        isRemainedTimeVisible = false   % Remained time text visible flag
        remainedTimePosition = 'right'  % Remained time text position
        
        progressPercent = 0             % Percent of progress
        elapsedSeconds = 0              % Elapsed time from start in seconds
        remainedSeconds = Inf           % Remained time (prediction)
        
    end % GetAccess public properties
    
    
    properties (Access = private)
        % Private Properties
        
        isStarted = false
        
        progressBuffer = ''
        progressString = ''
        prevProgressStringLength = 0
        
        timeStart = 0
        timeStamp = []
        pctStamp = []
        
    end % Private Properties
    
    
    %======================================================================
    methods (Access = public)
        % Public Methods
        
        %------------------------------------------------------------------
        function obj = ConsoleProgressBar()
            %ConsoleProgressBar Constructor
            %
            % Description:
            %   Creates an instance of ConsoleProgressBar
            %
            % Usage:
            %   obj = ConsoleProgressBar()
            %
            
            error(nargchk(0, 0, nargin));
        end
        
        %------------------------------------------------------------------
        function start(obj)
            %start Start progress bar
            %
            % Description:
            %   Initialize and start new progress bar
            %
            % Usage:
            %   obj.start()
            %   start(obj)
            %
            
            error(nargchk(1, 1, nargin));
            
            obj.stop();
            
            obj.isStarted = true;
            obj.update();
            obj.resetTime();
        end
        
        %------------------------------------------------------------------
        function stop(obj)
            %start Stop progress bar
            %
            % Description:
            %   Stop the current progress bar
            %
            % Usage:
            %   obj.stop()
            %   stop(obj)
            %
            
            error(nargchk(1, 1, nargin));
            
            obj.progressString = '';
            obj.prevProgressStringLength = 0;
            
            obj.isStarted = false;
            obj.shiftToNewLine();
        end
        
        %------------------------------------------------------------------
        function reset(obj)
            %reset Reset progress bar
            %
            % Description:
            %   Reset progress bar to the minimum value.
            %
            % Usage:
            %   obj.reset()
            %   reset(obj)
            %
            
            error(nargchk(1, 1, nargin));
            
            obj.setValue(obj.minimum);
            obj.resetTime();
        end
        
        %------------------------------------------------------------------
        function setLength(obj, len)
            %setLength Set progress bar length
            %
            % Description:
            %   Sets the length of progress bar
            %
            % Usage:
            %   obj.setLength(len)
            %   setLength(obj, len)
            %
            % Inputs:
            %   len -- length in characters
            %
            
            error(nargchk(2, 2, nargin));
            
            validateattributes(len, {'numeric'}, ...
                {'scalar', '>=', 10, '<=', 500}, ...
                mfilename('fullpath'), 'Progress Bar Length');
            
            obj.progressLength = len;
            obj.update();
        end
        
        %------------------------------------------------------------------
        function setLeftMargin(obj, margin)
            %setLeftMargin Set progress bar left margin
            %
            % Description:
            %   Sets margin of the left edge of the command window
            %
            % Usage:
            %   obj.setLeftMargin(margin)
            %   setLeftMargin(obj, margin)
            %
            % Inputs:
            %   margin -- margin in characters
            %
            
            error(nargchk(2, 2, nargin));
            
            validateattributes(margin, {'numeric'}, ...
                {'scalar', '>=', 0, '<=', 100}, ...
                mfilename('fullpath'), 'Left Margin');
            
            obj.leftMargin = margin;
            obj.update();
        end
        
        %------------------------------------------------------------------
        function setTopMargin(obj, margin)
            %setTopMargin Set progress bar top margin
            %
            % Description:
            %   Sets the number of lines shift progress bar down
            %
            % Usage:
            %   obj.setTopMargin(margin)
            %   setTopMargin(obj, margin)
            %
            % Inputs:
            %   margin -- margin in rows
            %
            
            error(nargchk(2, 2, nargin));
            
            validateattributes(margin, {'numeric'}, ...
                {'scalar', '>=', 0, '<=', 10}, ...
                mfilename('fullpath'), 'Top Margin');
            
            obj.topMargin = margin;
        end
        
        %------------------------------------------------------------------
        function setValue(obj, val)
            %setValue Set current progress value
            %
            % Description:
            %   Sets current progress value
            %
            % Usage:
            %   obj.setValue(val)
            %   setValue(obj, val)
            %
            % Inputs:
            %   val -- progres value in range [min max]
            %
            
            error(nargchk(2, 2, nargin));
            
            validateattributes(val, {'numeric'}, ...
                {'scalar', '>=', obj.minimum, '<=', obj.maximum}, ...
                mfilename('fullpath'), 'Progress Value');
            
            obj.value = val;
            obj.update();
            obj.measureTime()
        end
        
        %------------------------------------------------------------------
        function setMinimum(obj, minVal)
            %setMinimum Set minimum progress range value
            %
            % Description:
            %   Sets minimum value of progress range
            %
            % Usage:
            %   obj.setMinimum(minVal)
            %   setMinimum(obj, minVal)
            %
            % Inputs:
            %   minVal -- minimum value of progress range [min max]
            %
            
            error(nargchk(2, 2, nargin));
            
            validateattributes(minVal, {'numeric'}, ...
                {'scalar', '<', obj.maximum}, ...
                mfilename('fullpath'), 'Minimum Value');
            
            obj.minimum = minVal;
            
            if (obj.value < obj.minimum)
                obj.value = obj.minimum;
            end
            
            obj.update();
        end
        
        %------------------------------------------------------------------
        function setMaximum(obj, maxVal)
            %setMaximum Set maximum progress range value
            %
            % Description:
            %   Sets maximum value of progress range
            %
            % Usage:
            %   obj.setMaximum(maxVal)
            %   setMaximum(obj, maxVal)
            %
            % Inputs:
            %   maxVal -- maximum value of progress range [min max]
            %
            
            error(nargchk(2, 2, nargin));
            
            validateattributes(maxVal, {'numeric'}, ...
                {'scalar', '>', obj.minimum}, ...
                mfilename('fullpath'), 'Maximum Value');
            
            obj.maximum = maxVal;
            
            if (obj.value > obj.maximum)
                obj.value = obj.maximum;
            end
            
            obj.update();
        end
        
        %------------------------------------------------------------------
        function setText(obj, userText)
            %setText Set user text
            %
            % Description:
            %   Sets any text, which will be displayed
            %
            % Usage:
            %   obj.setText(text)
            %   setText(obj, text)
            %
            % Inputs:
            %   text -- text string
            %
            
            error(nargchk(2, 2, nargin));
            
            validateattributes(userText, {'char'}, {'row'}, ...
                mfilename('fullpath'), 'Text');
            
            % Remove \n and \r characters
            userText = regexprep(userText, sprintf('(\n|\r)'), '');
            obj.text = userText;
            
            obj.update();
        end
        
        %------------------------------------------------------------------
        function setTextVisible(obj, flag)
            %setTextVisible Set user text visible
            %
            % Description:
            %   Sets visible flag of user text
            %
            % Usage:
            %   obj.setTextVisible(flag)
            %   setTextVisible(obj, flag)
            %
            % Inputs:
            %   flag -- true - visible, false - invisible
            %
            
            error(nargchk(2, 2, nargin));
            
            validateattributes(flag, {'numeric', 'logical'}, ...
                {'scalar', '>=', 0, '<=', 1}, ...
                mfilename('fullpath'), 'Text Visible');
            
            obj.isTextVisible = flag;
            
            obj.update();
            obj.showProgressBar();
        end
        
        %------------------------------------------------------------------
        function setTextPosition(obj, pos)
            %setTextPosition Set user text position
            %
            % Description:
            %   Sets position of user text
            %
            % Usage:
            %   obj.setTextPosition(pos)
            %   setTextPosition(obj, pos)
            %
            % Inputs:
            %   pos -- text position: 'left' or 'right'
            %
            
            error(nargchk(2, 2, nargin));
            
            pos = validatestring(lower(pos), {'left', 'right'}, ...
                mfilename('fullpath'), 'Text Position');
            
            obj.textPosition = pos;
            obj.update();
        end
        
        %------------------------------------------------------------------
        function setPercentVisible(obj, flag)
            %setPercentVisible Set percent text visible
            %
            % Description:
            %   Sets visible flag of percent text
            %
            % Usage:
            %   obj.setPercentVisible(flag)
            %   setPercentVisible(obj, flag)
            %
            % Inputs:
            %   flag -- true - text visible, false - text invisible
            %
            
            error(nargchk(2, 2, nargin));
            
            validateattributes(flag, {'numeric', 'logical'}, ...
                {'scalar', '>=', 0, '<=', 1}, ...
                mfilename('fullpath'), 'Percent Text Visible');
            
            obj.isPercentVisible = flag;
            obj.update();
        end
        
        %------------------------------------------------------------------
        function setPercentPosition(obj, pos)
            %setPercentPosition Set percent text position
            %
            % Description:
            %   Sets position of percent text
            %
            % Usage:
            %   obj.setPercentPosition(pos)
            %   setPercentPosition(obj, pos)
            %
            % Inputs:
            %   pos -- text position: 'left' or 'right'
            %
            
            error(nargchk(2, 2, nargin));
            
            pos = validatestring(lower(pos), {'left', 'right'}, ...
                mfilename('fullpath'), 'Percent Text Position');
            
            obj.percentPosition = pos;
            obj.update();
        end
        
        %------------------------------------------------------------------
        function setElapsedTimeVisible(obj, flag)
            %setElapsedTimeVisible Set elapsed time text visible
            %
            % Description:
            %   Sets visible flag of elapsed time text
            %
            % Usage:
            %   obj.setElapsedTimeVisible(flag)
            %   setElapsedTimeVisible(obj, flag)
            %
            % Inputs:
            %   flag -- true - text visible, false - text invisible
            %
            
            error(nargchk(2, 2, nargin));
            
            validateattributes(flag, {'numeric', 'logical'}, ...
                {'scalar', '>=', 0, '<=', 1}, ...
                mfilename('fullpath'), 'Elapsed Time Text Visible');
            
            obj.isElapsedTimeVisible = flag;
            obj.update();
        end
        
        %------------------------------------------------------------------
        function setElapsedTimePosition(obj, pos)
            %setElapsedTimePosition Set elapsed time text position
            %
            % Description:
            %   Sets position of elapsed time text
            %
            % Usage:
            %   obj.setElapsedTimePosition(pos)
            %   setElapsedTimePosition(obj, pos)
            %
            % Inputs:
            %   pos -- text position: 'left' or 'right'
            %
            
            error(nargchk(2, 2, nargin));
            
            pos = validatestring(lower(pos), {'left', 'right'}, ...
                mfilename('fullpath'), 'Elapsed Time Text Position');
            
            obj.elapsedTimePosition = pos;
            obj.update();
        end
        
        %------------------------------------------------------------------
        function setRemainedTimeVisible(obj, flag)
            %setRemainedTimeVisible Set elapsed time text visible
            %
            % Description:
            %   Sets visible flag of remained time text
            %
            % Usage:
            %   obj.setRemainedTimeVisible(flag)
            %   setRemainedTimeVisible(obj, flag)
            %
            % Inputs:
            %   flag -- true - text visible, false - text invisible
            %
            
            error(nargchk(2, 2, nargin));
            
            validateattributes(flag, {'numeric', 'logical'}, ...
                {'scalar', '>=', 0, '<=', 1}, ...
                mfilename('fullpath'), 'Elapsed Time Text Visible');
            
            obj.isRemainedTimeVisible = flag;
            obj.update();
        end
        
        %------------------------------------------------------------------
        function setRemainedTimePosition(obj, pos)
            %setRemainedTimePosition Set remained time text position
            %
            % Description:
            %   Sets position of remained time text
            %
            % Usage:
            %   obj.setRemainedTimePosition(pos)
            %   setRemainedTimePosition(obj, pos)
            %
            % Inputs:
            %   pos -- text position: 'left' or 'right'
            %
            
            error(nargchk(2, 2, nargin));
            
            pos = validatestring(lower(pos), {'left', 'right'}, ...
                mfilename('fullpath'), 'Remained Time Position');
            
            obj.remainedTimePosition = pos;
            obj.update();
        end
        
        %------------------------------------------------------------------
        function str = getElapsedTimeStr(obj, format)
            %getElapsedTimeStr Get elapsed time formated string
            %
            % Description:
            %   Gets formated string of elapsed time from start progress
            %
            % Note:
            %   This method is slow.
            %
            % Usage:
            %   obj.getElapsedTimeStr(format)
            %   getElapsedTimeStr(obj, format)
            %
            % Inputs:
            %   format -- DATESTR supported format
            %
            % Example:
            %   str = obj.getElapsedTimeStr('HH:MM:SS.FFF')
            %
            
            error(nargchk(2, 2, nargin));
            
            elapsedTimeNum = datenum([0 0 0 0 0 obj.elapsedSeconds]);
            str = datestr(elapsedTimeNum, format);
        end
        
        %------------------------------------------------------------------
        function str = getRemainedTimeStr(obj, format)
            %getRemainedTimeStr Get remained time formated string
            %
            % Description:
            %   Gets formated string of remained time from progress
            %
            % Usage:
            %   obj.getRemainedTimeStr(format)
            %   getRemainedTimeStr(obj, format)
            %
            % Inputs:
            %   format -- DATESTR supported format
            %
            % Example:
            %   str = obj.getRemainedTimeStr('HH:MM:SS.FFF')
            %
            
            error(nargchk(2, 2, nargin));
            
            if ~isinf(obj.remainedSeconds)
                remainedTimeNum = datenum([0 0 0 0 0 obj.remainedSeconds]);
                str = datestr(remainedTimeNum, format);
            else
                str = '--/--';
            end
        end
        
    end % Public Methods
    
    
    %======================================================================
    methods (Access = private)
        % Private Methods
        
        %------------------------------------------------------------------
        function update(obj)
            %update
            
            if obj.isStarted
                obj.calculatePercent();
                obj.fillProgressBuffer();
                obj.updateProgressString();
                obj.showProgressBar();
            end
        end
        
        %------------------------------------------------------------------
        function calculatePercent(obj)
            %calculatePercent
            
            obj.progressPercent = ...
                ((obj.value - obj.minimum) / (obj.maximum - obj.minimum)) * 100;
        end
        
        %------------------------------------------------------------------
        function fillProgressBuffer(obj)
            %fillProgressBuffer
            
            obj.progressBuffer = repmat('.', [1 obj.progressLength]);
            
            filledPart = round(...
                obj.progressLength * obj.progressPercent / 100);
            
            if (filledPart > 0)
                obj.progressBuffer(1:filledPart-1) = '=';
                obj.progressBuffer(filledPart) = '>';
            end
            
            obj.progressBuffer = ['[', obj.progressBuffer, ']'];
        end
        
        %------------------------------------------------------------------
        function updateProgressString(obj)
            %updateProgressString
            
            obj.prevProgressStringLength = length(obj.progressString);
            obj.progressString = obj.progressBuffer;
            
            obj.addPercentString();
            obj.addTimeString();
            obj.addTextString();
            
            obj.progressString = [blanks(obj.leftMargin), obj.progressString];
        end
        
        %------------------------------------------------------------------
        function addPercentString(obj)
            %addPercentString
            
            if obj.isPercentVisible
                percentText = sprintf('%d%%', fix(obj.progressPercent));
                marginBlanks = blanks(4 - length(percentText));
                
                switch obj.percentPosition
                    case 'left',  obj.progressString = [marginBlanks, percentText, ' ', obj.progressString];
                    case 'right', obj.progressString = [obj.progressString, ' ', marginBlanks, percentText];
                end
            end
        end
        
        %------------------------------------------------------------------
        function addTimeString(obj)
            %addTimeString
            
            if all(strcmpi('left', {obj.elapsedTimePosition, obj.remainedTimePosition}))
                % 
                cat_fun = @(e,r,es,rs) [e, es, r, rs, obj.progressString];
                
            elseif all(strcmpi('right', {obj.elapsedTimePosition, obj.remainedTimePosition}))
                % 
                cat_fun = @(e,r,es,rs) [obj.progressString, es, e, rs, r];
                
            elseif all(strcmpi({'left', 'right'}, {obj.elapsedTimePosition, obj.remainedTimePosition}))
                % 
                cat_fun = @(e,r,es,rs) [e, es, obj.progressString, rs, r];
            
            elseif all(strcmpi({'right', 'left'}, {obj.elapsedTimePosition, obj.remainedTimePosition}))
                % 
                cat_fun = @(e,r,es,rs) [r, rs, obj.progressString, es, e];
            
            end
            
            %FIXME: slow
            
            elapsedTimeStr = '';
            remainedTimeStr = '';
            es = '';
            rs = '';
            
            to_vec = @(s) datevec(datenum([0 0 0 0 0 round(s)]));
            to_str = @(h,m,s) sprintf('%02d:%02d:%02d', h, m, s);
            
            if obj.isElapsedTimeVisible
                dv = to_vec(obj.elapsedSeconds);
                elapsedTimeStr = ['E ', to_str(dv(4), dv(5), dv(6))];
                es = ' ';
            end
            
            if obj.isRemainedTimeVisible
                if ~isinf(obj.remainedSeconds)
                    dv = to_vec(obj.remainedSeconds);
                    remainedTimeStr = ['R ', to_str(dv(4), dv(5), dv(6))];
                else
                    remainedTimeStr = 'R --/--';
                end
                rs = ' ';
            end
            
            obj.progressString = cat_fun(elapsedTimeStr, remainedTimeStr, es, rs);
        end
        
        %------------------------------------------------------------------
        function addTextString(obj)
            %addTextString
            
            if obj.isTextVisible
                switch obj.textPosition
                    case 'left',  obj.progressString = [obj.text, ' ', obj.progressString];
                    case 'right', obj.progressString = [obj.progressString, ' ', obj.text];
                end
            end
        end
        
        %------------------------------------------------------------------
        function showProgressBar(obj)
            %showProgressBar
            
            backspaceStr = sprintf('\b');
            backspaceStr = backspaceStr(ones(1, obj.prevProgressStringLength));
            
            fprintf('%s%s', backspaceStr, obj.progressString);
        end
        
        %------------------------------------------------------------------
        function shiftToNewLine(obj)
            %switchNewLine
            
            newLines = sprintf('\n');
            newLines = newLines(ones(1, obj.topMargin));
            
            fprintf(newLines);
        end
        
        %------------------------------------------------------------------
        function resetTime(obj)
            %resetTime
            
            obj.timeStart = tic;
            
            obj.elapsedSeconds = 0;
            obj.remainedSeconds = Inf;
            obj.timeStamp = [];
            obj.pctStamp = [];
        end
        
        %------------------------------------------------------------------
        function measureTime(obj)
            %measureTime
            
            obj.elapsedSeconds = toc(obj.timeStart);
            
            % Accumulation time and percent values (timespent)
            obj.timeStamp = horzcat(obj.timeStamp, obj.elapsedSeconds);
            obj.pctStamp = horzcat(obj.pctStamp, obj.progressPercent);
            
            % Prediction remained time (simple linear least square estimate method)
            if (length(obj.timeStamp) > 1)
                % y = k*x + b, b == 0
                xy = sum(obj.timeStamp .* obj.pctStamp);
                x2 = sum(obj.pctStamp .* obj.pctStamp);
                predictionSeconds = (xy / x2) * 100;
                
                obj.remainedSeconds = abs(predictionSeconds - obj.elapsedSeconds);
            end
        end
        
    end % Private Methods
    
end % ConsoleProgressBar

