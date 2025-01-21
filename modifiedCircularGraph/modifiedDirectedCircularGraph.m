classdef modifiedDirectedCircularGraph < handle
% CIRCULARGRAPH Plot an interactive circular graph to illustrate connections in a network.
%
%% Syntax
% circularGraph(X)
% circularGraph(X,'PropertyName',propertyvalue,...)
% h = circularGraph(...)
%
%% Description
% A 'circular graph' is a visualization of a network of nodes and their
% connections. The nodes are laid out along a circle, and the connections
% are drawn within the circle. Click on a node to make the connections that
% emanate from it more visible or less visible. Click on the 'Show All'
% button to make all nodes and their connections visible. Click on the
% 'Hide All' button to make all nodes and their connections less visible.
%
% Required input arguments.
% X : A symmetric matrix of numeric or logical values.
%
% Optional properties.
% Colormap : A N by 3 matrix of [r g b] triples, where N is the 
%            length(adjacenyMatrix).
% Label    : A cell array of N strings.
%%
% Copyright 2016 The MathWorks, Inc.
  properties
    Node = node(0,0); % Array of nodes
    ColorMap;         % Colormap
    Label;            % Cell array of strings
    % ShowButton;       % Turn all nodes on
    % HideButton;       % Turn all nodes off
  end
  
  methods
    function drawConnections(this, adjacencyMatrix, t, colorMap, coeff)
        % Find non-zero values of adjacencyMatrix and their indices
        [row, col, v] = find(adjacencyMatrix);

        % Calculate line widths based on values of adjacencyMatrix
        minLineWidth = 0.5;
        lineWidthCoef = 5;
        lineWidth = v ./ max(v);
        if sum(lineWidth) == numel(lineWidth) % all lines are the same width
            lineWidth = repmat(minLineWidth, numel(lineWidth), 1);
        else % lines of variable width
            lineWidth = lineWidthCoef * lineWidth + minLineWidth;
        end

        % Define arrowhead size
        % Arrowhead parameters
        arrowSize = 0.07; % Adjust this size for the arrowheads
        arrowBaseWidth = 0.75 * arrowSize; % Width of the arrowhead base

        % Draw connections
        for i = 1:length(v)
            if row(i) ~= col(i)
                u = [cos(t(row(i))); sin(t(row(i)))]; % from
                v = [cos(t(col(i))); sin(t(col(i)))]; % to
            % try
                if abs(row(i) - col(i)) - length(adjacencyMatrix) / 2 == 0
                    % Straight line for diametric points
                    this.Node(row(i)).Connection(end + 1) = line(...
                        [u(1); v(1)],...
                        [u(2); v(2)],...
                        'LineWidth', lineWidth(i),...
                        'Color', colorMap(round(coeff * length(adjacencyMatrix)/3), :),...
                        'PickableParts', 'none');
                else
                    % Arc for other connections
                    u  = [cos(t(row(i)));sin(t(row(i)))];
                    v  = [cos(t(col(i)));sin(t(col(i)))];
                    x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
                    y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
                    r  = sqrt(x0^2 + y0^2 - 1);
                    thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
                    thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
                    
                    if u(1) >= 0 && v(1) >= 0 
                      % ensure the arc is within the unit disk
                      theta = [linspace(max(thetaLim),pi,50),...
                               linspace(-pi,min(thetaLim),50)].';
                    else
                      theta = linspace(thetaLim(1),thetaLim(2)).';
                    end
                    % checker1 = r * cos(theta) + x0;
                    % checker2 = r * sin(theta) + y0;
                    xArc = r * cos(theta) + x0;
                    yArc = r * sin(theta) + y0;
                    this.Node(row(i)).Connection(end + 1) = line(...
                        r * cos(theta) + x0,...
                        r * sin(theta) + y0,...
                        'LineWidth', lineWidth(i),...
                        'Color', colorMap(round(coeff * length(adjacencyMatrix)/3), :),...
                        'PickableParts', 'none');
                   
                    % % Add arrowhead for the arc
                    % direction = [xArc(end) - xArc(end-1); yArc(end) - yArc(end-1)];
                    % direction = direction / norm(direction); % Normalize
                    % perp = [-direction(2); direction(1)]; % Perpendicular vector
                    % arrowTip = [xArc(end); yArc(end)]; % Arrowhead tip
                    % arrowBase1 = arrowTip - arrowSize * direction + arrowBaseWidth * perp;
                    % arrowBase2 = arrowTip - arrowSize * direction - arrowBaseWidth * perp;
                    % 
                    % patch(...
                    %     [arrowTip(1), arrowBase1(1), arrowBase2(1)], ...
                    %     [arrowTip(2), arrowBase1(2), arrowBase2(2)], ...
                    %     colorMap(round(coeff * length(adjacencyMatrix)/3), :), ...
                    %     'EdgeColor', 'none');
                    
                    % Draw arrowhead at the target node
                    % Calculate the angle of the connection for arrow orientation
                    angle = atan2(v(2), v(1)); 
                    arrowSize = 0.07; % Adjust arrow size
                    arrowWidth = 0.05; % Adjust arrow width
                    
                    % Arrowhead triangle coordinates
                    arrowX = [...
                        v(1), ...
                        v(1) - arrowSize * cos(angle) + arrowWidth * sin(angle), ...
                        v(1) - arrowSize * cos(angle) - arrowWidth * sin(angle)];
                    arrowY = [...
                        v(2), ...
                        v(2) - arrowSize * sin(angle) - arrowWidth * cos(angle), ...
                        v(2) - arrowSize * sin(angle) + arrowWidth * cos(angle)];
                    
                    % Create the arrowhead patch
                    patch(arrowX, arrowY, colorMap(round(coeff * length(adjacencyMatrix) / 3), :), ...
                        'EdgeColor', 'none', ...
                        'PickableParts', 'none');
                    
                end
                % % Add an arrowhead to indicate direction
                % arrowBase = v; % Endpoint of the line (arrow base)
                % arrowDirection = [arrowBase(1) - u(1), arrowBase(2) - u(2)];
                % arrowDirection = arrowDirection / norm(arrowDirection); % Normalize
                % 
                % % Define the arrowhead vertices
                % perpendicular = [-arrowDirection(2), arrowDirection(1)]; % Perpendicular vector
                % arrowPoint1 = arrowBase - arrowSize * arrowDirection + 0.5 * arrowSize * perpendicular;
                % arrowPoint2 = arrowBase - arrowSize * arrowDirection - 0.5 * arrowSize * perpendicular;
                % 
                % % Plot the arrowhead using a patch
                % patch(...
                %     'XData', [arrowBase(1), arrowPoint1(1), arrowPoint2(1)],...
                %     'YData', [arrowBase(2), arrowPoint1(2), arrowPoint2(2)],...
                %     'FaceColor', colorMap(length(adjacencyMatrix), :),...
                %     'EdgeColor', 'none');
                % 
                % scatter([arrowTip(1), arrowBase1(1), arrowBase2(1)], ...
                %     [arrowTip(2), arrowBase1(2), arrowBase2(2)], 'o');

            % catch
            %     warning('Skipping arc (%d, %d): invalid radius.', row(i), col(i));
            % end
            end
        end
    end
end

  
  methods
    function this = modifiedDirectedCircularGraph(adjacencyMatrix,varargin)
      % Constructor
      p = inputParser;
      
      defaultColorMap = parula(length(adjacencyMatrix));
      defaultLabel = cell(length(adjacencyMatrix));
      for i = 1:length(defaultLabel)
        defaultLabel{i} = num2str(i);
      end
      
      addRequired(p,'adjacencyMatrix',@(x)(isnumeric(x) || islogical(x)));
      addParameter(p, 'SecondAdjacencyMatrix', [], @(x) isempty(x) || (isnumeric(x) && isequal(size(x), size(adjacencyMatrix))));
      addParameter(p,'ColorMap',defaultColorMap,@(colormap)length(colormap) == length(adjacencyMatrix));
      addParameter(p,'Label'   ,defaultLabel   ,@iscell);
      
      parse(p,adjacencyMatrix,varargin{:});
      secondAdjacencyMatrix = p.Results.SecondAdjacencyMatrix;
      this.ColorMap = p.Results.ColorMap;
      this.Label    = p.Results.Label;
      
      % this.ShowButton = uicontrol(...
      %   'Style','pushbutton',...
      %   'Position',[0 40 80 40],...
      %   'String','Show All',...
      %   'Callback',@circularGraph.showNodes,...
      %   'UserData',this);
      % 
      % this.HideButton = uicontrol(...
      %   'Style','pushbutton',...
      %   'Position',[0 0 80 40],...
      %   'String','Hide All',...
      %   'Callback',@circularGraph.hideNodes,...
      %   'UserData',this);
      
      fig = gcf;
      set(fig,...
        'UserData',this,...
        'CloseRequestFcn',@modifiedDirectedCircularGraph.CloseRequestFcn);
      
      % Draw the nodes
      delete(this.Node);
      t = linspace(-pi,pi,length(adjacencyMatrix) + 1).'; % theta for each node
      extent = zeros(length(adjacencyMatrix),1);
      for i = 1:length(adjacencyMatrix)
        this.Node(i) = node(cos(t(i)),sin(t(i)));
        this.Node(i).Color = this.ColorMap(length(adjacencyMatrix),:); % set node color to be the fullest of either red or blue colormap 
        this.Node(i).Label = this.Label{i};
      end
      
      % % Find non-zero values of s and their indices
      % [row,col,v] = find(adjacencyMatrix);
      % 
      % % Calculate line widths based on values of s (stored in v).
      % minLineWidth  = 0.5;
      % lineWidthCoef = 5;
      % lineWidth = v./max(v);
      % if sum(lineWidth) == numel(lineWidth) % all lines are the same width.
      %   lineWidth = repmat(minLineWidth,numel(lineWidth),1);
      % else % lines of variable width.
      %   lineWidth = lineWidthCoef*lineWidth + minLineWidth;
      % end
      % 
      % % Draw connections on the Poincare hyperbolic disk.
      % %
      % % Equation of the circles on the disk:
      % % x^2 + y^2 
      % % + 2*(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1))*x 
      % % - 2*(u(1)-v(1))/(u(1)*v(2)-u(2)*v(1))*y + 1 = 0,
      % % where u and v are points on the boundary.
      % %
      % % Standard form of equation of a circle
      % % (x - x0)^2 + (y - y0)^2 = r^2
      % %
      % % Therefore we can identify
      % % x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
      % % y0 = (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
      % % r^2 = x0^2 + y0^2 - 1
      % 
      % for i = 1:length(v)
      %   if row(i) ~= col(i)
      %     if abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0 
      %       % points are diametric, so draw a straight line
      %       u = [cos(t(row(i)));sin(t(row(i)))];
      %       v = [cos(t(col(i)));sin(t(col(i)))];
      %       this.Node(row(i)).Connection(end+1) = line(...
      %         [u(1);v(1)],...
      %         [u(2);v(2)],...
      %         'LineWidth', lineWidth(i),...
      %         'Color', this.ColorMap(floor(length(v)/3),:),... % set color of 0.01 significance connections to be pale
      %         'PickableParts','none');
      %     else % points are not diametric, so draw an arc
      %       u  = [cos(t(row(i)));sin(t(row(i)))];
      %       v  = [cos(t(col(i)));sin(t(col(i)))];
      %       x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
      %       y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
      %       r  = sqrt(x0^2 + y0^2 - 1);
      %       thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
      %       thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
      % 
      %       if u(1) >= 0 && v(1) >= 0 
      %         % ensure the arc is within the unit disk
      %         theta = [linspace(max(thetaLim),pi,50),...
      %                  linspace(-pi,min(thetaLim),50)].';
      %       else
      %         theta = linspace(thetaLim(1),thetaLim(2)).';
      %       end
      % 
      %       this.Node(row(i)).Connection(end+1) = line(...
      %         r*cos(theta)+x0,...
      %         r*sin(theta)+y0,...
      %         'LineWidth', lineWidth(i),...
      %         'Color', this.ColorMap(floor(length(v)/3),:),... % set color of 0.01 significance connections to be pale
      %         'PickableParts','none');
      %     end
      %   end
      % end

      % Draw connections from the first adjacency matrix
      this.drawConnections(adjacencyMatrix, t, this.ColorMap, 1);

      % Draw connections from the second adjacency matrix (if provided)
      if ~isempty(secondAdjacencyMatrix)
          this.drawConnections(secondAdjacencyMatrix, t, this.ColorMap, 2);
      end
      
      axis image;
      ax = gca;
      for i = 1:length(adjacencyMatrix)
        extent(i) = this.Node(i).Extent;
      end
      extent = max(extent(:));
      ax.XLim = ax.XLim + extent*[-1 1];
      fudgeFactor = 1.75; % Not sure why this is necessary. Eyeballed it.
      ax.YLim = ax.YLim + fudgeFactor*extent*[-1 1];
      ax.Visible = 'off';
      ax.SortMethod = 'depth';
      
      fig = gcf;
      fig.Color = [1 1 1];
    end
    
  end
  
  methods (Static = true)
    function showNodes(this,~)
      % Callback for 'Show All' button
      n = this.UserData.Node;
      for i = 1:length(n)
        n(i).Visible = true;
      end
    end
    
    function hideNodes(this,~)
      % Callback for 'Hide All' button
      n = this.UserData.Node;
      for i = 1:length(n)
        n(i).Visible = false;
      end
    end
    
    function CloseRequestFcn(this,~)
      % Callback for figure CloseRequestFcn
      c = this.UserData;
      for i = 1:length(c.Node)
        delete(c.Node(i));
      end
      delete(gcf);
    end
    
  end
  
end