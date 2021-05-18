% Triangle Thresholding Algorithm
% pixelCounts is the image histogram
% side is either 'R' or 'Right' to find a threshold on the right side of the histogram peak,
%     or 'L' or 'Left' to find a threshold on the left side of the histogram peak.
% Mark Hayworth
function thresholdLevel = triangle_threshold(pixelCounts, side, showPlot)
%     This technique is due to Zack (Zack GW, Rogers WE, Latt SA (1977),
%     "Automatic measurement of sister chromatid exchange frequency",
%     J. Histochem. Cytochem. 25 (7): 741?53, )
%     A line is constructed between the maximum of the histogram at
%     (b) and the lowest (or highest depending on context) value (a) in the
%     histogram. The distance L normal to the line and between the line and
%     the histogram h[b] is computed for all values from a to b. The level
%     where the distance between the histogram and the line is maximal is the
%     threshold value (level). This technique is particularly effective
%     when the object pixels produce skewed unimodal histogram,
%     where it will find the "corner" of the histogram.
% Sample call:
% 	binOfThreshold = triangle_threshold(pixelCount, 'R', true);

try
	thresholdLevel = 1; % Initialize
	side = upper(side(1)); % Get just first letter and capitalize it.
	
	% Find maximum of histogram and its location along the x axis
	% Find x and y of the peak of the histogram.
	[yPeak, xPeak] = max(pixelCounts);
	
	% Find location of first and last non-zero values.
	% Values < (h / 10000) are considered zeros.
	nonZeroBins = find(pixelCounts > yPeak / 10000);
	if isempty(nonZeroBins)
		% Didn't find any.  Lower the threshold to absolute zero and try again.
		nonZeroBins = find(pixelCounts > 0);
		if isempty(nonZeroBins)
			warningMessage = sprintf('Could not find a threshold value.\nNo counts were found in the histogram');
			uiwait(warndlg(warningMessage));
			return;
		end
	end
	
	firstNonZeroBin = nonZeroBins(1);
	lastNonZeroBin = nonZeroBins(end);
	
	% Sometimes, through contrast stretching or whatever, there are gaps in the histogram.  Fill those.
	gaps = pixelCounts == 0;
	% Don't include flat parts on the left or right.
	gaps(1:firstNonZeroBin-1) = false;
	gaps(lastNonZeroBin+1 : end) = false;
	% 	pixelCounts = regionfill(pixelCounts, gaps); % Doesn't work because regionfill doesn't work with 1-D "images"
	xAll = 1 : length(pixelCounts);
	xGood = 1 : length(pixelCounts);
	xGood(gaps) = []; % Remove gaps
	yGood = pixelCounts;
	yGood(gaps) = [];
	% Fill gaps by interpolating in 1-D across the gaps.
	pixelCounts = interp1(xGood, yGood, xAll);
	
	% The max can occur over multiple bins (if histogram hump has a flat top),
	% so find the bin that is closest to the furthest non-zero bin on the side the user wants to threshold on.
	if side(1) == 'R'
		xPeak = find(pixelCounts == yPeak, 1, 'last');
		xm = lastNonZeroBin; % Point on far right side, at the end of the tail on the right.
		ym = pixelCounts(xm);
	else
		xPeak = find(pixelCounts == yPeak, 1, 'first');
		xm = firstNonZeroBin; % Point on far left side, at the end of the tail on the left.
		ym = pixelCounts(xm);
	end
	
	% Find the denominator for our point-to-line distance formula.
	% Reference: http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
	denominator = sqrt((xm - xPeak) ^ 2 + (ym - yPeak) ^ 2);
	
	% Get the points along the histogram
	xh = 1 : length(pixelCounts);
	yh = pixelCounts;
	
	if showPlot
		% Plot it
		hFig = figure;
		bar(xh, yh, 'BarWidth', 1);
		xlim([1, length(xh)]);
		grid on;
		% Maximize the window via undocumented Java call.
		% Reference: http://undocumentedmatlab.com/blog/minimize-maximize-figure-window
% 		MaximizeFigureWindow;
		% Plot a line from peak to end
		hold on;
		plot([xPeak, xm], [yPeak, ym], 'r-', 'LineWidth', 2);
		drawnow;
	end
	
	% Find the numerator for our point-to-line distance formula.
	numerator = abs((xm-xPeak)*(yPeak-yh) - (xPeak - xh) * (ym - yPeak));
	
	% Compute all the distances
	distances = numerator ./ denominator;
	
	% Sometimes because of the shape of the histogram, which bulges up, there are some bins that are above
	% the line from point m (the tail) to the peak.  This means the distance is being measured from the line
	% going upwards, above the hypoteneuse, instead of downwards.  This is not what we want.
	% We need to find those distances and zero them out (set to -infinity) so they won't be considered.
	slope = (yPeak-ym)/(xPeak-xm);
	yLine = slope * (xh - xm) + ym;
	% 	plot(xh(1:xPeak), yLine(1:xPeak), 'm-', 'LineWidth', 2)
	badDistanceIndexes = yh > yLine;
	distances(badDistanceIndexes) = -inf;
	
	% Zero out the distances on the side of the histogram that we are not considering.
	if side(1) == 'R'
		% Looking on the right (bright) side.  Zero out bins to the left of the peak.
		startBin = xPeak - 1;
		if startBin < 1
			startBin = 1;
		end
		distances(1:startBin) = 0;
		% Zero out bins to the right of the right-most point of the hypoteneuse.
		startBin = lastNonZeroBin + 1;
		if startBin > length(distances)
			startBin = length(distances);
		end
		distances(startBin:end) = 0;
	else
		% Looking on the left (dark) side.  Zero out bins to the right of the peak.
		startBin = xPeak + 1;
		if startBin > length(distances)
			startBin = length(distances);
		end
		distances(startBin:end) = 0;
		% Zero out bins to the left of the left-most point of the hypoteneuse.
		startBin = lastNonZeroBin - 1;
		if startBin < 1
			startBin = 1;
		end
		distances(startBin:end) = 0;
	end
	
	% Find the max distance.  Because we set bins above the hypoteneuse to -inf, they won't be found.
	[maxDistance, indexOfMaxDistance] = max(distances);
	thresholdLevel = indexOfMaxDistance;
	
	if showPlot
		% Find the slope and intercept of the peak to end line.
		coefficients = polyfit([xPeak, xm], [yPeak, ym], 1);
		m = coefficients(1); % Slope.
		b = coefficients(2); % Intercept.
		x0 = indexOfMaxDistance;
		y0 = yh(indexOfMaxDistance);
		xOnLine = (-b + m * x0 + y0) / (2*m);
		yOnLine = m * xOnLine + b;
		% Plot a line from peak to end
		hold on;
		plot([x0, xOnLine], [y0, yOnLine], 'r-', 'LineWidth', 2);
		plot([x0, x0], [0, y0], 'r-', 'LineWidth', 2);
		xlabel('Gray Level', 'FontSize', 20)
		ylabel('Pixel Count', 'FontSize', 20)
		title('Triangle Threshold Determination', 'FontSize', 20)
		str = sprintf('Threshold = %.1f Gray Levels', x0);
		text(x0 + 9, y0, str, 'FontSize', 20, 'Color', 'r', 'Rotation', 0)
		xl = xlim(); % Get range in x direction.
		yl = ylim(); % Get range in y direction.
		% Put up "Longest Line" text angled over the line, more or less.
		text(x0 + 6, y0 + 0.15*yl(2), 'Longest Line', 'FontSize', 20, 'Color', 'r', 'Rotation', 50)
		drawnow;
		msgboxh('This is the threshold determination');
		if xl(2) <= 255
			% Set up tick marks every 10 on the x axis.
			xticks(0 : 10 : xl(2));
		end
		close(hFig);
	end
	
catch ME
	% Some error happened if you get here.
	errorMessage = sprintf('Error in program %s, function %s(), at line %d.\n\nError Message:\n%s', ...
		mfilename, ME.stack(1).name, ME.stack(1).line, ME.message);
	warning(errorMessage);
end
return;  % triangle_threshold

