```
git clone https://github.com/arrmgt/maneuvers.git
```
```
cd processing26
% Note this is a matlab script, not a function
% run as >> do10process1 % no arguments
edit do10process1.m
% Edit file *_raw.nc name(s).
% Edit segment indices
%    Note that 10 Hz data is expected. 
%	 Need start and stop point numbers for each segment.
%    Use 20260408 as a model
% All files/segments will be concatenated
NOTES:
1.  Process runs twice for ship static and boom static measurements, respectively.  Plots automatically be generated for each.
2.  Global attributes (AWINDS.*) are in the following files:
    maneuvers_BOOM.txt (uses BOOM static)
	maneuvers_SHIP.txt (uses SHIP static)
3.  The *xlsx spreadsheets have tabulated results and confidence limits.
```

