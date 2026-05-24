function catchME(ME)
% Basic info
fprintf('Error id: %s\n', ME.identifier);
fprintf('Message: %s\n\n', ME.message);

% Stack (topmost frame is the location where the error occurred)
for k = 1:numel(ME.stack)
    fprintf('Frame %d: %s\n  File: %s\n  Line: %d\n\n', ...
        k, ME.stack(k).name, ME.stack(k).file, ME.stack(k).line);
end

% Or get a single-line top frame
if ~isempty(ME.stack)
    top = ME.stack(1);
    fprintf('Top location: %s (line %d) in %s\n', top.name, top.line, top.file);
end

% Full formatted report (no hyperlinks)
disp(getReport(ME,'extended','hyperlinks','off'));

