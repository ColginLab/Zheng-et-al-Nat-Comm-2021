

function plotI(Responses, MaxResponses)

%plots binary data
t = 1:length(Responses);
if(length(MaxResponses) == sum(MaxResponses == 1))
  hold on; [y, x] = find(Responses > 0);
  h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor','k');
  set(h, 'MarkerEdgeColor', 'k');
  hold on; [y, x] = find(Responses == 0);
  h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor', [0.75 0.75 0.75]);
  set(h, 'MarkerEdgeColor', 'k');
end
axis([1 t(end)  0 1.10]); box on;
  
  