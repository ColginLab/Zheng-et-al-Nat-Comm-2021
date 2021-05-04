 function [nol] = getnumberofcells_cz_v1(file)
nol = 0;
% file = 'cells-CA1.txt';
if exist(file,'file')
     
    tets = nan(numel(textread(file,'%1c%*[^\n]')),1); %#ok<REMFF1>
%     disp(pwd);
    fid = fopen(file,'r');
    if fid == -1
        msgbox('Could not open the input file','ERROR');
        nol = 0;
    else
        for i = 1:size(tets,1)
          tline = fgetl(fid);
          %disp(tline);
          if ~ischar(tline) 
              break 
          else
              if tline(4)=='_'
                tets(i) = str2double(tline(3));
              elseif tline(5) == '_'
                tets(i) = str2double(tline(3:4));
              end
          end          
        end

        fclose(fid);

        tets(isnan(tets)) = [];
        nol = size(tets,1);
    end
else
    nol = 0;
end