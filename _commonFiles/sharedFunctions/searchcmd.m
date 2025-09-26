function outcell=searchcmd(cmd,keyword,endkey,endkey2,endkey3)
% searches text for keywords and takes chars till stopps with endkey

idx = strfind(cmd,keyword);

if ~exist('endkey2','var')
       for iit=1:size(idx,2)
            idxt=0;
            while ~strcmp(cmd(idx(iit)+size(keyword,2)+idxt),endkey)
                outcell{iit,1}(idxt+1)=cmd(idx(iit)+size(keyword,2)+idxt);
                idxt=idxt+1;
            end
       end
elseif ~exist('endkey3','var')
       for iit=1:size(idx,2)
            idxt=0;
            while ~strcmp(cmd(idx(iit)+size(keyword,2)+idxt),endkey) && ~strcmp(cmd(idx(iit)+size(keyword,2)+idxt),endkey2)
                outcell{iit,1}(idxt+1)=cmd(idx(iit)+size(keyword,2)+idxt);
                idxt=idxt+1;
            end
       end
else
       for iit=1:size(idx,2)
            idxt=0;
            while ~strcmp(cmd(idx(iit)+size(keyword,2)+idxt),endkey) && ~strcmp(cmd(idx(iit)+size(keyword,2)+idxt),endkey2) && ~strcmp(cmd(idx(iit)+size(keyword,2)+idxt),endkey3)
                outcell{iit,1}(idxt+1)=cmd(idx(iit)+size(keyword,2)+idxt);
                idxt=idxt+1;
            end
       end
end