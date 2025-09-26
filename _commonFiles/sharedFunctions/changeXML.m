function changeXML(xmlpath, tag, value, times, start)
% Mark Simonlehner 13.05.2019

% saves an edited *.xml file 
% INPUT: xmlPath, Node name, value to change to, and how often this should
% be repeated

%Change log:
% Brian Horsak 30.04.2020: changed start node level and added start as optional input variable

% -------------------------------------------------------------------------
% Set the optinal variable to 1 when not assigned
if nargin == 4; start = 1; end
   
for i=start:times
DOMnode = xmlread(xmlpath);
thisListItem=DOMnode.getElementsByTagName(tag);
thisListItem=thisListItem.item(i-1);
thisListItem.setTextContent(value);
xmlwrite(xmlpath,DOMnode);
end