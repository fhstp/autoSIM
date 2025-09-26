function [area,areaTri] = heronsFormula(vertices,faces)
% This function computes the area of all triangles by Heron's formula.

% Input vars
ver = vertices;
tri = faces;

% Side 1
s1 = ver(tri(:,2),:) - ver(tri(:,1),:);
s1 = sqrt(sum(s1.^2,2));

% Side 2
s2 = ver(tri(:,3),:) - ver(tri(:,2),:);
s2 = sqrt(sum(s2.^2,2));

% Side 3
s3 = ver(tri(:,1),:) - ver(tri(:,3),:);
s3 = sqrt(sum(s3.^2,2));

% Half circumference
s = (s1+s2+s3)/2;

% Area
areaTri = sqrt(s.*(s-s1).*(s-s2).*(s-s3));
area = sum(areaTri);
end

