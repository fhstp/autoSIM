function nor = compute_triangle_normal(ver,tri)
%% Compute Triangle Normal from vertices
%==========================================================================
%Author: Colin Smith
%Date: 1/19/2015
%--------------------------------------------------------------------------
%Compute the triangle normals for each triangle using the ver and nor from 
%load_asc. Follows the procedure in TRIANGLE_NORMAL in contacts_v4.c.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Inputs
%------
%ver : [nVer, 3] 
%   Matrix of vertex coordinates [x1 y1 z1; x2 y2 z2; ...]
%tri : [nTri, 3]
%   Matrix of triangles defined by three vertices 
%   [ver11 ver12 ver13; ver21 ver22 ver23; ...]
%   
%Outputs
%-------
%nor : [nTri, 3]
%   Matrix of triangle normals
%   [N1x N1y N1z; N2x N2y N2z; ...]
%==========================================================================
[nTri, ~] = size(tri);

nor = zeros(nTri,3);

for i = 1:nTri
    %Get coordinates of each vertex
    ver1 = ver(tri(i,1),:);
    ver2 = ver(tri(i,2),:);
    ver3 = ver(tri(i,3),:);
    
    %Draw Vector along two sides
    e1 = ver3-ver1;
    e2 = ver2-ver1;
    
    %Take cross product of two sides and normalize to get unit normal vector
    Nx = e1(2)*e2(3)-e1(3)*e2(2); 
    Ny = e1(3)*e2(1)-e1(1)*e2(3);
    Nz = e1(1)*e2(2)-e1(2)*e2(1);
    
    mag = sqrt(Nx^2+Ny^2+Nz^2);
    nor(i,:) = [Nx/mag Ny/mag Nz/mag];
end


