function [V,nr] = MY_con2vert(A,b)
% ----------------------------------------------------
% CON2VERT - convert a convex set of constraint inequalities into the set
%            of vertices at the intersections of those inequalities;i.e.,
%            solve the "vertex enumeration" problem. Additionally,
%            identify redundant entries in the list of inequalities.
% 
% V = con2vert(A,b)
% [V,nr] = con2vert(A,b)
% 
% Converts the polytope (convex polygon, polyhedron, etc.) defined by the
% system of inequalities A*x <= b into a list of vertices V. Each ROW
% of V is a vertex. For n variables:
% A = m x n matrix, where m >= n (m constraints, n variables)
% b = m x 1 vector (m constraints)
% V = p x n matrix (p vertices, n variables)
% nr = list of the rows in A which are NOT redundant constraints
flg = 0;
c = A\b;
tol = 1e-07;
if ~all(abs(A*c - b) < tol)
    %obj1 = @(c) obj(c, {A,b});
    [c,~,ef] = fminsearch( @(x)obj(x, {A,b}),c);
%     [c,~,ef] = fminsearch(@obj,c,A,b);
    if ef ~= 1
        flg = 1;
    end
end
if flg ==1
    V = [];
%     error('Unable to locate a point within the interior of a feasible region.')
else
b = b - A*c;
D = A ./ repmat(b,[1 size(A,2)]);
[~,v2] = convhulln([D;zeros(1,size(D,2))]);
[k,v1] = convhulln(D);
%if v2 > v1
%    error('Non-bounding constraints detected. (Consider box constraints on variables.)')
%end
nr = unique(k(:));
G  = zeros(size(k,1),size(D,2));
for ix = 1:size(k,1)
    F = D(k(ix,:),:);
    G(ix,:)=F\ones(size(F,1),1);
end
V = G + repmat(c',[size(G,1),1]);
[~,I]=unique(num2str(V,6),'rows');
V=V(I,:);
end
return
function d = obj(c,param)
A = param{1};
b = param{2};
% ,A,b
% A=params{1};
% b=params{2};
d = A*c-b;
k=(d>=-1e-15);
d(k)=d(k)+1;
d = max([0;d]);
return
