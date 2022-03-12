function [detA] = detcell33(A)

% det(A) = a11*(a22*a33 - a23*a32) - a12*(a21*a33 - a23*a31) + a13*(a21*a32-a22*a31)
detA = + A{1,1}.*(A{2,2}.*A{3,3} - A{2,3}.*A{3,2}) ...
       - A{1,2}.*(A{2,1}.*A{3,3} - A{2,3}.*A{3,1}) ...
       + A{1,3}.*(A{2,1}.*A{3,2} - A{2,2}.*A{3,1});