function [v, angle] = myncc(s1,s2)

v = norm(s1(:).'*s2(:))/(norm(s1(:))*norm(s2(:)));
angle = acos(v);