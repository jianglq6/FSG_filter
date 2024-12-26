function [vertices, val, faces] = patch_pre(nx, nz, xvec, zvec, data)
%  
% velocity layout:
%  v---0---v
%  |   |   |
%  0---v---0
%  |   |   |
%  v---0---v
% that is, there's no value on integer-half
%  and half-integer points.
%
% organize them to small patches centered on integer-integer 
%   and half-half points, for patch showing
%
%  v   0   v  
%    /   \ 
%  0   v   0
%    \   /
%  v   0   v
%
%

% int-point
vertices = [];
val = [];
faces = [];
i = 1;
for z = 1 : nz
    if mod(z, 2) == 0
        xcir = 2:2:nx;
    else
        xcir = 1:2:nx;
    end
    for x = xcir
        vertices(i,:) =  [xvec(x), zvec(z)];
        val(i,:) = data(x, z);
        i = i+1;
    end
end

half_nx = nx/2;
i = 1;
for h = 2:nz-1
    for w = 1:half_nx-1
        id = w + (h -1) * half_nx;
        if (mod(h, 2) == 0)
            face = [id, id - half_nx + 1, id + 1, id + half_nx + 1];
        else
            face = [id, id - half_nx, id + 1, id + half_nx];
        end
        faces(i, :) = face;
        i = i+1;
    end
end

end





