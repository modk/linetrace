function [path, ds] = traceLine3D( p1, p2, discretization)
% Accurately trace a 3D line from point1 to point2 returning all intersection
% points with the standard cartesian grid and the lengths of the corresponding
% line segments.
% If a discretization factor is provided the points are scaled before
% calculating the lengths.
%
% Return values: 
%	path - list of intersection points beginning with point1 and
%		ending with point2
%	ds - length of each line segment 
%		(such that ds(i) = euclidian_dist(path(i)-path(i-1)), ds(1) = 0)

MAX_PATH_LEN = 1500;
EPS = 1e-6;
INFTY = 999999;

% check if discretization was given
if ~exist('discretization', 'var')
	discretization = [1 1 1];
end

dir = p2 - p1;
dir = dir / norm( dir);

% tan of slope alpha (projected on z)
ta = dir( 2) / dir( 1);

% tan slope beta (projected on y)
tb = dir( 3) / dir( 1);

% (projected on x)
tc = dir( 2) / dir( 3);

pts = zeros( MAX_PATH_LEN, 3);
ds = zeros( 1, MAX_PATH_LEN);
pts( 1, :) = p1;
ds( 1) = 0;
cnt = 1;

curpt = p1;

% determine directions (sgx > 0 means movement to the "right" etc.)
sgx = sign( dir( 1));
sgy = sign( dir( 2));
sgz = sign( dir( 3));

if abs( dir( 1)) < EPS
		sgx = 0;
end
if abs( dir( 2)) < EPS
		sgy = 0;
end
if abs( dir( 3)) < EPS
		sgz = 0;
end

while abs( norm( curpt - p2)) > EPS
		
		% abort if too many steps
		if cnt > MAX_PATH_LEN
				error( 'cnt > %d', MAX_PATH_LEN);
		end
		
		cnt = cnt + 1;
		
		% check if final voxel boundaries of final pt reached
		cc = 0;
		if sgx > 0
				if floor( curpt( 1) + EPS) == floor( p2( 1))
						cc = cc + 1;
				end
		elseif sgx < 0
				if ceil( curpt( 1) - EPS) == ceil( p2( 1))
						cc = cc + 1;
				end
		else
				cc = cc + 1;
		end
		if sgy > 0
				if floor( curpt( 2) + EPS) == floor( p2( 2))
						cc = cc + 1;
				end
		elseif sgy < 0
				if ceil( curpt( 2) - EPS) == ceil( p2( 2))
						cc = cc + 1;
				end
		else
				cc = cc + 1;
		end
		if sgz > 0
				if floor( curpt( 3) + EPS) == floor( p2( 3))
						cc = cc + 1;
				end
		elseif sgz < 0
				if ceil( curpt( 3) - EPS) == ceil( p2( 3))
						cc = cc + 1;
				end
		else
				cc = cc + 1;
		end
		
		if cc == 3
				% we are at the final voxel, exit
				pts( cnt, :) = p2;
				ds( cnt) = norm( curpt - p2);
				break;
		end

		% distances to next grid intersection along the standard
		% base vectors
		a = INFTY;
		if sgx > 0
				a = ceil( curpt( 1) + EPS) - curpt( 1);
		elseif sgx < 0
				a = floor( curpt( 1) - EPS) - curpt( 1);
		end
		b = INFTY;
		if sgy > 0
				b = ceil( curpt( 2) + EPS) - curpt( 2);
		elseif sgy < 0
				b = floor( curpt( 2) - EPS) - curpt( 2);
		end
		c = INFTY;
		if sgz > 0
				c = ceil( curpt( 3) + EPS) - curpt( 3);
		elseif sgz < 0
				c = floor( curpt( 3) - EPS) - curpt( 3);
		end
		
		% calculate axis sections
		y1 = a * ta;
		z1 = a * tb;
		 
		x2 = b / ta;
		z2 = b / tc;

		x3 = c / tb;
		y3 = c * tc;
					 
		v = [a y1 z1; x2 b z2; x3 y3 c];

		% look for the intersection with the smallest distance
		% by sorting
		[tmp, ix] = sortrows( abs( v));
		% first element in sorted indices is the nearest intersection
		ix = ix( 1);
		
		if ix == 1
				nextpt = curpt + [a y1 z1];
		elseif ix == 2
				nextpt = curpt + [x2 b z2];
		else
				nextpt = curpt + [x3 y3 c];
		end
		
		pts( cnt, :) = nextpt;
		ds( cnt) = norm( ( curpt - nextpt) .* discretization);
		curpt = nextpt;
end

path = pts( 1 : cnt, :);
ds = ds( 1 : cnt);

end
