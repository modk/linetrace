function [path, ds] = traceLine2D( p1, p2)
% Accurately trace a line from point1 to point2 returning all intersection
% points with the standard cartesian grid and the lengths of the corresponding
% line segments.
%
% Return values: 
%	path - list of intersection points beginning with point1 and
%		ending with point2
%	ds - length of each line segment 
%		(such that ds(i) = euclidian_dist(path(i)-path(i-1)), ds(1) = 0)

MAX_PATH_LEN = 1500;
EPS = 1e-6;
INFTY = 999999;

dir = p2 - p1;
dir = dir / norm( dir);

% tan of slope alpha
ta = dir( 2) / dir( 1);

pts = zeros( MAX_PATH_LEN, 2);
ds = zeros( 1, MAX_PATH_LEN);
pts( 1, :) = p1;
ds( 1) = 0;
cnt = 1;

curpt = p1;

% determine directions (sgx > 0 means movement to the "right" etc.)
sgx = sign( dir( 1));
sgy = sign( dir( 2));
if abs( dir( 1)) < EPS
	sgx = 0;
end
if abs( dir( 2)) < EPS
	sgy = 0;
end

while abs( norm( curpt - p2)) > EPS

	% abort if too many steps
	if cnt > MAX_PATH_LEN
		error( 'cnt > %d', MAX_PATH_LEN);
	end

	cnt = cnt + 1;
	
	% check if voxel boundaries of final pt reached
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
	if cc == 2
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
	else
		a = floor( curpt( 1) - EPS) - curpt( 1);
	end
	b = INFTY;
	if sgy > 0
		b = ceil( curpt( 2) + EPS) - curpt( 2);
	else
		b = floor( curpt( 2) - EPS) - curpt( 2);
	end
	
	% calculate axis sections
	x = b / ta;
	if abs( dir( 2)) < EPS
		x = INFTY * sgx;
	end		 
	y = ta * a;
	if abs( dir( 1)) < EPS
		y = 0;
	end

	if a > x
		% step to the right
		nextpt = curpt + [a y];
	else
		% step to the top
		nextpt = curpt + [x b];
	end

	pts( cnt, :) = nextpt;
	ds( cnt) = norm( curpt - nextpt);
	curpt = nextpt;
end

path = pts( 1 : cnt, :);
ds = ds( 1 : cnt);

end
