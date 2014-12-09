linetrace
=========

Simple and accurate 2D and 3D line tracing through a cartesian grid.
`[path, ds] = traceLine2D(p1, p2)` returns the intersection points (`path`)
of the line from `p1` to `p2` with the standard cartesian grid as well 
as the lengths of the line segments (`ds`) from one intersection point 
to another such that `eucdist(path(i+1) - path(i)) = ds(i)`.

`traceLine3D` works similarly in 3D. Both functions are not particularly
optimized for speed. The default `MAX_PATH_LEN = 1500` can be adapted to
larger grids and the `EPS = 1e-6` should be tuned per application.
