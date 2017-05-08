module InPolygon

export inpolygon

# Copyright 2000 softSurfer, 2012 Dan Sunday
# This code may be freely used and modified for any purpose
# providing that this copyright notice is included with it.
# SoftSurfer makes no warranty for this code, and cannot be held
# liable for any real or imagined damage resulting from its use.
# Users of this code must verify correctness for their application.
 
# a Point is defined by its coordinates {int x, y;}
#

type Point
	x::Real
	y::Real
end

"""
# isLeft(): tests if a point is Left|On|Right of an infinite line.
#    Input:  three points P0, P1, and P2
#    Return: >0 for P2 left of the line through P0 and P1
#            =0 for P2  on the line
#            <0 for P2  right of the line
#    See: Algorithm 1 "Area of Triangles and Polygons
"""
@inline function isLeft(P0::Point, P1::Point, P2::Point)::Int
    sign( (P1.x - P0.x) * (P2.y - P0.y) - (P2.x -  P0.x) * (P1.y - P0.y) )
end

"""
# cn_PnPoly(): crossing number test for a point in a polygon
#      Input:   P = a point,
#               V[] = vertex points of a polygon V[1::n] with V[n]=V[1]
#      Return:  0 = outside, 1 = inside
# This code is patterned after [Franklin, 2000]
"""
function cn_PnPoly(P::Point, V::Vector{Point})::Int
  cn = 0    # the  crossing number counter
  n = length(V)

  vi = V[n]
  viP = vi.y - P.y
  # loop through all edges of the polygon
  for vii in V  # edge from V[i]  to V[i+1]
    viiP = vii.y - P.y	
    if (viP <= 0) != (viiP <= 0)
      # compute  the actual edge-ray intersect x-coordinate
	  vt = viP / (viP - viiP)
      if P.x <  vi.x + vt * (vii.x - vi.x) # P.x < intersect
        cn += 1   # a valid crossing of y=P.y right of P.x
      end
    end
	vi = vii
	viP = viiP
  end
  cn    # 0 if even (out), and 1 if  odd (in)
end
#

"""
# wn_PnPoly(): winding number test for a point in a polygon
#      Input:   P = a point,
#               V[] = vertex points of a polygon V[1:n] with V[n]=V[1]
#      Return:  wn = the winding number (=0 only when P is outside)
"""
function wn_PnPoly(P::Point, V::Vector{Point})::Int
  wn = 0    # the  winding number counter
  n = length(V)

  vi = V[n]
  # loop through all edges of the polygon
  for vii in V   # edge from V[i]  to V[i+1]
    if vi.y <= P.y         # start y <= P.y
      if vii.y  > P.y      # an upward crossing
        if isLeft(vi, vii, P) > 0  # P left of  edge
          wn += 1             # have  a valid up intersect
        end
      end
    else                      # start y > P.y (no test needed)
      if vii.y  <= P.y     # a downward crossing
        if isLeft(vi, vii, P) < 0  # P right of  edge
          wn -= 1             # have  a valid down intersect
	    end
      end
	end
	vi = vii
  end
  return wn
end

## -*- texinfo -*-
## @deftypefn  {Function File} {@var{in} =} inpolygon (@var{x}, @var{y}, @var{xv}, @var{yv})
## @deftypefnx {Function File} {[@var{in}, @var{on}] =} inpolygon (@var{x}, @var{y}, @var{xv}, @var{yv})
##
## For a polygon defined by vertex points @code{(@var{xv}, @var{yv})}, return
## true if the points @code{(@var{x}, @var{y})} are inside (or on the boundary)
## of the polygon; Otherwise, return false.
##
## The input variables @var{x} and @var{y}, must have the same dimension.
##
## The optional output @var{on} returns true if the points are exactly on the
## polygon edge, and false otherwise.
## @seealso{delaunay}
## @end deftypefn

## Author: Frederick (Rick) A Niles <niles@rickniles.com>
## Created: 14 November 2006

## Vectorized by Søren Hauberg <soren@hauberg.org>
"""
## The method for determining if a point is in in a polygon is based on
## the algorithm shown on
## http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/
## and is credited to Randolph Franklin.
"""
function inpolygon{T <: Real}(x, y, xv::Vector{T}, yv::Vector{T})

  if size(x) != size(y)
    error("inpolygon: X and Y must be real arrays of the same size");
  elseif size(xv) != size(yv)
    error("inpolygon: XV and YV must be real vectors of the same size");
  end

  npol = length(xv);
  in = falses(x);
  on = falses(x);

  j = npol;
  for i = 1 : npol
    delta_xv = xv[j] - xv[i];
    delta_yv = yv[j] - yv[i];
    ## distance = [distance from (x,y) to edge] * length(edge)
	distance = delta_xv .* (y - yv[i]) - (x - xv[i]) .* delta_yv;

    ## is y between the y-values of edge i,j AND (x,y) on the left of the edge?
	t1 = ((yv[i] .<= y) & (y .< yv[j]))
    t2 = ((yv[j] .<= y) & (y .< yv[i]))
	t3 = (0 .< distance*delta_yv)

	idx1 = (((yv[i] .<= y) & (y .< yv[j])) |
            ((yv[j] .<= y) & (y .< yv[i]))  ) &
	       (0 .< distance*delta_yv);

		   println("$j $i: ($t1 | $t2 ) & $t3 = $idx1")
    in[idx1] = !in[idx1];

    ## Check if (x,y) are actually on the boundary of the polygon.
	idx2 = ((((yv[i] .<= y) & (y .<= yv[j])) | ((yv[j] .<= y) & (y .<= yv[i])))
		 & (((xv[i] .<= x) & (x .<= xv[j])) | ((xv[j] .<= x) & (x .<= xv[i])))
		 & ((0 .== distance) | (0 != delta_xv)));
    on[idx2] = true;

    j = i;
  end

  ## Matlab definition include both in polygon and on polygon points.
  in |= on;
  (in, on)

end


end # module
