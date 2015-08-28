# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 09:05:48 2015

@author: heistermann
"""

from osgeo import ogr, osr
import wradlib
import pylab as plt
import numpy as np
from matplotlib.path import Path
from matplotlib.collections import PolyCollection


def points_in_polygon(polygon, points, buffer=0.):
    """Select points inside polygon
    """
    mpath = Path( polygon )
    return  mpath.contains_points(points, radius=-buffer)


def subset_points(pts, bbox, buffer=0.):
    """Subset a large set of points by polygon bbox
    """
    x = pts[:,0]
    y = pts[:,1]
    return np.where(
            (x >= bbox["left"]  -buffer) & \
            (x <= bbox["right"] +buffer) & \
            (y >= bbox["bottom"]-buffer) & \
            (y <= bbox["top"]   +buffer) )[0]
        

def get_bbox(x,y, buffer=0.):
    """Return dictionary of bbox
    """
    return dict(left=np.min(x), 
                right=np.max(x), 
                bottom=np.min(y), 
                top=np.max(y))
    

if __name__ == '__main__':

    # Get coordinates
    grid_xy_radolan = wradlib.georef.get_radolan_grid(900, 900)
    x_radolan = grid_xy_radolan[:, :, 0]
    y_radolan = grid_xy_radolan[:, :, 1]
    
    # create radolan projection osr object
    proj_stereo = wradlib.georef.create_osr("dwd-radolan")

    # create Gauss Krueger zone 4 projection osr object
    proj_gk = osr.SpatialReference()
    proj_gk.ImportFromEPSG(31468)

    # transform radolan polar stereographic projection to gk
    xy = wradlib.georef.reproject(grid_xy_radolan,
                                  projection_source=proj_stereo,
                                  projection_target=proj_gk)

    # Open shapefile (already in GK4)
    shpfile = "mulde/Mulde.shp"
    dataset, inLayer = wradlib.io.open_shape(shpfile)
    
    # Reduce grid for enhancing performance
#    bbox = inLayer.GetExtent()
#    bbox = dict(left=bbox[0], right=bbox[1], bottom=bbox[2], top=bbox[3])
#    mask = (x >= bbox["left"]) & (x <= bbox["right"]) & (y >= bbox["bottom"]) & (y <= bbox["top"])
#    si, se = np.where(mask)
#    x2 = x[si.min():si.max() + 1, se.min():se.max() + 1]
#    y2 = y[si.min():si.max() + 1, se.min():se.max() + 1]
#    rwdata2 = rwdata[si.min():si.max() + 1, se.min():se.max() + 1]
   
    
    # Compute indices for sub-catchments
    cats, keys = wradlib.georef.get_shape_coordinates(inLayer, key='GWKZ')
    xy_ = np.reshape(xy, (-1,2))

    pips = []  # these are those which we consider inside or near
    bboxs = [] # this is our pre-selection based on a simple bbox (just for comparison)    
    for cat in cats:
        # Pre-selection to increase performance 
        ix = subset_points(xy_, get_bbox(cat[:,0],cat[:,1]), buffer=500.)
        ixix = ix[points_in_polygon(cat, xy_[ix,:], buffer=500.)]
        if len(ixix)==0:
            # For very small catchments: increase buffer size
            ix = subset_points(xy_, get_bbox(cat[:,0],cat[:,1]), buffer=1000.)
            ixix = ix[points_in_polygon(cat, xy_[ix,:], buffer=1000.)]            
        pips.append( ixix )
        bboxs.append( ix )
      
    # Read and prepare the actual data (RADOLAN)
    ifile = "raa01-rw_10000-1408030950-dwd---bin.gz"
    data, attrs = wradlib.io.read_RADOLAN_composite(ifile, missing=np.nan)
    sec = attrs['secondary']
    data.flat[sec] = np.nan

    # Compute average areal rainfall based on the indices
    avg = np.array([])
    for i, cat in enumerate(cats):
        if len(pips[i])>0:
            avg = np.append(avg, np.nanmean(data.ravel()[pips[i]]) )
        else:
            avg = np.append(avg, np.nan )
            
    # Check if some catchments still are NaN
    invalids = np.where(np.isnan(avg))[0]
    assert len(invalids)==0, "Attention: No average rainfall computed for %d catchments" % len(invalids)
              
    # Plot average rainfall and original data
    fig = plt.figure(figsize=(14,8))
    # Average rainfall sum
    ax = fig.add_subplot(121, aspect="equal")
    wradlib.vis.add_lines(ax, cats, color='black', lw=0.5)
    coll = PolyCollection(cats, array=avg, cmap=plt.cm.jet, edgecolors='none')
    ax.add_collection(coll)
    ax.autoscale()
    cb = plt.colorbar(coll, ax=ax, shrink=0.5)
    cb.set_label("(mm/h)")
    plt.xlabel("GK4 Easting")
    plt.ylabel("GK4 Northing")
    plt.title("Areal average rain sums")
    plt.draw()
    # Original RADOLAN data
    ax1 = fig.add_subplot(122, aspect="equal")
    pm = plt.pcolormesh(xy[:, :, 0], xy[:, :, 1], np.ma.masked_invalid(data), vmax=2.)
    wradlib.vis.add_lines(ax1, cats, color='black', lw=0.5)
    bbox = inLayer.GetExtent()
    plt.xlim(ax.get_xlim())
    plt.ylim(ax.get_ylim())
    cb = plt.colorbar(pm, ax=ax1, shrink=0.5)
    cb.set_label("(mm/h)")
    plt.xlabel("GK4 Easting")
    plt.ylabel("GK4 Northing")
    plt.title("Original RADOLAN rain sums")
    plt.draw()
    plt.tight_layout()
    plt.savefig("gebietsniederschlag.png")

    # Inspect behaviour of the algorithms for example catchments
#    i = 0     
#    plt.plot(cats[i][:,0],cats[i][:,1])
#    plt.plot(xy_[bboxs[i],0],xy_[bboxs[i],1], "bo")
#    plt.plot(xy_[pips[i] ,0],xy_[pips[i] ,1], "ro")

    # Export table with areal average rainfall
    with open('gebietsniederschlag.txt', 'w') as f:
        f.write("gmkz\trainfall\n")
        for i, item in enumerate(avg):
            f.write("%s\t%.2f\n" % (keys[i],item))

    
