"""Preliminary canopy proxy; does not measure total current vegetation cover.
Run with rasterio, shapely, pyproj and numpy. Inputs are retained OSM and CHMv2.
"""
from pathlib import Path
import gzip, json, re
import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from shapely.geometry import LineString, mapping, shape, box
from shapely.ops import transform, unary_union
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
fwd = Transformer.from_crs(4326,3826,always_xy=True).transform
rev = Transformer.from_crs(3826,4326,always_xy=True).transform
data=json.load(gzip.open(ROOT.parent/'maps/infrastructure/full-corridor/geometry.json.gz'))
ways=[w for w in data['elements'] if re.fullmatch('市民大道[一二三四五六七八]段',w.get('tags',{}).get('name','')) and w['tags'].get('highway')=='tertiary']
lines=[transform(fwd,LineString([(p['lon'],p['lat']) for p in w['geometry']])) for w in ways]
surface_west=min(l.bounds[0] for l in lines)
western=[]
for w in data['elements']:
    t=w.get('tags',{})
    if t.get('name') in ['市民大道高架道路','市民大道'] and t.get('highway')=='trunk':
        line=transform(fwd,LineString([(p['lon'],p['lat']) for p in w['geometry']]))
        clipped=line.intersection(box(-1e7,-1e7,surface_west,1e7))
        if not clipped.is_empty: western.append(clipped)
roads=unary_union(lines)
roads=unary_union([roads,*western])
x0,y0,x1,y1=roads.bounds
# Approximate the boulevard axis with the midpoint of outer carriageway
# intersections at 20 m eastings. This remains a REVIEW boundary, not survey data.
xs=np.arange(x0+0.1,x1,20.)
coords=[];missing=[];spreads=[]
for x in xs:
    cut=roads.intersection(LineString([(x,y0-10),(x,y1+10)]))
    points=list(cut.geoms) if hasattr(cut,'geoms') else [cut]
    ys=[p.y for p in points if p.geom_type=='Point' and not p.is_empty]
    if ys:
        coords.append((float(x),float((min(ys)+max(ys))/2)))
        spreads.append(max(ys)-min(ys))
    else: missing.append(float(x))
axis=LineString(coords).simplify(3)
aoi=axis.buffer(250,cap_style=2)
def save_geo(name,features):
    (ROOT/name).write_text(json.dumps({'type':'FeatureCollection','features':features},ensure_ascii=False))
save_geo('study-boundary-review.geojson',[{'type':'Feature','properties':{'status':'approximate; needs endpoint and route review','buffer_m':250,'cap':'flat'},'geometry':mapping(transform(rev,aoi))},{'type':'Feature','properties':{'role':'approximate centerline'},'geometry':mapping(transform(rev,axis))}])
save_geo('source-road-alignment.geojson',[{'type':'Feature','properties':{'osm_id':w['id'],'name':w['tags']['name']},'geometry':mapping(transform(rev,l))} for w,l in zip(ways,lines)])
left,bottom,right,top=aoi.bounds
left,bottom=np.floor([left,bottom]);right,top=np.ceil([right,top])
dst_transform=from_origin(left,top,1,1)
height,width=int(top-bottom),int(right-left)
dst=np.full((height,width),255,dtype=np.uint8)
valid=np.zeros_like(dst)
with rasterio.open(ROOT/'canopy-height-v2-context.tif') as src:
    for source,target in [(src.read(1),dst),(src.read_masks(1),valid)]:
        reproject(source,target,src_transform=src.transform,src_crs=src.crs,dst_transform=dst_transform,dst_crs='EPSG:3826',resampling=Resampling.nearest,dst_nodata=255 if target is dst else 0)
inside=rasterize([(aoi,1)],out_shape=dst.shape,transform=dst_transform).astype(bool)
observed=inside & (valid>0) & (dst!=255)
stats={'status':'preliminary historical canopy proxy, NOT current total green cover','crs':'EPSG:3826','resampling':'nearest; 1m target grid, not increased source detail','source_pixel_mercator_m':1.1943285669558463,'axis_length_m':axis.length,'vector_area_m2':aoi.area,'pixel_area_m2':int(inside.sum()),'valid_m2':int(observed.sum()),'unknown_m2':int((inside & ~observed).sum()),'source_road_count':len(ways),'missing_axis_samples':len(missing),'max_road_intersection_spread_m':max(spreads),'thresholds':{}}
for threshold in [1,2,5]:
    n=int((observed & (dst>=threshold)).sum())
    stats['thresholds'][str(threshold)]={'definition':f'height >= {threshold} m','area_m2':n,'percent_of_aoi':n/int(inside.sum())*100,'percent_of_valid':n/int(observed.sum())*100}
features=json.load(open(ROOT/'taipei-v2-imagery-metadata.geojson'))['features']
dates=[]
for f in features:
    intersection=transform(fwd,shape(f['geometry'])).intersection(aoi)
    if not intersection.is_empty:dates.append({'date':f['properties']['acq_date'],'intersection_m2':intersection.area})
stats['imagery_date_intersections']=dates
(ROOT/'preliminary-canopy-statistics.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2))
with rasterio.open(ROOT/'canopy-v2-study-review.tif','w',driver='GTiff',width=width,height=height,count=1,dtype='uint8',crs='EPSG:3826',transform=dst_transform,nodata=255,compress='deflate',tiled=True) as out:
    out.write(np.where(observed,dst,255).astype('uint8'),1)
print(json.dumps(stats,ensure_ascii=False,indent=2))
