from pathlib import Path
import json,numpy as np,rasterio
from rasterio.warp import transform_bounds,reproject,Resampling
from rasterio.windows import from_bounds
from rasterio.features import rasterize
from shapely.geometry import shape
from shapely.ops import transform
from pyproj import Transformer
P=Path(__file__).resolve().parents[1];O=P/'sentinel'
item=json.load(open(O/'S2B_51RUH_20260603_0_L2A.json'))
aoi=shape(json.load(open(P/'study-boundary-review.geojson'))['features'][0]['geometry'])
bands={};masks={};target=None
for key in ['red','nir','green','blue','scl']:
 asset=item['assets'][key]
 with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.tif',GDAL_HTTP_TIMEOUT='30'):
  with rasterio.open(asset['href']) as src:
   bb=transform_bounds(4326,src.crs,*aoi.buffer(.0005).bounds);w=from_bounds(*bb,src.transform).round_offsets().round_lengths();raw=src.read(1,window=w);tr=src.window_transform(w)
   if target is None:target=(tr,src.crs,raw.shape)
   dest=np.zeros(target[2],raw.dtype);reproject(raw,dest,src_transform=tr,src_crs=src.crs,dst_transform=target[0],dst_crs=target[1],resampling=Resampling.nearest)
   masks[key]=dest!=0
   if key!='scl':
    info=asset['raster:bands'][0]
    offset=0 if item['properties'].get('earthsearch:boa_offset_applied') is True else info.get('offset',0)
    bands[key]=dest.astype('float32')*info.get('scale',1)+offset
   else:bands[key]=dest
 print('Read',key,flush=True)
tr,crs,sh=target
poly=transform(Transformer.from_crs(4326,crs,always_xy=True).transform,aoi);inside=rasterize([(poly,1)],out_shape=sh,transform=tr).astype(bool)
valid=np.isin(bands['scl'],[4,5,6])&masks['red']&masks['nir']&(bands['red']>=0)&(bands['nir']>=0)&((bands['red']+bands['nir'])>0)
ndvi=np.full(sh,-9999.,'float32');np.divide(bands['nir']-bands['red'],bands['nir']+bands['red'],out=ndvi,where=valid)
profile=dict(driver='GTiff',height=sh[0],width=sh[1],count=1,dtype='float32',crs=crs,transform=tr,nodata=-9999,compress='deflate')
with rasterio.open(O/'ndvi-2026-06-03.tif','w',**profile) as out:out.write(ndvi,1)
stats={'scene_id':item['id'],'date':item['properties']['datetime'],'native_pixel_m':10,'scl_resolution_m':20,'scale_offset':'scale 0.0001; earthsearch:boa_offset_applied=true, so no second -0.1 offset','metric':'share of 10m pixel centres within AOI with NDVI above threshold; NOT fractional green cover','aoi_pixel_count':int(inside.sum()),'valid_pixel_count':int((inside&valid).sum()),'unknown_pixel_count':int((inside&~valid).sum()),'thresholds':{}}
for q in [.2,.3,.4,.5]:
 n=int((inside&valid&(ndvi>=q)).sum());stats['thresholds'][str(q)]={'pixel_count':n,'nominal_area_m2':n*100,'percent_of_aoi_pixels':n/inside.sum()*100,'percent_of_valid_pixels':n/(inside&valid).sum()*100}
(O/'ndvi-statistics.json').write_text(json.dumps(stats,indent=2));print(json.dumps(stats,indent=2))
from PIL import Image
rgb=np.dstack([bands[k] for k in ['red','green','blue']]);rgb=np.clip(rgb*3,0,1);Image.fromarray((rgb*255).astype('uint8')).save(O/'sentinel-true-color.png')
