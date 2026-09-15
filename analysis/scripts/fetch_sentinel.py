from pathlib import Path
import json,numpy as np,rasterio
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds
from rasterio.features import rasterize
from shapely.geometry import shape,mapping
from shapely.ops import transform
from pyproj import Transformer
P=Path(__file__).resolve().parents[1];O=P/'sentinel';O.mkdir(exist_ok=True)
items=json.load(open('/tmp/civic-sentinel-search.json'))['features']
selected=[f for f in items if any(d in f['id'] for d in ['20260531','20260802','20260829','20260603','20260703','20260414'])]
aoi=shape(json.load(open(P/'study-boundary-review.geojson'))['features'][0]['geometry'])
records=json.load(open(O/'scene-audit.json')) if (O/'scene-audit.json').exists() else []
for f in selected:
 if any(r['id']==f['id'] for r in records):continue
 asset=f['assets']['scl']
 with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.tif',GDAL_HTTP_TIMEOUT='30'):
  with rasterio.open(asset['href']) as src:
   bb=transform_bounds(4326,src.crs,*aoi.bounds);w=from_bounds(*bb,src.transform).round_offsets().round_lengths();a=src.read(1,window=w);tr=src.window_transform(w)
   poly=transform(Transformer.from_crs(4326,src.crs,always_xy=True).transform,aoi);inside=rasterize([(poly,1)],out_shape=a.shape,transform=tr).astype(bool)
   good=inside & np.isin(a,[4,5,6]);percent=float(good.sum()/inside.sum()*100)
   rec={'id':f['id'],'date':f['properties']['datetime'],'study_clear_percent':percent,'scl_counts':{str(i):int(((a==i)&inside).sum()) for i in np.unique(a[inside])}}
   records.append(rec);(O/(f['id']+'.json')).write_text(json.dumps(f))
 print(rec,flush=True)
(O/'scene-audit.json').write_text(json.dumps(records,indent=2))
