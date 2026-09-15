"""Small sequential study-area exports. Do not redistribute source imagery."""
from pathlib import Path
import json,math,time,urllib.parse,subprocess
from shapely.geometry import shape,box
from shapely.ops import transform
from pyproj import Transformer
P=Path(__file__).resolve().parents[1]
OUT=P/'private-ortho-review'; OUT.mkdir(exist_ok=True)
fwd=Transformer.from_crs(4326,3826,always_xy=True).transform
aoi=transform(fwd,shape(json.load(open(P/'study-boundary-review.geojson'))['features'][0]['geometry']))
l,b,r,t=aoi.bounds
base='https://www.historygis.udd.gov.taipei/arcgis/rest/services/Aerial/Ortho_2025/MapServer/export'
manifest=[]
for i,x in enumerate(range(math.floor(l/1000)*1000,math.ceil(r/1000)*1000,1000)):
    part=aoi.intersection(box(x,b-1,x+1000,t+1))
    if part.is_empty:continue
    bb=(x,math.floor(part.bounds[1]/10)*10,x+1000,math.ceil(part.bounds[3]/10)*10)
    wh=(1000,int(bb[3]-bb[1]))
    params={'bbox':','.join(map(str,bb)),'bboxSR':3826,'imageSR':3826,'size':','.join(map(str,wh)),'format':'png32','transparent':'true','layers':'show:0','f':'image'}
    url=base+'?'+urllib.parse.urlencode(params,safe=',:')
    name=f'{i+1:02d}.png'; path=OUT/name
    if not path.exists():
        subprocess.run(['curl','-sS','-L','--fail','--max-time','60',url,'-o',str(path)],check=True)
        time.sleep(3)
    manifest.append({'id':i+1,'file':name,'bbox_epsg3826':bb,'width':wh[0],'height':wh[1],'requested_pixel_m':1,'edition':2025,'capture_date':'not provided by service','url':url})
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print('Saved panel',i+1,flush=True)
