"""Build full-corridor multisource screening with explicit unresolved pixels."""
from pathlib import Path
import json,gzip,numpy as np,rasterio
from PIL import Image
from rasterio.features import rasterize,shapes
from rasterio.warp import reproject,Resampling
from rasterio.transform import from_origin
from shapely.geometry import shape,box,Polygon,LineString,mapping
from shapely.ops import transform,unary_union
from pyproj import Transformer
P=Path(__file__).resolve().parents[1];D=P/'private-ortho-review';O=P/'full-line-results';O.mkdir(exist_ok=True)
fwd=Transformer.from_crs(4326,3826,always_xy=True).transform;rev=Transformer.from_crs(3826,4326,always_xy=True).transform
gj=json.load(open(P/'study-boundary-review.geojson'));aoi=transform(fwd,shape(gj['features'][0]['geometry']))
source=json.load(gzip.open(P.parent/'maps/infrastructure/full-corridor/geometry.json.gz'))
bridges=[]
for w in source['elements']:
 t=w.get('tags',{})
 if t.get('bridge')=='yes' and '市民大道' in t.get('name','') and t.get('highway') in ['trunk','trunk_link'] and w.get('geometry'):
  bridges.append(transform(fwd,LineString([(p['lon'],p['lat']) for p in w['geometry']])).buffer(12))
bridge_screen=unary_union(bridges).intersection(aoi)
manifest=json.load(open(D/'manifest.json'));records=[];geoms={2:[],3:[]}
# Material cannot be confidently resolved in these sports/roof patches.
material_boxes={2:[(207,152,267,245)],3:[(806,565,951,661)],5:[(190,438,289,492)],7:[(319,620,373,715)],8:[(412,477,496,540),(819,278,925,343)],9:[(392,128,443,179)],10:[(833,0,964,60),(927,133,991,184)],12:[(480,137,552,189)],13:[(369,362,451,428)]}
for rec in manifest:
 n=rec['id'];rgba=np.array(Image.open(D/rec['file']).convert('RGBA'));rgb=rgba[:,:,:3];score=np.load(D/f'{n:02d}-fusion-score.npy');h,w=score.shape;x0,y0,x1,y1=rec['bbox_epsg3826'];tr=from_origin(x0,y1,1,1)
 inside=rasterize([(aoi,1)],out_shape=(h,w),transform=tr).astype(bool)
 bridge=rasterize([(bridge_screen,1)],out_shape=(h,w),transform=tr).astype(bool)
 nd=np.full((h,w),-9999.,'float32')
 with rasterio.open(P/'sentinel/ndvi-2026-06-03.tif') as src:reproject(rasterio.band(src,1),nd,src_transform=src.transform,src_crs=src.crs,dst_transform=tr,dst_crs='EPSG:3826',src_nodata=-9999,dst_nodata=-9999,resampling=Resampling.nearest)
 obscured=(rgba[:,:,3]==0)|(rgb.max(2)<48)|(nd==-9999)|bridge
 material=np.zeros((h,w),bool)
 for xx0,yy0,xx1,yy1 in material_boxes.get(n,[]):material[yy0:yy1,xx0:xx1]=True
 if n==10:
  coords_list=[[(0,426),(42,426),(109,495),(72,543),(0,543)],[(639,441),(701,441),(730,513),(848,560),(790,647),(790,690),(501,690),(552,611),(615,578)]]
  for coords in coords_list:obscured|=rasterize([(Polygon([(x0+x,y1-y) for x,y in coords]),1)],out_shape=(h,w),transform=tr).astype(bool)
 if n==11:obscured[-30:,:650]=True
 hard_unknown=inside&(obscured|material)
 cls=np.zeros((h,w),'uint8');cls[inside]=3
 usable=inside&~hard_unknown
 cls[usable&(score<=.35)]=1;cls[usable&(score>=.65)]=2
 # Main mapped areas distinguish uncertain classification from physical masking.
 counts={str(k):int((cls==k).sum()) for k in [1,2,3]};total=int(inside.sum())
 records.append({'panel':n,'bbox_epsg3826':rec['bbox_epsg3826'],'study_area_m2':total,'model_vegetation_m2':counts['2'],'model_nonvegetation_m2':counts['1'],'unresolved_m2':counts['3'],'physically_or_material_unresolved_m2':int(hard_unknown.sum()),'model_vegetation_pct_of_total':counts['2']/total*100,'unresolved_pct':counts['3']/total*100,'binary_0_5_vegetation_pct_of_observable':float(((score>=.5)&usable).sum()/max(usable.sum(),1)*100),'ortho_visual_review':'completed; full image viewed','status':'multisource screening; unresolved cells retained'})
 with rasterio.open(O/f'{n:02d}-classification.tif','w',driver='GTiff',width=w,height=h,count=1,dtype='uint8',crs='EPSG:3826',transform=tr,nodata=0,compress='deflate') as dst:dst.write(cls,1)
 for geom,val in shapes(cls,mask=(cls>=2),transform=tr):geoms[int(val)].append(shape(geom))
 colors=np.array([[248,248,243],[218,218,207],[39,111,70],[229,186,108]],dtype='uint8');Image.fromarray(colors[cls]).save(O/f'{n:02d}-classification.png')
 overlay=rgb.astype(float);overlay[cls==2]=overlay[cls==2]*.5+np.array([40,240,60])*.5;overlay[cls==3]=overlay[cls==3]*.7+np.array([245,166,50])*.3;overlay[~inside]=overlay[~inside]*.3+255*.7;Image.fromarray(overlay.astype('uint8')).save(D/f'{n:02d}-final-review.png')
total=sum(r['study_area_m2'] for r in records);veg=sum(r['model_vegetation_m2'] for r in records);unknown=sum(r['unresolved_m2'] for r in records)
summary={'scope':'full Civic Boulevard, both sides 250m; western branches retained separately','vector_area_m2':aoi.area,'raster_area_m2':total,'model_vegetation_m2':veg,'model_vegetation_pct_of_total':veg/total*100,'unresolved_m2':unknown,'unresolved_pct':unknown/total*100,'status':'full-corridor remote-sensing screening completed; not complete field verification or definitive total green-cover rate','epochs':{'ortho':'Taipei 2025 edition, exact exposure date not supplied','nir':'Sentinel-2 2026-06-03','height_reference_only':'2018-03-23 and 2018-08-07'},'thresholds':{'vegetation_score_min':.65,'nonvegetation_score_max':.35,'dark_rgb_max_below':48,'bridge_screen_buffer_each_side_m':12},'limitations':['classifier scores are not calibrated probabilities','bridge buffer is an occlusion screening approximation, not measured bridge footprint','1m output grid does not make 10m near-infrared data 1m resolution','do not treat unresolved cells as zero vegetation','no legal land-use vegetation exclusion mask used','not a single-date current census; roofs and ground overlap in plan view are counted once']}
(O/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));(O/'panels.json').write_text(json.dumps(records,indent=2))
features=[{'type':'Feature','properties':{'class':c,'label':'model_vegetation' if c==2 else 'unresolved'},'geometry':mapping(transform(rev,g))} for c,gs in geoms.items() for g in gs]
with gzip.open(O/'cover-screening.geojson.gz','wt',encoding='utf-8') as out:json.dump({'type':'FeatureCollection','features':features},out,separators=(',',':'))
# Map coordinates are metres in EPSG:3826; SVG preserves uniform scale.
l,b,r,t=aoi.bounds
def path(g):
 polys=list(g.geoms) if g.geom_type=='MultiPolygon' else [g];parts=[]
 for poly in polys:
  for ring in [poly.exterior,*poly.interiors]:parts.append('M'+'L'.join(f'{x-l:.1f},{t-y:.1f}' for x,y in ring.coords)+'Z')
 return ''.join(parts)
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {r-l:.1f} {t-b:.1f}" width="2400"><title>全線多源植被判讀；橙色為待核實，不是零綠覆</title><path fill="#dadacf" fill-rule="evenodd" d="{path(aoi)}"/>']
for c,color in [(3,'#e5ba6c'),(2,'#276f46')]:
 for g in geoms[c]:svg.append(f'<path fill="{color}" fill-rule="evenodd" d="{path(g)}"/>')
svg.append('</svg>');(O/'full-corridor-screening.svg').write_text(''.join(svg))
print(json.dumps(summary,ensure_ascii=False,indent=2));print('vector features',len(features))
