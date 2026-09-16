"""Retain complete road-enclosed blocks intersecting the existing 250m AOI.
OSM surface-road centreline blocks are a working boundary, not cadastral parcels.
No artificial closure is inserted across missing roads.
"""
from pathlib import Path
import json,gzip
from shapely.geometry import shape,LineString,mapping,box
from shapely.ops import transform,unary_union,polygonize
from pyproj import Transformer
P=Path(__file__).resolve().parents[1];O=P/'whole-block-results';O.mkdir(exist_ok=True)
fwd=Transformer.from_crs(4326,3826,always_xy=True).transform;rev=Transformer.from_crs(3826,4326,always_xy=True).transform
old=transform(fwd,shape(json.loads((P/'study-boundary-review.geojson').read_text())['features'][0]['geometry']))
src=json.load(gzip.open(P.parent/'maps/infrastructure/full-corridor/geometry.json.gz'))
roads=[]
for w in src['elements']:
 t=w.get('tags',{});h=t.get('highway');g=w.get('geometry',[])
 if h not in ['primary','secondary','tertiary','residential','unclassified','living_street','pedestrian','service','primary_link','secondary_link','tertiary_link']:continue
 if t.get('bridge')=='yes' or t.get('tunnel')=='yes' or t.get('layer','0') not in ['0',0] or t.get('area')=='yes':continue
 if t.get('service') in ['parking_aisle','driveway','drive-through'] or t.get('access')=='private':continue
 if len(g)>1:roads.append(transform(fwd,LineString([(v['lon'],v['lat']) for v in g])))
blocks=list(polygonize(unary_union(roads)))
selected=[b for b in blocks if b.intersection(old).area>0.001]
new=unary_union([old,*selected]);added=new.difference(old)
features=[{'type':'Feature','properties':{'id':i+1,'area_m2':round(b.area,2),'basis':'OSM surface road enclosed block'},'geometry':mapping(transform(rev,b))} for i,b in enumerate(selected)]
(O/'retained-blocks.geojson').write_text(json.dumps({'type':'FeatureCollection','features':features}))
(O/'study-boundary.geojson').write_text(json.dumps({'type':'FeatureCollection','features':[{'type':'Feature','properties':{'rule':'250m buffer union all intersected complete road-enclosed blocks','area_m2':new.area},'geometry':mapping(transform(rev,new))}]}))
# Preserve existing classifications; new land must remain explicitly unclassified.
classified=json.load(gzip.open(P/'full-line-results/cover-screening.geojson.gz'))
classified['features'].append({'type':'Feature','properties':{'class':4,'label':'expanded_area_not_yet_classified'},'geometry':mapping(transform(rev,added))})
with gzip.open(O/'vegetation.geojson.gz','wt') as f:json.dump(classified,f,separators=(',',':'))
l,b,r,t=new.bounds

def path(g):
 ps=list(g.geoms) if g.geom_type=='MultiPolygon' else [g];s=[]
 for p in ps:
  if p.geom_type!='Polygon':continue
  for ring in [p.exterior,*p.interiors]:s.append('M'+'L'.join(f'{x-l:.2f},{t-y:.2f}' for x,y in ring.coords)+'Z')
 return ''.join(s)
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="-20 -20 {r-l+40:.2f} {t-b+40:.2f}" width="2600"><title>市民大道 250m 選取完整街廓；紫灰為新增未判讀</title><g id="study-area"><path fill="#dadacf" fill-rule="evenodd" d="{path(new)}"/></g>']
for c,color in [(4,'#bcb4ca'),(3,'#e5ba6c'),(2,'#276f46')]:
 svg.append(f'<g id="class-{c}" fill="{color}" fill-rule="evenodd">')
 for ft in classified['features']:
  if ft['properties']['class']==c:svg.append(f'<path d="{path(transform(fwd,shape(ft["geometry"])))}"/>')
 svg.append('</g>')
svg.append('</svg>');(O/'whole-block-vegetation.svg').write_text(''.join(svg))
bound=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="-20 -20 {r-l+40:.2f} {t-b+40:.2f}" width="2600"><g fill="#e7e8df" stroke="#738372" stroke-width="1" fill-rule="evenodd">']
for bb in selected:bound.append(f'<path d="{path(bb)}"/>')
bound.append(f'</g><path d="{path(old)}" fill="none" stroke="#be824c" stroke-width="2" stroke-dasharray="8 6"/></svg>');(O/'block-boundary.svg').write_text(''.join(bound))
# Area checks are geometric, independent of raster sampling.
assert old.difference(new).area<.01
assert all(bb.difference(new).area<.01 for bb in selected)
s={'selected_block_count':len(selected),'old_area_m2':old.area,'new_area_m2':new.area,'added_unclassified_m2':added.area,'rule':'retain all complete OSM road-enclosed polygons intersecting 250m AOI by >0.001m² (numerical tolerance); union original buffer to retain roads','minimum_overlap_m2':0.001,'limitations':['OSM road-centreline block approximation, not cadastral or surveyed curb boundary','missing road geometry cannot be artificially closed','expanded area not classified; old 10.36% must not be reused for this new extent']}
(O/'summary.json').write_text(json.dumps(s,ensure_ascii=False,indent=2))
print(json.dumps(s,ensure_ascii=False,indent=2));print('Largest selected blocks m2',sorted([round(v.area) for v in selected],reverse=True)[:10])
