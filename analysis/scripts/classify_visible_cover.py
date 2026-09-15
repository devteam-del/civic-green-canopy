"""Interpret 2025 RGB imagery; output candidate vegetation and uncertainty.
Manual visual patch labels, never land-use masks. Not field-validated cover.
The classifier score is not a calibrated probability or confidence interval.
"""
from pathlib import Path
import json
import numpy as np
from PIL import Image
from scipy.ndimage import uniform_filter
from sklearn.ensemble import RandomForestClassifier
from shapely.geometry import shape,box,mapping
from shapely.ops import transform
from pyproj import Transformer
import rasterio
from rasterio.features import rasterize,shapes
from rasterio.transform import from_origin

P=Path(__file__).resolve().parents[1];D=P/'private-ortho-review'
OUT=P/'visible-cover';OUT.mkdir(exist_ok=True)
manifest=json.load(open(D/'manifest.json'))
# Each point denotes a 5x5 pixel patch, selected by visual inspection of
# original 1000px-wide official image panels. 1=vegetation, 0=other surface.
patches={
1:{1:[(545,271),(555,408),(271,619),(557,680),(97,18)],0:[(389,124),(238,388),(732,478),(742,652),(832,264)]},
2:{1:[(357,459),(407,482),(318,436),(578,602),(656,623),(703,608),(594,471)],0:[(756,469),(959,714),(867,471),(727,602),(239,179),(361,395),(580,79),(553,44)]},
3:{1:[(938,461),(942,490),(730,397),(488,409)],0:[(701,288),(595,316),(269,92),(660,225),(772,321),(425,313),(171,366),(540,436)]},
4:{1:[(96,382),(274,445),(343,490),(433,474),(211,257),(325,216),(622,13)],0:[(579,29),(624,763),(282,560),(813,112),(782,515),(667,585),(615,88),(565,445)]},
5:{1:[(251,150),(294,20),(407,340),(705,369),(716,359),(304,452)],0:[(220,458),(260,459),(597,260),(595,249),(839,114),(44,3),(862,532),(149,325)]},
6:{1:[(400,122),(398,222),(407,424),(373,561),(552,117),(561,101)],0:[(673,249),(852,278),(553,433),(832,48),(233,363),(637,103),(444,565)]},
7:{1:[(223,340),(651,352),(645,385),(747,715),(708,791),(408,672)],0:[(602,452),(677,495),(426,654),(534,840),(687,47),(514,258),(735,69),(294,123)]},
8:{1:[(597,263),(833,245),(963,244),(414,449),(493,549),(780,411)],0:[(460,510),(980,310),(562,313),(589,405),(455,201),(146,350),(213,345),(795,529)]},
9:{1:[(714,301),(657,345),(591,387),(870,284),(470,157)],0:[(37,19),(410,150),(195,135),(715,74),(864,134),(351,324),(627,443),(827,348)]},
10:{1:[(503,237),(570,257),(714,265),(850,205),(302,143),(937,165)],0:[(936,135),(889,28),(551,107),(141,53),(37,367),(192,295),(552,366),(272,497)]},
11:{1:[(745,599),(464,621),(617,378),(743,354),(481,116),(909,489)],0:[(362,283),(529,377),(726,275),(798,376),(899,18),(254,544),(194,243)]},
12:{1:[(262,299),(208,651),(493,639),(708,558),(471,179),(607,135)],0:[(551,390),(771,88),(308,311),(373,107),(62,42),(822,454),(498,265)]},
13:{1:[(289,228),(741,65),(89,505),(419,440),(555,42)],0:[(397,61),(111,92),(80,359),(507,336),(474,181),(126,276)]}
}

def features(rgb):
 a=rgb.astype('float32')/255.;r,g,b=a[:,:,0],a[:,:,1],a[:,:,2]
 den=np.maximum(r+g+b,0.01);lum=a.mean(2)
 mean=uniform_filter(lum,5);std=np.sqrt(np.maximum(uniform_filter(lum**2,5)-mean**2,0))
 return np.dstack([r,g,b,r/den,g/den,b/den,g-r,g-b,2*g-r-b,mean,std]).astype('float32')
X=[];Y=[]
for rec in manifest:
 rgb=np.array(Image.open(D/rec['file']).convert('RGB'));f=features(rgb)
 for label,points in patches[rec['id']].items():
  for x,y in points:
   if y>=rgb.shape[0] or x>=rgb.shape[1]:raise ValueError((rec['id'],x,y))
   patch=f[y-2:y+3,x-2:x+3].reshape(-1,f.shape[2]);X.append(patch);Y.extend([label]*len(patch))
model=RandomForestClassifier(n_estimators=160,min_samples_leaf=4,max_depth=14,class_weight='balanced',random_state=20260915,n_jobs=2)
model.fit(np.concatenate(X),Y)
(OUT/'manual-training-patches.json').write_text(json.dumps(patches,indent=2))
fwd=Transformer.from_crs(4326,3826,always_xy=True).transform
rev=Transformer.from_crs(3826,4326,always_xy=True).transform
aoi=transform(fwd,shape(json.load(open(P/'study-boundary-review.geojson'))['features'][0]['geometry']))
result=[]
for rec in manifest:
 rgba=np.array(Image.open(D/rec['file']).convert('RGBA'));rgb=rgba[:,:,:3];f=features(rgb)
 score=model.predict_proba(f.reshape(-1,f.shape[2]))[:,1].reshape(rgb.shape[:2])
 np.save(D/f'{rec["id"]:02d}-score.npy',score.astype('float32'))
 x0,y0,x1,y1=rec['bbox_epsg3826'];tr=from_origin(x0,y1,1,1)
 inside=rasterize([(aoi,1)],out_shape=score.shape,transform=tr).astype(bool)
 # Unobservable dark patches and intentionally pixelated source areas are
 # retained as uncertainty, not filled using the old height model.
 dark=rgb.max(2)<48
 pixelated=np.zeros_like(inside)
 if rec['id']==10:
  for coords in [[(0,426),(42,426),(109,495),(72,543),(0,543)],[(639,441),(701,441),(730,513),(848,560),(790,647),(790,690),(501,690),(552,611),(615,578)]]:
   from shapely.geometry import Polygon
   pp=Polygon([(x0+x,y1-y) for x,y in coords]);pixelated|=rasterize([(pp,1)],out_shape=score.shape,transform=tr).astype(bool)
 if rec['id']==11:pixelated[-30:,:650]=True
 # 0=outside;1=nonvegetation candidate;2=vegetation candidate;3=uncertain.
 cls=np.zeros(score.shape,'uint8');cls[inside]=3
 usable=inside & (rgba[:,:,3]>0) & ~dark & ~pixelated
 cls[usable & (score<=0.4)]=1;cls[usable & (score>=0.65)]=2
 Image.fromarray(cls).save(OUT/f'{rec["id"]:02d}-classes.png')
 with rasterio.open(OUT/f'{rec["id"]:02d}-classes.tif','w',driver='GTiff',width=cls.shape[1],height=cls.shape[0],count=1,dtype='uint8',crs='EPSG:3826',transform=tr,nodata=0,compress='deflate') as dst:dst.write(cls,1)
 area=int(inside.sum());counts={str(k):int((cls==k).sum()) for k in [1,2,3]}
 result.append({'id':rec['id'],'bbox_epsg3826':rec['bbox_epsg3826'],'area_m2':area,'class_m2':counts,'candidate_vegetation_pct':counts['2']/area*100 if area else None,'uncertain_pct':counts['3']/area*100 if area else None,'status':'full panel screened; classifier output pending independent validation'})
 # Review overlay stays local with source imagery; not redistributed.
 overlay=rgb.copy().astype(float);overlay[cls==2]=overlay[cls==2]*.4+np.array([60,245,65])*.6;overlay[cls==3]=overlay[cls==3]*.7+np.array([255,150,40])*.3;overlay[~inside]=overlay[~inside]*.3+255*.7
 Image.fromarray(overlay.astype('uint8')).save(D/f'{rec["id"]:02d}-review.png')
(OUT/'panel-statistics.json').write_text(json.dumps(result,indent=2))
print(json.dumps({'panels':len(result),'total_area_m2':sum(r['area_m2'] for r in result),'classes':{k:sum(r['class_m2'][k] for r in result) for k in ['1','2','3']}},indent=2))
