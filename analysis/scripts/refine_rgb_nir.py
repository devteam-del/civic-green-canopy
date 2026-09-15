"""Use explicitly reviewed reference points; spatial holdout by 1km panel.
The previous draft classifier is rejected and not used for final area figures.
"""
from pathlib import Path
import json,numpy as np
from PIL import Image
from scipy.ndimage import uniform_filter
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
P=Path(__file__).resolve().parents[1];D=P/'private-ortho-review';O=P/'visible-cover'
refs=json.load(open(O/'reviewed-reference-points.json'))
def feats(rgb):
 a=rgb.astype('float32')/255.;r,g,b=a[:,:,0],a[:,:,1],a[:,:,2];den=np.maximum(r+g+b,.01);lum=a.mean(2);m=uniform_filter(lum,5);sd=np.sqrt(np.maximum(uniform_filter(lum**2,5)-m*m,0))
 return np.dstack([r,g,b,r/den,g/den,b/den,g-r,g-b,2*g-r-b,m,sd]).astype('float32')
manifest=json.load(open(D/'manifest.json'));fs={}
import rasterio
from rasterio.warp import reproject,Resampling
from rasterio.transform import from_origin
with rasterio.open(P/'sentinel/ndvi-2026-06-03.tif') as src:
 for r in manifest:
  f=feats(np.array(Image.open(D/r['file']).convert('RGB')))
  d=np.full(f.shape[:2],-9999.,dtype='float32');x0,y0,x1,y1=r['bbox_epsg3826']
  reproject(rasterio.band(src,1),d,src_transform=src.transform,src_crs=src.crs,dst_transform=from_origin(x0,y1,1,1),dst_crs='EPSG:3826',src_nodata=-9999,dst_nodata=-9999,resampling=Resampling.nearest)
  fs[r['id']]=np.dstack([f,np.where(d==-9999,-1,d),d==-9999]).astype('float32')
X=[];y=[];groups=[]
for r in refs:
 if r['visual_label']==3:continue
 f=fs[r['panel']];x=r['x'];yy=r['y']
 for dy in [-1,0,1]:
  for dx in [-1,0,1]:
   if 0<=yy+dy<f.shape[0] and 0<=x+dx<f.shape[1]:X.append(f[yy+dy,x+dx]);y.append(r['visual_label']);groups.append(r['panel'])
X=np.array(X);y=np.array(y);groups=np.array(groups)
def fit(xx,yy):
 m=RandomForestClassifier(n_estimators=180,min_samples_leaf=3,max_depth=8,random_state=20260915,n_jobs=2)
 return m.fit(xx,yy)
for panel in range(1,14):
 m=fit(X[groups!=panel],y[groups!=panel])
 for r in refs:
  if r['panel']==panel and r['visual_label']!=3:
   r['fusion_holdout_score']=float(m.predict_proba(fs[panel][r['y'],r['x']][None,:])[0,1])
   r['fusion_holdout_label']=int(r['fusion_holdout_score']>=.5)
valid=[r for r in refs if r['visual_label']!=3]
cm=confusion_matrix([r['visual_label'] for r in valid],[r['fusion_holdout_label'] for r in valid],labels=[0,1])
result={'method':'leave-one-panel-out; no training pixels from held-out panel','reference':'visual interpretation of 2025 edition RGB; not field survey','sample_count':len(valid),'ambiguous_reference_count':len(refs)-len(valid),'unweighted_confusion_matrix_rows_reference_0_1':cm.tolist(),'unweighted_accuracy':float(np.trace(cm)/cm.sum()),'vegetation_precision':float(cm[1,1]/max(cm[:,1].sum(),1)),'vegetation_recall':float(cm[1,1]/max(cm[1,:].sum(),1)),'limitation':'stratified small visual reference sample; not an area-weighted independent field accuracy or calibrated uncertainty interval'}
(O/'fusion-validation.json').write_text(json.dumps(result,indent=2));(O/'reviewed-reference-points.json').write_text(json.dumps(refs,indent=2))
model=fit(X,y)
for rec in manifest:
 n=rec['id'];f=fs[n];score=model.predict_proba(f.reshape(-1,f.shape[2]))[:,1].reshape(f.shape[:2]);np.save(D/f'{n:02d}-fusion-score.npy',score.astype('float32'))
print(json.dumps(result,indent=2))
