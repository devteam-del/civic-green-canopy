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
manifest=json.load(open(D/'manifest.json'));fs={r['id']:feats(np.array(Image.open(D/r['file']).convert('RGB'))) for r in manifest}
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
   r['spatial_holdout_score']=float(m.predict_proba(fs[panel][r['y'],r['x']][None,:])[0,1])
   r['spatial_holdout_label']=int(r['spatial_holdout_score']>=.5)
valid=[r for r in refs if r['visual_label']!=3]
cm=confusion_matrix([r['visual_label'] for r in valid],[r['spatial_holdout_label'] for r in valid],labels=[0,1])
result={'method':'leave-one-panel-out; no training pixels from held-out panel','reference':'visual interpretation of 2025 edition RGB; not field survey','sample_count':len(valid),'ambiguous_reference_count':len(refs)-len(valid),'unweighted_confusion_matrix_rows_reference_0_1':cm.tolist(),'unweighted_accuracy':float(np.trace(cm)/cm.sum()),'vegetation_precision':float(cm[1,1]/max(cm[:,1].sum(),1)),'vegetation_recall':float(cm[1,1]/max(cm[1,:].sum(),1)),'limitation':'stratified small visual reference sample; not an area-weighted independent field accuracy or calibrated uncertainty interval'}
(O/'spatial-validation.json').write_text(json.dumps(result,indent=2));(O/'reviewed-reference-points.json').write_text(json.dumps(refs,indent=2))
model=fit(X,y)
for rec in manifest:
 n=rec['id'];f=fs[n];score=model.predict_proba(f.reshape(-1,f.shape[2]))[:,1].reshape(f.shape[:2]);np.save(D/f'{n:02d}-refined-score.npy',score.astype('float32'))
print(json.dumps(result,indent=2))
