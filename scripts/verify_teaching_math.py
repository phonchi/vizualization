"""Independent numerical checks of the mathematical corrections in the rewrite."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from gdms_toolkit.teaching import learning_catalog

ROOT=Path(__file__).resolve().parents[1]
rows=[]

def check(name,actual,expected,tol=1e-8):
    passed=bool(np.isfinite(actual) and abs(actual-expected)<=tol)
    rows.append(dict(name=name,actual=float(actual),expected=float(expected),tolerance=tol,passed=passed))
    if not passed:raise AssertionError(rows[-1])

# Supercritical offspring mean does not imply a divergent finite-window likelihood.
integral=.1+2*quad(lambda u:(1+u)**-2,0,.5)[0]
check('n=2 finite-window compensator',integral,.1+2/3)
check('finite-window log likelihood',np.log(.1)-integral,np.log(.1)-(.1+2/3))

for p in [.8,1.,1.2]:
    c,T=.2,10.
    analytic=np.log((T+c)/c) if p==1 else ((T+c)**(1-p)-c**(1-p))/(1-p)
    check(f'finite Omori integral p={p}',quad(lambda t:(t+c)**-p,0,T)[0],analytic)

# Normalized planar spatial kernel: include the polar Jacobian.
q,sigma=1.7,4.
check('spatial kernel normalization',quad(lambda r:2*np.pi*r*(q-1)/(np.pi*sigma)*(1+r*r/sigma)**-q,0,np.inf)[0],1.)

# Discrete magnitude likelihood score is zero at the stated geometric MLE.
dm,u=.1,.35
beta=np.log1p(dm/u)/dm
check('discrete magnitude MLE score per event',-u+dm/np.expm1(beta*dm),0.)
C=1+dm/u
N=100
bhat=np.log(C)/(dm*np.log(10))
d=np.sqrt(C/N)
lo=np.log((C+d)/(1+d))/(dm*np.log(10))
hi=np.log((C-d)/(1-d))/(dm*np.log(10))
assert lo<bhat<hi and N>C
rows.append(dict(name='asymmetric magnitude interval brackets estimate',passed=True,lower=lo,estimate=bhat,upper=hi))

# Gamma-Poisson predictive variance, independently sampled.
rng=np.random.default_rng(20260914)
a,b,h=5.,7.,10.
mus=rng.gamma(a,1/b,500000)
future=rng.poisson(h*mus)
expected=h*a/b+h*h*a/(b*b)
check('posterior predictive total variance',future.var(),expected,.15)

# Concavity permits a boundary optimum for a mixture.
# One event at lambda1=4, lambda2=1; both total integrated intensities equal 5.
weights=np.linspace(0,1,1001)
scores=np.log(weights*4+(1-weights)*1)-5
check('concave mixture maximum may be on boundary',weights[np.argmax(scores)],1.)

for kind in ['regular','poisson','cluster']:
    frame=learning_catalog(kind)
    assert len(frame)==80 and frame.day.is_monotonic_increasing
    assert frame.day.between(0,120,inclusive='left').all() and (frame.m>=3).all()
    assert frame.equals(learning_catalog(kind))
    rows.append(dict(name=f'{kind} shared illustrative catalog provenance',passed=True,events=len(frame)))

out=ROOT/'reference/notes/rewrite_20260911/math_checks.json'
out.write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(rows,ensure_ascii=False,indent=2))
