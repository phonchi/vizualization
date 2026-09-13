"""Bounded, independent checks of the teaching formulas and cached Italy inputs.

Run from the repository environment. This script never rebuilds forecast caches.
Evidence and export round-trip examples go to rewrite_20260913/continuation/.
"""

from _teaching_runtime import prepare_plotting_environment
prepare_plotting_environment()

from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
from math import factorial
from fractions import Fraction
from pathlib import Path
import json
import os

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.stats import gamma, norm, t

from gdms_toolkit import csep_teaching as csep, italy, italy_models as models
from gdms_toolkit.teaching import learning_catalog

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / os.environ.get('TEACHING_REPORT_DIR','reference/notes/refresh_20260913_exhibition/math')
ROWS = []
PROVENANCE = {}
SEED = 20260915


def check(name, actual, expected, tol=1e-8):
    actual, expected = float(actual), float(expected)
    passed = bool(np.isfinite(actual) and abs(actual - expected) <= tol)
    row = dict(name=name, actual=actual, expected=expected, tolerance=tol, passed=passed)
    ROWS.append(row)
    if not passed:
        raise AssertionError(row)


def require(name, condition, **details):
    row = dict(name=name, passed=bool(condition), **details)
    ROWS.append(row)
    if not condition:
        raise AssertionError(row)


def raises_value_error(name, callback):
    try:
        callback()
    except ValueError:
        require(name, True)
    else:
        require(name, False)


def digest(path):
    hasher = sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            hasher.update(block)
    return dict(path=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                bytes=path.stat().st_size, sha256=hasher.hexdigest())


def legacy_checks():
    integral = .1 + 2 * quad(lambda u: (1 + u) ** -2, 0, .5)[0]
    check('n=2 finite-window compensator', integral, .1 + 2/3)
    check('finite-window log likelihood', np.log(.1)-integral, np.log(.1)-(.1+2/3))
    for p in [.8, 1., 1.2]:
        c, horizon = .2, 10.
        analytic = (np.log((horizon+c)/c) if p == 1 else
                    ((horizon+c)**(1-p)-c**(1-p))/(1-p))
        check(f'finite Omori integral p={p}', quad(lambda u: (u+c)**-p, 0, horizon)[0], analytic)
    q, sigma = 1.7, 4.
    check('planar kernel normalization with polar Jacobian',
          quad(lambda r: 2*np.pi*r*(q-1)/(np.pi*sigma)*(1+r*r/sigma)**-q, 0, np.inf)[0], 1.)
    dm, u = .1, .35
    beta = np.log1p(dm/u)/dm
    check('discrete magnitude MLE score per event', -u+dm/np.expm1(beta*dm), 0.)
    C, n = 1+dm/u, 100
    bhat, d = np.log(C)/(dm*np.log(10)), np.sqrt(C/n)
    lo = np.log((C+d)/(1+d))/(dm*np.log(10))
    hi = np.log((C-d)/(1-d))/(dm*np.log(10))
    require('asymmetric magnitude interval brackets estimate', lo < bhat < hi and n > C,
            lower=lo, estimate=bhat, upper=hi)
    # Integrate conditional first/second moments: replaces the old 500,000 draws.
    a, b, h = 5., 7., 10.
    first = quad(lambda mu: h*mu*gamma.pdf(mu, a=a, scale=1/b), 0, np.inf)[0]
    second = quad(lambda mu: (h*mu+(h*mu)**2)*gamma.pdf(mu, a=a, scale=1/b), 0, np.inf)[0]
    check('Gamma-Poisson predictive total variance by quadrature', second-first**2,
          h*a/b+h*h*a/(b*b))
    weights = np.linspace(0, 1, 1001)
    check('concave mixture maximum may be on boundary',
          weights[np.argmax(np.log(weights*4+(1-weights))-5)], 1.)
    for kind in ['regular', 'poisson', 'cluster']:
        frame = learning_catalog(kind)
        require(f'{kind} illustrative catalog provenance',
                len(frame) == 80 and frame.day.is_monotonic_increasing
                and frame.day.between(0, 120, inclusive='left').all()
                and (frame.m >= 3).all() and frame.equals(learning_catalog(kind)), events=len(frame))


def score_checks():
    rate = np.array([0., .3, 1.7, 4.])
    count = np.array([0, 0, 2, 3])
    manual = np.array([0., -.3, -1.7+2*np.log(1.7)-np.log(2),
                       -4+3*np.log(4)-np.log(6)])
    check('POLL includes factorial at repeated events', np.max(abs(csep.poll(rate, count)-manual)), 0.)
    check('jPOLL equals independent full Poisson log likelihood', csep.jpoll(rate, count), manual.sum())
    require('positive event at zero rate has negative infinite score',
            np.isneginf(csep.poll(np.array([0.]), np.array([1]))[0]))
    mu, variance, obs = 4., 10., 6
    p, size = mu/variance, mu*mu/(variance-mu)
    check('NB parameterization mean', size*(1-p)/p, mu)
    check('NB parameterization variance', size*(1-p)/(p*p), variance)
    # Direct PMF recurrence, independent of scipy nbinom CDF/sf.
    masses = [p**size]
    for k in range(1, obs+1):
        masses.append(masses[-1]*(k-1+size)/k*(1-p))
    nb = csep.n_test_nb(np.array([mu]), np.array([obs]), variance)
    check('NB inclusive lower tail', nb['delta2'], sum(masses))
    check('NB inclusive upper tail', nb['delta1'], 1-sum(masses[:-1]))
    check('NB tails overlap by observed PMF', nb['delta1']+nb['delta2'], 1+masses[-1])
    raises_value_error('NB rejects variance not greater than mean',
                       lambda: csep.n_test_nb([mu], [obs], variance=mu))
    a, b, w = np.array([.5, 2., 3.]), np.array([1., 1.3, 2.]), np.array([2, 1, 3])
    gain = csep.information_gain(a, b, w)
    diff = sum(int(n)*np.log(x/y)-(x-y) for x, y, n in zip(a, b, w))
    check('IGPE equals independently expanded LL difference per event', gain['igpe'], diff/w.sum())
    check('IGPE agrees with jPOLL difference/N', gain['igpe'],
          (csep.jpoll(a, w)-csep.jpoll(b, w))/w.sum())
    event_diff = np.array([np.log(.5)]*2 + [np.log(2/1.3)] + [np.log(3/2)]*3)
    half = t.ppf(.975, len(event_diff)-1)*np.std(event_diff, ddof=1)/np.sqrt(len(event_diff))
    check('IGPE paired t interval halfwidth', (gain['ci'][1]-gain['ci'][0])/2, half)
    raises_value_error('IGPE rejects zero events', lambda: csep.information_gain(a, b, np.zeros(3, int)))
    for mixing in [0., .25, .6, 1.]:
        check(f'convex expected-count conservation weight={mixing}',
              (mixing*a+(1-mixing)*b).sum(), mixing*a.sum()+(1-mixing)*b.sum())
    # Same uniforms, independently allocated via cumulative intervals and dense counts.
    lam = np.array([0., .4, 1.2, 0., 2.4])
    ns = np.array([0, 1, 2, 9, 3, 6, 10, 1])
    sparse = csep._sim_jpoll(lam, ns, len(ns), np.random.default_rng(SEED))
    rng = np.random.default_rng(SEED)
    dense = []
    for n in ns:
        counts = np.zeros(lam.size, int)
        for uniform in rng.random(int(n)):
            index = next(j for j in range(lam.size) if uniform < sum(lam[:j+1])/sum(lam))
            counts[index] += 1
        value = -sum(lam)
        for intensity, multiplicity in zip(lam, counts):
            if multiplicity:
                value += int(multiplicity)*np.log(intensity)-np.log(factorial(int(multiplicity)))
        dense.append(value)
    check('sparse categorical jPOLL equals same-seed dense allocation', np.max(abs(sparse-dense)), 0.)
    zero_sim = csep._sim_jpoll(np.zeros(3), [0, 1], 2, np.random.default_rng(SEED))
    require('zero-rate simulation handles zero/positive counts', zero_sim[0] == 0 and np.isneginf(zero_sim[1]))


def kernel_checks():
    edges = italy.magnitude_edges()
    mass = models._gr_bin_mass(edges, italy.SPEC.mT)
    beta = italy.SPEC.b_value*np.log(10)
    expected = 1-10**(-italy.SPEC.b_value*(italy.SPEC.m_max-italy.SPEC.mT))
    check('SUP finite GR mass leaves upper tail unnormalized', mass.sum(), expected)
    check('SUP first GR bin by density quadrature', mass[0],
          quad(lambda m: beta*np.exp(-beta*(m-italy.SPEC.mT)), edges[0], edges[1])[0])
    P = models.PARAMS['EEPAS']
    magnitude, event_time = 3.7, 100.
    center = P['aT']+P['bT']*magnitude
    left, right = 10**(center-.3), 10**(center+.2)
    numeric = quad(lambda tau: np.exp(-.5*((np.log10(tau)-center)/P['sT'])**2)
                   /(tau*np.log(10)*P['sT']*np.sqrt(2*np.pi)), left, right)[0]
    actual = models.eepas_time_mass(np.array([event_time]), np.array([magnitude]),
                                   event_time+left, event_time+right)[0]
    check('EEPAS time CDF difference vs direct lognormal integral', actual, numeric)
    check('EEPAS zero-length time interval', models.eepas_time_mass(
        np.array([event_time]), np.array([magnitude]), event_time+left, event_time+left)[0], 0.)
    src = pd.DataFrame(dict(mb=[magnitude], x_km=[1.2], y_km=[-2.]))
    cells = pd.DataFrame(dict(x0=[-10.], x1=[7.], y0=[-6.], y1=[13.]))
    eta, mag, spatial = models.eepas_pieces(src, cells)
    sigma = P['sA']*10**(P['bA']*magnitude/2)
    xmass = quad(lambda x: np.exp(-.5*((x-1.2)/sigma)**2)/(sigma*np.sqrt(2*np.pi)), -10, 7)[0]
    ymass = quad(lambda y: np.exp(-.5*((y+2)/sigma)**2)/(sigma*np.sqrt(2*np.pi)), -6, 13)[0]
    check('EEPAS Gaussian rectangle vs separate density quadratures', spatial[0, 0], xmass*ymass)
    eta_expected = P['bM']*(1-P['mu'])*10**(-italy.SPEC.b_value*(P['aM']+(P['bM']-1)*magnitude))
    eta_expected *= np.exp(-.5*(beta*P['sM'])**2)
    check('EEPAS eta base-10 rate balance formula', eta[0], eta_expected)
    # The implemented magnitude integral deliberately uses five midpoint samples.
    def corrected_magnitude(m):
        delta = norm.cdf((m-P['aM']-P['bM']*italy.SPEC.m0-P['sM']**2*beta)/P['sM'])
        return norm.pdf((m-P['aM']-P['bM']*magnitude)/P['sM'])/P['sM']/delta
    midpoint = sum(corrected_magnitude(edges[0]+(i+.5)*.02) for i in range(5))*.02
    check('EEPAS magnitude five-midpoint implementation', mag[0, 0], midpoint)
    truth = quad(corrected_magnitude, edges[0], edges[1])[0]
    check('EEPAS magnitude midpoint accuracy in selected example', mag[0, 0], truth, 3e-4)



def model_concept_checks():
    """Small synthetic fixtures for new prose; not full forecasting implementations."""
    PROVENANCE['concept_fixtures'] = {
        'hawkes_scene': {'background': .3, 'amplitude': 1.5, 'decay_days': 1.3,
                         'event_days': [2., 5., 5.8], 'stationary_claim': False},
        'stationary_hawkes': {'background': .3, 'offspring_mean': .6,
                             'assumption': 'nonnegative integrable kernel, n < 1'},
        'stress_release': {'initial_state': .4, 'loading': .2, 'gamma': .7,
                           'reference_rate': .3, 'event_day': 2., 'release': .8,
                           'units': 'dimensionless synthetic state; days'},
        'score_distribution': {'counts': [0, 4], 'probabilities': [.75, .25],
                               'mean': 1., 'variance': 3.},
        'multiplicative_fixture': {'baseline': [2, 3, 5], 'conjugate': [0, 1, 2],
                                   'transformation': 'ln(1+z)', 'target_total': 10},
        'time_completeness': {
            'catalogue_lead_definition': 'target time minus catalogue start',
            'lag_definition': 'target time minus input cutoff; cutoff = issue minus delay',
            'endpoint_fixture': {'catalog_start': 0., 'issue': 100., 'delay': 10., 'target': 125., 'window_end': 150.},
            'input_magnitudes': [3., 4.], 'target_magnitudes': [5., 6.],
            'assumption': 'two-source quadrature fixture, no spatial boundaries; not full GR integration'},
        'eas_temporal_fixture': {'mainshock_rate': .2, 'direct_offspring_mean': .6,
                                 'decay_days': 2., 'forecast_days': 5.,
                                 'assumption': 'constant future mainshock rate, exponential first-generation kernel'},
    }
    amplitude, decay, background = 1.5, 1.3, .3
    check('Hawkes display exponential area is 1.95',
          quad(lambda lag: amplitude*np.exp(-lag/decay), 0, np.inf)[0], 1.95)
    require('Hawkes display fixture is not subcritical', amplitude*decay > 1)
    events = np.array([2., 5., 5.8])
    # Match the published scene expression, then compare scalar past-only sums.
    times = np.array([0., 2., 2.+1e-7, 5., 5.8, 6.])
    vector_rate = np.full(times.shape, background)
    for event in events:
        vector_rate += np.where(times > event,
                                amplitude*np.exp(-np.maximum(times-event, 0)/decay), 0.)
    scalar_rate = np.array([
        background+sum(amplitude*np.exp(-(now-event)/decay)
                       for event in events if event < now)
        for now in times])
    check('Hawkes scene vector and causal past-only sum agree',
          np.max(abs(vector_rate-scalar_rate)), 0.)
    check('Hawkes left rate at first event excludes self', vector_rate[1], background)
    check('Hawkes left rate at second event includes only earlier event',
          vector_rate[3], background+amplitude*np.exp(-3/decay))
    check('Hawkes left rate at third event includes two earlier events',
          vector_rate[4], background+amplitude*(np.exp(-3.8/decay)+np.exp(-.8/decay)))
    check('Hawkes immediate right contribution approaches kernel amplitude',
          vector_rate[2]-background, amplitude, 2e-7)
    # Separate subcritical example: never insert n=1.95 into this formula.
    n = .6
    by_generations = background*sum(n**generation for generation in range(100))
    mean_rate = background/(1-n)
    check('subcritical Hawkes mean vs summed immigration generations', mean_rate, by_generations)
    check('subcritical Hawkes mean satisfies stationary balance', mean_rate, background+n*mean_rate)
    require('stationary Hawkes fixture uses its own n below one', 0 <= n < 1)

    initial, loading, response, ref, event, release = .4, .2, .7, .3, 2., .8
    before_loading = ref*np.exp(response*initial)
    before_event = ref*np.exp(response*(initial+loading*event))
    after_event = ref*np.exp(response*(initial+loading*event-release))
    check('stress-release positive loading rate ratio', before_event/before_loading,
          np.exp(response*loading*event))
    require('stress-release positive loading raises rate', before_event > before_loading)
    check('stress-release event rate drop ratio', after_event/before_event, np.exp(-response*release))
    require('stress-release positive release lowers positive rate', 0 < after_event < before_event)
    after_state = initial+loading*event-release
    require('stress-release fixture does not reset to initial state', not np.isclose(after_state, initial))

    baseline = [Fraction(2), Fraction(3), Fraction(5)]
    multipliers = [Fraction(1), Fraction(2), Fraction(3)]
    unnormalized = [a*b for a, b in zip(baseline, multipliers)]
    require('multiplicative three-cell unnormalized exact total', sum(unnormalized) == 23)
    normalized = [value*Fraction(10, 23) for value in unnormalized]
    require('multiplicative three-cell exact normalized fractions',
            normalized == [Fraction(20, 23), Fraction(60, 23), Fraction(150, 23)]
            and sum(normalized) == 10)
    transformed = np.array([2., 3., 5.])*np.exp(np.log1p([0., 1., 2.]))
    check('multiplicative log-transform normalization vs exact fractions',
          np.max(abs(transformed*10/transformed.sum()-np.array([float(v) for v in normalized]))), 0.)

    # Deliberately non-Poisson: support {0, 4}, mean 1 and variance 3.
    values, probabilities = [0, 4], [Fraction(3, 4), Fraction(1, 4)]
    mu = sum(Fraction(value)*probability for value, probability in zip(values, probabilities))
    variance = sum(probability*(value-mu)**2 for value, probability in zip(values, probabilities))
    require('score fixture is overdispersed, non-Poisson with exact mean',
            mu == 1 and variance == 3)
    def score(rate, count):
        return count*np.log(rate)-rate-np.log(factorial(count))
    for candidate in [.25, .5, 1., 2., 4.]:
        difference = sum(float(probability)*(score(float(mu), value)-score(candidate, value))
                         for value, probability in zip(values, probabilities))
        analytic = candidate-float(mu)+float(mu)*np.log(float(mu)/candidate)
        check(f'non-Poisson exact expected score difference lambda={candidate}', difference, analytic)
        require(f'mean-consistent Poisson-form score lambda={candidate}', difference >= -1e-12)
    # Conditional mean optimality says nothing about Poisson tail calibration.

    catalog_start, issue, delay, target, window_end = 0., 100., 10., 125., 150.
    catalog_lead = target-catalog_start
    time_lag = target-(issue-delay)
    horizon = window_end-issue
    check('catalogue lead time is catalog start to target', catalog_lead, 125.)
    check('input delay is data cutoff to issue', issue-(issue-delay), 10.)
    check('time lag includes target offset after issue', time_lag, 35.)
    check('forecast horizon is issue to window end', horizon, 50.)
    check('oldest source lag endpoint recovers catalog start', target-catalog_lead, catalog_start)
    check('newest source lag endpoint recovers issue minus delay', target-time_lag, issue-delay)
    check('usable history length differs from catalogue lead with positive lag',
          catalog_lead-time_lag, issue-delay-catalog_start)

    P = models.PARAMS['EEPAS']
    source_mag = np.array([3., 4.])
    beta = italy.SPEC.b_value*np.log(10)
    lead = 10**(P['aT']+P['bT']*3.5)
    log_means = P['aT']+P['bT']*source_mag
    cdf_mass = models.eepas_time_mass(np.zeros(2), source_mag, 0., lead)
    # Independent integration in log10-time coordinates, not the erf code.
    integrated = np.array([
        quad(lambda u: np.exp(-.5*((u-center)/P['sT'])**2)/(P['sT']*np.sqrt(2*np.pi)),
             -np.inf, np.log10(lead))[0] for center in log_means])
    check('time completeness source CDF vs direct log-time quadrature',
          np.max(abs(cdf_mass-integrated)), 0., 1e-9)
    def completeness(target_magnitude, temporal_mass):
        weights = np.exp(-beta*source_mag)*np.exp(
            -.5*((target_magnitude-P['aM']-P['bM']*source_mag)/P['sM'])**2)
        # eta is independent of source magnitude for the fixed bM=1 NW fixture;
        # Gaussian normalization and other common factors cancel in the ratio.
        return np.dot(weights, temporal_mass)/weights.sum()
    require('completeness fixture uses constant eta condition bM=1', P['bM'] == 1.)
    for target_magnitude in [5., 6.]:
        fraction = completeness(target_magnitude, cdf_mass)
        check(f'time completeness weighted ratio M={target_magnitude}', fraction,
              completeness(target_magnitude, integrated))
        check(f'time completeness infinite-history ratio M={target_magnitude}',
              completeness(target_magnitude, np.ones(2)), 1.)
        require(f'time completeness lies in unit interval M={target_magnitude}', 0 <= fraction <= 1)
    require('specified EEPAS scaling gives lower completeness for larger target',
            completeness(6., cdf_mass) < completeness(5., cdf_mass))
    # Rate-balance identity for the two compensation endpoints; synthetic values.
    mix, fraction, full_rate = .18, .4, 2.
    finite_signal = (1-mix)*fraction*full_rate
    smooth = (mix+(1-mix)*(1-fraction))*full_rate+finite_signal
    signal = mix*full_rate+finite_signal/fraction
    check('EEPAS smooth compensation restores assumed average', smooth, full_rate)
    check('EEPAS signal compensation restores assumed average', signal, full_rate)

    main_rate, direct_mean, decay, horizon = .2, .6, 2., 5.
    offspring_rate = quad(lambda u: main_rate*direct_mean/decay*np.exp(-(horizon-u)/decay),
                          0, horizon)[0]
    check('EAS synthetic future-mainshock temporal convolution', offspring_rate,
          main_rate*direct_mean*(1-np.exp(-horizon/decay)))
    expected_children = quad(lambda u: main_rate*direct_mean*(1-np.exp(-(horizon-u)/decay)),
                             0, horizon)[0]
    check('EAS finite-window expected direct children by convolution',
          expected_children, main_rate*direct_mean*(horizon-decay*(1-np.exp(-horizon/decay))))
    require('EAS finite-window fixture does not count all infinite-time offspring',
            0 < expected_children < main_rate*horizon*direct_mean)
    PROVENANCE['concept_fixture_limits'] = [
        'These synthetic checks verify selected formulas, not Hawkes, stress-release or EAS forecast implementations.',
        'The n=1.95 display and n=0.6 stationary examples are distinct parameter settings.',
        'Two-magnitude completeness checks do not verify the continuous input-magnitude integral or spatial boundaries.',
        'Non-Poisson mean-score consistency does not establish Poisson-calibrated tests or confidence intervals.',
        'No new forecast cache or scientific fitted result is produced by these checks.',
    ]

def boundary_checks():
    # Include a gap between windows to distinguish start and end checks.
    windows = pd.DataFrame(dict(t1=[10., 30.], t2=[20., 40.]))
    targets = pd.DataFrame(dict(
        t_days=[10., 19.999, 20., 30., 40., 9., 31., 31., 31., 31., 31.],
        cell=[0, 0, 0, 1, 0, 0, -1, 177, 0, 0, 0],
        mb=[5., 5.1, 5., 7.4, 5., 5., 5., 5., 7.5, 4.9, 5.]))
    observed = italy.bin_targets(targets, windows)
    expected = np.zeros((2, 177, 25), int)
    expected[0, 0, 0] = 1
    expected[0, 0, 1] = 1
    expected[1, 1, 24] = 1
    expected[1, 0, 0] = 1
    require('bin_targets half-open bins, gaps, region and 7.5 exclusion', np.array_equal(observed, expected))
    adjacent = pd.DataFrame(dict(t1=[10., 20.], t2=[20., 30.]))
    event = pd.DataFrame(dict(t_days=[20.], cell=[0], mb=[5.]))
    counts = italy.bin_targets(event, adjacent)
    require('shared time boundary belongs to right window', counts[0].sum() == 0 and counts[1, 0, 0] == 1)


def cache_checks():
    cache_dir = ROOT / 'data/cache/italy_forecasts'
    files = [cache_dir/f'{model}_{period}.npy' for model in models.MODELS for period in ['learning', 'testing']]
    for path in files:
        require(f'precomputed cache exists: {path.name}', path.is_file())
    clean = ROOT / 'data/cache/italy/horus_clean.csv'
    regrid = ROOT / 'data/cache/italy/regrid_01deg.npz'
    require('clean HORUS exists; no download permitted', clean.is_file())
    require('regrid cache exists; no reconstruction permitted', regrid.is_file())
    PROVENANCE['caches'] = [digest(path) for path in files+[clean, regrid]]
    PROVENANCE['data_sources'] = [digest(path) for path in
                                  [italy.CELLS_MAT, italy.POLYGON_MAT, italy.RAW_CATALOG]
                                  if path.is_file()]
    testing = {}
    for path in files:
        array = np.load(path, allow_pickle=False)
        period = path.stem.split('_')[-1]
        shape = (88 if period == 'learning' else 40, 177, 25)
        require(f'{path.stem} shape/nonnegative/finite', array.shape == shape
                and np.isfinite(array).all() and (array >= 0).all(), shape=list(array.shape), total=float(array.sum()))
        if period == 'testing':
            testing[path.stem.split('_')[0]] = array
    cells = italy.testing_cells()
    require('177 positive projected areas', len(cells) == 177 and (cells.area_km2 > 0).all())
    # Use corners from projected boundaries rather than stored geographic centers.
    x, y = cells.x0.to_numpy(), cells.y1.to_numpy()
    lon, lat = italy.km_to_lonlat(x, y)
    xx, yy = italy.lonlat_to_km(lon, lat)
    check('EPSG:7794 km geographic round-trip', np.max(np.hypot(xx-x, yy-y)), 0., 1e-5)
    check('177 equal cell areas', np.max(abs(cells.area_km2.to_numpy()-1800.)), 0., 1e-7)
    cat = italy.experiment_catalog()
    learn, target = italy.target_events(cat, 'learning'), italy.target_events(cat, 'testing')
    check('HORUS 2024 learning targets', len(learn), 27)
    check('HORUS 2024 testing targets', len(target), 25)
    learning_windows = models.learning_windows()
    check('learning final boundary is 2012-01-01', learning_windows.t2.iloc[-1], italy.year_start_days(2012))
    check('learning windows cover exact learning period', (learning_windows.t2-learning_windows.t1).sum(),
          italy.year_start_days(2012)-italy.year_start_days(1990))
    windows = italy.forecast_windows()
    check('forecast window count', len(windows), 40)
    check('forecast total duration', (windows.t2-windows.t1).sum(), 3652.4, 1e-8)
    tail = italy.year_start_days(2022)-windows.t2.iloc[-1]
    check('uncovered testing-period tail in days', tail, .6, 1e-8)
    check('target count in uncovered tail', ((target.t_days >= windows.t2.iloc[-1])
          & (target.t_days < italy.year_start_days(2022))).sum(), 0)
    omega = italy.bin_targets(target, windows)
    check('all real testing targets retained by bins', omega.sum(), len(target))
    # Independent SUP total, including finite magnitude range and exact duration.
    duration = italy.year_start_days(2012)-italy.year_start_days(1990)
    mass = 1-10**(-italy.SPEC.b_value*(italy.SPEC.m_max-italy.SPEC.mT))
    total = len(learn)/duration*(windows.t2-windows.t1).sum()*mass
    check('cached SUP total vs independent duration and GR mass', testing['SUP'].sum(), total)
    W, grid = italy.regrid_matrix_01deg()
    require('regrid weights finite/nonnegative', np.isfinite(W).all() and (W >= 0).all())
    check('regrid every source column sums to one', np.max(abs(W.sum(axis=0)-1)), 0., 1e-12)
    native = testing['EEPAS'].sum(axis=0)
    mapped = W @ native
    check('regrid preserves every magnitude-bin total', np.max(abs(mapped.sum(axis=0)-native.sum(axis=0))), 0., 1e-10)
    # Verify actual writer roundoff, not only an independent format-string copy.
    for label, writer, expected_rows in [
        ('native', italy.write_csep_10col, 177*25),
        ('regridded', italy.write_csep_10col_regridded, len(grid)*25),
    ]:
        path = OUT/f'math_export_{label}.tsv'
        writer(native, path)
        exported = pd.read_csv(path, sep='\t')
        require(f'{label} CSEP export shape/flags', exported.shape == (expected_rows, 10)
                and (exported.FLAG == 1).all())
        # %.6e retains seven significant digits; the total error is bounded by
        # 5e-7 of the nonnegative total, with a small floating-point margin.
        check(f'{label} CSEP export roundoff total', exported.RATE.sum(), native.sum(),
              max(1e-10, native.sum()*5.1e-7))
        PROVENANCE.setdefault('verification_exports', []).append(digest(path))
    PROVENANCE['limitations'] = [
        'Forecast caches were read, not regenerated; this is not a complete independent model reconstruction.',
        'Regrid column normalization proves mass conservation, not geometric intersection accuracy.',
        'Midpoint magnitude accuracy is checked on one explicit EEPAS event, not every source event.',
        'Source and cache hashes identify this run; historical cache-generation code was not inferred from hashes.',
    ]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    PROVENANCE.update(
        seed=SEED, maximum_score_simulation_catalogues=8,
        spec=asdict(italy.SPEC), published_model_parameters=models.PARAMS,
        paper='Biondini, Rhoades & Gasperini (2023), doi:10.1093/gji/ggad123, Table 3 and Appendix A',
        cache_rebuilt=False,
    )
    source_paths = [ROOT/'scripts/verify_teaching_math.py', ROOT/'gdms_toolkit/csep_teaching.py',
                    ROOT/'gdms_toolkit/italy.py', ROOT/'gdms_toolkit/italy_models.py',
                    ROOT/'book/11_conditional_intensity.py',
                    ROOT/'book/appendix_a_point_process.md',
                    ROOT/'book/appendix_d_eepas.md',
                    ROOT/'book/appendix_e_testing.md',
                    ROOT/'book/appendix_f_hazard.md']
    paper_text = Path('/tmp/italy_eepas_source.txt')
    if paper_text.is_file():
        source_paths.append(paper_text)
    PROVENANCE['sources'] = [digest(path) for path in source_paths]
    error = None
    try:
        legacy_checks()
        score_checks()
        kernel_checks()
        boundary_checks()
        model_concept_checks()
        cache_checks()
        from unittest.mock import patch
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as temporary:
            absent=Path(temporary)/'EEPAS_testing.npy'
            with patch.object(models,'_cache_path',return_value=absent), patch.object(models,'eepas_forecast',side_effect=AssertionError('Unexpected model recomputation')):
                try:
                    models.get_forecast('EEPAS')
                except FileNotFoundError:
                    require('missing forecast cache never triggers implicit recomputation',True)
                else:
                    require('missing forecast cache never triggers implicit recomputation',False)
    except Exception as exc:
        error = dict(type=type(exc).__name__, message=str(exc))
        raise
    finally:
        report = dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
                      passed=error is None and all(row['passed'] for row in ROWS),
                      checks=ROWS, provenance=PROVENANCE, error=error)
        path = OUT/'math_checks.json'
        path.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
        print(json.dumps(dict(report=str(path), passed=report['passed'], checks=len(ROWS), error=error),
                         ensure_ascii=False))


if __name__ == '__main__':
    main()
