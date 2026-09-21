#!/usr/bin/env python3
"""
bdcz_oracle_v3.py — Canonical BDCZ oracle, unified with Paper 1's final oracle
=================================================================================
This module replaces bdcz_oracle_v2.py. It is derived directly from
Paper 1's canonical, submitted oracle
(Paper1_REVIEW_PACKAGE_v13/01_code_and_data/bdcz_oracle.py::oracle_evaluate),
which is the exact function used to generate QI-Bench-100 -- not the
older -log(F)-product/floyd_warshall proxy with an undisclosed fidelity
floor that an earlier Paper 2 draft had accidentally used instead (that
proxy exists in Paper 1's own codebase as oracle_evaluate_proxy, kept
there only for a speed-vs-accuracy comparison, and was never the oracle
that produced any reported Paper 1 or Paper 2 number).

Two differences from bdcz_oracle_v2.py, both correcting real
discrepancies against Paper 1's canonical physics:
  1. SAT_BASE_FID: 0.97 -> 0.83 (Micius-calibrated, Yin et al. 2017 /
     Khatri et al. 2021), with an added distance-dependent fidelity
     decay (SAT_FID_DECOHERENCE) that the earlier draft omitted entirely.
  2. FSO_BASE_FID: 0.98 -> 0.92 (matches the literature-cited value used
     throughout Paper 1).
  3. Path selection: exact Werner-optimal widest-path search
     (_werner_optimal_paths_from_source), not the -log(F)-product proxy
     with cost_fid = max(LINK_FID, 0.5) floor.

hub_boost remains removed (11 free parameters, not 12) for the reason
already documented in the paper: it was a dead parameter in every draft
of the decoder, never referenced in decode_topology.
"""

import heapq

import numpy as np

# ── Cities ──────────────────────────────────────────────────────────────────
CITIES = [
    (40.4168,-3.7038),(41.3851,2.1734),(38.7223,-9.1393),(48.8566,2.3522),
    (51.5074,-0.1278),(52.3676,4.9041),(50.8503,4.3517),(52.52,13.405),
    (48.1351,11.582),(53.5511,9.9937),(50.1109,8.6821),(47.3769,8.5417),
    (46.2044,6.1432),(48.2082,16.3738),(50.0755,14.4378),(52.2297,21.0122),
    (55.6761,12.5683),(59.3293,18.0686),(59.9139,10.7522),(60.1699,24.9384),
    (53.3498,-6.2603),(41.9028,12.4964),(45.4642,9.19),(37.9838,23.7275),
    (47.4979,19.0402),(44.4268,26.1025),(42.6977,23.3219),(55.9533,-3.1883),
    (52.0116,4.3571),(47.2692,11.4041),
    (40.7128,-74.006),(42.3601,-71.0589),(38.9072,-77.0369),(41.8781,-87.6298),
    (33.749,-84.388),(25.7617,-80.1918),(29.7604,-95.3698),(32.7767,-96.797),
    (39.7392,-104.9903),(34.0522,-118.2437),(37.7749,-122.4194),(47.6062,-122.3321),
    (49.2827,-123.1207),(43.6532,-79.3832),(45.5017,-73.5673),(19.4326,-99.1332),
    (35.9606,-84.2073),(41.712,-87.9799),(34.1478,-118.1445),(40.3573,-74.6672),
    (39.9042,116.4074),(31.2304,121.4737),(31.8206,117.2272),(36.6512,117.1201),
    (22.3193,114.1694),(22.5431,114.0579),(23.1291,113.2644),(30.5728,104.0668),
    (30.5928,114.3055),(34.3416,108.9398),(35.6762,139.6503),(34.6937,135.5023),
    (37.5665,126.978),(25.033,121.5654),(1.3521,103.8198),(13.7563,100.5018),
    (3.139,101.6869),(-6.2088,106.8456),(14.5995,120.9842),(21.0278,105.8342),
    (19.076,72.8777),(28.7041,77.1025),(12.9716,77.5946),(-33.8688,151.2093),
    (-37.8136,144.9631),
    (41.0082,28.9784),(25.2048,55.2708),(24.4539,54.3773),(25.2854,51.531),
    (24.7136,46.6753),(32.0853,34.7818),(31.7683,35.2137),(30.0444,31.2357),
    (35.6892,51.389),(6.5244,3.3792),(-1.2921,36.8219),(-26.2041,28.0473),
    (-33.9249,18.4241),(33.5731,-7.5898),(9.145,40.4897),
    (-23.5505,-46.6333),(-22.9068,-43.1729),(-34.6037,-58.3816),(-33.4489,-70.6693),
    (-12.0464,-77.0428),(4.711,-74.0721),(10.4806,-66.9036),(-0.1807,-78.4678),
    (-15.8267,-47.9218),(-34.9011,-56.1745),
]

city_lat = np.array([c[0] for c in CITIES], dtype=np.float64)
city_lon = np.array([c[1] for c in CITIES], dtype=np.float64)
city_region = np.array([0]*30 + [1]*20 + [2]*25 + [3]*15 + [4]*10, dtype=np.int64)
N_NODES = 100


def haversine(a, b, c, d):
    """Haversine distance in km between (lat a, lon b) and (lat c, lon d)."""
    R = 6371.0
    p1, p2 = np.radians(a), np.radians(c)
    dp, dl = np.radians(c - a), np.radians(d - b)
    x = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return 2 * R * np.arcsin(np.sqrt(np.clip(x, 0, 1)))


DIST = np.zeros((N_NODES, N_NODES), dtype=np.float64)
for _i in range(N_NODES):
    DIST[_i] = haversine(city_lat[_i], city_lon[_i], city_lat, city_lon)


# ── BDCZ Oracle Physical Parameters (canonical, matches Paper 1) ───────────
FIBER_ATTENUATION       = 0.20
FIBER_BASE_FID          = 0.99
FIBER_DECOHERENCE       = 1.5e-4
FIBER_DET_EFF           = 0.85
FIBER_SOURCE_RATE       = 1e7
FIBER_REPEATER_SPACING  = 50
FIBER_SWAP_SUCCESS      = 0.50

EARTH_RADIUS_KM         = 6371.0      # same value haversine() uses
SAT_ALTITUDE            = 500
SAT_TX_APERTURE         = 0.30
SAT_RX_APERTURE         = 1.00
SAT_WAVELENGTH          = 810e-9
SAT_ATM_TRANSMITTANCE   = 0.80
SAT_POINTING_LOSS_DB    = 3.0
SAT_BASE_FID            = 0.83        # Micius-calibrated (Yin et al. 2017 / Khatri et al. 2021)
SAT_FID_DECOHERENCE     = 5e-6        # 1/km, fidelity loss with slant range
SAT_DET_EFF             = 0.50
SAT_SOURCE_RATE         = 5.9e6
SAT_DUTY_CYCLE          = 0.005

FSO_BASE_FID            = 0.92        # matches Paper 1's literature-cited value
FSO_EXTINCTION          = 0.50
FSO_DET_EFF             = 0.70
FSO_SOURCE_RATE         = 1e7
FSO_MAX_DIST            = 150  # physical hard limit -- decoder search bounds (XU[5]) must not exceed this
# FIX (external review round 6, item 5): this fidelity-decoherence rate was
# a bare literal (1e-4) inside oracle_evaluate_fso below, undocumented and
# missing from the manuscript's physical-parameter table (Table 1) even
# though every other medium's decoherence rate is a named constant listed
# there. Named and added to the table.
FSO_FID_DECOHERENCE     = 1e-4  # 1/km, fidelity loss with distance

C_FIBER                 = 2e8


def fiber_link(dist_km):
    n_spans = max(1, int(np.ceil(dist_km / FIBER_REPEATER_SPACING)))
    span_len = dist_km / n_spans
    loss_db = FIBER_ATTENUATION * span_len
    eta_span = 10**(-loss_db / 10) * FIBER_DET_EFF
    F_span = FIBER_BASE_FID * np.exp(-FIBER_DECOHERENCE * span_len)
    F_e2e = F_span
    for _ in range(n_spans - 1):
        F_e2e = F_e2e * F_span + (1 - F_e2e) * (1 - F_span) / 3
    rate_span = FIBER_SOURCE_RATE * eta_span
    rate = rate_span * (FIBER_SWAP_SUCCESS ** (n_spans - 1))
    # FIX (physics audit, PHYSICS_ORACLE_AUDIT.md item 1d, confirmed real):
    # was `dist_km*1000/C_FIBER*(2*n_spans)`, which -- since
    # span_len = dist_km/n_spans -- multiplies the ALREADY-total distance
    # by n_spans again, giving latency quadratic in distance for fixed
    # repeater spacing. Neither the naive-sequential model (round-trip
    # propagation summed over spans collapses to 2*dist_km/v, independent
    # of n_spans once total distance is fixed) nor the BDCZ nested-repeater
    # model (Sangouard, Simon, de Riedmatten, Gisin, Rev. Mod. Phys. 83, 33
    # (2011), Eq. 12: the propagation-only elementary period L0/c DECREASES
    # with more nesting levels at fixed total distance) supports latency
    # increasing with n_spans. This also matches the manuscript's own
    # stated definition ("latency is the sum of per-link propagation
    # times", Introduction) exactly once the erroneous n_spans factor is
    # removed.
    latency_s = dist_km * 1000 / C_FIBER * 2
    return F_e2e, rate, latency_s


# Maximum ground separation of two stations that can see the SAME
# satellite simultaneously, at zero elevation angle (horizon-grazing):
# each station subtends a central half-angle arccos(R/(R+h)) from the
# sub-satellite point, so the great-circle separation cannot exceed
# twice that arc.
SAT_HORIZON_HALF_ANGLE = np.arccos(
    EARTH_RADIUS_KM / (EARTH_RADIUS_KM + SAT_ALTITUDE))          # rad
SAT_MAX_GROUND_SEPARATION = 2 * EARTH_RADIUS_KM * SAT_HORIZON_HALF_ANGLE
# = 4891.0 km at h = 500 km


def satellite_link(dist_km):
    """Single-satellite downlink pair between two ground stations d km apart.

    FIX (external review, Round 12, confirmed real and quantified): the
    previous implementation used the flat-Earth slant range
    sqrt((d/2)^2 + h^2) and imposed NO visibility limit, while the decoder
    was allowed to request satellite links out to max_d_sat = 15,000 km.
    Beyond a ground separation of SAT_MAX_GROUND_SEPARATION the Earth
    itself blocks the line of sight from a single satellite to both
    stations, so those links cannot exist. Measured on the previously
    released front, 14.6% of all satellite edges (69,800 of 476,827,
    spanning 382 distinct city pairs) violated this constraint. Three
    corrections are applied here:

      (1) Spherical slant range. With half-angle phi = d/(2R), the
          station-to-satellite distance is
              s = sqrt(R^2 + (R+h)^2 - 2 R (R+h) cos phi),
          which reduces to the old expression only in the small-phi
          (flat-Earth) limit and is otherwise strictly larger.
      (2) Hard visibility limit. If phi exceeds arccos(R/(R+h)) the link
          is infeasible and is reported as such, exactly as fso_link()
          already does beyond FSO_MAX_DIST. The decoder's upper bound on
          max_d_sat is capped to match (see XU below), so the search
          cannot propose links the model treats as impossible.
      (3) Consistent two-arm efficiency. The source sits on the satellite
          and sends one photon down each arm; both must be detected for a
          pair to be heralded, so the coincidence efficiency is the
          product over the two arms. The previous expression squared only
          the atmospheric transmittance while applying the geometric,
          pointing and detector terms once, which is internally
          inconsistent regardless of whether those terms are read as
          per-arm or aggregate. Here every per-arm factor is squared.

    Latency is likewise corrected: the two downlinks are simultaneous and
    the arms are equal by the midpoint geometry, so the pair is delivered
    after ONE transit time s/c, not the sum 2s/c of both arms' transits.
    The classical heralding channel between the two ground stations, over
    terrestrial fibre, is unchanged.
    """
    phi = dist_km / (2.0 * EARTH_RADIUS_KM)            # central half-angle, rad
    if phi > SAT_HORIZON_HALF_ANGLE:
        return 0.0, 0.0, 1e6                            # no line of sight

    Rs = EARTH_RADIUS_KM + SAT_ALTITUDE
    slant_km = np.sqrt(EARTH_RADIUS_KM**2 + Rs**2
                       - 2.0 * EARTH_RADIUS_KM * Rs * np.cos(phi))
    slant_range = slant_km * 1000.0                     # m

    theta_div = SAT_WAVELENGTH / (np.pi * SAT_TX_APERTURE)
    spot_size = theta_div * slant_range
    eta_geo = min(1.0, (SAT_RX_APERTURE / spot_size)**2)
    pointing_loss = 10**(-SAT_POINTING_LOSS_DB / 10)
    eta_arm = eta_geo * SAT_ATM_TRANSMITTANCE * pointing_loss * SAT_DET_EFF
    eta_total = eta_arm**2                              # both arms must fire

    F_e2e = SAT_BASE_FID * np.exp(-SAT_FID_DECOHERENCE * dist_km)
    rate = SAT_SOURCE_RATE * eta_total * SAT_DUTY_CYCLE
    #   (1) slant_range/c: one downlink transit. Both photons are emitted
    #       simultaneously and the slant ranges are equal by the
    #       midpoint-satellite geometry, so the pair is delivered after a
    #       single transit time, not the sum of both arms'.
    #   (2) dist_km*1000/C_FIBER: the classical confirmation channel
    #       between the two ground stations, needed to correlate detection
    #       events and herald success -- assumed routed over the same
    #       terrestrial fibre network this benchmark already assumes
    #       exists between any two cities, rather than relayed back
    #       through the satellite.
    latency_s = slant_range / 3e8 + dist_km * 1000 / C_FIBER
    return F_e2e, rate, latency_s


def fso_link(dist_km):
    if dist_km > FSO_MAX_DIST:
        return 0.0, 0.0, 1e6
    loss_db = FSO_EXTINCTION * dist_km
    eta = 10**(-loss_db / 10) * FSO_DET_EFF
    F_e2e = FSO_BASE_FID * np.exp(-FSO_FID_DECOHERENCE * dist_km)
    rate = FSO_SOURCE_RATE * eta
    latency_s = dist_km * 1000 / 3e8 * 2
    return F_e2e, rate, latency_s


LINK_FID  = np.zeros((N_NODES, N_NODES, 3), dtype=np.float64)
LINK_RATE = np.zeros((N_NODES, N_NODES, 3), dtype=np.float64)
LINK_LAT  = np.zeros((N_NODES, N_NODES, 3), dtype=np.float64)

for _i in range(N_NODES):
    for _j in range(_i + 1, N_NODES):
        for _t, _fn in enumerate([fiber_link, satellite_link, fso_link]):
            _f, _r, _l = _fn(DIST[_i, _j])
            LINK_FID[_i, _j, _t] = LINK_FID[_j, _i, _t] = _f
            LINK_RATE[_i, _j, _t] = LINK_RATE[_j, _i, _t] = _r
            LINK_LAT[_i, _j, _t] = LINK_LAT[_j, _i, _t] = _l

BEST_MEDIUM = np.argmax(LINK_FID, axis=2).astype(np.int64)


# ── Variable bounds: 11 free parameters (hub_boost removed) ────────────────
#   0  k_fibre           [3,  25]
#   1  k_sat             [2,  25]
#   2  k_fso             [0,  10]
#   3  max_d_fibre  km   [300, 3000]
#   4  max_d_sat    km   [1000, 15000]
#   5  max_d_fso    km   [10,  150]   -- FIX (Round 3 review): was [10,200],
#      allowing the optimiser to search FSO distances the oracle itself
#      treats as physically invalid (fso_link returns 0 fidelity/rate above
#      FSO_MAX_DIST=150km). Bounded to match the physical model exactly.
#   6  p_fibre           [0.10, 1.0]
#   7  p_sat             [0.10, 1.0]
#   8  p_fso             [0.10, 1.0]
#   9  interregion_bonus [0.0, 1.0]
#  10  seed_offset       [0,   999]
# FIX (external review, Round 12): max_d_sat's upper bound was 15,000 km,
# far beyond the ground separation at which a single satellite can see both
# stations (SAT_MAX_GROUND_SEPARATION = 4891 km at h = 500 km), so the
# search could and did propose satellite links with no line of sight. The
# bound is now capped to the physical limit, exactly as max_d_fso was
# bounded to FSO_MAX_DIST in Round 3 for the same reason: the search space
# must not contain points the physical model treats as impossible.
XL = np.array([ 3,  2,  0,  300,  1000,  10, 0.10, 0.10, 0.10, 0.0,   0],
               dtype=np.float64)
XU = np.array([25, 25, 10, 3000, SAT_MAX_GROUND_SEPARATION, 150,
               1.00, 1.00, 1.00, 1.0, 999], dtype=np.float64)
N_VAR = 11


def decode_topology(params, rng=None, enable_fibre=True, enable_satellite=True,
                     enable_fso=True):
    """Decode an 11-element parameter vector into a network topology.
    Returns (eu, ev, et) or None if fewer than 5 edges result.

    enable_fibre/enable_satellite/enable_fso: FIX (Round 3 review, confirmed
    real bug) -- the ablation study needs to genuinely disable a medium
    (e.g. 'satellite-only' must produce ZERO fibre edges), but decode_topology
    internally floors k_fibre/k_sat at 1 (`max(1, ...)`) and p_fibre/p_sat/
    p_fso at 0.05 (`np.clip(..., 0.05, 1.0)`) for the Full Decoder's own
    legitimate search-space reasons (a fibre/satellite hub count or edge
    probability of exactly 0 is not a meaningful point in the Full Decoder's
    continuous search space). Setting an ablation variant's params[0]/[6]
    etc. to 0.0 therefore does NOT actually zero out that medium -- it still
    floors to k=1, p=0.05, so "satellite-only" could still emit fibre edges.
    These three flags provide a genuine, unambiguous hard disable that does
    not depend on what value the (now-irrelevant) decoder parameters take,
    without changing the Full Decoder's own bounds/behaviour when all three
    default to True (verified: the flags do not alter any code path when all
    are True, so every pre-existing caller -- main comparison, classical
    baselines, binary baseline -- is unaffected by this change)."""
    k_fibre          = max(1, int(round(params[0])))
    k_sat            = max(1, int(round(params[1])))
    k_fso            = max(0, int(round(params[2])))
    max_d_fibre      = params[3]
    max_d_sat        = params[4]
    max_d_fso        = params[5]
    p_fibre          = np.clip(params[6], 0.05, 1.0)
    p_sat            = np.clip(params[7], 0.05, 1.0)
    p_fso            = np.clip(params[8], 0.05, 1.0)
    interregion_bonus= np.clip(params[9], 0.0, 1.0)
    seed_offset      = int(round(params[10])) % 1000

    if rng is None:
        rng = np.random.default_rng(seed_offset)

    edges = set()
    for u in range(N_NODES):
        d_row = DIST[u]

        if enable_fibre:
            fibre_cands = [(d_row[v], v) for v in range(N_NODES)
                           if u != v and d_row[v] <= max_d_fibre]
            fibre_cands.sort()
            for _, v in fibre_cands[:k_fibre]:
                p = p_fibre
                if city_region[u] != city_region[v]:
                    p *= (1 + interregion_bonus)
                if rng.random() < min(p, 1.0):
                    edges.add((min(u, v), max(u, v), 0))

        if enable_satellite:
            sat_cands = [(d_row[v], v) for v in range(N_NODES)
                         if u != v and 500 <= d_row[v] <= max_d_sat]
            sat_cands.sort()
            for _, v in sat_cands[:k_sat]:
                p = p_sat
                if city_region[u] != city_region[v]:
                    p = min(1.0, p * (1 + interregion_bonus))
                if rng.random() < p:
                    edges.add((min(u, v), max(u, v), 1))

        if enable_fso and k_fso > 0:
            fso_cands = [(d_row[v], v) for v in range(N_NODES)
                         if u != v and d_row[v] <= max_d_fso]
            fso_cands.sort()
            for _, v in fso_cands[:k_fso]:
                if rng.random() < p_fso:
                    edges.add((min(u, v), max(u, v), 2))

    edges = sorted(edges)
    if len(edges) < 5:
        return None

    eu = np.array([e[0] for e in edges], dtype=np.int64)
    ev = np.array([e[1] for e in edges], dtype=np.int64)
    et = np.array([e[2] for e in edges], dtype=np.int64)
    return eu, ev, et


def _werner_compose(f1, f2):
    return f1 * f2 + (1 - f1) * (1 - f2) / 3


def _build_best_medium_adjacency(eu, ev, et):
    """adj[u] = list of (v, F_uv, rate_uv, lat_uv), collapsing parallel
    media between the same city pair to the higher-fidelity one for
    path search (matches Paper 1's canonical oracle exactly).

    FIX (external review round 6, confirmed real): an edge whose chosen
    medium is physically infeasible at that city pair's distance (e.g.\
    FSO assigned beyond FSO_MAX_DIST) has LINK_FID==0 and LINK_LAT set to
    the corresponding `*_link()` function's 1e6 infeasibility sentinel,
    not a real latency. The parametric decoder's search bounds never
    propose such a pair (verified: max latency across all 150 main-run
    checkpoints is 0.146s), but the unconstrained categorical baseline
    (every one of the 4,950 pairs is a free {none,fibre,satellite,fso}
    choice) can and did propose one -- and because the resulting edge's
    fidelity is not used to reject it here, on a very sparse topology the
    widest-path search was occasionally forced to route through it as
    the only available path, leaking the 1e6-second sentinel into a
    real, averaged, reported latency (one confirmed instance: NSGA-III
    categorical seed=48, mean latency 3696.99s, traced to exactly this).
    A medium with zero fidelity is not a usable repeater link, so it is
    excluded from the path-search graph entirely, matching the physical
    model's own semantics."""
    best = {}
    for k in range(len(eu)):
        u, v, t = int(eu[k]), int(ev[k]), int(et[k])
        f = LINK_FID[u, v, t]
        if f <= 0.0:
            continue
        key = (min(u, v), max(u, v))
        if key not in best or f > best[key][0]:
            # FIX (external review, Round 12): this used to floor the
            # per-hop rate at 1e-3 s^-1. The floor was not a numerical
            # guard -- no reported quantity divides by the rate, and every
            # value involved is representable in double precision -- but it
            # silently inflated long-haul fibre rates by up to twelve
            # orders of magnitude (at 3,000 km the model gives 1.5e-12 and
            # the floor reported 1e-3). Removed; the true rate is used.
            best[key] = (f, LINK_RATE[u, v, t], LINK_LAT[u, v, t])
    adj = [[] for _ in range(N_NODES)]
    for (u, v), (f, r, l) in best.items():
        adj[u].append((v, f, r, l))
        adj[v].append((u, f, r, l))
    return adj


def _werner_optimal_paths_from_source(adj, source):
    """Maximum Werner-composed-fidelity path search (Dijkstra-like
    relaxation) -- not a classical widest/bottleneck path (Werner
    composition is neither a min nor a simple product)."""
    N = N_NODES
    best_F = np.zeros(N)
    best_F[source] = 1.0
    pred = -np.ones(N, dtype=np.int64)
    visited = np.zeros(N, dtype=bool)
    heap = [(-1.0, source)]
    while heap:
        neg_f, u = heapq.heappop(heap)
        if visited[u]:
            continue
        visited[u] = True
        f_u = -neg_f
        for v, f_uv, _r, _l in adj[u]:
            if visited[v]:
                continue
            cand = _werner_compose(f_u, f_uv)
            if cand > best_F[v]:
                best_F[v] = cand
                pred[v] = u
                heapq.heappush(heap, (-cand, v))
    return best_F, pred


def oracle_evaluate(eu, ev, et):
    """Evaluate a topology using exact Werner-optimal path selection.
    Returns (avg_fid, avg_rate_log, avg_lat, coverage), where avg_rate_log
    is the mean of log1p(rate) over reachable pairs (not the log of the
    mean rate -- Eq. 2 in the manuscript is written to match this)."""
    N = N_NODES
    FID_THRESHOLD = 0.70
    total_pairs = N * (N - 1) // 2

    adj = _build_best_medium_adjacency(eu, ev, et)
    adj_lookup = {(u, v): (r, l) for u in range(N) for v, f, r, l in adj[u]}

    n_pairs = 0
    sum_fid = 0.0
    sum_rate_log = 0.0
    sum_lat = 0.0
    n_covered = 0

    for src in range(N):
        best_F, pred = _werner_optimal_paths_from_source(adj, src)
        for tgt in range(src + 1, N):
            if best_F[tgt] <= 0.0 or pred[tgt] < 0:
                continue
            n_pairs += 1
            F_e2e = best_F[tgt]

            path_nodes = []
            cur = tgt
            while cur != src and cur >= 0:
                path_nodes.append(cur)
                cur = pred[cur]
            path_nodes.append(src)
            path_nodes.reverse()

            hop_rates = []
            hop_lats = []
            for a, b in zip(path_nodes[:-1], path_nodes[1:]):
                r, l = adj_lookup.get((a, b), adj_lookup.get((b, a), (1e-3, 1e18)))
                hop_rates.append(r)
                hop_lats.append(l)

            n_hops = len(hop_rates)
            bottleneck_rate = min(hop_rates) if hop_rates else 1e-3
            rate_e2e = bottleneck_rate * (FIBER_SWAP_SUCCESS ** max(n_hops - 1, 0))
            lat_e2e = sum(hop_lats)

            sum_fid += F_e2e
            sum_rate_log += np.log1p(rate_e2e)
            sum_lat += lat_e2e
            if F_e2e >= FID_THRESHOLD:
                n_covered += 1

    if n_pairs == 0:
        return 0.5, 0.0, 1.0, 0.0

    return sum_fid/n_pairs, sum_rate_log/n_pairs, sum_lat/n_pairs, n_covered/total_pairs


def eval_single(params, enable_fibre=True, enable_satellite=True, enable_fso=True):
    """Evaluate a single 11-element parameter vector -> (4,) objective array."""
    result = decode_topology(params, enable_fibre=enable_fibre,
                              enable_satellite=enable_satellite, enable_fso=enable_fso)
    if result is None:
        return np.array([0.5, 0.0, 1.0, 0.0])
    eu, ev, et = result
    return np.array(oracle_evaluate(eu, ev, et))


def dedup_front(F, X, decimals=6):
    """Deduplicate a non-dominated front by OBJECTIVE VECTOR (rounded to
    `decimals`), keeping the first occurrence's decision vector.

    FIX (external review round 4, confirmed real and material): a single
    MOEA/D run's own returned non-dominated population routinely contains
    >80% exact-duplicate objective vectors (verified directly: seed=42,
    80 non-dominated individuals, only 11 unique objective vectors) --
    the population converges onto a handful of true optima and keeps
    re-finding them. Computing HV/IGD/spacing/|F*| on the raw,
    duplicate-laden front is wrong for three of those four metrics: HV is
    a set measure and is provably invariant to exact duplicates (adding
    the same point twice does not change the dominated volume), but IGD's
    reference front implicitly overweights duplicate-heavy regions,
    spacing's nearest-neighbour distances are deflated by near-zero
    distances between duplicates (making a degenerately-clustered front
    look artificially "even"), and |F*| as the raw population size
    overstates how many genuinely distinct solutions were found. Call this
    immediately after `find_non_dominated`, before computing any metric or
    persisting a run's front, so every downstream consumer (this script,
    the ablation script, global-front pooling) works from deduplicated
    data by construction rather than needing to remember to dedup later."""
    if len(F) == 0:
        return F, X
    _, idx = np.unique(np.round(F, decimals), axis=0, return_index=True)
    idx = np.sort(idx)
    return F[idx], X[idx]
