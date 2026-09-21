#!/usr/bin/env python3
"""
multi_instance_oracle_v3.py — external review round 7 (author-approved
full scope): a parameterised version of bdcz_oracle_v3's decoder/oracle
for city subsets OTHER than the main 100-city benchmark, to test whether
the algorithm ranking, medium-necessity finding, and effective
dimensionality are stable across problem instances (a recurring "single
instance" objection across three rounds of review).

Reuses bdcz_oracle_v3's physical constants and pure per-link formulas
(fiber_link, satellite_link, fso_link, _werner_compose -- these depend
only on dist_km, not on which/how-many cities exist) directly, so the
physics is identical to the main benchmark's; only the city set (and
therefore the distance matrix, region labels, and haversine-derived
adjacency) differs. decode_topology/oracle_evaluate/eval_single are
reimplemented here as closures over a per-instance (DIST, city_region,
LINK_FID, LINK_RATE, LINK_LAT) tuple, rather than module globals as in
bdcz_oracle_v3.py, since a single process needs several DIFFERENT
instances simultaneously.

New instances (all real cities -- verified subsets/regions of the
existing, already-real 100-city list, no fabricated coordinates):
  - global_n50: every other city (indices 0,2,4,...) from the main
    100-city list, preserving global geographic spread at half the size.
  - regional_eu_na: the 30 Europe + 20 North America cities only (two
    regions, cross-continental but not globally dispersed) -- 50 cities.
  - metro_europe: the 30 Europe-region cities only -- a genuinely denser,
    single-region instance (shorter typical distances) than the main
    global-sparse benchmark.

Round 8: a genuinely independent N=200 instance (`global_n200`) is now
included, sourced from GeoNames' public `cities15000` export (real
cities, real coordinates, verified -- see cities.csv rows 100-199 and
their `coordinate_source` column for the exact provenance of each of
the 100 new cities; the original 100 are unchanged). This was
previously left out on the grounds that fabricating coordinates to
reach N=200 would violate this project's no-fabrication standard --
the resolution is sourcing genuinely real additional cities, not
relaxing that standard.
"""
import csv
import heapq
from pathlib import Path

import numpy as np

from bdcz_oracle_v3 import (
    fiber_link, satellite_link, fso_link, _werner_compose,
    XL, XU, N_VAR, dedup_front, city_lat, city_lon, city_region as _MAIN_REGION,
)

CITIES_CSV = Path(__file__).resolve().parent.parent / 'cities.csv'
_REGION_NAME_TO_CODE = {'Europe': 0, 'North America': 1, 'Asia-Pacific': 2,
                        'Middle East/Africa': 3, 'South America': 4}


def _load_extended_200_city_table():
    """Reads cities.csv (200 rows: the original 100 + 100 new real
    cities sourced from GeoNames, Round 8) into lat/lon/region arrays,
    independent of bdcz_oracle_v3's fixed 100-city module globals."""
    with open(CITIES_CSV) as f:
        rows = list(csv.DictReader(f))
    rows.sort(key=lambda r: int(r['index']))
    lat = np.array([float(r['latitude']) for r in rows])
    lon = np.array([float(r['longitude']) for r in rows])
    region = np.array([_REGION_NAME_TO_CODE[r['region']] for r in rows], dtype=np.int64)
    return lat, lon, region


_city_lat_200, _city_lon_200, _city_region_200 = _load_extended_200_city_table()

assert np.allclose(_city_lat_200[:100], city_lat), \
    'cities.csv rows 0-99 latitude must match bdcz_oracle_v3.city_lat exactly'
assert np.allclose(_city_lon_200[:100], city_lon), \
    'cities.csv rows 0-99 longitude must match bdcz_oracle_v3.city_lon exactly'
assert np.array_equal(_city_region_200[:100], _MAIN_REGION), \
    'cities.csv rows 0-99 region must match bdcz_oracle_v3.city_region exactly'


def haversine(a, b, c, d):
    R = 6371.0
    p1, p2 = np.radians(a), np.radians(c)
    dp, dl = np.radians(c - a), np.radians(d - b)
    x = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return 2 * R * np.arcsin(np.sqrt(np.clip(x, 0, 1)))


INSTANCES = {
    'global_n50': {'indices': list(range(0, 100, 2)), 'description': '50 cities, every other of the main 100 (preserves global spread)'},
    'regional_eu_na': {'indices': list(range(0, 50)), 'description': '50 cities: Europe (30) + North America (20) only'},
    'metro_europe': {'indices': list(range(0, 30)), 'description': '30 cities: Europe region only (denser, single-region)'},
    'global_n200': {'indices': list(range(0, 200)), 'description': '200 cities: the original 100 plus 100 additional real '
                    'cities sourced from GeoNames cities15000 (Round 8), a genuinely independent scale-up rather than a '
                    'subset of the main benchmark', 'use_200_table': True},
}


class Instance:
    def __init__(self, name):
        cfg = INSTANCES[name]
        idx = np.array(cfg['indices'])
        self.name = name
        self.description = cfg['description']
        self.n_nodes = len(idx)
        if cfg.get('use_200_table'):
            self.city_lat = _city_lat_200[idx]
            self.city_lon = _city_lon_200[idx]
            self.city_region = _city_region_200[idx]
        else:
            self.city_lat = city_lat[idx]
            self.city_lon = city_lon[idx]
            self.city_region = _MAIN_REGION[idx]
        self.dist = np.zeros((self.n_nodes, self.n_nodes))
        for i in range(self.n_nodes):
            self.dist[i] = haversine(self.city_lat[i], self.city_lon[i], self.city_lat, self.city_lon)

        self.link_fid = np.zeros((self.n_nodes, self.n_nodes, 3))
        self.link_rate = np.zeros((self.n_nodes, self.n_nodes, 3))
        self.link_lat = np.zeros((self.n_nodes, self.n_nodes, 3))
        for i in range(self.n_nodes):
            for j in range(i + 1, self.n_nodes):
                for t, fn in enumerate([fiber_link, satellite_link, fso_link]):
                    f, r, l = fn(self.dist[i, j])
                    self.link_fid[i, j, t] = self.link_fid[j, i, t] = f
                    self.link_rate[i, j, t] = self.link_rate[j, i, t] = r
                    self.link_lat[i, j, t] = self.link_lat[j, i, t] = l

    def _build_adjacency(self, eu, ev, et):
        best = {}
        for k in range(len(eu)):
            u, v, t = int(eu[k]), int(ev[k]), int(et[k])
            f = self.link_fid[u, v, t]
            if f <= 0.0:
                continue
            key = (min(u, v), max(u, v))
            if key not in best or f > best[key][0]:
                best[key] = (f, max(self.link_rate[u, v, t], 1e-3), self.link_lat[u, v, t])
        adj = [[] for _ in range(self.n_nodes)]
        for (u, v), (f, r, l) in best.items():
            adj[u].append((v, f, r, l))
            adj[v].append((u, f, r, l))
        return adj

    def _widest_paths(self, adj, source):
        N = self.n_nodes
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

    def oracle_evaluate(self, eu, ev, et):
        N = self.n_nodes
        FID_THRESHOLD = 0.70
        total_pairs = N * (N - 1) // 2
        adj = self._build_adjacency(eu, ev, et)
        adj_lookup = {(u, v): (r, l) for u in range(N) for v, f, r, l in adj[u]}

        n_pairs = 0
        sum_fid = sum_rate_log = sum_lat = 0.0
        n_covered = 0
        for src in range(N):
            best_F, pred = self._widest_paths(adj, src)
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
                hop_rates, hop_lats = [], []
                for a, b in zip(path_nodes[:-1], path_nodes[1:]):
                    r, l = adj_lookup.get((a, b), adj_lookup.get((b, a), (1e-3, 1e18)))
                    hop_rates.append(r)
                    hop_lats.append(l)
                n_hops = len(hop_rates)
                bottleneck_rate = min(hop_rates) if hop_rates else 1e-3
                rate_e2e = bottleneck_rate * (0.50 ** max(n_hops - 1, 0))
                lat_e2e = sum(hop_lats)
                sum_fid += F_e2e
                sum_rate_log += np.log1p(rate_e2e)
                sum_lat += lat_e2e
                if F_e2e >= FID_THRESHOLD:
                    n_covered += 1
        if n_pairs == 0:
            return 0.5, 0.0, 1.0, 0.0
        return sum_fid / n_pairs, sum_rate_log / n_pairs, sum_lat / n_pairs, n_covered / total_pairs

    def decode_topology(self, params, rng=None, enable_fibre=True, enable_satellite=True, enable_fso=True):
        N = self.n_nodes
        k_fibre = max(1, int(round(params[0])))
        k_sat = max(1, int(round(params[1])))
        k_fso = max(0, int(round(params[2])))
        max_d_fibre, max_d_sat, max_d_fso = params[3], params[4], params[5]
        p_fibre = np.clip(params[6], 0.05, 1.0)
        p_sat = np.clip(params[7], 0.05, 1.0)
        p_fso = np.clip(params[8], 0.05, 1.0)
        interregion_bonus = np.clip(params[9], 0.0, 1.0)
        seed_offset = int(round(params[10])) % 1000
        if rng is None:
            rng = np.random.default_rng(seed_offset)

        edges = set()
        for u in range(N):
            d_row = self.dist[u]
            if enable_fibre:
                cands = sorted((d_row[v], v) for v in range(N) if u != v and d_row[v] <= max_d_fibre)
                for _, v in cands[:k_fibre]:
                    p = p_fibre
                    if self.city_region[u] != self.city_region[v]:
                        p *= (1 + interregion_bonus)
                    if rng.random() < min(p, 1.0):
                        edges.add((min(u, v), max(u, v), 0))
            if enable_satellite:
                cands = sorted((d_row[v], v) for v in range(N) if u != v and 500 <= d_row[v] <= max_d_sat)
                for _, v in cands[:k_sat]:
                    p = p_sat
                    if self.city_region[u] != self.city_region[v]:
                        p = min(1.0, p * (1 + interregion_bonus))
                    if rng.random() < p:
                        edges.add((min(u, v), max(u, v), 1))
            if enable_fso and k_fso > 0:
                cands = sorted((d_row[v], v) for v in range(N) if u != v and d_row[v] <= max_d_fso)
                for _, v in cands[:k_fso]:
                    if rng.random() < p_fso:
                        edges.add((min(u, v), max(u, v), 2))
        edges = sorted(edges)
        if len(edges) < 5:
            return None
        eu = np.array([e[0] for e in edges], dtype=np.int64)
        ev = np.array([e[1] for e in edges], dtype=np.int64)
        et = np.array([e[2] for e in edges], dtype=np.int64)
        return eu, ev, et

    def eval_single(self, params, enable_fibre=True, enable_satellite=True, enable_fso=True):
        result = self.decode_topology(params, enable_fibre=enable_fibre,
                                       enable_satellite=enable_satellite, enable_fso=enable_fso)
        if result is None:
            return np.array([0.5, 0.0, 1.0, 0.0])
        eu, ev, et = result
        return np.array(self.oracle_evaluate(eu, ev, et))
