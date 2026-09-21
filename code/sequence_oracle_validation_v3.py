#!/usr/bin/env python3
"""
sequence_oracle_validation_v3.py -- Round 7 (SeQUeNCe external oracle validation)
=================================================================================
Scope, honestly stated up front
--------------------------------
This script does NOT attempt to numerically reproduce this paper's analytical
BDCZ oracle (bdcz_oracle_v3.py) inside SeQUeNCe (Argonne's discrete-event
quantum-network simulator, pip package `sequence`, v0.8.5). That would
require either (a) forcing SeQUeNCe's entanglement-generation protocol to
use this paper's exact Werner-fidelity recursion, which its API does not
expose at the elementary-link level (see below), or (b) recalibrating
SeQUeNCe's own physical constants to match ours, which risks silently
introducing a mismatch that *looks* like agreement but isn't -- exactly the
failure mode "Nunca reportar resultados sin verificar" warns against.

Instead, this script runs SeQUeNCe's own, independently-implemented
entanglement-generation protocol (`EntanglementGenerationA`, midpoint-heralded,
the textbook DLCZ/BDCZ elementary-link scheme -- the same physical scheme
this paper's `fiber_link()` models at the single-span level) at several fibre
distances drawn from this paper's own parameter range, and checks that the
*direction* of the distance-dependence agrees: does entanglement-generation
rate fall as distance grows, in both models?

What is and isn't comparable
-----------------------------
Inspecting SeQUeNCe's source (`sequence/components/optical_channel.py`):

    loss = 1 - 10 ** (self.distance * self.attenuation / -10)

is structurally IDENTICAL to this paper's per-span transmittance,

    eta_span = 10**(-loss_db / 10) * FIBER_DET_EFF,  loss_db = FIBER_ATTENUATION * span_len

(both are the standard dB-attenuation-vs-distance law; SeQUeNCe's
`attenuation` is dB/**metre**, this paper's `FIBER_ATTENUATION` is dB/km --
converted below by a factor of 1000). Both models therefore predict
monotonically increasing photon loss, hence monotonically decreasing
elementary-link success rate, as fibre distance increases. That shared
prediction is what this script tests empirically inside SeQUeNCe's own
discrete-event kernel.

Inspecting `sequence/entanglement_management/generation/generation_base.py`,
the fidelity of the state produced by a successful `EntanglementGenerationA`
round is set as `self.fidelity = memory.raw_fidelity` -- a fixed,
distance-independent constant supplied by the memory, not a function of
channel loss or distance. SeQUeNCe's minimal elementary-link protocol
therefore has no built-in mechanism analogous to this paper's
`F_span = FIBER_BASE_FID * exp(-FIBER_DECOHERENCE * span_len)` distance-linked
fidelity decay, short of implementing custom noise/decoherence modelling
that is not part of the standard demo protocol. Consequently: this script
validates the RATE dimension of the analytical oracle against an independent
simulator's own physics; it does NOT independently validate the FIDELITY
dimension, and does not claim to.

Method
------
For a range of single-span fibre distances (1-50 km, i.e. within this
paper's own `FIBER_REPEATER_SPACING = 50` km single-span regime, so no
repeater chaining or entanglement swapping is needed for a fair single-span
comparison), build the exact two-router-plus-BSM-node topology from
SeQUeNCe's official `demo_for_beginners/two_node_eg.ipynb` reference example
(midpoint-heralded entanglement generation between quantum routers `r1`/`r2`
via BSM node `m1`), run the simulation for a fixed wall-clock sim time, and
record how many of the (up to 50, one per quantum memory) entanglement
attempts per node succeed. Repeat for several seeds per distance and report
the median attempted-entanglement count (a direct proxy for rate, since sim
time and memory count are held fixed across distances).

Output: sequence_oracle_validation_v3.json
"""

import argparse
import json
import statistics
import sys
from pathlib import Path

from sequence.kernel.timeline import Timeline
from sequence.topology.node import QuantumRouter, BSMNode
from sequence.components.optical_channel import QuantumChannel, ClassicalChannel
from sequence.resource_management.rule_manager import Rule
from sequence.entanglement_management.generation import EntanglementGenerationA
from sequence.constants import MILLISECOND

HERE = Path(__file__).resolve().parent
OUT_PATH = HERE.parent / "results" / "sequence_oracle_validation_v3.json"

# This paper's canonical fibre attenuation (dB/km) -- see bdcz_oracle_v3.py.
FIBER_ATTENUATION_DB_PER_KM = 0.20
FIBER_ATTENUATION_DB_PER_M = FIBER_ATTENUATION_DB_PER_KM / 1000.0

# Single-span regime only (bdcz_oracle_v3.FIBER_REPEATER_SPACING = 50 km):
# no repeater chaining/swapping needed for a fair single-span comparison.
DISTANCES_KM = [1, 5, 10, 20, 30, 40, 50]
SEEDS = list(range(10))
# Chosen short enough that the entangled-memory count does not saturate at
# N_MEMORIES for any tested distance (a 2000 ms budget saturates at 50/50
# memories for every distance up to 30 km, which flattens the very
# dose-response curve this check exists to show); verified empirically
# before fixing this value.
SIM_TIME_MS = 50
CC_DELAY_MS = 1.0
N_MEMORIES = 50


def eg_rule_condition(memory_info, manager, args):
    if memory_info.state == "RAW":
        return [memory_info]
    return []


def eg_rule_action1(memories_info, args):
    def eg_req_func(protocols, args):
        for protocol in protocols:
            if isinstance(protocol, EntanglementGenerationA):
                return protocol
        return None

    memories = [info.memory for info in memories_info]
    memory = memories[0]
    protocol = EntanglementGenerationA.create(None, "EGA." + memory.name, "m1", "r2", memory)
    protocol.primary = True
    return [protocol, ["r2"], [eg_req_func], [None]]


def eg_rule_action2(memories_info, args):
    memories = [info.memory for info in memories_info]
    memory = memories[0]
    protocol = EntanglementGenerationA.create(None, "EGA." + memory.name, "m1", "r1", memory)
    return [protocol, [None], [None], [None]]


def run_once(dist_km, seed):
    """Run the two-node midpoint-heralded EG demo at a given fibre distance
    and seed; return the number of memories that reached a nonzero
    entangle_time within SIM_TIME_MS (a rate proxy, since sim time and
    memory count are fixed across distances)."""
    PS_PER_MS = 1e9
    M_PER_KM = 1e3

    cc_delay_ps = CC_DELAY_MS * PS_PER_MS
    qc_dist_m = dist_km * M_PER_KM

    tl = Timeline(SIM_TIME_MS * PS_PER_MS)

    r1 = QuantumRouter("r1", tl, N_MEMORIES)
    r2 = QuantumRouter("r2", tl, N_MEMORIES)
    m1 = BSMNode("m1", tl, ["r1", "r2"])

    r1.set_seed(seed)
    r2.set_seed(seed + 1000)
    m1.set_seed(seed + 2000)

    nodes = [r1, r2, m1]
    for node1 in nodes:
        for node2 in nodes:
            if node1 is node2:
                continue
            cc = ClassicalChannel("_".join(["cc", node1.name, node2.name]), tl, 1e3, delay=cc_delay_ps)
            cc.set_ends(node1, node2.name)

    qc1 = QuantumChannel("qc_r1_m1", tl, FIBER_ATTENUATION_DB_PER_M, qc_dist_m)
    qc1.set_ends(r1, m1.name)
    qc2 = QuantumChannel("qc_r2_m1", tl, FIBER_ATTENUATION_DB_PER_M, qc_dist_m)
    qc2.set_ends(r2, m1.name)

    tl.init()
    rule1 = Rule(10, eg_rule_action1, eg_rule_condition, None, None)
    r1.resource_manager.load(rule1)
    rule2 = Rule(10, eg_rule_action2, eg_rule_condition, None, None)
    r2.resource_manager.load(rule2)

    tl.run()

    entangle_times = []
    for info in r1.resource_manager.memory_manager:
        if info.entangle_time > 0:
            entangle_times.append(info.entangle_time / MILLISECOND)
    return len(entangle_times), sorted(entangle_times)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--distances", type=int, nargs="+", default=DISTANCES_KM)
    parser.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
    args = parser.parse_args()

    per_distance = {}
    for dist_km in args.distances:
        counts = []
        for seed in args.seeds:
            n_entangled, _times = run_once(dist_km, seed)
            counts.append(n_entangled)
            print(f"dist={dist_km:3d} km  seed={seed}  entangled={n_entangled}/{N_MEMORIES}", flush=True)
        per_distance[dist_km] = {
            "counts_per_seed": counts,
            "median_count": statistics.median(counts),
            "mean_count": statistics.fmean(counts),
        }

    medians = [per_distance[d]["median_count"] for d in args.distances]
    # Direction check: is the median entangled-memory count non-increasing
    # as distance grows? (Allow ties -- SeQUeNCe's rate is a stochastic,
    # count-saturating quantity, not required to be strictly monotone with
    # only 5 seeds per point.)
    n_violations = sum(
        1 for i in range(1, len(medians)) if medians[i] > medians[i - 1]
    )
    is_monotone_nonincreasing = n_violations == 0

    analytical_rates = []
    for dist_km in args.distances:
        loss_db = FIBER_ATTENUATION_DB_PER_KM * dist_km
        eta = 10 ** (-loss_db / 10)
        analytical_rates.append(eta)
    analytical_monotone = all(
        analytical_rates[i] <= analytical_rates[i - 1] for i in range(1, len(analytical_rates))
    )

    def _rank(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        ranks = [0.0] * len(values)
        for r, i in enumerate(order):
            ranks[i] = r
        return ranks

    rank_med = _rank(medians)
    rank_ana = _rank(analytical_rates)
    n = len(medians)
    if n > 1:
        d2_sum = sum((rank_med[i] - rank_ana[i]) ** 2 for i in range(n))
        spearman_rho = 1 - (6 * d2_sum) / (n * (n**2 - 1))
    else:
        spearman_rho = float("nan")

    result = {
        "description": (
            "Directional cross-check of BDCZ oracle's fibre transmittance-vs-distance "
            "law against SeQUeNCe's independent discrete-event simulation of the same "
            "midpoint-heralded (DLCZ/BDCZ) elementary-link protocol. Validates RATE "
            "direction only -- see module docstring for why fidelity is out of scope "
            "for this specific cross-check."
        ),
        "sequence_version": "0.8.5",
        "fiber_attenuation_db_per_km": FIBER_ATTENUATION_DB_PER_KM,
        "sim_time_ms": SIM_TIME_MS,
        "n_memories_per_node": N_MEMORIES,
        "seeds": args.seeds,
        "distances_km": args.distances,
        "per_distance": per_distance,
        "median_counts_by_distance": medians,
        "n_monotonicity_violations": n_violations,
        "sequence_rate_is_monotone_nonincreasing": is_monotone_nonincreasing,
        "analytical_eta_by_distance": analytical_rates,
        "analytical_rate_is_monotone_nonincreasing": analytical_monotone,
        "directional_agreement": is_monotone_nonincreasing and analytical_monotone,
        "spearman_rho_sequence_vs_analytical": spearman_rho,
    }

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=2)

    print(json.dumps(result, indent=2))
    if not result["directional_agreement"]:
        print("WARNING: directional agreement check FAILED.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
