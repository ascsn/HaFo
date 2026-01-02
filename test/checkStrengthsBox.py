#!/usr/bin/env python3
# fit_pair_strength.py
# Minimal driver: fits pairing strength to a target gap by running the Fortran code
# and reading the two gap values from "for011".
#
# - No bracketing, no reuse of previous fits.
# - Secant-style updates; optional step clamp and min/max clamps.
# - Default: fit to the mean of the two gaps on the "Z,N,r_C,E,epair,gap,sumuv=" line.

import argparse, os, re, subprocess, sys, time
from typing import Optional, Tuple

# ---- tiny helpers ----

def parse_fortran_float(x: str) -> float:
    return float(x.replace("D","e").replace("d","e"))

LINE_RE = re.compile(
    r"Z,N,r_C,E,epair,gap,sumuv=\s*"
    r"(\d+)\s+"                 # Z
    r"(\d+)\s+"                 # N
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # r_C
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # E
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # epair_1
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # epair_2
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # gap_1
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # gap_2
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # sumuv_1
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)"      # sumuv_2
)

def pick_gap_scalar(g1: float, g2: float, mode: str) -> float:
    if mode == "proton":  return g1
    if mode == "neutron": return g2
    if mode == "both_mean": return 0.5*(g1+g2)
    if mode == "both_max":  return max(g1, g2)
    if mode == "both_min":  return min(g1, g2)
    if mode == "both_rss":  return ((g1*g1 + g2*g2)/2.0)**0.5
    raise ValueError(f"Unknown gap mode: {mode}")

# ---- I/O with your solver ----

def write_for005(
    *,
    ngrid: int,
    hfbwidth: str,
    ipair: int,
    pair_strength: float,
    icut: int,
):
    # Mirrors your echo'd lines (kept exactly so your solver behaves the same)
    with open("for005", "w") as f:
        f.write("1,,    IWUNIT\n")
        f.write("'SV-bas HFB'\n")
        f.write(f"{ngrid},0.3,              NGRID,RSTEP\n")
        f.write("-1,1E-8,1E-9,200,20000,3,-1, IWR,DABR/E,ITPR/MAX,IFPL\n")
        f.write("0.3,100.,0.2,-1,,,,,               X0DAMP,E0DAMP,ADDNEW,WALL\n")
        f.write("1,0,0,-1,1,3,0,0,0,10,,, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk\n")
        f.write("1,,,                       ifform\n")
        f.write("0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12\n")
        f.write("0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12\n")
        f.write(" -1879.639,   313.749,   112.677, 12527.376,   124.634,      t0,t1,t2,t3,t4\n")
        f.write("0.25854594    ,-.38168900    ,-2.8236379    ,0.12322898    ,   34.11167,       x0,x1,x2,x3,b4p\n")
        f.write("   0.00000   0.20113   0.00000,,,      riftm,rho0pr,powerp\n")
        f.write("     0.000,  0.000000,  0.000000,           tyuk,xyuk,ayuk\n")
        # POWER, ACOUL, IPAIR, PAIRFP/N, ILNRESP, ..., CUT, HFBWIDTH
        f.write(f"  0.300000,1.0,{ipair},10,{pair_strength},0,0.0,-{icut},{hfbwidth},,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,...\n")
        f.write("50,70,,,7,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh\n")
        f.write(" 0         IFGSC\n")
        f.write("-1,,    IWUNIT\n")

def run_solver(exe_path: str) -> None:
    # Equivalent of: ../hafom_huge.exe < for005 > /dev/null
    with open("for005", "rb") as fin, open(os.devnull, "wb") as devnull:
        subprocess.run([exe_path], stdin=fin, stdout=devnull, stderr=devnull, check=True)

def read_gap_from_for011(
    *,
    expect_Z: Optional[int],
    expect_N: Optional[int],
    gap_mode: str,
    swap_gap_order: bool,
    for011_path: str = "for011",
) -> Tuple[float, float]:
    if not os.path.exists(for011_path):
        raise RuntimeError("for011 not found after solver run.")
    chosen = None
    with open(for011_path, "r", errors="ignore") as f:
        for line in f:
            m = LINE_RE.search(line)
            if m:
                Z = int(m.group(1)); N = int(m.group(2))
                # choose the matching isotope if specified; otherwise keep last occurrence
                if expect_Z is not None and expect_N is not None:
                    if Z == expect_Z and N == expect_N:
                        chosen = m
                else:
                    chosen = m
    if chosen is None:
        raise RuntimeError("Could not find 'Z,N,r_C,E,epair,gap,sumuv=' line in for011.")
    g1 = parse_fortran_float(chosen.group(7))
    g2 = parse_fortran_float(chosen.group(8))
    if swap_gap_order:
        g1, g2 = g2, g1
    return pick_gap_scalar(g1, g2, gap_mode), (g1 if gap_mode!="neutron" else g2)

# ---- simple secant optimizer (no bracketing) ----

def fit_strength_secant(
    *,
    exe_path: str,
    ngrid: int,
    hfbwidth: str,
    ipair: int,
    icut: int,
    start_strength: float,
    second_strength: Optional[float],
    target_gap: float,
    gap_mode: str,
    expect_Z: Optional[int],
    expect_N: Optional[int],
    swap_gap_order: bool,
    tolerance: float = 1e-3,
    max_iter: int = 25,
    min_strength: Optional[float] = None,
    max_strength: Optional[float] = None,
    step_clip: Optional[float] = None,
) -> Tuple[float, float]:
    def eval_gap(s: float) -> float:
        if min_strength is not None: s = max(min_strength, s)
        if max_strength is not None: s = min(max_strength, s)
        write_for005(ngrid=ngrid, hfbwidth=hfbwidth, ipair=ipair, pair_strength=s, icut=icut)
        run_solver(exe_path)
        g_scalar, _ = read_gap_from_for011(expect_Z=expect_Z, expect_N=expect_N,
                                           gap_mode=gap_mode, swap_gap_order=swap_gap_order)
        return g_scalar

    s0 = start_strength
    s1 = s0 + (second_strength if second_strength is not None else max(5.0, 0.05*abs(s0)))
    f0 = eval_gap(s0) - target_gap
    if abs(f0) <= tolerance: return s0, f0 + target_gap
    f1 = eval_gap(s1) - target_gap
    if abs(f1) <= tolerance: return s1, f1 + target_gap

    for _ in range(max_iter):
        denom = (f1 - f0)
        step = (f1 * (s1 - s0) / denom) if abs(denom) > 1e-12 else (0.1*(s1 - s0) if s1 != s0 else max(5.0, 0.05*abs(s1)))
        s2 = s1 - step
        if step_clip is not None:
            s2 = s1 + max(-step_clip, min(step_clip, -step))
        if min_strength is not None: s2 = max(min_strength, s2)
        if max_strength is not None: s2 = min(max_strength, s2)

        f2 = eval_gap(s2) - target_gap
        if abs(f2) <= tolerance or abs(s2 - s1) <= max(1e-6, 1e-6*abs(s2)):
            return s2, f2 + target_gap
        s0, f0, s1, f1 = s1, f1, s2, f2

    return s1, f1 + target_gap  # best we have

# ---- CLI ----

def main():
    ap = argparse.ArgumentParser(description="Fit pairing strength to target gap by running the Fortran code and reading for011.")
    ap.add_argument("--exe-path", default="../hafom_huge.exe", help="Path to your Fortran executable")
    ap.add_argument("--hfbwidth", required=True, help="e.g. 1.0d0")
    ap.add_argument("--start-strength", type=parse_fortran_float, required=True, help="Initial guess for pairing strength (e.g. 385.3)")
    ap.add_argument("--target-gap", type=float, default=1.61)
    ap.add_argument("--ngrid", type=int, default=71, help="Grid at which to fit (default 71)")
    ap.add_argument("--ipair", type=int, default=8, choices=[7,8], help="7=BCS, 8=HFB (default 8)")
    ap.add_argument("--icut", type=int, default=25)

    ap.add_argument("--gap-mode", default="both_mean", choices=["proton","neutron","both_mean","both_rss","both_max","both_min"])
    ap.add_argument("--expect-Z", type=int, default=None)
    ap.add_argument("--expect-N", type=int, default=None)
    ap.add_argument("--swap-gap-order", action="store_true", help="If your code prints neutron first")

    ap.add_argument("--second-strength", type=parse_fortran_float, default=None, help="Optional second seed for secant (default: start+5% or +5)")
    ap.add_argument("--tolerance", type=float, default=1e-3)
    ap.add_argument("--max-iter", type=int, default=25)
    ap.add_argument("--min-strength", type=parse_fortran_float, default=None)
    ap.add_argument("--max-strength", type=parse_fortran_float, default=None)
    ap.add_argument("--step-clip", type=parse_fortran_float, default=None, help="Limit per-iteration change in strength")

    args = ap.parse_args()

    try:
        s_opt, g_ach = fit_strength_secant(
            exe_path=args.exe_path,
            ngrid=args.ngrid,
            hfbwidth=args.hfbwidth,
            ipair=args.ipair,
            icut=args.icut,
            start_strength=args.start_strength,
            second_strength=args.second_strength,
            target_gap=args.target_gap,
            gap_mode=args.gap_mode,
            expect_Z=args.expect_Z,
            expect_N=args.expect_N,
            swap_gap_order=args.swap_gap_order,
            tolerance=args.tolerance,
            max_iter=args.max_iter,
            min_strength=args.min_strength,
            max_strength=args.max_strength,
            step_clip=args.step_clip,
        )
        print(f"Optimized pairing strength: {s_opt:.6f}")
        print(f"Achieved gap ({args.gap_mode}): {g_ach:.6f}  (target={args.target_gap})")
        print(f"Settings: hfbwidth={args.hfbwidth}, ngrid={args.ngrid}, ipair={args.ipair}, icut={args.icut}")
    except subprocess.CalledProcessError:
        sys.exit("ERROR: Fortran executable failed to run (check --exe-path and working directory).")
    except Exception as e:
        sys.exit(f"ERROR: {e}")

if __name__ == "__main__":
    main()
