#!/usr/bin/env python3
# sweep_ngrid_collect.py
# Use a given pairing strength to run ngrid=51..129 and collect results from for011.

import argparse, os, re, subprocess, sys, csv
from typing import Optional, Dict, Any

# Rf = [0, 1, 2, 3]
# strengths = [456.824428, 1529.696198, 10996.526316, 30235.012020]

def parse_fortran_float(x: str) -> float:
    return float(x.replace("D","e").replace("d","e"))

LINE_RE = re.compile(
    r"Z,N,r_C,E,epair,gap,sumuv=\s*"
    r"(\d+)\s+"                 # 1: Z
    r"(\d+)\s+"                 # 2: N
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # 3: r_C
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # 4: E
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # 5: epair_1
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # 6: epair_2
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # 7: gap_1
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # 8: gap_2
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"   # 9: sumuv_1
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)"      # 10: sumuv_2
)

def write_for005(*, ngrid: int, hfbwidth: str, ipair: int, pair_strength: float, icut: int) -> None:
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
    with open("for005", "rb") as fin, open(os.devnull, "wb") as devnull:
        subprocess.run([exe_path], stdin=fin, stdout=devnull, stderr=devnull, check=True)

def parse_for011(expect_Z: Optional[int], expect_N: Optional[int], for011_path: str="for011") -> Dict[str, Any]:
    if not os.path.exists(for011_path):
        raise RuntimeError("for011 not found after solver run.")
    chosen = None
    with open(for011_path, "r", errors="ignore") as f:
        for line in f:
            m = LINE_RE.search(line)
            if m:
                if expect_Z is not None and expect_N is not None:
                    Z, N = int(m.group(1)), int(m.group(2))
                    if Z == expect_Z and N == expect_N:
                        chosen = m
                else:
                    chosen = m  # keep last match
    if chosen is None:
        raise RuntimeError("No 'Z,N,r_C,E,epair,gap,sumuv=' line found.")
    vals = {
        "Z": int(chosen.group(1)),
        "N": int(chosen.group(2)),
        "r_C": parse_fortran_float(chosen.group(3)),
        "E": parse_fortran_float(chosen.group(4)),
        "epair1": parse_fortran_float(chosen.group(5)),
        "epair2": parse_fortran_float(chosen.group(6)),
        "gap1": parse_fortran_float(chosen.group(7)),
        "gap2": parse_fortran_float(chosen.group(8)),
        "sumuv1": parse_fortran_float(chosen.group(9)),
        "sumuv2": parse_fortran_float(chosen.group(10)),
        "raw_line": chosen.group(0).strip(),
    }
    return vals

def main():
    ap = argparse.ArgumentParser(description="Sweep ngrid=51..129 using a fixed pairing strength and collect for011 results.")
    ap.add_argument("--exe-path", default="../hafom_huge.exe", help="Path to your Fortran executable")
    ap.add_argument("--hfbwidth", required=True, help="e.g. 1.0d0")
    ap.add_argument("--pair-strength", required=True, type=parse_fortran_float, help="Optimized pairing strength to use")
    ap.add_argument("--ipair", type=int, default=8, choices=[7,8], help="7=BCS, 8=HFB (default 8)")
    ap.add_argument("--icut", type=int, default=25)
    ap.add_argument("--ngrid-start", type=int, default=17)
    ap.add_argument("--ngrid-stop", type=int, default=129)
    ap.add_argument("--ngrid-step", type=int, default=17)
    # optional bookkeeping names (match your shell style if you want)
    ap.add_argument("--edf", default="SVbas")
    ap.add_argument("--nucleus", default="unknown")
    ap.add_argument("--outdir", default="collect_cat_box_paper_new")
    ap.add_argument("--expect-Z", type=int, default=None, help="If multiple lines exist in for011, prefer this Z")
    ap.add_argument("--expect-N", type=int, default=None, help="If multiple lines exist in for011, prefer this N")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    cond = "bcs" if args.ipair == 7 else "hfb"

    dat_name = f"collect_cat_{args.edf}_{args.nucleus}_{args.hfbwidth}_{cond}_Cut{args.icut}.dat"
    csv_name = f"sweep_results_test_{cond}_Cut{args.icut}.csv"

    # CSV header
    with open(csv_name, "w", newline="") as fout:
        w = csv.writer(fout)
        w.writerow(["ngrid","hfbwidth","ipair","pair_strength",
                    "Z","N","r_C","E","epair1","epair2","gap1","gap2","sumuv1","sumuv2"])

    # sweep
    for ngrid in range(args.ngrid_start, args.ngrid_stop+1, args.ngrid_step):
        try:
            write_for005(ngrid=ngrid, hfbwidth=args.hfbwidth, ipair=args.ipair,
                         pair_strength=args.pair_strength, icut=args.icut)
            run_solver(args.exe_path)
            vals = parse_for011(args.expect_Z, args.expect_N)
            process = subprocess.run(['cp', 'for011', f"for011_Rf3/for011_{ngrid}.txt"], check=True)

            # append to CSV
            with open(csv_name, "a", newline="") as fout:
                w = csv.writer(fout)
                w.writerow([ngrid, args.hfbwidth, args.ipair, f"{args.pair_strength:.6f}",
                            vals["Z"], vals["N"], f"{vals['r_C']:.6f}", f"{vals['E']:.6f}",
                            f"{vals['epair1']:.6f}", f"{vals['epair2']:.6f}",
                            f"{vals['gap1']:.6f}", f"{vals['gap2']:.6f}",
                            f"{vals['sumuv1']:.6f}", f"{vals['sumuv2']:.6f}"])

            # append a human-readable .dat with the raw line (like grep)
            with open(dat_name, "a") as fdat:
                fdat.write(f"ngrid={ngrid},hfbwidth={args.hfbwidth},pair.strength={args.pair_strength:.6f}\n")
                fdat.write(vals["raw_line"] + "\n")

            print(f"[OK] ngrid={ngrid}: gaps=({vals['gap1']:.6f},{vals['gap2']:.6f})")
            os.replace("for011", os.path.join(args.outdir, f"for011_{ngrid}_{args.hfbwidth}.txt"))

        except subprocess.CalledProcessError:
            print(f"[ERR] ngrid={ngrid}: executable failed", file=sys.stderr)
        except Exception as e:
            print(f"[ERR] ngrid={ngrid}: {e}", file=sys.stderr)

    # move .dat to outdir (like your shell script)
    try:
        os.replace(dat_name, os.path.join(args.outdir, dat_name))
        print(f"[DONE] moved {dat_name} -> {args.outdir}/")
        print(f"[DONE] CSV: {csv_name}")
    except Exception as e:
        print(f"[WARN] could not move {dat_name}: {e}")
        print(f"[DONE] CSV: {csv_name}")

if __name__ == "__main__":
    main()
