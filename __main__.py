"""Command-line entry point:  python -m ld_maker ...

Example
-------
    python -m ld_maker --reference modules/run.gro \
        --nx 12 --ny 12 --n-trio 200 \
        --composition POPC:176,DOPE:74,SAPI:20 \
        --core-thickness 40 --water 15 \
        --out planar_ld
"""

import argparse

from .builder import build_from_reference
from .bilayer import bilayer_to_ld, replace_lipid
from .topology import write_top


def _parse_composition(text):
    comp = {}
    for item in text.split(","):
        name, _, weight = item.partition(":")
        comp[name.strip()] = float(weight) if weight else 1.0
    return comp


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="ld_maker",
        description="Build a planar lipid-droplet (PL / TRIO core / PL) starting structure.",
    )
    p.add_argument("--bilayer", default=None,
                   help="BILAYER MODE: path to a solvated phospholipid bilayer .gro; "
                        "inserts a TRIO core between its leaflets (reuses its water/ions)")
    p.add_argument("--replace", default=None,
                   help="REPLACE MODE (with --bilayer): substitute a leaflet lipid, "
                        "e.g. POPC:CHYO:0.5 to swap half the POPC for CHYO in place")
    p.add_argument("--core", default=None,
                   help="bilayer mode: mixed core composition, e.g. TRIO:120,CHYO:30 "
                        "(overrides --n-trio; needs --core-source containing every species)")
    p.add_argument("--core-source", default=None,
                   help="bilayer mode: structure to take core templates from "
                        "(default: packaged TRIO+CHYO templates; override to use your own)")
    p.add_argument("--trio-source", default=None,
                   help="bilayer mode: alias of --core-source (kept for compatibility)")
    p.add_argument("--packing-fraction", type=float, default=0.35,
                   help="bilayer mode: initial core fill (looser = safer minimization)")
    p.add_argument("--reference", default=None,
                   help="grid mode: structure with one copy of each residue (e.g. modules/run.gro)")
    p.add_argument("--nx", type=int, default=12, help="lipids per leaflet along x")
    p.add_argument("--ny", type=int, default=12, help="lipids per leaflet along y")
    p.add_argument("--n-trio", type=int, default=200, help="TRIO molecules in the core")
    p.add_argument("--composition", default="POPC:176,DOPE:74,SAPI:20",
                   help="phospholipid ratios, e.g. POPC:176,DOPE:74,SAPI:20")
    p.add_argument("--core-thickness", type=float, default=40.0, help="core slab thickness (A)")
    p.add_argument("--water", type=float, default=15.0, help="empty pad above/below for solvation (A)")
    p.add_argument("--apl", type=float, default=65.0, help="target area per phospholipid (A^2)")
    p.add_argument("--d-min", type=float, default=1.6, help="min inter-molecular contact target (A)")
    p.add_argument("--no-relax", action="store_true", help="skip clash relaxation (faster, clashy)")
    p.add_argument("--seed", type=int, default=0, help="RNG seed")
    p.add_argument("--out", default="planar_ld", help="output basename (writes .gro and .top)")
    p.add_argument("--vmd", action="store_true", help="also write a <out>.vmd view script")
    args = p.parse_args(argv)

    gro, top = f"{args.out}.gro", f"{args.out}.top"

    if args.bilayer and args.replace:
        # --- replace mode: substitute a leaflet lipid in place ---
        reps = {}
        for spec in args.replace.split(";"):
            old, new, frac = spec.split(":")
            reps[old.strip()] = (new.strip(), float(frac))
        print(f"[ld_maker] replace mode: {reps} in {args.bilayer}")
        u = replace_lipid(args.bilayer, reps, out_gro=gro, out_top=top, seed=args.seed)
    elif args.bilayer:
        # --- bilayer -> LD: insert a neutral-lipid core between the leaflets ---
        core_t = args.core_thickness if args.core_thickness != 40.0 else None
        if args.core:
            core = {n: int(w) for n, w in _parse_composition(args.core).items()}
        else:
            core = {"TRIO": args.n_trio}
        core_source = args.core_source or args.trio_source
        print(f"[ld_maker] bilayer mode: inserting core {core} into {args.bilayer}")
        u = bilayer_to_ld(
            args.bilayer, core=core, out_gro=gro, out_top=top,
            core_source=core_source, core_thickness=core_t,
            packing_fraction=args.packing_fraction, seed=args.seed,
            relax=not args.no_relax, d_min=args.d_min,
        )
    else:
        # --- from-scratch grid build ---
        if not args.reference:
            p.error("grid mode needs --reference (or use --bilayer for bilayer mode)")
        comp = _parse_composition(args.composition)
        print(f"[ld_maker] grid mode: {args.nx}x{args.ny} leaflets, {args.n_trio} TRIO, comp={comp}")
        u = build_from_reference(
            args.reference, nx=args.nx, ny=args.ny, n_trio=args.n_trio,
            pl_composition=comp, core_thickness=args.core_thickness,
            water_thickness=args.water, apl=args.apl, seed=args.seed,
            relax=not args.no_relax, d_min=args.d_min,
        )
        u.atoms.write(gro)
        write_top(u, top)

    print(f"[ld_maker] wrote {gro} ({u.atoms.n_atoms} atoms) and {top}")
    print(f"[ld_maker] box = {u.dimensions[:3].round(1)} A -- minimise in GROMACS before use")

    if args.vmd:
        from .vmd import write_vmd
        vmd_path = f"{args.out}.vmd"
        write_vmd(gro, vmd_path)
        print(f"[ld_maker] wrote {vmd_path}  (vmd -e {vmd_path})")


if __name__ == "__main__":
    main()
