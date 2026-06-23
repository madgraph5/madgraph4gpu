#!/usr/bin/env python3
"""Shared constants/helpers for the MadtRex CI test (run_madtrex.py) and the
dev-only asset regeneration script (generate_madtrex_assets.py)."""
import csv
import os
import re
from pathlib import Path

ALLOWED_PROCESSES = [ "ee_mumu", "gg_tt", "gg_tt01g", "gg_ttg", "gg_ttgg", "gg_ttggg", "gq_ttq", "heft_gg_bb", "nobm_pp_ttW", "pp_tt012j", "smeft_gg_tttt", "susy_gg_t1t1", "susy_gg_tt" ]

PROCESSES_NON_TRIVIAL_MODELS = {
    "heft_gg_bb": "heft",
    "smeft_gg_tttt": "SMEFTsim_topU3l_MwScheme_UFO-massless",
    "susy_gg_t1t1": "MSSM_SLHA2",
    "susy_gg_tt": "MSSM_SLHA2",
}

# Run name used for the (real or restored) event sample: must match MadGraph's
# default first-run naming, since the asset is restored without a results database.
RUN_NAME = "run_01"

ISEED = 489

BASELINE_DIR = Path("CODEGEN") / "PLUGIN" / "CUDACPP_SA_OUTPUT" / "test" / "MadtRex_baseline"

def baseline_csv_path(home: Path, process: str) -> Path:
    return home / BASELINE_DIR / f"{process}_rwgt.csv"

def baseline_lhe_path(home: Path, process: str) -> Path:
    return home / BASELINE_DIR / f"{process}_events.lhe.gz"

def model_import_lines(process: str) -> str:
    """Lines needed to pre-convert/import a non-trivial model before reweighting,
    so that if it is still Python 2, it is converted before the procedure takes place."""
    if process not in PROCESSES_NON_TRIVIAL_MODELS:
        return ""
    model = PROCESSES_NON_TRIVIAL_MODELS[process]
    return f"set auto_convert_model True\nimport model {model}\n"

def is_executable(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)

def mg5amcnlo_dir(home: Path) -> Path:
    """Location of the MG5aMC/mg5amcnlo checkout relative to HOME (epochX/cudacpp)."""
    return home / ".." / ".." / "MG5aMC" / "mg5amcnlo"

def check_cudacpp_plugin_present(home: Path) -> str:
    """Check that PLUGIN/CUDACPP_OUTPUT exists (directory or symlink to one) inside
    the mg5amcnlo checkout: it is required for MadtRex reweighting to work.
    Returns an error message if missing, or an empty string if the check passes."""
    plugin_dir = mg5amcnlo_dir(home) / "PLUGIN" / "CUDACPP_OUTPUT"
    if not plugin_dir.is_dir():
        return (
            f"ERROR: CUDACPP_OUTPUT plugin not found at:\n  {plugin_dir}\n"
            f"It is required for MadtRex reweighting. Create it, e.g. with:\n"
            f"  cd {plugin_dir.parent} && ln -s ../../MG5aMC_PLUGIN/CUDACPP_OUTPUT ./"
        )
    return ""

def set_mg5_path(me5_configuration_path: Path, mg5_path: Path) -> None:
    """Set (uncomment/overwrite) the mg5_path entry in me5_configuration.txt so that
    bin/madevent can find the mg5amcnlo checkout needed for MadtRex reweighting,
    without going through bin/mg5_aMC's 'launch' command."""
    text = me5_configuration_path.read_text(encoding="utf-8")
    text = re.sub(r"^#?\s*mg5_path\s*=.*$", f"mg5_path = {mg5_path}", text, flags=re.MULTILINE)
    me5_configuration_path.write_text(text, encoding="utf-8")

def write_rwgt_card(path: Path) -> None:
    if path.exists():
        return
    content = """launch\nset sminputs 1 scan:[j for j in range(100,200,10)]\n"""
    path.write_text(content, encoding="utf-8")

def load_csv(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f, fieldnames=["RWGT", "VALUE", "ERROR"])
        for row in reader:
            yield float(row["VALUE"]), float(row["ERROR"])

def dump_logs(stdout_log, stderr_log):
    print("Dumping run logs...")
    print("==== STDOUT ====")
    with open(stdout_log, "r") as file:
        print(file.read())
    print("\n\n==== STDERR ====")
    with open(stderr_log, "r") as file:
        print(file.read())
    print("================")

def compare_csv(baseline_csv: Path, madtrex_csv: Path) -> bool:
    all_ok = True
    for i, ((v_base, _), (v_mad, _)) in enumerate(zip(load_csv(baseline_csv), load_csv(madtrex_csv)), start=1):
        diff = abs(v_base - v_mad)
        tol = 0.05 * v_mad
        if diff >= tol:
            print(f"Error: Row {i}: |{v_base} - {v_mad}| = {diff} >= {tol}")
            all_ok = False
    return all_ok
