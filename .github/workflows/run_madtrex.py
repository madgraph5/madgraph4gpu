#!/usr/bin/env python3
import argparse
import os
import sys
import shutil
import subprocess
from pathlib import Path
import csv
ALLOWED_PROCESS_DICTIONARY = {
    "ee_mumu": {"model": "sm", "process": "e+ e- > mu+ mu-"},
    "gg_tt": {"model": "sm", "process": "g g > t t~"},
    "gg_tt01g": {"model": "sm", "process": "g g > t t~\nadd process g g > t t~ g"},
    "gg_ttg": {"model": "sm", "process": "g g > t t~ g"},
    "gg_ttgg": {"model": "sm", "process": "g g > t t~ g g"},
    "gg_ttggg": {"model": "sm", "process": "g g > t t~ g g g"},
    "gq_ttq": {"model": "sm", "process": "g q > t t~ q"},
    "heft_gg_bb": {"model": "heft", "process": "g g > b b~"},
    "nobm_pp_ttW": {"model": "sm-no_b_mass", "process": "p p > t t~ w+"},
    "pp_tt012j": {"model": "sm", "process": "p p > t t~\nadd process p p > t t~ j\nadd process p p > t t~ j j"},
    "smeft_gg_tttt": {"model": "SMEFTsim_topU3l_MwScheme_UFO-massless", "process": "g g > t t~ t t~"},
    "susy_gg_t1t1": {"model": "MSSM_SLHA2", "process": "g g > t1 t1~"},
    "susy_gg_tt": {"model": "MSSM_SLHA2", "process": "g g > t t~"},
}


def generate_dat_content(process: str, rwgt_card_path: Path) -> str:
    if not process in ALLOWED_PROCESS_DICTIONARY:
        raise ValueError(f"Process '{process}' is not in the allowed processes.")
    proc_info = ALLOWED_PROCESS_DICTIONARY[process]
    proc = proc_info.get("process", "unknown_process")
    model = proc_info.get("model", "unknown_model")
    if proc == "unknown_process" or model == "unknown_model":
        raise ValueError(f"Process '{process}' is not properly defined in the dictionary.")
    run_card = f"import model {model}\n"
    run_card += "define q = u c d s u~ c~ d~ s~\n"
    run_card += f"generate {proc}\n"
    run_card += f"output {process}.rw\n"
    run_card += "launch\n"
    run_card += "reweight=madtrex\n"
    run_card += "0\n"
    run_card += "set nevents 10000\n"
    run_card += "set iseed 489\n"
    run_card += f"{rwgt_card_path}\n"
    run_card += "0\n"
    return run_card

def is_executable(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)

def write_rwgt_card(path: Path) -> None:
    if path.exists():
        return
    content = """launch\nset sminputs 1 scan:[j for j in range(100,200,10)]\n"""
    path.write_text(content, encoding="utf-8")


def load_csv(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            yield float(row["VALUE"]), float(row["ERROR"])

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run mg5_aMC for a validated PROCESS and move reweighting results."
    )
    parser.add_argument(
        "mg5_rel_path",
        help="Relative path (from current working directory HOME) to the mg5_aMC executable, e.g. SOMETHING/SOMETHING/bin/mg5_aMC",
    )
    parser.add_argument(
        "process",
        help="Label for the process (must be in the allowed list)."
    )
    args = parser.parse_args()

    # Treat current working directory as HOME
    HOME = Path.cwd()

    # Resolve mg5 executable path from HOME
    mg5_path = (HOME / args.mg5_rel_path).resolve()
    if not mg5_path.exists():
        print(f"ERROR: mg5_aMC not found at: {mg5_path}", file=sys.stderr)
        return 1
    if not is_executable(mg5_path):
        print(f"ERROR: Not executable: {mg5_path}", file=sys.stderr)
        return 1

    process = args.process.strip()
    if process not in ALLOWED_PROCESS_DICTIONARY:
        print(
            f"ERROR: PROCESS '{process}' is not in the allowed list.\n"
            f"Allowed: {sorted(ALLOWED_PROCESS_DICTIONARY.keys())}",
            file=sys.stderr,
        )
        return 1

    # Check that baseline rwgt.csv exists
    baseline_csv = HOME / "baseline" / f"{process}_rwgt.csv"
    if not baseline_csv.exists():
        print(
            f"ERROR: Baseline rwgt.csv not found at:\n  {baseline_csv}\n"
            f"Ensure the baseline file exists before running.",
            file=sys.stderr,
        )
        return 1

    # Write rwgt_card.dat if not exists
    rwgt_card_path = HOME / "rwgt_card.dat"
    try:
        write_rwgt_card(rwgt_card_path)
    except Exception as e:
        print(f"ERROR: Failed to write rwgt_card.dat: {e}", file=sys.stderr)
        return 1

    # Write PROCESS.dat to HOME
    dat_path = HOME / f"{process}.dat"
    try:
        dat_content = generate_dat_content(process, rwgt_card_path)
        dat_path.write_text(dat_content, encoding="utf-8")
    except Exception as e:
        print(f"ERROR: Failed to write {dat_path}: {e}", file=sys.stderr)
        return 1

    # Run mg5_aMC with PROCESS.dat as argument, wait for completion
    LOGS = HOME / "logs"
    LOGS.mkdir(exist_ok=True)
    stdout_log = LOGS / f"mg5_{process}.stdout.log"
    stderr_log = LOGS / f"mg5_{process}.stderr.log"
    print(f"Launching: {mg5_path} {dat_path}")
    try:
        with stdout_log.open("wb") as out, stderr_log.open("wb") as err:
            result = subprocess.run(
                [str(mg5_path), str(dat_path)],
                cwd=str(HOME),
                stdout=out,
                stderr=err,
                check=False,
            )
        if result.returncode != 0:
            print(
                f"ERROR: mg5_aMC exited with code {result.returncode}. "
                f"See logs:\n  stdout: {stdout_log}\n  stderr: {stderr_log}",
                file=sys.stderr,
            )
            return result.returncode
        else:
            print(f"mg5_aMC finished. Logs:\n  stdout: {stdout_log}\n  stderr: {stderr_log}")
    except FileNotFoundError:
        print(f"ERROR: Failed to launch mg5_aMC at {mg5_path}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"ERROR: mg5_aMC run failed: {e}", file=sys.stderr)
        return 1

    # Remove process.dat
    dat_path.unlink(missing_ok=True)

    # Move rwgt_results.csv → HOME/baseline/PROCESS_rwgt.csv
    madtrex_csv = HOME / f"{process}.rw" / "rw_me" / "SubProcesses" / "rwgt_results.csv"
    if not madtrex_csv.exists():
        print(
            f"ERROR: Expected results not found at:\n  {madtrex_csv}\n"
            f"Ensure the run produced rwgt_results.csv.",
            file=sys.stderr,
        )
        return 1

    for i, ((v_base, e_base), (v_mad, e_mad)) in enumerate(zip(load_csv(baseline_csv), load_csv(madtrex_csv)), start=1):
        diff = abs(v_base - v_mad)
        tol = min(e_base, e_mad)

        if diff >= tol:
            print(f"Error: Row {i}: |{v_base} - {v_mad}| = {diff} >= {tol}")
            all_ok = False

    if not all_ok:
        print(f"Some checks failed for process {process}.", file=sys.stderr)
        sys.exit(1)
    print(f"All checks passed for process {process}.")

    return 0

if __name__ == "__main__":
    sys.exit(main())
