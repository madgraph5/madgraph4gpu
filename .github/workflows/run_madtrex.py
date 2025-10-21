#!/usr/bin/env python3
import os
import sys
import subprocess
from pathlib import Path
import csv

ALLOWED_PROCESSES = [ "ee_mumu", "gg_tt", "gg_tt01g", "gg_ttg", "gg_ttgg", "gg_ttggg", "gq_ttq", "heft_gg_bb", "nobm_pp_ttW", "pp_tt012j", "smeft_gg_tttt", "susy_gg_t1t1", "susy_gg_tt" ]
MADGRAPH_CLI = Path.cwd() / ".." / ".." / "MG5aMC" / "mg5amcnlo" / "bin" / "mg5_aMC"

def generate_dat_content(process_dir: str, rwgt_card_path: Path) -> str:
    run_card = f"launch {process_dir}\n"
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
        reader = csv.DictReader(f, fieldnames=["RWGT", "VALUE", "ERROR"])
        for row in reader:
            yield float(row["VALUE"]), float(row["ERROR"])

def main() -> int:
    # Name of the directory of the process to test
    process_dir = sys.argv[1]
    # Label for the process (must be in the allowed list)
    process = process_dir.replace(".mad", "")

    # Treat current working directory as HOME
    HOME = Path.cwd()

    process_path = (HOME / process_dir).resolve()
    if not process_path.exists():
        print(f"ERROR: Process {process} not found at: {process_path}", file=sys.stderr)
        return 1

    if process not in ALLOWED_PROCESSES:
        print(
            f"ERROR: PROCESS '{process}' is not in the allowed list.\n"
            f"Allowed: {sorted(ALLOWED_PROCESSES)}",
            file=sys.stderr,
        )
        return 1

    # Check that baseline rwgt.csv exists
    baseline_csv = HOME / "CODEGEN" / "PLUGIN" / "CUDACPP_SA_OUTPUT" / "test" / "MadtRex_baseline" / f"{process}_rwgt.csv"
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
        dat_content = generate_dat_content(process_dir, rwgt_card_path)
        dat_path.write_text(dat_content, encoding="utf-8")
    except Exception as e:
        print(f"ERROR: Failed to write {dat_path}: {e}", file=sys.stderr)
        return 1

    # Run mg5_aMC with PROCESS.dat as argument, wait for completion
    LOGS = HOME / "logs"
    LOGS.mkdir(exist_ok=True)
    stdout_log = LOGS / f"mg5_{process}.stdout.log"
    stderr_log = LOGS / f"mg5_{process}.stderr.log"
    print(f"Launching: {MADGRAPH_CLI} {dat_path}")
    try:
        with stdout_log.open("wb") as out, stderr_log.open("wb") as err:
            result = subprocess.run(
                [str(MADGRAPH_CLI), str(dat_path)],
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
        print(f"ERROR: Failed to launch {MADGRAPH_CLI}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"ERROR: mg5_aMC run failed: {e}", file=sys.stderr)
        return 1

    # Remove process.dat
    dat_path.unlink(missing_ok=True)

    # Move rwgt_results.csv → HOME/baseline/PROCESS_rwgt.csv
    madtrex_csv = process_path / "rw_me" / "SubProcesses" / "rwgt_results.csv"
    if not madtrex_csv.exists():
        print(
            f"ERROR: Expected results not found at:\n  {madtrex_csv}\n"
            f"Ensure the run produced rwgt_results.csv.",
            file=sys.stderr,
        )
        return 1

    all_ok = True
    for i, ((v_base, _), (v_mad, _)) in enumerate(zip(load_csv(baseline_csv), load_csv(madtrex_csv)), start=1):
        diff = abs(v_base - v_mad)
        tol = 0.05 * v_mad

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
