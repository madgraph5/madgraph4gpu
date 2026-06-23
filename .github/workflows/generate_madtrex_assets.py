#!/usr/bin/env python3
"""Script to (re)generate the MadtRex CI test assets for one process:
 - an LHE event sample (test/MadtRex_baseline/<process>_events.lhe.gz)
 - the corresponding reweight baseline (test/MadtRex_baseline/<process>_rwgt.csv)

This is NOT run in CI. The idea is to run it locally every now and then
(whenever the generated *.mad code or something crucial changes) and commit the
resulting files, the same way epochX/cudacpp/CODEGEN/allGenerateAndCompare.sh
is used to regenerate reference outputs.

Usage (from epochX/cudacpp):
  ../../.github/workflows/generate_madtrex_assets.py gg_tt.mad
"""
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import madtrex_common as common

MADGRAPH_CLI = Path.cwd() / ".." / ".." / "MG5aMC" / "mg5amcnlo" / "bin" / "mg5_aMC"

def generate_dat_content(process_dir: Path, rwgt_card_path: Path, process: str) -> str:
    dat = common.model_import_lines(process)
    dat += f"launch {process_dir}\n"
    dat += "reweight=madtrex\n"
    dat += "set nevents 10000\n"
    dat += f"set iseed {common.ISEED}\n"
    dat += f"{rwgt_card_path}\n"
    return dat

def main() -> int:
    process_dir = Path(sys.argv[1])
    process = process_dir.name.replace(".mad", "")

    HOME = Path.cwd()
    process_path = (HOME / process_dir).resolve()
    if not process_path.exists():
        print(f"ERROR: Process {process} not found at: {process_path}", file=sys.stderr)
        return 1

    if process not in common.ALLOWED_PROCESSES:
        print(
            f"ERROR: PROCESS '{process}' is not in the allowed list.\n"
            f"Allowed: {sorted(common.ALLOWED_PROCESSES)}",
            file=sys.stderr,
        )
        return 1

    rwgt_card_path = HOME / "rwgt_card.dat"
    common.write_rwgt_card(rwgt_card_path)

    dat_path = HOME / f"{process}.dat"
    dat_path.write_text(generate_dat_content(process_dir, rwgt_card_path, process), encoding="utf-8")

    LOGS = HOME / "logs"
    LOGS.mkdir(exist_ok=True)
    stdout_log = LOGS / f"mg5_{process}.stdout.log"
    stderr_log = LOGS / f"mg5_{process}.stderr.log"
    print(f"Launching: {MADGRAPH_CLI} {dat_path}")
    with stdout_log.open("wb") as out, stderr_log.open("wb") as err:
        result = subprocess.run(
            [str(MADGRAPH_CLI), str(dat_path)],
            cwd=str(HOME),
            stdout=out,
            stderr=err,
            check=False,
        )
    if result.returncode != 0:
        print(f"ERROR: mg5_aMC exited with code {result.returncode}.", file=sys.stderr)
        common.dump_logs(stdout_log, stderr_log)
        return result.returncode
    print(f"mg5_aMC finished. Logs:\n  stdout: {stdout_log}\n  stderr: {stderr_log}")

    dat_path.unlink(missing_ok=True)

    # Locate the generated run's LHE file (the only run in a freshly generated process dir)
    run_dirs = sorted((process_path / "Events").glob("run_*"))
    if not run_dirs:
        print(f"ERROR: No Events/run_* directory found under {process_path}", file=sys.stderr)
        common.dump_logs(stdout_log, stderr_log)
        return 1
    if len(run_dirs) > 1:
        print(f"WARNING: Multiple runs found, using the last one: {[d.name for d in run_dirs]}")
    run_dir = run_dirs[-1]
    lhe_src = run_dir / "unweighted_events.lhe.gz"
    if not lhe_src.exists():
        print(f"ERROR: Expected LHE file not found at: {lhe_src}", file=sys.stderr)
        common.dump_logs(stdout_log, stderr_log)
        return 1

    madtrex_csv = process_path / "rw_me" / "SubProcesses" / "rwgt_results.csv"
    if not madtrex_csv.exists():
        print(f"ERROR: Expected results not found at: {madtrex_csv}", file=sys.stderr)
        common.dump_logs(stdout_log, stderr_log)
        return 1

    lhe_dst = common.baseline_lhe_path(HOME, process)
    csv_dst = common.baseline_csv_path(HOME, process)
    lhe_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(lhe_src, lhe_dst)
    shutil.copyfile(madtrex_csv, csv_dst)
    print(f"Updated asset: {lhe_dst}")
    print(f"Updated asset: {csv_dst}")
    print("Please review and commit these files.")

    return 0

if __name__ == "__main__":
    sys.exit(main())
