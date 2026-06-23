#!/usr/bin/env python3
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import madtrex_common as common

def generate_dat_content(process: str) -> str:
    dat = common.model_import_lines(process)
    dat += f"reweight {common.RUN_NAME} -from_cards --plugin=madtrex\n"
    return dat

def restore_lhe(lhe_asset: Path, process_path: Path) -> Path:
    """Decompress the committed LHE asset into Events/<run_name>/unweighted_events.lhe.gz,
    the location MadGraph expects for an existing run (no run database needed: the banner
    is read straight from the LHE file itself)."""
    events_dir = process_path / "Events" / common.RUN_NAME
    events_dir.mkdir(parents=True, exist_ok=True)
    lhe_dst = events_dir / "unweighted_events.lhe.gz"
    shutil.copyfile(lhe_asset, lhe_dst)
    return lhe_dst

def main() -> int:
    # Name of the directory of the process to test
    process_dir = Path(sys.argv[1])
    # Label for the process (must be in the allowed list)
    process = process_dir.name.replace(".mad", "")

    # Treat current working directory as HOME
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

    # Check that baseline rwgt.csv exists
    baseline_csv = common.baseline_csv_path(HOME, process)
    if not baseline_csv.exists():
        print(
            f"ERROR: Baseline rwgt.csv not found at:\n  {baseline_csv}\n"
            f"Ensure the baseline file exists before running.",
            file=sys.stderr,
        )
        return 1

    # Check that the pre-generated LHE asset exists, and restore it into the run directory
    lhe_asset = common.baseline_lhe_path(HOME, process)
    if not lhe_asset.exists():
        print(
            f"ERROR: LHE asset not found at:\n  {lhe_asset}\n"
            f"Generate it with generate_madtrex_assets.py and commit it before running.",
            file=sys.stderr,
        )
        return 1
    restore_lhe(lhe_asset, process_path)

    # Write reweight_card.dat directly into the process Cards directory (consumed via -from_cards)
    rwgt_card_path = process_path / "Cards" / "reweight_card.dat"
    common.write_rwgt_card(rwgt_card_path)

    # Write PROCESS.dat to HOME
    dat_path = HOME / f"{process}.dat"
    dat_path.write_text(generate_dat_content(process), encoding="utf-8")

    # Check that the CUDACPP_OUTPUT plugin is present: required for MadtRex reweighting
    error = common.check_cudacpp_plugin_present(HOME)
    if error:
        print(error, file=sys.stderr)
        return 1

    # Run bin/madevent with PROCESS.dat as argument, wait for completion
    madevent_bin = process_path / "bin" / "madevent"
    LOGS = HOME / "logs"
    LOGS.mkdir(exist_ok=True)
    stdout_log = LOGS / f"mg5_{process}.stdout.log"
    stderr_log = LOGS / f"mg5_{process}.stderr.log"
    print(f"Launching: {madevent_bin} {dat_path}")
    try:
        with stdout_log.open("wb") as out, stderr_log.open("wb") as err:
            result = subprocess.run(
                [str(madevent_bin), str(dat_path)],
                cwd=str(process_path),
                stdout=out,
                stderr=err,
                check=False,
            )
        if result.returncode != 0:
            print(
                f"ERROR: bin/madevent exited with code {result.returncode}. "
                f"See logs:\n  stdout: {stdout_log}\n  stderr: {stderr_log}",
                file=sys.stderr,
            )
            common.dump_logs(stdout_log, stderr_log)
            return result.returncode
        else:
            print(f"bin/madevent finished. Logs:\n  stdout: {stdout_log}\n  stderr: {stderr_log}")
    except FileNotFoundError:
        print(f"ERROR: Failed to launch {madevent_bin}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"ERROR: bin/madevent run failed: {e}", file=sys.stderr)
        return 1

    # Remove process.dat
    dat_path.unlink(missing_ok=True)

    # Get rwgt_results.csv results file
    madtrex_csv = process_path / "rw_me" / "SubProcesses" / "rwgt_results.csv"
    if not madtrex_csv.exists():
        print(
            f"ERROR: Expected results not found at:\n  {madtrex_csv}\n"
            f"Ensure the run produced rwgt_results.csv.",
            file=sys.stderr,
        )
        common.dump_logs(stdout_log, stderr_log)
        return 1

    if not common.compare_csv(baseline_csv, madtrex_csv):
        print(f"Some checks failed for process {process}.", file=sys.stderr)
        sys.exit(1)
    print(f"All checks passed for process {process}.")

    return 0

if __name__ == "__main__":
    sys.exit(main())
