"""Run every check script and report.  Usage:

    python run_all.py            # scripts 01-05 and the quick version of 06 (about 1 minute)
    python run_all.py --full     # also the full timing table (adds about 1 minute)

Each script's output is also saved to results/<script>.txt.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["01_examples.py", "02_exhaustive_checks.py", "03_weighted_loops.py",
           "04_dominance_and_open_questions.py", "05_symbolic_identities.py",
           "06_timing_table1.py"]

full = "--full" in sys.argv
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
failed = []
for s in SCRIPTS:
    args = [sys.executable, s]
    if s.startswith("06") and not full:
        args.append("--quick")
    print(f"=== {s} {' '.join(args[2:])}")
    p = subprocess.run(args, cwd=os.path.join(HERE, "scripts"), capture_output=True, text=True)
    out = p.stdout + p.stderr
    name = s.replace(".py", "_quick.txt" if (s.startswith("06") and not full) else ".txt")
    with open(os.path.join(HERE, "results", name), "w") as f:
        f.write(out)
    print(out.strip().splitlines()[-1] if out.strip() else "(no output)")
    if p.returncode != 0:
        failed.append(s)

print("\nFAILED: " + ", ".join(failed) if failed else "\nAll scripts passed.")
sys.exit(1 if failed else 0)
