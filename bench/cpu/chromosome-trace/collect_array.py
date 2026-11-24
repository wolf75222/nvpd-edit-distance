"""Collect independent Slurm blocks without discarding their build or host identity."""
import csv
import shutil
import sys
from pathlib import Path

from summarize import summarize


def collect(source, destination):
    blocks = sorted(source.glob("block-*"))
    if len(blocks) != 5:
        raise ValueError("five Slurm tasks required")
    destination.mkdir(exist_ok=True,parents=True)
    trials, origins = [], []
    for number, path in enumerate(blocks,1):
        with (path/"trials.csv").open() as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames
            obs = list(reader)
        if {r["method"] for r in obs} != {"split","biwfa","wfa2","edlib"} or len(obs) != 4:
            raise ValueError(f"incomplete comparator block {path}")
        if any(int(r["block"]) != number or r["case"] != "chr1" or r["mode"] != "trace" for r in obs):
            raise ValueError(f"unexpected input in {path}")
        trials.extend(obs)
        env = (path/"environment.txt").read_text().splitlines()
        origins.append((number,env[0],(path/"binary.sha256").read_text().split()[0]))
        for filename in ("sources.sha256","inputs.sha256"):
            if (path/filename).read_bytes() != (blocks[0]/filename).read_bytes():
                raise ValueError(f"{filename} differs between tasks")
    with (destination/"trials.csv").open("w") as f:
        writer = csv.DictWriter(f,fields)
        writer.writeheader()
        writer.writerows(trials)
    with (destination/"blocks.csv").open("w") as f:
        wr = csv.writer(f)
        wr.writerow(("block","host","binary_sha256"))
        wr.writerows(origins)
    for filename in ("sources.sha256","inputs.sha256"):
        shutil.copyfile(blocks[0]/filename,destination/filename)
    (destination/"environment.txt").write_text("\n".join((blocks[0]/"environment.txt").read_text().splitlines()[1:])+"\n")
    summarize(destination)


if __name__ == "__main__":
    collect(*(Path(s) for s in sys.argv[1:]))
