# OC43 RdRp binder design —  runbook

Everything below were the issues I faced..🫠

!! Read the **Gotchas** section before changing anything

---

## 0. Environment / GPUs

GPU allocation is pinned inside `run_env.sh`, **not** on the command line:

```bash
grep -n "CUDA_VISIBLE_DEVICES" run_env.sh
# line 13: export CUDA_VISIBLE_DEVICES=2,3,4,5
```

To change which cards are used, edit that line. Then **match the job count** to
the number of GPUs in every run command:

```
++gen_njobs=4 ++eval_njobs=4      # 4 GPUs listed 
```

If `njobs` > number of GPUs the jobs fight over cards and slow down. If
`njobs` < GPUs you waste cards.

---

## 1. Prepare the target PDB

Put the structure in `data/`. If cropping from the full protein, **preserve the
original residue numbering** — do not renumber to 1:

```python
START, END = 461, 928

with open('data/oc43.pdb') as fin, open('data/oc43_palmcore_ext.pdb','w') as fout:
    for line in fin:
        if line[:4] in ('ATOM','HETA'):
            if START <= int(line[22:26]) <= END:
                fout.write(line)          # line copied verbatim -> numbering kept
        elif line[:3] in ('TER','END'):
            fout.write(line)
```

Verify:

```bash
grep "^ATOM" data/oc43_palmcore_ext.pdb | grep " CA " | head -1   # should show 461
grep "^ATOM" data/oc43_palmcore_ext.pdb | grep " CA " | tail -1   # should show 928
# and confirm the hotspots are actually there, with the right residue types
for r in 570 573 574 578 579; do
  echo -n "Residue $r: "
  grep "^ATOM" data/oc43_palmcore_ext.pdb | grep " CA " | awk -v r=$r '$6==r {print $4}'
done
# expect: LYS LYS SER THR ARG
```

---

## 2. Add the target config

```bash
# ALWAYS check first — running cat >> twice creates a duplicate YAML key
# and the whole file stops parsing ("found duplicate key")
grep -n "^  oc43_palmcore" configs/targets/targets_dict.yaml
```

```bash
cat >> configs/targets/targets_dict.yaml << 'YAMLEOF'

  oc43_palmcore_ext_v2:
    target_input: "A461-928"
    target_path: "/home/seda/proteina-complexa/Proteina-Complexa/data/oc43_palmcore_ext.pdb"
    source: custom_targets
    target_filename: oc43_palmcore_ext
    hotspot_residues: ["A485","A487","A488","A490","A492","A493","A494","A495","A496","A570","A573","A574","A576","A577","A578","A579","A580","A603","A604","A678","A679","A681","A684"]
    binder_length: [80, 170]
    pdb_id: oc43
YAMLEOF

grep -A 8 "oc43_palmcore_ext_v2:" configs/targets/targets_dict.yaml
```

`target_input` must match the actual residue range in the PDB.
`target_filename` is the PDB file stem — it is **not** the task name; several
task entries can share one PDB and differ only in `hotspot_residues`.

If you did create a duplicate:

```bash
grep -n "oc43_palmcore_ext_v2" configs/targets/targets_dict.yaml   # find both line numbers
sed -i '506,517d' configs/targets/targets_dict.yaml                # delete the later block
```

---

## 3. Validate

```bash
./run_env.sh complexa validate design configs/search_binder_local_pipeline.yaml
```

---

## 4. Pilot (only when something changed)

Needed after changing binder length, hotspot residues, or the crop.
Not needed for another batch of the same configuration.

```bash
./run_env.sh complexa design configs/search_binder_local_pipeline.yaml \
  ++run_name=oc43_palmcore_ext_v2_pilot \
  ++generation.task_name=oc43_palmcore_ext_v2 \
  ++generation.dataloader.dataset.nres.nsamples=2 \
  ++generation.dataloader.batch_size=2 \
  ++generation.search.algorithm=beam-search \
  ++generation.search.beam_search.beam_width=4 \
  ++generation.search.beam_search.n_branch=4 \
  ++gen_njobs=4 \
  ++eval_njobs=4 \
  ++metric.binder_folding_method=colabdesign \
  ++metric.num_redesign_seqs=1
```

`run_name` and `task_name` must be consistent, because the output directory is
`{task_name}_{run_name}` and every check below reads that path.

---

## 5. Full batch

```bash
./run_env.sh complexa design configs/search_binder_local_pipeline.yaml \
  ++run_name=oc43_palmcore_ext_v2_full4 \
  ++generation.task_name=oc43_palmcore_ext_v2 \
  ++generation.dataloader.dataset.nres.nsamples=15 \
  ++generation.dataloader.batch_size=2 \
  ++generation.search.algorithm=beam-search \
  ++generation.search.beam_search.beam_width=4 \
  ++generation.search.beam_search.n_branch=4 \
  ++gen_njobs=4 \
  ++eval_njobs=4 \
  ++seed=11 \
  ++metric.binder_folding_method=colabdesign \
  ++metric.num_redesign_seqs=4
```

A missing backslash silently drops all following overrides.

For each new batch: change `run_name` **and** `seed`. Keep `task_name` fixed.
Seeds used so far: 5 (full), 7 (full2), 9 (full3).


### Note on the nsamples

Raising `n_branch` or `beam_width` is not the way to get more designs:
`beam_width=8` runs out of GPU memory, and raising `n_branch` multiplies the
reward compute without increasing how many trajectories survive each step.
Use `nsamples`, or repeat batches with new seeds.

---

## 6. Checks after the run

Set the run label once:

```bash
RUN=oc43_palmcore_ext_v2_oc43_palmcore_ext_v2_full4
DIR=evaluation_results/search_binder_local_pipeline_$RUN
```

### 6a. AF2 metrics

```bash
python3 -c "
import pandas as pd, glob
files = glob.glob('$DIR/binder_results_*.csv')
df = pd.concat([pd.read_csv(f) for f in files])
print('N:', len(df))
print('pLDDT mean:', round(df['self_complex_pLDDT'].mean(),3))
print('min_ipAE mean:', round(df['self_complex_min_ipAE'].mean(),3))
"
```

### 6b. Full-target clash filter

Not part of Complexa. Generation only ever sees the crop, so designs can
collide with residues the model never saw. The contig argument **must match**
`target_input`:

```bash
./run_env.sh python oc43_v5_clashfilter.py "$DIR" "A461-928" 2>&1 | tail -3
```

### 6c. Core hotspot contact — the metric that actually matters

**The refolded PDBs renumber chain A to 1..N.** Hotspot indices must be
remapped: `crop_index = original_residue - crop_start + 1`.

For the A461-928 crop: A570->110, A573->113, A574->114, A578->118, A579->119.
For the old A461-885 crop the mapping is identical (same start).
For a multi-segment crop it is **not** a simple offset — count segment lengths.

```bash
cat > check_contact.py << 'PYEOF'
import glob, sys, collections
import numpy as np

BASE = sys.argv[1]
HOTSPOTS = [110, 113, 114, 118, 119]   # A570/573/574/578/579, A461-928 crop
CONTACT_DIST = 8.0

def ca(pdb):
    t, b = {}, []
    for l in open(pdb):
        if l[:4] == 'ATOM' and l[12:16].strip() == 'CA':
            x = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            if l[21] == 'A': t[int(l[22:26])] = x
            elif l[21] == 'B': b.append(x)
    return t, np.array(b)

res = []
for p in glob.glob(f'{BASE}/job_*/AF2/*_self_seq_0_model1.pdb'):
    t, b = ca(p)
    if len(b) == 0: continue
    hs = np.array([t[h] for h in HOTSPOTS if h in t])
    if len(hs) == 0: continue
    d = np.sqrt(((hs[:,None,:] - b[None,:,:])**2).sum(-1)).min(axis=1)
    res.append(int((d < CONTACT_DIST).sum()))

c = collections.Counter(res)
print(f'N: {len(res)}   mean: {np.mean(res):.2f}/5')
for k in range(6): print(f'  {k}/5: {c.get(k,0)}')
print(f'>=4/5: {sum(r>=4 for r in res)}   >=3/5: {sum(r>=3 for r in res)}')
PYEOF

python3 check_contact.py "$DIR"
```


---

## 7. Pooling batches and picking candidates

Merges several runs, computes contact, ranks by pLDDT then min_ipAE, and writes 
engaged_candidates.txt (>=4/5) and shortlist_candidates.txt (5/5 with good metrics).

**Add each new run to the RUNS dict at the top.**

```bash
cat > pool_runs.py << 'PYEOF'
import glob, os, shutil, collections
import numpy as np, pandas as pd

ROOT = "evaluation_results/search_binder_local_pipeline_oc43_palmcore_ext_v2_oc43_palmcore_ext_v2_"
RUNS = ["full", "full2", "full3"]          # add each new batch here
HOT = [110, 113, 114, 118, 119]            # A570/573/574/578/579 in the A461-928 crop
MIN_CONTACT, MIN_PLDDT, MAX_IPAE = 4, 0.83, 0.15
OUT = "shortlist"

def contact(pdb):
    t, b = {}, []
    for l in open(pdb):
        if l[:4] == "ATOM" and l[12:16].strip() == "CA":
            x = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            (t.__setitem__(int(l[22:26]), x) if l[21] == "A" else
             b.append(x) if l[21] == "B" else None)
    if not b: return None
    hs = np.array([t[h] for h in HOT if h in t])
    d = np.sqrt(((hs[:, None] - np.array(b)[None]) ** 2).sum(-1)).min(1)
    return int((d < 8).sum()), round(float(d.min()), 2), len(b)

rows = []
for run in RUNS:
    base = ROOT + run
    bf = glob.glob(f"{base}/binder_results_*.csv")
    if not bf:
        print(f"  ! {run} missing"); continue
    df = pd.concat(map(pd.read_csv, bf))
    df["id"] = df.pdb_path.str.split("/").str[-2]
    mf = glob.glob(f"{base}/monomer_results_*.csv")
    seq = {}
    if mf:
        m = pd.concat(map(pd.read_csv, mf))
        seq = dict(zip(m.pdb_path.str.split("/").str[-2], m._res_mpnn_best_sequence))
    for r in df.itertuples():
        p = glob.glob(f"{base}/{r.id}/AF2/*_self_seq_0_model1.pdb")
        c = contact(p[0]) if p else None
        if not c: continue
        rows.append(dict(design_id=r.id, run=run, contact=c[0], min_hs_dist=c[1],
                         binder_len=c[2], pLDDT=r.self_complex_pLDDT,
                         min_ipAE=r.self_complex_min_ipAE,
                         sequence=seq.get(r.id, ""), pdb=p[0]))
    print(f"  {run}: {sum(x['run'] == run for x in rows)} designs")

res = pd.DataFrame(rows).sort_values(["pLDDT", "min_ipAE"], ascending=[False, True])
res.to_csv("pooled_designs.csv", index=False)
print(f"\npooled {len(res)} -> pooled_designs.csv")
print("contact:", dict(sorted(collections.Counter(res.contact).items())))

short = res[(res.contact >= MIN_CONTACT) & (res.pLDDT >= MIN_PLDDT) &
            (res.min_ipAE <= MAX_IPAE)].copy()
print(f"\nshortlist (>={MIN_CONTACT}/5, pLDDT>={MIN_PLDDT}, ipAE<={MAX_IPAE}): "
      f"{len(short)}   unique seqs: {short.sequence.nunique()}   "
      f"lengths: {sorted(set(short.binder_len))}")
if not len(short):
    raise SystemExit("  empty -- loosen the thresholds at the top")

shutil.rmtree(OUT, ignore_errors=True); os.makedirs(OUT)
short["pdb_file"] = [f"cand{i:02d}_{r.run}_{os.path.basename(r.pdb)}"
                     for i, r in enumerate(short.itertuples(), 1)]
for r in short.itertuples():
    shutil.copy(r.pdb, f"{OUT}/{r.pdb_file}")
short.drop(columns="pdb").to_csv(f"{OUT}/candidates.csv", index=False)

print(short[["run", "contact", "pLDDT", "min_ipAE", "min_hs_dist",
             "binder_len"]].to_string(index=False))
print(f"\n-> {OUT}/ : {len(short)} PDBs + candidates.csv")
PYEOF

python3 pool_runs.py
```
This writes: 
1. `pooled_designs.csv` - every design from every batch, with contact count, nearest-hotspot distance,
   binder length, pLDDT, min_ipAE and sequence.
2. `shortlist/` - the PDBs that pass all thresholds, renamed `cand01_...` in rank order with the batch
   name embedded, plus `candidates.csv` holding design_id, run, binder length, contact count,
   nearest-hotspot distance, pLDDT, min_ipAE, the ProteinMPNN sequence, and the filename in that folder.

`candidates.csv` contains **every design** that passes, not just the ones printed to screen.

Thresholds are the three constants on one line near the top. `MIN_CONTACT = 4` is the default. 

If the shortlist comes back empty the script exits with a message loosen `MIN_PLDDT` or `MAX_IPAE`.

Check the `unique seqs` and `lengths` numbers it prints. It might come from a single beam lineage,
so sequence diversity matters more than the raw count.

`shortlist/` is deleted and rebuilt on every run, so rerunning after adding a batch is safe.

Download it:

```
tar czf shortlist.tar.gz shortlist/
```

---

## 8. Where things live

| Path | Contents |
|---|---|
| `configs/targets/targets_dict.yaml` | all target definitions |
| `configs/pipeline/binder/binder_generate.yaml` | search + reward defaults |
| `data/*.pdb` | target structures and crops |
| `inference/<task>_<run>/` | raw generated PDBs, rewards CSV, timing |
| `evaluation_results/<task>_<run>/` | refolded structures + all metrics CSVs |
| `logs/design_pipeline_<task>_<run>_<timestamp>/` | per-stage logs |

Refolded structures used for every check:
`evaluation_results/.../job_*/AF2/*_self_seq_0_model1.pdb`

---

## Gotchas

**Missing trailing backslash** — silently drops every override after it.
Always check the printed config banner at the start of the run.

**Duplicate YAML key** — `cat >>` run twice breaks the whole file. `grep -n`
before appending.

**nsteps and step_checkpoints are coupled** — `step_checkpoints[-1]` must equal
`nsteps`. Lowering nsteps alone fails with:
`step_checkpoints[-1] (400) must equal nsteps (200)`. To halve the steps:

```
++generation.args.nsteps=200 "++generation.search.step_checkpoints=[0,50,100,150,200]"
```

**CUDA out of memory** — `beam_width=8` does not fit. Stay at 4, or lower
`++generation.dataloader.batch_size`.

**Reward files only appear when generate finishes** — you cannot run `filter`
on a partially finished generate; it fails with "No reward files found!".

**Crop boundaries create blind spots** — every crop hides part of the protein
from the generator, and binders will bury themselves into the missing volume.
The A461-885 crop gave the best AF2 metrics of any run (pLDDT 0.936,
min_ipAE 0.053, contact 4.60/5) but 0/10 designs survived the clash filter,
because the binder occupied the space the C-domain actually fills. Always run
the clash filter against the full structure before trusting any metric.

**Sequence adjacency is not spatial adjacency** — the old six-segment crop
handed the model junctions where residues 11-32 A apart in the real protein
were presented as covalently bonded (A497->A560: 10.7 A, A603->A629: 29.6 A,
A760->A916: 31.9 A, versus 3.8 A for a real Ca-Ca bond). The surface it
optimised against did not exist.

---

## Reward composition (for reference)

```
af2folding        i_pae: -1.0          native Complexa reward, all other terms 0
hotspot_contact   fraction within 8 A  custom addition, scored on the AF2-refolded pose
hotspot_proximity 1/(1 + max(d-8,0)/8) custom addition; the hard 8 A cutoff alone was
                                       flat (0 for 132/136 designs) so beam search had
                                       no gradient to climb
```

`af2folding/i_pae` is Complexa's own; the other two were added for this
project. Worth stating explicitly in any writeup.
