# Alzheimer's (AChE) in-silico screening kit — data → model → app

A minimal, working, three-part pipeline: two Colab notebooks that build and
validate a QSAR potency model from real ChEMBL data, and a Streamlit app that
serves it. Nothing here needs a local Python install.

## Why this works on an old macOS (10.13.2 / High Sierra)

Every heavy step — data download, descriptor calculation, model training,
model hosting — runs on Google's or Streamlit's servers, not on your machine.
Your Mac only needs a browser to:
1. Open Colab and click "Run all."
2. Upload two small files (`cleaned_data.csv`, then `model.pkl` +
   `feature_names.json`) to GitHub, via github.com's web upload — no `git`
   command line needed.
3. Point Streamlit Community Cloud at your GitHub repo.

The one caveat: High Sierra (2017) can no longer run the very latest Chrome
or Safari. If colab.research.google.com or streamlit.io feel broken or won't
load, install a current version of **Firefox** — it maintains support for
older macOS versions longer than Chrome or Safari do — and use that instead.

## Files in this kit

| File | Purpose | Where it runs |
|---|---|---|
| `01_data_collection_cleaning.ipynb` | Pulls acetylcholinesterase (ChEMBL target `CHEMBL220`) bioactivity data, cleans it, computes descriptors | Google Colab |
| `02_model_training_validation.ipynb` | Trains a RandomForest QSAR model, validates with cross-validation + held-out test set | Google Colab |
| `03_docking_receptor_prep.ipynb` | One-time: prepares `receptor.pdbqt` + `box_config.json` for the Docking tab | Google Colab |
| `04_md_stability_summary.ipynb` | Lightweight conformational-rigidity proxy for your top hits (see note below) | Google Colab |
| `app.py` | Streamlit app: **Screening**, **Docking**, and **MD summary** tabs sharing one results table | Streamlit Community Cloud |
| `requirements.txt` | Exact package versions for the app | Streamlit Community Cloud |

## How the three tabs fit together

The app now has three tabs, all reading from and writing to **one shared
table** (`st.session_state["compound_table"]`) — dock a compound after
screening it and its docking score just becomes a new column next to its
potency score, no re-uploading anything.

- **Screening** — potency (QSAR model) + druglikeness + CNS-friendliness. Same as before.
- **Docking** — runs real **AutoDock Vina** docking live in the app. This is
  light enough to actually work on free hosting. It needs a prepared receptor:
  run `03_docking_receptor_prep.ipynb` once (defaults to PDB `4EY7`, human AChE
  with donepezil bound) and commit the two files it produces —
  `receptor.pdbqt` and `box_config.json` — into the same GitHub folder as
  `app.py`. Users can also upload their own receptor/box files per session
  instead.
- **MD summary** — full protein-ligand molecular dynamics (explicit solvent,
  nanosecond-scale runs) is too heavy for free hosting, so this tab is
  intentionally an *upload*, not a live compute step. Run
  `04_md_stability_summary.ipynb` in Colab on your top docking hits — it
  measures how tightly each compound's low-energy conformers cluster, a fast
  rigidity proxy — then upload the resulting `md_summary.csv` here to merge
  it in. Treat it as a first-pass filter before investing in real MD on
  whatever survives.

## Alzheimer's targets you can swap in

`01_data_collection_cleaning.ipynb` includes a `find_target()` lookup cell so you
can search ChEMBL by name and confirm the exact ID before using it — never
hardcode a guessed ID. Two are already confirmed and ready to use directly:

| Target | Role in Alzheimer's | ChEMBL target ID |
|---|---|---|
| Acetylcholinesterase (AChE) | Target of all approved AD drugs (donepezil, rivastigmine, galantamine) | `CHEMBL220` (default in this kit) |
| Butyrylcholinesterase (BChE) | Secondary cholinesterase, rises in importance as AD progresses | `CHEMBL1914` |
| Beta-secretase 1 (BACE1) | Cleaves APP to produce amyloid-beta | look up via `find_target("secretase 1")` |
| Glycogen synthase kinase-3 beta (GSK-3β) | Drives tau hyperphosphorylation | look up via `find_target("Glycogen synthase kinase-3 beta")` |

To switch, change `TARGET_CHEMBL_ID` in notebook 1's config cell, then re-run
both notebooks — nothing else needs to change.

## Step-by-step

### 1. Generate the dataset (Colab)
1. Go to [colab.research.google.com](https://colab.research.google.com), sign in with any free Google account.
2. `File → Upload notebook` → select `01_data_collection_cleaning.ipynb`.
3. `Runtime → Run all`. This takes a few minutes and ends with `cleaned_data.csv` downloading to your Mac automatically.

### 2. Train and validate the model (Colab)
1. Open a fresh Colab tab, upload `02_model_training_validation.ipynb`.
2. `Runtime → Run all`. When prompted, upload the `cleaned_data.csv` from step 1.
3. Read the printed R² / RMSE / MAE and look at the two plots — this is your validation evidence, keep it for your records.
4. The notebook ends by downloading `model.pkl` and `feature_names.json` to your Mac.

### 3. (Optional but recommended) Prepare the docking receptor
1. Open `03_docking_receptor_prep.ipynb` in Colab, `Runtime → Run all`.
2. Downloads `receptor.pdbqt` and `box_config.json` to your Mac.

### 4. Assemble the GitHub repo
1. Create a new repository on [github.com](https://github.com) (e.g. `alzheimers-screening-app`).
2. Use the web "Add file → Upload files" button (no terminal needed) to upload, all in the repo root:
   - `app.py`
   - `requirements.txt`
   - `model.pkl`
   - `feature_names.json`
   - `receptor.pdbqt` and `box_config.json` (from step 3, if you did it — otherwise users can upload their own in the Docking tab each session)
3. Commit.

### 5. Deploy on Streamlit
1. Go to [share.streamlit.io](https://share.streamlit.io), sign in with GitHub.
2. "New app" → pick your repo → set main file to `app.py` → Deploy.
3. You'll get a public URL. Anyone can now paste SMILES, dock them, and upload an MD summary — no local setup on their end either.

### A note on `vina` and old macOS
`vina` (the AutoDock Vina Python bindings) ships pre-built wheels for Linux,
which is what Streamlit Community Cloud runs on — so this installs fine
server-side regardless of your Mac's OS version. You never install it
locally; it's only listed in `requirements.txt` for Streamlit's servers to
install.

## Scaling this up later

This kit intentionally ships a lightweight, honest MVP: one target (EGFR),
one property (potency), plus simple druglikeness rules — because that's what
reliably deploys on free-tier compute today. When you're ready to add
docking, MD, and full ADMET as more pages in the same app, the architecture
from our earlier discussion (a shared compound table that every module reads
from and writes back to) is the natural next step; AutoDock Vina and OpenMM
are the free, pip-installable tools for that, though full docking/MD is
compute-heavy enough that you'll likely want to run those jobs from Colab and
feed the results back into this Streamlit app rather than running them live
on Streamlit's free tier.
