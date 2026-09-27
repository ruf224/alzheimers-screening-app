"""
Alzheimer's (AChE) in-silico screening app — single-platform front end.

Three tabs, one shared results table (st.session_state["compound_table"]):
  - Screening: potency prediction (model.pkl) + druglikeness / CNS flags
  - Docking:   live AutoDock Vina docking against a prepared receptor
  - MD summary: upload the lightweight stability proxy from
                04_md_stability_summary.ipynb (real MD is too heavy to run
                live on free hosting -- see that notebook for why)

Run locally:      streamlit run app.py
Deploy:            push this folder to GitHub, then deploy on
                    https://share.streamlit.io pointing at app.py

For the Docking tab to work out of the box, commit receptor.pdbqt and
box_config.json (from 03_docking_receptor_prep.ipynb) into this same folder --
otherwise users can upload their own each session.
"""

import json
import os

import joblib
import numpy as np
import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski

st.set_page_config(page_title="Alzheimer's in-silico screen", layout="wide")

FEATURES = ["MolWt", "LogP", "TPSA", "HBD", "HBA", "RotatableBonds", "RingCount", "AromaticRings"]

if "compound_table" not in st.session_state:
    st.session_state["compound_table"] = pd.DataFrame()


# ---------- shared helpers ----------

@st.cache_resource
def load_model():
    model = joblib.load("model.pkl")
    with open("feature_names.json") as f:
        feature_order = json.load(f)
    return model, feature_order


def compute_descriptors(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "RotatableBonds": Descriptors.NumRotatableBonds(mol),
        "RingCount": Descriptors.RingCount(mol),
        "AromaticRings": Lipinski.NumAromaticRings(mol),
    }


def druglikeness_flags(d: dict):
    lipinski_pass = (d["MolWt"] <= 500) and (d["LogP"] <= 5) and (d["HBD"] <= 5) and (d["HBA"] <= 10)
    veber_pass = (d["RotatableBonds"] <= 10) and (d["TPSA"] <= 140)
    cns_mpo_friendly = (d["MolWt"] <= 450) and (d["TPSA"] <= 90) and (d["LogP"] <= 5)
    return lipinski_pass, veber_pass, cns_mpo_friendly


def verdict(pIC50: float, lipinski_pass: bool, veber_pass: bool, cns_friendly: bool) -> str:
    if pIC50 >= 7 and lipinski_pass and veber_pass and cns_friendly:
        return "Proceed to wet lab"
    if pIC50 >= 6 and (lipinski_pass or veber_pass):
        return "Borderline — review manually"
    return "Deprioritize"


def merge_into_shared_table(new_df: pd.DataFrame, key: str = "SMILES"):
    existing = st.session_state["compound_table"]
    if existing is None or existing.empty:
        st.session_state["compound_table"] = new_df
    else:
        st.session_state["compound_table"] = existing.merge(new_df, on=key, how="outer")


st.title("Alzheimer's (AChE) in-silico screening")
st.caption(
    "Screen, dock, and review MD stability for the same compound list — "
    "each tab adds columns to one shared results table."
)

tab_screen, tab_dock, tab_md = st.tabs(["Screening", "Docking", "MD summary"])


# ---------- Tab 1: Screening ----------

with tab_screen:
    st.subheader("Potency and druglikeness screening")

    try:
        model, feature_order = load_model()
    except FileNotFoundError:
        model, feature_order = None, None
        st.error(
            "model.pkl / feature_names.json not found. Run 02_model_training_validation.ipynb "
            "first and place both files in this same folder."
        )

    smiles_text = st.text_area(
        "SMILES (one per line)",
        height=150,
        placeholder="CCOc1ccc2nc(S(N)(=O)=O)sc2c1\nCC(=O)Nc1ccc(O)cc1",
        key="screen_smiles",
    )

    if st.button("Run screen", type="primary", disabled=model is None) and smiles_text.strip():
        rows = []
        for smi in [s.strip() for s in smiles_text.splitlines() if s.strip()]:
            d = compute_descriptors(smi)
            if d is None:
                rows.append({"SMILES": smi, "error": "invalid SMILES"})
                continue
            x = np.array([[d[f] for f in feature_order]])
            pIC50 = float(model.predict(x)[0])
            ic50_nM = 10 ** (9 - pIC50)
            lipinski_pass, veber_pass, cns_friendly = druglikeness_flags(d)
            rows.append(
                {
                    "SMILES": smi,
                    "predicted_pIC50": round(pIC50, 2),
                    "predicted_IC50_nM": round(ic50_nM, 1),
                    "MolWt": round(d["MolWt"], 1),
                    "LogP": round(d["LogP"], 2),
                    "TPSA": round(d["TPSA"], 1),
                    "Lipinski_pass": lipinski_pass,
                    "Veber_pass": veber_pass,
                    "CNS_friendly": cns_friendly,
                    "Verdict": verdict(pIC50, lipinski_pass, veber_pass, cns_friendly),
                }
            )

        results = pd.DataFrame(rows)
        merge_into_shared_table(results)
        st.dataframe(results, use_container_width=True)


# ---------- Tab 2: Docking ----------

def smiles_to_ligand_pdbqt(smiles: str):
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol)

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    if not mol_setups:
        return None
    pdbqt_string, is_ok, _err = PDBQTWriterLegacy.write_string(mol_setups[0])
    return pdbqt_string if is_ok else None


def dock_compound(smiles: str, receptor_path: str, box_center, box_size, exhaustiveness: int = 8):
    from vina import Vina

    ligand_pdbqt = smiles_to_ligand_pdbqt(smiles)
    if ligand_pdbqt is None:
        return None

    v = Vina(sf_name="vina")
    v.set_receptor(receptor_path)
    v.set_ligand_from_string(ligand_pdbqt)
    v.compute_vina_maps(center=box_center, box_size=box_size)
    v.dock(exhaustiveness=exhaustiveness, n_poses=5)
    energies = v.energies(n_poses=1)
    return float(energies[0][0])  # kcal/mol, more negative = stronger predicted binding


with tab_dock:
    st.subheader("Docking (AutoDock Vina)")
    st.caption(
        "Docks each compound against a prepared receptor. Generate receptor.pdbqt "
        "and box_config.json once with 03_docking_receptor_prep.ipynb and commit "
        "them alongside app.py, or upload your own below."
    )

    receptor_file = st.file_uploader("Receptor (.pdbqt)", type=["pdbqt"])
    box_file = st.file_uploader("Box config (.json)", type=["json"])

    receptor_path = None
    if receptor_file is not None:
        receptor_path = "uploaded_receptor.pdbqt"
        with open(receptor_path, "wb") as f:
            f.write(receptor_file.read())
    elif os.path.exists("receptor.pdbqt"):
        receptor_path = "receptor.pdbqt"

    box_cfg = None
    if box_file is not None:
        box_cfg = json.load(box_file)
    elif os.path.exists("box_config.json"):
        box_cfg = json.load(open("box_config.json"))

    dock_smiles_text = st.text_area(
        "SMILES to dock (defaults to the Screening tab's list if left blank)",
        height=100,
        key="dock_smiles",
    )

    if st.button("Dock compounds", type="primary"):
        if receptor_path is None or box_cfg is None:
            st.error(
                "Missing receptor.pdbqt or box_config.json. Run "
                "03_docking_receptor_prep.ipynb and upload the two files above, "
                "or commit them to this app's folder."
            )
        else:
            box_center = [box_cfg["center_x"], box_cfg["center_y"], box_cfg["center_z"]]
            box_size = [box_cfg["size_x"], box_cfg["size_y"], box_cfg["size_z"]]

            smiles_list = [s.strip() for s in dock_smiles_text.splitlines() if s.strip()]
            if not smiles_list:
                table = st.session_state["compound_table"]
                if table is not None and not table.empty and "SMILES" in table.columns:
                    smiles_list = table["SMILES"].dropna().unique().tolist()

            if not smiles_list:
                st.warning("No compounds to dock -- paste SMILES here or run the Screening tab first.")
            else:
                results = []
                progress = st.progress(0.0)
                for i, smi in enumerate(smiles_list):
                    score = dock_compound(smi, receptor_path, box_center, box_size)
                    results.append({"SMILES": smi, "docking_score_kcal_mol": score})
                    progress.progress((i + 1) / len(smiles_list))

                dock_df = pd.DataFrame(results)
                merge_into_shared_table(dock_df)
                st.dataframe(dock_df, use_container_width=True)


# ---------- Tab 3: MD summary ----------

with tab_md:
    st.subheader("MD stability summary")
    st.caption(
        "Full molecular dynamics is too compute-heavy to run live on free hosting. "
        "Run 04_md_stability_summary.ipynb in Colab on your top hits (a fast "
        "conformational-rigidity proxy, not a substitute for real complex MD), "
        "then upload its output CSV here to merge it into the shared table."
    )

    md_file = st.file_uploader("MD summary CSV", type=["csv"], key="md_upload")
    if md_file is not None:
        md_df = pd.read_csv(md_file)
        if "SMILES" not in md_df.columns:
            st.error("This CSV needs a SMILES column to merge -- check 04_md_stability_summary.ipynb's output.")
        else:
            merge_into_shared_table(md_df)
            st.dataframe(md_df, use_container_width=True)
            st.success("Merged into the shared results table below.")


# ---------- Shared results table (all tabs feed this) ----------

st.divider()
st.subheader("Combined results")
combined = st.session_state["compound_table"]
if combined is not None and not combined.empty:
    st.dataframe(combined, use_container_width=True)
    st.download_button(
        "Download combined results as CSV",
        data=combined.to_csv(index=False).encode("utf-8"),
        file_name="combined_results.csv",
        mime="text/csv",
    )
else:
    st.info("Run the Screening or Docking tab to start building the shared results table.")

with st.expander("What this app does and doesn't cover"):
    st.markdown(
        """
- **Screening:** predicted AChE potency (QSAR model) plus druglikeness and
  CNS-permeability flags.
- **Docking:** real AutoDock Vina docking against a prepared receptor --
  computed live, in-browser, no other tools needed.
- **MD summary:** a lightweight ligand-rigidity proxy computed offline in
  Colab and merged in here. True protein-ligand molecular dynamics (explicit
  solvent, nanosecond-scale production runs) is intentionally out of scope
  for a free, hosted app -- treat this tab's verdict as a first-pass filter,
  not a final answer.
"""
    )
