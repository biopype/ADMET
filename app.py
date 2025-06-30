import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors

st.set_page_config(page_title="Drug-Likeness & ADMET Predictor", page_icon="💊")

st.title("💊 Drug-Likeness & ADMET Predictor")
st.markdown("""
Enter a **SMILES** string to evaluate drug-likeness based on:
- Lipinski's Rule of 5
- Basic physicochemical properties
""")

smiles = st.text_input("Enter SMILES string (e.g., CC(=O)Oc1ccccc1C(=O)O):")

if smiles:
    try:
        mol = Chem.MolFromSmiles(smiles)

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        st.subheader("Molecular Properties:")
        st.write(f"**Molecular Weight**: {mw:.2f} g/mol")
        st.write(f"**LogP (lipophilicity)**: {logp:.2f}")
        st.write(f"**H-bond Donors**: {hbd}")
        st.write(f"**H-bond Acceptors**: {hba}")

        lipinski_pass = mw < 500 and logp < 5 and hbd <= 5 and hba <= 10

        st.subheader("Lipinski's Rule of 5:")
        if lipinski_pass:
            st.success("✅ Passes Lipinski's Rule")
        else:
            st.error("❌ Fails Lipinski's Rule")

    except:
        st.error("Invalid SMILES string. Please check the input.")
