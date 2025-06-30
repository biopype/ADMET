import streamlit as st
import requests

st.set_page_config(page_title="Drug-Likeness & ADMET Predictor", page_icon="💊")

st.title("💊 Drug-Likeness & ADMET Predictor (API-Based)")
st.markdown("""
Enter a **SMILES** string to evaluate drug-likeness using:
- Lipinski's Rule of 5 (calculated manually)
- ADMET properties via **ADMETlab 2.0 API**
""")

smiles = st.text_input("Enter SMILES string:")

def get_admetlab_data(smiles):
    url = "https://admetmesh.scbdd.com/service/predict"
    payload = {"smiles": smiles}
    try:
        response = requests.post(url, json=payload)
        if response.status_code == 200:
            return response.json()
        else:
            return None
    except:
        return None

if smiles:
    with st.spinner("Fetching data from ADMETlab..."):
        data = get_admetlab_data(smiles)

    if data and 'data' in data:
        props = data['data'][0]

        mw = float(props.get("MW", 0))
        logp = float(props.get("MLOGP", 0))
        hbd = int(props.get("nHDon", 0))
        hba = int(props.get("nHAcc", 0))

        st.subheader("Molecular Properties")
        st.write(f"**Molecular Weight**: {mw} g/mol")
        st.write(f"**LogP (MLogP)**: {logp}")
        st.write(f"**H-bond Donors**: {hbd}")
        st.write(f"**H-bond Acceptors**: {hba}")

        st.subheader("Lipinski’s Rule of 5")
        lipinski = mw < 500 and logp < 5 and hbd <= 5 and hba <= 10
        if lipinski:
            st.success("✅ Passes Lipinski's Rule")
        else:
            st.error("❌ Fails Lipinski's Rule")

        st.subheader("Selected ADMET Properties")
        st.write(f"**Human Intestinal Absorption**: {props.get('HIA', 'N/A')}")
        st.write(f"**CYP1A2 Inhibitor**: {props.get('CYP1A2_inh', 'N/A')}")
        st.write(f"**hERG Blocker**: {props.get('hERG', 'N/A')}")
        st.write(f"**Ames Toxicity**: {props.get('AMES', 'N/A')}")
        st.write(f"**Carcinogenicity**: {props.get('Carcinogens', 'N/A')}")
        st.write(f"**LD50 (rat, oral)**: {props.get('LD50', 'N/A')}")

    else:
        st.error("Failed to retrieve data. Make sure the SMILES is valid.")
