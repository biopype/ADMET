# 💊 Drug-Likeness & ADMET Predictor

A simple Streamlit web app to evaluate molecular drug-likeness using SMILES strings based on:

- Lipinski's Rule of 5
- Molecular descriptors (MW, LogP, HBD, HBA)

## 🚀 Try it Live
[Live App Link](https://your-username.streamlit.app) ← replace after deployment

## 🧪 Example SMILES
- Aspirin: `CC(=O)Oc1ccccc1C(=O)O`
- Ibuprofen: `CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O`

## 🛠 Run Locally

```bash
git clone https://github.com/your-username/drug-likeness-app.git
cd drug-likeness-app
pip install -r requirements.txt
streamlit run app.py