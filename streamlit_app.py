import streamlit as st

# Try to import RDKit with proper error handling
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    st.error("⚠️ RDKit is not available. Please install it or use the alternative approach below.")

st.set_page_config(page_title="Drug-Likeness & ADMET Predictor", page_icon="💊")
st.title("💊 Drug-Likeness & ADMET Predictor")

if not RDKIT_AVAILABLE:
    st.markdown("""
    ## RDKit Installation Issue
    
    RDKit is required for this application but couldn't be imported. This is common in cloud deployments.
    
    **Solutions:**
    1. Create a `requirements.txt` file with: `rdkit-pypi`
    2. Or use conda with `environment.yml` file
    3. Consider using alternative chemistry libraries
    """)
    
    st.info("The app will still show the interface below, but calculations won't work without RDKit.")

st.markdown("""
Enter a **SMILES** string to evaluate drug-likeness based on:
- Lipinski's Rule of 5
- Basic physicochemical properties
""")

smiles = st.text_input("Enter SMILES string (e.g., CC(=O)Oc1ccccc1C(=O)O):")

if smiles:
    if not RDKIT_AVAILABLE:
        st.error("Cannot process SMILES: RDKit is not available.")
        st.info("Please check the installation instructions above.")
        st.stop()
    
    try:
        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            st.error("Invalid SMILES string. Please check the input.")
            st.stop()
        
        # Calculate molecular properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        # Display results
        st.subheader("Molecular Properties:")
        
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Molecular Weight", f"{mw:.2f} g/mol")
            st.metric("H-bond Donors", hbd)
        
        with col2:
            st.metric("LogP (lipophilicity)", f"{logp:.2f}")
            st.metric("H-bond Acceptors", hba)
        
        # Lipinski's Rule of 5 evaluation
        st.subheader("Lipinski's Rule of 5:")
        
        # Individual rule checks
        rules = {
            "Molecular Weight < 500 Da": (mw < 500, mw),
            "LogP < 5": (logp < 5, logp),
            "H-bond Donors ≤ 5": (hbd <= 5, hbd),
            "H-bond Acceptors ≤ 10": (hba <= 10, hba)
        }
        
        passed_rules = 0
        for rule, (passed, value) in rules.items():
            if passed:
                st.success(f"✅ {rule} (Value: {value:.2f})")
                passed_rules += 1
            else:
                st.error(f"❌ {rule} (Value: {value:.2f})")
        
        # Overall assessment
        st.subheader("Overall Assessment:")
        if passed_rules == 4:
            st.success("🎉 **Drug-like**: Passes all Lipinski's Rules!")
        elif passed_rules >= 3:
            st.warning(f"⚠️ **Borderline**: Passes {passed_rules}/4 rules")
        else:
            st.error(f"❌ **Not drug-like**: Only passes {passed_rules}/4 rules")
            
    except Exception as e:
        st.error(f"Error processing SMILES: {str(e)}")
        st.info("Please check that you've entered a valid SMILES string.")

# Add some example SMILES for testing
st.subheader("Example SMILES to try:")
examples = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1"
}

cols = st.columns(len(examples))
for i, (name, smile) in enumerate(examples.items()):
    with cols[i]:
        if st.button(f"{name}", key=f"btn_{i}"):
            st.rerun()