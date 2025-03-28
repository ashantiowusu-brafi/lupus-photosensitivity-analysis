# app.py
import streamlit as st
import pandas as pd
import plotly.express as px
from gseapy import enrichr
from utils.analysis import run_differential_expression

st.set_page_config(page_title="Lupus Photosensitivity Explorer", layout="wide")

# Load preprocessed data (subset for speed)
@st.cache_data
def load_data():
    return pd.read_csv("data/processed/lupus_expression_subset.csv")

df = load_data()

# Sidebar filters
st.sidebar.header("Filters")
group1 = st.sidebar.selectbox("Group 1", options=df['group'].unique())
group2 = st.sidebar.selectbox("Group 2", options=df['group'].unique())

# Run differential expression analysis
if st.sidebar.button("Analyze"):
    deg_results = run_differential_expression(df, group1, group2)
    
    # Show volcano plot
    fig = px.scatter(
        deg_results, 
        x="log2_fold_change", 
        y="-log10(padj)", 
        hover_name="gene",
        title=f"DEGs: {group1} vs {group2}"
    )
    st.plotly_chart(fig, use_container_width=True)
    
    # Pathway enrichment
    enriched = enrichr(
        gene_list=deg_results[deg_results["padj"] < 0.05]["gene"].tolist(),
        gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023']
    )
    st.dataframe(enriched.results)