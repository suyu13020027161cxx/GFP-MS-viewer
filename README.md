# 🔬 Proteomics Viewer

An interactive web app for MS proteomics data visualization.

## Features
- 🌋 **Volcano plots** — with gene search & highlight, fold-change/p-value filters, SVG export
- ⭕ **Venn diagrams** — up to 5 sets, custom colors, SVG export
- 🔥 **Heatmaps** — top-N by variance, gene highlighting, multiple color scales
- 📋 **Protein table** — search, filter, CSV download

---

## Deploy to Streamlit Cloud (free, shareable link)

1. **Push to GitHub**
   ```
   git init
   git add app.py requirements.txt
   git commit -m "Initial commit"
   gh repo create proteomics-viewer --public --push
   ```

2. **Deploy on Streamlit Cloud**
   - Go to https://share.streamlit.io
   - Sign in with GitHub
   - Click **New app** → select your repo → `app.py`
   - Click **Deploy** → you get a public shareable link!

---

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

Then open http://localhost:8501 in your browser.

---

## Data format

Upload an **Excel (.xlsx)** file where:
- Each **sheet** = one condition/dataset (e.g. LC3A, LC3B, GBRP...)
- Columns include: **Gene name**, **log2 fold-change**, **p-value**, and any intensity columns
