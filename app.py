import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from venn import venn
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import io
import base64

# ─── Page config ────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Proteomics Viewer",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─── Custom CSS — force full light theme regardless of Streamlit settings ────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600&display=swap');

/* Force white background on every Streamlit container */
html, body                          { background: #ffffff !important; color: #1a1a2e !important; }
.stApp                              { background: #ffffff !important; color: #1a1a2e !important; }
.stApp > header                     { background: #ffffff !important; }
.main .block-container              { background: #ffffff !important; padding-top: 2rem; }
section[data-testid="stSidebar"]    { background: #f5f7fa !important; }
section[data-testid="stSidebar"] *  { color: #1a1a2e !important; }

/* Widgets */
div[data-baseweb="select"] > div,
div[data-baseweb="input"] > div,
div[data-baseweb="textarea"] > div  { background: #ffffff !important; color: #1a1a2e !important; border-color: #dde1f0 !important; }
.stSlider label, .stCheckbox label,
.stSelectbox label, .stTextArea label,
.stTextInput label, .stMultiSelect label,
.stRadio label                      { color: #1a1a2e !important; }
.stMarkdown, .stText, p, li         { color: #1a1a2e !important; }

/* Typography */
h1, h2, h3 {
    font-family: 'IBM Plex Mono', monospace !important;
    letter-spacing: -0.5px;
    color: #1a1a2e !important;
}
html, body, [class*="css"] { font-family: 'IBM Plex Sans', sans-serif; }

/* Custom components */
.metric-card {
    background: #f5f7fa;
    border: 1px solid #dde1f0;
    border-radius: 8px;
    padding: 1rem 1.5rem;
    text-align: center;
}
.metric-card .value {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 2rem;
    font-weight: 600;
    color: #3a6bc4;
}
.metric-card .label {
    font-size: 0.8rem;
    color: #888;
    text-transform: uppercase;
    letter-spacing: 1px;
}
.section-header {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.75rem;
    text-transform: uppercase;
    letter-spacing: 2px;
    color: #3a6bc4;
    border-bottom: 1px solid #dde1f0;
    padding-bottom: 0.5rem;
    margin-bottom: 1rem;
}
</style>
""", unsafe_allow_html=True)

# ─── Sidebar ─────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## 🔬 Proteomics Viewer")
    st.markdown("---")
    page = st.radio(
        "Navigation",
        ["🌋 Volcano Plot", "⭕ Venn Diagram", "🔥 Heatmap", "📋 Protein Table"],
        label_visibility="collapsed"
    )
    st.markdown("---")
    st.markdown("### Upload Data")
    uploaded_file = st.file_uploader(
        "Upload Excel file (.xlsx)",
        type=["xlsx"],
        help="Each sheet = one dataset/condition"
    )
    st.markdown("---")
    st.caption("Proteomics Viewer v1.0")

# ─── Helpers ─────────────────────────────────────────────────────────────────
COLORS = {
    'LC3a':   '#9B59B6',
    'LC3b':   '#E8A0C8',
    'GBRP':   '#B5A800',
    'GBRPL1': '#2D7D2D',
    'GBRPL2': '#00AAAA',
    'default': '#7eb8f7',
}

def fig_to_bytes(fig, fmt="svg"):
    buf = io.BytesIO()
    fig.write_image(buf, format=fmt, scale=3)
    buf.seek(0)
    return buf.read()

def mpl_to_bytes(fig, fmt="svg"):
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, bbox_inches="tight", dpi=300, facecolor=fig.get_facecolor())
    buf.seek(0)
    return buf.read()

def load_sheets(file):
    xl = pd.ExcelFile(file)
    return xl.sheet_names, {s: xl.parse(s) for s in xl.sheet_names}

def guess_col(cols, keywords):
    """Return the first column whose name contains any keyword (case-insensitive)."""
    for kw in keywords:
        for c in cols:
            if kw.lower() in c.lower():
                return c
    return cols[0]

# ─── No file fallback ────────────────────────────────────────────────────────
def no_file_warning():
    st.info("📂 Upload an Excel file in the sidebar to get started.", icon="📂")

# ═══════════════════════════════════════════════════════════════════════════════
# PAGE: VOLCANO PLOT
# ═══════════════════════════════════════════════════════════════════════════════
if page == "🌋 Volcano Plot":
    st.markdown("## Volcano Plot")
    st.markdown('<p class="section-header">Differential expression analysis</p>', unsafe_allow_html=True)

    if not uploaded_file:
        no_file_warning()
    else:
        sheet_names, sheets = load_sheets(uploaded_file)
        col_settings, col_plot = st.columns([1, 3])

        with col_settings:
            # ── Dataset & columns ──────────────────────────────────────────
            st.markdown("**Dataset**")
            sheet = st.selectbox("Select dataset", sheet_names)
            df    = sheets[sheet].copy()
            cols  = df.columns.tolist()

            # Smart auto-detection of columns
            gene_col  = st.selectbox("Gene / Protein column",
                                     cols,
                                     index=cols.index(guess_col(cols, ["gene","protein","name","id","accession"])))
            fc_col    = st.selectbox("Fold-change column (log2FC)",
                                     cols,
                                     index=cols.index(guess_col(cols, ["logfc","log2fc","fc","fold","ratio","diff"])))
            pval_col  = st.selectbox("P-value column",
                                     cols,
                                     index=cols.index(guess_col(cols, ["pval","p.val","p_val","p-val","adj","qval","fdr","pvalue"])))
            use_log_p = st.checkbox("Apply –log10 to p-value", value=True,
                                    help="Check if your p-value column contains raw p-values (not already –log10 transformed)")

            # ── Thresholds ─────────────────────────────────────────────────
            st.markdown("**Thresholds**")
            fc_thresh   = st.slider("Fold-change cutoff (|log2FC|)", 0.0, 5.0, 1.0, 0.1)
            pval_thresh = st.slider("p-value cutoff", 0.001, 0.1, 0.05, 0.001, format="%.3f")

            # ── Dot colors ─────────────────────────────────────────────────
            st.markdown("**Dot colors**")
            c1, c2, c3, c4 = st.columns(4)
            with c1: color_up   = st.color_picker("Up",   "#e05c5c")
            with c2: color_down = st.color_picker("Down", "#5c9ee0")
            with c3: color_ns   = st.color_picker("NS",   "#888899")
            with c4: color_hi   = st.color_picker("Hit",  "#f5c518")

            dot_size    = st.slider("Dot size", 2, 15, 5)
            dot_opacity = st.slider("Dot opacity", 0.1, 1.0, 0.75, 0.05)

            # ── Axis ranges ────────────────────────────────────────────────
            st.markdown("**Axis ranges**")
            auto_axes = st.checkbox("Auto axes", value=True)
            if not auto_axes:
                ax_c1, ax_c2 = st.columns(2)
                with ax_c1:
                    x_min = st.number_input("X min", value=-10.0, step=0.5)
                    y_min = st.number_input("Y min", value=0.0,   step=0.5)
                with ax_c2:
                    x_max = st.number_input("X max", value=10.0,  step=0.5)
                    y_max = st.number_input("Y max", value=5.0,   step=0.5)

            # ── Gene search with categories ────────────────────────────────
            st.markdown("**Gene search & highlight by category**")
            st.caption("Format: `GENE_NAME:CATEGORY` or just `GENE_NAME`")
            st.caption("Example: `TP53:tumor suppressors` or `MDM2:regulators`")
            
            search_input   = st.text_area("Genes to highlight (one per line)",
                                          placeholder="TP53:tumor suppressors\nMDM2:regulators\nBRCA1:tumor suppressors", 
                                          height=100)
            
            # Parse gene input with categories
            highlight_genes_dict = {}  # {gene: category}
            for line in search_input.split("\n"):
                line = line.strip()
                if line:
                    if ":" in line:
                        gene, category = line.split(":", 1)
                        highlight_genes_dict[gene.strip()] = category.strip()
                    else:
                        highlight_genes_dict[line] = "Highlighted"
            
            # Get unique categories and assign colors
            categories = list(set(highlight_genes_dict.values()))
            category_colors = {}
            default_category_colors = ["#e05c5c", "#5c9ee0", "#f5c518", "#2dd4bf", "#8b5cf6", "#ec4899", "#f97316"]
            for i, cat in enumerate(sorted(categories)):
                category_colors[cat] = default_category_colors[i % len(default_category_colors)]
            
            # Show category colors
            if categories:
                st.markdown("**Category colors:**")
                for cat in sorted(categories):
                    st.markdown(f"<span style='color:{category_colors[cat]};'>●</span> {cat}", 
                              unsafe_allow_html=True)
            
            label_hits      = st.checkbox("Show gene labels on plot", value=True)

        # ── Plot ───────────────────────────────────────────────────────────
        with col_plot:
            try:
                df = df.copy()
                df[fc_col]   = pd.to_numeric(df[fc_col],   errors="coerce")
                df[pval_col] = pd.to_numeric(df[pval_col], errors="coerce")
                df = df.dropna(subset=[fc_col, pval_col])

                if use_log_p:
                    # guard against zero / negative p-values
                    df = df[df[pval_col] > 0]
                    df["_y"] = -np.log10(df[pval_col])
                    y_label  = "–log10(p-value)"
                    sig_p    = -np.log10(pval_thresh)
                else:
                    df["_y"] = df[pval_col]
                    y_label  = "p-value"
                    sig_p    = pval_thresh

                def classify(fc, y):
                    if abs(fc) >= fc_thresh and y >= sig_p:
                        return "Significant Up" if fc > 0 else "Significant Down"
                    return "Not significant"

                df["_class"] = df.apply(lambda r: classify(r[fc_col], r["_y"]), axis=1)

                color_map = {
                    "Significant Up":   color_up,
                    "Significant Down": color_down,
                    "Not significant":  color_ns,
                }

                fig = go.Figure()

                for cls, color in color_map.items():
                    sub = df[df["_class"] == cls]
                    fig.add_trace(go.Scatter(
                        x=sub[fc_col], y=sub["_y"],
                        mode="markers",
                        name=cls,
                        marker=dict(color=color, size=dot_size, opacity=dot_opacity),
                        text=sub[gene_col],
                        hovertemplate="<b>%{text}</b><br>log2FC: %{x:.3f}<br>" + y_label + ": %{y:.3f}<extra></extra>"
                    ))

                # Highlighted genes by category (case-insensitive match, drawn on top)
                if highlight_genes_dict:
                    for category, color in category_colors.items():
                        # Find genes in this category
                        genes_in_cat = [g for g, c in highlight_genes_dict.items() if c == category]
                        genes_lower = [g.lower() for g in genes_in_cat]
                        hi = df[df[gene_col].astype(str).str.lower().isin(genes_lower)]
                        
                        if not hi.empty:
                            # Smart text positioning: alternate between top and bottom to avoid overlap
                            text_positions = []
                            for idx in range(len(hi)):
                                if idx % 2 == 0:
                                    text_positions.append("top center")
                                else:
                                    text_positions.append("bottom center")
                            
                            mode = "markers+text" if label_hits else "markers"
                            fig.add_trace(go.Scatter(
                                x=hi[fc_col], y=hi["_y"],
                                mode=mode,
                                name=category,
                                marker=dict(color=color, size=dot_size + 6, symbol="circle",
                                            line=dict(color="white", width=1.5)),
                                text=hi[gene_col],
                                textposition=text_positions[0] if len(text_positions) == 1 else "top center",
                                textfont=dict(color=color, size=11, family="Arial"),
                                hovertemplate="<b>%{text}</b><br>" + category + "<br>log2FC: %{x:.3f}<br>" + y_label + ": %{y:.3f}<extra></extra>"
                            ))
                    
                    if len(highlight_genes_dict) > 0 and len(df[df[gene_col].astype(str).str.lower().isin([g.lower() for g in highlight_genes_dict.keys()])]) == 0:
                        st.warning("No matching genes found — check spelling or gene column selection.")

                # Threshold lines
                fig.add_hline(y=sig_p,        line_dash="dash", line_color="#aaa", opacity=0.5)
                fig.add_vline(x= fc_thresh,   line_dash="dash", line_color="#aaa", opacity=0.5)
                fig.add_vline(x=-fc_thresh,   line_dash="dash", line_color="#aaa", opacity=0.5)

                layout_kwargs = dict(
                    template="plotly_white",
                    paper_bgcolor="#ffffff",
                    plot_bgcolor="#f5f7fa",
                    font=dict(family="Arial", color="#1a1a2e"),
                    xaxis_title="log2 Fold Change",
                    yaxis_title=y_label,
                    legend=dict(bgcolor="#ffffff", bordercolor="#dde1f0", borderwidth=1),
                    height=580,
                    margin=dict(l=50, r=20, t=30, b=50),
                )
                if not auto_axes:
                    layout_kwargs["xaxis"] = dict(range=[x_min, x_max])
                    layout_kwargs["yaxis"] = dict(range=[y_min, y_max])

                fig.update_layout(**layout_kwargs)
                st.plotly_chart(fig, use_container_width=True)

                # ── Summary metrics ────────────────────────────────────────
                m1, m2, m3, m4 = st.columns(4)
                up   = (df["_class"] == "Significant Up").sum()
                down = (df["_class"] == "Significant Down").sum()
                total_hi = len(df[df[gene_col].astype(str).str.lower().isin([g.lower() for g in highlight_genes_dict.keys()])]) if highlight_genes_dict else 0
                
                for col_, val, lbl in zip(
                    [m1, m2, m3, m4],
                    [len(df), up, down, total_hi],
                    ["Total proteins", "Up-regulated", "Down-regulated", "Highlighted"]
                ):
                    col_.markdown(f"""
                    <div class="metric-card">
                        <div class="value">{val}</div>
                        <div class="label">{lbl}</div>
                    </div>""", unsafe_allow_html=True)

                # ── Export buttons ─────────────────────────────────────────
                st.markdown("####")
                st.markdown("**Export plot**")
                try:
                    ex1, ex2, ex3 = st.columns(3)
                    with ex1:
                        st.download_button("⬇ SVG",
                            data=fig_to_bytes(fig, "svg"),
                            file_name="volcano.svg", mime="image/svg+xml")
                    with ex2:
                        st.download_button("⬇ PNG",
                            data=fig_to_bytes(fig, "png"),
                            file_name="volcano.png", mime="image/png")
                    with ex3:
                        st.download_button("⬇ TIFF",
                            data=fig_to_bytes(fig, "pdf"),   # kaleido doesn't do tiff; PDF is vector like tiff
                            file_name="volcano.pdf", mime="application/pdf")
                except Exception:
                    st.caption("Install `kaleido` (`pip install kaleido`) to enable image export.")

                # ── Highlighted gene table ─────────────────────────────────
                if highlight_genes_dict:
                    hi_all = pd.DataFrame()
                    for gene, category in highlight_genes_dict.items():
                        hi_gene = df[df[gene_col].astype(str).str.lower() == gene.lower()]
                        if not hi_gene.empty:
                            hi_gene = hi_gene.copy()
                            hi_gene["Gene Category"] = category
                            hi_all = pd.concat([hi_all, hi_gene], ignore_index=True)
                    
                    if not hi_all.empty:
                        st.markdown("**Highlighted gene details**")
                        st.dataframe(hi_all[[gene_col, fc_col, pval_col, "_class", "Gene Category"]].rename(
                            columns={"_class": "Expression Class"}
                        ), use_container_width=True)

            except Exception as e:
                st.error(f"Error generating plot: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# PAGE: VENN DIAGRAM
# ═══════════════════════════════════════════════════════════════════════════════
elif page == "⭕ Venn Diagram":
    st.markdown("## Venn Diagram")
    st.markdown('<p class="section-header">Set overlap analysis</p>', unsafe_allow_html=True)

    if not uploaded_file:
        no_file_warning()
    else:
        sheet_names, sheets = load_sheets(uploaded_file)

        col_s, col_p = st.columns([1, 3])
        with col_s:
            st.markdown("**Settings**")
            selected_sheets = st.multiselect(
                "Select 2–5 datasets", sheet_names,
                default=sheet_names[:min(5, len(sheet_names))]
            )
            if selected_sheets:
                sample_df = sheets[selected_sheets[0]]
                gene_col  = st.selectbox("Gene/protein column", sample_df.columns.tolist())

                st.markdown("**Colors**")
                color_list = []
                defaults = ['#9B59B6','#E8A0C8','#B5A800','#2D7D2D','#00AAAA']
                for i, s in enumerate(selected_sheets):
                    c = st.color_picker(s, defaults[i % len(defaults)])
                    color_list.append(c)

        with col_p:
            if selected_sheets and 2 <= len(selected_sheets) <= 5:
                try:
                    sets = {s: set(sheets[s][gene_col].dropna()) for s in selected_sheets}

                    fig_v, ax = plt.subplots(figsize=(8, 8))
                    fig_v.patch.set_facecolor("#ffffff")
                    ax.set_facecolor("#ffffff")

                    v = venn(sets, ax=ax, fmt="{size}")

                    for i, patch in enumerate(ax.patches):
                        patch.set_facecolor(color_list[i % len(color_list)])
                        patch.set_alpha(0.5)

                    for text in ax.texts:
                        text.set_color("#1a1a2e")
                        text.set_fontsize(10)

                    ax.set_title("Overlapping proteins", color="#1a1a2e",
                                 fontsize=14, fontfamily="monospace", pad=15)

                    st.pyplot(fig_v)
                    ve1, ve2, ve3 = st.columns(3)
                    with ve1:
                        st.download_button("⬇ SVG",  mpl_to_bytes(fig_v, "svg"),  "venn.svg",  "image/svg+xml")
                    with ve2:
                        st.download_button("⬇ PNG",  mpl_to_bytes(fig_v, "png"),  "venn.png",  "image/png")
                    with ve3:
                        st.download_button("⬇ TIFF", mpl_to_bytes(fig_v, "tiff"), "venn.tiff", "image/tiff")

                    # Overlap table
                    st.markdown("**Unique to each set**")
                    rows = []
                    for s in selected_sheets:
                        others = set().union(*[sets[o] for o in selected_sheets if o != s])
                        unique = sets[s] - others
                        rows.append({"Dataset": s, "Unique proteins": len(unique),
                                     "Examples": ", ".join(list(unique)[:5])})
                    st.dataframe(pd.DataFrame(rows), use_container_width=True)

                except Exception as e:
                    st.error(f"Error: {e}")
            elif selected_sheets:
                st.warning("Please select between 2 and 5 datasets.")


# ═══════════════════════════════════════════════════════════════════════════════
# PAGE: HEATMAP
# ═══════════════════════════════════════════════════════════════════════════════
elif page == "🔥 Heatmap":
    st.markdown("## Heatmap")
    st.markdown('<p class="section-header">Expression across conditions</p>', unsafe_allow_html=True)

    if not uploaded_file:
        no_file_warning()
    else:
        sheet_names, sheets = load_sheets(uploaded_file)

        col_s, col_p = st.columns([1, 3])
        with col_s:
            st.markdown("**Settings**")
            sheet    = st.selectbox("Select sheet", sheet_names)
            df       = sheets[sheet].copy()
            gene_col = st.selectbox("Gene column", df.columns.tolist())
            num_cols = df.select_dtypes(include=np.number).columns.tolist()
            val_cols = st.multiselect("Value columns", num_cols, default=num_cols[:min(6, len(num_cols))])
            top_n    = st.slider("Top N proteins (by variance)", 10, 200, 50)
            colorscale = st.selectbox("Color scale", ["RdBu_r", "Viridis", "Plasma", "Cividis", "Turbo"])

            search_hm = st.text_area("Highlight genes", placeholder="TP53\nBRCA1", height=100)
            highlight_hm = [g.strip() for g in search_hm.split("\n") if g.strip()]

        with col_p:
            if val_cols:
                try:
                    hm = df[[gene_col] + val_cols].dropna()
                    hm = hm.set_index(gene_col)

                    if highlight_hm:
                        hi_df  = hm[hm.index.isin(highlight_hm)]
                        rest   = hm[~hm.index.isin(highlight_hm)]
                        rest   = rest.loc[rest.var(axis=1).nlargest(top_n).index]
                        hm     = pd.concat([hi_df, rest])
                    else:
                        hm = hm.loc[hm.var(axis=1).nlargest(top_n).index]

                    fig = px.imshow(
                        hm,
                        color_continuous_scale=colorscale,
                        aspect="auto",
                        labels=dict(color="Intensity"),
                    )
                    fig.update_layout(
                        template="plotly_white",
                        paper_bgcolor="#ffffff",
                        plot_bgcolor="#f5f7fa",
                        font=dict(family="Arial", color="#1a1a2e"),
                        height=650,
                        margin=dict(l=120, r=20, t=30, b=40),
                    )
                    st.plotly_chart(fig, use_container_width=True)

                    try:
                        he1, he2, he3 = st.columns(3)
                        with he1:
                            st.download_button("⬇ SVG", fig_to_bytes(fig, "svg"), "heatmap.svg", "image/svg+xml")
                        with he2:
                            st.download_button("⬇ PNG", fig_to_bytes(fig, "png"), "heatmap.png", "image/png")
                        with he3:
                            st.download_button("⬇ PDF", fig_to_bytes(fig, "pdf"), "heatmap.pdf", "application/pdf")
                    except Exception:
                        st.caption("Install `kaleido` to enable image export.")

                except Exception as e:
                    st.error(f"Error: {e}")
            else:
                st.info("Select at least one value column.")


# ═══════════════════════════════════════════════════════════════════════════════
# PAGE: PROTEIN TABLE
# ═══════════════════════════════════════════════════════════════════════════════
elif page == "📋 Protein Table":
    st.markdown("## Protein Table")
    st.markdown('<p class="section-header">Browse and filter your dataset</p>', unsafe_allow_html=True)

    if not uploaded_file:
        no_file_warning()
    else:
        sheet_names, sheets = load_sheets(uploaded_file)

        sheet = st.selectbox("Select dataset", sheet_names)
        df    = sheets[sheet].copy()

        col_f1, col_f2, col_f3 = st.columns(3)
        with col_f1:
            search_gene = st.text_input("🔍 Search gene/protein name", "")
        with col_f2:
            num_cols = df.select_dtypes(include=np.number).columns.tolist()
            if num_cols:
                fc_col_t = st.selectbox("Filter by column", ["None"] + num_cols)
        with col_f3:
            if num_cols and fc_col_t != "None":
                min_v = float(df[fc_col_t].min())
                max_v = float(df[fc_col_t].max())
                range_v = st.slider("Value range", min_v, max_v, (min_v, max_v))

        filtered = df.copy()
        if search_gene:
            mask = filtered.apply(
                lambda col: col.astype(str).str.contains(search_gene, case=False, na=False)
            ).any(axis=1)
            filtered = filtered[mask]

        if num_cols and fc_col_t != "None":
            filtered = filtered[
                (filtered[fc_col_t] >= range_v[0]) & (filtered[fc_col_t] <= range_v[1])
            ]

        st.markdown(f"**{len(filtered):,} proteins** shown")
        st.dataframe(filtered, use_container_width=True, height=500)

        csv = filtered.to_csv(index=False).encode()
        st.download_button("⬇ Download filtered table (.csv)", csv,
                           file_name="filtered_proteins.csv", mime="text/csv")
