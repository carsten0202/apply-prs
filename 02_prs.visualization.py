import os
# os.chdir(r"C:\Study\Europe\Torben\small\2020 08 CVD Target scores")
# os.chdir(r"C:\Study\Europe\Torben\small\2020 08 CVD Target scores\testing CAD")
os.chdir(r"C:\Study\Europe\Torben\small\2020 08 CVD Target scores\olink")
import numpy as np
import seaborn as sns
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
import matplotlib
import matplotlib.pyplot as plt
import ppscore as pps
import pandas as pd
sns.color_palette()
# Loading PRS scores
df = pd.DataFrame()
for fname in os.listdir("."):
    if not fname.endswith(".profile"):
        continue
    phenotype = fname.split(".")[0]
    _df = pd.read_table(fname, sep="\s+").set_index("IID")
    if df.shape[0] == 0:
        # initialize indices
        df = _df[[]].copy()
    # TODO: this will use the *first* file's index
    # this works with current pipeline, because all files are generated on the same genetic input ->
    # thus index is the same everywhere, but this may potentially cause errors
    df[phenotype] = _df["SCORESUM"] / _df["CNT2"]

# Loading phenotypes
cardiotraits = ["z_BMI", "chol_total", "chol_hdl", "chol_ldl", "triglycerides",
                "waist", "hip", "bp_sys_avg", "bp_dia_avg", "z_bp_sys", "z_bp_dia",
                "glucose", "insulin", "HOMA"]
phenotypes_df = pd.read_excel(
    r'C:\Study\Europe\Torben\small\2020 08 CVD Target scores\three_scores_without_rare_variants.xlsx').set_index("IID")
phenotypes_df = phenotypes_df[["age"] + cardiotraits]

# Merging
total_df = pd.merge(df, phenotypes_df, left_index=True, right_index=True)
print(total_df.shape)

# Selecting non-obese
samples = [
    f"66-{sid.strip()}" for sid in open(r"C:\Study\Europe\Torben\small\2020 08 CVD Target scores\HealthySamples.txt", 'r')]
# we want to exclude the following samples from the PRS selection process, to avoid overfitting
olinked_samples = [
    f"66-{sid.strip()}" for sid in open(r"C:\Study\Europe\Torben\small\2020 08 CVD Target scores\Olink_blood_sample_IDs.txt", 'r')]

translate = {
    "olink-PGS000217": "ADM_PRS",
    "olink-PGS000218": "AGRP_PRS",
    "olink-PGS000221": "CCL3_PRS",
    "olink-PGS000230": "CXCL1_PRS",
    "olink-PGS000233": "Dkk_1_PRS",
    "olink-PGS000239": "FGF_23_PRS",
    "olink-PGS000240": "FS_PRS",
    "olink-PGS000244": "GH_PRS",
    "olink-PGS000245": "HB_EGF_PRS",
    "olink-PGS000248": "HSP_27_PRS",
    "olink-PGS000249": "IL_18_PRS",
    "olink-PGS000250": "IL_1ra_PRS",
    "olink-PGS000251": "IL_27_PRS",
    "olink-PGS000252": "IL_6_PRS",
    "olink-PGS000255": "IL_16_PRS",
    "olink-PGS000256": "ITGB1BP2_PRS",
    "olink-PGS000257": "KIM_1_PRS",
    "olink-PGS000259": "LEP_PRS",
    "olink-PGS000260": "LOX_1_PRS",
    "olink-PGS000266": "MMP_12_PRS",
    "olink-PGS000268": "MMP_7_PRS",
    "olink-PGS000272": "PAPPA_PRS",
    "olink-PGS000273": "PAR_1_PRS",
    "olink-PGS000274": "PDGF_subunit_B_PRS",
    "olink-PGS000277": "PSGL_1_PRS",
    "olink-PGS000278": "PTX3_PRS",
    "olink-PGS000279": "RAGE_PRS",
    "olink-PGS000280": "REN_PRS",
    "olink-PGS000282": "SCF_PRS",
    "olink-PGS000285": "TF_PRS",
    "olink-PGS000286": "TM_PRS",
    "olink-PGS000288": "TNF_R2_PRS",
    "olink-PGS000291": "TRAIL_R2_PRS",
    "olink-PGS000295": "VEGF_D_PRS",
}
total_df.rename(columns=translate, inplace=True)

# there are 1863 population samples, for whom genetics is available - len(set(total_df.index).intersection(samples))
# of them, 198 have olink data for proteins
print("from here, use code below if you want only pop analysis, otherwise scroll to the bottom")
assert False

# Selecting non-obese - only healthy, have no olink data (to avoid an overfit)
total_df = total_df.loc[total_df.index.intersection(samples).difference(olinked_samples)]  # 1665 samples
total_df_sample = total_df
print(total_df_sample.shape)

# Using PPS
matrix_df = pps.matrix(total_df_sample)[['x', 'y', 'ppscore']].pivot(columns='x', index='y', values='ppscore')
fig = plt.figure(figsize=(16, 16))
ax = sns.heatmap(matrix_df, vmin=0, vmax=1, cmap="Blues", linewidths=0.5, annot=True)
axlim = total_df_sample.shape[1]
ax.set_ylim(0, axlim)  # a patch for the heatmap in sns - it was broken in matplotlib 3.1.1
plt.savefig("pps.png")
plt.close()
# No findings

# Correlations?
fig = plt.figure(figsize=(16, 16))
# ax = sns.heatmap(total_df_sample.corr()[["CAD-PGS000013"]].abs(), vmin=0, vmax=1, cmap="Blues", annot=True)
ax = sns.heatmap(total_df_sample.corr().abs(), vmin=0, vmax=1, cmap="Blues", annot=True)
axlim = total_df_sample.shape[1]
ax.set_ylim(0, axlim)  # a patch for the heatmap in sns - it was broken in matplotlib 3.1.1
plt.savefig("corr.png")
plt.close()

# For olink PRS
# # Correlations?
# fig = plt.figure(figsize=(32, 32))
# # ax = sns.heatmap(total_df_sample.corr()[["CAD-PGS000013"]].abs(), vmin=0, vmax=1, cmap="Blues", annot=True)
# ax = sns.heatmap(total_df_sample.corr(), vmin=-1, vmax=1, cmap="vlag", annot=True)
# axlim = total_df_sample.shape[1]
# ax.set_ylim(0, axlim)  # a patch for the heatmap in sns - it was broken in matplotlib 3.1.1
# plt.savefig("corr.png")
# plt.close()

# Analyzing correlations
studied_prs = ['CAD-PGS000012', 'CAD-PGS000013', 'CAD-PGS000018', 'CAD-PGS000116', 'CAD-PGS000296']
negative_control_prs = ['prostatecancer-PGS000333', 'psychiatristvisit-PGS000141']
prs_corrs = total_df_sample.corr().abs().loc[studied_prs+negative_control_prs, studied_prs+negative_control_prs]
print("Average cross-correlation between the PRS")
print("We recommend to exclude the PRS that are unlike others (i.e. outliers with small cross-correlation), unless they perform extremely well")
for prs in studied_prs:
    avg_corr = (prs_corrs.loc[prs, studied_prs].sum() - 1) / (len(studied_prs) - 1)  # correcting for 1.0 correlation in the table
    print(f"{prs}\t{avg_corr:.3f}")
print("Here are the negative controls cross-correlations for comparison")
for prs in negative_control_prs:
    avg_corr = prs_corrs.loc[prs, studied_prs].sum() / len(studied_prs)
    print(f"{prs}\t{avg_corr:.3f}")
check_phenotypes = cardiotraits
prs_phen_corrs = total_df_sample.corr().abs().loc[studied_prs+negative_control_prs, check_phenotypes]
avg_control_performance = prs_phen_corrs.loc[negative_control_prs].mean(axis=0)
print("Performance-wise: how much each score outperforms negative control PRS average in each phenotype")
print("We recommend picking the best performing one, if exists; "
      "otherwise pick a well-correlated with other PRS (so that it does not matter if you picked wrongly) with good performance")
print(prs_phen_corrs.loc[studied_prs] // avg_control_performance)


# -----
# This code was used to figure out that 500 samples will be enough to select a PRS to use, on average
# for N in (100, 300, 500, 1000):
#     corrs = []
#     for i in range(20):
#         total_df_sample = total_df.sample(n=N)
#         corrs.append(total_df_sample.corr().abs().values)
#     arr = np.array(corrs)
#     means = pd.DataFrame(arr.mean(axis=0), index=total_df.columns, columns=total_df.columns)
#     stds = pd.DataFrame(arr.std(axis=0), index=total_df.columns, columns=total_df.columns)
#     annotations = np.array([f"{means.values[i][j]:.3f}\n+-{stds.values[i][j]:.3f}" for i in range(total_df.shape[1]) for j in range(total_df.shape[1])]).reshape((total_df.shape[1], total_df.shape[1]))
#     fig = plt.figure(figsize=(16, 16))
#     sns.heatmap(means, vmin=0, vmax=1, cmap="Blues", annot=annotations, fmt='')
#     plt.savefig(f"corr-{N}.png")
#     plt.close()

# =======

total_df["population"] = total_df.index.isin(samples)
total_df["obese"] = ~total_df["population"]

olink_df = pd.read_table(r"C:\Study\Europe\Torben\small\2020 08 CVD Target scores\olink\CVDII_INF-panel_NPX_20210205.txt", sep="\t").set_index("blood_sample_ID")
olink_df.index = olink_df.index.map(lambda x: f"66-{x}")
correction_table = {
    # WARNING: this table only fixes inconsistencies that we needed, not all of them!
    # fix the upstream data delivery to ensure correct name matching in the future!
    # .y .x issues
    "CXCL1.y": "CXCL1",
    "SCF.y": "SCF",
    "IL18.y": "IL_18",  # subject to this and next issue
    "CCL3.y": "CCL3",
    "FGF_21.x": "FGF_21",
    # separation inconsistency (some interleukins are called "ILXX", others are "IL_XX")
    "IL10": "IL_10",
    "IL16": "IL_16",
    "IL1RL2": "IL_1RL2",
    "IL6": "IL_6",
    "IL7": "IL_7",
    "IL8": "IL_8",
    "KIM1": "KIM_1",
    "MMP12": "MMP_12",
    "MMP7": "MMP_7",
    "VEGFA": "VEGF_A",
    "VEGFD": "VEGF_D",
}
olink_df.rename(columns=correction_table, inplace=True)

# gPCs for the plots
gpcs_df = pd.read_table(r"C:\Study\Europe\Torben\small\2020 08 CVD Target scores\olink\TARGET_mergedPRS_LDpred_v2020.11.20.txt",
                        sep="\t", usecols=["IID", "pc1", "pc2", "pc3", "pc4"]).set_index("IID")

# gender info - TODO better
gender_df = pd.read_excel(r"C:\Study\Europe\Torben\small\2020 08 CVD Target scores\Phenotypes.xlsx",
                          usecols=["blood_sample_ID","gender"]).dropna().set_index("blood_sample_ID")
gender_df.index = gender_df.index.map(lambda x: f"66-{x}")

_full_df1 = pd.merge(total_df, olink_df, how="left", left_index=True, right_index=True)
_full_df2 = pd.merge(_full_df1, gender_df, how="left", left_index=True, right_index=True)
full_df = pd.merge(_full_df2, gpcs_df, how="left", left_index=True, right_index=True)

# QC naming
missing_prs = []
for prs_name in translate.values():
    prot_name = prs_name[:-4]
    if prot_name not in full_df:
        missing_prs.append(prot_name)
assert missing_prs == ['TNF_R2'], "Oops! Some of the olink PRS are not available in data - inconsistent naming?"

proteins = list(sorted(set(x[:-4] for x in translate.values()).difference(missing_prs)))
bonf_adj = 1. / len(proteins)

joint_patch_coords = {
    'CXCL1': (0.2, 1.15),
    'PTX3': (0.25, 2.25),
    'SCF': (0.28, 3.5),
    'CCL3': (0.31, 4.5),
    'TM': (0.31, 6),
}
pop_patch_coords = {
    'MMP_12': (0.2, 6.9),
    'CCL3': (0.28, 3.2),
    'SCF': (0.15, 2.8),
}


for run, _df in (
        ("Joint", full_df),
        ("Population", full_df[full_df.population]),
        ("Obese", full_df[full_df.obese]),
):
    full_df_std = (_df - _df.mean()) / _df.std()

    results = pd.DataFrame(index=proteins, columns=["pvalue", "-log10(pvalue)", "Beta", "R2_adj", "Significant"])
    for p in proteins:
        model = smf.ols(formula=f"{p} ~ 1 + {p}_PRS + "
                                f"C(gender) + age + "
                                f"pc1 + pc2 + pc3 + pc4",
                        data=full_df_std,
                        missing='drop').fit()
        results.at[p, "pvalue"] = model.pvalues[f"{p}_PRS"]
        results.at[p, "-log10(pvalue)"] = -np.log10(model.pvalues[f"{p}_PRS"])
        results.at[p, "Beta"] = model.params[f"{p}_PRS"]
        results.at[p, "R2_adj"] = model.rsquared_adj
        results.at[p, "Significant"] = True if model.pvalues[f"{p}_PRS"] < 0.05 * bonf_adj else False

    colorbar_limits = (-results["Beta"].abs().max() * 1.1,
                        results["Beta"].abs().max() * 1.1)
    ax = sns.scatterplot(
        x="Beta", y="-log10(pvalue)",
        hue="Beta", palette="vlag", hue_norm=colorbar_limits,
        edgecolor="black",
        size="Significant", sizes={True: 250, False: 15},
        legend=False,
        data=results
    )
    # Patching: adding labels to significant dots
    significant_data = results[results["Significant"]][["Beta", "-log10(pvalue)"]].T.to_dict()
    print(significant_data)
    for name, d in significant_data.items():
        if run == "Joint":
            y = d['-log10(pvalue)'] - 0.3
            if name in joint_patch_coords:
                x, y = joint_patch_coords[name]
            else:
                x = d['Beta'] - 0.02 * len(name)
        elif run == "Population":
            y = d['-log10(pvalue)'] - 0.2
            if name in pop_patch_coords:
                x, y = pop_patch_coords[name]
            else:
                x = d['Beta'] - 0.02 * len(name)
        elif run == "Obese":
            y = d['-log10(pvalue)'] - 0.2
            if name == "IL_27":
                x, y = 0.43, 6.2
            elif y > 5:
                x = d['Beta'] - 0.025 * len(name)
            else:
                x = d['Beta'] + 0.05
        ax.text(x=x, y=y, s=name).set_backgroundcolor("#ffffff")
    # Patching: adding a synthetic colorbar to match existing hue_norm
    norm = plt.Normalize(*colorbar_limits)
    sm = plt.cm.ScalarMappable(cmap="vlag", norm=norm)
    sm.set_array([])
    ax.figure.colorbar(sm)
    ax.set_title(run)
    plt.show()
    # TODO: do in R - ask Sara or Evelina
    is_fdr_significant, adj_pv = fdrcorrection(results["pvalue"])
