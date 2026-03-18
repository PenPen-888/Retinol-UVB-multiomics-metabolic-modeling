import cobra
import riptide
import os
import pandas as pd
import gurobipy
print(gurobipy.gurobi.version())
model = cobra.io.read_sbml_model(r"E:\human_gems\Human-GEM\model\Human-GEM.xml")
print(len(model.reactions))
print(len(model.metabolites))
print(len(model.genes))
# 读取整个HFF.csv表格
HFF_df = pd.read_csv('HFF.csv')

UVBcols = ["HFF_UVB1", "HFF_UVB2", "HFF_UVB3"]
UAcols = ["HFF_UA1", "HFF_UA2", "HFF_UA3"]
HFF_mean = pd.DataFrame({
    "gene_id": HFF_df["gene_id"],

    "UVBmean": HFF_df[UVBcols].mean(axis=1),
    "UAmean": HFF_df[UAcols].mean(axis=1)
})
HFF_mean.to_csv("HFF_mean.tsv", sep='\t', index=False, header=True)

cobra.Configuration().solver = 'gurobi'
output_dir = "riptide_output_mean"
os.makedirs(output_dir, exist_ok=True)


expression_dict = riptide.read_transcription_file("HFF_mean.tsv", header=True, norm=True)


replicates = expression_dict['replicates']
gene_list = [k for k in expression_dict.keys() if k != 'replicates']


all_models = {}

for i in range(replicates):
    sample_expression = {gene: expression_dict[gene][i] for gene in gene_list}
    sample_name = f"sample_{i + 1}"
    print(f"\n>>> Processing: {sample_name}")

    try:
        contextualized_model = riptide.maxfit(model.copy(), transcriptome=sample_expression, silent=False)
        if contextualized_model is None:
            raise ValueError("RIPTiDe return None")

        all_models[sample_name] = contextualized_model
        output_path = os.path.join(output_dir, f"{sample_name}.tsv")
        riptide.save_output(riptide_obj=contextualized_model, path=output_path)
        print(f"✔ Saved in : {output_path}")

    except Exception as e:
        print(f"❌ Error: {sample_name}  {e}")