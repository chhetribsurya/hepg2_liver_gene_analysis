import pandas as pd
import numpy as np
import pybedtools
import pickle
import scipy 
import os
import re
import glob
import scipy.stats as stats
import statsmodels.stats.multitest as smm
from rpy2.robjects.packages import importr  #stats = importr('stats')
from rpy2.robjects.vectors import FloatVector
from os.path import join
from os.path import basename
from os.path import splitext


main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_factor_count_with_ideas_bins"
if not os.path.exists(main_dir):
	os.makedirs(main_dir)

piechart_file_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/piecharts_ideas_total"

""" All TF containing dir """
#dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/test_analysis/SL*"
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
all_tf_file_list = glob.glob(dir_path) # represents all tf file list


"""Read multiple excel sheets from the same file"""
excelfile = "/gpfs/gpfs1/home/schhetri/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx"
df_xls = pd.ExcelFile(excelfile).parse("unique_TFs_DBF_CR_annotation")
df_xls["Category_split_join"] = df_xls["Category"].apply(lambda x: "_".join(str(x).split("/"))) 
uniq_tf_list =  df_xls["Target"].tolist()

file_list = []
for each_tf in uniq_tf_list:
	glob_pattern =  join(piechart_file_dir, each_tf + "_intersectbed_counts_with_ideas.txt" )
	glob_pattern = re.sub(r'\[', '[[]', glob_pattern)
	glob_pattern = re.sub(r'(?<!\[)\]', '[]]', glob_pattern)
	tf_ideas_file = glob.glob(glob_pattern)
	file_list.append(tf_ideas_file)

# flatten the list:
final_file_list = sum(file_list, [])

concat_df_dict = {}
tf_name_list =  []
for each_tf_file in final_file_list:
	tf_name = re.compile(r"(.*)_intersectbed_counts_with_ideas.txt").findall(basename(each_tf_file))[0]
	tf_df = pd.read_csv(each_tf_file, sep="\t", header =None)
	tf_df.columns = ["ideas_state", "peak_count"]
	#concat_df_list.append(tf_df)
	tf_name_list.append(tf_name)
	concat_df_dict[tf_name] = tf_df

combined_tf_df = pd.concat(concat_df_dict).reset_index()
combined_tf_df = combined_tf_df.iloc[:, [0,2,3]]
combined_tf_df = combined_tf_df.rename(columns={"level_0":"tf_name"})
combined_tf_df["bool_peak_count"] = np.where(combined_tf_df["peak_count"] > 0, 1, 0)
combined_tf_df.to_csv(join(piechart_file_dir, "final_tf_ideas_piechart_combined.bed"), sep="\t", header=True, index=False)

""" Generate a barplot data """
tf_count_with_ideas_bin = combined_tf_df.groupby(["ideas_state"]).apply(lambda x : x["bool_peak_count"].sum())
tf_count_with_ideas_bin = tf_count_with_ideas_bin.reset_index().rename(columns={0: "tf_counts"})
tf_count_with_ideas_bin.to_csv(join(piechart_file_dir, "final_tf_ideas_piechart_combined_barplot_data.bed"), sep="\t", header=True, index=False)

""" Generate a heatmap data for the TF binding in each bin"""
combined_tf_df_heatmap_data = combined_tf_df.pivot(index="tf_name", columns="ideas_state", values= "bool_peak_count")
combined_tf_df_heatmap_data.to_csv(join(piechart_file_dir, "final_tf_ideas_piechart_combined_heatmap_data.bed"), sep="\t", header=True, index=True)












