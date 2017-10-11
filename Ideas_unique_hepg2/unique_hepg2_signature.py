### Merging of the files:


In [12]: pwd
Out[12]: u'/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks/test_analysis'


In [13]: file1='SL120595.filt.nodup.srt.SE_VS_SL146889.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak_FOXO1'

In [14]: file3="SL88833.filt.nodup.srt.SE_VS_SL88834.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak_RAD21[FLAG]"

In [15]: file2="SL151607.filt.nodup.srt.SE_VS_SL151608.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak_FOXA3[FLAG]"

In [16]: files=[file1, file2, file3]

In [17]: pybed = pybedtools.BedTool(files[0])

In [18]: concat_pybed = pybed.cat(*files[1:])

In [20]: concat_pybed.count()
Out[20]: 75156


### Other alternative:
##########################
In [22]: file_list = glob.glob("./SL*")

In [38]: concat_list = []

In [39]: for each in file_list:
    read_file = pd.read_csv(each, sep="\t", header=None)
    concat_list.append(read_file)
   ....:

In [40]: concat_df = pd.concat(concat_list, ignore_index=True)
n [50]: final_df = concat_df.loc[:,[0,1,2]]

In [52]: final_df.columns = ["chrom", "start", "end"]

In [53]: final_df_sorted = final_df.sort_values(["chrom", "start"])


In [60]: pybed_df = pybedtools.BedTool.from_dataframe(final_df_sorted)

In [61]: pybed_merge = pybed_df.merge()

In [64]: pybed_merge.count()
Out[64]: 75156

In [65]: concat_pybed == pybed_merge
Out[65]: True
