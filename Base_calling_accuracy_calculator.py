#!/usr/bin/env python
# coding: utf-8

#################################################################################################
#
#  Python script to calculate the base calling accuracy based on Phred Score for bases or reads count
#
#  Required for fastqc_per_sequence_quality_scores_plot.tsv: 
#  export from data of Per Sequence Quality Scores in MultiQC, 
#  or other matrix (bases count matrix or reads count matrix) with same strcture, 
#  for reads count matrix, all reads should be same length
#
################################################################################################


import numpy as np
import pandas as pd


df=pd.read_csv('./fastqc_per_sequence_quality_scores_plot.tsv', header=0, sep='\t')
df.fillna(0, inplace=True)
colraw=df.columns.tolist()
colraw.remove('Mean Sequence Quality (Phred Score)')
dfTraw=df[colraw].copy()
dfTraw.loc['Base calling accuracy (%)']=dfTraw.apply(lambda x: x.sum(), axis=0)
dfTraw['Total']= dfTraw.apply(lambda x: x.sum(), axis=1)
dfTraw=dfTraw[-1:]
dfEraw= df.apply(lambda x: x[colraw]*(10**(x['Mean Sequence Quality (Phred Score)']*(-0.1))), axis=1)
dfEraw.loc['Base calling accuracy (%)']=dfEraw.apply(lambda x: x.sum(), axis=0)
dfEraw['Total']= dfEraw.apply(lambda x: x.sum(), axis=1)
dfEraw=dfEraw[-1:]
dfraw=100*(dfTraw-dfEraw)/dfTraw
dfraw.to_csv('./Base calling accuracy_raw.txt', header=True, index=True, sep='\t', na_rep='')
