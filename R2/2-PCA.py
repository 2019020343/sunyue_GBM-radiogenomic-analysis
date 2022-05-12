import pandas as pd
import numpy as np
import os
from pandas.core.frame import DataFrame
from sklearn.preprocessing import scale
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans



os.chdir(u'E:\\Mask_RCNN_master\\pyradiomics\\test')
ex14_GBM = pd.read_table('F-S_matrix.txt', header=0, index_col=0, sep="\t")

os.chdir(u'E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\PCA')
ex14_LGG = pd.read_table('F-S_matrix.txt', header=0, index_col=0, sep="\t")



ex14_GBM.index
# cluster
ex141_GBM=ex14_GBM._values
df_GBM = ex141_GBM.transpose()
df_GBM = DataFrame(df_GBM)
# plt.figure(figsize=(8, 6))
# sns.heatmap(data=df_GBM.corr(), annot=True)
# PCA 通常用中心标准化，也就是都转化成 Z 分数的形式
data_GBM = scale(df_GBM)

ex14_LGG.index
# cluster
ex141_LGG=ex14_LGG._values
df_LGG = ex141_LGG.transpose()
df_LGG = DataFrame(df_LGG)
# plt.figure(figsize=(8, 6))
# sns.heatmap(data=df_LGG.corr(), annot=True)
# PCA 通常用中心标准化，也就是都转化成 Z 分数的形式
data_LGG = scale(df_LGG)

from sklearn.decomposition import PCA

# 说明：
## 1. 第一次的 n_components 参数最好设置得大一些（保留的主成份）
## 2. 观察 explained_variance_ratio_ 取值变化，即每个主成分能够解释原始数据变异的百分比
# pca = PCA(n_components=len(data)) # 直接与变量个数相同的主成分
# pca.fit(data)


# pca.explained_variance_ratio_
# plt.plot(np.cumsum(pca.explained_variance_ratio_), linewidth=3)
# plt.xlabel('成份数')
# plt.ylabel('累积解释方差')
# plt.grid(True)

# 可知两个主成分就已经足够了
pca = PCA(n_components=10) # 直接与变量个数相同的主成分
pca.fit(data_GBM)
pca.explained_variance_ratio_
new_data_GBM = pca.fit_transform(data_GBM) #  fit_transform 表示将生成降维后的数据
# 查看规模差别
print("原始数据集规模:   ", data_GBM.shape)
print("降维后的数据集规模:", new_data_GBM.shape)

## 1. 输出GBM数据的10个主成分特征值 62*10 （样本*主成分）
new_data_GBM = DataFrame(new_data_GBM)
new_data_GBM.insert(0, "samples", ex14_GBM.columns)
col_name = np.append(np.array(['samples']), ['pca_1', 'pca_2', 'pca_3', 'pca_4', 'pca_5', 'pca_6', 'pca_7', 'pca_8', 'pca_9',  'pca_10'])
#col_name = np.append(np.array(['samples']), ['pca_1', 'pca_2', 'pca_3', 'pca_4', 'pca_5', 'pca_6', 'pca_7'])
out_path ="E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA1.xlsx"
new_data_GBM.to_excel(out_path, header=col_name, index=0)
## 2. 输出GBM数据原特征值的系数 10*86 （主成分*原特征）
components_data = DataFrame(pca.components_)
out_path ="E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA1_score.xlsx"
components_data.to_excel(out_path, header=False, index=0)


## 3. 输出LGG数据的标准化的86个特征值 44*86 （样本*原特征）
data = DataFrame(data_LGG)
out_path ="E:\\Mask_RCNN_master\\pyradiomics\\LGGtest\\PCA_bz1.xlsx"
data.to_excel(out_path, header=False, index=0)
## 3.1 输出GBM数据的标准化的86个特征值 62*86 （样本*原特征）
data1 = DataFrame(data_GBM)
out_path ="E:\\Mask_RCNN_master\\pyradiomics\\GBMtest\\PCA_bz.xlsx"
data1.to_excel(out_path, header=False, index=0)

## 4. 输出LGG数据的是个主成分特征值  44*10（样本*主成分）
bz_LGG = DataFrame(data_LGG)
wright = DataFrame(pca.components_)
LGG_result = np.dot(bz_LGG, wright.transpose())
LGG_result1 = DataFrame(LGG_result)
LGG_result1.insert(0, "samples", ex14_LGG.columns)
col_name = np.append(np.array(['samples']), ['pca_1', 'pca_2', 'pca_3', 'pca_4', 'pca_5', 'pca_6', 'pca_7', 'pca_8', 'pca_9',  'pca_10'])
out_path ="E:\\Mask_RCNN_master\\pyradiomics\\test\\PCA.xlsx"
LGG_result1.to_excel(out_path, header=col_name, index=0)










