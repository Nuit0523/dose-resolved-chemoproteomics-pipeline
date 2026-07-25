# ============================================================================
# CurveCurator Configuration File
# Direct Labeling Titration Experiment
# ============================================================================

# ----------------------------------------------------------------------------
# [Meta] - 样本特定信息
# ----------------------------------------------------------------------------
[Meta]
id = "titration_experiment_001"
description = "Direct labeling titration experiment with concentration gradient 20-640 μM"
condition = "drug_treatment"
treatment_time = 24  # 处理时间（小时），根据实际情况修改

# ----------------------------------------------------------------------------
# [Experiment] - 实验设计信息
# ----------------------------------------------------------------------------
[Experiment]
# 实验ID列表（对应你的6个浓度点）
experiments = ["1", "2", "3", "4", "5", "6"]

# 药物浓度（与experiments顺序对应）
# 将20 μM设为参考点(0.0)，其他浓度相对于20 μM
doses = [0.0, 40.0, 80.0, 160.0, 320.0, 640.0]

# 浓度单位换算（微摩尔）
dose_scale = "1e-6"

# 浓度基础单位
dose_unit = "M"

# Control实验（使用最低浓度20 μM作为参考）
control_experiment = "1"

# 如果是蛋白质组学数据，取消下面注释并填写
# measurement_type = "TMT"  # 或 "LFQ", "DIA", "OTHER"
# data_type = "PROTEIN"     # 或 "PEPTIDE", "OTHER"
# search_engine = "MAXQUANT"  # 或 "DIANN", "PD", "MSFRAGGER", "OTHER"
# search_engine_version = "2.0.0"

# ----------------------------------------------------------------------------
# [Paths] - 输入输出文件路径（相对于TOML文件位置）
# ----------------------------------------------------------------------------
[Paths]
# 输入数据文件（必需）
input_file = "./your_data.csv"  # 修改为你的实际文件名

# 输出文件（可选，不指定则使用默认名称）
curves_file = "./output_curves.tsv"
# decoys_file = "./output_decoys.tsv"  # FDR模式时使用
# fdr_file = "./output_fdr.tsv"        # FDR模式时使用
# normalization_file = "./normalization_factors.tsv"
# mad_file = "./mad_analysis.tsv"
# dashboard = "./interactive_dashboard.html"

# ----------------------------------------------------------------------------
# [Processing] - 数据预处理参数
# ----------------------------------------------------------------------------
[Processing]
# 并行计算核心数
available_cores = 4  # 根据你的计算机配置调整

# 缺失值填补
imputation = true
imputation_pct = 0.005  # 使用0.5%分位数填补

# 缺失值容忍度
max_missing = 2  # 每条曲线最多允许2个缺失值（不包括control）
max_imputation = 2  # 最多允许2个填补值

# 全局标准化（推荐用于蛋白质组学数据）
normalization = true

# 比值范围限制（可选）
# ratio_range = [0.0, 10.0]  # 限制比值在0-10之间

# ----------------------------------------------------------------------------
# [Curve Fit] - 曲线拟合参数
# ----------------------------------------------------------------------------
[Curve Fit]
# 拟合类型
type = "OLS"  # "OLS" 或 "MLE"

# 拟合速度
speed = "standard"  # "fast", "standard", "exhaustive", "basinhopping"

# 最大迭代次数
max_iterations = 1000

# 权重（可选，默认所有点权重相同）
# weights = [1, 1, 1, 1, 1, 1]

# 固定参数（可选，根据需要取消注释）
# slope = 1.0
# front = 1.0
# back = 0.0

# Fold change计算相对于control
control_fold_change = true

# 插值（可选，增加拟合鲁棒性但会略微降低p值）
interpolation = false

# ----------------------------------------------------------------------------
# [F Statistic] - 统计检验和显著性阈值
# ----------------------------------------------------------------------------
[F Statistic]
# 显著性阈值
alpha = 0.05

# Fold change阈值（log2）
fc_lim = 1.0  # 相当于2倍变化

# F分布参数（使用优化的自由度）
optimized_dofs = true

# 质量过滤（可选）
# quality_min = 0.0

# 多重检验校正方法
mtc_method = "fdr_bh"  # Benjamini-Hochberg FDR校正
# 其他选项: "bonferroni", "fdr_by", "holm", 等

# NOT分类阈值
not_rmse_limit = 0.1
# not_p_limit = 0.05  # 可选的额外p值过滤

# FDR估计（如果使用--fdr模式）
# decoy_ratio = 1.0

# ----------------------------------------------------------------------------
# [Dashboard] - 可视化参数（可选）
# ----------------------------------------------------------------------------
[Dashboard]
# 后端渲染引擎
backend = "webgl"  # "webgl", "svg", "canvas"




import pandas as pd

# ============================================================================
# 数据格式转换脚本 - 适配CurveCurator
# ============================================================================

def convert_to_curvecurator_format(input_file, output_file):
    """
    将原始数据转换为CurveCurator要求的格式
    
    Parameters:
    -----------
    input_file : str
        输入CSV文件路径
    output_file : str
        输出TXT文件路径（tab分隔）
    """
    
    # 读取原始数据
    df = pd.read_csv(input_file)
    
    print("原始数据形状:", df.shape)
    print("原始列名:", df.columns.tolist())
    print("\n前5行数据:")
    print(df.head())
    
    # ========================================================================
    # 方案1: 如果你的数据是这样的结构：
    # | Sample_ID | 20uM | 40uM | 80uM | 160uM | 320uM | 640uM |
    # ========================================================================
    
    # 假设第一列是样本名称，后面6列是浓度数据
    # 根据你的实际列名修改下面的列名列表
    
    # 选项A: 如果列名包含浓度信息（如 "20uM", "40uM"等）
    concentration_columns = [col for col in df.columns if any(c in str(col) for c in ['20', '40', '80', '160', '320', '640'])]
    
    # 选项B: 如果列名是Sample_1, Sample_2等，按顺序对应浓度
    # concentration_columns = [col for col in df.columns if 'Sample' in str(col)]
    
    # 选项C: 直接指定列名（最保险）
    # concentration_columns = ['20uM', '40uM', '80uM', '160uM', '320uM', '640uM']
    
    if len(concentration_columns) != 6:
        print(f"\n警告: 检测到{len(concentration_columns)}个浓度列，期望6个")
        print("检测到的列:", concentration_columns)
        print("\n请手动指定浓度列名，修改代码中的concentration_columns变量")
        return
    
    # 获取样本名称列（假设是第一列）
    name_column = df.columns[0]
    
    # 创建输出数据框
    output_df = pd.DataFrame()
    output_df['Name'] = df[name_column]
    
    # 重命名浓度列为Raw_1到Raw_6
    for i, col in enumerate(concentration_columns, start=1):
        output_df[f'Raw_{i}'] = df[col]
    
    # ========================================================================
    # 数据质量检查
    # ========================================================================
    print("\n" + "="*60)
    print("数据质量检查")
    print("="*60)
    
    # 检查缺失值
    missing_counts = output_df.isnull().sum()
    print("\n各列缺失值数量:")
    print(missing_counts)
    
    # 检查每行缺失值数量
    row_missing = output_df.iloc[:, 1:].isnull().sum(axis=1)
    print(f"\n缺失值统计:")
    print(f"  - 无缺失值的样本: {(row_missing == 0).sum()}")
    print(f"  - 有1个缺失值: {(row_missing == 1).sum()}")
    print(f"  - 有2个缺失值: {(row_missing == 2).sum()}")
    print(f"  - 有>2个缺失值: {(row_missing > 2).sum()}")
    
    if (row_missing > 2).sum() > 0:
        print("\n警告: 有样本缺失值>2个，可能会被CurveCurator过滤")
        print("缺失值>2的样本:")
        print(output_df[row_missing > 2]['Name'].tolist())
    
    # 检查负值
    numeric_cols = [f'Raw_{i}' for i in range(1, 7)]
    negative_mask = (output_df[numeric_cols] < 0).any(axis=1)
    if negative_mask.sum() > 0:
        print(f"\n警告: 发现{negative_mask.sum()}个样本含有负值")
        print("含负值的样本:")
        print(output_df[negative_mask]['Name'].tolist())
    
    # 数据范围统计
    print("\n数据范围统计:")
    print(output_df[numeric_cols].describe())
    
    # ========================================================================
    # 保存输出文件
    # ========================================================================
    # 保存为tab分隔的txt文件
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n转换完成!")
    print(f"输出文件: {output_file}")
    print(f"样本数量: {len(output_df)}")
    print(f"\n输出数据预览:")
    print(output_df.head(10))
    
    # 同时保存一份CSV格式方便检查
    csv_output = output_file.replace('.txt', '_check.csv')
    output_df.to_csv(csv_output, index=False)
    print(f"\n同时保存CSV格式用于检查: {csv_output}")
    
    return output_df


# ============================================================================
# 方案2: 如果你的数据是长格式（每行一个样本-浓度组合）
# ============================================================================

def convert_long_format(input_file, output_file):
    """
    转换长格式数据
    假设数据格式为: | Sample_ID | Concentration | Value |
    """
    df = pd.read_csv(input_file)
    
    # 透视表转换
    pivot_df = df.pivot(index='Sample_ID', columns='Concentration', values='Value')
    
    # 确保浓度按顺序排列
    concentrations = [20, 40, 80, 160, 320, 640]
    pivot_df = pivot_df[concentrations]
    
    # 重命名列
    output_df = pd.DataFrame()
    output_df['Name'] = pivot_df.index
    for i in range(1, 7):
        output_df[f'Raw_{i}'] = pivot_df.iloc[:, i-1].values
    
    # 保存
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"转换完成: {output_file}")
    
    return output_df


# ============================================================================
# 主程序
# ============================================================================

if __name__ == "__main__":
    
    # 设置输入输出文件路径
    input_file = "your_data.csv"  # 修改为你的实际文件名
    output_file = "curvecurator_input.txt"
    
    # 执行转换
    try:
        result_df = convert_to_curvecurator_format(input_file, output_file)
        
        print("\n" + "="*60)
        print("转换成功！")
        print("="*60)
        print("\n下一步:")
        print("1. 检查输出文件格式是否正确")
        print("2. 使用生成的TOML配置文件")
        print("3. 运行CurveCurator:")
        print("   curvecurator --config your_config.toml")
        
    except Exception as e:
        print(f"\n错误: {e}")
        print("\n请检查:")
        print("1. 输入文件路径是否正确")
        print("2. 文件格式是否为CSV")
        print("3. 列名是否正确")
        print("\n如需帮助，请提供:")
        print("- 原始数据的前几行")
        print("- 完整的列名列表")

import pandas as pd
import numpy as np

# ============================================================================
# 步骤2: 转换数据为CurveCurator格式
# ============================================================================

# 文件路径
input_file = r"E:\python\curvecurator\curvecurator_input.txt"
output_file = r"E:\python\curvecurator\curvecurator_input.txt"

# 读取数据
print("正在读取数据...")
df = pd.read_csv(input_file)
print(f"读取完成: {df.shape[0]} 行")

# ============================================================================
# 数据处理
# ============================================================================

print("\n开始数据转换...")

# 1. 对每个site和浓度，计算3个重复的平均值
print("  - 计算每个浓度的平均值（合并3个重复）...")
df_avg = df.groupby(['site_id', 'concentration', 'Protein.Id', 'gene_symbol', 
                      'prot_description', 'Site.Position', 'Motif', 
                      'Final_Protein_Class']).agg({
    'intensity': 'mean',  # 取平均值
    'log2_intensity': 'mean'
}).reset_index()

print(f"  - 合并后: {df_avg.shape[0]} 行")

# 2. 透视表：将浓度从行转为列
print("  - 转换为宽格式...")
pivot_df = df_avg.pivot(
    index=['site_id', 'Protein.Id', 'gene_symbol', 'prot_description', 
           'Site.Position', 'Motif', 'Final_Protein_Class'],
    columns='concentration',
    values='intensity'
).reset_index()

# 3. 重命名列
print("  - 重命名列...")
# 浓度列按顺序对应 Raw_1 到 Raw_6
concentrations = [20, 40, 80, 160, 320, 640]
column_mapping = {conc: f'Raw_{i}' for i, conc in enumerate(concentrations, 1)}

# 重命名浓度列
pivot_df.rename(columns=column_mapping, inplace=True)

# 4. 创建Name列（唯一标识符）
print("  - 创建样本标识符...")
# 格式: gene_symbol_siteID_position
pivot_df['Name'] = (pivot_df['gene_symbol'].astype(str) + '_' + 
                     pivot_df['site_id'].astype(str) + '_' + 
                     pivot_df['Site.Position'].astype(str))

# 5. 选择并排序列
output_columns = ['Name'] + [f'Raw_{i}' for i in range(1, 7)]
output_df = pivot_df[output_columns]

# ============================================================================
# 数据质量检查
# ============================================================================

print("\n" + "="*80)
print("数据质量检查")
print("="*80)

# 检查缺失值
missing_counts = output_df.isnull().sum()
print("\n各列缺失值数量:")
for col in output_df.columns:
    if missing_counts[col] > 0:
        print(f"  {col}: {missing_counts[col]}")

if missing_counts.sum() == 0:
    print("  ✓ 无缺失值")

# 检查每行缺失值数量
row_missing = output_df.iloc[:, 1:].isnull().sum(axis=1)
print(f"\n缺失值统计:")
print(f"  - 无缺失值的样本: {(row_missing == 0).sum()}")
print(f"  - 有1个缺失值: {(row_missing == 1).sum()}")
print(f"  - 有2个缺失值: {(row_missing == 2).sum()}")
print(f"  - 有>2个缺失值: {(row_missing > 2).sum()}")

if (row_missing > 2).sum() > 0:
    print("\n  ⚠ 警告: 有样本缺失值>2个，可能会被CurveCurator过滤")

# 检查负值和零值
numeric_cols = [f'Raw_{i}' for i in range(1, 7)]
negative_mask = (output_df[numeric_cols] < 0).any(axis=1)
zero_mask = (output_df[numeric_cols] == 0).any(axis=1)

if negative_mask.sum() > 0:
    print(f"\n  ⚠ 警告: 发现 {negative_mask.sum()} 个样本含有负值")
if zero_mask.sum() > 0:
    print(f"\n  ⚠ 警告: 发现 {zero_mask.sum()} 个样本含有零值")

# 数据范围统计
print("\n数据范围统计:")
print(output_df[numeric_cols].describe())

# ============================================================================
# 保存输出文件
# ============================================================================

print("\n" + "="*80)
print("保存文件")
print("="*80)

# 保存为tab分隔的txt文件（CurveCurator要求）
output_df.to_csv(output_file, sep='\t', index=False)
print(f"✓ 主输出文件: {output_file}")
print(f"  格式: Tab分隔")
print(f"  样本数量: {len(output_df)}")

# 同时保存CSV格式方便检查
csv_output = output_file.replace('.txt', '_check.csv')
output_df.to_csv(csv_output, index=False)
print(f"✓ 检查文件: {csv_output}")
print(f"  格式: CSV")

# 保存一个包含额外信息的完整文件
full_output = output_file.replace('.txt', '_full.csv')
full_df = pivot_df[['Name', 'site_id', 'Protein.Id', 'gene_symbol', 
                     'prot_description', 'Site.Position', 'Motif', 
                     'Final_Protein_Class'] + [f'Raw_{i}' for i in range(1, 7)]]
full_df.to_csv(full_output, index=False)
print(f"✓ 完整信息文件: {full_output}")
print(f"  包含: Protein ID, Gene symbol, 描述等")

# ============================================================================
# 显示结果预览
# ============================================================================

print("\n" + "="*80)
print("输出数据预览（前10行）")
print("="*80)
print(output_df.head(10).to_string())

print("\n" + "="*80)
print("转换完成！")
print("="*80)
print("\n下一步:")
print("1. 检查文件: " + csv_output)
print("2. 确认数据格式正确")
print("3. 准备运行CurveCurator")
print("\n生成的文件:")
print(f"  - {output_file} (CurveCurator输入)")
print(f"  - {csv_output} (Excel检查)")
print(f"  - {full_output} (完整信息)")

# ============================================================================
# 步骤3: 生成CurveCurator配置文件
# ============================================================================

toml_content = """# ============================================================================
# CurveCurator Configuration File
# TMT Direct Labeling Titration Experiment
# ============================================================================

# ----------------------------------------------------------------------------
# [Meta] - 实验元数据
# ----------------------------------------------------------------------------
[Meta]
id = "TMT4_direct_labeling_wxr35"
description = "Direct labeling desthiobiotin pipeline - Site level analysis"
condition = "wxr35_treatment"
warhead = "wxr35"
data_level = "site"

# ----------------------------------------------------------------------------
# [Experiment] - 实验设计
# ----------------------------------------------------------------------------
[Experiment]
# 实验ID（对应6个浓度）
experiments = ["1", "2", "3", "4", "5", "6"]

# 药物浓度（微摩尔）
# Raw_1=20, Raw_2=40, Raw_3=80, Raw_4=160, Raw_5=320, Raw_6=640
doses = [20.0, 40.0, 80.0, 160.0, 320.0, 640.0]

# 浓度单位
dose_scale = "1e-6"  # 微摩尔转为摩尔
dose_unit = "M"

# Control实验（使用最低浓度20 μM作为参考）
control_experiment = "1"

# TMT蛋白质组学数据
measurement_type = "TMT"
data_type = "SITE"  # 磷酸化位点水平
search_engine = "OTHER"
search_engine_version = "custom_pipeline"

# ----------------------------------------------------------------------------
# [Paths] - 文件路径
# ----------------------------------------------------------------------------
[Paths]
# 输入数据文件
input_file = "./curvecurator_input.txt"

# 输出文件
curves_file = "./output_curves.tsv"
normalization_file = "./normalization_factors.tsv"
mad_file = "./mad_analysis.tsv"
dashboard = "./interactive_dashboard.html"

# 如果使用FDR模式，取消下面注释
# decoys_file = "./output_decoys.tsv"
# fdr_file = "./output_fdr.tsv"

# ----------------------------------------------------------------------------
# [Processing] - 数据预处理
# ----------------------------------------------------------------------------
[Processing]
# 并行计算
available_cores = 4  # 根据你的CPU调整

# 缺失值处理
imputation = true
imputation_pct = 0.005  # 使用0.5%分位数填补

# 缺失值容忍度
max_missing = 2  # 每条曲线最多允许2个缺失值（不包括control）
max_imputation = 2  # 最多允许2个填补值

# 全局标准化（推荐用于TMT数据）
normalization = true

# 比值范围限制（可选）
# ratio_range = [0.0, 10.0]

# ----------------------------------------------------------------------------
# [Curve Fit] - 曲线拟合参数
# ----------------------------------------------------------------------------
[Curve Fit]
# 拟合类型
type = "OLS"  # Ordinary Least Squares

# 拟合速度
speed = "standard"  # "fast", "standard", "exhaustive"

# 最大迭代次数
max_iterations = 1000

# Fold change相对于control
control_fold_change = true

# 插值（可选，增加鲁棒性）
interpolation = false

# 固定参数（可选）
# slope = 1.0
# front = 1.0
# back = 0.0

# ----------------------------------------------------------------------------
# [F Statistic] - 统计检验
# ----------------------------------------------------------------------------
[F Statistic]
# 显著性阈值
alpha = 0.05

# Fold change阈值（log2）
fc_lim = 1.0  # 相当于2倍变化

# F分布自由度优化
optimized_dofs = true

# 多重检验校正
mtc_method = "fdr_bh"  # Benjamini-Hochberg FDR

# NOT分类阈值
not_rmse_limit = 0.1

# 质量过滤（可选）
# quality_min = 0.0

# ----------------------------------------------------------------------------
# [Dashboard] - 可视化
# ----------------------------------------------------------------------------
[Dashboard]
# 渲染引擎
backend = "webgl"  # "webgl", "svg", "canvas"

# 其他可视化选项
# plot_width = 800
# plot_height = 600
"""

# 保存TOML文件
toml_file = r"E:\python\curvecurator\curvecurator_config.toml"
with open(toml_file, 'w', encoding='utf-8') as f:
    f.write(toml_content)

print("="*80)
print("TOML配置文件已生成")
print("="*80)
print(f"文件位置: {toml_file}")
print("\n配置摘要:")
print("  - 浓度梯度: 20, 40, 80, 160, 320, 640 μM")
print("  - Control: 20 μM (Raw_1)")
print("  - 数据类型: TMT site-level")
print("  - 标准化: 开启")
print("  - 缺失值填补: 开启")
print("  - FDR校正: Benjamini-Hochberg")
print("\n" + "="*80)
print("下一步: 运行CurveCurator")
print("="*80)
print("\n命令:")
print(f"cd E:\\python\\curvecurator")
print(f"curvecurator --config curvecurator_config.toml")
print("\n或者如果需要详细日志:")
print(f"curvecurator --config curvecurator_config.toml --verbose")