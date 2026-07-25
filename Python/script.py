import pandas as pd
import numpy as np

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

    concentration_columns = [
        col
        for col in df.columns
        if any(c in str(col) for c in ["20", "40", "80", "160", "320", "640"])
    ]

    # 选项B: 如果列名是Sample_1, Sample_2等，按顺序对应浓度
    # concentration_columns = [col for col in df.columns if 'Sample' in str(col)]

    # 选项C: 直接指定列名（最保险）
    # concentration_columns = ['20uM', '40uM', '80uM', '160uM', '320uM', '640uM']

    if len(concentration_columns) != 6:
        print(f"\n警告: 检测到{len(concentration_columns)}个浓度列，期望6个")
        print("检测到的列:", concentration_columns)
        print("\n请手动指定浓度列名，修改代码中的concentration_columns变量")
        return None

    # 获取样本名称列（假设是第一列）
    name_column = df.columns[0]

    # 创建输出数据框
    output_df = pd.DataFrame()
    output_df["Name"] = df[name_column]

    # 重命名浓度列为Raw_1到Raw_6
    for i, col in enumerate(concentration_columns, start=1):
        output_df[f"Raw_{i}"] = df[col]

    # ========================================================================
    # 数据质量检查
    # ========================================================================
    print("\n" + "=" * 60)
    print("数据质量检查")
    print("=" * 60)

    missing_counts = output_df.isnull().sum()
    print("\n各列缺失值数量:")
    print(missing_counts)

    row_missing = output_df.iloc[:, 1:].isnull().sum(axis=1)
    print(f"\n缺失值统计:")
    print(f"  - 无缺失值的样本: {(row_missing == 0).sum()}")
    print(f"  - 有1个缺失值: {(row_missing == 1).sum()}")
    print(f"  - 有2个缺失值: {(row_missing == 2).sum()}")
    print(f"  - 有>2个缺失值: {(row_missing > 2).sum()}")

    if (row_missing > 2).sum() > 0:
        print("\n警告: 有样本缺失值>2个，可能会被CurveCurator过滤")
        print("缺失值>2的样本:")
        print(output_df[row_missing > 2]["Name"].tolist())

    numeric_cols = [f"Raw_{i}" for i in range(1, 7)]
    negative_mask = (output_df[numeric_cols] < 0).any(axis=1)
    if negative_mask.sum() > 0:
        print(f"\n警告: 发现{negative_mask.sum()}个样本含有负值")
        print("含负值的样本:")
        print(output_df[negative_mask]["Name"].tolist())

    print("\n数据范围统计:")
    print(output_df[numeric_cols].describe())

    output_df.to_csv(output_file, sep="\t", index=False)
    print(f"\n转换完成!")
    print(f"输出文件: {output_file}")
    print(f"样本数量: {len(output_df)}")
    print(f"\n输出数据预览:")
    print(output_df.head(10))

    csv_output = output_file.replace(".txt", "_check.csv")
    output_df.to_csv(csv_output, index=False)
    print(f"\n同时保存CSV格式用于检查: {csv_output}")

    return output_df


# ============================================================================
# 方案2: 长格式转换
# ============================================================================


def convert_long_format(input_file, output_file):
    """
    转换长格式数据
    假设数据格式为: | Sample_ID | Concentration | Value |
    """
    df = pd.read_csv(input_file)

    pivot_df = df.pivot(index="Sample_ID", columns="Concentration", values="Value")

    concentrations = [20, 40, 80, 160, 320, 640]
    pivot_df = pivot_df[concentrations]

    output_df = pd.DataFrame()
    output_df["Name"] = pivot_df.index
    for i in range(1, 7):
        output_df[f"Raw_{i}"] = pivot_df.iloc[:, i - 1].values

    output_df.to_csv(output_file, sep="\t", index=False)
    print(f"转换完成: {output_file}")

    return output_df


# ============================================================================
# 步骤2: site-level 数据处理
# ============================================================================


def process_site_level_data(input_file, output_file):
    """
    对 site-level 数据进行聚合、透视、导出 CurveCurator 所需格式
    """
    print("正在读取数据...")
    df = pd.read_csv(input_file)
    print(f"读取完成: {df.shape[0]} 行")

    print("\n开始数据转换...")

    # 1. 对每个site和浓度，计算3个重复的平均值
    print("  - 计算每个浓度的平均值（合并3个重复）...")
    df_avg = (
        df.groupby(
            [
                "site_id",
                "concentration",
                "Protein.Id",
                "gene_symbol",
                "prot_description",
                "Site.Position",
                "Motif",
                "Final_Protein_Class",
            ]
        )
        .agg(
            {
                "intensity": "mean",
                "log2_intensity": "mean",
            }
        )
        .reset_index()
    )

    print(f"  - 合并后: {df_avg.shape[0]} 行")

    # 2. 透视表：将浓度从行转为列
    print("  - 转换为宽格式...")
    pivot_df = df_avg.pivot(
        index=[
            "site_id",
            "Protein.Id",
            "gene_symbol",
            "prot_description",
            "Site.Position",
            "Motif",
            "Final_Protein_Class",
        ],
        columns="concentration",
        values="intensity",
    ).reset_index()

    # 3. 重命名列
    print("  - 重命名列...")
    concentrations = [20, 40, 80, 160, 320, 640]
    column_mapping = {conc: f"Raw_{i}" for i, conc in enumerate(concentrations, 1)}
    pivot_df.rename(columns=column_mapping, inplace=True)

    # 4. 创建Name列
    print("  - 创建样本标识符...")
    pivot_df["Name"] = (
        pivot_df["gene_symbol"].astype(str)
        + "_"
        + pivot_df["site_id"].astype(str)
        + "_"
        + pivot_df["Site.Position"].astype(str)
    )

    # 5. 选择输出列
    output_columns = ["Name"] + [f"Raw_{i}" for i in range(1, 7)]
    output_df = pivot_df[output_columns]

    # ========================================================================
    # 数据质量检查
    # ========================================================================
    print("\n" + "=" * 80)
    print("数据质量检查")
    print("=" * 80)

    missing_counts = output_df.isnull().sum()
    print("\n各列缺失值数量:")
    for col in output_df.columns:
        if missing_counts[col] > 0:
            print(f"  {col}: {missing_counts[col]}")

    if missing_counts.sum() == 0:
        print("  ✓ 无缺失值")

    row_missing = output_df.iloc[:, 1:].isnull().sum(axis=1)
    print(f"\n缺失值统计:")
    print(f"  - 无缺失值的样本: {(row_missing == 0).sum()}")
    print(f"  - 有1个缺失值: {(row_missing == 1).sum()}")
    print(f"  - 有2个缺失值: {(row_missing == 2).sum()}")
    print(f"  - 有>2个缺失值: {(row_missing > 2).sum()}")

    if (row_missing > 2).sum() > 0:
        print("\n  ⚠ 警告: 有样本缺失值>2个，可能会被CurveCurator过滤")

    numeric_cols = [f"Raw_{i}" for i in range(1, 7)]
    negative_mask = (output_df[numeric_cols] < 0).any(axis=1)
    zero_mask = (output_df[numeric_cols] == 0).any(axis=1)

    if negative_mask.sum() > 0:
        print(f"\n  ⚠ 警告: 发现 {negative_mask.sum()} 个样本含有负值")
    if zero_mask.sum() > 0:
        print(f"\n  ⚠ 警告: 发现 {zero_mask.sum()} 个样本含有零值")

    print("\n数据范围统计:")
    print(output_df[numeric_cols].describe())

    # ========================================================================
    # 保存输出文件
    # ========================================================================
    print("\n" + "=" * 80)
    print("保存文件")
    print("=" * 80)

    output_df.to_csv(output_file, sep="\t", index=False)
    print(f"✓ 主输出文件: {output_file}")
    print(f"  格式: Tab分隔")
    print(f"  样本数量: {len(output_df)}")

    csv_output = output_file.replace(".txt", "_check.csv")
    output_df.to_csv(csv_output, index=False)
    print(f"✓ 检查文件: {csv_output}")
    print(f"  格式: CSV")

    full_output = output_file.replace(".txt", "_full.csv")
    full_df = pivot_df[
        [
            "Name",
            "site_id",
            "Protein.Id",
            "gene_symbol",
            "prot_description",
            "Site.Position",
            "Motif",
            "Final_Protein_Class",
        ]
        + [f"Raw_{i}" for i in range(1, 7)]
    ]
    full_df.to_csv(full_output, index=False)
    print(f"✓ 完整信息文件: {full_output}")
    print(f"  包含: Protein ID, Gene symbol, 描述等")

    print("\n" + "=" * 80)
    print("输出数据预览（前10行）")
    print("=" * 80)
    print(output_df.head(10).to_string())

    print("\n" + "=" * 80)
    print("转换完成！")
    print("=" * 80)
    print("\n下一步:")
    print("1. 检查文件: " + csv_output)
    print("2. 确认数据格式正确")
    print("3. 准备运行CurveCurator")
    print("\n生成的文件:")
    print(f"  - {output_file} (CurveCurator输入)")
    print(f"  - {csv_output} (Excel检查)")
    print(f"  - {full_output} (完整信息)")

    return output_df, pivot_df, csv_output, full_output


# ============================================================================
# 步骤3: 生成 CurveCurator 配置文件
# ============================================================================


def generate_toml_config(toml_file):
    toml_content = """# ============================================================================
# CurveCurator Configuration File
# TMT Direct Labeling Titration Experiment
# ============================================================================

# ----------------------------------------------------------------------------
# [Meta] - 实验元数据
# ----------------------------------------------------------------------------
['Meta']
id = "TMT4_direct_labeling_wxr35"
description = "Direct labeling desthiobiotin pipeline - Site level analysis"
condition = "wxr35_treatment"
warhead = "wxr35"
data_level = "site"

# ----------------------------------------------------------------------------
# [Experiment] - 实验设计
# ----------------------------------------------------------------------------
['Experiment']
experiments = ["1", "2", "3", "4", "5", "6"]

# Raw_1=20, Raw_2=40, Raw_3=80, Raw_4=160, Raw_5=320, Raw_6=640
doses = [20.0, 40.0, 80.0, 160.0, 320.0, 640.0]

dose_scale = "1e-6"
dose_unit = "M"

control_experiment = "1"

measurement_type = "TMT"
data_type = "SITE"
search_engine = "OTHER"
search_engine_version = "custom_pipeline"

# ----------------------------------------------------------------------------
# [Paths] - 文件路径
# ----------------------------------------------------------------------------
['Paths']
input_file = "./curvecurator_input.txt"
curves_file = "./output_curves.tsv"
normalization_file = "./normalization_factors.tsv"
mad_file = "./mad_analysis.tsv"
dashboard = "./interactive_dashboard.html"

# decoys_file = "./output_decoys.tsv"
# fdr_file = "./output_fdr.tsv"

# ----------------------------------------------------------------------------
# [Processing] - 数据预处理
# ----------------------------------------------------------------------------
['Processing']
available_cores = 4

imputation = true
imputation_pct = 0.005

max_missing = 2
max_imputation = 2

normalization = true

# ratio_range = [0.0, 10.0]

# ----------------------------------------------------------------------------
# [Curve Fit] - 曲线拟合参数
# ----------------------------------------------------------------------------
['Curve Fit']
type = "OLS"
speed = "standard"
max_iterations = 1000
control_fold_change = true
interpolation = false

# slope = 1.0
# front = 1.0
# back = 0.0

# ----------------------------------------------------------------------------
# [F Statistic] - 统计检验
# ----------------------------------------------------------------------------
['F Statistic']
alpha = 0.05
fc_lim = 1.0
optimized_dofs = true
mtc_method = "fdr_bh"
not_rmse_limit = 0.1

# quality_min = 0.0

# ----------------------------------------------------------------------------
# [Dashboard] - 可视化
# ----------------------------------------------------------------------------
['Dashboard']
backend = "webgl"

# plot_width = 800
# plot_height = 600
"""

    with open(toml_file, "w", encoding="utf-8") as f:
        f.write(toml_content)

    print("=" * 80)
    print("TOML配置文件已生成")
    print("=" * 80)
    print(f"文件位置: {toml_file}")
    print("\n配置摘要:")
    print("  - 浓度梯度: 20, 40, 80, 160, 320, 640 μM")
    print("  - Control: 20 μM (Raw_1)")
    print("  - 数据类型: TMT site-level")
    print("  - 标准化: 开启")
    print("  - 缺失值填补: 开启")
    print("  - FDR校正: Benjamini-Hochberg")
    print("\n" + "=" * 80)
    print("下一步: 运行CurveCurator")
    print("=" * 80)
    print("\n命令:")
    print("cd E:\\python\\curvecurator")
    print("curvecurator --config curvecurator_config.toml")
    print("\n或者如果需要详细日志:")
    print("curvecurator --config curvecurator_config.toml --verbose")


# ============================================================================
# 主程序入口
# ============================================================================


import sys


def main(input_file, output_file, toml_file):
    try:
        process_site_level_data(input_file, output_file)
        generate_toml_config(toml_file)

    except Exception as e:
        print(f"\n错误: {e}")
        print("\n请检查:")
        print("1. 输入文件路径是否正确")
        print("2. 文件格式是否正确")
        print("3. 列名是否正确")
        print("\n如需帮助，请提供:")
        print("- 原始数据的前几行")
        print("- 完整的列名列表")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("用法:")
        print("python script.py <input_file> <output_file> <toml_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    toml_file = sys.argv[3]

    main(input_file, output_file, toml_file)