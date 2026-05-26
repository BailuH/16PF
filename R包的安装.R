# 1. 设置 CRAN 镜像（避免交互式选择，建议加到 ~/.Rprofile）
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# 2. 安装步骤3需要的核心包（一次搞定）
install.packages(c(
  "psych",           # EFA 全能包：平行分析、VSS、MAP、碎石图
  "EFA.dimensions",  # Hull method
  "EGAnet",          # 探索图分析 EGA
  "polycor",         # Polychoric 相关矩阵
  "GPArotation",     # 旋转方法补充
  "lavaan",          # 如果用 ML 提取，需要拟合指数
  "ggplot2",         # 可视化
  "dplyr"            # 数据处理（R 版的 pandas）
))
