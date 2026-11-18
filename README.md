# Fig1
library(tidyverse)
library(patchwork)
library(ggbreak)
library(scales)
library(grid)
library(gridExtra)

# 2. 读取和处理数据
# ----------------------------------------------------------------------------
# 病毒科数据
virus_data <- data.frame(
  Family = c("Asfarviridae", "Corticoviridae", "Inoviridae", "Iridoviridae", 
             "Lavidaviridae", "Mimiviridae", "Phycodnaviridae", "Tectiviridae", "Unclassified"),
  N1 = c(0.4517, 119.9007, 3.0724, 922.3715, 2622.2934, 4283.0017, 17004.088, 0.1344, 5961834.083),
  N2 = c(0, 155.5145, 67.5185, 992.5232, 2747.3746, 4422.3179, 17350.4801, 0.1283, 5959931.497),
  N3 = c(0, 63.0825, 2.0279, 1158.0371, 3414.5848, 4890.9867, 19064.0558, 0.308, 5961694.223),
  Re = c(0, 0.1663, 56.4683, 31.6899, 1.2761, 470.1775, 169.3644, 0, 5988413.457),
  R1 = c(0, 5.7406, 219.244, 85.4376, 5.3135, 915.9295, 696.2958, 0, 5965385.211),
  R2 = c(0, 3803.7706, 1928.4635, 58.6277, 0.7426, 328.8376, 114.2519, 0, 5864039.253),
  R3 = c(0, 4972.1758, 1549.3035, 28.3134, 0.6228, 110.3821, 41.8608, 0, 5854346.599),
  We = c(385.8279, 1174.5795, 517.1175, 2056.2139, 775.0657, 1069.7999, 839.5661, 1655.2059, 5883381.22),
  W1 = c(192.5888, 470.5351, 3797.1967, 2261.9927, 320.5657, 633.3007, 665.1353, 0, 5352282.171),
  W2 = c(203.08, 232.825, 683.6177, 13022.2018, 286.2027, 1206.3287, 5107.3589, 0, 5708845.309),
  W3 = c(19.647, 115.595, 24185.6037, 18479.929, 68.5909, 1959.8854, 1947.9552, 0, 5588555.781)
)

# Shannon多样性数据
shannon_data <- data.frame(
  Sample = c("N1", "N2", "N3", "Re", "R1", "R2", "R3", "We", "W1", "W2", "W3"),
  Shannon = c(9.295591, 9.306378, 9.194032, 8.871925, 8.807752, 8.488196, 6.765096, 8.776596, 7.904375, 8.954394, 8.374754),
  Pielou = c(0.6211176, 0.6220485, 0.6151404, 0.594189, 0.5879481, 0.5643033, 0.4550897, 0.5996193, 0.5359294, 0.599072, 0.5649255),
  Group = c("N", "N", "N", "R", "R", "R", "R", "W", "W", "W", "W")
)

# 3. 设置样本顺序
# ----------------------------------------------------------------------------
sample_order <- c("N1", "N2", "N3", "Re", "R1", "R2", "R3", "We", "W1", "W2", "W3")
shannon_data$Sample <- factor(shannon_data$Sample, levels = sample_order)

# 4. 计算病毒科相对丰度
# ----------------------------------------------------------------------------
virus_relative <- virus_data %>%
  pivot_longer(cols = -Family, 
               names_to = "Sample", 
               values_to = "Abundance") %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = (Abundance / sum(Abundance, na.rm = TRUE)) * 100) %>%
  ungroup()

virus_relative$Sample <- factor(virus_relative$Sample, levels = sample_order)

# 5. 准备病毒科堆叠数据（带断点）
# ----------------------------------------------------------------------------
virus_complete <- virus_relative %>%
  group_by(Sample) %>%
  arrange(Relative_Abundance) %>%  # 按丰度升序排列
  mutate(
    Cum_Position = cumsum(Relative_Abundance),
    Start_Position = Cum_Position - Relative_Abundance
  ) %>%
  ungroup()

# 6. 设置颜色
# ----------------------------------------------------------------------------
specified_colors <- c("#fb6a4b", "#fcbca1", "#fdbb84", "#9fcbe2", 
                      "#bbbddc", "#9e9bc1", "#4292c7", "#acde8b", "#cdcdcd")

# 设置病毒科顺序
family_order <- c("Asfarviridae", "Corticoviridae", "Inoviridae", "Iridoviridae", 
                  "Lavidaviridae", "Mimiviridae", "Phycodnaviridae", "Tectiviridae", "Unclassified")

virus_complete$Family <- factor(virus_complete$Family, levels = family_order)

# 创建颜色映射
family_colors <- setNames(specified_colors, family_order)

# 7. 创建带断点的病毒科柱状图（左侧）
# ----------------------------------------------------------------------------
p_virus <- ggplot(virus_complete, aes(y = Sample)) +
  # 绘制所有病毒的堆叠柱状图
  geom_rect(
    aes(
      xmin = Start_Position,
      xmax = Cum_Position,
      ymin = as.numeric(Sample) - 0.35,
      ymax = as.numeric(Sample) + 0.35,
      fill = Family
    )
  ) +
  scale_fill_manual(values = family_colors, name = "病毒科") +
  # 添加断点：在1%和99.5%之间
  scale_x_break(c(1, 99.5), scales = 0.8) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1, 99.5, 99.75, 100),
    labels = c("0%", "0.5%", "1%", "99.5%", "99.75%", "100%"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_discrete(
    limits = rev(levels(virus_relative$Sample)),  # 反转顺序，使N1在顶部
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "病毒科相对丰度分布",
    x = "相对丰度 (%)",
    y = "样本"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "none",  # 移除图例
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.2),
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")
  ) +
  # 添加分组标签
  annotate("text", x = -0.05, y = "N2", label = "N组", 
           size = 4, fontface = "bold", color = "blue", hjust = 1) +
  annotate("text", x = -0.05, y = "R2", label = "R组", 
           size = 4, fontface = "bold", color = "darkgreen", hjust = 1) +
  annotate("text", x = -0.05, y = "W2", label = "W组", 
           size = 4, fontface = "bold", color = "darkred", hjust = 1) +
  coord_cartesian(xlim = c(0, 100), clip = "off")

# 8. 创建Shannon多样性柱状图（右侧）
# ----------------------------------------------------------------------------
p_shannon <- ggplot(shannon_data, aes(y = Sample, x = Shannon)) +
  geom_bar(stat = "identity", fill = "#FF8C00", width = 0.7) +
  geom_text(aes(label = round(Shannon, 2)), hjust = -0.2, size = 3.5) +
  scale_x_continuous(
    limits = c(0, 10),
    breaks = seq(0, 10, 2),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_y_discrete(
    limits = rev(levels(shannon_data$Sample)),  # 与病毒科图保持一致
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Shannon多样性指数",
    x = "Shannon多样性指数",
    y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.major.x = element_line(color = "gray90", size = 0.2),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.2),
    panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm")
  )
