#!/usr/bin/env Rscript
# =============================================================
# FEISTY 출력 비교: Reference (COBALT) vs NEMO-PISCES
# - quarto/FESITY methods.qmd Figure 4 코드 기반
# - Mollweide 투영 + viridis log10 스케일
# =============================================================

library(tidyverse)
library(sf)
library(viridis)
library(rnaturalearth)
library(patchwork)

# =============================================================
# 파일 경로
# =============================================================
BASE_DIR <- "/data01/labdisk/sungjin/NEMO-FEISTY/04.run"
OUT_DIR  <- "/data01/labdisk/sungjin/NEMO-FEISTY/05.figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Reference: quarto/data/ 에서 복사
load(file.path(BASE_DIR, "data", "Ref_fish_biomass.RData"))   # → out
df_ref <- out

# NEMO: RData
load(file.path(BASE_DIR, "data", "NEMO_fish_biomass.RData"))  # → out
df_nemo <- out

cat("Reference:", nrow(df_ref), "points\n")
cat("NEMO:     ", nrow(df_nemo), "points\n")

# =============================================================
# 데이터 전처리
# =============================================================
# lon 보정 (0~360 → -180~180)
df_ref$lon  <- ifelse(df_ref$lon > 179.5,  df_ref$lon - 360,  df_ref$lon)
df_nemo$lon <- ifelse(df_nemo$lon > 179.5, df_nemo$lon - 360, df_nemo$lon)

# Total biomass 계산
fish_cols <- c("totB_smpel", "totB_mesopel", "totB_largepel", "totB_midwpred", "totB_dem")
df_ref$totB_all  <- rowSums(df_ref[, fish_cols])
df_nemo$totB_all <- rowSums(df_nemo[, fish_cols])

# 세계 지도
world <- ne_countries(scale = 110, returnclass = "sf")

# =============================================================
# sf 격자 생성 함수
# =============================================================
make_grid <- function(df) {
  sf_pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  grid <- st_as_sf(
    sf_pts %>% st_make_grid(cellsize = 1, what = "polygons")
  ) %>%
    st_join(sf_pts, join = st_intersects, left = TRUE)
  return(grid)
}

# =============================================================
# 공통 플롯 테마
# =============================================================
map_theme <- theme(
  plot.background   = element_blank(),
  panel.background  = element_blank(),
  axis.text.y       = element_blank(),
  axis.text.x       = element_blank(),
  axis.ticks        = element_blank(),
  axis.title.y      = element_blank(),
  axis.title.x      = element_blank(),
  panel.border      = element_blank(),
  legend.position   = "bottom",
  legend.title      = element_text(size = 8),
  plot.title        = element_text(size = 10, hjust = 0.5)
)

# =============================================================
# 단일 변수 맵 생성 함수
# =============================================================
plot_map <- function(grid, var, title_text) {
  ggplot(grid) +
    geom_sf(aes(fill = log10(.data[[var]] + 0.01)), colour = NA) +
    scale_fill_viridis(
      name   = bquote(atop("Biomass g m"^-2, "")),
      labels = c("<0.01", "0.1", "1", "10", "100"),
      breaks = c(-2, -1, 0, 1, 2),
      limits = c(-2.5, 2.5),
      oob    = scales::squish
    ) +
    geom_sf(data = world, col = "grey50", fill = "grey80") +
    coord_sf(crs = "+proj=moll") +
    ggtitle(title_text) +
    map_theme
}

# =============================================================
# [1] Total biomass: Reference vs NEMO
# =============================================================
cat("\n[1] Total biomass comparison map...\n")

grid_ref  <- make_grid(df_ref)
grid_nemo <- make_grid(df_nemo)

p1 <- plot_map(grid_ref,  "totB_all", "Reference (COBALT)")
p2 <- plot_map(grid_nemo, "totB_all", "NEMO-PISCES")

p_total <- p1 + p2 +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "feisty_total_biomass_comparison.png"),
       p_total, width = 8, height = 9, dpi = 150)
cat("  저장:", file.path(OUT_DIR, "feisty_total_biomass_comparison.png"), "\n")

# =============================================================
# [2] 그룹별 비교 맵
# =============================================================
cat("\n[2] Fish group comparison maps...\n")

group_info <- data.frame(
  var   = fish_cols,
  label = c("Small pelagic", "Mesopelagic", "Large pelagic",
            "Midwater predator", "Demersal"),
  stringsAsFactors = FALSE
)

plot_list <- list()
for (i in seq_len(nrow(group_info))) {
  var   <- group_info$var[i]
  label <- group_info$label[i]

  p_ref  <- plot_map(grid_ref,  var, paste0("Ref: ", label))
  p_nemo <- plot_map(grid_nemo, var, paste0("NEMO: ", label))

  plot_list[[2 * i - 1]] <- p_ref
  plot_list[[2 * i]]     <- p_nemo
}

p_groups <- wrap_plots(plot_list, ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "feisty_groups_comparison.png"),
       p_groups, width = 14, height = 20, dpi = 150)
cat("  저장:", file.path(OUT_DIR, "feisty_groups_comparison.png"), "\n")

# =============================================================
# [3] 바 차트: 그룹별 평균 비교
# =============================================================
cat("\n[3] Bar chart...\n")

df_bar <- data.frame(
  group = rep(group_info$label, 2),
  source = rep(c("Reference", "NEMO-PISCES"), each = nrow(group_info)),
  mean_biomass = c(
    sapply(fish_cols, function(v) mean(df_ref[[v]], na.rm = TRUE)),
    sapply(fish_cols, function(v) mean(df_nemo[[v]], na.rm = TRUE))
  )
)
df_bar$group <- factor(df_bar$group, levels = group_info$label)

p_bar <- ggplot(df_bar, aes(x = group, y = mean_biomass, fill = source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c("Reference" = "steelblue", "NEMO-PISCES" = "tomato")) +
  labs(y = expression("Mean biomass [g m"^-2*"]"),
       title = "Mean Fish Biomass by Functional Group",
       fill = NULL) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        legend.position = "top")

ggsave(file.path(OUT_DIR, "feisty_bar_comparison.png"),
       p_bar, width = 9, height = 5, dpi = 150)
cat("  저장:", file.path(OUT_DIR, "feisty_bar_comparison.png"), "\n")

# =============================================================
# [4] 위도별 평균 (Zonal mean)
# =============================================================
cat("\n[4] Zonal mean...\n")

zonal_ref <- df_ref %>%
  mutate(lat_bin = round(lat)) %>%
  group_by(lat_bin) %>%
  summarise(across(all_of(c(fish_cols, "totB_all")), ~mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(source = "Reference")

zonal_nemo <- df_nemo %>%
  mutate(lat_bin = round(lat)) %>%
  group_by(lat_bin) %>%
  summarise(across(all_of(c(fish_cols, "totB_all")), ~mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(source = "NEMO-PISCES")

zonal <- bind_rows(zonal_ref, zonal_nemo)

all_vars <- c(fish_cols, "totB_all")
all_labels <- c(group_info$label, "Total biomass")

zonal_plots <- list()
for (i in seq_along(all_vars)) {
  zonal_plots[[i]] <- ggplot(zonal, aes(x = lat_bin, y = .data[[all_vars[i]]],
                                         colour = source)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = c("Reference" = "steelblue", "NEMO-PISCES" = "tomato")) +
    labs(title = all_labels[i], x = "Latitude", y = expression("g m"^-2)) +
    xlim(-90, 90) +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = "none")
}

p_zonal <- wrap_plots(zonal_plots, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "feisty_zonal_comparison.png"),
       p_zonal, width = 14, height = 8, dpi = 150)
cat("  저장:", file.path(OUT_DIR, "feisty_zonal_comparison.png"), "\n")

# =============================================================
# [5] 요약 통계
# =============================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("%18s | %10s %10s | %10s %10s | %6s\n",
            "Group", "Ref mean", "Ref max", "NEMO mean", "NEMO max", "Ratio"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in seq_along(all_vars)) {
  rv <- df_ref[[all_vars[i]]]
  nv <- df_nemo[[all_vars[i]]]
  ratio <- if (mean(rv, na.rm = TRUE) > 0) mean(nv, na.rm = TRUE) / mean(rv, na.rm = TRUE) else 0
  cat(sprintf("%18s | %10.3f %10.3f | %10.3f %10.3f | %5.2fx\n",
              all_labels[i],
              mean(rv, na.rm = TRUE), max(rv, na.rm = TRUE),
              mean(nv, na.rm = TRUE), max(nv, na.rm = TRUE),
              ratio))
}
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\n완료! 그림 →", OUT_DIR, "\n")
