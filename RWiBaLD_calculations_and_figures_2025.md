# Range Weighted Branch Length Difference (RWiBaLD) Calculations &
Figure Generation
Nunzio Knerr
2026-02-26

- [<span class="toc-section-number">1</span>
  Introduction](#introduction)
- [<span class="toc-section-number">2</span> Load
  packages](#load-packages)
- [<span class="toc-section-number">3</span> The RWiBaLD
  metric](#the-rwibald-metric)
- [<span class="toc-section-number">4</span> Specify biodiverse results
  data](#specify-biodiverse-results-data)
- [<span class="toc-section-number">5</span> Load data & calculate
  RWiBaLD](#load-data--calculate-rwibald)
- [<span class="toc-section-number">6</span> Elbow point calculation
  function](#elbow-point-calculation-function)
- [<span class="toc-section-number">7</span> Calculate branches of
  interest](#calculate-branches-of-interest)
- [<span class="toc-section-number">8</span> Generate figure 2A -
  RWiBaLD branches of
  interest](#generate-figure-2a---rwibald-branches-of-interest)
- [<span class="toc-section-number">9</span> Calculate RWiBaLD branch
  categories (neo-endemic, meso-endemic,
  paleo-endemic)](#calculate-rwibald-branch-categories-neo-endemic-meso-endemic-paleo-endemic)
- [<span class="toc-section-number">10</span> Load RWiBaLD results &
  tree data by cells](#load-rwibald-results--tree-data-by-cells)
- [<span class="toc-section-number">11</span> Generate figure 2B -
  ranked RWiBaLD score by RWiBaLD
  score](#generate-figure-2b---ranked-rwibald-score-by-rwibald-score)
- [<span class="toc-section-number">12</span> Generate figure 2C -
  ranked BaLD by BaLD (non range
  weighted)](#generate-figure-2c---ranked-bald-by-bald-non-range-weighted)
- [<span class="toc-section-number">13</span> Generate multi-panel
  figure 2ABC](#generate-multi-panel-figure-2abc)
- [<span class="toc-section-number">14</span> Generate Figure
  3A](#generate-figure-3a)
- [<span class="toc-section-number">15</span> Generate Figure
  3B](#generate-figure-3b)
- [<span class="toc-section-number">16</span> Generate multi-panel
  figure 3AB](#generate-multi-panel-figure-3ab)
- [<span class="toc-section-number">17</span> Read tree file & generate
  supplementary figure
  1](#read-tree-file--generate-supplementary-figure-1)
- [<span class="toc-section-number">18</span> Generate tree in 2 halves
  supplementary figure 1A and
  1B](#generate-tree-in-2-halves-supplementary-figure-1a-and-1b)
- [<span class="toc-section-number">19</span> Load data & generate
  histograms for all
  cells](#load-data--generate-histograms-for-all-cells)
- [<span class="toc-section-number">20</span> Load data & generate 12
  histograms for figure
  4A-L](#load-data--generate-12-histograms-for-figure-4a-l)
- [<span class="toc-section-number">21</span> Generate multi-panel
  histogram figure 4A-L](#generate-multi-panel-histogram-figure-4a-l)
- [<span class="toc-section-number">22</span> Load results data &
  generate data files for biome
  maps](#load-results-data--generate-data-files-for-biome-maps)
- [<span class="toc-section-number">23</span> Two tailed relative
  phylogenetic diversity (RPD) & CANAPE
  functions](#two-tailed-relative-phylogenetic-diversity-rpd--canape-functions)
- [<span class="toc-section-number">24</span> Load CANAPE &
  randomisation results](#load-canape--randomisation-results)
- [<span class="toc-section-number">25</span> Generate CANAPE map figure
  5A](#generate-canape-map-figure-5a)
- [<span class="toc-section-number">26</span> Generate biome map with
  neo-endemics figure
  5B](#generate-biome-map-with-neo-endemics-figure-5b)
- [<span class="toc-section-number">27</span> Generate biome map with
  meso-endemics figure
  5C](#generate-biome-map-with-meso-endemics-figure-5c)
- [<span class="toc-section-number">28</span> Generate biome map with
  paleo-endemics figure
  5D](#generate-biome-map-with-paleo-endemics-figure-5d)
- [<span class="toc-section-number">29</span> Generate 4 up map figure
  5ABCD](#generate-4-up-map-figure-5abcd)
- [<span class="toc-section-number">30</span> Table of RWiBaLD
  results supplementary table
  1](#table-of-rwibald-results-supplementary-table-1)

------------------------------------------------------------------------

## Introduction

This document contains code for calculating and plotting Range Weighted
Branch Length Difference (RWiBaLD). It uses input from
[Biodiverse](https://github.com/shawnlaffan/biodiverse) analyses of
phylogenetic data. The method is related to the [Categorical Analysis of
Neo And Paleo-Endemisim
(CANAPE)](https://www.nature.com/articles/ncomms5473) published in 2014
and follows up on that paper. Consequently, the [acacia dataset from the
CANAPE paper](https://doi.org/10.5061/dryad.dv4qk) is used to illustrate
the method.

## Load packages

Here we load the packages for use in subsequent code.

<details class="code-fold">
<summary>Show the code</summary>

``` r
library(sf)
library(ggplot2)
library(Cairo)
library(extrafont)
library(dplyr)
library(ggrepel)
library(stringr)
library(patchwork)
library(forcats)
library(ggtree)
library(ape)
library(knitr)
library(gt)
library(colorBlindness)
```

</details>

## The RWiBaLD metric

The steps for calculating RWiBaLD are as follows:

1.  For each branch on the phylogeny (terminal or internal), calculate
    its *range weighted branch length difference* (RWiBaLD) score as the
    difference between its length on the RWoT and the RWcT. These scores
    can be used directly as a continuous measure of neo- and
    paleo-endemism (negative and positive values, respectively).
    Classification into neo-endemic, paleo-endemic, and meso-endemic
    requires further processing using steps 2-4:  

2.  Identify the “highly endemic branches,” i.e., those with range sizes
    below a threshold. We calculate this threshold by identifying the
    “elbow” of the distribution, using the *maximum Euclidean distance*
    method described by [Ramer
    (1972)](https://doi.org/10.1016/S0146-664X(72)80017-0). This is done
    by ranking all branches according to the inverse of their range size
    (functionally the same as their lengths on the RWcT), plotting a
    straight line from the first to last points on that curve, then
    identifying the point on the curve that produces the longest
    perpendicular line drawn to it from the straight line (Fig 2A).
    Branches greater than or equal to the threshold are considered the
    highly endemic branches of interest.

3.  The same elbow statistic is then applied to the RWiBaLD scores from
    step 1. The positive and negative differences are processed
    separately (divided at the RWiBaLD = 0 value), resulting in two
    thresholds (Fig 2B).

4.  The highly endemic branches identified in step 2 are then classified
    using the thresholds from step 3. Highly endemic branches with
    negative differences less than or equal to the threshold are
    classified as neo-endemic, those with positive differences greater
    than or equal to the threshold are classified as paleo-endemic, and
    those between these two categories are classified as meso-endemic
    (Fig 2B). 

## Specify biodiverse results data

Specify the tabular data files exported from a biodiverse analysis. Both
the observed data & the equal branch length data from a range weighted
tree.

<details class="code-fold">
<summary>Show the code</summary>

``` r
# Set the data locations etc.
data_dir <- "Acacia_biodiverse_exports/"
observed_data_csv <- paste0(data_dir, "Acacia_RWT_tabular_export.csv")
equal_branch_length_data_csv  <- paste0(data_dir, "Acacia_RWT_EQBL_tabular_export.csv")
```

</details>

## Load data & calculate RWiBaLD

Load the data in, merge the tables and calculate the RWiBaLD score for
each branch, then write out the results to a file.

<details class="code-fold">
<summary>Show the code</summary>

``` r
# Load data
observed_data <- read.table(observed_data_csv, header=T,sep=",")
equal_branch_length_data <- read.table(equal_branch_length_data_csv, header=T,sep=",")

# Calculate RWiBaLD Statistic etc.
rwibald_results <- observed_data |>
  left_join(equal_branch_length_data, by = 'NAME') |>
  rename(branch_length_observed_tree = LENGTH.x, branch_length_comparison_tree = LENGTH.y) |> 
  mutate(rank_comparison_tree = rank(desc(branch_length_comparison_tree), ties.method = "min")) |>
  mutate(sorted_asc_bl_comparison_tree = rank(branch_length_comparison_tree, ties.method = "first")) |>
  mutate(rwibald_score = (branch_length_observed_tree - branch_length_comparison_tree)) |>
  mutate(rwibald_score_rank = rank(rwibald_score, ties.method = "first"))

# View(rwibald_results)
output_dir <- "quarto_outputs/"

# Write data to file
write.csv(rwibald_results, paste0(output_dir, "Acacia_RWiBaLD_results.csv"), row.names=FALSE)
```

</details>

## Elbow point calculation function

We propose a function to calculate the elbow point of a curve in order
to determine branches of interest and identify the different categories
of RWiBaLD (e.g., neo, meso, paleo).

It takes a numerical vector `data` as input and calculates the point
where the rate of decrease in a measure of fit (e.g., sum of squared
distances) starts to slow down. The function first orders the values of
the data and calculates vectors from the first point to all other
points, normalizing the vector from the first to the last point. It then
computes the Euclidean distance of each point to this line and returns
the y-coordinate (value) of the point with the maximum distance from the
line, considered the “elbow” point threshold.

The procedure is described step by step below:

1.  **Function Definition and Input:**

    - The function `get_elbow` takes one argument, `data`, which is a
      numerical vector.

2.  **Initial Setup:**

    - `n_pts` is assigned the length of `data`.

    - `x_coords` is a sequence from 1 to `n_pts`.

    - `y_coords` is the values of `data`, sorted in ascending order.

3.  **First Point Coordinates:**

    - `x1` and `y1` are the coordinates of the first point in the
      sequence.

4.  **Normalize the Line Vector:**

    - `x_vec` and `y_vec` are vectors from the first point to all other
      points.

    - `x_vec_max` and `y_vec_max` are the vectors from the first point
      to the last point.

    - `normaliser` is the length of the line vector.

    - Normalize `x_vec_max` and `y_vec_max` by dividing by `normaliser`.

5.  **Vectors from First Point:**

    - `v_x` and `v_y` are vectors from the first point to all other
      points.

    - `scalar_prod` calculates the scalar projection of `v_x` and `v_y`
      onto the normalized line vector.

6.  **Calculate Distance to Line:**

    - `vec_to_line_x` and `vec_to_line_y` are the components of the
      vectors from each point to the line.

    - `dist_to_line` is the Euclidean distance of each point to the
      line.

    - `i_max` is the index of the maximum distance (the elbow)

7.  **Return the Elbow Point Threshold:**

    - The function returns the coordinates (x,y) of the point with the
      maximum distance from the line, considered the “elbow” point
      threshold.

In summary, the function identifies the point in the data set that is
farthest away from a line drawn between the first and last data points,
which often corresponds to a significant change. In case of ties we take
the point closest to the minimum.

see: Ramer, U. (1972). *An iterative procedure for the polygonal
approximation of plane curves.* **Computer Graphics and Image
Processing**, 1, 244–256.
<https://doi.org/10.1016/S0146-664X(72)80017-0>

<details class="code-fold">
<summary>Show the code</summary>

``` r
get_elbow <- function(data){
    
   n_pts = length(data)

   x_coords = 1:n_pts
   y_coords = data[order(data)]

   # First points
   x1 = x_coords[1]
   y1 = y_coords[1]

   # Normalize the line vector
   x_vec = x_coords - x1
   y_vec = y_coords - y1
   x_vec_max = x_vec[n_pts]
   y_vec_max = y_vec[n_pts]

   normaliser = sqrt(x_vec_max^2 + y_vec_max^2)

   x_vec_max = x_vec_max / normaliser
   y_vec_max = y_vec_max / normaliser

   # Vectors from first point
   v_x = x_coords - x1
   v_y = y_coords - y1

   scalar_prod = v_x * x_vec_max + v_y * y_vec_max

   vec_to_line_x = v_x - scalar_prod * x_vec_max
   vec_to_line_y = v_y - scalar_prod * y_vec_max

   # Distance to line is the norm
   dist_to_line = sqrt(vec_to_line_x^2 + vec_to_line_y^2)
   
   # y value at the point of maximum distance from the hypotenuse
   #return(y_coords[which(dist_to_line==max(dist_to_line))])
   
   # Find index of the maximum distance (the elbow)
   i_max <- which.max(dist_to_line)
  
   # Return both coordinates
   return(list(x = x_coords[i_max], y = y_coords[i_max]))
}
```

</details>

## Calculate branches of interest

Now we take the data and determine the key branches of interest using
the elbow point statistic defined in the function above.

<details class="code-fold">
<summary>Show the code</summary>

``` r
# For range weighted tree
df <- rwibald_results[,c("NAME", "rwibald_score", "rwibald_score_rank",
                         "sorted_asc_bl_comparison_tree",
                         "branch_length_comparison_tree",
                         "branch_length_observed_tree")] %>%
      arrange(., branch_length_comparison_tree)

# Call the get_elbow function to find the threshold
elbow_thresh <- get_elbow(df$branch_length_comparison_tree)

# Define key branches
key_branches <- df$branch_length_comparison_tree >= elbow_thresh$y

# Count them
number_of_key_branches <- sum(key_branches)

# Add a column to the results marking the branches of interest as 'sig' 
rwibald_results_key <- rwibald_results %>%
  arrange(., branch_length_comparison_tree) %>%
  mutate(bl_comparison_tree_rwibald_elbow_key = ifelse(key_branches, 'key', 'not_key'))

# Merge data
rwibald_results_all <- merge(rwibald_results, rwibald_results_key[, c("NAME","bl_comparison_tree_rwibald_elbow_key")], by=c('NAME','NAME'),all.x=T)

#View(rwibald_results_all)

# Print number of key branches
print(paste0("Number of key branches: ", number_of_key_branches))
```

</details>

    [1] "Number of key branches: 135"

## Generate figure 2A - RWiBaLD branches of interest

Now we plot the ranked branch lengths on the comparison tree showing the
cut-off line for the key RWiBaLD branches of interest.

<details class="code-fold">
<summary>Show the code</summary>

``` r
# Check the total number of branches
#length(df$branch_length_comparison_tree)

# Add x_marks every 100 
x_marks <- c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, length(df$branch_length_comparison_tree))

branch_length_by_rank <- ggplot(
  df,
  aes(x = as.numeric(rownames(df)), y = sort(branch_length_comparison_tree))
) +
  geom_point() +
  geom_vline(
    xintercept = (length(df$branch_length_comparison_tree) - number_of_key_branches + 1),
    color = "red", linetype = "dashed"
  ) +
  geom_segment(x = 0, y = 0, xend = 1015, yend = 0.0090413023,
    color = "orange", linewidth = 0.5, inherit.aes = FALSE) +
  geom_segment(x = 881, y = 0.001004589, xend = 755.55, yend = 0.00675, color = "orange", linewidth = 0.5, inherit.aes = FALSE) +  
  annotate("text", x = 724, y = 0.0052, label = "longest\nperpendicular\nline",
           color = "orange") +
  xlab("Branch Length Rank (shortest-longest)") +
  ylab("Branch Length\nComparison Tree (RWcT)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = x_marks, labels = x_marks) +
  annotate("text", x = 945, y = 0.0072,
           label = paste0(number_of_key_branches, " key\n branches\n of interest"),
           vjust = 1, hjust = 0.5)


print(branch_length_by_rank)
```

</details>

![](RWiBaLD_calculations_and_figures_2025_files/figure-commonmark/unnamed-chunk-6-1.png)

<details class="code-fold">
<summary>Show the code</summary>

``` r
output_dir <- "quarto_outputs/"

#cvdPlot(branch_length_by_rank)

ggsave(paste0(output_dir, "figures/Figure2_A.png"), plot = branch_length_by_rank,  scale = 1,
        width = 3000,
        height = 1354,
        units = "px",
        dpi = 300)
```

</details>

## Calculate RWiBaLD branch categories (neo-endemic, meso-endemic, paleo-endemic)

Next we take the key branches of interest and define which of them are
neo-endemic, meso-endemic or paleo-endemic. We use the same elbow
statistic defined above but only consider the branches of key interest
calculated above as valid. We split the data at 0 with those equal or
below the elbow point threshold in the negative data being neo-endemics.
Those in the positive dataset equal or above the elbow point threshold
being classified as paleo-endemic, and rest are classified as
meso-endemic. We have described three categories, but of course one
could break the continuous distribution into further categories if
required.

<details class="code-fold">
<summary>Show the code</summary>

``` r
# Load the data
data <- rwibald_results_all[,c("NAME", "bl_comparison_tree_rwibald_elbow_key", "rwibald_score", "rwibald_score_rank")] %>%
  arrange(., rwibald_score_rank)

# Subset the left and right dataframes
subset_left_df <- data[data$rwibald_score <= 0, ]
subset_right_df <- data[data$rwibald_score >= 0, ]

# Create a list to store the results for the table
results <- list()

# Records in negative data
negative_records <- nrow(subset_left_df)
results$Negative_Records <- negative_records

# Neo cutoff
left_point_x <- get_elbow(subset_left_df$rwibald_score)
rwibald_score_rank_at_left_point <- subset_left_df[subset_left_df$rwibald_score == left_point_x$y, "rwibald_score_rank"]
results$Neo_Cutoff <- rwibald_score_rank_at_left_point

# Plot for left half of curve
l <- ggplot(subset_left_df, aes(x = rwibald_score_rank, y = rwibald_score)) +
  geom_point() +
  geom_vline(xintercept = rwibald_score_rank_at_left_point, color = "red", linetype = "dashed") +
  ggtitle(label= "Neo cutoff") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

print(l)
```

</details>

![](RWiBaLD_calculations_and_figures_2025_files/figure-commonmark/unnamed-chunk-7-1.png)

<details class="code-fold">
<summary>Show the code</summary>

``` r
ggsave(paste0(output_dir, "figures/Figure2B_L.png"), plot = l,  scale = 1, 
       width = 2048, 
       height = 1024, 
       units = "px", 
       dpi = 300)

# Records in positive data
positive_records <- nrow(subset_right_df)
results$Positive_Records <- positive_records

# Paleo cutoff
right_point_x <- get_elbow(subset_right_df$rwibald_score)
rwibald_score_rank_at_right_point <- subset_right_df[subset_right_df$rwibald_score == right_point_x$y, "rwibald_score_rank"]
results$Paleo_Cutoff <- rwibald_score_rank_at_right_point

# Plot for right half of curve
r <- ggplot(subset_right_df, aes(x = rwibald_score_rank, y = rwibald_score)) +
  geom_point() +
  geom_vline(xintercept = rwibald_score_rank_at_right_point, color = "red", linetype = "dashed") +
  ggtitle(label= "Paleo cutoff") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
print(r)
```

</details>

![](RWiBaLD_calculations_and_figures_2025_files/figure-commonmark/unnamed-chunk-7-2.png)

<details class="code-fold">
<summary>Show the code</summary>

``` r
ggsave(paste0(output_dir, "figures/Figure2B_R.png"), plot = r,  scale = 1,
        width = 2048,
        height = 1024,
        units = "px",
        dpi = 300)

# Create a column with the results
rwibald_results_key_type <- rwibald_results_all %>%
  filter(bl_comparison_tree_rwibald_elbow_key == "key") %>%
  mutate(rwibald_type = case_when(rwibald_score_rank <= rwibald_score_rank_at_left_point ~ "neo-endemic", 
                                  rwibald_score_rank >= rwibald_score_rank_at_right_point ~ "paleo-endemic",
                                  TRUE ~ "meso-endemic")) %>%
  select(NAME, rwibald_type)

# ---- Table output: gt for HTML/PDF, kable for GFM ----
summary_results_df <- as.data.frame(t(sapply(results, c)))

if (knitr::is_html_output() || knitr::is_latex_output()) {

  results_gt <- gt::gt(data = summary_results_df) |>
    gt::tab_header(title = "Summary of Records and Cutoffs") |>
    gt::fmt_number(columns = dplyr::everything(), decimals = 0) |>
    gt::tab_style(
      style = gt::cell_text(align = "center"),
      locations = gt::cells_body()
    )

  results_gt

} else {

  knitr::kable(
    summary_results_df,
    format = "pipe",
    align  = "c"
  )

}

# results_gt <- gt(data = summary_results_df)
# 
# # Display the table
# results_gt %>%
#   tab_header(title = "Summary of Records and Cutoffs") %>%
#   fmt_number(columns = everything(), decimals = 0) %>%
#   tab_style(style = cell_text(align = "center"),
#             locations = cells_body()
#             )
# View(rwibald_results_key_type)

# Merge the data
rwibald_results_all <- merge(rwibald_results_all, rwibald_results_key_type[, c("NAME","rwibald_type")], by=c('NAME','NAME'),all.x=T)

# View(rwibald_results_all)
output_dir <- "quarto_outputs/"

# Write data to file
write.csv(rwibald_results_all,paste0(output_dir, "Acacia_RWiBaLD_results_all.csv") ,row.names=FALSE)
```

</details>

| Negative_Records | Neo_Cutoff | Positive_Records | Paleo_Cutoff |
|:----------------:|:----------:|:----------------:|:------------:|
|       686        |     71     |       330        |     991      |

## Load RWiBaLD results & tree data by cells

Now we load the data for the RWiBaLD results and the range weighted tree
data in tabular format and then replace the counts in the matrix with
the branch lengths for that group so we can generate the histograms for
each grid cell location. We also add the CANAPE scores to the data table
for comparison and calculate the ranges of each branch.

<details class="code-fold">
<summary>Show the code</summary>

``` r
data_dir <- "Acacia_biodiverse_exports/"
output_dir <- "quarto_outputs/"

# Data with internal branches (spatial anaylsis from BD) x grid cells
all_tree_data_csv <- paste0(data_dir, "Acacia_PD_Included_Node_List.csv")
RWiBaLD_results_csv  <- paste0(output_dir, "Acacia_RWiBaLD_results_all.csv")

all_tree_data <- read.table(all_tree_data_csv, header=T, sep=",", check.names = FALSE )
RWiBaLD_results <- read.table(RWiBaLD_results_csv, header=T, sep=",")

#View(all_tree_data)
#View(RWiBaLD_results)
valueToUse <- "rwibald_score"

# Iterate over the columns in Table 1
for (col in colnames(all_tree_data)[-c(1:3)]) {
  # get the corresponding "valueToUse" column specified above from Table 2
  branch_length <- RWiBaLD_results[RWiBaLD_results$NAME == col, paste0(valueToUse)]
  # replace the numbers in Table 1 with the "valueToUse" value
  all_tree_data[, col] <- ifelse(is.na(all_tree_data[, col]), NA, as.numeric(branch_length))
}

# View(all_tree_data)

# Write data to file
write.csv(all_tree_data, paste0(output_dir, "Acacia_PD_Included_Node_List_", valueToUse, ".csv"), row.names=FALSE)

# Get CANAPE data from biodiverse as well
canape_csv <- paste0(data_dir, "Acacia_Rand1_CANAPE_Export.csv")
canape_results <- read.table(canape_csv, header=T,sep=",")

# View(canape_results)

canape_group_data <- canape_results %>%
                     left_join(all_tree_data %>%
                               select(-c(Axis_0, Axis_1)), by = "ELEMENT"
                               )

# Transpose the dataframe
data <- t(canape_group_data)

# View(data)

# Set the first row as column names
transposed_df <- setNames(data.frame(data[-1,]), data[1,])

# View(transposed_df)

# Merge CANAPE and RWiBaLD data
canape_group_data_rwibald <- transposed_df %>%
  mutate(range_cell_count = rowSums(!is.na(.))) %>%
  tibble::rownames_to_column(var = "NAME") %>%
  left_join(select(rwibald_results_all, NAME, rwibald_type), by = 'NAME') %>%
  mutate(rwibald_type = ifelse(is.na(rwibald_type), "other", rwibald_type)) 

# Set row names
rownames(canape_group_data_rwibald) <- canape_group_data_rwibald$NAME 

# View(canape_group_data_rwibald)

# Add range data to results and create three new columns
rwibald_results_all_with_range <- rwibald_results_all %>%
  left_join(select(canape_group_data_rwibald, NAME, range_cell_count), by = "NAME") %>%
  mutate(rwibald_type = replace(rwibald_type, is.na(rwibald_type), "other"),
         non_range_weighted_branch_length_observed_tree = branch_length_observed_tree * range_cell_count,
         non_range_weighted_branch_length_comparison_tree = branch_length_comparison_tree * range_cell_count,
         diff_non_range_weighted_observed_to_comparison_bl = non_range_weighted_branch_length_observed_tree - non_range_weighted_branch_length_comparison_tree
        ) %>%
  mutate(rank_BaLD = rank(diff_non_range_weighted_observed_to_comparison_bl, ties.method = "first"))


# Write data to file
write.csv(rwibald_results_all_with_range, paste0(output_dir, "Acacia_RWiBaLD_results_all_with_range.csv"), row.names=FALSE)

# View(rwibald_results_all_with_range)
# View(rwibald_results_all)
```

</details>

## Generate figure 2B - ranked RWiBaLD score by RWiBaLD score

Here we plot the ranked RWiBaLD score (x) by the RWiBaLD score (y)
(Figure 2B), colouring the different RWiBaLD Categories: neo-endemic
(red), meso-endemic (#FFD851) and paleo-endemic (royalblue1).

<details class="code-fold">
<summary>Show the code</summary>

``` r
# Find the neo_cutoff
neo_cutoff <- rwibald_results_all_with_range %>%
  filter(rwibald_type == "neo-endemic") %>%
  summarize(highest_rank = max(rwibald_score_rank, na.rm = TRUE)) %>%
  pull(highest_rank)

# Find the paleo_cutoff
paleo_cutoff <- rwibald_results_all_with_range %>%
  filter(rwibald_type == "paleo-endemic") %>%
  summarize(lowest_rank = min(rwibald_score_rank, na.rm = TRUE)) %>%
  pull(lowest_rank)

# Get only the rwibald branches of key interest
rwibald_results_all_key_only <- rwibald_results_all_with_range %>% 
                                filter(bl_comparison_tree_rwibald_elbow_key == "key")

# View(rwibald_results_all_key_only)
# View(rwibald_results_all_with_range)

# Plot rwibald rank X rwibald score
rwibald_key_plot <- ggplot() +
    geom_point(data=rwibald_results_all_with_range, aes(x = rwibald_score_rank, y = rwibald_score), fill="transparent", size = 3, colour= "grey77", alpha = 0.4, pch = 21) +
  geom_point(data = rwibald_results_all_key_only, aes(x = rwibald_score_rank, y = rwibald_score, fill = rwibald_type, color = rwibald_type), size = 3, alpha = 0.9, pch = 21) +
    scale_fill_manual(name = "RWiBaLD Category", values = c("paleo-endemic" = "royalblue1", "neo-endemic" = "red", "meso-endemic" = "#FFD851"), labels = c("neo-endemic","meso-endemic", "paleo-endemic"), limits = c("neo-endemic","meso-endemic", "paleo-endemic"),  na.value = "transparent") +
  scale_color_manual(name = "RWiBaLD Category", values = c("paleo-endemic" = "royalblue4", "neo-endemic" = "red4", "meso-endemic" = "#D0A100"), labels = c("neo-endemic", "meso-endemic", "paleo-endemic"), limits = c("neo-endemic", "meso-endemic", "paleo-endemic"), na.value = "transparent") +
  geom_vline(xintercept = neo_cutoff, color = "red", linetype = "dashed") +
  geom_vline(xintercept = paleo_cutoff, color = "red", linetype = "dashed") +
  xlab("RWiBaLD Score Rank (shortest-longest)") +
  ylab("Range Weighted Branch\nLength Difference (RWiBaLD) Score") +
 #annotate("segment", x = results$Negative_Records, xend = results$Negative_Records, 
  #          y = max(rwibald_results_all_with_range$rwibald_score) * 1.1, yend = 0, 
  #          arrow = arrow(length = unit(0.2, "cm")), color = "grey74") +
  geom_vline(xintercept = results$Negative_Records, color = "grey74", linetype = "solid") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom", legend.text = element_text(size = 12))

print(rwibald_key_plot)
```

</details>

![](RWiBaLD_calculations_and_figures_2025_files/figure-commonmark/unnamed-chunk-9-1.png)

<details class="code-fold">
<summary>Show the code</summary>

``` r
#cvdPlot(rwibald_key_plot)

ggsave("quarto_outputs/figures/Figure2_B.png", plot = rwibald_key_plot,  scale = 1,
        width = 3000,
        height = 1354,
        units = "px",
        dpi = 300)
```

</details>

## Generate figure 2C - ranked BaLD by BaLD (non range weighted)

Illustrating the distribution of branch lengths on the observed tree
with no range-weighting; the Y-axis is the difference between length on
the observed tree and length on the comparison tree for each branch; the
X-axis is ranked branch length difference (shortest - longest).  The
branches of interest are colored by RWiBaLD category as in 2B

<details class="code-fold">
<summary>Show the code</summary>

``` r
#View(rwibald_results_all_with_range)

#View(rwibald_results_all_key_only)

# Plot non_range_weighted score
rwibald_key_plot_non_range_weighted <- ggplot() +
    #geom_point(data=rwibald_results_all_with_range, aes(x = reorder(non_range_weighted_branch_length_observed_tree, non_range_weighted_branch_length_observed_tree), y = diff_non_range_weighted_observed_to_comparison_bl), fill="transparent", size = 3, colour= "grey77", alpha = 0.4, pch = 21) +
  geom_point(data=rwibald_results_all_with_range, aes(x = rank_BaLD, y = diff_non_range_weighted_observed_to_comparison_bl), fill="transparent", size = 3, colour= "grey77", alpha = 0.4, pch = 21) +
  geom_point(data = rwibald_results_all_key_only, aes(x =  rank_BaLD, y = diff_non_range_weighted_observed_to_comparison_bl, fill = rwibald_type, color = rwibald_type), size = 3, alpha = 0.9, pch = 21) +
  scale_fill_manual(name = "RWiBaLD Category", values = c("paleo-endemic" = "royalblue1", "neo-endemic" = "red", "meso-endemic" = "#FFD851"), labels = c("neo-endemic","meso-endemic", "paleo-endemic"), limits = c("neo-endemic","meso-endemic", "paleo-endemic"),  na.value = "transparent") +
  scale_color_manual(name = "RWiBaLD Category", values = c("paleo-endemic" = "royalblue4", "neo-endemic" = "red4", "meso-endemic" = "#D0A100"), labels = c("neo-endemic", "meso-endemic", "paleo-endemic"), limits = c("neo-endemic", "meso-endemic", "paleo-endemic"), na.value = "transparent") +
  #geom_vline(xintercept = neo_cutoff, color = "red", linetype = "dashed") +
  #geom_vline(xintercept = paleo_cutoff, color = "red", linetype = "dashed") +
  xlab("Unweighted Branch Length Difference (BaLD) Score Rank (shortest−longest)") +
  ylab("Branch Length\nDifference (BaLD) Score") +
  geom_vline(xintercept = results$Negative_Records, color = "grey74", linetype = "solid") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none", legend.text = element_text(size = 12), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8))

print(rwibald_key_plot_non_range_weighted)
```

</details>

![](RWiBaLD_calculations_and_figures_2025_files/figure-commonmark/unnamed-chunk-10-1.png)

<details class="code-fold">
<summary>Show the code</summary>

``` r
ggsave("quarto_outputs/figures/Figure2_C.png", plot = rwibald_key_plot_non_range_weighted,  scale = 1,
        width = 3000,
        height = 1354,
        units = "px",
        dpi = 300)
```

</details>

## Generate multi-panel figure 2ABC

Here we generate a 2 up figure using patchwork

<details class="code-fold">
<summary>Show the code</summary>

``` r
patchwork <- branch_length_by_rank / rwibald_key_plot / rwibald_key_plot_non_range_weighted + plot_annotation(
  #title = 'Branches of Interest for Range Weighted Branch Length Difference (RWiBaLD)\n&\nRWiBaLD Scores With Categories Calculated Using Elbow Statistic',
  theme = theme(plot.title = element_text(hjust = 0.5)),
  subtitle = '',
  caption = '',
tag_levels = 'A') + 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 12, hjust = 0, vjust = 0)) 

#print(patchwork)
#4062
ggsave("quarto_outputs/figures/Figure2_ABC.png", patchwork, width = 3000, height = 4062, units = "px")
```

</details>

<img src="quarto_outputs/figures/Figure2_ABC.png"
data-fig-align="center" />

## Generate Figure 3A

Here we create a plot of RWiBaLD score (x) by branch range (cells)
(Figure 3A), coloured by RWiBaLD Categories: neo-endemic (red),
meso-endemic (green) and paleo-endemic (blue).

<details class="code-fold">
<summary>Show the code</summary>

``` r
# Plot RWiBaLD rank X branch range
rwibald_key_plot_range <- ggplot() +
    geom_point(data=rwibald_results_all_with_range, aes(x = rwibald_score_rank, y = range_cell_count), fill="transparent", size = 2, colour= "grey77", alpha = 0.4, pch = 21) +
  geom_point(data = rwibald_results_all_key_only, aes(x = rwibald_score_rank, y = range_cell_count, fill = rwibald_type, color = rwibald_type), size = 2, alpha = 0.9, pch = 21) +
  scale_fill_manual(name = "RWiBaLD Category", values = c("paleo-endemic" = "royalblue1","neo-endemic" = "red", "meso-endemic" = "#FFD851"), labels = c("neo-endemic","meso-endemic", "paleo-endemic"), limits = c("neo-endemic","meso-endemic", "paleo-endemic"),  na.value = "transparent") +
  scale_color_manual(name = "RWiBaLD Category", values = c("paleo-endemic" = "royalblue4", "neo-endemic" = "red4", "meso-endemic" = "#D0A100"), labels = c("neo-endemic", "meso-endemic", "paleo-endemic"), limits = c("neo-endemic", "meso-endemic", "paleo-endemic"), na.value = "transparent") +
  geom_vline(xintercept = neo_cutoff, color = "red", linetype = "dashed") +
  geom_vline(xintercept = paleo_cutoff, color = "red", linetype = "dashed") +
  geom_vline(xintercept = results$Negative_Records, color = "grey74", linetype = "solid") +
  xlab("RWiBaLD Score Rank (shortest-longest)") +
  ylab("Branch Range (cells)") +
  theme_classic() +
  theme(plot.title = element_text(hjust =1),
        legend.position="none", 
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)), 
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)))

print(rwibald_key_plot_range)
```

</details>

![](RWiBaLD_calculations_and_figures_2025_files/figure-commonmark/unnamed-chunk-12-1.png)

<details class="code-fold">
<summary>Show the code</summary>

``` r
ggsave("quarto_outputs/figures/Figure3_A.png", plot = rwibald_key_plot_range,  scale = 1,
        width = 2048,
        height = 1024,
        units = "px",
        dpi = 300)
```

</details>

## Generate Figure 3B

Now the same data as figure 3A but first removing all internal branches
before plotting the data. Note: the neo / meso /paleo cuttofs are
calculated on the whole dataset and the filtering to terminals done
after that.

<details class="code-fold">
<summary>Show the code</summary>

``` r
# For plotting purposes get only the rwibald branches of key interest and terminals only
rwibald_results_all_key_only_terminals <- rwibald_results_all_with_range %>%
  filter(!grepl("^[0-9]", NAME)) %>%
  mutate(rwibald_score_rank_terminals = rank(rwibald_score, ties.method = "first")) %>% 
    filter(bl_comparison_tree_rwibald_elbow_key == "key") 

rwibald_results_all_with_range_terminals <- rwibald_results_all_with_range %>% 
    filter(!grepl("^[0-9]", NAME)) %>%
  mutate(rwibald_score_rank_terminals = rank(rwibald_score, ties.method = "first"))

# Find the neo_cutoff
neo_cutoff_terms <- rwibald_results_all_with_range_terminals %>%
  filter(rwibald_type == "neo-endemic") %>%
  summarize(highest_rank = max(rwibald_score_rank_terminals, na.rm = TRUE)) %>%
  pull(highest_rank)

# Find the paleo_cutoff
paleo_cutoff_terms <- rwibald_results_all_with_range_terminals %>%
  filter(rwibald_type == "paleo-endemic") %>%
  summarize(lowest_rank = min(rwibald_score_rank_terminals, na.rm = TRUE)) %>%
  pull(lowest_rank)

# Subset the left and right dataframes
subset_left_df_terms <- rwibald_results_all_with_range_terminals[rwibald_results_all_with_range_terminals$rwibald_score <= 0, ]

# Records in negative data
negative_records_terms <- nrow(subset_left_df_terms)

# View(rwibald_results_all_key_only_terminals)
# View(rwibald_results_all_with_range_terminals)

# View(rwibald_results_all_with_range)

# Plot rwibald rank x rwibald score
# rwibald_key_plot_terminals <- ggplot() +
#     geom_point(data=rwibald_results_all_with_range_terminals, aes(x = rwibald_score_rank_terminals, y = rwibald_score), fill="transparent", size =3, colour= "grey77", alpha = 0.4, pch=21) +
#   geom_point(data = rwibald_results_all_key_only_terminals, aes(x = rwibald_score_rank_terminals, y = rwibald_score, fill = rwibald_type, color = rwibald_type), size =3, alpha = 0.9, pch=21) +
#   scale_fill_manual(name = "RWiBaLD Category", values = c("paleo-endemic" = "royalblue1", "neo-endemic" = "red", "meso-endemic" = "#FFD851"), labels = c("neo-endemic","meso-endemic", "paleo-endemic"), limits = c("neo-endemic","meso-endemic", "paleo-endemic"),  na.value = "transparent") +
#   scale_color_manual(name = "RWiBaLD Category", values = c("paleo-endemic" = "royalblue4", "neo-endemic" = "red4", "meso-endemic" = "#D0A100"), labels = c("neo-endemic", "meso-endemic", "paleo-endemic"), limits = c("neo-endemic", "meso-endemic", "paleo-endemic"), na.value = "transparent") +
#   geom_vline(xintercept = neo_cutoff, color = "red", linetype = "dashed") +
#   geom_vline(xintercept = paleo_cutoff, color = "red", linetype = "dashed") +
#   xlab("RWiBaLD Score Rank Terminals Only (shortest-longest)") +
#   ylab("Range Weighted Branch\nLength Difference (RWiBaLD) Score") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5), legend.position="bottom", legend.text = element_text(size = 12))
# 
# print(rwibald_key_plot_terminals)
# 
# ggsave("quarto_outputs/figures/Figure3_B_terminals_only.png", plot = rwibald_key_plot_terminals,  scale = 1,
#         width = 2048,
#         height = 1024,
#         units = "px",
#         dpi = 300)

# Plot rwibald rank x branch range
rwibald_key_plot_range_terminals <- ggplot() +
    geom_point(data=rwibald_results_all_with_range_terminals, aes(x = rwibald_score_rank_terminals, y = range_cell_count), fill="transparent", size =3, colour= "grey77", alpha = 0.4, pch=21) +
  geom_point(data = rwibald_results_all_key_only_terminals, aes(x = rwibald_score_rank_terminals, y = range_cell_count, fill = rwibald_type, color = rwibald_type), size =3, alpha = 0.9, pch=21) +
  scale_fill_manual(name = "RWiBaLD Category", values = c("paleo-endemic" = "royalblue1","neo-endemic" = "red", "meso-endemic" = "#FFD851"), labels = c("neo-endemic","meso-endemic", "paleo-endemic"), limits = c("neo-endemic","meso-endemic", "paleo-endemic"),  na.value = "transparent") +
  scale_color_manual(name = "RWiBaLD Category", values = c("paleo-endemic" = "royalblue4", "neo-endemic" = "red4", "meso-endemic" = "#D0A100"), labels = c("neo-endemic", "meso-endemic", "paleo-endemic"), limits = c("neo-endemic", "meso-endemic", "paleo-endemic"), na.value = "transparent") +
  geom_vline(xintercept = neo_cutoff_terms, color = "red", linetype = "dashed") +
  geom_vline(xintercept = paleo_cutoff_terms, color = "red", linetype = "dashed") +
    geom_vline(xintercept = negative_records_terms, color = "grey74", linetype = "solid") +
  xlab("RWiBaLD Score Rank Terminals Only (shortest-longest)") +
  ylab("Branch Range (cells)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom", legend.text = element_text(size=rel(1)))

print(rwibald_key_plot_range_terminals)
```

</details>

![](RWiBaLD_calculations_and_figures_2025_files/figure-commonmark/unnamed-chunk-13-1.png)

<details class="code-fold">
<summary>Show the code</summary>

``` r
ggsave("quarto_outputs/figures/Figure3_B.png", plot = rwibald_key_plot_range_terminals,  scale = 1,
        width = 2048,
        height = 1024,
        units = "px",
        dpi = 300)
```

</details>

## Generate multi-panel figure 3AB

<details class="code-fold">
<summary>Show the code</summary>

``` r
patchwork <- rwibald_key_plot_range / rwibald_key_plot_range_terminals + plot_annotation(
  theme = theme(plot.title = element_text(hjust = 0.5)),
  subtitle = '',
  caption = '',
tag_levels = 'A') + 
   theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 12, hjust = 0, vjust = 0)) 

#print(patchwork)

ggsave("quarto_outputs/figures/Figure3_AB.png", patchwork, width = 3000, height = 3600, units = "px")
```

</details>

![](quarto_outputs/figures/Figure3_AB.png)

## Read tree file & generate supplementary figure 1

Now we load the range weighted tree file exported from biodiverse and
colour code it by the RWiBaLD categories.

<details class="code-fold">
<summary>Show the code</summary>

``` r
# if need be check and install ggtree using this code
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#  
# BiocManager::install("ggtree")

data_dir <- "Acacia_biodiverse_exports/"

all_tree_data_nwk <- paste0(data_dir, "acacia_tree_TRIMMED1_EQ1_RW1.nwk")

myTree <- read.tree(file=all_tree_data_nwk)

rwibald_results_all <- rwibald_results_all %>% 
  mutate(rwibald_type = replace(rwibald_type,is.na(rwibald_type),"other"))

# View(rwibald_results_all)

rwibald_results_all$tree_cols <- as.factor(rwibald_results_all$rwibald_type)

# Remove quotes from labels so they match dataframe
myTree$tip.label <- gsub("'","",myTree$tip.label)
myTree$node.label <- gsub("'","",myTree$node.label)

# Create ggtree object
g <- ggtree(myTree) %<+% rwibald_results_all +
      aes(color=tree_cols) +
      scale_color_manual(values = c("paleo-endemic" = "royalblue1",
                                 "neo-endemic" = "red",
                                 "meso-endemic" = "#FFD851")) +
     # geom_text2(aes(label = node), 
     #         hjust =0,    # adjust horizontal position
     #         vjust = 0,     # adjust vertical position
     #         size = 1,        # adjust text size
     #         color = "blue")  + # adjust text color
      theme(legend.position="none") +
       geom_tiplab(as_ylab=FALSE, size=1)

# Display the plot
#print(g)

ggsave("quarto_outputs/figures/SupFigure1_full_tree.png", g, width = 2048,
        height = 4024, units = "px")
```

</details>

## Generate tree in 2 halves supplementary figure 1A and 1B

<details class="code-fold">
<summary>Show the code</summary>

``` r
data_dir <- "Acacia_biodiverse_exports/"
all_tree_data_nwk <- paste0(data_dir, "acacia_tree_TRIMMED1_EQ1_RW1.nwk")

myTree <- read.tree(file=all_tree_data_nwk)

rwibald_results_all <- rwibald_results_all %>% 
  mutate(rwibald_type = replace(rwibald_type, is.na(rwibald_type), "other"))

rwibald_results_all$tree_cols <- as.factor(rwibald_results_all$rwibald_type)

# Remove quotes from labels so they match the dataframe
myTree$tip.label <- gsub("'", "", myTree$tip.label)
myTree$node.label <- gsub("'", "", myTree$node.label)

# Calculate the halfway point
half <- length(myTree$tip.label) / 2

# Sort the tip labels for consistency
tip_labels <- myTree$tip.label
first_half_tips <- tip_labels[1:half]
second_half_tips <- tip_labels[(half + 1):length(tip_labels)]

# Identify MRCA nodes for each half
node_first_half <- getMRCA(myTree, first_half_tips)
node_second_half <- getMRCA(myTree, second_half_tips)


# Create the base ggtree object
g <- ggtree(myTree) %<+% rwibald_results_all +
      aes(color = tree_cols) +
      scale_color_manual(values = c("paleo-endemic" = "royalblue1",
                                    "neo-endemic" = "red",
                                    "meso-endemic" = "#FFD851",
                                    "other" = "gray66")) +
      theme(legend.position = "none") +
      geom_tiplab(as_ylab = FALSE, size=1)

#plot(g)

# Identify MRCA node
mrca_node <- MRCA(g, second_half_tips)

second_half_tree <- viewClade(tree_view = g, mrca_node) +
    geom_point2(
      aes(subset = (node == mrca_node)),
    shape = 21,
    size = 2,
    fill = 'gray26',
    nudge_x = 20,
    position = position_nudge(y = -88)) 


# +
#   theme(panel.ontop = FALSE,
#          panel.spacing.x = unit(c(1, 1, 1, 1), "cm"),
#          panel.background = element_rect(color = "green", fill="transparent", linewidth = 10),
#          plot.margin = unit(c(1, 1, 1, 1), "cm"))


#plot(second_half_tree)
 
# Save the figure showing the second half expanded and the first half collapsed
ggsave("quarto_outputs/figures/SupFigure_1A.png", 
        second_half_tree, width = 2048, height = 4024, units = "px")

#plot(g)
# Collapse the second half of the tree
g_half_collapsed <- collapse(g, node = node_second_half) +
  # Add a point to indicate the collapsed node for reference
  geom_point2(aes(subset=(node == node_second_half)), shape=21, size=2, fill='gray26')

#plot(g_half_collapsed)

# Save the figure showing the first half expanded and the second half collapsed
ggsave("quarto_outputs/figures/SupFigure_1B.png", 
       g_half_collapsed, width = 2048, height = 4024, units = "px")
```

</details>

![](quarto_outputs/figures/SupFigure_1A.png)

![](quarto_outputs/figures/SupFigure_1B.png)

## Load data & generate histograms for all cells

Next we load the data and generate histograms. This code generates a
histogram for every cell in the dataset. For convenience the file names
contain the CANAPE code and whether the histogram contains branches of
key interest.

<details class="code-fold">
<summary>Show the code</summary>

``` r
#function to change decimal places displayed on plot
scaleFUN <- function(x) sprintf("%.8g", x)

# Function to replace '-' with 'm' and ':' with '_'
clean_column_names <- function(df) {
  new_names <- colnames(df) %>%
    gsub("-", "neg", .) %>%
    gsub(":", "_", .)
  colnames(df) <- new_names
  return(df)
}

# Clean column names using the function
canape_group_data_rwibald <- clean_column_names(canape_group_data_rwibald)

#View(canape_group_data_rwibald)

col_scheme <- c("paleo-endemic" = "royalblue1","neo-endemic" = "red", "meso-endemic" = "#FFD851", "other" = "lightgoldenrodyellow")
#legend_labels <- c("neo-endemic"="Neo Endemic","paleo-endemic"="paleo Endemic", "other" = "other", "meso-endemic" = "Meso Endemic")

#data <- transposed_df
#View(canape_key_group_data_rwibald)

#Generate all histograms
for (i in 2:(ncol(canape_group_data_rwibald)-2)) {# skip first column and last 2 as they are NAME,cell_count,rwibald_type
#for (i in 74:74) {
  #i <- 74
  data <- canape_group_data_rwibald[-(1:6),] # remove the rows not needed
  col_name <- colnames(data)[i] # get the column name to generate the histogram of
  colourByCol <- "rwibald_type"  # set the column to use for colouring the histogram
  numberOfBins <- 41 # set number of bins to use 
  hist_data <- na.omit(as.numeric(data[[col_name]])) # extract that column data for ploting
  hist_range <- range(hist_data,na.rm=TRUE) # calculate the range of all values in the column ie. min/max value

#calculate the total range of the whole data set to get a consistent x scale on all histograms.    
dataForRange <- data %>%
  select(-first(colnames(.)), -last(colnames(.))) %>%
  mutate(across(everything(), as.numeric)) 
totalRange <- range(dataForRange, na.rm=TRUE)#lowest and highedt values in the dataset
totalRangeGap <- max(dataForRange, na.rm = TRUE) - min(dataForRange, na.rm = TRUE) # distance between range above

  range <- max(hist_data) - min(hist_data)
  binSize <- range/numberOfBins
  farthest_number <- max(abs(hist_range)) # calculate the number farthest from 0
  x_min <- -(farthest_number)-2*(binSize) # set xmin so zero is centred
  x_max <- (farthest_number)+2*(binSize) # set xmax so zero is centred
 
  #check to see if any of the branches in the histogram are RWIBALD significant to add to the file names
  result <- data %>%
  select(!!sym(col_name), !!sym(colourByCol)) %>%  
  filter(!is.na(!!sym(col_name))) %>%
  summarise(result = ifelse(any(!!sym(colourByCol) %in% c("paleo-endemic", "neo-endemic", "meso-endemic")), "key", "zero-key")) %>%
  pull(result)
  
  #View(result)
  # setup the file names 
  CANAPE_CODE <- canape_group_data_rwibald["CANAPE_CODE", col_name]
  cell <- str_replace_all(col_name, "[:]", "_")
  filename <- paste0("quarto_outputs/histograms/", cell, "_CC_", CANAPE_CODE,"_", result, "_", numberOfBins,"bins.png")
  
  #print(filename)
  # get a subset of dataframe for the histogram
  plot_data <- data %>%
  select(!!sym(col_name), !!sym(colourByCol)) %>%  
  filter(!is.na(!!sym(col_name))) %>%
  arrange(!!sym(colourByCol))
  
  #View(plot_data)

 plot <- ggplot(plot_data, aes(x = as.numeric(.data[[col_name]]), fill = forcats::fct_rev(rwibald_type))) +
  geom_histogram(color = "black", linewidth= 0.2, bins = numberOfBins) +
  scale_fill_manual(values = col_scheme) +
  xlim(x_min, x_max) +
  coord_cartesian(ylim = c(0, 15)) +
  labs(x = "RWiBaLD Score", y = "Frequency", fill = "Branch Category") +
  theme_bw()
 
 plot <- plot +
  stat_bin(
    aes(label = after_stat(if_else (condition = count>15, as.character(count), ""))),
    bins= as.numeric(numberOfBins),
    position=position_stack(vjust=0.1),
    pad=TRUE,
    geom = "text",
    color = "black",
    size = 2,
    y = 14
  )
 
 CairoPNG(width = 1024, height = 600, file = filename, canvas="white", bg = "white", units="px", dpi=96, title = "") 
  print(plot)
 dev.off()
}
```

</details>

## Load data & generate 12 histograms for figure 4A-L

This code filters out only the 12 histograms included in Figure 4 and
generates the ggplot2 objects of them.

<details class="code-fold">
<summary>Show the code</summary>

``` r
data <- t(canape_group_data)
#View(data)

#Transpose the dataframe and set the first row as column names
transposed_df <- setNames(data.frame(data[-1,]), data[1,])

#View(transposed_df)

#cells to keep
cols_to_keep <- c("NAME", "-1825000:-2975000", "75000:-2525000", "625000:-3475000","-1425000:-3225000","1775000:-3725000", "-1425000:-3475000", "-1225000:-3775000", "-1425000:-3275000", "1225000:-1275000", "-1175000:-3275000", "-1075000:-2825000", "-1075000:-2375000")

#Calculate RWiBaLD significants and Statistics etc.
canape_group_data_rwibald <- transposed_df %>%
  tibble::rownames_to_column(var = "NAME") %>%
  select(all_of(cols_to_keep)) %>%
  left_join(select(rwibald_results_all, NAME, rwibald_type), by = 'NAME') %>%
  mutate(rwibald_type = ifelse(is.na(rwibald_type), "other", rwibald_type))   

# Check if the column is already a factor
 if (!is.factor(canape_group_data_rwibald$rwibald_type)) {
   # If not a factor, convert and specify levels
   canape_group_data_rwibald$rwibald_type <- factor(canape_group_data_rwibald$rwibald_type, levels = c("neo-endemic", "meso-endemic", "paleo-endemic", "other"))
 }

#View(canape_group_data_rwibald)

# Function to replace '-' with 'm' and ':' with '_'
clean_column_names <- function(df) {
  new_names <- colnames(df) %>%
    gsub("-", "neg", .) %>%
    gsub(":", "_", .)
  colnames(df) <- new_names
  return(df)
}

# Clean column names using the function
canape_group_data_rwibald <- clean_column_names(canape_group_data_rwibald)

rownames(canape_group_data_rwibald) <- canape_group_data_rwibald$NAME 

#View(canape_group_data_rwibald)

# Write data to file
write.csv(canape_group_data_rwibald, paste0("quarto_outputs/transposed_df_figure_only.csv"), row.names=FALSE)

#function to change decimal places displayed on plot
scaleFUN <- function(x) sprintf("%.8g", x)

#View(canape_group_data_rwibald)

col_scheme <- c("paleo-endemic" = "royalblue1",
                "neo-endemic" = "red",
                "meso-endemic" = "#FFD851",
                "other" = "lightgoldenrodyellow"
                )

legend_labels <- c("neo-endemic" = "neo-endemic",
                   "paleo-endemic" = "paleo-endemic",
                   "meso-endemic" = "meso-endemic",
                   "other" = "other"
                   )

legend_order <- c("neo-endemic", "meso-endemic", "paleo-endemic", "other")

#data <- transposed_df
#View(canape_group_data_rwibald)

# Create an empty list to store plots
plot_list <- list()

#Generate histograms
for (i in 2:(ncol(canape_group_data_rwibald)-1)) {# skip first column and last 2 as they are NAME,rwibald_type
  #i <- 2
  local({ # have to make this local to allow multiple plots in patchwork
    i <- i
  data <- canape_group_data_rwibald[-(1:6),] # remove the rows not needed
  col_name <- colnames(data)[i] # get the column name to generate the histogram of
  colourByCol <- "rwibald_type"  # set the column to use for colouring the histogram
  numberOfBins <- 41 # set number of bins to use 
  
  # have to show only one legend as the 'collect' feature in patchwork does not collect the legends as desired
  if(i == 3){
  showLegend <- TRUE
  } else {
  showLegend <- FALSE
  }
 
  hist_data <- na.omit(as.numeric(data[[col_name]])) # extract that column data for ploting
  hist_range <- range(hist_data,na.rm=TRUE) # calculate the range of all values in the column ie. min/max value

  #calculate the total range of the whole data set to get a consistent x scale on all histograms.    
  dataForRange <- data %>%
  select(-first(colnames(.)), -last(colnames(.))) %>%
  mutate(across(everything(), as.numeric)) 
  totalRange <- range(dataForRange, na.rm=TRUE)#lowest and highest values in the dataset
  totalRangeGap <- max(dataForRange, na.rm = TRUE) - min(dataForRange, na.rm = TRUE) # distance between range above

  range <- max(hist_data) - min(hist_data)
  binSize <- range/numberOfBins
  farthest_number <- max(abs(hist_range)) # calculate the number farthest from 0
  x_min <- -(farthest_number)-2*(binSize) # set xmin so zero is centred
  x_max <- (farthest_number)+2*(binSize) # set xmax so zero is centred
 
  # get a subset of dataframe for the histogram
  plot_data <- data %>%
  mutate(!!sym(col_name) := as.numeric(!!sym(col_name))) %>% 
  select(!!sym(col_name), !!sym(colourByCol)) %>%  
  filter(!is.na(!!sym(col_name))) %>%
  arrange(!!sym(colourByCol))
  
  #View(plot_data)
  #str(plot_data)
 
  plot <- ggplot(plot_data, aes(x = .data[[col_name]], fill = forcats::fct_rev(rwibald_type))) +
  geom_histogram(color = "black", linewidth = 0.2, bins = numberOfBins,  show.legend = showLegend) +
  scale_fill_manual(values = col_scheme, labels = legend_labels, breaks = legend_order, drop = FALSE, na.value = "transparent") +
  #xlim(x_min, x_max) +
  xlim(-0.011, 0.011) +  
  coord_cartesian(ylim = c(0, 13)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  labs(x = "RWiBaLD Score", y = "Frequency", fill = "RWiBaLD Category") +
  theme_bw() 

 plot <- plot +
  geom_rect(
    aes(xmin = -0.0005, xmax = 0.0005, ymin = 12.2, ymax = 12.8),
    color = "white",
    fill = "white"
  ) + 
   stat_bin(
    aes(label = after_stat(if_else (condition = count>12, as.character(count), ""))),
    bins= as.numeric(numberOfBins),
    position=position_stack(vjust=0.1),
    pad=TRUE,
    geom = "text",
    color = "black",
    size = 2.5,
    y = 12.5
  ) 

  #print(plot)
 
  # Dynamic variable name for the plot
  plot_name <- paste0("plot_", i - 1)  # Subtracting 1 to start numbering from 1
  assign(plot_name, plot)

  # Add the plot to the list
  plot_list[[i - 1]] <<- get(plot_name)
  #rm(plot)
  })
 
}

 #print(plot_list[[1]])

#save the plotlist for later on
saveRDS(plot_list, "quarto_outputs/plotlist.rds")
```

</details>

## Generate multi-panel histogram figure 4A-L

This code compiles all 12 histograms generated above into one figure.

<details class="code-fold">
<summary>Show the code</summary>

``` r
plot_list <- readRDS("quarto_outputs/plotlist.rds")

patchwork <-  plot_list[[1]] + xlab(NULL) +
              guides(x = "none") +
              plot_list[[2]] + xlab(NULL) + ylab(NULL) +
              guides(x = "none", y = "none") +
              plot_list[[3]] + xlab(NULL) + ylab(NULL) +
              guides(x = "none", y = "none") +
              plot_list[[4]] + xlab(NULL) + ylab(NULL) +
              guides(x = "none", y = "none") +
              plot_list[[5]] + xlab(NULL) + guides(x = "none") +
              plot_list[[6]] + xlab(NULL) + ylab(NULL) +
              guides(x = "none", y = "none") +
              plot_list[[7]] + xlab(NULL) + ylab(NULL) +
              guides(x = "none", y = "none") +
              plot_list[[8]] + xlab(NULL) + ylab(NULL) +
              guides(x = "none", y = "none") +
              plot_list[[9]] +
              plot_list[[10]] + ylab(NULL) + guides(y = "none") +
              plot_list[[11]] + ylab(NULL) + guides(y = "none") +
              plot_list[[12]] + ylab(NULL) + guides(y = "none") +
              plot_annotation(
              subtitle = '                     CANAPE-Neo                                CANAPE-paleo                                                       CANAPE-Mixed',
              caption = '',
              tag_levels = 'A') +
              plot_layout(ncol = 4, guides = 'collect') &
              theme(legend.position = "bottom", plot.tag.position = c(0.9, 0.9),
                    plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

#print(patchwork)

ggsave("quarto_outputs/figures/Figure4_A.png", patchwork, width = 3000, height = 2400, units = "px")
```

</details>

<img src="quarto_outputs/figures/figure4_A_inkscape.png"
data-fig-align="center" />

## Load results data & generate data files for biome maps

Now we take the RWiBaLD results and split the taxa in the original
specimen data from the [Categorical Analysis of Neo And Paleo-Endemisim
(CANAPE)](https://www.nature.com/articles/ncomms5473) paper according to
these RWiBaLD categories so we can overlay them on the biome map
published by [Crisp et al.](http://dx.doi.org/10.1098/rstb.2004.1528) in
2004.

<details class="code-fold">
<summary>Show the code</summary>

``` r
original_paper_specimen_data_file <- "doi_10_5061_dryad_dv4qk__v20150514/Point_distribution_Australian_Phylogenetic_Diversity_Acacia.csv"

original_paper_specimen_data <- read.table(original_paper_specimen_data_file, header=T, sep=",")

#View(original_paper_specimen_data)
#View(rwibald_results_all_with_range)

# Summarize the data in the rwibaled_type column
summary_df <- rwibald_results_all_with_range %>%
  group_by(rwibald_type) %>%
  summarise(count = n())  # Count of each rwibaled_type

# ---- Table output: gt for HTML/PDF, kable for GFM ----
if (knitr::is_html_output() || knitr::is_latex_output()) {

  summary_table <- summary_df |>
    gt::gt() |>
    gt::tab_header(
      title = "Summary of rwibald_type",
      subtitle = "Count of each type"
    )

  summary_table

} else {

  knitr::kable(
    summary_df,
    format = "pipe",
    align  = "c"
  )

}

# Print the table
#summary_table

# Extract the 'NAME' values where 'rwibald_type' is 'meso-endemic'
meso_names <- rwibald_results_all_with_range[rwibald_results_all_with_range$rwibald_type == "meso-endemic", "NAME"]

# Filter the 'original_paper_specimen_data' dataframe
meso_taxa <- original_paper_specimen_data[original_paper_specimen_data$Species %in% meso_names, ]

# Write data to file
write.csv(meso_taxa, "quarto_outputs/Meso_endemics.csv" ,row.names=FALSE)

# Extract the 'NAME' values where 'rwibald_type' is 'neo-endemic'
neo_names <- rwibald_results_all_with_range[rwibald_results_all_with_range$rwibald_type == "neo-endemic", "NAME"]

# Filter the 'original_paper_specimen_data' dataframe
neo_taxa <- original_paper_specimen_data[original_paper_specimen_data$Species %in% neo_names, ]

# Write data to file
write.csv(neo_taxa, "quarto_outputs/Neo_endemics.csv" ,row.names=FALSE)

# Extract the 'NAME' values where 'rwibald_type' is 'paleo-endemic'
paleo_names <- rwibald_results_all_with_range[rwibald_results_all_with_range$rwibald_type == "paleo-endemic", "NAME"]

# Filter the 'original_paper_specimen_data' dataframe
paleo_taxa <- original_paper_specimen_data[original_paper_specimen_data$Species %in% paleo_names, ]

# Write data to file
write.csv(paleo_taxa, "quarto_outputs/Paleo_endemics.csv" ,row.names=FALSE)
```

</details>

| rwibald_type  | count |
|:-------------:|:-----:|
| meso-endemic  |  59   |
|  neo-endemic  |  55   |
|     other     |  880  |
| paleo-endemic |  21   |

## Two tailed relative phylogenetic diversity (RPD) & CANAPE functions

This is the function for calculating Categorical Analysis of Neo- and
Paleo-endemism (CANAPE) from ranked P scores of Biodiverse analysis.

<details class="code-fold">
<summary>Show the code</summary>

``` r
#Standard 2 tailed test for RPD
significance_fun <- function(x){
  if (x >= 0.99) {
    return("Very Highly Sig")
  } else if (x >= 0.975){
    return ("Highly Sig")
  } else if (x <= 0.01){
    return ("Very Sig Low")
  } else if (x <= 0.025){
    return ("Sig Low")
  } else {
    return("Not Sig")
  }
}

#two pass test for RPE
# x=P_PE_WE_P, y=P_PHYLO_RPE_NULL2, z=P_PHYLO_RPE2
significance_super_fun <- function(x, y, z){
  if (x > 0.95 || y > 0.95) {
    if (z <= 0.025){
      return ("Neo")
    } else if (z >= 0.975){
      return ("paleo")
    } else if (x >= 0.99 || y >= 0.99){
      return ("Super")
    } else {
      return("Mixed")
    }
  } else {
    return("Not Sig")
  }
}
```

</details>

## Load CANAPE & randomisation results

Here we load the re-created biodiverse results from the original paper,
calculate CANAPE and re-generate the CANAPE map. Note: these results may
differ slightly from the original paper as we re-ran the randomisation
in biodiverse.

<details class="code-fold">
<summary>Show the code</summary>

``` r
data_dir <- "Acacia_biodiverse_exports/"

# The CANAPE results file calculated in Biodiverse directly and exported
biodiverse_canape_results_file <- paste0(data_dir, "Acacia_Rand1_CANAPE_Export.csv")

# The spatial & randomisation results files calculated in Biodiverse and exported
biodiverse_observed_data_file  <- paste0(data_dir, "Acacia_SPATIAL_RESULTS_Export.csv")
biodiverse_rand_results_file <- paste0(data_dir, "Acacia_Rand1_SPATIAL_RESULTS_Export.csv")

biodiverse_canape_results <- read.table(biodiverse_canape_results_file, header=T,sep=",", check.names = FALSE )

biodiverse_observed_spatial_results <- read.table(biodiverse_observed_data_file, header=T,sep=",")
biodiverse_rand_spatial_results <- read.table(biodiverse_rand_results_file, header=T,sep=",")

biodiverse_results_concatenated <- cbind(biodiverse_observed_spatial_results, biodiverse_rand_spatial_results)

#View(biodiverse_results_concatenated)

###############################################
#Create new columns in dataframe and 
#populate them using the functions above
###############################################

targets <- c("PHYLO_RPD2", "PD_P", "PE_WE_P", "PD_P_per_taxon", "PHYLO_RPE2")

for (name in targets) {
  colname <- paste0("P_", name)  #  prepend the P_ since we want the proportions, saves some typing above
  new_colname = paste0(colname, "_SIG")
  trait_index <- match (colname, colnames(biodiverse_results_concatenated))
  # Apply the function to every row of column with index "trait_index"
  #  and generate a column in the dataframe showing significant cells
  if (!is.na(trait_index)) {
    biodiverse_results_concatenated[[new_colname]] <- apply (biodiverse_results_concatenated[trait_index],  MARGIN=c(1), significance_fun)
  } else {
    print (paste("Cannot find index", colname, "in data frame"))
  }
}

biodiverse_results_concatenated$P_PHYLO_RPE2_CANAPE_SIG <- sapply(
  1:nrow(biodiverse_results_concatenated),
  function(x) significance_super_fun(
    biodiverse_results_concatenated$P_PE_WE_P[x],
    biodiverse_results_concatenated$P_PHYLO_RPE_NULL2[x],
    biodiverse_results_concatenated$P_PHYLO_RPE2[x]
  )
)

#View(biodiverse_results_concatenated)
```

</details>

## Generate CANAPE map figure 5A

Here we generate the CANAPE map, indicating cell locations from the
histograms in figure 4A-L. CANAPE identifies geographic concentrations
of high PE, and gives a summary classification of the type of endemism
dominating in a location. The RWiBaLD histograms (figure 4A-L) identify
the specific branches that contribute the most to PE in a given cell,
and what type of endemism they represent.

<details class="code-fold">
<summary>Show the code</summary>

``` r
#loadfonts()

myFont <- choose_font(c("HelvLight", "Arial", "sans"), quiet = TRUE) #load a font if available
   
map_text <- "Categorical Analysis of Neo- And Paleo- Endemism"
sigplot <- "P_PHYLO_RPE2_CANAPE_SIG"
col_scheme <- c("paleo" = "royalblue1","Not Sig" = "snow2", "Neo" = "red", "Super" = "#9D00FF", "Mixed"= "#CB7FFF")
legend_order <-c("Neo","paleo", "Not Sig", "Mixed", "Super")
legend_labels <- c("Neo"="Neo","paleo"="Paleo", "Not Sig"="Not significant", "Mixed"="Mixed", "Super"="Super")

biodiverse_results_concatenated[, sigplot] <- factor(biodiverse_results_concatenated[, sigplot], levels=legend_order)
Axis_0 <- "Axis_0"
Axis_1 <- "Axis_1"   
   
map_shape_file <- paste0("shape_files/coastline_albers.shp")
map_data   <- st_read(map_shape_file)  

cols_to_keep <- c("-1825000:-2975000", "75000:-2525000", "625000:-3475000","-1425000:-3225000","1775000:-3725000", "-1425000:-3475000", "-1225000:-3775000", "-1425000:-3275000", "1225000:-1275000", "-1175000:-3275000", "-1075000:-2825000", "-1075000:-2375000")

# Initialize empty vectors for x and y
x <- c()
y <- c()

# Split each element of cols_to_keep and append to x and y for plotting on the map
for (val in cols_to_keep) {
  parts <- strsplit(val, ":")[[1]]
  x <- c(x, as.numeric(parts[1]))
  y <- c(y, as.numeric(parts[2]))
}

# Print the result
#print(x)
#print(y)

labels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")

# Create data frame
histograms <- data.frame(x = x, y = y, label = labels)

# Plot the sf object with ggplot2, coloring by the BIOME column
map_plot_CANAPE  <- ggplot() +
  geom_tile(data=biodiverse_results_concatenated, aes_string(x=Axis_0, y=Axis_1, fill=sigplot)) +
  geom_sf(data = map_data, colour = "grey77" , fill="transparent") +
  scale_fill_manual(values = col_scheme,  labels=legend_labels, name="CANAPE", guide = guide_legend(direction = "horizontal", title.position = "bottom", title.hjust=0.5, title.vjust=0.5, label.position="bottom", label.hjust = 0.5, label.vjust = 0.1, lineheight=0.5))+
  geom_text_repel(data=histograms, aes(x = x, y = y, label = labels), fontface = "bold", size = 2.5, nudge_x = c(-200000, 0, 0, -100000, 250000, -300000, 200000, 100000, 200000, 200000, 200000, 200000), nudge_y = c(200000, 250000, 350000, 300000, -100000, 0, -200000, 200000, 200000, 200000, 200000, 200000)) +
   annotate("rect", xmin = -750000, xmax = -250000, ymin = -4500000, ymax = -4550000, fill = "black", colour = "black", alpha = 1, linewidth = 0.1) +
    annotate("rect", xmin = -250000, xmax = 250000, ymin = -4500000, ymax = -4550000, fill = "white", colour = "black", alpha = 1, linewidth = 0.1) +
    annotate("text", label = "0", x = -750000, y = -4650000, size=rel(2),  face = 'plain', family = myFont) +
    annotate("text", label = "500", x = -250000, y = -4650000, size=rel(2),  face = 'plain', family = myFont) +
    annotate("text", label = "1000", x = 250000, y = -4650000, size=rel(2),  face = 'plain', family = myFont) +
    annotate("text", label = "km", x = 500000, y = -4650000, size=rel(2),  face = 'plain', family = myFont) +
   theme(text = element_text(family = myFont),
         strip.background = element_blank(),
         title = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont),
         axis.line=element_blank(),axis.text.x=element_blank(),
         axis.text.y=element_blank(),axis.ticks=element_blank(),
         axis.title.x=element_blank(), axis.title.y=element_blank(),
         legend.position="none",
         legend.direction='horizontal',
         legend.text = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont),
         panel.grid = element_blank(),
         panel.background=element_blank(),#element_rect(colour = "black", fill="white", size = 1),
         panel.border = element_blank(),
         plot.background=element_blank(),#element_rect(colour = "black", fill="white", size = 1),
         plot.margin=unit(c(0,0,0,0),"line"))

print(map_plot_CANAPE)
```

</details>

![](RWiBaLD_calculations_and_figures_2025_files/figure-commonmark/unnamed-chunk-23-1.png)

<details class="code-fold">
<summary>Show the code</summary>

``` r
ggsave('quarto_outputs/figures/Figure5_A.png', map_plot_CANAPE, width = 1024, height = 1024, units = "px", bg = "white")
```

</details>

    Reading layer `coastline_albers' from data source 
      `C:\E\gitRepos\rwibald\shape_files\coastline_albers.shp' using driver `ESRI Shapefile'
    Simple feature collection with 44 features and 4 fields
    Geometry type: POLYGON
    Dimension:     XY
    Bounding box:  xmin: -1887741 ymin: -4840771 xmax: 2121462 ymax: -1029671
    Projected CRS: GDA94 / Australian Albers

## Generate biome map with neo-endemics figure 5B

Here we show the specimen data for all of the neo-endemic taxa
categorised by RWiBaLD overlaid on the biome map published by [Crisp et
al.](http://dx.doi.org/10.1098/rstb.2004.1528) in 2004.

<details class="code-fold">
<summary>Show the code</summary>

``` r
#loadfonts()

#myFont <- choose_font(c("HelvLight", "Arial", "sans"), quiet = TRUE) #load a font if available

col_scheme <- c("aseasonal wet" = "cornflowerblue",
                "eremean" = "lightyellow1",
                "monsoonal tropics" = "powderblue",
                "southeastern temperate" = "moccasin",
                "southwestern temperate"= "#fbb4ae")
legend_order <-c("aseasonal wet", 
                 "eremean", 
                 "monsoonal tropics",
                 "southeastern temperate", 
                 "southwestern temperate")
legend_labels <- c("aseasonal wet" = "aseasonal wet",
                   "eremean" = "eremean",
                   "monsoonal tropics" = "monsoonal tropics",
                   "southeastern temperate"="southeastern temperate",
                   "southwestern temperate"="southwestern temperate")

map_shape_file1 <- paste0("shape_files/coastline_albers.shp")
map_data   <- st_read(map_shape_file1)   
#geom_sf(data = map_data, colour = "grey77" , fill="transparent") +

map_shape_file <- paste0("shape_files/biomes_crisp_3577/biomes_crisp_3577.shp")
# Read the shapefile as an sf object
biomes_sf <- st_read(map_shape_file)

csv_file <- paste0("quarto_outputs/Neo_endemics.csv") 
# Read the CSV file into a data frame
csv_data <- read.csv(csv_file)
# Convert data frame to sf object
csv_sf <- st_as_sf(csv_data, coords = c("x_metres_EPSG_3577_Albers_Equal_Area", "y_metres_EPSG_3577_Albers_Equal_Area"), crs = st_crs(biomes_sf))


# Plot the sf object with ggplot2, coloring by the BIOME column
map_plot_Neo  <- ggplot() +
  geom_sf(data = biomes_sf, aes(fill = BIOME), colour = "grey77") +
  geom_sf(data = map_data, colour = "grey77" , fill="transparent") +
  geom_sf(data = csv_sf, fill= "red", color = "red4", shape = 21) +
 scale_fill_manual(values = col_scheme,  labels=legend_labels, name="Biomes", guide = guide_legend(direction = "horizontal", title.position = "bottom", title.hjust=0.5, title.vjust=0.5, label.position="bottom", label.hjust = 0.5, label.vjust = 0.1, lineheight=0.5))+
   theme(text = element_text(family = myFont),
         strip.background = element_blank(),
         title = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont),
         axis.line=element_blank(),axis.text.x=element_blank(),
         axis.text.y=element_blank(),axis.ticks=element_blank(),
         axis.title.x=element_blank(), axis.title.y=element_blank(),
         legend.position="none",
         legend.direction='horizontal',
         legend.text = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont),
         legend.key = element_rect(color = "black", size = 0),
         panel.grid = element_blank(),
         panel.background=element_blank(),#element_rect(colour = "black", fill="white", size = 1),
         panel.border = element_blank(),
         plot.background=element_blank(),#element_rect(colour = "black", fill="white", size = 1),
         plot.margin=unit(c(0,0,0,0),"line"))


print(map_plot_Neo)
```

</details>

![](RWiBaLD_calculations_and_figures_2025_files/figure-commonmark/unnamed-chunk-24-1.png)

<details class="code-fold">
<summary>Show the code</summary>

``` r
ggsave("quarto_outputs/figures/Figure5_B.png", map_plot_Neo, width = 3000, height = 2600, units = "px")
```

</details>

    Reading layer `coastline_albers' from data source 
      `C:\E\gitRepos\rwibald\shape_files\coastline_albers.shp' using driver `ESRI Shapefile'
    Simple feature collection with 44 features and 4 fields
    Geometry type: POLYGON
    Dimension:     XY
    Bounding box:  xmin: -1887741 ymin: -4840771 xmax: 2121462 ymax: -1029671
    Projected CRS: GDA94 / Australian Albers
    Reading layer `biomes_crisp_3577' from data source 
      `C:\E\gitRepos\rwibald\shape_files\biomes_crisp_3577\biomes_crisp_3577.shp' 
      using driver `ESRI Shapefile'
    Simple feature collection with 5 features and 3 fields
    Geometry type: MULTIPOLYGON
    Dimension:     XY
    Bounding box:  xmin: -1887414 ymin: -4840588 xmax: 2121442 ymax: -1087816
    Projected CRS: GDA94 / Australian Albers

## Generate biome map with meso-endemics figure 5C

Here we show the specimen data for all of the meso-endemic taxa
categorised by RWiBaLD overlaid on the biome map published by [Crisp et
al.](http://dx.doi.org/10.1098/rstb.2004.1528) in 2004.

<details class="code-fold">
<summary>Show the code</summary>

``` r
#loadfonts()

#myFont <- choose_font(c("HelvLight", "Arial", "sans"), quiet = TRUE) #load a font if available

map_text <- "Categorical Analysis of Neo- And Paleo- Endemism"

col_scheme <- c("aseasonal wet" = "cornflowerblue",
                "eremean" = "lightyellow1",
                "monsoonal tropics" = "powderblue",
                "southeastern temperate" = "moccasin",
                "southwestern temperate"= "#fbb4ae")
legend_order <-c("aseasonal wet", 
                 "eremean", 
                 "monsoonal tropics",
                 "southeastern temperate", 
                 "southwestern temperate")
legend_labels <- c("aseasonal wet" = "aseasonal wet",
                   "eremean" = "eremean",
                   "monsoonal tropics" = "monsoonal tropics",
                   "southeastern temperate"="southeastern temperate",
                   "southwestern temperate"="southwestern temperate")

map_shape_file1 <- paste0("shape_files/coastline_albers.shp")
map_data   <- st_read(map_shape_file1)   
#geom_sf(data = map_data, colour = "grey77" , fill="transparent") +

map_shape_file <- paste0("shape_files/biomes_crisp_3577/biomes_crisp_3577.shp")
# Read the shapefile as an sf object
biomes_sf <- st_read(map_shape_file)

csv_file <- paste0("quarto_outputs/Meso_endemics.csv") 
# Read the CSV file into a data frame
csv_data <- read.csv(csv_file)
# Convert data frame to sf object
csv_sf <- st_as_sf(csv_data, coords = c("x_metres_EPSG_3577_Albers_Equal_Area", "y_metres_EPSG_3577_Albers_Equal_Area"), crs = st_crs(biomes_sf))


# Plot the sf object with ggplot2, coloring by the BIOME column
map_plot_Meso  <- ggplot() +
  geom_sf(data = biomes_sf, aes(fill = BIOME), color = "grey77") +
  geom_sf(data = map_data, colour = "grey77" , fill="transparent") +
  geom_sf(data = csv_sf, fill= "#FFD851", color = "#D0A100", shape = 21) +
  scale_fill_manual(values = col_scheme,  labels=legend_labels, name="Biomes", guide = guide_legend(direction = "horizontal", title.position = "bottom", title.hjust=0.5, title.vjust=0.5, label.position="bottom", label.hjust = 0.5, label.vjust = 0.1, lineheight=0.5))+
   theme(text = element_text(family = myFont),
         strip.background = element_blank(),
         title = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont),
         axis.line=element_blank(),axis.text.x=element_blank(),
         axis.text.y=element_blank(),axis.ticks=element_blank(),
         axis.title.x=element_blank(), axis.title.y=element_blank(),
         legend.position="none",
         legend.direction='horizontal',
         legend.text = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont),
         legend.key = element_rect(color = "black", size = 0),
         panel.grid = element_blank(),
         panel.background=element_blank(),#element_rect(colour = "black", fill="white", size = 1),
         panel.border = element_blank(),
         plot.background=element_blank(),#element_rect(colour = "black", fill="white", size = 1),
         plot.margin=unit(c(0,0,0,0),"line"))

print(map_plot_Meso)
```

</details>

![](RWiBaLD_calculations_and_figures_2025_files/figure-commonmark/unnamed-chunk-25-1.png)

<details class="code-fold">
<summary>Show the code</summary>

``` r
#cvdPlot(map_plot_Meso)

ggsave("quarto_outputs/figures/Figure5_C.png", map_plot_Meso, width = 3000, height = 2600, units = "px")
```

</details>

    Reading layer `coastline_albers' from data source 
      `C:\E\gitRepos\rwibald\shape_files\coastline_albers.shp' using driver `ESRI Shapefile'
    Simple feature collection with 44 features and 4 fields
    Geometry type: POLYGON
    Dimension:     XY
    Bounding box:  xmin: -1887741 ymin: -4840771 xmax: 2121462 ymax: -1029671
    Projected CRS: GDA94 / Australian Albers
    Reading layer `biomes_crisp_3577' from data source 
      `C:\E\gitRepos\rwibald\shape_files\biomes_crisp_3577\biomes_crisp_3577.shp' 
      using driver `ESRI Shapefile'
    Simple feature collection with 5 features and 3 fields
    Geometry type: MULTIPOLYGON
    Dimension:     XY
    Bounding box:  xmin: -1887414 ymin: -4840588 xmax: 2121442 ymax: -1087816
    Projected CRS: GDA94 / Australian Albers

## Generate biome map with paleo-endemics figure 5D

Here we show the specimen data for all of the paleo-endemic taxa
categorised by RWiBaLD overlaid on the biome map published by [Crisp et
al.](http://dx.doi.org/10.1098/rstb.2004.1528) in 2004.

<details class="code-fold">
<summary>Show the code</summary>

``` r
#loadfonts()

#myFont <- choose_font(c("HelvLight", "Arial", "sans"), quiet = TRUE) #load a font if available

map_text <- "Categorical Analysis of Neo- And Paleo- Endemism"

col_scheme <- c("aseasonal wet" = "cornflowerblue",
                "eremean" = "lightyellow1",
                "monsoonal tropics" = "powderblue",
                "southeastern temperate" = "moccasin",
                "southwestern temperate"= "#fbb4ae")
legend_order <-c("aseasonal wet", 
                 "eremean", 
                 "monsoonal tropics",
                 "southeastern temperate", 
                 "southwestern temperate")
legend_labels <- c("aseasonal wet" = "aseasonal wet",
                   "eremean" = "eremean",
                   "monsoonal tropics" = "monsoonal tropics",
                   "southeastern temperate"="southeastern temperate",
                   "southwestern temperate"="southwestern temperate")

map_shape_file1 <- paste0("shape_files/coastline_albers.shp")
map_data  <- st_read(map_shape_file1)   
#geom_sf(data = map_data, colour = "grey77" , fill="transparent") +

map_shape_file <- paste0("shape_files/biomes_crisp_3577/biomes_crisp_3577.shp")

# Read the shapefile as an sf object
biomes_sf <- st_read(map_shape_file)

csv_file <- paste0("quarto_outputs/Paleo_endemics.csv") 

# Read the CSV file into a data frame
csv_data <- read.csv(csv_file)

# Convert data frame to sf object
csv_sf <- st_as_sf(csv_data, coords = c("x_metres_EPSG_3577_Albers_Equal_Area", "y_metres_EPSG_3577_Albers_Equal_Area"), crs = st_crs(biomes_sf))

# Plot the sf object with ggplot2, coloring by the BIOME column
map_plot_Paleo  <- ggplot() +
  geom_sf(data = biomes_sf, aes(fill = BIOME), color = "grey77") +
  geom_sf(data = map_data, colour = "grey77" , fill="transparent") +
  geom_sf(data = csv_sf, fill= "royalblue3", color = "darkblue", shape = 21) +
  scale_fill_manual(values = col_scheme,  labels=legend_labels, name="Biomes", guide = guide_legend(direction = "horizontal", title.position = "bottom", title.hjust=0.5, title.vjust=0.5, label.position="bottom", label.hjust = 0.5, label.vjust = 0.1, lineheight=0.5))+
    theme(text = element_text(family = myFont),
         strip.background = element_blank(),
         title = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont),
         axis.line=element_blank(),axis.text.x=element_blank(),
         axis.text.y=element_blank(),axis.ticks=element_blank(),
         axis.title.x=element_blank(), axis.title.y=element_blank(),
         legend.position="none",
         legend.direction='horizontal',
         legend.text = element_text(colour = 'black', angle = 0, size=rel(0.5), face = 'plain', family = myFont),
         legend.key = element_rect(color = "black", size = 0),
         panel.grid = element_blank(),
         panel.background=element_blank(),#element_rect(colour = "black", fill="white", size = 1),
         panel.border = element_blank(),
         plot.background=element_blank(),#element_rect(colour = "black", fill="white", size = 1),
         plot.margin=unit(c(0,0,0,0),"line"))
  

print(map_plot_Paleo)
```

</details>

![](RWiBaLD_calculations_and_figures_2025_files/figure-commonmark/unnamed-chunk-26-1.png)

<details class="code-fold">
<summary>Show the code</summary>

``` r
ggsave("quarto_outputs/figures/Figure5_D.png", map_plot_Paleo, width = 3000, height = 2600, units = "px")
```

</details>

    Reading layer `coastline_albers' from data source 
      `C:\E\gitRepos\rwibald\shape_files\coastline_albers.shp' using driver `ESRI Shapefile'
    Simple feature collection with 44 features and 4 fields
    Geometry type: POLYGON
    Dimension:     XY
    Bounding box:  xmin: -1887741 ymin: -4840771 xmax: 2121462 ymax: -1029671
    Projected CRS: GDA94 / Australian Albers
    Reading layer `biomes_crisp_3577' from data source 
      `C:\E\gitRepos\rwibald\shape_files\biomes_crisp_3577\biomes_crisp_3577.shp' 
      using driver `ESRI Shapefile'
    Simple feature collection with 5 features and 3 fields
    Geometry type: MULTIPOLYGON
    Dimension:     XY
    Bounding box:  xmin: -1887414 ymin: -4840588 xmax: 2121442 ymax: -1087816
    Projected CRS: GDA94 / Australian Albers

## Generate 4 up map figure 5ABCD

Here we compile the maps into one figure 5A-D

<details class="code-fold">
<summary>Show the code</summary>

``` r
patchwork_map <- map_plot_CANAPE +
              map_plot_Neo + 
              map_plot_Meso + 
              map_plot_Paleo

patchwork_map <- patchwork_map +
              plot_annotation(subtitle = '',
                              caption = '',
                              tag_levels = 'A') +
              plot_layout(ncol = 2, guides = 'collect') &
              theme(legend.position = "bottom", 
                    legend.title=element_text(size=rel(0.6)),
                    legend.text=element_text(size=rel(0.5)),
                    legend.title.position = "bottom",
                    legend.label.position= "bottom",
                    legend.key = element_rect(color = "black", size = 0.5),
                    legend.key.height = unit(0.5, "cm"),
                    legend.key.width = unit(0.5, "cm"),
                    plot.tag.position = c(0.9, 0.9),
                    plot.tag = element_text(size = 12, hjust = 0, vjust = 0))
             
#print(patchwork_map)

#cvdPlot(patchwork_map)

ggsave('quarto_outputs/figures/Figure5_ABCD.png', patchwork_map, width = 3000, height = 2600, units = "px")
```

</details>

<img src="quarto_outputs/figures/Figure5_ABCD.png"
data-fig-align="center" />

## Table of RWiBaLD results supplementary table 1

RWiBaLD results for each branch on the phylogeny of Australian Acacia.

<details class="code-fold">
<summary>Show the code</summary>

``` r
RWiBaLD_results_csv  <- paste0("quarto_outputs/Acacia_RWiBaLD_results_all_with_range.csv")

rwibald_results_all_with_range <- read.table(RWiBaLD_results_csv, header=T,sep=",")

#View(rwibald_results_all_with_range)

#colnames(rwibald_results_all_with_range)

gt_table <- rwibald_results_all_with_range %>%
            select("NAME", "branch_length_comparison_tree", "branch_length_observed_tree",
                   "rwibald_score", "rwibald_type", "range_cell_count")

# ---- Table output: gt for HTML/PDF, kable for GFM ----
if (knitr::is_html_output() || knitr::is_latex_output()) {

  rwibald_results_all_with_range_table <- gt::gt(gt_table) |>
    gt::tab_options(
      table.width = gt::pct(100),
      table.layout = "auto",
      table.align = "left",
      table.margin.left = gt::px(5),
      table.margin.right = gt::px(5),
      table.font.size = gt::px(8),
      column_labels.font.size = gt::px(10),
      heading.align = "center",
      heading.title.font.size = gt::px(12),
      quarto.use_bootstrap = TRUE
    ) |>
    gt::opt_row_striping() |>
    gt::tab_header(title = gt::md("RWiBaLD Results")) |>
    gt::tab_style(
      style = gt::cell_text(align = "center"),
      locations = gt::cells_body()
    )

  rwibald_results_all_with_range_table

} else {

  # GFM-safe fallback (pipe table)
  knitr::kable(
    gt_table,
    format = "pipe",
    align  = "c"
  )

}


# print table using gt
# rwibald_results_all_with_range_table <- gt(gt_table) %>%
#   tab_options(
#     table.width = pct(100),
#     table.layout = "auto",
#     table.align = "left",
#     table.margin.left = px(5),
#     table.margin.right = px(5),
#     table.font.size = px(8),
#     column_labels.font.size = px(10),
#     heading.align = "center",
#     heading.title.font.size = px(12),
#     quarto.use_bootstrap = TRUE
#   ) %>%
#   opt_row_striping() %>%
#   tab_header(
#     title = md("RWiBaLD Results")
#   ) %>%
#   tab_style(
#     style = cell_text(align = "center"),
#     locations = cells_body()
#   )

#rwibald_results_all_with_range_table
```

</details>

| NAME | branch_length_comparison_tree | branch_length_observed_tree | rwibald_score | rwibald_type | range_cell_count |
|:--:|:--:|:--:|:--:|:--:|:--:|
| 1\_\_\_ | 0.0000878 | 0.0002235 | 0.0001357 | other | 103 |
| 10\_\_\_ | 0.0000895 | 0.0000358 | -0.0000537 | other | 101 |
| 100\_\_\_ | 0.0000411 | 0.0000072 | -0.0000339 | other | 220 |
| 101\_\_\_ | 0.0000237 | 0.0000336 | 0.0000099 | other | 382 |
| 102\_\_\_ | 0.0000232 | 0.0000188 | -0.0000043 | other | 390 |
| 103\_\_\_ | 0.0000214 | 0.0000470 | 0.0000256 | other | 423 |
| 104\_\_\_ | 0.0002318 | 0.0009137 | 0.0006819 | other | 39 |
| 105\_\_\_ | 0.0006028 | 0.0016249 | 0.0010222 | other | 15 |
| 106\_\_\_ | 0.0001532 | 0.0001057 | -0.0000475 | other | 59 |
| 107\_\_\_ | 0.0000994 | 0.0003546 | 0.0002552 | other | 91 |
| 108\_\_\_ | 0.0000983 | 0.0000343 | -0.0000640 | other | 92 |
| 109\_\_\_ | 0.0000845 | 0.0001130 | 0.0000285 | other | 107 |
| 11\_\_\_ | 0.0000878 | 0.0000290 | -0.0000588 | other | 103 |
| 110\_\_\_ | 0.0000189 | 0.0000090 | -0.0000099 | other | 478 |
| 111\_\_\_ | 0.0007534 | 0.0013431 | 0.0005896 | other | 12 |
| 112\_\_\_ | 0.0001615 | 0.0000434 | -0.0001180 | other | 56 |
| 113\_\_\_ | 0.0001051 | 0.0001208 | 0.0000157 | other | 86 |
| 114\_\_\_ | 0.0000878 | 0.0000000 | -0.0000878 | other | 103 |
| 115\_\_\_ | 0.0000786 | 0.0000099 | -0.0000687 | other | 115 |
| 116\_\_\_ | 0.0000718 | 0.0000135 | -0.0000583 | other | 126 |
| 117\_\_\_ | 0.0000701 | 0.0000408 | -0.0000293 | other | 129 |
| 118\_\_\_ | 0.0000701 | 0.0000562 | -0.0000139 | other | 129 |
| 119\_\_\_ | 0.0000695 | 0.0000178 | -0.0000518 | other | 130 |
| 12\_\_\_ | 0.0000134 | 0.0000110 | -0.0000024 | other | 677 |
| 120\_\_\_ | 0.0000160 | 0.0000078 | -0.0000081 | other | 566 |
| 121\_\_\_ | 0.0000150 | 0.0000147 | -0.0000004 | other | 601 |
| 122\_\_\_ | 0.0000099 | 0.0000041 | -0.0000058 | other | 915 |
| 123\_\_\_ | 0.0000099 | 0.0000014 | -0.0000084 | other | 917 |
| 124\_\_\_ | 0.0003118 | 0.0003328 | 0.0000211 | other | 29 |
| 125\_\_\_ | 0.0002511 | 0.0003888 | 0.0001377 | other | 36 |
| 126\_\_\_ | 0.0002055 | 0.0000457 | -0.0001597 | other | 44 |
| 127\_\_\_ | 0.0002825 | 0.0000691 | -0.0002134 | other | 32 |
| 128\_\_\_ | 0.0002825 | 0.0002412 | -0.0000413 | other | 32 |
| 129\_\_\_ | 0.0005023 | 0.0002914 | -0.0002109 | other | 18 |
| 13\_\_\_ | 0.0000131 | 0.0000056 | -0.0000075 | other | 690 |
| 130\_\_\_ | 0.0001130 | 0.0000390 | -0.0000740 | other | 80 |
| 131\_\_\_ | 0.0001005 | 0.0000416 | -0.0000589 | other | 90 |
| 132\_\_\_ | 0.0000786 | 0.0001279 | 0.0000492 | other | 115 |
| 133\_\_\_ | 0.0000766 | 0.0000395 | -0.0000371 | other | 118 |
| 134\_\_\_ | 0.0000386 | 0.0000018 | -0.0000368 | other | 234 |
| 135\_\_\_ | 0.0000435 | 0.0000426 | -0.0000008 | other | 208 |
| 136\_\_\_ | 0.0000655 | 0.0000444 | -0.0000211 | other | 138 |
| 137\_\_\_ | 0.0000281 | 0.0000524 | 0.0000243 | other | 322 |
| 138\_\_\_ | 0.0000254 | 0.0000742 | 0.0000488 | other | 356 |
| 139\_\_\_ | 0.0004759 | 0.0000652 | -0.0004106 | other | 19 |
| 14\_\_\_ | 0.0000115 | 0.0000044 | -0.0000071 | other | 789 |
| 140\_\_\_ | 0.0003767 | 0.0002540 | -0.0001227 | other | 24 |
| 141\_\_\_ | 0.0003477 | 0.0005309 | 0.0001831 | other | 26 |
| 142\_\_\_ | 0.0002917 | 0.0004352 | 0.0001436 | other | 31 |
| 143\_\_\_ | 0.0002055 | 0.0000266 | -0.0001789 | other | 44 |
| 144\_\_\_ | 0.0001507 | 0.0000000 | -0.0001507 | other | 60 |
| 145\_\_\_ | 0.0001239 | 0.0000595 | -0.0000644 | other | 73 |
| 146\_\_\_ | 0.0001222 | 0.0000387 | -0.0000834 | other | 74 |
| 147\_\_\_ | 0.0000815 | 0.0000912 | 0.0000098 | other | 111 |
| 148\_\_\_ | 0.0000194 | 0.0000087 | -0.0000106 | other | 467 |
| 149\_\_\_ | 0.0000136 | 0.0000072 | -0.0000064 | other | 664 |
| 15\_\_\_ | 0.0000109 | 0.0000250 | 0.0000141 | other | 829 |
| 150\_\_\_ | 0.0000130 | 0.0000037 | -0.0000093 | other | 694 |
| 151\_\_\_ | 0.0003931 | 0.0001186 | -0.0002745 | other | 23 |
| 152\_\_\_ | 0.0002318 | 0.0000684 | -0.0001634 | other | 39 |
| 153\_\_\_ | 0.0001586 | 0.0001557 | -0.0000029 | other | 57 |
| 154\_\_\_ | 0.0000701 | 0.0000063 | -0.0000638 | other | 129 |
| 155\_\_\_ | 0.0000655 | 0.0000117 | -0.0000538 | other | 138 |
| 156\_\_\_ | 0.0000650 | 0.0000442 | -0.0000208 | other | 139 |
| 157\_\_\_ | 0.0000363 | 0.0000161 | -0.0000202 | other | 249 |
| 158\_\_\_ | 0.0000576 | 0.0000524 | -0.0000052 | other | 157 |
| 159\_\_\_ | 0.0000459 | 0.0000186 | -0.0000273 | other | 197 |
| 16\_\_\_ | 0.0000690 | 0.0000939 | 0.0000249 | other | 131 |
| 160\_\_\_ | 0.0000400 | 0.0000212 | -0.0000188 | other | 226 |
| 161\_\_\_ | 0.0000262 | 0.0000000 | -0.0000262 | other | 345 |
| 162\_\_\_ | 0.0001370 | 0.0000577 | -0.0000793 | other | 66 |
| 163\_\_\_ | 0.0001027 | 0.0000169 | -0.0000858 | other | 88 |
| 164\_\_\_ | 0.0000800 | 0.0000395 | -0.0000406 | other | 113 |
| 165\_\_\_ | 0.0000532 | 0.0000000 | -0.0000532 | other | 170 |
| 166\_\_\_ | 0.0000422 | 0.0000034 | -0.0000388 | other | 214 |
| 167\_\_\_ | 0.0000913 | 0.0000250 | -0.0000664 | other | 99 |
| 168\_\_\_ | 0.0002825 | 0.0000485 | -0.0002341 | other | 32 |
| 169\_\_\_ | 0.0001674 | 0.0000278 | -0.0001396 | other | 54 |
| 17\_\_\_ | 0.0000619 | 0.0000891 | 0.0000271 | other | 146 |
| 170\_\_\_ | 0.0001206 | 0.0000380 | -0.0000826 | other | 75 |
| 171\_\_\_ | 0.0001206 | 0.0000569 | -0.0000636 | other | 75 |
| 172\_\_\_ | 0.0001206 | 0.0001042 | -0.0000163 | other | 75 |
| 173\_\_\_ | 0.0001174 | 0.0000460 | -0.0000715 | other | 77 |
| 174\_\_\_ | 0.0001064 | 0.0000016 | -0.0001047 | other | 85 |
| 175\_\_\_ | 0.0000878 | 0.0000067 | -0.0000811 | other | 103 |
| 176\_\_\_ | 0.0000508 | 0.0000000 | -0.0000508 | other | 178 |
| 177\_\_\_ | 0.0000497 | 0.0000165 | -0.0000332 | other | 182 |
| 178\_\_\_ | 0.0000363 | 0.0000000 | -0.0000363 | other | 249 |
| 179\_\_\_ | 0.0000282 | 0.0000123 | -0.0000158 | other | 321 |
| 18\_\_\_ | 0.0002379 | 0.0000277 | -0.0002102 | other | 38 |
| 180\_\_\_ | 0.0000223 | 0.0000021 | -0.0000202 | other | 406 |
| 181\_\_\_ | 0.0000158 | 0.0000098 | -0.0000060 | other | 571 |
| 182\_\_\_ | 0.0000138 | 0.0000022 | -0.0000116 | other | 654 |
| 183\_\_\_ | 0.0000580 | 0.0000649 | 0.0000069 | other | 156 |
| 184\_\_\_ | 0.0000541 | 0.0000312 | -0.0000229 | other | 167 |
| 185\_\_\_ | 0.0000837 | 0.0000606 | -0.0000231 | other | 108 |
| 186\_\_\_ | 0.0000807 | 0.0000000 | -0.0000807 | other | 112 |
| 187\_\_\_ | 0.0000793 | 0.0000233 | -0.0000560 | other | 114 |
| 188\_\_\_ | 0.0000786 | 0.0000388 | -0.0000398 | other | 115 |
| 189\_\_\_ | 0.0006028 | 0.0006558 | 0.0000530 | other | 15 |
| 19\_\_\_ | 0.0000983 | 0.0000612 | -0.0000371 | other | 92 |
| 190\_\_\_ | 0.0002659 | 0.0001347 | -0.0001313 | other | 34 |
| 191\_\_\_ | 0.0002444 | 0.0000000 | -0.0002443 | other | 37 |
| 192\_\_\_ | 0.0000773 | 0.0000261 | -0.0000512 | other | 117 |
| 193\_\_\_ | 0.0000675 | 0.0000413 | -0.0000262 | other | 134 |
| 194\_\_\_ | 0.0005318 | 0.0000749 | -0.0004569 | other | 17 |
| 195\_\_\_ | 0.0001924 | 0.0001879 | -0.0000045 | other | 47 |
| 196\_\_\_ | 0.0000690 | 0.0000417 | -0.0000273 | other | 131 |
| 197\_\_\_ | 0.0000655 | 0.0000062 | -0.0000593 | other | 138 |
| 198\_\_\_ | 0.0000655 | 0.0000248 | -0.0000407 | other | 138 |
| 199\_\_\_ | 0.0000619 | 0.0000358 | -0.0000262 | other | 146 |
| 2\_\_\_ | 0.0000853 | 0.0001408 | 0.0000555 | other | 106 |
| 20\_\_\_ | 0.0000437 | 0.0000500 | 0.0000063 | other | 207 |
| 200\_\_\_ | 0.0000569 | 0.0000473 | -0.0000096 | other | 159 |
| 201\_\_\_ | 0.0000385 | 0.0000117 | -0.0000267 | other | 235 |
| 202\_\_\_ | 0.0000363 | 0.0000048 | -0.0000315 | other | 249 |
| 203\_\_\_ | 0.0000349 | 0.0000034 | -0.0000316 | other | 259 |
| 204\_\_\_ | 0.0004521 | 0.0000809 | -0.0003712 | other | 20 |
| 205\_\_\_ | 0.0002009 | 0.0000407 | -0.0001602 | other | 45 |
| 206\_\_\_ | 0.0001586 | 0.0000272 | -0.0001314 | other | 57 |
| 207\_\_\_ | 0.0000468 | 0.0000039 | -0.0000430 | other | 193 |
| 208\_\_\_ | 0.0000464 | 0.0000018 | -0.0000446 | other | 195 |
| 209\_\_\_ | 0.0000402 | 0.0000037 | -0.0000365 | other | 225 |
| 21\_\_\_ | 0.0000437 | 0.0000505 | 0.0000068 | other | 207 |
| 210\_\_\_ | 0.0000349 | 0.0000032 | -0.0000318 | other | 259 |
| 211\_\_\_ | 0.0004759 | 0.0003479 | -0.0001280 | other | 19 |
| 212\_\_\_ | 0.0010046 | 0.0001252 | -0.0008794 | neo-endemic | 9 |
| 213\_\_\_ | 0.0003767 | 0.0000000 | -0.0003767 | other | 24 |
| 214\_\_\_ | 0.0002917 | 0.0001446 | -0.0001471 | other | 31 |
| 215\_\_\_ | 0.0002825 | 0.0000228 | -0.0002597 | other | 32 |
| 216\_\_\_ | 0.0002205 | 0.0000183 | -0.0002022 | other | 41 |
| 217\_\_\_ | 0.0001773 | 0.0000000 | -0.0001773 | other | 51 |
| 218\_\_\_ | 0.0001644 | 0.0000000 | -0.0001644 | other | 55 |
| 219\_\_\_ | 0.0001559 | 0.0000000 | -0.0001559 | other | 58 |
| 22\_\_\_ | 0.0000318 | 0.0000297 | -0.0000021 | other | 284 |
| 220\_\_\_ | 0.0001256 | 0.0000000 | -0.0001256 | other | 72 |
| 221\_\_\_ | 0.0000344 | 0.0000099 | -0.0000245 | other | 263 |
| 222\_\_\_ | 0.0000690 | 0.0000000 | -0.0000690 | other | 131 |
| 223\_\_\_ | 0.0000641 | 0.0000136 | -0.0000505 | other | 141 |
| 224\_\_\_ | 0.0002103 | 0.0000000 | -0.0002102 | other | 43 |
| 225\_\_\_ | 0.0003477 | 0.0000569 | -0.0002908 | other | 26 |
| 226\_\_\_ | 0.0002511 | 0.0000478 | -0.0002033 | other | 36 |
| 227\_\_\_ | 0.0002511 | 0.0000000 | -0.0002511 | other | 36 |
| 228\_\_\_ | 0.0002444 | 0.0000000 | -0.0002443 | other | 37 |
| 229\_\_\_ | 0.0002444 | 0.0000597 | -0.0001846 | other | 37 |
| 23\_\_\_ | 0.0000306 | 0.0001227 | 0.0000920 | other | 295 |
| 230\_\_\_ | 0.0000535 | 0.0000000 | -0.0000535 | other | 169 |
| 231\_\_\_ | 0.0000505 | 0.0000124 | -0.0000381 | other | 179 |
| 232\_\_\_ | 0.0000483 | 0.0000000 | -0.0000483 | other | 187 |
| 233\_\_\_ | 0.0001064 | 0.0000000 | -0.0001064 | other | 85 |
| 234\_\_\_ | 0.0001370 | 0.0000354 | -0.0001016 | other | 66 |
| 235\_\_\_ | 0.0001330 | 0.0000469 | -0.0000861 | other | 68 |
| 236\_\_\_ | 0.0001884 | 0.0000000 | -0.0001883 | other | 48 |
| 237\_\_\_ | 0.0002917 | 0.0000724 | -0.0002193 | other | 31 |
| 238\_\_\_ | 0.0004110 | 0.0000487 | -0.0003622 | other | 22 |
| 239\_\_\_ | 0.0001808 | 0.0000164 | -0.0001644 | other | 50 |
| 24\_\_\_ | 0.0000215 | 0.0000094 | -0.0000120 | other | 421 |
| 240\_\_\_ | 0.0000994 | 0.0000087 | -0.0000906 | other | 91 |
| 241\_\_\_ | 0.0000972 | 0.0000000 | -0.0000972 | other | 93 |
| 242\_\_\_ | 0.0005651 | 0.0001624 | -0.0004027 | other | 16 |
| 243\_\_\_ | 0.0000481 | 0.0000186 | -0.0000294 | other | 188 |
| 244\_\_\_ | 0.0000437 | 0.0000071 | -0.0000365 | other | 207 |
| 245\_\_\_ | 0.0000861 | 0.0000000 | -0.0000861 | other | 105 |
| 246\_\_\_ | 0.0000800 | 0.0000067 | -0.0000733 | other | 113 |
| 247\_\_\_ | 0.0000523 | 0.0000044 | -0.0000479 | other | 173 |
| 248\_\_\_ | 0.0000279 | 0.0000000 | -0.0000279 | other | 324 |
| 249\_\_\_ | 0.0000274 | 0.0000025 | -0.0000249 | other | 330 |
| 25\_\_\_ | 0.0000189 | 0.0000046 | -0.0000143 | other | 478 |
| 250\_\_\_ | 0.0000263 | 0.0000000 | -0.0000263 | other | 344 |
| 251\_\_\_ | 0.0000249 | 0.0000023 | -0.0000226 | other | 363 |
| 252\_\_\_ | 0.0000248 | 0.0000000 | -0.0000248 | other | 365 |
| 253\_\_\_ | 0.0000240 | 0.0000024 | -0.0000216 | other | 376 |
| 254\_\_\_ | 0.0000220 | 0.0000045 | -0.0000175 | other | 411 |
| 255\_\_\_ | 0.0000219 | 0.0000022 | -0.0000197 | other | 413 |
| 256\_\_\_ | 0.0000214 | 0.0000000 | -0.0000214 | other | 422 |
| 257\_\_\_ | 0.0000198 | 0.0000060 | -0.0000138 | other | 456 |
| 258\_\_\_ | 0.0000195 | 0.0000000 | -0.0000195 | other | 464 |
| 259\_\_\_ | 0.0000183 | 0.0000038 | -0.0000146 | other | 493 |
| 26\_\_\_ | 0.0000185 | 0.0000045 | -0.0000140 | other | 489 |
| 260\_\_\_ | 0.0000178 | 0.0000014 | -0.0000164 | other | 507 |
| 261\_\_\_ | 0.0000178 | 0.0000052 | -0.0000126 | other | 508 |
| 262\_\_\_ | 0.0000162 | 0.0000027 | -0.0000135 | other | 558 |
| 263\_\_\_ | 0.0000157 | 0.0000029 | -0.0000128 | other | 576 |
| 264\_\_\_ | 0.0000094 | 0.0000087 | -0.0000007 | other | 962 |
| 265\_\_\_ | 0.0000089 | 0.0000009 | -0.0000080 | other | 1015 |
| 266\_\_\_ | 0.0000089 | 0.0000212 | 0.0000123 | other | 1016 |
| 267\_\_\_ | 0.0001089 | 0.0001555 | 0.0000465 | other | 83 |
| 268\_\_\_ | 0.0018083 | 0.0000001 | -0.0018081 | neo-endemic | 5 |
| 269\_\_\_ | 0.0004110 | 0.0000285 | -0.0003824 | other | 22 |
| 27\_\_\_ | 0.0000123 | 0.0000141 | 0.0000018 | other | 734 |
| 270\_\_\_ | 0.0000461 | 0.0000078 | -0.0000383 | other | 196 |
| 271\_\_\_ | 0.0000402 | 0.0000042 | -0.0000360 | other | 225 |
| 272\_\_\_ | 0.0000400 | 0.0000000 | -0.0000400 | other | 226 |
| 273\_\_\_ | 0.0000391 | 0.0000131 | -0.0000260 | other | 231 |
| 274\_\_\_ | 0.0000388 | 0.0000049 | -0.0000339 | other | 233 |
| 275\_\_\_ | 0.0000385 | 0.0000173 | -0.0000212 | other | 235 |
| 276\_\_\_ | 0.0000374 | 0.0000067 | -0.0000307 | other | 242 |
| 277\_\_\_ | 0.0000290 | 0.0000075 | -0.0000215 | other | 312 |
| 278\_\_\_ | 0.0002103 | 0.0001287 | -0.0000816 | other | 43 |
| 279\_\_\_ | 0.0000837 | 0.0000000 | -0.0000837 | other | 108 |
| 28\_\_\_ | 0.0000111 | 0.0000227 | 0.0000116 | other | 814 |
| 280\_\_\_ | 0.0015069 | 0.0006490 | -0.0008578 | neo-endemic | 6 |
| 281\_\_\_ | 0.0006458 | 0.0022568 | 0.0016110 | other | 14 |
| 282\_\_\_ | 0.0000535 | 0.0000103 | -0.0000432 | other | 169 |
| 283\_\_\_ | 0.0000371 | 0.0000066 | -0.0000305 | other | 244 |
| 284\_\_\_ | 0.0000526 | 0.0000537 | 0.0000011 | other | 172 |
| 285\_\_\_ | 0.0001413 | 0.0000082 | -0.0001331 | other | 64 |
| 286\_\_\_ | 0.0000439 | 0.0000260 | -0.0000179 | other | 206 |
| 287\_\_\_ | 0.0000134 | 0.0000147 | 0.0000013 | other | 675 |
| 288\_\_\_ | 0.0005023 | 0.0002488 | -0.0002535 | other | 18 |
| 289\_\_\_ | 0.0000426 | 0.0000255 | -0.0000172 | other | 212 |
| 29\_\_\_ | 0.0000655 | 0.0000618 | -0.0000038 | other | 138 |
| 290\_\_\_ | 0.0000120 | 0.0000032 | -0.0000088 | other | 755 |
| 291\_\_\_ | 0.0000116 | 0.0000073 | -0.0000043 | other | 779 |
| 292\_\_\_ | 0.0000741 | 0.0000503 | -0.0000238 | other | 122 |
| 293\_\_\_ | 0.0001391 | 0.0000121 | -0.0001270 | other | 65 |
| 294\_\_\_ | 0.0000607 | 0.0000419 | -0.0000187 | other | 149 |
| 295\_\_\_ | 0.0000595 | 0.0000000 | -0.0000595 | other | 152 |
| 296\_\_\_ | 0.0000562 | 0.0000000 | -0.0000562 | other | 161 |
| 297\_\_\_ | 0.0000545 | 0.0000255 | -0.0000289 | other | 166 |
| 298\_\_\_ | 0.0000445 | 0.0000415 | -0.0000030 | other | 203 |
| 299\_\_\_ | 0.0000261 | 0.0000117 | -0.0000143 | other | 347 |
| 3\_\_\_ | 0.0000845 | 0.0000228 | -0.0000617 | other | 107 |
| 30\_\_\_ | 0.0000718 | 0.0002195 | 0.0001477 | other | 126 |
| 300\_\_\_ | 0.0000514 | 0.0000636 | 0.0000122 | other | 176 |
| 301\_\_\_ | 0.0000326 | 0.0000242 | -0.0000084 | other | 277 |
| 302\_\_\_ | 0.0001206 | 0.0000000 | -0.0001205 | other | 75 |
| 303\_\_\_ | 0.0000576 | 0.0000126 | -0.0000450 | other | 157 |
| 304\_\_\_ | 0.0000466 | 0.0000205 | -0.0000261 | other | 194 |
| 305\_\_\_ | 0.0000491 | 0.0000165 | -0.0000327 | other | 184 |
| 306\_\_\_ | 0.0003118 | 0.0000392 | -0.0002726 | other | 29 |
| 307\_\_\_ | 0.0003229 | 0.0001402 | -0.0001827 | other | 28 |
| 308\_\_\_ | 0.0001706 | 0.0000991 | -0.0000715 | other | 53 |
| 309\_\_\_ | 0.0001116 | 0.0000103 | -0.0001013 | other | 81 |
| 31\_\_\_ | 0.0000362 | 0.0000675 | 0.0000314 | other | 250 |
| 310\_\_\_ | 0.0000735 | 0.0000304 | -0.0000431 | other | 123 |
| 311\_\_\_ | 0.0000356 | 0.0000059 | -0.0000297 | other | 254 |
| 312\_\_\_ | 0.0000305 | 0.0000222 | -0.0000083 | other | 296 |
| 313\_\_\_ | 0.0000291 | 0.0000071 | -0.0000220 | other | 311 |
| 314\_\_\_ | 0.0000229 | 0.0000314 | 0.0000086 | other | 395 |
| 315\_\_\_ | 0.0000156 | 0.0000020 | -0.0000136 | other | 580 |
| 316\_\_\_ | 0.0000145 | 0.0000063 | -0.0000082 | other | 624 |
| 317\_\_\_ | 0.0000113 | 0.0000077 | -0.0000037 | other | 797 |
| 318\_\_\_ | 0.0000112 | 0.0000083 | -0.0000029 | other | 807 |
| 319\_\_\_ | 0.0000111 | 0.0000008 | -0.0000104 | other | 812 |
| 32\_\_\_ | 0.0002659 | 0.0003391 | 0.0000732 | other | 34 |
| 320\_\_\_ | 0.0000093 | 0.0000017 | -0.0000076 | other | 970 |
| 321\_\_\_ | 0.0004759 | 0.0004871 | 0.0000113 | other | 19 |
| 322\_\_\_ | 0.0004110 | 0.0002386 | -0.0001723 | other | 22 |
| 323\_\_\_ | 0.0002444 | 0.0000311 | -0.0002133 | other | 37 |
| 324\_\_\_ | 0.0001292 | 0.0001215 | -0.0000076 | other | 70 |
| 325\_\_\_ | 0.0000273 | 0.0000199 | -0.0000074 | other | 331 |
| 326\_\_\_ | 0.0000272 | 0.0000304 | 0.0000032 | other | 332 |
| 327\_\_\_ | 0.0000261 | 0.0000071 | -0.0000189 | other | 347 |
| 328\_\_\_ | 0.0000256 | 0.0000096 | -0.0000160 | other | 353 |
| 329\_\_\_ | 0.0000845 | 0.0000200 | -0.0000645 | other | 107 |
| 33\_\_\_ | 0.0000837 | 0.0001185 | 0.0000348 | other | 108 |
| 330\_\_\_ | 0.0000822 | 0.0000482 | -0.0000340 | other | 110 |
| 331\_\_\_ | 0.0000779 | 0.0000073 | -0.0000707 | other | 116 |
| 332\_\_\_ | 0.0000632 | 0.0000053 | -0.0000579 | other | 143 |
| 333\_\_\_ | 0.0001808 | 0.0000203 | -0.0001605 | other | 50 |
| 334\_\_\_ | 0.0000459 | 0.0000000 | -0.0000459 | other | 197 |
| 335\_\_\_ | 0.0003014 | 0.0000430 | -0.0002584 | other | 30 |
| 336\_\_\_ | 0.0002825 | 0.0000226 | -0.0002599 | other | 32 |
| 337\_\_\_ | 0.0000558 | 0.0000134 | -0.0000424 | other | 162 |
| 338\_\_\_ | 0.0000514 | 0.0000041 | -0.0000472 | other | 176 |
| 339\_\_\_ | 0.0000330 | 0.0000084 | -0.0000246 | other | 274 |
| 34\_\_\_ | 0.0000641 | 0.0000895 | 0.0000254 | other | 141 |
| 340\_\_\_ | 0.0000267 | 0.0000132 | -0.0000135 | other | 339 |
| 341\_\_\_ | 0.0000248 | 0.0000065 | -0.0000182 | other | 365 |
| 342\_\_\_ | 0.0000187 | 0.0000037 | -0.0000150 | other | 484 |
| 343\_\_\_ | 0.0000180 | 0.0000031 | -0.0000150 | other | 502 |
| 344\_\_\_ | 0.0000179 | 0.0000088 | -0.0000091 | other | 506 |
| 345\_\_\_ | 0.0000072 | 0.0000033 | -0.0000039 | other | 1250 |
| 346\_\_\_ | 0.0000071 | 0.0000052 | -0.0000020 | other | 1273 |
| 347\_\_\_ | 0.0000071 | 0.0000011 | -0.0000060 | other | 1274 |
| 348\_\_\_ | 0.0000060 | 0.0000010 | -0.0000051 | other | 1501 |
| 349\_\_\_ | 0.0000059 | 0.0000030 | -0.0000029 | other | 1545 |
| 35\_\_\_ | 0.0000246 | 0.0000400 | 0.0000155 | other | 368 |
| 350\_\_\_ | 0.0000058 | 0.0000017 | -0.0000041 | other | 1571 |
| 351\_\_\_ | 0.0000057 | 0.0000024 | -0.0000033 | other | 1588 |
| 352\_\_\_ | 0.0001370 | 0.0003417 | 0.0002047 | other | 66 |
| 353\_\_\_ | 0.0000476 | 0.0001998 | 0.0001522 | other | 190 |
| 354\_\_\_ | 0.0003477 | 0.0000816 | -0.0002662 | other | 26 |
| 355\_\_\_ | 0.0002583 | 0.0002895 | 0.0000311 | other | 35 |
| 356\_\_\_ | 0.0000402 | 0.0000188 | -0.0000214 | other | 225 |
| 357\_\_\_ | 0.0000348 | 0.0000084 | -0.0000264 | other | 260 |
| 358\_\_\_ | 0.0000296 | 0.0001035 | 0.0000738 | other | 305 |
| 359\_\_\_ | 0.0001144 | 0.0000000 | -0.0001144 | other | 79 |
| 36\_\_\_ | 0.0018083 | 0.0008962 | -0.0009120 | neo-endemic | 5 |
| 360\_\_\_ | 0.0000302 | 0.0000210 | -0.0000092 | other | 299 |
| 361\_\_\_ | 0.0000246 | 0.0000089 | -0.0000157 | other | 368 |
| 362\_\_\_ | 0.0000244 | 0.0000125 | -0.0000120 | other | 370 |
| 363\_\_\_ | 0.0000231 | 0.0000237 | 0.0000006 | other | 391 |
| 364\_\_\_ | 0.0000228 | 0.0000070 | -0.0000158 | other | 397 |
| 365\_\_\_ | 0.0001435 | 0.0000794 | -0.0000641 | other | 63 |
| 366\_\_\_ | 0.0001190 | 0.0000215 | -0.0000975 | other | 76 |
| 367\_\_\_ | 0.0002055 | 0.0000000 | -0.0002055 | other | 44 |
| 368\_\_\_ | 0.0000753 | 0.0000404 | -0.0000349 | other | 120 |
| 369\_\_\_ | 0.0000680 | 0.0000063 | -0.0000617 | other | 133 |
| 37\_\_\_ | 0.0011302 | 0.0019571 | 0.0008269 | meso-endemic | 8 |
| 370\_\_\_ | 0.0004521 | 0.0000000 | -0.0004520 | other | 20 |
| 371\_\_\_ | 0.0009041 | 0.0002487 | -0.0006554 | other | 10 |
| 372\_\_\_ | 0.0001966 | 0.0000463 | -0.0001503 | other | 46 |
| 373\_\_\_ | 0.0001370 | 0.0000600 | -0.0000770 | other | 66 |
| 374\_\_\_ | 0.0001005 | 0.0000069 | -0.0000936 | other | 90 |
| 375\_\_\_ | 0.0001103 | 0.0000195 | -0.0000907 | other | 82 |
| 376\_\_\_ | 0.0001089 | 0.0000376 | -0.0000714 | other | 83 |
| 377\_\_\_ | 0.0000942 | 0.0000294 | -0.0000647 | other | 96 |
| 378\_\_\_ | 0.0000587 | 0.0000759 | 0.0000172 | other | 154 |
| 379\_\_\_ | 0.0000368 | 0.0000217 | -0.0000150 | other | 246 |
| 38\_\_\_ | 0.0000196 | 0.0000109 | -0.0000087 | other | 462 |
| 380\_\_\_ | 0.0000316 | 0.0000121 | -0.0000195 | other | 286 |
| 381\_\_\_ | 0.0000305 | 0.0000083 | -0.0000222 | other | 296 |
| 382\_\_\_ | 0.0000293 | 0.0000109 | -0.0000184 | other | 309 |
| 383\_\_\_ | 0.0000223 | 0.0000088 | -0.0000134 | other | 406 |
| 384\_\_\_ | 0.0000180 | 0.0000017 | -0.0000163 | other | 502 |
| 385\_\_\_ | 0.0000753 | 0.0000483 | -0.0000270 | other | 120 |
| 386\_\_\_ | 0.0000517 | 0.0000358 | -0.0000159 | other | 175 |
| 387\_\_\_ | 0.0000318 | 0.0000080 | -0.0000239 | other | 284 |
| 388\_\_\_ | 0.0000267 | 0.0000429 | 0.0000161 | other | 338 |
| 389\_\_\_ | 0.0005651 | 0.0002012 | -0.0003638 | other | 16 |
| 39\_\_\_ | 0.0000779 | 0.0000194 | -0.0000585 | other | 116 |
| 390\_\_\_ | 0.0000431 | 0.0000585 | 0.0000154 | other | 210 |
| 391\_\_\_ | 0.0000424 | 0.0000284 | -0.0000140 | other | 213 |
| 392\_\_\_ | 0.0000422 | 0.0000041 | -0.0000382 | other | 214 |
| 393\_\_\_ | 0.0000407 | 0.0000045 | -0.0000362 | other | 222 |
| 394\_\_\_ | 0.0000272 | 0.0000035 | -0.0000237 | other | 332 |
| 395\_\_\_ | 0.0000203 | 0.0000020 | -0.0000183 | other | 445 |
| 396\_\_\_ | 0.0000150 | 0.0000020 | -0.0000130 | other | 601 |
| 397\_\_\_ | 0.0000122 | 0.0000039 | -0.0000084 | other | 740 |
| 398\_\_\_ | 0.0000094 | 0.0000020 | -0.0000074 | other | 957 |
| 399\_\_\_ | 0.0000087 | 0.0000042 | -0.0000045 | other | 1034 |
| 4\_\_\_ | 0.0000815 | 0.0000782 | -0.0000032 | other | 111 |
| 40\_\_\_ | 0.0001159 | 0.0000952 | -0.0000207 | other | 78 |
| 400\_\_\_ | 0.0000085 | 0.0000009 | -0.0000076 | other | 1060 |
| 401\_\_\_ | 0.0012916 | 0.0001536 | -0.0011380 | neo-endemic | 7 |
| 402\_\_\_ | 0.0006028 | 0.0006600 | 0.0000572 | other | 15 |
| 403\_\_\_ | 0.0001130 | 0.0000254 | -0.0000877 | other | 80 |
| 404\_\_\_ | 0.0002659 | 0.0001328 | -0.0001331 | other | 34 |
| 405\_\_\_ | 0.0001016 | 0.0000189 | -0.0000826 | other | 89 |
| 406\_\_\_ | 0.0000186 | 0.0000016 | -0.0000170 | other | 486 |
| 407\_\_\_ | 0.0000180 | 0.0000150 | -0.0000029 | other | 503 |
| 408\_\_\_ | 0.0000109 | 0.0000041 | -0.0000067 | other | 832 |
| 409\_\_\_ | 0.0000426 | 0.0000174 | -0.0000253 | other | 212 |
| 41\_\_\_ | 0.0000132 | 0.0000000 | -0.0000132 | other | 683 |
| 410\_\_\_ | 0.0000349 | 0.0000035 | -0.0000314 | other | 259 |
| 411\_\_\_ | 0.0000317 | 0.0000032 | -0.0000286 | other | 285 |
| 412\_\_\_ | 0.0000305 | 0.0000054 | -0.0000252 | other | 296 |
| 413\_\_\_ | 0.0000103 | 0.0000054 | -0.0000049 | other | 877 |
| 414\_\_\_ | 0.0000097 | 0.0000052 | -0.0000045 | other | 932 |
| 415\_\_\_ | 0.0000094 | 0.0000031 | -0.0000063 | other | 966 |
| 416\_\_\_ | 0.0000082 | 0.0000010 | -0.0000073 | other | 1100 |
| 417\_\_\_ | 0.0000476 | 0.0000181 | -0.0000295 | other | 190 |
| 418\_\_\_ | 0.0000450 | 0.0000080 | -0.0000370 | other | 201 |
| 419\_\_\_ | 0.0000431 | 0.0000424 | -0.0000007 | other | 210 |
| 42\_\_\_ | 0.0001808 | 0.0000151 | -0.0001657 | other | 50 |
| 420\_\_\_ | 0.0000314 | 0.0000041 | -0.0000273 | other | 288 |
| 421\_\_\_ | 0.0000236 | 0.0000137 | -0.0000099 | other | 383 |
| 422\_\_\_ | 0.0000583 | 0.0000070 | -0.0000513 | other | 155 |
| 423\_\_\_ | 0.0000395 | 0.0000096 | -0.0000299 | other | 229 |
| 424\_\_\_ | 0.0000170 | 0.0000039 | -0.0000131 | other | 532 |
| 425\_\_\_ | 0.0000140 | 0.0000048 | -0.0000092 | other | 644 |
| 426\_\_\_ | 0.0010046 | 0.0016691 | 0.0006645 | meso-endemic | 9 |
| 427\_\_\_ | 0.0006955 | 0.0000174 | -0.0006781 | other | 13 |
| 428\_\_\_ | 0.0000741 | 0.0001887 | 0.0001146 | other | 122 |
| 429\_\_\_ | 0.0003118 | 0.0002183 | -0.0000935 | other | 29 |
| 43\_\_\_ | 0.0000128 | 0.0000011 | -0.0000118 | other | 706 |
| 430\_\_\_ | 0.0001586 | 0.0001106 | -0.0000480 | other | 57 |
| 431\_\_\_ | 0.0005318 | 0.0001050 | -0.0004268 | other | 17 |
| 432\_\_\_ | 0.0001103 | 0.0000377 | -0.0000725 | other | 82 |
| 433\_\_\_ | 0.0000685 | 0.0000347 | -0.0000338 | other | 132 |
| 434\_\_\_ | 0.0000407 | 0.0000046 | -0.0000361 | other | 222 |
| 435\_\_\_ | 0.0000206 | 0.0000100 | -0.0000105 | other | 439 |
| 436\_\_\_ | 0.0000090 | 0.0000044 | -0.0000045 | other | 1007 |
| 437\_\_\_ | 0.0000085 | 0.0000010 | -0.0000076 | other | 1058 |
| 438\_\_\_ | 0.0000056 | 0.0000000 | -0.0000056 | other | 1618 |
| 439\_\_\_ | 0.0000056 | 0.0000004 | -0.0000051 | other | 1618 |
| 44\_\_\_ | 0.0000107 | 0.0000031 | -0.0000076 | other | 843 |
| 440\_\_\_ | 0.0000055 | 0.0000032 | -0.0000023 | other | 1640 |
| 441\_\_\_ | 0.0000655 | 0.0000664 | 0.0000009 | other | 138 |
| 442\_\_\_ | 0.0000701 | 0.0000000 | -0.0000701 | other | 129 |
| 443\_\_\_ | 0.0000932 | 0.0000000 | -0.0000932 | other | 97 |
| 444\_\_\_ | 0.0000878 | 0.0000380 | -0.0000498 | other | 103 |
| 445\_\_\_ | 0.0000441 | 0.0000107 | -0.0000334 | other | 205 |
| 446\_\_\_ | 0.0000368 | 0.0000410 | 0.0000043 | other | 246 |
| 447\_\_\_ | 0.0000276 | 0.0000082 | -0.0000194 | other | 328 |
| 448\_\_\_ | 0.0000520 | 0.0000406 | -0.0000113 | other | 174 |
| 449\_\_\_ | 0.0000407 | 0.0000081 | -0.0000326 | other | 222 |
| 45\_\_\_ | 0.0000102 | 0.0000046 | -0.0000055 | other | 890 |
| 450\_\_\_ | 0.0000292 | 0.0000082 | -0.0000210 | other | 310 |
| 451\_\_\_ | 0.0000267 | 0.0000068 | -0.0000199 | other | 338 |
| 452\_\_\_ | 0.0000962 | 0.0000383 | -0.0000579 | other | 94 |
| 453\_\_\_ | 0.0000624 | 0.0000196 | -0.0000427 | other | 145 |
| 454\_\_\_ | 0.0000591 | 0.0000160 | -0.0000431 | other | 153 |
| 455\_\_\_ | 0.0002103 | 0.0002624 | 0.0000522 | other | 43 |
| 456\_\_\_ | 0.0000665 | 0.0000186 | -0.0000479 | other | 136 |
| 457\_\_\_ | 0.0000353 | 0.0000069 | -0.0000284 | other | 256 |
| 458\_\_\_ | 0.0000296 | 0.0000015 | -0.0000282 | other | 305 |
| 459\_\_\_ | 0.0000288 | 0.0000049 | -0.0000239 | other | 314 |
| 46\_\_\_ | 0.0000099 | 0.0000009 | -0.0000090 | other | 917 |
| 460\_\_\_ | 0.0000173 | 0.0000066 | -0.0000106 | other | 524 |
| 461\_\_\_ | 0.0000123 | 0.0000027 | -0.0000095 | other | 738 |
| 462\_\_\_ | 0.0000360 | 0.0000000 | -0.0000360 | other | 251 |
| 463\_\_\_ | 0.0000240 | 0.0000111 | -0.0000129 | other | 377 |
| 464\_\_\_ | 0.0000235 | 0.0000181 | -0.0000054 | other | 385 |
| 465\_\_\_ | 0.0000233 | 0.0000000 | -0.0000233 | other | 388 |
| 466\_\_\_ | 0.0000206 | 0.0000020 | -0.0000186 | other | 438 |
| 467\_\_\_ | 0.0000437 | 0.0000339 | -0.0000098 | other | 207 |
| 468\_\_\_ | 0.0000380 | 0.0000028 | -0.0000352 | other | 238 |
| 469\_\_\_ | 0.0000248 | 0.0000172 | -0.0000076 | other | 364 |
| 47\_\_\_ | 0.0000076 | 0.0000108 | 0.0000032 | other | 1191 |
| 470\_\_\_ | 0.0002511 | 0.0001074 | -0.0001437 | other | 36 |
| 471\_\_\_ | 0.0000665 | 0.0000176 | -0.0000488 | other | 136 |
| 472\_\_\_ | 0.0001051 | 0.0000501 | -0.0000550 | other | 86 |
| 473\_\_\_ | 0.0000520 | 0.0000077 | -0.0000442 | other | 174 |
| 474\_\_\_ | 0.0000932 | 0.0000353 | -0.0000579 | other | 97 |
| 475\_\_\_ | 0.0000471 | 0.0000040 | -0.0000431 | other | 192 |
| 476\_\_\_ | 0.0000378 | 0.0000033 | -0.0000345 | other | 239 |
| 477\_\_\_ | 0.0000248 | 0.0000021 | -0.0000226 | other | 365 |
| 478\_\_\_ | 0.0000214 | 0.0000059 | -0.0000156 | other | 422 |
| 479\_\_\_ | 0.0000212 | 0.0000000 | -0.0000212 | other | 427 |
| 48\_\_\_ | 0.0000076 | 0.0000263 | 0.0000187 | other | 1191 |
| 480\_\_\_ | 0.0000200 | 0.0000019 | -0.0000181 | other | 452 |
| 481\_\_\_ | 0.0000117 | 0.0000010 | -0.0000106 | other | 776 |
| 482\_\_\_ | 0.0000097 | 0.0000021 | -0.0000077 | other | 931 |
| 483\_\_\_ | 0.0000095 | 0.0000032 | -0.0000064 | other | 951 |
| 484\_\_\_ | 0.0000090 | 0.0000042 | -0.0000048 | other | 1003 |
| 485\_\_\_ | 0.0000090 | 0.0000026 | -0.0000064 | other | 1006 |
| 486\_\_\_ | 0.0000077 | 0.0000034 | -0.0000042 | other | 1177 |
| 487\_\_\_ | 0.0000076 | 0.0000025 | -0.0000051 | other | 1194 |
| 488\_\_\_ | 0.0000040 | 0.0000009 | -0.0000031 | other | 2259 |
| 489\_\_\_ | 0.0000040 | 0.0000025 | -0.0000015 | other | 2263 |
| 49\_\_\_ | 0.0000076 | 0.0000026 | -0.0000050 | other | 1192 |
| 490\_\_\_ | 0.0000038 | 0.0000004 | -0.0000034 | other | 2390 |
| 491\_\_\_ | 0.0000038 | 0.0000012 | -0.0000025 | other | 2393 |
| 492\_\_\_ | 0.0000038 | 0.0000004 | -0.0000034 | other | 2393 |
| 493\_\_\_ | 0.0000038 | 0.0000004 | -0.0000033 | other | 2397 |
| 494\_\_\_ | 0.0000038 | 0.0000008 | -0.0000030 | other | 2397 |
| 495\_\_\_ | 0.0000037 | 0.0000011 | -0.0000026 | other | 2438 |
| 496\_\_\_ | 0.0000037 | 0.0000012 | -0.0000025 | other | 2457 |
| 497\_\_\_ | 0.0000037 | 0.0000007 | -0.0000030 | other | 2469 |
| 498\_\_\_ | 0.0000031 | 0.0000037 | 0.0000006 | other | 2888 |
| 499\_\_\_ | 0.0000031 | 0.0000013 | -0.0000018 | other | 2907 |
| 5\_\_\_ | 0.0002511 | 0.0003377 | 0.0000866 | other | 36 |
| 50\_\_\_ | 0.0000070 | 0.0000044 | -0.0000026 | other | 1294 |
| 500\_\_\_ | 0.0000031 | 0.0000008 | -0.0000023 | other | 2927 |
| 501\_\_\_ | 0.0000031 | 0.0000000 | -0.0000031 | other | 2956 |
| 502\_\_\_ | 0.0000031 | 0.0000043 | 0.0000013 | other | 2957 |
| 503\_\_\_ | 0.0000030 | 0.0000021 | -0.0000010 | other | 2971 |
| 504\_\_\_ | 0.0000030 | 0.0000014 | -0.0000016 | other | 3001 |
| 505\_\_\_ | 0.0000030 | 0.0000018 | -0.0000012 | other | 3015 |
| 506\_\_\_ | 0.0000030 | 0.0000006 | -0.0000024 | other | 3016 |
| 507\_\_\_ | 0.0000000 | 0.0000000 | 0.0000000 | other | NA |
| 51\_\_\_ | 0.0002583 | 0.0004403 | 0.0001820 | other | 35 |
| 52\_\_\_ | 0.0002055 | 0.0005069 | 0.0003014 | other | 44 |
| 53\_\_\_ | 0.0000753 | 0.0000360 | -0.0000393 | other | 120 |
| 54\_\_\_ | 0.0001310 | 0.0000254 | -0.0001057 | other | 69 |
| 55\_\_\_ | 0.0000665 | 0.0000184 | -0.0000481 | other | 136 |
| 56\_\_\_ | 0.0000395 | 0.0000141 | -0.0000253 | other | 229 |
| 57\_\_\_ | 0.0000395 | 0.0000539 | 0.0000144 | other | 229 |
| 58\_\_\_ | 0.0000391 | 0.0000578 | 0.0000187 | other | 231 |
| 59\_\_\_ | 0.0000329 | 0.0000291 | -0.0000038 | other | 275 |
| 6\_\_\_ | 0.0000580 | 0.0001104 | 0.0000524 | other | 156 |
| 60\_\_\_ | 0.0001884 | 0.0002229 | 0.0000345 | other | 48 |
| 61\_\_\_ | 0.0001586 | 0.0009409 | 0.0007823 | other | 57 |
| 62\_\_\_ | 0.0003931 | 0.0002258 | -0.0001673 | other | 23 |
| 63\_\_\_ | 0.0003014 | 0.0000417 | -0.0002596 | other | 30 |
| 64\_\_\_ | 0.0002825 | 0.0001212 | -0.0001613 | other | 32 |
| 65\_\_\_ | 0.0001845 | 0.0000865 | -0.0000980 | other | 49 |
| 66\_\_\_ | 0.0001773 | 0.0004721 | 0.0002948 | other | 51 |
| 67\_\_\_ | 0.0001532 | 0.0002036 | 0.0000503 | other | 59 |
| 68\_\_\_ | 0.0015069 | 0.0004977 | -0.0010091 | neo-endemic | 6 |
| 69\_\_\_ | 0.0002318 | 0.0001788 | -0.0000530 | other | 39 |
| 7\_\_\_ | 0.0000520 | 0.0000478 | -0.0000042 | other | 174 |
| 70\_\_\_ | 0.0001051 | 0.0000472 | -0.0000579 | other | 86 |
| 71\_\_\_ | 0.0001051 | 0.0004661 | 0.0003610 | other | 86 |
| 72\_\_\_ | 0.0005023 | 0.0014626 | 0.0009603 | other | 18 |
| 73\_\_\_ | 0.0003349 | 0.0005882 | 0.0002533 | other | 27 |
| 74\_\_\_ | 0.0006458 | 0.0000001 | -0.0006458 | other | 14 |
| 75\_\_\_ | 0.0000395 | 0.0001462 | 0.0001067 | other | 229 |
| 76\_\_\_ | 0.0000375 | 0.0000800 | 0.0000425 | other | 241 |
| 77\_\_\_ | 0.0000369 | 0.0000095 | -0.0000274 | other | 245 |
| 78\_\_\_ | 0.0000368 | 0.0000349 | -0.0000019 | other | 246 |
| 79\_\_\_ | 0.0000357 | 0.0000000 | -0.0000357 | other | 253 |
| 8\_\_\_ | 0.0000457 | 0.0000169 | -0.0000288 | other | 198 |
| 80\_\_\_ | 0.0000330 | 0.0000319 | -0.0000011 | other | 274 |
| 81\_\_\_ | 0.0000328 | 0.0000229 | -0.0000099 | other | 276 |
| 82\_\_\_ | 0.0000322 | 0.0000261 | -0.0000061 | other | 281 |
| 83\_\_\_ | 0.0000318 | 0.0000102 | -0.0000216 | other | 284 |
| 84\_\_\_ | 0.0000171 | 0.0000079 | -0.0000092 | other | 529 |
| 85\_\_\_ | 0.0000057 | 0.0000076 | 0.0000019 | other | 1588 |
| 86\_\_\_ | 0.0008219 | 0.0009148 | 0.0000929 | other | 11 |
| 87\_\_\_ | 0.0000222 | 0.0000144 | -0.0000078 | other | 408 |
| 88\_\_\_ | 0.0000214 | 0.0000275 | 0.0000061 | other | 422 |
| 89\_\_\_ | 0.0000211 | 0.0000456 | 0.0000245 | other | 429 |
| 9\_\_\_ | 0.0004305 | 0.0004711 | 0.0000405 | other | 21 |
| 90\_\_\_ | 0.0000637 | 0.0000160 | -0.0000476 | other | 142 |
| 91\_\_\_ | 0.0000637 | 0.0003152 | 0.0002515 | other | 142 |
| 92\_\_\_ | 0.0000580 | 0.0001907 | 0.0001328 | other | 156 |
| 93\_\_\_ | 0.0000807 | 0.0000086 | -0.0000722 | other | 112 |
| 94\_\_\_ | 0.0000735 | 0.0000879 | 0.0000144 | other | 123 |
| 95\_\_\_ | 0.0000324 | 0.0000573 | 0.0000249 | other | 279 |
| 96\_\_\_ | 0.0000232 | 0.0000247 | 0.0000015 | other | 389 |
| 97\_\_\_ | 0.0000119 | 0.0000114 | -0.0000005 | other | 762 |
| 98\_\_\_ | 0.0008219 | 0.0047984 | 0.0039764 | other | 11 |
| 99\_\_\_ | 0.0000439 | 0.0000033 | -0.0000406 | other | 206 |
| abbreviata | 0.0015069 | 0.0014361 | -0.0000708 | meso-endemic | 6 |
| acanthaster | 0.0009041 | 0.0025008 | 0.0015967 | other | 10 |
| acanthoclada | 0.0000913 | 0.0001844 | 0.0000931 | other | 99 |
| acinacea | 0.0000723 | 0.0000722 | -0.0000001 | other | 125 |
| aciphylla | 0.0007534 | 0.0003225 | -0.0004310 | other | 12 |
| acoma | 0.0006458 | 0.0003465 | -0.0002993 | other | 14 |
| acradenia | 0.0000514 | 0.0000562 | 0.0000048 | other | 176 |
| acrionastes | 0.0018083 | 0.0007818 | -0.0010264 | neo-endemic | 5 |
| acuaria | 0.0002825 | 0.0000000 | -0.0002825 | other | 32 |
| aculeatissima | 0.0002659 | 0.0007569 | 0.0004910 | other | 34 |
| acuminata | 0.0000655 | 0.0000189 | -0.0000467 | other | 138 |
| acutata | 0.0003617 | 0.0004275 | 0.0000659 | other | 25 |
| adinophylla | 0.0090413 | 0.0178638 | 0.0088225 | paleo-endemic | 1 |
| adnata | 0.0015069 | 0.0004797 | -0.0010272 | neo-endemic | 6 |
| adoxa | 0.0000786 | 0.0000182 | -0.0000604 | other | 115 |
| adsurgens | 0.0000583 | 0.0000105 | -0.0000478 | other | 155 |
| adunca | 0.0005651 | 0.0007894 | 0.0002243 | other | 16 |
| aemula | 0.0005023 | 0.0003673 | -0.0001350 | other | 18 |
| aestivalis | 0.0003767 | 0.0002934 | -0.0000833 | other | 24 |
| alata | 0.0001808 | 0.0004010 | 0.0002201 | other | 50 |
| alcockii | 0.0012916 | 0.0016062 | 0.0003146 | meso-endemic | 7 |
| alexandri | 0.0022603 | 0.0091523 | 0.0068919 | paleo-endemic | 4 |
| alpina | 0.0005318 | 0.0005692 | 0.0000374 | other | 17 |
| amblygona | 0.0001039 | 0.0002060 | 0.0001021 | other | 87 |
| amblyophylla | 0.0045207 | 0.0140873 | 0.0095667 | paleo-endemic | 2 |
| ammobia | 0.0008219 | 0.0007911 | -0.0000309 | other | 11 |
| ampliceps | 0.0001116 | 0.0000526 | -0.0000591 | other | 81 |
| amyctica | 0.0030138 | 0.0010406 | -0.0019732 | neo-endemic | 3 |
| anasilla | 0.0015069 | 0.0034281 | 0.0019212 | meso-endemic | 6 |
| anastema | 0.0015069 | 0.0022936 | 0.0007867 | meso-endemic | 6 |
| anaticeps | 0.0006955 | 0.0021430 | 0.0014475 | other | 13 |
| anceps | 0.0002511 | 0.0002016 | -0.0000496 | other | 36 |
| ancistrocarpa | 0.0000435 | 0.0000770 | 0.0000335 | other | 208 |
| ancistrophylla | 0.0001239 | 0.0000257 | -0.0000981 | other | 73 |
| andrewsii | 0.0002205 | 0.0004744 | 0.0002539 | other | 41 |
| aneura | 0.0000109 | 0.0000042 | -0.0000067 | other | 826 |
| anfractuosa | 0.0005651 | 0.0005546 | -0.0000105 | other | 16 |
| angusta | 0.0005651 | 0.0011561 | 0.0005910 | other | 16 |
| anthochaera | 0.0004759 | 0.0001535 | -0.0003223 | other | 19 |
| aphanoclada | 0.0022603 | 0.0069506 | 0.0046902 | paleo-endemic | 4 |
| aphylla | 0.0030138 | 0.0120579 | 0.0090441 | paleo-endemic | 3 |
| applanata | 0.0002825 | 0.0021956 | 0.0019130 | other | 32 |
| aprepta | 0.0002825 | 0.0000000 | -0.0002825 | other | 32 |
| aprica | 0.0090413 | 0.0022933 | -0.0067480 | neo-endemic | 1 |
| araneosa | 0.0030138 | 0.0009975 | -0.0020163 | neo-endemic | 3 |
| arcuatilis | 0.0015069 | 0.0001846 | -0.0013223 | neo-endemic | 6 |
| areolata | 0.0012916 | 0.0000001 | -0.0012915 | neo-endemic | 7 |
| argutifolia | 0.0030138 | 0.0050869 | 0.0020731 | meso-endemic | 3 |
| argyraea | 0.0001706 | 0.0003021 | 0.0001315 | other | 53 |
| argyrodendron | 0.0003767 | 0.0000544 | -0.0003223 | other | 24 |
| argyrophylla | 0.0002825 | 0.0000265 | -0.0002561 | other | 32 |
| arida | 0.0002583 | 0.0000791 | -0.0001793 | other | 35 |
| armitii | 0.0006028 | 0.0000911 | -0.0005116 | other | 15 |
| arrecta | 0.0006458 | 0.0003389 | -0.0003069 | other | 14 |
| ashbyae | 0.0009041 | 0.0036321 | 0.0027280 | other | 10 |
| aspera | 0.0003014 | 0.0000881 | -0.0002133 | other | 30 |
| asperulacea | 0.0001966 | 0.0001089 | -0.0000877 | other | 46 |
| assimilis | 0.0001051 | 0.0000353 | -0.0000698 | other | 86 |
| ataxiphylla | 0.0010046 | 0.0037297 | 0.0027251 | meso-endemic | 9 |
| atkinsiana | 0.0004759 | 0.0006311 | 0.0001552 | other | 19 |
| atrox | 0.0090413 | 0.0057638 | -0.0032775 | neo-endemic | 1 |
| attenuata | 0.0008219 | 0.0008138 | -0.0000081 | other | 11 |
| aulacocarpa | 0.0001051 | 0.0000508 | -0.0000543 | other | 86 |
| aulacophylla | 0.0005023 | 0.0002951 | -0.0002072 | other | 18 |
| auratiflora | 0.0030138 | 0.0127717 | 0.0097579 | paleo-endemic | 3 |
| aureocrinita | 0.0010046 | 0.0003679 | -0.0006366 | neo-endemic | 9 |
| auricoma | 0.0015069 | 0.0029563 | 0.0014494 | meso-endemic | 6 |
| auriculiformis | 0.0000837 | 0.0002506 | 0.0001669 | other | 108 |
| auronitens | 0.0007534 | 0.0010848 | 0.0003313 | other | 12 |
| ausfeldii | 0.0006028 | 0.0000000 | -0.0006027 | other | 15 |
| axillaris | 0.0008219 | 0.0025477 | 0.0017257 | other | 11 |
| ayersiana | 0.0000760 | 0.0000531 | -0.0000229 | other | 119 |
| baeuerlenii | 0.0006028 | 0.0006646 | 0.0000618 | other | 15 |
| baileyana | 0.0007534 | 0.0001348 | -0.0006186 | other | 12 |
| bakeri | 0.0006028 | 0.0027102 | 0.0021075 | other | 15 |
| balsamea | 0.0004759 | 0.0001139 | -0.0003620 | other | 19 |
| bancroftiorum | 0.0001144 | 0.0000000 | -0.0001144 | other | 79 |
| barattensis | 0.0030138 | 0.0071249 | 0.0041111 | paleo-endemic | 3 |
| barringtonensis | 0.0012916 | 0.0015173 | 0.0002257 | meso-endemic | 7 |
| bartleana | 0.0045207 | 0.0000004 | -0.0045203 | neo-endemic | 2 |
| basedowii | 0.0003229 | 0.0006167 | 0.0002938 | other | 28 |
| baueri | 0.0003617 | 0.0009259 | 0.0005642 | other | 25 |
| baxteri | 0.0010046 | 0.0038221 | 0.0028175 | meso-endemic | 9 |
| beadleana | 0.0045207 | 0.0137569 | 0.0092362 | paleo-endemic | 2 |
| beauverdiana | 0.0003477 | 0.0000405 | -0.0003072 | other | 26 |
| beckleri | 0.0001674 | 0.0001375 | -0.0000299 | other | 54 |
| benthamii | 0.0011302 | 0.0003384 | -0.0007917 | neo-endemic | 8 |
| betchei | 0.0007534 | 0.0002232 | -0.0005302 | other | 12 |
| bidentata | 0.0001674 | 0.0004496 | 0.0002822 | other | 54 |
| bifaria | 0.0015069 | 0.0007277 | -0.0007791 | neo-endemic | 6 |
| biflora | 0.0006955 | 0.0011955 | 0.0005000 | other | 13 |
| binata | 0.0007534 | 0.0004723 | -0.0002812 | other | 12 |
| binervata | 0.0002318 | 0.0000207 | -0.0002111 | other | 39 |
| binervia | 0.0003349 | 0.0002210 | -0.0001139 | other | 27 |
| bivenosa | 0.0000459 | 0.0000755 | 0.0000296 | other | 197 |
| blakei | 0.0001027 | 0.0006606 | 0.0005578 | other | 88 |
| blakelyi | 0.0003014 | 0.0004671 | 0.0001658 | other | 30 |
| blayana | 0.0030138 | 0.0005415 | -0.0024723 | neo-endemic | 3 |
| boormanii | 0.0004521 | 0.0005758 | 0.0001237 | other | 20 |
| brachybotrya | 0.0000706 | 0.0000147 | -0.0000559 | other | 128 |
| brachyclada | 0.0003617 | 0.0011760 | 0.0008143 | other | 25 |
| brachyphylla | 0.0005023 | 0.0008906 | 0.0003883 | other | 18 |
| brachypoda | 0.0018083 | 0.0022440 | 0.0004358 | meso-endemic | 5 |
| brachystachya | 0.0000380 | 0.0000252 | -0.0000128 | other | 238 |
| bracteolata | 0.0012916 | 0.0014367 | 0.0001450 | meso-endemic | 7 |
| brassii | 0.0005318 | 0.0011589 | 0.0006271 | other | 17 |
| brockii | 0.0010046 | 0.0029205 | 0.0019159 | meso-endemic | 9 |
| bromilowiana | 0.0010046 | 0.0017710 | 0.0007664 | meso-endemic | 9 |
| browniana | 0.0002318 | 0.0004626 | 0.0002308 | other | 39 |
| bulgaensis | 0.0022603 | 0.0006533 | -0.0016070 | neo-endemic | 4 |
| burkittii | 0.0000435 | 0.0000241 | -0.0000194 | other | 208 |
| buxifolia | 0.0000723 | 0.0001774 | 0.0001051 | other | 125 |
| caesiella | 0.0004521 | 0.0011153 | 0.0006632 | other | 20 |
| calamifolia | 0.0000913 | 0.0000522 | -0.0000392 | other | 99 |
| calantha | 0.0011302 | 0.0011742 | 0.0000441 | meso-endemic | 8 |
| calcicola | 0.0000869 | 0.0000082 | -0.0000787 | other | 104 |
| camptoclada | 0.0002825 | 0.0004585 | 0.0001760 | other | 32 |
| campylophylla | 0.0012916 | 0.0027230 | 0.0014314 | meso-endemic | 7 |
| cana | 0.0003617 | 0.0001285 | -0.0002331 | other | 25 |
| cangaiensis | 0.0045207 | 0.0004788 | -0.0040418 | neo-endemic | 2 |
| cardiophylla | 0.0005651 | 0.0002138 | -0.0003513 | other | 16 |
| carneorum | 0.0003229 | 0.0007079 | 0.0003850 | other | 28 |
| carnosula | 0.0022603 | 0.0053800 | 0.0031197 | paleo-endemic | 4 |
| caroleae | 0.0001190 | 0.0000498 | -0.0000692 | other | 76 |
| catenulata | 0.0001532 | 0.0000827 | -0.0000705 | other | 59 |
| cavealis | 0.0015069 | 0.0036569 | 0.0021500 | meso-endemic | 6 |
| cedroides | 0.0018083 | 0.0046625 | 0.0028543 | meso-endemic | 5 |
| celastrifolia | 0.0004305 | 0.0008603 | 0.0004297 | other | 21 |
| celsa | 0.0006028 | 0.0004034 | -0.0001993 | other | 15 |
| chartacea | 0.0005318 | 0.0005367 | 0.0000049 | other | 17 |
| cheelii | 0.0002659 | 0.0001810 | -0.0000850 | other | 34 |
| chinchillensis | 0.0010046 | 0.0017117 | 0.0007072 | meso-endemic | 9 |
| chisholmii | 0.0001349 | 0.0000000 | -0.0001349 | other | 67 |
| chrysotricha | 0.0045207 | 0.0052073 | 0.0006866 | meso-endemic | 2 |
| citrinoviridis | 0.0001773 | 0.0000205 | -0.0001568 | other | 51 |
| clelandii | 0.0002205 | 0.0001819 | -0.0000387 | other | 41 |
| clunies | 0.0045207 | 0.0018181 | -0.0027026 | neo-endemic | 2 |
| clydonophora | 0.0015069 | 0.0031881 | 0.0016812 | meso-endemic | 6 |
| cochlearis | 0.0001706 | 0.0002671 | 0.0000965 | other | 53 |
| cognata | 0.0005651 | 0.0001803 | -0.0003848 | other | 16 |
| colei | 0.0000468 | 0.0000000 | -0.0000468 | other | 193 |
| colletioides | 0.0000583 | 0.0001280 | 0.0000696 | other | 155 |
| complanata | 0.0001130 | 0.0002001 | 0.0000871 | other | 80 |
| concurrens | 0.0001884 | 0.0000278 | -0.0001605 | other | 48 |
| conferta | 0.0000895 | 0.0000645 | -0.0000250 | other | 101 |
| congesta | 0.0006028 | 0.0020491 | 0.0014463 | other | 15 |
| consobrina | 0.0006955 | 0.0012820 | 0.0005865 | other | 13 |
| conspersa | 0.0003014 | 0.0004142 | 0.0001128 | other | 30 |
| constablei | 0.0090413 | 0.0362315 | 0.0271902 | paleo-endemic | 1 |
| continua | 0.0001330 | 0.0001970 | 0.0000641 | other | 68 |
| coolgardiensis | 0.0001016 | 0.0001628 | 0.0000612 | other | 89 |
| coriacea | 0.0000471 | 0.0000000 | -0.0000471 | other | 192 |
| costata | 0.0008219 | 0.0021725 | 0.0013506 | other | 11 |
| courtii | 0.0090413 | 0.0060647 | -0.0029766 | neo-endemic | 1 |
| covenyi | 0.0015069 | 0.0031430 | 0.0016361 | meso-endemic | 6 |
| cowleana | 0.0000615 | 0.0001050 | 0.0000435 | other | 147 |
| craspedocarpa | 0.0002009 | 0.0000791 | -0.0001218 | other | 45 |
| crassa | 0.0000942 | 0.0000260 | -0.0000682 | other | 96 |
| crassicarpa | 0.0001482 | 0.0000000 | -0.0001482 | other | 61 |
| crassiuscula | 0.0006028 | 0.0006188 | 0.0000160 | other | 15 |
| cremiflora | 0.0006458 | 0.0005289 | -0.0001169 | other | 14 |
| cultriformis | 0.0002318 | 0.0003305 | 0.0000986 | other | 39 |
| cupularis | 0.0000800 | 0.0000341 | -0.0000459 | other | 113 |
| curranii | 0.0006955 | 0.0012426 | 0.0005471 | other | 13 |
| curvata | 0.0011302 | 0.0023828 | 0.0012527 | meso-endemic | 8 |
| cuspidifolia | 0.0003118 | 0.0021724 | 0.0018606 | other | 29 |
| cuthbertsonii | 0.0000545 | 0.0001963 | 0.0001418 | other | 166 |
| cyclops | 0.0001370 | 0.0006521 | 0.0005151 | other | 66 |
| cyperophylla | 0.0000983 | 0.0000736 | -0.0000247 | other | 92 |
| dangarensis | 0.0045207 | 0.0074512 | 0.0029306 | paleo-endemic | 2 |
| dawsonii | 0.0003014 | 0.0008916 | 0.0005903 | other | 30 |
| dealbata | 0.0000576 | 0.0000109 | -0.0000467 | other | 157 |
| deanei | 0.0000483 | 0.0000087 | -0.0000396 | other | 187 |
| debilis | 0.0006028 | 0.0005540 | -0.0000488 | other | 15 |
| declinata | 0.0022603 | 0.0019148 | -0.0003455 | meso-endemic | 4 |
| decora | 0.0000365 | 0.0000597 | 0.0000233 | other | 248 |
| decurrens | 0.0004521 | 0.0000000 | -0.0004520 | other | 20 |
| delibrata | 0.0004521 | 0.0003779 | -0.0000742 | other | 20 |
| delphina | 0.0007534 | 0.0043143 | 0.0035609 | other | 12 |
| dempsteri | 0.0005023 | 0.0018092 | 0.0013069 | other | 18 |
| denticulosa | 0.0011302 | 0.0002010 | -0.0009291 | neo-endemic | 8 |
| desmondii | 0.0012916 | 0.0003279 | -0.0009637 | neo-endemic | 7 |
| diallaga | 0.0045207 | 0.0003270 | -0.0041936 | neo-endemic | 2 |
| dictyoneura | 0.0045207 | 0.0033072 | -0.0012135 | neo-endemic | 2 |
| dictyophleba | 0.0000478 | 0.0000606 | 0.0000128 | other | 189 |
| didyma | 0.0045207 | 0.0000004 | -0.0045203 | neo-endemic | 2 |
| difficilis | 0.0000760 | 0.0000223 | -0.0000537 | other | 119 |
| dimidiata | 0.0001159 | 0.0002613 | 0.0001454 | other | 78 |
| diphylla | 0.0022603 | 0.0023750 | 0.0001146 | meso-endemic | 4 |
| disparrima | 0.0001064 | 0.0000730 | -0.0000334 | other | 85 |
| distans | 0.0006028 | 0.0005797 | -0.0000230 | other | 15 |
| dodonaeifolia | 0.0004521 | 0.0006013 | 0.0001492 | other | 20 |
| dolichophylla | 0.0022603 | 0.0108795 | 0.0086192 | paleo-endemic | 4 |
| doratoxylon | 0.0000822 | 0.0001727 | 0.0000905 | other | 110 |
| dorothea | 0.0009041 | 0.0009913 | 0.0000871 | other | 10 |
| drepanocarpa | 0.0000962 | 0.0000214 | -0.0000748 | other | 94 |
| drepanophylla | 0.0011302 | 0.0031551 | 0.0020249 | meso-endemic | 8 |
| drummondii | 0.0002444 | 0.0002884 | 0.0000441 | other | 37 |
| dunnii | 0.0002659 | 0.0003924 | 0.0001265 | other | 34 |
| elachantha | 0.0000587 | 0.0000000 | -0.0000587 | other | 154 |
| elata | 0.0006028 | 0.0017073 | 0.0011046 | other | 15 |
| elongata | 0.0002444 | 0.0001557 | -0.0000886 | other | 37 |
| empelioclada | 0.0018083 | 0.0021011 | 0.0002929 | meso-endemic | 5 |
| enervia | 0.0001966 | 0.0001332 | -0.0000633 | other | 46 |
| enterocarpa | 0.0005651 | 0.0003797 | -0.0001854 | other | 16 |
| epacantha | 0.0030138 | 0.0054094 | 0.0023956 | meso-endemic | 3 |
| eriopoda | 0.0001586 | 0.0007485 | 0.0005899 | other | 57 |
| errabunda | 0.0015069 | 0.0016904 | 0.0001835 | meso-endemic | 6 |
| estrophiolata | 0.0000853 | 0.0001650 | 0.0000797 | other | 106 |
| euthycarpa | 0.0001144 | 0.0000952 | -0.0000193 | other | 79 |
| excelsa | 0.0000478 | 0.0000184 | -0.0000295 | other | 189 |
| excentrica | 0.0009041 | 0.0007147 | -0.0001894 | other | 10 |
| exilis | 0.0015069 | 0.0024212 | 0.0009143 | meso-endemic | 6 |
| extensa | 0.0003349 | 0.0014162 | 0.0010813 | other | 27 |
| fagonioides | 0.0018083 | 0.0032413 | 0.0014330 | meso-endemic | 5 |
| falcata | 0.0001016 | 0.0000000 | -0.0001016 | other | 89 |
| falciformis | 0.0000766 | 0.0000482 | -0.0000284 | other | 118 |
| faucium | 0.0018083 | 0.0007260 | -0.0010823 | neo-endemic | 5 |
| fecunda | 0.0022603 | 0.0021429 | -0.0001174 | meso-endemic | 4 |
| filicifolia | 0.0002205 | 0.0000837 | -0.0001368 | other | 41 |
| fimbriata | 0.0001027 | 0.0000174 | -0.0000854 | other | 88 |
| flexifolia | 0.0001924 | 0.0001762 | -0.0000162 | other | 47 |
| floribunda | 0.0001076 | 0.0001948 | 0.0000872 | other | 84 |
| fragilis | 0.0002103 | 0.0001453 | -0.0000650 | other | 43 |
| fulva | 0.0018083 | 0.0013245 | -0.0004837 | meso-endemic | 5 |
| gardneri | 0.0009041 | 0.0003984 | -0.0005057 | other | 10 |
| gelasina | 0.0022603 | 0.0025579 | 0.0002976 | meso-endemic | 4 |
| genistifolia | 0.0000779 | 0.0001179 | 0.0000399 | other | 116 |
| georginae | 0.0000895 | 0.0001435 | 0.0000540 | other | 101 |
| gilbertii | 0.0006955 | 0.0005579 | -0.0001376 | other | 13 |
| gillii | 0.0022603 | 0.0025167 | 0.0002563 | meso-endemic | 4 |
| gittinsii | 0.0009041 | 0.0017413 | 0.0008372 | other | 10 |
| gladiiformis | 0.0002825 | 0.0005709 | 0.0002884 | other | 32 |
| glaucissima | 0.0018083 | 0.0023141 | 0.0005058 | meso-endemic | 5 |
| glaucocarpa | 0.0002055 | 0.0000506 | -0.0001549 | other | 44 |
| glaucoptera | 0.0003014 | 0.0001167 | -0.0001846 | other | 30 |
| gnidium | 0.0005318 | 0.0003442 | -0.0001877 | other | 17 |
| gonocarpa | 0.0001532 | 0.0008435 | 0.0006903 | other | 59 |
| gonoclada | 0.0000675 | 0.0000718 | 0.0000043 | other | 134 |
| gonophylla | 0.0003931 | 0.0010913 | 0.0006982 | other | 23 |
| gracillima | 0.0006955 | 0.0012706 | 0.0005751 | other | 13 |
| grandifolia | 0.0018083 | 0.0007016 | -0.0011067 | neo-endemic | 5 |
| grasbyi | 0.0001458 | 0.0000436 | -0.0001022 | other | 62 |
| gregorii | 0.0006955 | 0.0035898 | 0.0028943 | other | 13 |
| guinetii | 0.0045207 | 0.0009093 | -0.0036114 | neo-endemic | 2 |
| hakeoides | 0.0000437 | 0.0000421 | -0.0000016 | other | 207 |
| halliana | 0.0001458 | 0.0002543 | 0.0001084 | other | 62 |
| hamersleyensis | 0.0004521 | 0.0006617 | 0.0002096 | other | 20 |
| hammondii | 0.0000650 | 0.0000279 | -0.0000372 | other | 139 |
| harpophylla | 0.0000558 | 0.0000133 | -0.0000425 | other | 162 |
| harveyi | 0.0006028 | 0.0002449 | -0.0003578 | other | 15 |
| hastulata | 0.0006458 | 0.0000001 | -0.0006458 | other | 14 |
| havilandiorum | 0.0001391 | 0.0001968 | 0.0000577 | other | 65 |
| hemiteles | 0.0001174 | 0.0000407 | -0.0000767 | other | 77 |
| hemsleyi | 0.0001027 | 0.0005079 | 0.0004052 | other | 88 |
| heterochroa | 0.0011302 | 0.0026994 | 0.0015692 | meso-endemic | 8 |
| heteroclita | 0.0003931 | 0.0001467 | -0.0002464 | other | 23 |
| hexaneura | 0.0015069 | 0.0010806 | -0.0004263 | meso-endemic | 6 |
| hilliana | 0.0000523 | 0.0000215 | -0.0000307 | other | 173 |
| holosericea | 0.0000261 | 0.0000000 | -0.0000261 | other | 346 |
| hopperiana | 0.0010046 | 0.0000001 | -0.0010045 | neo-endemic | 9 |
| howittii | 0.0006458 | 0.0001393 | -0.0005065 | other | 14 |
| huegelii | 0.0004110 | 0.0013705 | 0.0009595 | other | 22 |
| hylonoma | 0.0045207 | 0.0501322 | 0.0456116 | paleo-endemic | 2 |
| hypermeces | 0.0090413 | 0.0107441 | 0.0017028 | meso-endemic | 1 |
| hystrix | 0.0007534 | 0.0024795 | 0.0017260 | other | 12 |
| imbricata | 0.0018083 | 0.0013854 | -0.0004229 | meso-endemic | 5 |
| implexa | 0.0000464 | 0.0001266 | 0.0000802 | other | 195 |
| inaequilatera | 0.0001051 | 0.0001812 | 0.0000761 | other | 86 |
| inceana | 0.0003229 | 0.0004896 | 0.0001667 | other | 28 |
| incrassata | 0.0010046 | 0.0057008 | 0.0046962 | paleo-endemic | 9 |
| ingramii | 0.0018083 | 0.0017079 | -0.0001004 | meso-endemic | 5 |
| irrorata | 0.0001130 | 0.0000188 | -0.0000942 | other | 80 |
| iteaphylla | 0.0002659 | 0.0005605 | 0.0002946 | other | 34 |
| ixiophylla | 0.0001435 | 0.0000977 | -0.0000458 | other | 63 |
| ixodes | 0.0002583 | 0.0006563 | 0.0003980 | other | 35 |
| jacksonioides | 0.0007534 | 0.0045560 | 0.0038025 | other | 12 |
| jamesiana | 0.0003767 | 0.0006492 | 0.0002724 | other | 24 |
| jennerae | 0.0001089 | 0.0001573 | 0.0000484 | other | 83 |
| jensenii | 0.0004759 | 0.0017114 | 0.0012356 | other | 19 |
| jibberdingensis | 0.0003617 | 0.0000899 | -0.0002718 | other | 25 |
| jonesii | 0.0015069 | 0.0002428 | -0.0012641 | neo-endemic | 6 |
| jucunda | 0.0002740 | 0.0000815 | -0.0001925 | other | 33 |
| julifera | 0.0000641 | 0.0000879 | 0.0000238 | other | 141 |
| karina | 0.0030138 | 0.0007507 | -0.0022631 | neo-endemic | 3 |
| kempeana | 0.0000353 | 0.0000095 | -0.0000258 | other | 256 |
| kybeanensis | 0.0007534 | 0.0000001 | -0.0007534 | other | 12 |
| lamprocarpa | 0.0000869 | 0.0000230 | -0.0000640 | other | 104 |
| lasiocalyx | 0.0001159 | 0.0002463 | 0.0001303 | other | 78 |
| latescens | 0.0001739 | 0.0003171 | 0.0001432 | other | 52 |
| latipes | 0.0001644 | 0.0000000 | -0.0001644 | other | 55 |
| latisepala | 0.0011302 | 0.0038085 | 0.0026784 | meso-endemic | 8 |
| latzii | 0.0012916 | 0.0003228 | -0.0009688 | neo-endemic | 7 |
| leiocalyx | 0.0000541 | 0.0000277 | -0.0000265 | other | 167 |
| leioderma | 0.0009041 | 0.0011651 | 0.0002610 | other | 10 |
| leptocarpa | 0.0000569 | 0.0000900 | 0.0000331 | other | 159 |
| leptoneura | 0.0022603 | 0.0128159 | 0.0105556 | paleo-endemic | 4 |
| leucoclada | 0.0001644 | 0.0001453 | -0.0000190 | other | 55 |
| leucolobia | 0.0004110 | 0.0004408 | 0.0000298 | other | 22 |
| ligulata | 0.0000132 | 0.0000093 | -0.0000039 | other | 683 |
| linearifolia | 0.0004110 | 0.0002358 | -0.0001752 | other | 22 |
| lineata | 0.0000913 | 0.0000128 | -0.0000785 | other | 99 |
| lineolata | 0.0002153 | 0.0003270 | 0.0001117 | other | 42 |
| linifolia | 0.0006458 | 0.0003512 | -0.0002946 | other | 14 |
| loderi | 0.0002153 | 0.0001508 | -0.0000645 | other | 42 |
| longifolia | 0.0000558 | 0.0000291 | -0.0000267 | other | 162 |
| longispicata | 0.0001089 | 0.0005300 | 0.0004211 | other | 83 |
| longispinea | 0.0001674 | 0.0006547 | 0.0004873 | other | 54 |
| longissima | 0.0002659 | 0.0002123 | -0.0000536 | other | 34 |
| loroloba | 0.0007534 | 0.0000001 | -0.0007534 | other | 12 |
| lycopodiifolia | 0.0001051 | 0.0001541 | 0.0000490 | other | 86 |
| lysiphloia | 0.0000445 | 0.0000161 | -0.0000284 | other | 203 |
| mabellae | 0.0006955 | 0.0003002 | -0.0003953 | other | 13 |
| macdonnellensis | 0.0001644 | 0.0004269 | 0.0002625 | other | 55 |
| mackeyana | 0.0003349 | 0.0001928 | -0.0001420 | other | 27 |
| macnuttiana | 0.0030138 | 0.0015428 | -0.0014710 | neo-endemic | 3 |
| maconochieana | 0.0007534 | 0.0011529 | 0.0003995 | other | 12 |
| macradenia | 0.0002009 | 0.0003951 | 0.0001941 | other | 45 |
| maitlandii | 0.0000461 | 0.0001546 | 0.0001085 | other | 196 |
| mangium | 0.0004305 | 0.0000000 | -0.0004305 | other | 21 |
| maranoensis | 0.0018083 | 0.0002107 | -0.0015976 | neo-endemic | 5 |
| marramamba | 0.0003477 | 0.0007539 | 0.0004062 | other | 26 |
| masliniana | 0.0004110 | 0.0005799 | 0.0001689 | other | 22 |
| mearnsii | 0.0000861 | 0.0000082 | -0.0000779 | other | 105 |
| meisneri | 0.0008219 | 0.0000764 | -0.0007456 | other | 11 |
| melanoxylon | 0.0000368 | 0.0000368 | 0.0000001 | other | 246 |
| melleodora | 0.0000348 | 0.0000430 | 0.0000083 | other | 260 |
| melvillei | 0.0000773 | 0.0000000 | -0.0000773 | other | 117 |
| menzelii | 0.0012916 | 0.0021945 | 0.0009029 | meso-endemic | 7 |
| microbotrya | 0.0001206 | 0.0000000 | -0.0001205 | other | 75 |
| microsperma | 0.0003477 | 0.0000000 | -0.0003477 | other | 26 |
| midgleyi | 0.0005651 | 0.0003780 | -0.0001871 | other | 16 |
| minyura | 0.0000932 | 0.0000000 | -0.0000932 | other | 97 |
| mitchellii | 0.0003118 | 0.0005187 | 0.0002070 | other | 29 |
| mollifolia | 0.0005023 | 0.0002514 | -0.0002509 | other | 18 |
| montana | 0.0000913 | 0.0000929 | 0.0000016 | other | 99 |
| monticola | 0.0000415 | 0.0001022 | 0.0000607 | other | 218 |
| mountfordiae | 0.0009041 | 0.0018211 | 0.0009170 | other | 10 |
| mucronata | 0.0001016 | 0.0000560 | -0.0000456 | other | 89 |
| muelleriana | 0.0003617 | 0.0002709 | -0.0000907 | other | 25 |
| multisiliqua | 0.0000695 | 0.0000241 | -0.0000455 | other | 130 |
| multispicata | 0.0001174 | 0.0000414 | -0.0000760 | other | 77 |
| murrayana | 0.0000255 | 0.0000172 | -0.0000082 | other | 355 |
| myrtifolia | 0.0000422 | 0.0001265 | 0.0000843 | other | 214 |
| nanodealbata | 0.0006458 | 0.0002465 | -0.0003993 | other | 14 |
| nematophylla | 0.0005023 | 0.0005499 | 0.0000476 | other | 18 |
| neriifolia | 0.0000895 | 0.0000642 | -0.0000253 | other | 101 |
| neurocarpa | 0.0001190 | 0.0000000 | -0.0001190 | other | 76 |
| neurophylla | 0.0001330 | 0.0002234 | 0.0000904 | other | 68 |
| nigricans | 0.0012916 | 0.0005550 | -0.0007366 | neo-endemic | 7 |
| notabilis | 0.0001391 | 0.0000662 | -0.0000729 | other | 65 |
| nuperrima | 0.0001310 | 0.0003761 | 0.0002450 | other | 69 |
| nyssophylla | 0.0000583 | 0.0000337 | -0.0000247 | other | 155 |
| obliquinervia | 0.0001586 | 0.0001016 | -0.0000570 | other | 57 |
| obtecta | 0.0008219 | 0.0048354 | 0.0040134 | other | 11 |
| obtusata | 0.0005023 | 0.0007313 | 0.0002290 | other | 18 |
| obtusifolia | 0.0002205 | 0.0003418 | 0.0001212 | other | 41 |
| oldfieldii | 0.0009041 | 0.0002881 | -0.0006160 | other | 10 |
| olgana | 0.0004305 | 0.0003791 | -0.0000514 | other | 21 |
| olsenii | 0.0045207 | 0.0000004 | -0.0045203 | neo-endemic | 2 |
| omalophylla | 0.0000800 | 0.0000000 | -0.0000800 | other | 113 |
| oncinocarpa | 0.0001330 | 0.0003339 | 0.0002009 | other | 68 |
| oncinophylla | 0.0008219 | 0.0000731 | -0.0007488 | other | 11 |
| orites | 0.0022603 | 0.0026259 | 0.0003656 | meso-endemic | 4 |
| orthocarpa | 0.0000895 | 0.0001359 | 0.0000464 | other | 101 |
| oshanesii | 0.0003617 | 0.0001790 | -0.0001827 | other | 25 |
| oswaldii | 0.0000146 | 0.0000097 | -0.0000049 | other | 620 |
| oxyclada | 0.0015069 | 0.0008782 | -0.0006287 | neo-endemic | 6 |
| pachyacra | 0.0001391 | 0.0000505 | -0.0000886 | other | 65 |
| pachycarpa | 0.0004759 | 0.0004835 | 0.0000076 | other | 19 |
| papyrocarpa | 0.0000807 | 0.0001343 | 0.0000536 | other | 112 |
| paradoxa | 0.0000478 | 0.0000792 | 0.0000314 | other | 189 |
| paraneura | 0.0000735 | 0.0000000 | -0.0000735 | other | 123 |
| parramattensis | 0.0002917 | 0.0000483 | -0.0002434 | other | 31 |
| parvipinnula | 0.0005651 | 0.0004284 | -0.0001367 | other | 16 |
| patagiata | 0.0004305 | 0.0005535 | 0.0001230 | other | 21 |
| pedina | 0.0030138 | 0.0000002 | -0.0030135 | neo-endemic | 3 |
| pedleyi | 0.0022603 | 0.0007846 | -0.0014757 | neo-endemic | 4 |
| pellita | 0.0002009 | 0.0000726 | -0.0001283 | other | 45 |
| pendula | 0.0000628 | 0.0000605 | -0.0000023 | other | 144 |
| penninervis | 0.0000576 | 0.0000083 | -0.0000493 | other | 157 |
| pentadenia | 0.0007534 | 0.0013640 | 0.0006106 | other | 12 |
| perangusta | 0.0018083 | 0.0002317 | -0.0015766 | neo-endemic | 5 |
| perryi | 0.0003349 | 0.0002962 | -0.0000386 | other | 27 |
| petraea | 0.0002511 | 0.0001253 | -0.0001259 | other | 36 |
| peuce | 0.0005023 | 0.0021979 | 0.0016956 | other | 18 |
| phlebopetala | 0.0011302 | 0.0024182 | 0.0012881 | meso-endemic | 8 |
| phlebophylla | 0.0030138 | 0.0040545 | 0.0010407 | meso-endemic | 3 |
| pickardii | 0.0012916 | 0.0023513 | 0.0010597 | meso-endemic | 7 |
| platycarpa | 0.0000483 | 0.0000273 | -0.0000210 | other | 187 |
| plectocarpa | 0.0000624 | 0.0000082 | -0.0000541 | other | 145 |
| podalyriifolia | 0.0001924 | 0.0000997 | -0.0000927 | other | 47 |
| polybotrya | 0.0002511 | 0.0000584 | -0.0001927 | other | 36 |
| porcata | 0.0045207 | 0.0131490 | 0.0086284 | paleo-endemic | 2 |
| praelongata | 0.0006458 | 0.0012620 | 0.0006162 | other | 14 |
| prainii | 0.0001016 | 0.0000357 | -0.0000659 | other | 89 |
| pravifolia | 0.0001330 | 0.0001335 | 0.0000005 | other | 68 |
| pravissima | 0.0002444 | 0.0002324 | -0.0000119 | other | 37 |
| producta | 0.0003229 | 0.0005649 | 0.0002420 | other | 28 |
| proiantha | 0.0015069 | 0.0010904 | -0.0004165 | meso-endemic | 6 |
| prominens | 0.0009041 | 0.0000001 | -0.0009041 | other | 10 |
| pruinocarpa | 0.0000779 | 0.0003362 | 0.0002583 | other | 116 |
| pruinosa | 0.0003931 | 0.0000372 | -0.0003559 | other | 23 |
| pterocaulon | 0.0045207 | 0.0099526 | 0.0054320 | paleo-endemic | 2 |
| ptychophylla | 0.0003617 | 0.0002377 | -0.0001239 | other | 25 |
| pubicosta | 0.0018083 | 0.0071960 | 0.0053877 | paleo-endemic | 5 |
| pubifolia | 0.0015069 | 0.0000001 | -0.0015068 | neo-endemic | 6 |
| pulchella | 0.0001089 | 0.0007290 | 0.0006201 | other | 83 |
| pustula | 0.0002825 | 0.0000429 | -0.0002396 | other | 32 |
| pycnantha | 0.0000580 | 0.0000108 | -0.0000472 | other | 156 |
| pycnostachya | 0.0015069 | 0.0022238 | 0.0007169 | meso-endemic | 6 |
| pygmaea | 0.0045207 | 0.0238021 | 0.0192815 | paleo-endemic | 2 |
| pyrifolia | 0.0000800 | 0.0000632 | -0.0000168 | other | 113 |
| ramulosa | 0.0000199 | 0.0000129 | -0.0000070 | other | 454 |
| redolens | 0.0006458 | 0.0012879 | 0.0006421 | other | 14 |
| repanda | 0.0045207 | 0.0007566 | -0.0037640 | neo-endemic | 2 |
| retinervis | 0.0006028 | 0.0000000 | -0.0006027 | other | 15 |
| retinodes | 0.0002055 | 0.0002181 | 0.0000126 | other | 44 |
| retivenea | 0.0000685 | 0.0000470 | -0.0000215 | other | 132 |
| rhamphophylla | 0.0030138 | 0.0041766 | 0.0011628 | meso-endemic | 3 |
| rhetinocarpa | 0.0006955 | 0.0002245 | -0.0004710 | other | 13 |
| rhigiophylla | 0.0005651 | 0.0001190 | -0.0004461 | other | 16 |
| rhodophloia | 0.0000845 | 0.0000311 | -0.0000534 | other | 107 |
| rigens | 0.0000431 | 0.0000564 | 0.0000134 | other | 210 |
| rivalis | 0.0005023 | 0.0003403 | -0.0001619 | other | 18 |
| rostellifera | 0.0001808 | 0.0000449 | -0.0001359 | other | 50 |
| rubida | 0.0000983 | 0.0001803 | 0.0000820 | other | 92 |
| ryaniana | 0.0018083 | 0.0018537 | 0.0000455 | meso-endemic | 5 |
| sabulosa | 0.0004305 | 0.0000000 | -0.0004305 | other | 21 |
| saliciformis | 0.0010046 | 0.0001850 | -0.0008195 | neo-endemic | 9 |
| salicina | 0.0000237 | 0.0000218 | -0.0000019 | other | 382 |
| saligna | 0.0000766 | 0.0001361 | 0.0000595 | other | 118 |
| saxatilis | 0.0011302 | 0.0026264 | 0.0014963 | meso-endemic | 8 |
| schinoides | 0.0011302 | 0.0000946 | -0.0010356 | neo-endemic | 8 |
| scirpifolia | 0.0004759 | 0.0000000 | -0.0004758 | other | 19 |
| sclerophylla | 0.0000822 | 0.0001333 | 0.0000511 | other | 110 |
| semicircinalis | 0.0022603 | 0.0085854 | 0.0063251 | paleo-endemic | 4 |
| semilunata | 0.0005023 | 0.0011956 | 0.0006933 | other | 18 |
| semitrullata | 0.0011302 | 0.0033723 | 0.0022422 | meso-endemic | 8 |
| sericoflora | 0.0005651 | 0.0006895 | 0.0001244 | other | 16 |
| sericophylla | 0.0000464 | 0.0000389 | -0.0000075 | other | 195 |
| sessilispica | 0.0003349 | 0.0003606 | 0.0000257 | other | 27 |
| shuttleworthii | 0.0006955 | 0.0016896 | 0.0009942 | other | 13 |
| sibilans | 0.0004305 | 0.0000777 | -0.0003528 | other | 21 |
| sibina | 0.0002379 | 0.0002373 | -0.0000006 | other | 38 |
| siculiformis | 0.0001586 | 0.0030597 | 0.0029011 | other | 57 |
| silvestris | 0.0005651 | 0.0001294 | -0.0004356 | other | 16 |
| simsii | 0.0000786 | 0.0000417 | -0.0000370 | other | 115 |
| simulans | 0.0022603 | 0.0075278 | 0.0052675 | paleo-endemic | 4 |
| sparsiflora | 0.0001391 | 0.0003417 | 0.0002026 | other | 65 |
| spathulifolia | 0.0003767 | 0.0013792 | 0.0010025 | other | 24 |
| spectabilis | 0.0001076 | 0.0003073 | 0.0001997 | other | 84 |
| spinescens | 0.0001089 | 0.0007008 | 0.0005919 | other | 83 |
| spirorbis | 0.0005318 | 0.0002273 | -0.0003046 | other | 17 |
| spongolitica | 0.0012916 | 0.0004427 | -0.0008489 | neo-endemic | 7 |
| stanleyi | 0.0030138 | 0.0000002 | -0.0030135 | neo-endemic | 3 |
| stenophylla | 0.0000321 | 0.0002240 | 0.0001919 | other | 282 |
| stigmatophylla | 0.0002205 | 0.0005408 | 0.0003203 | other | 41 |
| stipuligera | 0.0000660 | 0.0000889 | 0.0000229 | other | 137 |
| storyi | 0.0018083 | 0.0015320 | -0.0002763 | meso-endemic | 5 |
| stowardii | 0.0000680 | 0.0000085 | -0.0000595 | other | 133 |
| striatifolia | 0.0006955 | 0.0002125 | -0.0004830 | other | 13 |
| strongylophylla | 0.0002103 | 0.0003992 | 0.0001889 | other | 43 |
| suaveolens | 0.0000741 | 0.0001462 | 0.0000720 | other | 122 |
| subrigida | 0.0012916 | 0.0031112 | 0.0018196 | meso-endemic | 7 |
| subsessilis | 0.0010046 | 0.0005491 | -0.0004554 | meso-endemic | 9 |
| subtessarogona | 0.0006028 | 0.0002441 | -0.0003586 | other | 15 |
| subulata | 0.0003617 | 0.0002445 | -0.0001172 | other | 25 |
| sulcaticaulis | 0.0090413 | 0.0069265 | -0.0021148 | neo-endemic | 1 |
| synchronicia | 0.0001005 | 0.0002668 | 0.0001663 | other | 90 |
| tarculensis | 0.0003931 | 0.0000000 | -0.0003931 | other | 23 |
| telmica | 0.0030138 | 0.0046789 | 0.0016652 | meso-endemic | 3 |
| tenuinervis | 0.0011302 | 0.0004150 | -0.0007151 | neo-endemic | 8 |
| tenuispica | 0.0008219 | 0.0018347 | 0.0010127 | other | 11 |
| tenuissima | 0.0000417 | 0.0000734 | 0.0000318 | other | 217 |
| tephrina | 0.0001310 | 0.0000892 | -0.0000419 | other | 69 |
| terminalis | 0.0001016 | 0.0004265 | 0.0003249 | other | 89 |
| tessellata | 0.0012916 | 0.0006443 | -0.0006474 | neo-endemic | 7 |
| tetragonophylla | 0.0000193 | 0.0000367 | 0.0000173 | other | 468 |
| thomsonii | 0.0002740 | 0.0000233 | -0.0002507 | other | 33 |
| tindaleae | 0.0006458 | 0.0002756 | -0.0003702 | other | 14 |
| torulosa | 0.0000452 | 0.0000396 | -0.0000056 | other | 200 |
| trachycarpa | 0.0002260 | 0.0000000 | -0.0002260 | other | 40 |
| trachyphloia | 0.0012916 | 0.0005363 | -0.0007553 | neo-endemic | 7 |
| translucens | 0.0001239 | 0.0000734 | -0.0000504 | other | 73 |
| triptera | 0.0001330 | 0.0001211 | -0.0000118 | other | 68 |
| triquetra | 0.0003477 | 0.0001481 | -0.0001997 | other | 26 |
| tropica | 0.0001966 | 0.0004125 | 0.0002160 | other | 46 |
| tumida | 0.0000443 | 0.0000636 | 0.0000192 | other | 204 |
| tysonii | 0.0003477 | 0.0002221 | -0.0001256 | other | 26 |
| umbellata | 0.0000807 | 0.0001547 | 0.0000739 | other | 112 |
| umbraculiformis | 0.0003477 | 0.0000000 | -0.0003477 | other | 26 |
| uncinata | 0.0009041 | 0.0010621 | 0.0001580 | other | 10 |
| undoolyana | 0.0022603 | 0.0003809 | -0.0018795 | neo-endemic | 4 |
| validinervia | 0.0001586 | 0.0005343 | 0.0003757 | other | 57 |
| venulosa | 0.0003014 | 0.0001217 | -0.0001797 | other | 30 |
| verniciflua | 0.0000637 | 0.0000151 | -0.0000486 | other | 142 |
| verricula | 0.0002009 | 0.0002819 | 0.0000810 | other | 45 |
| vestita | 0.0005023 | 0.0003078 | -0.0001945 | other | 18 |
| victoriae | 0.0000151 | 0.0000427 | 0.0000276 | other | 598 |
| viscidula | 0.0002103 | 0.0012487 | 0.0010384 | other | 43 |
| wanyu | 0.0003014 | 0.0000650 | -0.0002364 | other | 30 |
| wattsiana | 0.0008219 | 0.0000001 | -0.0008219 | other | 11 |
| wilhelmiana | 0.0001076 | 0.0000775 | -0.0000301 | other | 84 |
| woodmaniorum | 0.0030138 | 0.0050042 | 0.0019904 | meso-endemic | 3 |
| xanthina | 0.0010046 | 0.0004608 | -0.0005438 | neo-endemic | 9 |
| xiphophylla | 0.0001924 | 0.0002994 | 0.0001070 | other | 47 |
| yirrkallensis | 0.0003118 | 0.0005226 | 0.0002109 | other | 29 |
