# “OLD-GROWTH” IS NOT ALWAYS COMPLEX: ASSESSING LIDAR-DERIVED FOREST STRUCTURAL COMPLEXITY ACROSS SUCCESSIONAL STAGES ON VANCOUVER ISLAND

**Authors:**\
Luizmar de Assis Barros<sup>a</sup>, Karen Price<sup>b</sup>, Camile Sothe<sup>c</sup>, Chris Johnson<sup>a</sup>, Juan Pablo Ramírez-Delgado<sup>a</sup>, Xavier Llano<sup>a</sup>, Michelle Venter<sup>a</sup>, Oscar Venter<sup>a</sup>

a *University of Northern British Columbia, 3333 University Way, Prince George, V2N 4Z9, British Columbia, Canada*\
b *Independent Researcher, 1355 Malkow Road, Smithers, BC V0J 2N7, Canada*\
c *Planet Labs PBC, San Francisco, 695571, California, USA*

---

## Abstract

Forest structural complexity (FSC) summarizes the three-dimensional (3D) arrangement of forest elements that directly influence both species’ habitat availability and canopy-atmosphere interactions (e.g., light interception, moisture, and gas exchange). While often considered synonymous with forest succession, factors such as site quality, topography, climate, and disturbance history also shape FSC. In this context, we investigated whether old-growth forests were inherently characterized by high levels of structural complexity, as often assumed. We assessed six LiDAR-derived FSC indices (canopy rugosity, canopy entropy, foliage diversity, median absolute deviation of height, coefficient of variation of height, fractal dimensions) against 13 old-growth structural attributes (e.g., stand age, biomass, vertical complexity, top height) measured in 343 plots across Vancouver Island, British Columbia, Canada. Random forest models with 5-fold cross-validation were used to test the performance and redundancies of the FSC indices. We then compared FSC across age, forest maturity index (IMAT; aggregation of all old-growth structural attributes, except age), and site productivity (site index; expected tree height at 50 years old). Canopy entropy, canopy rugosity, and median absolute deviation of height were highly redundant. However, canopy entropy had the strongest statistical relationship with old-growth structural attributes (R2 =0.57 ± 0.08), with top tree height, live tree volume, max diameter at breast height (DBH), and age as key predictors. We found an asymptotic relationship between FSC and the top height and biomass predictors: FSC plateaued beyond thresholds of approximately 45 m for height and 900 m³ ha⁻¹ for live tree volume. For DBH and age, the relationship with FSC was quadratic: FSC increased with DBH up to 250 cm and slightly declined thereafter. When assessing FSC across forest stand age, FSC peaked in late-mature forests (~210 years) but declined in older stands, particularly in lower-productivity stands (site index < 20m); old-growth forests on productive sites retained higher FSC (> 250 years and site index > 20m). The IMAT index had a stronger relationship with FSC than age but also showed an apparent plateau around 0.75 on a scale from 0 to 1.5. Our findings suggest that FSC is more linked to biomass and productivity than forest succession, making it valuable for describing forest development but insufficient to define old-growth.

---

## Index Terms
Forest structural complexity; Canopy Rugosity; Old-growth forests; LiDAR indices; Site productivity

---

## Main Result

![Forest Structural Complexity](image.jpeg)
**Figure 1** Scatter plot of FSC vs a) log transformed stand age and b) forest maturity index (IMAT). 

# R Scripts Description:
**- 1_BC_tree_data.R**\
  Pre-processing of British Columbia tree-level field inventory data and calculation of old-growth structural attributes\
**- 2_plot_level_comp_metrics.R**\
  Calculation of LiDAR-derived Forest Structural Complexity (FSC) indices at plot level\
**- 3_comp_met_analysis.R**\
  Random forest analysis of field-measured old-growth structures vs LiDAR FSC indices\
**- 4_Anova_analysis_Full.R**\
  Comparison of stand age groups, Maturity clusters, and Productivity vs FSC

# Data Summaries:
**data_Full.csv**: Field and LiDAR metrics unfiltered\
**data_75.csv**: Field and LiDAR metrics filtered with a voxel size of 0.4 m, resulting in a drop of 25% on the point cloud density\
**data_50.csv**: Field and LiDAR metrics filtered with a voxel size of 0.75 m, resulting in a drop of 50% on the point cloud density
