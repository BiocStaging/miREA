# miREA 0.99.0

* Initial Bioconductor submission.
* Four node-based and five edge-based miRNA-oriented enrichment methods:
  `TG_ORA`, `TG_Score`, `MiR_ORA`, `MiR_Score`, `Edge_ORA`, `Edge_Score1D`,
  `Edge_Score2D`, `Edge_Topology`, `Edge_Network`.
* `get_all_input_data()` to prepare method-specific input data from raw
  expression and differential-expression tables.
* `plot_summary()`, `plot_heatmap()`, and `plot_heatmap_sankey()` for
  visualizing and comparing enrichment results.
