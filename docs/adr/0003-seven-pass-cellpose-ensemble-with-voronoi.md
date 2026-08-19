# 3. 7-Pass Cellpose and CPSAM Ensemble with Voronoi Expansion

Single-pass deep learning segmentation fails to capture heterogeneously sized cancer and immune cells in dense tissue sections. We decided to run an ensemble across 3 Cyto3 diameter scales, Hematoxylin nuclear decomposition, and 3 CPSAM segmentation passes, followed by fast overlap-filtered merging and Voronoi membrane expansion to recover complete whole-cell morphologies without artificial cell collisions.
