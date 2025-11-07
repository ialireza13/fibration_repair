# Fibration Repair Algorithm
by David Phillips, Alireza Hashemi


Usage:
```
prohibit_file = "sample_data/collapsed_prohibited_edges.txt"
a, b = 1, 1

EdgesRemoved, EdgesAdded, G_result = repair_network("sample_data/cons_6_colors_collapsed.txt", "sample_data/collapsed_varshney.graph.txt", f"sample_data/collapsed_consensous_6_colors_o_", a, b, prohibit_file_path=None)
```