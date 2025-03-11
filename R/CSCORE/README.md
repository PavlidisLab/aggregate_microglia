This directory contains scripts to do the same scRNA-seq processing as done in
the TRsc project, but to keep raw counts (instead of CPM normalized) and feed
the resulting data into CSCORE. As before, input datasets are grouped by whether
they originated from Cell X Gene (and thus are standardized and a single script
could be used), or sourced from GEO.

All experiments were ran using the following commands:

```
bash batch_cellxgene_datasets.sh cellxgene_ids.tsv
bash batch_geo_datasets.sh geo_ids.tsv
```