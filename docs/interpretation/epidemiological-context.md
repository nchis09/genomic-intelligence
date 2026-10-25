# Interpreting epidemiological context

!!! info "Section in progress"
    This page will explain where the epidemiological layer comes from
    (WHO FluNet for influenza, HDX / EMPRES-i for Ebola and avian
    influenza), how `--epi_min_year` / `--epi_max_year` shape the historical
    comparison, and how R(t) at sampling is estimated.

## Key idea

The transmission and geographic panels are only as good as the epi window
fetched at run time. Historical questions need a window wide enough to
contain history — see
[Running the pipeline](../getting-started/running-the-pipeline.md#the-example-run).
