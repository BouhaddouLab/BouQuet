# Built-in tables

| File | What it is | Used by |
|---|---|---|
| `Regulatory_sites_2025-07-02.tsv.gz` | PhosphoSitePlus regulatory sites (downloaded 2025-07-02) | `get_psp_regulatory_sites()`, notebook D4 |
| `human_kinases.txt` | List of human kinases (column `Kinase`) | notebook D4 |
| `kinaseSiteList_BachmanGyoriSorger2019.csv.gz` | Kinase-substrate list (Bachman, Gyori & Sorger 2019) | `loadKinaseData()` |

Get the full path to any of these with `bouquet_file("human_kinases.txt")`.

**Updating PhosphoSitePlus:** download `Regulatory_sites.gz` from phosphosite.org (free login), save it here as
`Regulatory_sites_<YYYY-MM-DD>.tsv.gz` (delete the old one), and reinstall the package. A newer copy in a project's
`resources/` folder is always used before this one.

PhosphoSitePlus(R) is licensed under CC BY-NC-SA 3.0 (non-commercial use). Please acknowledge PhosphoSitePlus(R),
www.phosphosite.org, in publications.
