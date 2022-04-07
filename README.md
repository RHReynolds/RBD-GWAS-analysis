
<!-- README.md is generated from README.Rmd. Please edit that file -->
# Background

This repository contains scripts used to perform a colocalisation analysis using the REM sleep behavior disorder (RBD) GWAS and two eQTL datasets:

-   [eQTLGen](https://pubmed.ncbi.nlm.nih.gov/34475573/)
-   [PsychENCODE](https://www.ncbi.nlm.nih.gov/pubmed/30545857)

# Citation

If you use any of the code or data from this repository, please cite our [paper](https://www.medrxiv.org/content/10.1101/2021.09.08.21254232v2). Further, if you use any of the software used within this repository (e.g. `coloc`, `colochelpR`, etc.) please make sure to cite the software appropriately.

# License

The code in this repository is released under an MIT license. This repository is distributed in the hope that it will be useful to the wider community, but without any warranty of any kind. Please see the [LICENSE](LICENSE) file for more details.

# Code contents

Scripts and results of analyses have been described in the following workflows:

1.  [RBD\_coloc.Rmd](./docs/RBD_coloc.Rmd): `Rmd` detailing coloc set up and analysis of results.
2.  [RBD\_coloc\_tissue\_cell\_specificity.Rmd](./docs/RBD_coloc_tissue_cell_specificity.Rmd): `Rmd` detailing tissue- and cell-type-specificity of genes identified by coloc.
3.  [RBD\_manuscript\_figures](./docs/RBD_manuscript_figures.Rmd): `Rmd` detailing code used to produce manuscript figures pertaining to coloc analyses.

All can be view interactively at: <https://rhreynolds.github.io/RBD-GWAS-analysis/>

Within this repository you will otherwise find:

<table>
<colgroup>
<col width="11%" />
<col width="88%" />
</colgroup>
<thead>
<tr class="header">
<th>Directory</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><a href="docs" class="uri">docs</a></td>
<td>Contains all <code>.Rmd</code>s and their corresponding <code>.html</code>s describing analyses performed for this project. These can be view interactively at: <a href="https://rhreynolds.github.io/RBD-GWAS-analysis/" class="uri">https://rhreynolds.github.io/RBD-GWAS-analysis/</a></td>
</tr>
<tr class="even">
<td><a href="logs" class="uri">logs</a></td>
<td>For any scripts that were run outside of an <code>.Rmd</code> (e.g. scripts from the <a href="scripts" class="uri">scripts</a> directory), a log file was recorded and can be accessed here.</td>
</tr>
<tr class="odd">
<td><a href="manuscript" class="uri">manuscript</a></td>
<td>Figures and tables produced for the manuscript.</td>
</tr>
<tr class="even">
<td><a href="results" class="uri">results</a></td>
<td>Results from all analyses.</td>
</tr>
<tr class="odd">
<td><a href="scripts" class="uri">scripts</a></td>
<td>Contains analysis scripts. Each script contains a one-line description and is also referenced in its corresponding <code>.Rmd</code>.</td>
</tr>
<tr class="even">
<td><a href="R" class="uri">R</a></td>
<td>Various functions called in <a href="docs" class="uri">docs</a> and <a href="scripts" class="uri">scripts</a>.</td>
</tr>
</tbody>
</table>
