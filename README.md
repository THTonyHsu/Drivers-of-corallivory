# Data from: Corallivore assemblage, feeding behaviour, and cover of *Acropora* corals drive fish corallivory on tropical reefs

------------------------------------------------------------------------

This README.md file was generated on 25 Mar 2026 by **Tsai-Hsuan Tony Hsu** ([tsai-hsuan.hsu\@sydney.edu.au](mailto:tsai-hsuan.hsu@sydney.edu.au) / [tonytsaihsuanhsu\@gmail.com](mailto:tonytsaihsuanhsu@gmail.com)). Render [CorallivoryDrivers_Hsu_et_al_2026.qmd](CorallivoryDrivers_Hsu_et_al_2026.qmd) to generate HTML and extract R code.

1.  **Author Information**

    [Tsai-Hsuan Tony Hsu](https://orcid.org/0000-0002-9362-9769) <sup>1</sup>, [Andrew S. Hoey](https://orcid.org/0000-0002-4261-5594) <sup>2</sup>, [Renata Ferrari](https://orcid.org/0000-0002-5056-1178) <sup>3</sup>, [Maren Toor](https://orcid.org/0009-0004-9458-877X) <sup>3</sup>, [Sophie Gordon](https://orcid.org/0000-0003-1775-7221) <sup>3</sup>, [Tiny Remmers](https://orcid.org/0000-0002-0984-7940) <sup>2,3</sup>, [James Smallhorn-West](https://orcid.org/0009-0000-6459-9508) <sup>3</sup>, [Agustina Piccaluga](https://orcid.org/0009-0004-6571-2241) <sup>1</sup>, [Will F. Figueira](https://orcid.org/0000-0001-9472-8710) <sup>1</sup>

    1.  School of Life and Environmental Science, The University of Sydney, Sydney, NSW 2006, Australia
    2.  College of Science and Engineering, James Cook University, Townsville, QLD 4814, Australia
    3.  Australian Institute of Marine Science, Townsville, QLD 4810, Australia

2.  **Date of data collection**: Jan-May 2023

3.  **Geographic location of data collection**: Lizard Island, Davies Reef, and Heron Island in the northern, central and southern regions of the Great Barrier Reef (GBR), respectively, and Masig Island in the Torres Strait (TS) in Australia

4.  **Funding sources that supported the collection of the data**: This research was supported by the Reef Restoration and Adaptation Program, which is funded by a partnership between the Australian Government, Reef Trust, the Great Barrier Reef Foundation, and The Australian Institute of Marine Science, and by the Taiwan–University of Sydney Jointly Funded PhD Scholarship (USYD and Ministry of Education, Taiwan) awarded to T.-H. T. Hsu. The research was conducted under the Great Barrier Reef Marine Park permit G21/44774.1.

5.  **Recommended citation for this dataset:** Hsu T-H.T, Hoey A.S, Warren C.R, Marrable I, Figueira W.F (2026) Data from “Corallivore assemblage, feeding behaviour, and *Acropora* corals drive fish corallivory on tropical reefs.”. [Data set]. Sydney eScholarship. [doi](doi)

[Data and R script are also available through the GitHub repository <https://github.com/THTonyHsu/Drivers-of-corallivory.git> in order to replicate our analyses]

### Methodological information

A detailed description of data acquisition and processing can be found in the published manuscript in the *Ecosphere* (doixx).

### Description of data sturcture

1.  **Description of dataset**

    These data were generated to understand drivers of fish corallivory among coral reefs across \>1,800 km in northeast Australia

2.  **File List**

    **File 1**: Hsuetal_dataset_benthic_MFG.csv

    List of benthic taxa to morpho-functional groups

    **File 2**: Hsuetal_dataset_benthic_taxa.csv

    Benthic communites of each plot, characterized by per cent cover of 72 benthic taxa

    **File 3**: Hsuetal_dataset_Complexity.csv

    Structural complexity metrics for each plot

    **File 4**: Hsuetal_dataset_Corallivory.csv

    Feeding activities of obligate corallivores annotated in EventMeasure

    **File 5**: Hsuetal_dataset_FishUVC.csv

    Fish assemblages surveyed by Underwater Visual Census (UVC)

    **File 6**: Hsuetal_dataset_Obligate_corallivore.csv

    List of obligate corallivores as defined by [Cole et al. (2008)](https://doi.org/10.1111/j.1467-2979.2008.00290.x)

    **File 7**: Hsuetal_dataset_prey_MFG.csv

    List of coral prey to morpho-functional groups

### Data-specific information

##### **File 1: Hsuetal_dataset_benthic_MFG.csv**

1.  Number of variables/columns: 2

2.  Number of cases/rows: 73

3.  Missing data codes

    None

4.  Variable List

    -   Column A - taxa: benthic taxa
    -   Column B - MFG: morpho-functional group

5.  Abbreviations used

    -   for MFG:
        -   Abiotic: abiotic substrate
        -   Acro.arb: arborescent *Acropora*
        -   Acro.bran: other branching *Acropora*
        -   Acro.cory: corymbose *Acropora*
        -   Acro.tab: tabular *Acropora*
        -   alg.ca: calcareous algae
        -   alg.oth: other algae
        -   CCA: crustise coralline algae
        -   Foliose.lett: lettuce-like foliose hard coral
        -   Foliose.plat: plate-like foliose hard coral
        -   Freeliving: free living hard coral (e.g., Fugiidae)
        -   HC.other.bran: other branching hard coral
        -   HC.other.mas: other massive hard coral
        -   oth.ben: other benthos
        -   oth.mob: other mobile organisms
        -   Porites.mas: massive *Porites*
        -   SC.bran: branching soft coral
        -   SC.fan: fan-shaped soft coral
        -   SC.mas: massive soft coral
        -   SC.oth: other soft coral
        -   SC.whip: whip-like soft coral
        -   TurfL turf algae
        -   UnID: unidentifiable

##### **File 2: Hsuetal_dataset_benthic_taxa.csv**

1.  Number of variables/columns: 79

2.  Number of cases/rows: 65

3.  Missing data codes

    None

4.  Variable List

    -   Column A - Region
    -   Column B - Reef
    -   Column C - Habitat
    -   Column D - Position
    -   Column E - Site
    -   Column F - Plot
    -   Column G - Site_Plot
    -   Column H-CA: 72 benthic taxa

5.  Abbreviations used

    -   for Region:
        -   GBR: Great Barrier Reef
        -   TS: Torres Strait
    -   for Reef:
        -   CBHE: Heron Island, southern GBR
        -   OCDA: Davies Reef, central GBR
        -   ONLI: Lizard Island, northern GBR
        -   TSMA: Masig Island, TS
    -   for Habitat:
        -   BA: Back-reef
        -   FR: Fore-reef

##### **File 3: Hsuetal_dataset_Complexity.csv**

1.  Number of variables/columns: 6

2.  Number of cases/rows: 65

3.  Missing data codes

    None

4.  Variable List

    -   Column A - Reef
    -   Column B - Habitat
    -   Column C - Position
    -   Column D - Site
    -   Column E - Plot
    -   Column F - Rugosity_plane: the ratio of 3D surface area over the 2D surface area projected onto the plane of best fit

5.  Abbreviations used

    -   for Reef:
        -   CBHE: Heron Island, southern GBR
        -   OCDA: Davies Reef, central GBR
        -   ONLI: Lizard Island, northern GBR
        -   TSMA: Masig Island, TS
    -   for Habitat:
        -   BA: Back-reef
        -   FR: Fore-reef

##### **File 4: Hsuetal_dataset_Corallivory.csv**

1.  Number of variables/columns: 15

2.  Number of cases/rows: 6504

3.  Missing data codes

    None

4.  Variable List

    -   Column A - Reef
    -   Column B - Habitat
    -   Column C - Position
    -   Column D - Site
    -   Column E - Plot
    -   Column F - Time
    -   Column G - SN: Genus_species
    -   Column H - Family
    -   Column I - Genus
    -   Column J - Species
    -   Column K - Number
    -   Column L - Activity: either Passing or Feeding
    -   Column M - TL: total length (cm)
    -   Column N - Comment
    -   Column O - ID: ID to distinguish fish individulas for each species at each plot

5.  Abbreviations used

    -   for Reef:
        -   CBHE: Heron Island, southern GBR
        -   OCDA: Davies Reef, central GBR
        -   ONLI: Lizard Island, northern GBR
        -   TSMA: Masig Island, TS
    -   for Habitat:
        -   BA: Back-reef
        -   FR: Fore-reef
    -   for Comment:
        -   time in: the time each individual entered the field of view
        -   time out: the time each individual left the field of view

##### **File 5: Hsuetal_dataset_FishUVC.csv**

1.  Number of variables/columns: 13

2.  Number of cases/rows: 923

3.  Missing data codes: NA

4.  Variable List

    -   Column A - Reef
    -   Column B - Habitat
    -   Column C - Position
    -   Column D - Site
    -   Column E - Family
    -   Column F - Genus
    -   Column G - Species
    -   Column H - Valid_as: scientific name
    -   Column I - Trophic_group
    -   Column J - Feeding_group
    -   Column K - TR: transect
    -   Column L - Density (individuals per 100 m<sup>-2<sup>)
    -   Column M - Biomass (g per 100 m<sup>-2<sup>)

5.  Abbreviations used

    -   for Reef:
        -   CBHE: Heron Island, southern GBR
        -   OCDA: Davies Reef, central GBR
        -   ONLI: Lizard Island, northern GBR
        -   TSMA: Masig Island, TS
    -   for Habitat:
        -   BA: Back-reef
        -   FR: Fore-reef
    -   for TR:
        -   T1: Transect 1, spaning plots 1 & 2
        -   T2: Transect 2, spaning plots 3 & 4

##### **File 6: Hsuetal_dataset_Obligate_corallivore.csv**

1.  Number of variables/columns: 4

2.  Number of cases/rows: 30

3.  Missing data codes

    None

4.  Variable List

    -   Column A - Family
    -   Column B - Genus
    -   Column C - Species
    -   Column D - SN: Genus_species
    -   Column E - Vaild_as: scientific name

##### **File 7: Hsuetal_dataset_prey_MFG.csv**

1.  Number of variables/columns: 2

2.  Number of cases/rows: 29

3.  Missing data codes: NA

4.  Variable List

    -   Column A - Prey: coral preys form which corallivores took bites
    -   Column B - MFG: morpho-functional group

5.  Abbreviations used

    -   for MFG:
        -   Acropora.tab: Tabular *Acropora*
        -   Encrusting: Encrusting hard coral
        -   Acropora.Branching.Other
        -   Massive: Massive hard coral
