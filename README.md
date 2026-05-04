# kidney-bipartite-matching

A comparative empirical analysis of bipartite graph matching algorithms applied to kidney transplant donor-recipient allocation. Implements and evaluates three algorithms of increasing sophistication — Greedy Matching, Kuhn's Algorithm, and the Jonker-Volgenant Algorithm — on randomly generated compatibility graphs modeled after US/European kidney transplant demographics.

**Course project — Data Structures and Algorithms, Rutgers – New Brunswick**  
*Advaith "Ad" Subramanian Sahasranamam*

---

## Overview

Kidney transplant allocation is modeled as a bipartite graph G = (U, V, E) where:
- **U** — donors
- **V** — recipients
- **E** — compatible donor-recipient pairs

Compatibility is determined by four medical constraints:
1. ABO blood type compatibility
2. Rhesus (Rh) factor compatibility
3. HLA tissue typing across three loci (HLA-A, HLA-B, HLA-DR)
4. Crossmatch probability derived from recipient PRA score

---

## Algorithms

| Algorithm | Objective |
|---|---|
| Greedy Matching | Baseline — fast, locally optimal |
| Kuhn's Algorithm | Maximum cardinality matching |
| Jonker-Volgenant | Minimum cost optimal assignment |

---

## Project Structure

```
kidney-bipartite-matching/
│
├── graph/
│   ├── attributes.py        # Blood type, Rh, HLA, PRA generation
│   ├── compatibility.py     # Compatibility checks and cost functions
│   ├── builder.py           # Graph construction + matching algorithms
│   └── experiment.py        # All three experiments + plotting
│
└── README.md
```

---

## Experiments

Three experiments are conducted:

**1. Scalability** — graph size n varies from 100 to 500 (step 5), fixed rej_chance = 0.2, 3 trials per size.

**2. Donor : Recipient Ratio** — recipients fixed at 500, donor ratio varies from 0.1 to 0.9, 5 trials per ratio.

**3. Rejection Threshold** — graph size fixed at 500, rej_chance varies from 0.05 to 0.5, 5 trials per threshold.

Each experiment measures:
- Runtime (wall-clock)
- Average matching size
- Average HLA score (0–6)
- Average PRA score (0–100)
- Match percentage (ratio and threshold experiments)

---

## Installation

```bash
git clone https://github.com/<AdSahas>/kidney-bipartite-matching.git
cd kidney-bipartite-matching
pip install -r requirements.txt
```

---

## Usage

Run all three experiments in sequence:

```bash
cd graph
python experiment.py
```

Results are saved as CSVs in the working directory. Plots are saved as PNGs.

---

## Dependencies

- Python 3.14.3
- `numpy`
- `pandas`
- `matplotlib`
- `timeit`

---

## Key Results

All results and findings are documented in the report attached to this repository. 

## References

[1] UNOS, "Organ Procurement and Transplantation Network (OPTN) Data," U.S. Department of Health & Human Services, 2023. [Online]. Available: https://optn.transplant.hrsa.gov/data/

[2] University of Rochester Medical Center, "Organ Transplant Rejection," Life Sciences Learning Center, 2009. [Online]. Available: https://www.urmc.rochester.edu/MediaLibraries/URMCMedia/life-sciences-learning-center/documents/STUDENTRejection7-23-09.pdf

[3] A. E. Roth, T. Sonmez, and M. U. Unver, "Kidney exchange," Quarterly Journal of Economics, vol. 119, no. 2, pp. 457–488, 2004.

[4] C. Süsal and G. Opelz, "Current role of human leukocyte antigen matching in kidney transplantation," Current Opinion in Organ Transplantation, vol. 18, no. 4, pp. 438–444, 2013.

[5] UNOS, "How we match organs," 2023. [Online]. Available: https://unos.org/transplant/how-we-match-organs/

[6] S. A. Zenios, "Optimal control of a paired-kidney exchange program," Management Science, vol. 48, no. 3, pp. 328–342, 2002.

[7] H. W. Kuhn, "The Hungarian method for the assignment problem," Naval Research Logistics Quarterly, vol. 2, no. 1–2, pp. 83–97, 1955.

[8] R. Jonker and A. Volgenant, "A shortest augmenting path algorithm for dense and sparse linear assignment problems," Computing, vol. 38, no. 4, pp. 325–340, 1987.

[9] G. Garratty, S. A. Glynn, and R. McEntire, "ABO and Rh(D) phenotype frequencies of different racial/ethnic groups in the United States," Transfusion, vol. 44, no. 5, pp. 703–706, 2004.
