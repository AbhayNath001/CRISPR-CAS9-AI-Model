# DNA Re-engineering

A research-oriented deep learning and bioinformatics project for DNA mutation analysis using CRISPR sequence datasets.

This repository includes:
- Sequence mutation comparison utilities
- Deep learning models for mutation detection/classification
- Statistical testing workflows
- Visualization scripts for CRISPR workflows and sequence clustering

## Repository Structure

- `main.py`: Baseline deep learning workflow using encoded DNA sequences, mutation detection summary, and Tanimoto similarity.
- `sequencing.py`: CNN + BiLSTM mutation classification workflow with per-sample visualization.
- `metrics.py`: Enhanced per-position mutation detection model with attention, weighted loss, ROC analysis, and model checkpointing.
- `metric.py`: McNemar test-based statistical comparison per disease category.
- `plot.py`: Hierarchical clustering and dendrogram generation from sequence k-mer profiles.
- `editing.py`: CRISPR-Cas9 workflow visualization including target, cut site, and HDR template rendering.
- `more.ipynb`: Notebook for additional experimentation.
- `CRISPR_Mutation_Dataset_Provisional_200.xlsx`: Primary dataset used by scripts.

## Requirements

- Python 3.12 or higher recommended
- biopython==1.87
- numpy==2.4.4
- pandas==3.0.2
- python-dateutil==2.9.0.post0
- six==1.17.0
- tzdata==2026.1
- tensorboard==2.16.2
- tensorboard-data-server==0.7.2
- tensorflow==2.16.1
- tensorflow-hub==0.16.1
- tensorflow-intel==2.16.1
- scikit-learn==1.8.0
- scikit-misc==0.5.2
- scikit-plot==0.3.7
- scipy==1.16.3

Project dependencies are pinned in `requirements.txt`.

## Installation

1. Clone the repository.

```bash
git clone https://github.com/AbhayNath001/CRISPR-CAS9-AI-Model
cd CRISPR-CAS9-AI-Model
pip install -r requirements.txt
```

2. Create and activate a virtual environment.

Windows (PowerShell):

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

Linux/macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

3. Install dependencies.

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

## Dataset Format

Most scripts expect the Excel file `CRISPR_Mutation_Dataset_Provisional_200.xlsx` with at least these columns:
- `Reference DNA`
- `Patient DNA`
- `Disease`

`metric.py` also expects `Sheet1` and assumes 27 unique diseases.

## How To Use

Run scripts from the repository root.

### 1) Baseline Mutation Analysis

```bash
python main.py
```

What it does:
- Loads reference/patient DNA pairs
- Finds mismatch positions
- Computes Tanimoto similarity
- Trains a simple Keras model and prints predicted mutation impact

### 2) Sequence Modeling + Mutation Visualization

```bash
python sequencing.py
```

What it does:
- Trains Conv1D + BiLSTM model
- Reports accuracy, precision, recall, F1-score
- Visualizes reference vs mutated sequences

## Notes and Troubleshooting

- If TensorFlow installation fails, verify Python version (3.12 usually works best for TensorFlow 2.16.x).
- If script errors mention missing columns, verify dataset column names exactly match expected names.

## Reproducibility Tips

- Keep dataset files unchanged from expected schema.
- Use the pinned versions from `requirements.txt`.
- For stable comparisons, keep random seed settings in scripts unchanged.

## 📬 Contact

<table>
<tr>
<td width="65%">

### 👤 Research Profile

**Name:** Abhay Nath |
**Affiliation:** AI-ML Research Scientist

**📧 Email:** [dr.abhaynath001@gmail.com](mailto:dr.abhaynath001@gmail.com)
**🔗 GitHub:** https://github.com/AbhayNath001
**🎓 LinkedIn:** https://www.linkedin.com/in/abhaynath001/ 

**🧪 Research Areas:**

* Artificial Intelligence
* Machine Learning
* Deep Learning
* Computer Vision
* Generative AI
* MedTech

</td>

<td width="35%" align="center">

<img src="https://raw.githubusercontent.com/AbhayNath001/AbhayNath001/main/Passport%20Size%20Image%20(1).jpg" width="160" style="border:2px solid #e0e0e0; border-radius:12px; padding:6px;" />

<br><br>

<i>Abhay Nath | Kolkata, West Bengal, India</i>

</td>
</tr>
</table>

---
