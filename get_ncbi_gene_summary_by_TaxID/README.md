

# get_ncbi_gene_summary_by_TaxID

A robust, resume-capable Python tool to batch fetch gene summaries from NCBI Entrez for specific species (TaxIDs).

## 🚀 Why this script? (Common Pitfalls vs. This Solution)

Many bioinformatics scripts for fetching NCBI data are written as "one-off" snippets. They often fail when processing large datasets (like the Human genome with ~25k+ genes) due to network instability or API limitations.


| Common Issues in "Quick Scripts" | **How This Script Solves Them** |
| :-- | :-- |
| **Fragile** (Crashes on network timeout) | **Robust**: Implements automatic retries with exponential backoff for network errors. |
| **No Checkpoints** (Fail at 99% = Start over) | **Resume-Capable**: Checks for existing outputs. If you stop it, it resumes exactly where it left off. |
| **Corrupted Files** (Partial writes on crash) | **Atomic Writes**: Writes to `.tmp` first and only renames to `.csv.gz` upon 100% success. No half-written files. |
| **IP Bans** (Ignoring rate limits) | **Compliant**: strictly follows NCBI's 3 reqs/s (or 10 reqs/s with API Key) limits. |
| **URL Length Errors** (Requesting too many IDs) | **Smart Batching**: Automatically chunks requests (e.g., 300 genes/batch) to prevent HTTP 414 errors. |
| **Disk Space** | **Compressed**: Directly saves as `.csv.gz`, saving ~80% disk space compared to raw CSV. |

## ✨ Features

* **Batch Processing**: Reads a list of TaxIDs and processes them one by one.
* **Gene Discovery**: Automatically finds all valid `GeneID`s for a given TaxID (filtering out dead records).
* **Data Extraction**: Fetches Symbol, Description, Summary, and Aliases.
* **Pandas Integration**: Outputs structured data ready for analysis.
* **NCBI API Key Support**: Supports API keys for faster downloading (up to 10 requests/second).


## 🛠️ Installation

1. **Clone the repository**:

```bash
git clone https://github.com/houhaiyang/Tools.git
cd get_ncbi_gene_summary_by_TaxID
```

2. **Install dependencies**:
This script relies on `biopython` for API interaction and `pandas` for data handling.

```bash
pip install biopython pandas
```


## 📂 Usage

### 1. Prepare Input

Create a text file (default: `data/ncbi_summary/TaxID-list.txt`) containing one TaxID per line:

```text
9606
10090
7955
```


### 2. Run the Script

**Basic Usage:**

```bash
python get_ncbi_gene_summary_by_TaxID.py --email your_email@example.com
```

**With API Key (Recommended for speed):**
If you have an NCBI API Key, the script runs ~3x faster.

```bash
python get_ncbi_gene_summary_by_TaxID.py \
    --email your_email@example.com \
    --api_key your_api_key_string
```

**Custom Paths:**

```bash
python get_ncbi_gene_summary_by_TaxID.py \
    --input my_taxids.txt \
    --output_dir my_results/
```


### 3. Output Format

The script generates compressed CSV files in the output directory (e.g., `9606.csv.gz`).


| Column | Description |
| :-- | :-- |
| `GeneID` | NCBI Entrez Gene ID (Unique integer) |
| `Symbol` | Official Gene Symbol (e.g., TP53) |
| `Summary` | The curated gene function summary |
| `Description` | Full gene name |
| `OtherDesignations` | Aliases and other names |

## 📁 Directory Structure

```text
.
├── get_ncbi_gene_summary_by_TaxID.py
├── README.md
└── data
    └── ncbi_summary
        ├── TaxID-list.txt          # Input file
        └── species                 # Output directory
            ├── 9606.csv.gz         # Human
            ├── 10090.csv.gz        # Mouse
            └── ...
```


## ⚠️ Note on Rate Limits

* **Without API Key**: Max 3 requests/second.
* **With API Key**: Max 10 requests/second.
* The script automatically adjusts delays (`sleep`) to stay safe, but using an API key is highly recommended for large datasets (like Human or Mouse).


## 👨‍💻 Author

**Hou Haiyang**

* **Email**: [houhaiyang1@genomics.cn](mailto:houhaiyang1@genomics.cn)
* **Affiliation**: UCAS AI \& BGI Research AI Lab, Beijing
* **Date**: 2025-12-18


## 📄 License

This project is open source and available under the [MIT License](LICENSE).

