"""
Author: houhaiyang
Email: houhaiyang1@genomics.cn
Affiliation: UCAS AI & BGI Research AI Lab, Beijing
Date: 2025-12-18

Description:
This script fetches gene summaries from the NCBI database for a list of TaxIDs.
It handles batch processing, rate limiting, and saves the output as compressed CSV files (.csv.gz).
Features:
- Resume from interruption (checks for existing .csv.gz files)
- Atomic writes (writes to .tmp.csv.gz then renames)
- Robust GeneID extraction
"""

import os
import time
import argparse
import pandas as pd
from Bio import Entrez
from urllib.error import URLError, HTTPError
import socket

# ---------------- 配置区域 ----------------
DEFAULT_INPUT_FILE = "data/ncbi_summary/TaxID-list.txt"
DEFAULT_OUTPUT_DIR = "data/ncbi_summary/species/"
BATCH_SIZE = 300
REQUEST_DELAY = 0.34
RETRIES = 3

def setup_args():
    parser = argparse.ArgumentParser(description="Fetch NCBI Gene Summaries by TaxID")
    parser.add_argument("--input", default=DEFAULT_INPUT_FILE, help="Path to TaxID list file")
    parser.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR, help="Directory to save CSV.GZ files")
    parser.add_argument("--email", default="houhaiyang1@genomics.cn", help="Email required by NCBI Entrez")
    parser.add_argument("--api_key", default=None, help="NCBI API Key (optional, allows faster requests)")
    return parser.parse_args()

def safe_request(func, *args, **kwargs):
    for attempt in range(RETRIES):
        try:
            return func(*args, **kwargs)
        except (URLError, HTTPError, socket.timeout) as e:
            print(f"    [Warning] Network error: {e}. Retrying ({attempt + 1}/{RETRIES})...")
            time.sleep(2 * (attempt + 1))
        except Exception as e:
            print(f"    [Error] Unexpected error: {e}")
            break
    return None

def get_all_gene_ids(taxid):
    print(f"  Fetching Gene IDs for TaxID {taxid}...")
    search_term = f"txid{taxid}[Organism] AND alive[prop]"
    handle = safe_request(Entrez.esearch, db="gene", term=search_term, retmax=1)
    if not handle: return []

    try:
        record = Entrez.read(handle)
        count = int(record["Count"])
        print(f"    Found {count} genes.")
    except Exception as e:
        print(f"    [Error] Failed to parse search result: {e}")
        return []

    all_ids = []
    chunk_size = 10000
    for start in range(0, count, chunk_size):
        handle = safe_request(Entrez.esearch, db="gene", term=search_term, retstart=start, retmax=chunk_size)
        if handle:
            try:
                record = Entrez.read(handle)
                all_ids.extend(record["IdList"])
            except Exception:
                print("    [Error] Failed to read ID batch.")
        time.sleep(REQUEST_DELAY)
    return all_ids

def get_gene_summaries_batch(gene_ids):
    if not gene_ids: return []

    id_str = ",".join(gene_ids)
    handle = safe_request(Entrez.esummary, db="gene", id=id_str)

    results = []
    if handle:
        try:
            record = Entrez.read(handle)
            if 'DocumentSummarySet' in record and 'DocumentSummary' in record['DocumentSummarySet']:
                summaries = record['DocumentSummarySet']['DocumentSummary']
            else:
                summaries = record

            for gene in summaries:
                # ---------------- ID 提取核心修复 ----------------
                # 1. 尝试直接获取字典键 'Id' (最常见情况)
                g_id = gene.get("Id", "")

                # 2. 如果失败，尝试从 XML 属性获取 (Biopython 解析器特性)
                if not g_id and hasattr(gene, "attributes"):
                    g_id = gene.attributes.get("uid", "")

                # 3. 强制转为字符串并去除空格
                g_id = str(g_id).strip()

                # 4. 如果仍为空，标记为 Unknown，避免错位
                if not g_id:
                    g_id = "Unknown"
                # -----------------------------------------------

                item = {
                    "GeneID": g_id,
                    "Symbol": gene.get("Name", ""),
                    "Description": gene.get("Description", ""),
                    "Summary": gene.get("Summary", ""),
                    "OtherDesignations": gene.get("OtherDesignations", "")
                }
                results.append(item)
        except Exception as e:
            print(f"    [Error] Failed to parse summary batch: {e}")

    return results

def process_species(taxid, output_dir):
    taxid = str(taxid).strip()
    if not taxid: return

    final_output_path = os.path.join(output_dir, f"{taxid}.csv.gz")

    if os.path.exists(final_output_path):
        print(f"[Skip] TaxID {taxid}: Output file already exists.")
        return

    print(f"[Processing] TaxID {taxid}...")

    gene_ids = get_all_gene_ids(taxid)
    if not gene_ids:
        print(f"  [Warning] No genes found for {taxid}.")
        return

    all_data = []
    total_genes = len(gene_ids)
    print(f"  Downloading summaries for {total_genes} genes...")

    for i in range(0, total_genes, BATCH_SIZE):
        batch_ids = gene_ids[i : i + BATCH_SIZE]
        batch_data = get_gene_summaries_batch(batch_ids)
        all_data.extend(batch_data)

        percent = min(100, int((i + BATCH_SIZE) / total_genes * 100))
        print(f"    Progress: {percent}% ({len(all_data)}/{total_genes})", end="\r")
        time.sleep(REQUEST_DELAY)

    print(f"\n  Processing complete. Saving data...")

    if all_data:
        df = pd.DataFrame(all_data)

        # 确保 GeneID 列存在，如果是 NaN 则填 "Unknown"
        if "GeneID" in df.columns:
            df["GeneID"] = df["GeneID"].fillna("Unknown").astype(str)

        # 将其他列的 NaN 填充为空字符串
        df.fillna("", inplace=True)

        # 字段排序
        df = df[["GeneID", "Symbol", "Summary", "Description", "OtherDesignations"]]

        # ---------------- 文件保存逻辑修改 ----------------
        # 按照要求修改为：xxx.tmp.csv.gz
        temp_filename = f"{taxid}.tmp.csv.gz"
        temp_path = os.path.join(output_dir, temp_filename)

        try:
            # compression='gzip' 会自动处理 .gz 后缀
            df.to_csv(temp_path, index=False, encoding='utf-8', compression='gzip')

            # 原子操作：重命名为最终文件 xxx.csv.gz
            os.rename(temp_path, final_output_path)
            print(f"  [Success] Saved to {final_output_path}")

        except Exception as e:
            print(f"  [Error] Failed to save file: {e}")
            if os.path.exists(temp_path):
                os.remove(temp_path)
        # -----------------------------------------------
    else:
        print(f"  [Warning] No summary data retrieved for {taxid}.")

def main():
    args = setup_args()
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
        global REQUEST_DELAY
        REQUEST_DELAY = 0.11

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        return

    with open(args.input, 'r') as f:
        taxids = [line.strip() for line in f if line.strip()]

    print(f"Loaded {len(taxids)} TaxIDs from {args.input}")
    print("-" * 50)

    for taxid in taxids:
        try:
            process_species(taxid, args.output_dir)
        except KeyboardInterrupt:
            print("\n[Aborted] Script stopped by user.")
            break
        except Exception as e:
            print(f"[Error] Critical failure for TaxID {taxid}: {e}")
        print("-" * 50)

if __name__ == "__main__":
    main()
