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
- Atomic writes (writes to .tmp then renames)
- Automatic rate limiting for NCBI Entrez API
"""

import os
import time
import argparse
import pandas as pd
from Bio import Entrez
from urllib.error import URLError, HTTPError
import socket

# ---------------- 配置区域 ----------------
# 默认输入输出路径
DEFAULT_INPUT_FILE = "data/ncbi_summary/TaxID-list.txt"
DEFAULT_OUTPUT_DIR = "data/ncbi_summary/species/"

# NCBI API 限制配置
# BATCH_SIZE: 每次 esummary 请求的基因数量（建议 200-500）
BATCH_SIZE = 300
# REQUEST_DELAY: 基础延迟，保证无 Key 时不超过 3次/秒
REQUEST_DELAY = 0.34
RETRIES = 3  # 网络错误重试次数


def setup_args():
    parser = argparse.ArgumentParser(description="Fetch NCBI Gene Summaries by TaxID (Saved as .csv.gz)")
    parser.add_argument("--input", default=DEFAULT_INPUT_FILE, help="Path to TaxID list file")
    parser.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR, help="Directory to save CSV.GZ files")
    parser.add_argument("--email", default="houhaiyang1@genomics.cn", help="Email required by NCBI Entrez")
    parser.add_argument("--api_key", default=None, help="NCBI API Key (optional, allows faster requests)")
    return parser.parse_args()


def safe_request(func, *args, **kwargs):
    """带重试机制的 API 请求封装"""
    for attempt in range(RETRIES):
        try:
            return func(*args, **kwargs)
        except (URLError, HTTPError, socket.timeout) as e:
            print(f"    [Warning] Network error: {e}. Retrying ({attempt + 1}/{RETRIES})...")
            time.sleep(2 * (attempt + 1))  # 指数退避
        except Exception as e:
            print(f"    [Error] Unexpected error: {e}")
            break
    return None


def get_all_gene_ids(taxid):
    """获取指定 TaxID 下的所有 Gene ID"""
    print(f"  Fetching Gene IDs for TaxID {taxid}...")

    # 1. 获取总数
    # alive[prop] 排除已撤销的基因
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

    # 2. 分批获取所有 ID
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
    """批量获取基因摘要信息"""
    if not gene_ids:
        return []

    id_str = ",".join(gene_ids)
    handle = safe_request(Entrez.esummary, db="gene", id=id_str)

    results = []
    if handle:
        try:
            # esummary 返回的是 document summaries 列表
            record = Entrez.read(handle)
            # 处理 DocumentSummarySet 结构兼容性
            if 'DocumentSummarySet' in record and 'DocumentSummary' in record['DocumentSummarySet']:
                summaries = record['DocumentSummarySet']['DocumentSummary']
            else:
                summaries = record

            for gene in summaries:
                # 提取关键字段
                item = {
                    "GeneID": gene.get("Id", ""),
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
    """处理单个物种的主逻辑"""
    taxid = str(taxid).strip()
    if not taxid: return

    # 修改：文件名后缀改为 .csv.gz
    final_output_path = os.path.join(output_dir, f"{taxid}.csv.gz")

    # 增量检查：如果文件已存在，跳过
    if os.path.exists(final_output_path):
        print(f"[Skip] TaxID {taxid}: Output file already exists ({os.path.basename(final_output_path)}).")
        return

    print(f"[Processing] TaxID {taxid}...")

    # 1. 获取 ID 列表
    gene_ids = get_all_gene_ids(taxid)
    if not gene_ids:
        print(f"  [Warning] No genes found or error fetching IDs for {taxid}.")
        return

    # 2. 分批获取摘要
    all_data = []
    total_genes = len(gene_ids)
    print(f"  Downloading summaries for {total_genes} genes...")

    for i in range(0, total_genes, BATCH_SIZE):
        batch_ids = gene_ids[i: i + BATCH_SIZE]
        batch_data = get_gene_summaries_batch(batch_ids)
        all_data.extend(batch_data)

        # 简单的进度显示
        percent = min(100, int((i + BATCH_SIZE) / total_genes * 100))
        print(f"    Progress: {percent}% ({len(all_data)}/{total_genes})", end="\r")

        time.sleep(REQUEST_DELAY)

    print(f"\n  Processing complete. Saving data...")

    # 3. 保存数据（原子写入：先写临时文件，再重命名）
    if all_data:
        df = pd.DataFrame(all_data)
        # 字段筛选
        df = df[["GeneID", "Symbol", "Summary", "Description", "OtherDesignations"]]

        temp_path = final_output_path + ".tmp"
        try:
            # 修改：添加 compression='gzip' 参数
            df.to_csv(temp_path, index=False, encoding='utf-8', compression='gzip')
            os.rename(temp_path, final_output_path)
            print(f"  [Success] Saved to {final_output_path}")
        except Exception as e:
            print(f"  [Error] Failed to save file: {e}")
            if os.path.exists(temp_path):
                os.remove(temp_path)
    else:
        print(f"  [Warning] No summary data retrieved for {taxid}.")


def main():
    args = setup_args()

    # 设置 Entrez 全局参数
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
        global REQUEST_DELAY
        REQUEST_DELAY = 0.11  # 如果有 Key，可以加快速度 (约9次/秒)

    # 确保输出目录存在
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        print(f"Created output directory: {args.output_dir}")

    # 读取输入文件
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        return

    with open(args.input, 'r') as f:
        taxids = [line.strip() for line in f if line.strip()]

    print(f"Loaded {len(taxids)} TaxIDs from {args.input}")
    print("-" * 50)

    # 逐个处理
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
