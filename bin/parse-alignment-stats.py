#!/usr/bin/env python3
import csv
import argparse



def parse_samtools_stats_file(stats_file):
    """Parse a single samtools stats file and return key stats including total_bases."""
    total_reads = mapped_reads = mapped_bases = total_bases = None

    with open(stats_file) as f:
        for line in f:
            if line.startswith("SN"):
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                key, value = parts[1].strip(":"), parts[2]
                if key == "sequences":
                    total_reads = int(value)
                elif key == "reads mapped":
                    mapped_reads = int(value)
                elif key == "bases mapped (cigar)":
                    mapped_bases = int(value)
                elif key == "total length":
                    total_bases = int(value)

    if total_reads is None or mapped_reads is None or mapped_bases is None or total_bases is None:
        raise ValueError(f"Could not parse required fields from {stats_file}")

    return total_reads, mapped_reads, mapped_bases, total_bases

def main():
    parser = argparse.ArgumentParser(description="Summarize samtools stats files into CSV")
    parser.add_argument("-statsfile", help="Single samtools stats file to parse")
    parser.add_argument("-g", "--genome_size", type=float, default=3.1e9,
                        help="Haploid genome size (default=3.1e9 bp for human)")
    parser.add_argument("-sample_id", help="Sample ID for single stats file")
    
    args = parser.parse_args()

    if args.statsfile:
        if not args.sample_id:
            raise ValueError("When using -statsfile, you must also provide -sample_id")

        
        total_reads, mapped_reads, mapped_bases, total_bases = parse_samtools_stats_file(args.statsfile)
        mapping_rate = mapped_reads / total_reads
        mean_coverage = mapped_bases / args.genome_size
        theoretical_coverage = total_bases / args.genome_size

        

        # 4. Write data to CSV file
    fname = f"{args.sample_id}.aligned.bam.stats.csv"
    with open(fname, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        
        # Define the header row (updated to include the 12X proportion)
        header = [
            "sample_id", 
            
            "total_reads", 
            "total_bases", 
            "mapped_reads", 
            "total_bases_mapped", # mapped_bases variable
            "mapping_rate", 
            "mean_coverage", 
            "theoretical_coverage"
        ]

        writer.writerow(header)
        # NEW LINE ADDED: Writing the actual data row
        data_row = [
            args.sample_id, 
       
            total_reads, 
            total_bases, 
            mapped_reads, 
            mapped_bases, 
            mapping_rate, 
            mean_coverage, 
            theoretical_coverage,
           
        ]
        writer.writerow(data_row)
    print(f"Wrote summary to {fname}")



if __name__ == "__main__":
    main()
