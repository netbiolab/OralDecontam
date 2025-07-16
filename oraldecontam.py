import os
import subprocess
import argparse
import shutil
import gzip


def check_kraken2_db(db_path, downloader):
    if os.path.exists(db_path) and os.path.isdir(db_path):
        print(f"Kraken2 DB found at {db_path}")
        return True
    else:
        print(f"Kraken2 DB not found at {db_path}. Attempting to download...")
        # Try wget, curl, aria2 in that order
        kraken_url = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz"
        archive_name = "k2_standard.tar.gz"

        download_cmds = {
            "wget": ["wget", "-O", archive_name, kraken_url],
            "curl": ["curl", "-L", "-o", archive_name, kraken_url],
            "aria2": ["aria2c", "-o", archive_name, kraken_url],
        }

    if downloader:
        if downloader not in download_cmds:
            raise ValueError(f"Unknown downloader: {downloader}")
        if shutil.which(download_cmds[downloader][0]) is None:
            raise RuntimeError(f"{downloader} is not installed.")
        print(f"Using specified downloader: {downloader}")
        subprocess.run(download_cmds[downloader], check=True)
    else:
        print("Downloader was not specified. Resorting to automatic selection.")
        for name, cmd in download_cmds.items():
            if shutil.which(cmd[0]) is None:
                continue
            try:
                print(f"Attempting download with {name}...")
                subprocess.run(cmd, check=True)
                print(f"Download successful with {name}.")
                return
            except subprocess.CalledProcessError:
                print(f"{name} failed. Trying next...")
        raise RuntimeError("Failed to download Kraken2 DB with wget, curl, or aria2.")

    subprocess.run(["mkdir", "-p", db_path], check=True)
    subprocess.run(["tar", "-xzvf", archive_name, "-C", db_path, "--strip-components=1"], check=True)
    print("Kraken2 DB download complete.")

def run_kraken2(forward_path, reverse_path, db_path, classified_prefix, unclassified_prefix):
    # Construct Kraken2 command for paired reads
    command = [
        "kraken2",
        "--db", db_path,
        "--report", "kraken2_report.txt",
        "--classified-out", f"{classified_prefix}#.fastq",
        "--unclassified-out", f"{unclassified_prefix}#.fastq",
        "--output", "kraken2_output.txt",
        "--paired",              # indicate paired-end mode
        forward_path,            # forward reads file
        reverse_path             # reverse reads file
    ]
    print("Running Kraken2 for classification (paired-end mode)...")
    subprocess.run(command, check=True)
    print("Classification complete. Outputs saved.")

def separate_reads(kraken_output_path, forward_path, reverse_path, output_prefix, output_dir=None):
    # Define taxonomy markers indicating human reads (GRCh38/Homo sapiens)
    human_markers = {'d__Eukaryota', 'p__Chordata', 'c__Mammalia', 
                     'o__Primates', 'f__Hominidae', 'g__Homo', 
                     's__Homosapiens_GrCH38.fna', 'GRCH38.fna'}
    non_human_ids = set()
    root_ids = set()
    
    # Parse the Kraken2 output file to populate non_human_ids and root_ids
    with open(kraken_output_path, 'r') as kfile:
        for line in kfile:
            line = line.strip()
            if not line:
                continue
            # Check for root classification in the line
            if 'root (taxid 1)' in line:
                # Mark this read ID as ambiguous (classified only to root)
                # Kraken2 output format: <status>\t<read-id>\t<taxid>\t...<classification info>
                # The read ID is the second column:
                root_ids.add(line.split('\t')[1])
            # Check if any human-specific taxonomy marker is present
            human_flag = False
            for marker in human_markers:
                if marker in line:
                    human_flag = True
                    break
            if not human_flag:
                # This read is not classified as human
                non_human_ids.add(line.split('\t')[1])
    print("Finished parsing Kraken2 output. Starting FASTQ separation...")
    
    # Determine output directory and base prefix for file naming
    out_dir = output_dir or "."  # current directory if not specified
    base = output_prefix
    # Make sure the directory exists (create if needed)
    os.makedirs(out_dir, exist_ok=True)
    
    # Open output files for writing (six files: human, non-human, ambiguous for R1 and R2)
    human_out_1 = open(os.path.join(out_dir, f"{base}.grch38_1.fastq"), 'w')
    human_out_2 = open(os.path.join(out_dir, f"{base}.grch38_2.fastq"), 'w')
    non_human_out_1 = open(os.path.join(out_dir, f"{base}.non-human_1.fastq"), 'w')
    non_human_out_2 = open(os.path.join(out_dir, f"{base}.non-human_2.fastq"), 'w')
    ambiguous_out_1 = open(os.path.join(out_dir, f"{base}.ambiguous_1.fastq"), 'w')
    ambiguous_out_2 = open(os.path.join(out_dir, f"{base}.ambiguous_2.fastq"), 'w')
    
    # Helper to open possibly gzipped files
    def open_fastq(path):
        return gzip.open(path, 'rt') if path.endswith(".gz") else open(path, 'r')
    
    # Function to process one FASTQ file (either forward or reverse)
    def split_fastq_file(input_path, is_forward=True):
        out_human = human_out_1 if is_forward else human_out_2
        out_non_human = non_human_out_1 if is_forward else non_human_out_2
        out_ambiguous = ambiguous_out_1 if is_forward else ambiguous_out_2
        with open_fastq(input_path) as infile:
            while True:
                # Read one FASTQ record (4 lines)
                header = infile.readline()
                if not header:
                    break  # EOF reached
                seq = infile.readline().rstrip("\n")
                plus = infile.readline().rstrip("\n")
                qual = infile.readline().rstrip("\n")
                # Extract the read ID from header (remove '@' and any trailing space or pair indicator)
                read_id = header[1:].strip().split()[0]
                # Determine category and write to appropriate files
                if read_id in non_human_ids:
                    out_non_human.write(f"{header}{seq}\n{plus}\n{qual}\n")
                else:
                    out_human.write(f"{header}{seq}\n{plus}\n{qual}\n")
                if read_id in root_ids:
                    out_ambiguous.write(f"{header}{seq}\n{plus}\n{qual}\n")
        # (No explicit close needed for infile due to context manager)
    
    # Process both forward and reverse FASTQ files
    split_fastq_file(forward_path, is_forward=True)
    split_fastq_file(reverse_path, is_forward=False)
    
    # Close all output files
    human_out_1.close(); human_out_2.close()
    non_human_out_1.close(); non_human_out_2.close()
    ambiguous_out_1.close(); ambiguous_out_2.close()
    print(f"Read separation complete. Outputs written to directory: {out_dir}")


def main():
    parser = argparse.ArgumentParser(description="Classify FASTQ reads with Kraken2 and split into categories.")
    parser.add_argument("--forward", required=True, help="Path to forward reads FASTQ (R1)")
    parser.add_argument("--reverse", required=True, help="Path to reverse reads FASTQ (R2)")
    parser.add_argument("--db", required=True, help="Kraken2 database path")
    parser.add_argument("--classified", default="classified", help="Prefix for classified reads output")
    parser.add_argument("--unclassified", default="unclassified", help="Prefix for unclassified reads output")
    parser.add_argument("--prefix", default=None, help="Prefix for separated output FASTQs (e.g., sample name)")
    parser.add_argument("--read_dir", default=None, help="Base directory for input FASTQ files")
    parser.add_argument("--downloader", default=None, help="Downloader to use for DB (wget, curl, aria2)")
    args = parser.parse_args()

    # Resolve input file paths with read_dir if provided
    forward_path = args.forward
    reverse_path = args.reverse
    if args.read_dir:
        forward_path = os.path.join(args.read_dir, args.forward)
        reverse_path = os.path.join(args.read_dir, args.reverse)
    
    # Determine output prefix for separated files
    if args.prefix:
        sample_prefix = args.prefix
    else:
        # Derive from forward filename (e.g., remove trailing _1 or R1 and extension)
        fname = os.path.basename(forward_path)
        if fname.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
            fname = fname.split(".fastq")[0].split(".fq")[0]  # remove extension
        # remove common R1 designators
        sample_prefix = fname.replace("_1", "").replace("_R1", "")
    
    # Ensure Kraken2 DB is present (download if not)
    check_kraken2_db(args.db, args.downloader)
    # Run Kraken2 classification on paired reads
    run_kraken2(forward_path, reverse_path, args.db, args.classified, args.unclassified)
    # Perform read separation into human/non-human/ambiguous categories
    separate_reads("kraken2_output.txt", forward_path, reverse_path, sample_prefix)


if __name__ == "__main__":
    main()

