import multiprocessing
import gzip
from typing import List
import argparse
import xopen

def process_chunk(reads1: List[List[str]], reads2: List[List[str]], output_file1: List[str], lock: multiprocessing.Lock, *args, **kwargs) -> None:
    """
    Process a chunk of lines and write it to the output files.

    Args:
        lines1 (List[List[str]]): List of lines from the first input file.
        lines2 (List[List[str]]): List of lines from the second input file.
        output_file1 (str): Path to the first output file.
        output_file2 (str): Path to the second output file.
        lock (multiprocessing.Lock): Lock object for ensuring thread safety during file writes.
        *args: Additional positional arguments.
        **kwargs: Additional keyword arguments.

    Returns:
        None
    """
    for read1,read2 in zip(reads1, reads2):

    chunk1 = ''.join(lines1)
    chunk2 = ''.join(lines2)
    with lock:
        with open(output_file1, 'ab') as f1, open(output_file2, 'ab') as f2:
            f1.write(chunk1.encode())
            f2.write(chunk2.encode())

def parse_file_paths(file_paths_str: str) -> List[str]:
    """
    Parse input or output file paths separated by commas.

    Args:
        file_paths_str (str): String containing file paths separated by commas.

    Returns:
        List[str]: List of file paths.
    """
    return file_paths_str.split(',')

def process_compressed_files(input_files: List[str], output_files: List[str], barcode_dict_per_round, chunk_size: int = 4, *args, **kwargs) -> None:
    """
    Process two compressed files in parallel, reading line by line and writing chunks to corresponding output files.

    Args:
        input_files (List[str]): List of paths to input compressed files.
        output_files (List[str]): List of paths to output compressed files.
        chunk_size (int, optional): Number of lines to process as a chunk. Defaults to 4.
        *args: Additional positional arguments.
        **kwargs: Additional keyword arguments.

    Returns:
        None
    """
    lock = multiprocessing.Lock()
    process_pool = multiprocessing.Pool(len(input_files))

    with xopen.xopen(input_files[0], 'rt') as f1, xopen.xopen(input_files[1], 'rt') as f2:
        lines1 = []
        lines2 = []
        for line1, line2 in zip(f1, f2):
            lines1.append([next(f1) for _ in range(4)])
            lines2.append([next(f2) for _ in range(4)])
            if len(lines1) == chunk_size and len(lines2) == chunk_size:
                process_pool.apply_async(process_chunk, args=(lines1, lines2, output_files, lock) + args, kwds=kwargs)
                lines1 = []
                lines2 = []
        if lines1 and lines2:
            process_pool.apply_async(process_chunk, args=(lines1, lines2, output_files, lock) + args, kwds=kwargs)

    process_pool.close()
    process_pool.join()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process compressed files in parallel.')
    parser.add_argument('--input-files', required=True, help='Paths to input compressed files separated by commas')
    parser.add_argument('--output-files', required=True, help='Paths to output compressed files separated by commas')
    parser.add_argument('--exclusion-list', required=True, help='Path to the inclusion list file.')
    parser.add_argument('--chunk-size', type=int, default=10000, help='Number of lines to process as a chunk')
    parser.add_argument('-p', type=int, default=1, help='Number of cpus to use for processing the files')

    args = parser.parse_args()

    input_files = parse_file_paths(args.input_files)
    output_files = parse_file_paths(args.output_files)
    exclusion_list = args.exclusion_list
    chunk_size = args.chunk_size
    cpus = args.p

    barcode_dict_per_round = create_barcode_dicts_from_file(exclusion_list)
    process_compressed_files(input_files, output_files, barcode_dict_per_round, chunk_size)
