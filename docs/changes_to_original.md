## Changes to original repo (as of 20.07.2026)

### 1) Handling of renamed columns in output of dorado summary (dorado v2.0.x)
Several columns of the dorado summary outputs have been renamed between v0.9.x and v2.0.x, 
causing errors in `filter_join_final_snake.R`:
* `alignment_mapq` to `alignment_mapping_quality`
* `barcode` to `alias`, this now doesn't include the kit prefix!
* `alignment_length` was removed

The names were kept the same in the output of the R script for compatibility.

### 2) Performance improvements for the adapter filtering step
#### 2.1) Cython implementations for hamming and deletion scanning
Profiling the code revealed that the chokepoint of the script is the hamming scanning.  
The original `ham(a, b)` used `zip()` to calculate hamming distance, writing the procedure as a for loop already increases speed ~ x2.  
Far better improvement is achieved by implementing the function in Cython, which allows to take adavantage of C/C++'s speed while integrating seamlessly to Python code.  
The behaviour of `cy_stringscan.hamming_scan()` is almost the same as the original functionality, with two changes:
* Tuples for storing hits now have a `bytes` object at index 2 instead of the `string`, i.e.: `(3, 23, b'mm')`.
* Tuples are not appended one by one to the `hits` list, instead the function internally appends to a `std::vector` and later gives the references to a Python `list` preallocated to the correct size, avoiding the costly `append()` operation. The list is the return value.  

A lesser improvement is achieved by reimplementing the deletion scanning in Cython as well.
While the `set` lookup should be faster ($O(1)$) than iterating through all drop variants ($O(m)$, with $m$ = length of adapter sequence), it appears that the iteration's low overhead outcompetes the set lookup for small enough $m$. An (unused) implementation of the original `set` lookup is given using C++'s `std::unordered_map`, should longer adapter sequences be used at some point. Returning the hits is handled as described above for the hamming scan.

#### 2.2) Multiprocessing implementation for the script
More speedup is achieved by implementing the script to make use of more threads. The general layout consits of three stages:
* **Reader Process**:  
Reads a BAM file and gets **virtual offsets** of the reads' beginnings. The offsets are put in a queue shared with the worker processes. Since writing and reding from the queue is time-consuming, the offsets are sent for batches of reads, the number of reads in a batch is defined by parameter `CHUNK_SIZE` (currently set to 2000).  
An (unused) implementation in pure Python using `pysam` is given. However, it is not fast enough to supply 8 or more worker processes.
A functionally equivalent, but faster implementation is written in Cython, directly referencing `htslib`, using a lower level function to avoid unnecessary parsing of the read since only its offset in the file is needed.  

    A note regarding `ssize_t bgzf_read(BGZF* fp, void* data, size_t length)`:  
    Reading from the BAM file requires a memory buffer where the read bytes are written to. This buffer is currently set to 4096 kiB. If the read is detected to be longer than the buffer, the program will raise a `RuntimeError` to prevent exiting by segmentation fault. The longest read length observed was ~ 1500 kB, the usual size is ~ 2 kB. These 'freak reads' seem to originate from pores outputting garbage signal, as the reads consist of repetitive patterns. Occurences were > 128 kiB every 200k reads, > 256 kiB every 500k reads. 


* **Worker Processes**:  
Write accepted reads in their own particular output BAM files, which are afterwards merged with `samtools cat`. The extracted reads for the TSV file are put onto a writer queue as they are small. Writing to the queue happens when a new batch of reads is fetches, i.e. after `CHUNK_SIZE` reads.

* **Writer Process**:  
Fetches read extractions from the output queue and writes to a TSV file.

The number of processes is set by the `--n-threads <N>` command line option. The script will start one reader, one writer and N-2 worker threads.

**Important Note:** In contrast to the single threaded script, the multiprocessing version **does not** guarantee that the reads will be written in the same order that they were read, thus it cannot be used for sorted or paired-read BAM files.

Measured performance on Intel Xeon W5 (on 640k reads, keep in mind that e.g. 4 threads -> 2 workers):  
|  | single process <br>(original) | single process <br>(Cython) | 4 threads | 6 threads | 8 threads | 10 threads | 12 threads |   
| - | --- | --- | --- | --- | --- | --- | --- |
| 10^3 reads/s | 2.2 | 12.3 | 21.3 | 35.5 | 45.4 | 52 | 55 |
| rel. | | x1 | x1.7 | x2.9 | x3.6 | x4.2 | x4.5 |


