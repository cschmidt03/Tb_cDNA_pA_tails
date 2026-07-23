#cython: language_level=3, boundscheck=False, wraparound=False
import cython 
from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libc.stdio cimport stderr, fprintf, fflush
from libcpp.vector cimport vector
from multiprocessing import Queue

cdef extern from "htslib/bgzf.h":
    ctypedef struct BGZF

    cdef BGZF* bgzf_open(const char* path, const char* mode)
    cdef ssize_t bgzf_read(BGZF* fp, void* data, ssize_t length)
    cdef uint64_t cy_bgzf_tell(BGZF* fp)
    
cdef extern from *:
    '''
    #include <htslib/bgzf.h>
    static inline uint64_t cy_bgzf_tell(BGZF* fp) {return bgzf_tell(fp);}
    '''

DEF DISCARD_BUFFER_SIZE = 4 * 1024 * 1024 #4 MiB, will get allocated in .bss segment not stack!
cdef uint8_t[DISCARD_BUFFER_SIZE] discard_buffer
cdef uint32_t size
ctypedef (uint64_t, uint64_t) offset_pair

def get_virtual_offsets(bytes path, int chunk_size, input_queue):
    cdef BGZF* bamf = bgzf_open(path, b"r")
    cdef uint32_t n_ref
    cdef vector[uint64_t] offsets
    cdef offset_pair send_offsets = (0,0)

    bgzf_read(bamf, discard_buffer, 4) # reads "BAM\1" 
    bgzf_read(bamf, &size, 4) # reads header_length
    bgzf_read(bamf, discard_buffer, size) # skips header_length bytes
    bgzf_read(bamf, &n_ref, 4) #reads n_refsequences
    
    cdef uint32_t idx = 0
    while idx < n_ref:
        bgzf_read(bamf, &size, 4) # reads length_ref_sequence_naem
        bgzf_read(bamf, discard_buffer, size) # skips length_ref_sequence_naem bytes
        bgzf_read(bamf, discard_buffer, 4) # skips 4 bytes (len_ref_sequence)
        idx += 1
    #bamf should now point to the beginning of the alignments!

    cdef int read_n = 0
    send_offsets[0] = cy_bgzf_tell(bamf)
    send_offsets[1] = cy_bgzf_tell(bamf)

    #bzgf_read returns 0 on EOF, -1 on ERROR, this writes the size of the next read to size
    while bgzf_read(bamf, &size, 4) > 0: 
        if read_n % chunk_size == 0 and read_n > 0:
            send_offsets[1] = cy_bgzf_tell(bamf) - 4
            input_queue.put(send_offsets)
            send_offsets[0] = send_offsets[1]
        if size > DISCARD_BUFFER_SIZE:
            raise RuntimeError(
                f"Read size {size} exceeds discard buffer size {DISCARD_BUFFER_SIZE}")
        bgzf_read(bamf, discard_buffer, size)
        read_n += 1

    print(read_n, " reads read, sending the last offset!")
    send_offsets[0] = send_offsets[1]
    send_offsets[1] = cy_bgzf_tell(bamf)
    input_queue.put(send_offsets)
