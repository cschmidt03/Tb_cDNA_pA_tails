#cython: language_level=3, boundscheck=False, wraparound=False
import cython

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.unordered_set cimport unordered_set

ctypedef (int, int, char*) hit_tuple

cdef int hamming_distance(char* a, char* b, Py_ssize_t L) :
    cdef int distance = 0
    cdef int i = 0
    while i < L:
        distance += a[i] != b[i]
        i += 1
    return distance

def hamming_scan(str seq, str pat, int mm_max):
    cdef vector[hit_tuple] hits
    cdef hit_tuple hit = (0, 0, "mm")

    cdef Py_ssize_t n = len(seq)
    cdef Py_ssize_t L = len(pat)
    cdef bytes b_seq = seq.encode('utf-8')
    cdef bytes b_pat = pat.encode('utf-8')
    cdef char* c_seq = b_seq
    cdef char* c_pat = b_pat

    cdef Py_ssize_t i = 0
    cdef size_t j = 0

    while i < n - L + 1:
        if hamming_distance(c_seq + i, c_pat, L) <= mm_max:
            hit[0] = i
            hit[1] = i+L
            hits.push_back(hit)
        i += 1

    cdef list return_hits = [None] * hits.size()

    while j < hits.size(): 
        return_hits[j] = hits[j]
        j += 1

    return return_hits

def deletion_scan_set(str seq, str pat):
    cdef string s_seq = seq.encode('utf-8')
    cdef string s_pat = pat.encode('utf-8')
    cdef string del_var
    cdef unordered_set[string] del_variants

    cdef size_t n = s_seq.size()
    cdef size_t L = s_pat.size()
    cdef size_t i = 0
    cdef size_t j = 0

    cdef char last_drop = '\0'
    while i < L:
        if s_pat[i] != last_drop: #if this char is the same as the last char removed, the variants are the same!
            del_variants.insert(s_pat.substr(0, i) + s_pat.substr(i+1))
        last_drop = s_pat[i]
        i += 1

    cdef vector[hit_tuple] hits
    cdef hit_tuple hit = (0, 0, "del")

    while j < n - L + 2:
        if del_variants.find(s_seq.substr(j, L-1)) != del_variants.end():
            hit[0] = j
            hit[1] = j + L - 1
            hits.push_back(hit)
        j += 1

    cdef list return_hits = [None] * hits.size()
    j = 0
    while j < hits.size(): 
        return_hits[j] = hits[j]
        j += 1
        
    return return_hits

def deletion_scan_vector(str seq, str pat):
    cdef string s_seq = seq.encode('utf-8')
    cdef string s_pat = pat.encode('utf-8')
    cdef string del_var

    cdef size_t n = s_seq.size()
    cdef size_t L = s_pat.size()
    cdef size_t i = 0
    cdef size_t j = 0

    cdef vector[string] del_variants
    del_variants.reserve(L)

    cdef char last_drop = '\0'
    while i < L:
        if s_pat[i] != last_drop: #if this char is the same as the last char removed, the variants are the same!
            del_variants.push_back(s_pat.substr(0, i) + s_pat.substr(i+1))
        last_drop = s_pat[i]
        i += 1

    cdef vector[hit_tuple] hits
    cdef hit_tuple hit = (0, 0, "del")

    cdef int x = 0
    for v in del_variants:
        j = 0
        while j < n - L + 2:
            if s_seq.compare(j, L-1, v) == 0:
                hit[0] = j
                hit[1] = j + L - 1
                hits.push_back(hit)
            j += 1
        x += 1

    cdef list return_hits = [None] * hits.size()
    j = 0
    while j < hits.size(): 
        return_hits[j] = hits[j]
        j += 1
        
    return return_hits
