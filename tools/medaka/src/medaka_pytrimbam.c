#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "kvec.h"
#include "medaka_common.h"
#include "medaka_trimbam.h"
#include "medaka_pytrimbam.h"


seq_set PY_retrieve_trimmed_reads(
    const char *region, const char *bam_file, size_t num_dtypes, char *dtypes[],
    const char tag_name[2], const int tag_value, const bool keep_missing, 
    const bool partial, const char *read_group) {

    trimmed_reads t_reads = retrieve_trimmed_reads(
        region, bam_file, num_dtypes, dtypes,
        tag_name, tag_value, keep_missing, partial, read_group);
 
    seq_set reads = xalloc(1, sizeof(*reads), "seq_set");
    reads->n_seqs = t_reads->sequences.n;
    reads->seqs = t_reads->sequences.a;
    reads->is_rev = t_reads->is_reverse.a;
    reads->t_reads = (void*) t_reads;
    return reads;
}


/** Cleans up the result of PY_retrieve_trimmed_reads
 *
 * @param seqs as seq_set structure.
 *
 * @returns void
 *
 */
void PY_destroy_reads(seq_set reads) {
    trimmed_reads t_reads = (trimmed_reads) reads->t_reads;
    destroy_trimmed_reads(t_reads);
    free(reads);
}
