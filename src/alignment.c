/*
 * =============================================================================
 * FILE: alignment.c
 * MEMBER: 4 — Needleman-Wunsch Sequence Alignment
 * =============================================================================
 *
 * PURPOSE:
 *   Implements global sequence alignment using the Needleman-Wunsch algorithm.
 *   This module is called AFTER the search phase (suffix tree / trie search)
 *   when an exact match is NOT found. It computes a similarity score and
 *   prints the aligned sequences showing matches, mismatches, and gaps.
 *
 * HOW IT FITS IN THE PROJECT FLOW:
 *   1. Load DNA sequences (utils.c / trie.c)
 *   2. Fast filter via Hash Table (hashtable.c)
 *   3. Search with Suffix Tree (suffix_tree.c)
 *   4. *** YOU ARE HERE — alignment when exact match fails ***
 *   5. Score is passed to Skip List for ranking (skiplist.c)
 *
 * =============================================================================
 * NOTE FOR main.c (do NOT change main.c — comments only):
 *
 *   Currently, main.c uses a simple character-by-character similarity function
 *   called calculate_similarity() defined inside main.c itself. That function
 *   does NOT use Needleman-Wunsch.
 *
 *   To fully integrate this module, the team lead (main.c owner) should:
 *     1. Add:   #include "alignment.h"   at the top of main.c
 *     2. In case 5 (exact search), when sequence is NOT found, call:
 *            needleman_wunsch(dataset[i], query)
 *            print_alignment(dataset[i], query)
 *        for each dataset sequence to show approximate matches.
 *     3. Optionally replace calculate_similarity() in case 7 with
 *        needleman_wunsch() for more accurate scoring.
 *
 *   Until that integration is done, this file compiles cleanly on its own
 *   and all functions work correctly as standalone calls.
 * =============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alignment.h"

/* ── Scoring constants for Needleman-Wunsch ── */
#define MATCH_SCORE     1   /* Score when two bases match (A==A, T==T, etc.) */
#define MISMATCH_SCORE -1   /* Penalty when bases differ  (A vs T, etc.)     */
#define GAP_PENALTY    -2   /* Penalty for inserting a gap '-'               */

/* Maximum length for a single DNA sequence */
#define MAX_LEN 1000

/* ── Helper: return the maximum of three integers ── */
static int max3(int a, int b, int c) {
    if (a >= b && a >= c) return a;
    if (b >= a && b >= c) return b;
    return c;
}

/*
 * needleman_wunsch(seq1, seq2)
 * ----------------------------
 * Performs GLOBAL alignment of two DNA sequences using dynamic programming.
 *
 * HOW IT WORKS:
 *   1. A 2D DP table (dp[i][j]) is filled where each cell stores the best
 *      alignment score for the first i characters of seq1 and first j
 *      characters of seq2.
 *   2. Three choices at each cell:
 *        a) Diagonal  → match or mismatch  (dp[i-1][j-1] + score)
 *        b) Left      → gap in seq1        (dp[i][j-1]   + GAP_PENALTY)
 *        c) Up        → gap in seq2        (dp[i-1][j]   + GAP_PENALTY)
 *   3. The bottom-right cell dp[len1][len2] is the final alignment score.
 *
 * RETURNS:
 *   Integer alignment score (higher = more similar).
 *   If both sequences are identical, score == length of sequence.
 *   Negative scores indicate very different sequences.
 */
int needleman_wunsch(const char* seq1, const char* seq2) {

    if (!seq1 || !seq2) {
        printf("[alignment] ERROR: NULL sequence passed.\n");
        return 0;
    }

    int len1 = (int)strlen(seq1);
    int len2 = (int)strlen(seq2);

    if (len1 == 0 || len2 == 0) {
        printf("[alignment] WARNING: One of the sequences is empty.\n");
        return 0;
    }

    if (len1 >= MAX_LEN || len2 >= MAX_LEN) {
        printf("[alignment] WARNING: Sequence too long (max %d). Truncating.\n", MAX_LEN - 1);
        len1 = len1 < MAX_LEN - 1 ? len1 : MAX_LEN - 1;
        len2 = len2 < MAX_LEN - 1 ? len2 : MAX_LEN - 1;
    }

    /*
     * Allocate the DP table dynamically.
     * dp[i][j] = best alignment score for seq1[0..i-1] vs seq2[0..j-1]
     * Size: (len1+1) x (len2+1)
     */
    int** dp = (int**)malloc((len1 + 1) * sizeof(int*));
    if (!dp) {
        printf("[alignment] ERROR: Memory allocation failed.\n");
        return 0;
    }
    for (int i = 0; i <= len1; i++) {
        dp[i] = (int*)malloc((len2 + 1) * sizeof(int));
        if (!dp[i]) {
            printf("[alignment] ERROR: Memory allocation failed at row %d.\n", i);
            /* Free already-allocated rows */
            for (int k = 0; k < i; k++) free(dp[k]);
            free(dp);
            return 0;
        }
    }

    /* ── Step 1: Initialize first row and column with gap penalties ── */
    for (int i = 0; i <= len1; i++) dp[i][0] = i * GAP_PENALTY;
    for (int j = 0; j <= len2; j++) dp[0][j] = j * GAP_PENALTY;

    /* ── Step 2: Fill the DP table ── */
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {

            /* Score for aligning seq1[i-1] with seq2[j-1] */
            int diag_score = (seq1[i - 1] == seq2[j - 1]) ? MATCH_SCORE : MISMATCH_SCORE;

            int from_diagonal = dp[i - 1][j - 1] + diag_score;  /* match/mismatch */
            int from_left     = dp[i][j - 1]     + GAP_PENALTY;  /* gap in seq1   */
            int from_up       = dp[i - 1][j]     + GAP_PENALTY;  /* gap in seq2   */

            dp[i][j] = max3(from_diagonal, from_left, from_up);
        }
    }

    int final_score = dp[len1][len2];

    /* ── Free memory ── */
    for (int i = 0; i <= len1; i++) free(dp[i]);
    free(dp);

    return final_score;
}


/*
 * print_alignment(seq1, seq2)
 * ---------------------------
 * Traces back through the DP table to reconstruct and print the alignment.
 *
 * OUTPUT FORMAT (3 lines):
 *   Seq1:  A C G T - A
 *          | | |   * |      ← '|' = match, '*' = mismatch, ' ' = gap
 *   Seq2:  A C G - T A
 *
 * HOW TRACEBACK WORKS:
 *   Starting from dp[len1][len2], move:
 *     - Diagonal → if this cell came from a match/mismatch
 *     - Left     → gap was inserted in seq1
 *     - Up       → gap was inserted in seq2
 *   Since we build aligned strings backwards, we reverse them at the end.
 */
void print_alignment(const char* seq1, const char* seq2) {

    if (!seq1 || !seq2) {
        printf("[alignment] ERROR: NULL sequence passed.\n");
        return;
    }

    int len1 = (int)strlen(seq1);
    int len2 = (int)strlen(seq2);

    if (len1 == 0 || len2 == 0) {
        printf("[alignment] WARNING: Cannot print alignment — empty sequence.\n");
        return;
    }

    if (len1 >= MAX_LEN || len2 >= MAX_LEN) {
        len1 = len1 < MAX_LEN - 1 ? len1 : MAX_LEN - 1;
        len2 = len2 < MAX_LEN - 1 ? len2 : MAX_LEN - 1;
    }

    /* ── Rebuild DP table for traceback ── */
    int** dp = (int**)malloc((len1 + 1) * sizeof(int*));
    if (!dp) {
        printf("[alignment] ERROR: Memory allocation failed.\n");
        return;
    }
    for (int i = 0; i <= len1; i++) {
        dp[i] = (int*)malloc((len2 + 1) * sizeof(int));
        if (!dp[i]) {
            for (int k = 0; k < i; k++) free(dp[k]);
            free(dp);
            return;
        }
    }

    for (int i = 0; i <= len1; i++) dp[i][0] = i * GAP_PENALTY;
    for (int j = 0; j <= len2; j++) dp[0][j] = j * GAP_PENALTY;

    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            int diag_score = (seq1[i - 1] == seq2[j - 1]) ? MATCH_SCORE : MISMATCH_SCORE;
            dp[i][j] = max3(
                dp[i - 1][j - 1] + diag_score,
                dp[i][j - 1]     + GAP_PENALTY,
                dp[i - 1][j]     + GAP_PENALTY
            );
        }
    }

    /*
     * ── Traceback ──
     * Aligned strings can be at most len1 + len2 characters long
     * (worst case: all gaps). We build them in reverse.
     */
    int max_aligned = len1 + len2 + 1;
    char* aligned1  = (char*)malloc(max_aligned * sizeof(char));
    char* aligned2  = (char*)malloc(max_aligned * sizeof(char));
    char* match_str = (char*)malloc(max_aligned * sizeof(char));

    if (!aligned1 || !aligned2 || !match_str) {
        printf("[alignment] ERROR: Memory allocation failed during traceback.\n");
        free(aligned1); free(aligned2); free(match_str);
        for (int i = 0; i <= len1; i++) free(dp[i]);
        free(dp);
        return;
    }

    int idx = 0;
    int i = len1, j = len2;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            int diag_score = (seq1[i - 1] == seq2[j - 1]) ? MATCH_SCORE : MISMATCH_SCORE;

            if (dp[i][j] == dp[i - 1][j - 1] + diag_score) {
                /* Came from diagonal → match or mismatch */
                aligned1[idx]  = seq1[i - 1];
                aligned2[idx]  = seq2[j - 1];
                match_str[idx] = (seq1[i - 1] == seq2[j - 1]) ? '|' : '*';
                i--; j--;
            } else if (dp[i][j] == dp[i][j - 1] + GAP_PENALTY) {
                /* Came from left → gap in seq1 */
                aligned1[idx]  = '-';
                aligned2[idx]  = seq2[j - 1];
                match_str[idx] = ' ';
                j--;
            } else {
                /* Came from up → gap in seq2 */
                aligned1[idx]  = seq1[i - 1];
                aligned2[idx]  = '-';
                match_str[idx] = ' ';
                i--;
            }
        } else if (i > 0) {
            /* seq2 exhausted — fill remaining with gaps */
            aligned1[idx]  = seq1[i - 1];
            aligned2[idx]  = '-';
            match_str[idx] = ' ';
            i--;
        } else {
            /* seq1 exhausted — fill remaining with gaps */
            aligned1[idx]  = '-';
            aligned2[idx]  = seq2[j - 1];
            match_str[idx] = ' ';
            j--;
        }
        idx++;
    }

    aligned1[idx]  = '\0';
    aligned2[idx]  = '\0';
    match_str[idx] = '\0';

    /* Reverse all three strings (they were built backwards) */
    int left = 0, right = idx - 1;
    while (left < right) {
        char tmp;
        tmp = aligned1[left];  aligned1[left]  = aligned1[right];  aligned1[right]  = tmp;
        tmp = aligned2[left];  aligned2[left]  = aligned2[right];  aligned2[right]  = tmp;
        tmp = match_str[left]; match_str[left] = match_str[right]; match_str[right] = tmp;
        left++; right--;
    }

    /* ── Print the alignment ── */
    int score = needleman_wunsch(seq1, seq2);
    printf("\n--- Needleman-Wunsch Alignment ---\n");
    printf("Score : %d\n", score);
    printf("Seq1  : %s\n", aligned1);
    printf("        %s\n", match_str);   /* '|' match, '*' mismatch, ' ' gap */
    printf("Seq2  : %s\n", aligned2);
    printf("Legend: | = match   * = mismatch   (space) = gap\n");
    printf("----------------------------------\n");

    /* ── Free all memory ── */
    free(aligned1);
    free(aligned2);
    free(match_str);
    for (int k = 0; k <= len1; k++) free(dp[k]);
    free(dp);
}
