#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "utils.h"

int validate_sequence(const char* seq) {
    for (int i = 0; seq[i] != '\0'; i++) {
        char c = toupper(seq[i]);
        if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
            return 0; // invalid
        }
    }
    return 1;
}

void normalize_sequence(char* seq) {
    for (int i = 0; seq[i] != '\0'; i++) {
        seq[i] = toupper(seq[i]);
    }
}

void load_sample_data(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file.\n");
        return;
    }

    char line[256];
    while (fgets(line, sizeof(line), file)) {
        printf("%s", line); // later: pass to trie/suffix tree
    }

    fclose(file);
}