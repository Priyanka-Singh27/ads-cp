#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "trie.h"
#include "suffix_tree.h"
#include "skiplist.h"
#include "alignment.h"    
#include "hashtable.h"    

#define MAX_SEQ 1000
#define MAX_DATASET 100
#define KMER 4

TrieNode* trie_root = NULL;

char dataset[MAX_DATASET][MAX_SEQ];
char species[MAX_DATASET][50];
char labels[MAX_DATASET][50];
int dataset_size = 0;

void show_menu() {
    printf("\n========================================\n");
    printf("        DNA MATCHING SYSTEM\n");
    printf("========================================\n");
    printf("1. Load DNA Dataset\n");
    printf("2. Enter DNA Query\n");
    printf("4. Exit\n");
    printf("5. Exact Match (Trie)\n");
    printf("6. Pattern Search (Suffix Tree)\n");
    printf("7. Rank & Analyze Matches\n");
    printf("========================================\n");
    printf("Enter choice: ");
}

/* ================= LOAD ================= */
void load_and_store(const char* filename, const char* species_name) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening %s\n", filename);
        return;
    }

    char line[MAX_SEQ];
    char current_label[50];

    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = 0;

        if (strlen(line) == 0) continue;

        if (line[0] == '>') {
            strcpy(current_label, line + 1);
            continue;
        }

        normalize_sequence(line);
        if (!validate_sequence(line)) continue;

        strcpy(dataset[dataset_size], line);
        strcpy(species[dataset_size], species_name);
        strcpy(labels[dataset_size], current_label);

        insert_sequence(trie_root, line);

        int len = strlen(line);
        if (len >= KMER) {
            for (int i = 0; i <= len - KMER; i++) {
                char kmer[KMER + 1];
                strncpy(kmer, &line[i], KMER);
                kmer[KMER] = '\0';
                insert_kmer(kmer);
            }
        }

        dataset_size++;
    }

    fclose(file);
}

/* ================= MAIN ================= */
int main() {
    int choice;
    char query[MAX_SEQ];

    while (1) {
        show_menu();
        scanf("%d", &choice);
        getchar();

        switch (choice) {

        case 1:
            printf("\nLoading dataset...\n");

            dataset_size = 0;

            if (trie_root) free_trie(trie_root);
            trie_root = create_trie();
            free_table();

            load_and_store("data/human.txt", "Human");
            load_and_store("data/chimpanzee.txt", "Chimpanzee");
            load_and_store("data/mouse.txt", "Mouse");
            load_and_store("data/virus_strain_a.txt", "Virus A");
            load_and_store("data/virus_strain_b.txt", "Virus B");

            printf("Dataset loaded: %d sequences\n", dataset_size);
            break;

        case 2:
            while (1) {
                printf("\nEnter DNA query: ");
                fgets(query, MAX_SEQ, stdin);
                query[strcspn(query, "\n")] = 0;

                normalize_sequence(query);

                if (!validate_sequence(query)) {
                    printf("Invalid DNA sequence. Try again.\n");
                } else {
                    printf("Normalized Query: %s\n", query);
                    break;
                }
            }
            break;

        case 4:
            printf("Exiting...\n");
            return 0;

        /* ================= TRIE ================= */
        case 5:
            if (!trie_root) {
                printf("Load dataset first.\n");
                break;
            }

            int found = 0;

            for (int i = 0; i < dataset_size; i++) {
                if (strcmp(dataset[i], query) == 0) {
                    printf("Exact sequence FOUND in: %s (%s)\n",
                           species[i], labels[i]);
                    found = 1;
                }
            }

            if (!found) {
                printf("Exact sequence NOT found.\n");

                printf("\nRunning Needleman-Wunsch alignment (Top 3):\n");

                int count = 0;
                for (int i = 0; i < dataset_size && count < 3; i++) {
                    if (strlen(dataset[i]) == 0) continue;

                    printf("\n[%s - %s]\n", species[i], labels[i]);
                    print_alignment(dataset[i], query);
                    count++;
                }
            }
            break;

        /* ================= SUFFIX TREE ================= */
        case 6: {
            if (dataset_size == 0) {
                printf("Load dataset first.\n");
                break;
            }

            int found = 0;

            printf("Searching pattern across dataset...\n");

            for (int i = 0; i < dataset_size; i++) {
                build_suffix_tree(dataset[i]);

                if (search_pattern(query)) {
                    printf("Pattern found in: %s (%s)\n",
                           species[i], labels[i]);

                    // 🔥 ALSO ALIGN
                    print_alignment(dataset[i], query);

                    found = 1;
                }

                free_suffix_tree();
            }

            if (!found) {
                printf("Pattern NOT found in dataset.\n");
            }

            break;
        }

        /* ================= RANKING ================= */
        case 7: {
            if (dataset_size == 0) {
                printf("Load dataset first.\n");
                break;
            }

            free_skiplist();
            init_skiplist();

            printf("\nAnalyzing DNA Sequence...\n");
            printf("----------------------------------------\n");

            int qlen = strlen(query);

            // 🔥 Hash filter FIXED
            if (qlen >= KMER) {
                char kmer[KMER + 1];
                strncpy(kmer, query, KMER);
                kmer[KMER] = '\0';

                if (!search_kmer(kmer)) {
                    printf("No similar sequences found (hash filter).\n");
                    break;
                }
            } else {
                printf("Query too short, skipping hash filter...\n");
            }

            int inserted = 0;

            for (int i = 0; i < dataset_size; i++) {
                if (strlen(dataset[i]) == 0) continue;

                int score = needleman_wunsch(dataset[i], query);
                insert_skiplist(dataset[i], species[i], score);
                inserted++;
            }

            if (inserted == 0) {
                printf("No matches found.\n");
                break;
            }

            printf("\nTOP MATCHES\n");
            printf("----------------------------------------\n");
            display_top_matches(5);

            break;
        }

        default:
            printf("Invalid choice.\n");
        }
    }

    return 0;
}