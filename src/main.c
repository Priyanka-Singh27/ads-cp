#include <stdio.h>
#include <string.h>
#include "utils.h"

#define MAX_SEQ 1000

void show_menu() {
    printf("\n==== DNA Matching System ====\n");
    printf("1. Load Sample Data\n");
    printf("2. Enter DNA Sequence\n");
    printf("3. Validate Sequence\n");
    printf("4. Exit\n");
    printf("Choose an option: ");
}

int main() {
    int choice;
    char sequence[MAX_SEQ];

    while (1) {
        show_menu();
        scanf("%d", &choice);
        getchar(); // consume newline

        switch (choice) {

            case 1:
                printf("\nLoading sample data...\n");
                load_sample_data("data/human.txt");
                load_sample_data("data/chimpanzee.txt");
                load_sample_data("data/mouse.txt");
                load_sample_data("data/virus_strain_a.txt");
                load_sample_data("data/virus_strain_b.txt");
                break;

            case 2:
                printf("\nEnter DNA sequence: ");
                fgets(sequence, MAX_SEQ, stdin);
                sequence[strcspn(sequence, "\n")] = 0;

                normalize_sequence(sequence);
                printf("Normalized Sequence: %s\n", sequence);
                break;

            case 3:
                if (validate_sequence(sequence)) {
                    printf("Valid DNA sequence.\n");
                } else {
                    printf("Invalid DNA sequence.\n");
                }
                break;

            case 4:
                printf("Exiting...\n");
                return 0;

            default:
                printf("Invalid choice.\n");
        }
    }

    return 0;
}