#ifndef HASHTABLE_H
#define HASHTABLE_H

#define TABLE_SIZE 1000

typedef struct HashNode {
    char* kmer;
    struct HashNode* next;
} HashNode;

void insert_kmer(const char* kmer);
int search_kmer(const char* kmer);
void free_table();

#endif