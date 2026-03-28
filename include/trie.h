#ifndef TRIE_H
#define TRIE_H

typedef struct TrieNode {
    struct TrieNode* children[4];
    int isEnd;
} TrieNode;

TrieNode* create_trie();
void insert_sequence(TrieNode* root, const char* seq);
int search_sequence(TrieNode* root, const char* seq);
void free_trie(TrieNode* root);

#endif