/******************************************************************
  Copyright 2006 by Michael Farrar.  All rights reserved.
  This program may not be sold or incorporated into a commercial product,
  in whole or in part, without written consent of Michael Farrar.  For 
  further information regarding permission for use or reproduction, please 
  contact: Michael Farrar at farrar.michael@gmail.com.
*******************************************************************/

/*
  Written by Michael Farrar, 2006.
  Please send bug reports and/or suggestions to farrar.michael@gmail.com.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "swsse2.h"
#include "fastalib.h"

#define READ_BUFFER_SIZE (300000000)
#define SEQ_NAME_SIZE    (128)


int sequences=0;
int residues=0;
extern int thread_num;
BUFFER_NODE* buffer_list=NULL;
BUFFER_NODE* last_buffer=NULL;
FASTA_LIB **openLib (char *file, int pad)
{
    FILE *fp;

    // FASTA_LIB *lib;
    int totalSize;
    char* globalBuffer = (char *) malloc (READ_BUFFER_SIZE);
    
    if (!globalBuffer) {
        fprintf (stderr, "Unable to allocate memory for read buffer\n");
        exit (-1);
    }

    if ((fp = fopen (file, "r")) == NULL) {
        fprintf (stderr, "Unable to open file %s\n", file);
        exit (-1);
    }

    totalSize = (int) fread (globalBuffer, sizeof (char), READ_BUFFER_SIZE, fp);
    if (totalSize == 0 && !feof (fp)) {
        fprintf (stderr, "Error (%d) reading fasta file\n", ferror (fp));
        exit (-1);
    }

    buffer_list=malloc(sizeof(BUFFER_NODE));
    buffer_list->buffer_ptr=globalBuffer;
    buffer_list->size=totalSize;
    buffer_list->next=NULL;
    buffer_list->direction=1;
    // lib = (FASTA_LIB *) malloc (sizeof (FASTA_LIB));
    // if (!lib) {
    //     fprintf (stderr, "Unable to allocate memory for library header\n");
    //     exit (-1);
    // }

    FASTA_LIB** buffer_local=malloc(sizeof(FASTA_LIB*)*thread_num);
    int chunkSize=totalSize/thread_num;
    for(int i=0;i<thread_num;i++){
        buffer_local[i]=malloc(sizeof(FASTA_LIB));
        buffer_local[i]->readBuffer=globalBuffer+i*chunkSize;
        buffer_local[i]->pos=0;
        buffer_local[i]->size=chunkSize;
        buffer_local[i]->fp=fp;
        buffer_local[i]->tid=i;
        buffer_local[i]->belong=buffer_list;
        buffer_local[i]->residues=0;
        buffer_local[i]->sequences=0;
    }

    // lib->seqBuffer = (unsigned char *) malloc (MAX_SEQ_LENGTH);
    // if (!lib->seqBuffer) {
    //     fprintf (stderr, "Unable to allocate memory for sequence\n");
    //     exit (-1);
    // }

    // lib->seqName = (char *) malloc (SEQ_NAME_SIZE);
    // if (!lib->seqName) {
    //     fprintf (stderr, "Unable to allocate memory for sequence name\n");
    //     exit (-1);
    // }

    
    // BUFFER_NODE* node=malloc(sizeof(BUFFER_NODE));
    // node->buffer_ptr=lib->readBuffer;
    // node->size=lib->size;
    
    // buffer_local=malloc(sizeof(BUFFER_LOCAL*)*thread_num);
    // int chunkSize=lib->size/thread_num;
    // for(int i=0;i<thread_num;i++){
    //     buffer_local[i]=malloc(sizeof(BUFFER_LOCAL));
    //     buffer_local[i]->start=lib->readBuffer+i*chunkSize;
    //     buffer_local[i]->belong=node;
    //     buffer_local[i]->size=chunkSize;
    // }
    // //最后一块的结尾
    // buffer_local[thread_num-1]->size=lib->size-(thread_num-1)*chunkSize;


    // lib->pos = 0;

    // lib->fp = fp;

    // lib->sequences = 0;
    // lib->residues = 0;

    // lib->pad;

    return buffer_local;
}

QUERY_LIB *openQueryLib (char *file, int pad)
{
    FILE *fp;

    QUERY_LIB *lib;
    
    if ((fp = fopen (file, "r")) == NULL) {
        fprintf (stderr, "Unable to open file %s\n", file);
        exit (-1);
    }

    lib = (QUERY_LIB *) malloc (sizeof (QUERY_LIB));
    if (!lib) {
        fprintf (stderr, "Unable to allocate memory for library header\n");
        exit (-1);
    }

    lib->readBuffer = (char *) malloc (1024*10*3);
    if (!lib->readBuffer) {
        fprintf (stderr, "Unable to allocate memory for read buffer\n");
        exit (-1);
    }

    lib->seqBuffer = (unsigned char *) malloc (MAX_SEQ_LENGTH);
    if (!lib->seqBuffer) {
        fprintf (stderr, "Unable to allocate memory for sequence\n");
        exit (-1);
    }

    lib->seqName = (char *) malloc (SEQ_NAME_SIZE);
    if (!lib->seqName) {
        fprintf (stderr, "Unable to allocate memory for sequence name\n");
        exit (-1);
    }

    lib->size = (int) fread (lib->readBuffer, sizeof (char), 1024*10*3, fp);
    if (lib->size == 0 && !feof (fp)) {
        fprintf (stderr, "Error (%d) reading fasta file\n", ferror (fp));
        exit (-1);
    }
    // BUFFER_NODE* node=malloc(sizeof(BUFFER_NODE));
    // node->buffer_ptr=lib->readBuffer;
    // node

    lib->pos = 0;

    lib->fp = fp;

    lib->sequences = 0;
    lib->residues = 0;

    lib->pad;

    return lib;
}

static int
readNextBlock (FASTA_LIB *lib)
{
    FILE *fp = lib->fp;
    
    static int direction=0;
    #pragma omp critical
    {
        struct timeb s1,e1;
        double startTimeSec=0;
    double endTimeSec=0;
        ftime(&s1);
        if(lib->belong->next==NULL&&!feof(fp)){
            int totalSize=0;
            char* globalBuffer=malloc(READ_BUFFER_SIZE);
            totalSize = fread (globalBuffer, sizeof (char), READ_BUFFER_SIZE, fp);
            if (lib->size == 0 && !feof (fp)) {
                fprintf (stderr, "Error (%d) reading fasta file\n", ferror (fp));
                exit (-1);
            }
            BUFFER_NODE * new_node=malloc(sizeof(BUFFER_NODE));
            new_node->buffer_ptr=globalBuffer;
            new_node->size=totalSize;
            new_node->direction=direction;
            direction=!direction;
            new_node->next=NULL;
            // if(totalSize>=READ_BUFFER_SIZE){
            //     lib->belong->next=new_node;

            // }else{
            //     last_buffer=new_node;
            //     lib->belong=new_node;
            //     lib->pos=0;
            //     lib->size=totalSize;
            //     lib->readBuffer=globalBuffer;
            //     self_flag=1;
            // }
            
            lib->belong->next=new_node;
        }
        ftime(&e1);
        TPRINT(lib->tid,s1,e1,"asd:");
        if(lib->belong->next)
            printf("%d",lib->belong->next->size);
    }
    if(lib->belong->next==NULL){
        lib->size=0;
        return 0;
    }
    lib->belong=lib->belong->next;
    int chunkSize=lib->belong->size/thread_num;
    if(lib->belong->direction)
        lib->readBuffer=lib->belong->buffer_ptr+lib->tid*chunkSize;
    else
        lib->readBuffer=lib->belong->buffer_ptr+(thread_num-1-lib->tid)*chunkSize;
    lib->pos=0;
    if((direction&&lib->tid==thread_num-1)||(!direction&&lib->tid==0))
        lib->size=lib->belong->size-(thread_num-1)*chunkSize;
    else
        lib->size=chunkSize;
    //不用处理最后一个线程的空间大小问题
    
    

    return lib->size;
}
static int
readNextQueryBlock (QUERY_LIB *lib)
{
    FILE *fp = lib->fp;
    size_t size;
    
    size = fread (lib->readBuffer, sizeof (char), 1024*10*3,fp);
    if (lib->size == 0 && !feof (fp)) {
        fprintf (stderr, "Error (%d) reading fasta file\n", ferror (fp));
        exit (-1);
    }

    lib->pos = 0;
    lib->size = (int) size;

    return lib->size;
}

unsigned char *
nextQuerySeq (QUERY_LIB *lib, int *length)
{
    int inx;
    int size;
    int done;
    int len;

    char *name = lib->seqName;
    unsigned char *seq = lib->seqBuffer;

    /* check if we are at the end of the library */
    if (lib->size == 0) {
        *length = 0;
        return NULL;
    }

    if (lib->pos == lib->size) {
        readNextQueryBlock (lib);
    }

    inx = lib->pos;

    /* check for the start of a sequence */
    if (lib->readBuffer[inx] != '>') {
        fprintf (stderr, "Error parsing fasta file expecting > found %c\n",
            lib->readBuffer[inx]);
        exit (-1);
    }

    ++inx;

    /* read in the sequence name */
    len = 0;
    done = 0;
    do {
        if (inx >= lib->size) {
            size = readNextQueryBlock (lib);
            if (size == 0) {
                *length = 0;
                return NULL;
            }
            inx = lib->pos;
        } else if (lib->readBuffer[inx] == '\n') {
            *name = '\0';
            done = 1;
        } else if (len < SEQ_NAME_SIZE - 1) {
            *name++ = lib->readBuffer[inx];
            len++;
        }
        ++inx;
    } while (!done);

    lib->pos = inx;

    /* read in the sequence */
    len = 0;
    done = 0;
    do {
        if (inx >= lib->size) {//处理超出了buffer的情况
            size = readNextQueryBlock (lib);
            if (size == 0) {
                *seq = '\0';
                done = 1;
            }
            inx = 0;
        } else if (isspace(lib->readBuffer[inx])) {
            ++inx;
        } else if (lib->readBuffer[inx] == '>') {//读到了下一条记录的开始，以>标记
            *seq = '\0';
            done = 1;
        } else if (len >= MAX_SEQ_LENGTH) {
            fprintf (stderr, "Sequence %s exceeds maximum length\n", 
                lib->seqName);
            exit (-1);
        } else {
            int value = AMINO_ACID_VALUE[lib->readBuffer[inx]];
            if (value == -1) {
                fprintf (stderr, "Unknown amino acid %c in sequence %s\n",
                    lib->readBuffer[inx], lib->seqName);
                exit (-1);
            }
            *seq++ = (char) value;
            inx++;
            len++;
        }
    } while (!done);

    lib->pos = inx;
    *length = len;

    lib->sequences++;
    lib->residues += len;

    /*  check if we need to pad the sequence to a multiple of 16  */
    if (lib->pad) {
        inx = 16 - (len % 16);
        while (inx--) {
            *seq++ = ALPHA_SIZE;
        }
        *seq = '\0';
    }

    return lib->seqBuffer;
}

unsigned char *
nextSeq (LIB_LOCAL *lib_local,FASTA_LIB *lib, int *length)
{
    int inx;
    int size;
    int done;
    int len;

    char *name = lib_local->seqName;
    unsigned char *seq = lib_local->seqBuffer;

    /* check if we are at the end of the library */
    if (lib->size == 0) {
        *length = 0;
        return NULL;
    }

    if (lib->pos == lib->size) {
        readNextBlock (lib);
    }

    inx = lib->pos;

    /* check for the start of a sequence */
    // if (lib->readBuffer[inx] != '>') {
    //     fprintf (stderr, "Error parsing fasta file expecting > found %c\n",
    //         lib->readBuffer[inx]);
    //     exit (-1);
    // }
    while(lib->readBuffer[inx] != '>'&&inx<lib->belong->size)
        inx++;
    if(inx>=lib->belong->size)
        return NULL;
    ++inx;

    /* read in the sequence name */
    len = 0;
    done = 0;
    do {
        if (inx >= lib->size) {
            //分配得到的是最后一块
            if((lib->tid==thread_num-1&&lib->belong->direction)||(lib->tid==0&&!lib->belong->direction)){
                size = readNextBlock (lib);
                if (size == 0) {
                    *length = 0;
                    return NULL;
                }
                inx = lib->pos;
            }else{
                while(lib->readBuffer[inx] != '\n'){
                    if(len>=SEQ_NAME_SIZE - 1){
                        inx++;
                        continue;
                    }
                    *name++ = lib->readBuffer[inx++];
                    len++;
                }
                *name='\0';
                done=1;
            }
        } else if (lib->readBuffer[inx] == '\n') {
            *name = '\0';
            done = 1;
        } else if (len < SEQ_NAME_SIZE - 1) {
            *name++ = lib->readBuffer[inx];
            len++;
        }
        ++inx;
    } while (!done);

    lib->pos = inx;

    /* read in the sequence */
    len = 0;
    done = 0;
    do {
        if (inx >= lib->size) {
            //处理超出了buffer的情况
            if((lib->tid==thread_num-1&&lib->belong->direction)||(lib->tid==0&&!lib->belong->direction)){
                size = readNextBlock (lib);
                if (size == 0) {
                    *seq = '\0';
                    done = 1;
                }
                inx = 0;
            }else{
                while(lib->readBuffer[inx] != '>'&&lib->readBuffer[inx] != '\0'){
                    if(isspace(lib->readBuffer[inx])){
                        inx++;
                        continue;
                    }
                    int value = AMINO_ACID_VALUE[lib->readBuffer[inx]];
                    if (value == -1) {
                        fprintf (stderr, "Unknown amino acid %c in sequence %s\n",
                            lib->readBuffer[inx], lib_local->seqName);
                        exit (-1);
                    }
                    *seq++ = (char) value;
                    inx++;
                    len++;
                }
                *seq='\0';
                done=1;
            }
        } else if (isspace(lib->readBuffer[inx])) {
            ++inx;
        } else if (lib->readBuffer[inx] == '>') {//读到了下一条记录的开始，以>标记
            *seq = '\0';
            done = 1;
        } else if (len >= MAX_SEQ_LENGTH) {
            fprintf (stderr, "Sequence %s exceeds maximum length\n", 
                lib_local->seqName);
            exit (-1);
        } else {
            int value = AMINO_ACID_VALUE[lib->readBuffer[inx]];
            if (value == -1) {
                fprintf (stderr, "Unknown amino acid %c in sequence %s\n",
                    lib->readBuffer[inx], lib_local->seqName);
                exit (-1);
            }
            *seq++ = (char) value;
            inx++;
            len++;
        }
    } while (!done);


    //处理超出读取的范围的问题
    if(inx>=lib->size){
        int size_next=readNextBlock(lib);
        inx=0;
    }
    lib->pos = inx;
    *length = len;

    lib->sequences++;
    lib->residues += len;

    /*  check if we need to pad the sequence to a multiple of 16  */
    if (lib->pad) {
        inx = 16 - (len % 16);
        while (inx--) {
            *seq++ = ALPHA_SIZE;
        }
        *seq = '\0';
    }

    return lib_local->seqBuffer;
}

void closeQueryLib (QUERY_LIB *lib)
{
    fclose (lib->fp);
    
    free (lib->readBuffer);
    // free (lib->seqBuffer);
    // free (lib->seqName);

    free (lib);
}
void closeLib (FASTA_LIB **lib)
{
    fclose (lib[0]->fp);
    while(buffer_list){
        if(buffer_list->buffer_ptr)
            free(buffer_list->buffer_ptr);
        BUFFER_NODE* n=buffer_list;
        buffer_list=buffer_list->next;
        free(n);
    }
    if(last_buffer){
        if(last_buffer->buffer_ptr)
            free(last_buffer->buffer_ptr);
        free(last_buffer);
    }
    free(lib);
    
    // free (lib->readBuffer);
    // // free (lib->seqBuffer);
    // // free (lib->seqName);

    // free (lib);
}
