
#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include "libnspike.h"
#include "testmacros.h"

void open_test_files(FILE **in, FILE **out)
{
    int ret_val = 0;
    ret_val = system("openssl rand -out rand 1024000");
    fail_if(((ret_val >> 8) & 0xff) != 0);
    
    int bytes_read = 1;
    *in = fopen("rand", "r");
    fail_if(in == NULL);
    *out = fopen("rand2", "w");
    fail_if(out == NULL);
}

void compare_and_close_test_files(FILE **in, FILE **out)
{
    int ret_val = 0;
    fclose(*in);
    fclose(*out);
    ret_val = system("diff rand rand2");
    fail_if(((ret_val >> 8) & 0xff) != 0);
    unlink("rand");
    unlink("rand2");
}


START_TEST(test_create)
{
    RingBuffer *rb = RingBuffer_create(1000);
    fail_unless(rb != NULL);
    RingBuffer_destroy(rb);
}
END_TEST

START_TEST(test_destroy)
{
    RingBuffer *rb = RingBuffer_create(1000);
    RingBuffer_destroy(rb);
}
END_TEST

START_TEST(test_simple_overflow)
{
    char crap[1024];
    RingBuffer *rb = RingBuffer_create(1023);
    fail_if(RingBuffer_push_data(rb, crap, sizeof(crap)));
    RingBuffer_destroy(rb);
}
END_TEST

START_TEST(test_simple_push)
{
    char crap[1024];
    RingBuffer *rb = RingBuffer_create(1024);
    fail_if(RingBuffer_push_data(rb, crap, sizeof(crap)) == false);
    RingBuffer_destroy(rb);
}
END_TEST

START_TEST(test_overflow_on_wrap)
{
    char crap[1024];
    RingBuffer *rb = RingBuffer_create(512);
    printf("in buf: %d (0)\n",RingBuffer_data_available(rb));
    /* 512 free */
    printf("512 free -> push 256\n");
    fail_if(RingBuffer_push_data(rb, crap, 256) == false);
    printf("in buf: %d (256)\n",RingBuffer_data_available(rb));
    /* 256 free */
    printf("256 free -> get 256\n");
    fail_if(RingBuffer_get_data(rb, crap, 256) != 256);
    printf("in buf: %d (0)\n",RingBuffer_data_available(rb));
    /* 512 free */
    printf("512 free -> push 512\n");
    fail_if(RingBuffer_push_data(rb, crap, 512) == false);
    printf("in buf: %d (512)\n",RingBuffer_data_available(rb));
    /* 0 free */
    printf("0 free -> push 1\n");
    fail_if(RingBuffer_push_data(rb, crap, 1));
    RingBuffer_destroy(rb);
}
END_TEST

START_TEST(test_copy)
{
    FILE *in = NULL;
    FILE *out = NULL;

    char inbuf[2048];
    char outbuf[2048];
    int bytes_read = 1;
    
    open_test_files(&in, &out);
 
    RingBuffer *rb = RingBuffer_create(15000);
    int i = 0;
    int bytes_from_buf = 0;
    while(bytes_read > 0 || RingBuffer_data_available(rb) > 0)
    {
 //      printf("data available: %d\n", RingBuffer_data_available(rb));
       // fwrite(&inbuf, 1, bytes_read, out);
       bytes_read = fread(&inbuf, 1, 2048, in);
       if (bytes_read > 0)
        RingBuffer_push_data(rb, &inbuf, bytes_read);
        
       if (i > 4)
       {
        bytes_from_buf = RingBuffer_get_data(rb, &outbuf, 2048);
   //     printf("bytes read from buf: %d\n", bytes_from_buf);
        fwrite(&outbuf, 1, bytes_from_buf, out);
       } 
  //     printf("data available end: %d\n", RingBuffer_data_available(rb));
       i++;
    }
      
    RingBuffer_destroy(rb);
    
    compare_and_close_test_files(&in, &out);
}
END_TEST

START_TEST(test_copy_aligned)
{
    FILE *in = NULL;
    FILE *out = NULL;

    
    char inbuf[2048];
    char outbuf[4096];
    int bytes_read = 1;
    
    RingBuffer *rb = RingBuffer_create(32768);
    int i = 0;
    int bytes_from_buf = 0;
    open_test_files(&in, &out);
    while(bytes_read > 0 || RingBuffer_data_available(rb) > 0)
    {
       printf("data available: %d\n", RingBuffer_data_available(rb));
       // fwrite(&inbuf, 1, bytes_read, out);
       bytes_read = fread(&inbuf, 1, 2048, in);
       if (bytes_read > 0)
        RingBuffer_push_data(rb, &inbuf, bytes_read);
        
       if (i > 6)
       {
        bytes_from_buf = RingBuffer_get_data(rb, &outbuf, 2070);
        printf("bytes read from buf: %d\n", bytes_from_buf);
        fwrite(&outbuf, 1, bytes_from_buf, out);
       } 
       i++;
    }
      
    RingBuffer_destroy(rb);

    compare_and_close_test_files(&in, &out);
}
END_TEST

Suite *ringbuffer_suite(void)
{
    Suite *s = suite_create("RingBuffer");

    TCase *tc_core = tcase_create("core tests");
    tcase_add_test(tc_core, test_create);
    tcase_add_test(tc_core, test_destroy);
    tcase_add_test(tc_core, test_simple_overflow);
    tcase_add_test(tc_core, test_simple_push);
    tcase_add_test(tc_core, test_overflow_on_wrap);
    tcase_add_test(tc_core, test_copy);
    tcase_add_test(tc_core, test_copy_aligned);
    suite_add_tcase(s, tc_core);

    return s;
}


