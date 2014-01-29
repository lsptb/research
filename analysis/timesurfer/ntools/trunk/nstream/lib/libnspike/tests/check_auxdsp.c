
#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include "libnspike.h"
#include "testmacros.h"

#define AUX_DSP_IP "10.1.2.11"

START_TEST(test_create)
{
    AuxDSP *adsp = AuxDSP_create(0);
    fail_unless(adsp != NULL);
    AuxDSP_destroy(adsp);
}
END_TEST

START_TEST(test_destroy)
{
    AuxDSP *adsp = AuxDSP_create(0);
    AuxDSP_destroy(adsp);
}
END_TEST

START_TEST(test_connect)
{
    AuxDSP *adsp = AuxDSP_create(0);
    fail_unless(AuxDSP_connect(adsp, AUX_DSP_IP));
    AuxDSP_destroy(adsp);
}
END_TEST

START_TEST(test_getset_sample_rate)
{
    /* first verify default is 30k */
    AuxDSP *adsp = AuxDSP_create(0);
    fail_unless(AuxDSP_get_sample_rate(adsp) == 30000);
    AuxDSP_destroy(adsp);
   
    /* then verify that we can change it */ 
    adsp = AuxDSP_create(0);
    fail_unless(AuxDSP_set_sample_rate(adsp, 16000));
    fail_unless(AuxDSP_get_sample_rate(adsp) == 16000);
    AuxDSP_destroy(adsp);
}
END_TEST

START_TEST(test_initialize)
{
    AuxDSP *adsp = AuxDSP_create(0);
    fail_unless(AuxDSP_connect(adsp, AUX_DSP_IP));
    fail_unless(AuxDSP_initialize(adsp));
    AuxDSP_destroy(adsp);
}
END_TEST

START_TEST(test_set_channel)
{
    AuxDSP *adsp = AuxDSP_create(0);
    Channel *ch = Channel_create();
    fail_unless(AuxDSP_set_channel(adsp, ch));
    fail_unless(AuxDSP_get_channel_count(adsp) == 1);
    AuxDSP_destroy(adsp);
}
END_TEST

START_TEST(test_start_acquire)
{
    AuxDSP *adsp = AuxDSP_create(0);
    Channel *ch = Channel_create();
    fail_unless(AuxDSP_set_channel(adsp, ch));
    fail_unless(AuxDSP_connect(adsp, AUX_DSP_IP));
    fail_unless(AuxDSP_initialize(adsp));
    fail_unless(AuxDSP_start_acquire(adsp));
    AuxDSP_destroy(adsp);
}
END_TEST

START_TEST(test_stop_acquire)
{
    AuxDSP *adsp = AuxDSP_create(0);
    Channel *ch = Channel_create();
    fail_unless(AuxDSP_set_channel(adsp, ch));
    fail_unless(AuxDSP_connect(adsp, AUX_DSP_IP));
    fail_unless(AuxDSP_initialize(adsp));
    fail_unless(AuxDSP_stop_acquire(adsp));
    AuxDSP_destroy(adsp);
}
END_TEST

START_TEST(test_get_channel_by_num)
{
    AuxDSP *adsp = AuxDSP_create(0);
    Channel *ch = Channel_create();
    Channel_set_number(ch, 50);
    Channel_set_hw_number(ch, 150);
    fail_unless(AuxDSP_set_channel(adsp, ch));
    fail_unless(AuxDSP_get_channel_count(adsp) == 1);
    fail_unless(ch == AuxDSP_get_channel_by_number(adsp, 50));
    AuxDSP_destroy(adsp);
}
END_TEST

START_TEST(test_auxdsp_get_num)
{
    AuxDSP *adsp = AuxDSP_create(16);
    fail_unless(AuxDSP_get_number(adsp) == 16);
    AuxDSP_destroy(adsp);
}
END_TEST

Suite *auxdsp_suite(void)
{
    Suite *s = suite_create("AuxDSP");

    TCase *tc_core = tcase_create("core tests");
    suite_add_tcase(s, tc_core);
    tcase_add_test(tc_core, test_create);
    tcase_add_test(tc_core, test_destroy);
    tcase_add_test(tc_core, test_connect);
    tcase_add_test(tc_core, test_getset_sample_rate);
    tcase_add_test(tc_core, test_set_channel);
    tcase_add_test(tc_core, test_initialize);
    tcase_add_test(tc_core, test_start_acquire);
    tcase_add_test(tc_core, test_stop_acquire);
    tcase_add_test(tc_core, test_get_channel_by_num);
    tcase_add_test(tc_core, test_auxdsp_get_num);

    return s;
}


