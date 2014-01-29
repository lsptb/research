
#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include "libnspike.h"
#include "testmacros.h"

START_TEST(test_create)
{
    Channel *ch = Channel_create();
    fail_unless(ch != NULL);
    Channel_destroy(ch);
}
END_TEST

START_TEST(test_destroy)
{
    Channel *ch = Channel_create();
    Channel_destroy(ch);
}
END_TEST

START_TEST(test_getset_number)
{
    Channel *ch = Channel_create();
    Channel_set_number(ch, 100);
    fail_unless(Channel_get_number(ch) == 100);
    Channel_destroy(ch);
}
END_TEST

START_TEST(test_getset_cal)
{
    Channel *ch = Channel_create();
    Channel_set_calibration_factor(ch, 100);
    fail_unless(Channel_get_calibration_factor(ch) == 100);
    Channel_destroy(ch);
}
END_TEST

START_TEST(test_getset_hw_number)
{
    Channel *ch = Channel_create();
    Channel_set_hw_number(ch, 100);
    fail_unless(Channel_get_hw_number(ch) == 100);
    Channel_destroy(ch);
}
END_TEST

START_TEST(test_getset_ref_channel)
{
    Channel *ch = Channel_create();
    Channel *ref_ch = Channel_create();
    Channel_set_hw_number(ref_ch, 100);
    Channel_set_reference_channel(ch, ref_ch);
    ref_ch = NULL;
    ref_ch = Channel_get_reference_channel(ch);
    fail_unless(Channel_get_hw_number(ref_ch) == 100);
    Channel_destroy(ref_ch);

    /* now test using a NULL for ground */
    Channel_set_reference_channel(ch, NULL);
    ref_ch = NULL;
    ref_ch = Channel_get_reference_channel(ch);
    fail_unless(ref_ch == NULL); 
    Channel_destroy(ch);
}
END_TEST

START_TEST(test_getset_highpass)
{
    Channel *ch = Channel_create();
    Channel_set_highpass_filter(ch, 100);
    fail_unless(Channel_get_highpass_filter(ch) == 100);
    Channel_destroy(ch);
}
END_TEST

START_TEST(test_getset_lowpass)
{
    Channel *ch = Channel_create();
    Channel_set_lowpass_filter(ch, 100);
    fail_unless(Channel_get_lowpass_filter(ch) == 100);
    Channel_destroy(ch);
}
END_TEST

START_TEST(test_getset_sample_rate)
{
    /* first check for reasonable default */
    Channel *ch = Channel_create();
    fail_unless(Channel_get_sample_rate(ch) == 30000);
    Channel_destroy(ch);
   
    /* then try setting it */ 
    ch = Channel_create();
    Channel_set_sample_rate(ch, 16000);
    fail_unless(Channel_get_sample_rate(ch) == 16000);
    Channel_destroy(ch);
}
END_TEST

START_TEST(test_getset_auxdsp)
{
    AuxDSP *adsp = AuxDSP_create(0);
    Channel *ch = Channel_create();

    Channel_set_auxdsp(ch, adsp);
    fail_unless(adsp == Channel_get_auxdsp(ch));
}
END_TEST

START_TEST(test_getset_slotnum)
{
    Channel *ch = Channel_create();
    Channel_set_dsp_slot_num(ch, 100);
    fail_unless(Channel_get_dsp_slot_num(ch) == 100);
    Channel_destroy(ch);
}
END_TEST


Suite *channel_suite(void)
{
    Suite *s = suite_create("Channel");

    TCase *tc_core = tcase_create("core tests");
    suite_add_tcase(s, tc_core);
    tcase_add_test(tc_core, test_create);
    tcase_add_test(tc_core, test_destroy);
    tcase_add_test(tc_core, test_getset_number);
    tcase_add_test(tc_core, test_getset_hw_number);
    tcase_add_test(tc_core, test_getset_cal);
    tcase_add_test(tc_core, test_getset_ref_channel);
    tcase_add_test(tc_core, test_getset_highpass);
    tcase_add_test(tc_core, test_getset_lowpass);
    tcase_add_test(tc_core, test_getset_sample_rate);
    tcase_add_test(tc_core, test_getset_auxdsp);
    tcase_add_test(tc_core, test_getset_slotnum);

    return s;
}


