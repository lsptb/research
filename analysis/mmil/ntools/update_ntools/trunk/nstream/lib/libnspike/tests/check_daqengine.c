
#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include "libnspike.h"
#include "testmacros.h"

#define AUX_DSP_IP "10.1.2.11"

START_TEST(test_create)
{
    DAQEngine *daq = DAQEngine_create();
    fail_unless(daq != NULL);
    DAQEngine_destroy(daq);
}
END_TEST

START_TEST(test_destroy)
{
    DAQEngine *daq = DAQEngine_create();
    DAQEngine_destroy(daq);
}
END_TEST

#if 0
START_TEST(test_connect)
{
    DAQEngine *daq = DAQEngine_create(0);
    fail_unless(DAQEngine_connect(daq, AUX_DSP_IP));
    DAQEngine_destroy(daq);
}
END_TEST

START_TEST(test_getset_sample_rate)
{
    /* first verify default is 30k */
    DAQEngine *daq = DAQEngine_create(0);
    fail_unless(DAQEngine_get_sample_rate(daq) == 30000);
    DAQEngine_destroy(daq);
   
    /* then verify that we can change it */ 
    daq = DAQEngine_create(0);
    fail_unless(DAQEngine_set_sample_rate(daq, 16000));
    fail_unless(DAQEngine_get_sample_rate(daq) == 16000);
    DAQEngine_destroy(daq);
}
END_TEST

START_TEST(test_initialize)
{
    DAQEngine *daq = DAQEngine_create();
    fail_unless(DAQEngine_connect(daq, AUX_DSP_IP));
    fail_unless(DAQEngine_initialize(daq));
    DAQEngine_destroy(daq);
}
END_TEST

START_TEST(test_add_channel)
{
    DAQEngine *daq = DAQEngine_create(0);
    Channel *ch = Channel_create();
    fail_unless(DAQEngine_add_channel(daq, ch));
    fail_unless(DAQEngine_get_channel_count(daq) == 1);
    DAQEngine_destroy(daq);
}
END_TEST

START_TEST(test_start_acquire)
{
    DAQEngine *daq = DAQEngine_create(0);
    Channel *ch = Channel_create();
    fail_unless(DAQEngine_add_channel(daq, ch));
    fail_unless(DAQEngine_connect(daq, AUX_DSP_IP));
    fail_unless(DAQEngine_initialize(daq));
    fail_unless(DAQEngine_start_acquire(daq));
    DAQEngine_destroy(daq);
}
END_TEST

START_TEST(test_stop_acquire)
{
    DAQEngine *daq = DAQEngine_create(0);
    Channel *ch = Channel_create();
    fail_unless(DAQEngine_add_channel(daq, ch));
    fail_unless(DAQEngine_connect(daq, AUX_DSP_IP));
    fail_unless(DAQEngine_initialize(daq));
    fail_unless(DAQEngine_stop_acquire(daq));
    DAQEngine_destroy(daq);
}
END_TEST

START_TEST(test_get_channel_by_num)
{
    DAQEngine *daq = DAQEngine_create(0);
    Channel *ch = Channel_create();
    Channel_set_number(ch, 50);
    Channel_set_hw_number(ch, 150);
    fail_unless(DAQEngine_add_channel(daq, ch));
    fail_unless(DAQEngine_get_channel_count(daq) == 1);
    fail_unless(ch == DAQEngine_get_channel_by_number(daq, 50));
    DAQEngine_destroy(daq);
}
END_TEST

START_TEST(test_auxdsp_get_num)
{
    DAQEngine *daq = DAQEngine_create(16);
    fail_unless(DAQEngine_get_number(daq) == 16);
    DAQEngine_destroy(daq);
}
END_TEST

#endif

Suite *daqengine_suite(void)
{
    Suite *s = suite_create("DAQEngine");

    TCase *tc_core = tcase_create("core tests");
    suite_add_tcase(s, tc_core);
    tcase_add_test(tc_core, test_create);
    tcase_add_test(tc_core, test_destroy);

    return s;
}


