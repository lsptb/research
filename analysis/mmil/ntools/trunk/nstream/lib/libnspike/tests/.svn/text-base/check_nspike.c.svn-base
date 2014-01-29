
#include <stdlib.h>
#include <check.h>
#include "libnspike.h"
#include "testmacros.h"

START_TEST(test_create)
{
    NSpike *nsp = NSpike_create();
    fail_unless(nsp != NULL);
    NSpike_set_master_address(nsp, "10.0.0.1");
    NSpike_set_dioport_count(nsp, 32);
    char tmp[256];
    NSpike_get_master_address(nsp, tmp, 256);
    FAIL_UNLESS_STRING_EQUAL("10.0.0.1", tmp);
    fail_unless(NSpike_get_num_dioports(nsp) == 32);
}
END_TEST

START_TEST(test_destroy)
{
    NSpike *nsp = NSpike_create();
    NSpike_destroy(nsp);
}
END_TEST

START_TEST(test_bad_ip)
{
    NSpike *nsp = NSpike_create();
    fail_unless(nsp != NULL);
    NSpike_set_master_address(nsp, "10.0.x.0");
    char tmp[256];
    NSpike_get_error(nsp, tmp, 256);
    FAIL_UNLESS_STRING_EQUAL("Master DSP: Unable to parse IP address", tmp);

    fail_unless(!NSpike_get_master_address(nsp, tmp, 256));
    NSpike_get_error(nsp, tmp, 256);
    FAIL_UNLESS_STRING_EQUAL("Master IP address not set", tmp);
}
END_TEST

START_TEST(test_configure)
{
    char tmp[256];

    NSpike *nsp;

    nsp = NSpike_create(); 
    NSpike_set_dioport_count(nsp, 32);
    fail_unless(NSpike_configure(nsp) == false);
    NSpike_get_error(nsp, tmp, 256);
    FAIL_UNLESS_STRING_EQUAL("Cannot configure hardware until master_address is set", tmp);

    NSpike_destroy(nsp);

    nsp = NSpike_create(); 
    NSpike_set_master_address(nsp, "10.0.0.1");
    fail_unless(NSpike_configure(nsp) == false);
    NSpike_get_error(nsp, tmp, 256);
    FAIL_UNLESS_STRING_EQUAL("Cannot configure hardware until dioports is set", tmp);

    NSpike_destroy(nsp);

    nsp = NSpike_create();
    NSpike_set_master_address(nsp, "10.1.2.10");
    NSpike_set_dioport_count(nsp, 32);
    fail_unless(NSpike_configure(nsp) == true);
    NSpike_destroy(nsp);

}
END_TEST

START_TEST(test_set_channel)
{
    NSpike *nsp;
    Channel *ch;

    ch = Channel_create();
    Channel_set_sample_rate(ch, 16000);

    /* add channel no dsps */
    nsp = NSpike_create();
    fail_if(NSpike_set_channel(nsp, ch));
    NSpike_destroy(nsp);

    /* add channel 2 dsps wrong sample rate */
    nsp = NSpike_create();
    NSpike_add_auxdsp(nsp, "10.1.2.11", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.12", 30000);
    fail_if(NSpike_set_channel(nsp, ch));
    NSpike_destroy(nsp);

    /* add channel 2 dsps one with right sample rate */
    nsp = NSpike_create();
    NSpike_add_auxdsp(nsp, "10.1.2.11", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.12", 16000);
    fail_unless(NSpike_set_channel(nsp, ch));
    NSpike_destroy(nsp);
}
END_TEST

START_TEST(test_get_num_auxdsps)
{
    NSpike *nsp = NSpike_create();
    NSpike_add_auxdsp(nsp, "10.1.2.11", 30000);
    NSpike_add_auxdsp(nsp, "10.1.2.12", 30000);
    fail_unless(NSpike_get_num_auxdsps(nsp) == 2);
    NSpike_destroy(nsp);
}
END_TEST

START_TEST(test_get_channel_by_num)
{
    Channel *ch = Channel_create();
    Channel_set_number(ch, 50);

    NSpike *nsp = NSpike_create();
    NSpike_add_auxdsp(nsp, "10.1.2.12", 30000);
    fail_unless(NSpike_set_channel(nsp, ch));
    fail_unless(ch == NSpike_get_channel_by_number(nsp, 50));
    NSpike_destroy(nsp);
}
END_TEST

START_TEST(test_change_channel_params)
{
    Channel *ch = Channel_create();
    Channel_set_number(ch, 1);

    NSpike *nsp = NSpike_create();
    NSpike_add_auxdsp(nsp, "10.1.2.11", 30000);
    fail_unless(NSpike_set_channel(nsp, ch));
    fail_unless(ch == NSpike_get_channel_by_number(nsp, 1));
    Channel_set_highpass_filter(ch, 3000);
    fail_unless(NSpike_apply_channel_changes(nsp, ch));
    NSpike_destroy(nsp);
}
END_TEST

Suite *nspike_suite(void)
{
    Suite *s = suite_create("NSpike");

    TCase *tc_core = tcase_create("core tests");
    tcase_add_test(tc_core, test_create);
    tcase_add_test(tc_core, test_destroy);
    tcase_add_test(tc_core, test_bad_ip);
    tcase_add_test(tc_core, test_configure);
    tcase_add_test(tc_core, test_set_channel);
    tcase_add_test(tc_core, test_get_num_auxdsps);
    tcase_add_test(tc_core, test_get_channel_by_num);
    suite_add_tcase(s, tc_core);

    return s;
}


