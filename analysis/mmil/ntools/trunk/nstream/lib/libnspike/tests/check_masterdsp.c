
#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include "libnspike.h"
#include "testmacros.h"
#include "fakehw.h"

#define MASTER_DSP_IP "10.1.2.1"

/* this is the conversation that occurs when masterdsp_initialize is called
 * it has been placed here to make the tests of specific functions that
 * require that masterdsp_initialize be called shorter */
void setup_masterdsp_and_fakehw(fakehw *f, MasterDSP *mdsp)
{
    fail_unless(MasterDSP_set_dio_port_count(mdsp, 4));
    fail_unless(MasterDSP_set_dio_mode(mdsp, 0, true));
    fail_unless(MasterDSP_set_dio_mode(mdsp, 1, false));
    fail_unless(MasterDSP_set_dio_mode(mdsp, 2, true));
    fail_unless(MasterDSP_set_dio_mode(mdsp, 3, true));

    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd100, 0x0200, 0x0032}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd300, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd400, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xcc00, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd200, 0x0200, 0x0032}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd500, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd600, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xcd00, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd900, 0x0200, 0x0d00}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xda00, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xdc00, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xdd00, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xe200, 0x0200, 0xffff}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xe300, 0x0200, 0xffff}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xe400, 0x0200, 0xffff}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xe500, 0x0200, 0xffff}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xe800, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xed00, 0x0200, 0x0100}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
}

START_TEST(test_create)
{
    MasterDSP *mdsp = MasterDSP_create();
    fail_unless(mdsp != NULL);
    MasterDSP_destroy(mdsp);
}
END_TEST

START_TEST(test_destroy)
{
    MasterDSP *mdsp = MasterDSP_create();
    MasterDSP_destroy(mdsp);
}
END_TEST

START_TEST(test_connect)
{
    MasterDSP *mdsp = MasterDSP_create();
    fail_unless(MasterDSP_connect(mdsp, MASTER_DSP_IP));
    MasterDSP_destroy(mdsp);
}
END_TEST

START_TEST(test_initialize)
{
       
    fakehw *f = fakehw_create("10.1.2.1", NSPIKE_MESSAGE_PORT);
    MasterDSP *mdsp = MasterDSP_create();
    
    fakehw_start_seq(f,"MasterDSP_initialize()");
    setup_masterdsp_and_fakehw(f, mdsp);
    fakehw_run(f);

    fail_unless(MasterDSP_connect(mdsp, MASTER_DSP_IP));
    fail_unless(MasterDSP_initialize(mdsp));    
    MasterDSP_destroy(mdsp);
    fail_unless(fakehw_check(f), f->failreason);
    fakehw_destroy(f);
}
END_TEST

START_TEST(test_dac_gain)
{
    fakehw *f = fakehw_create("10.1.2.1", NSPIKE_MESSAGE_PORT);
    MasterDSP *mdsp = MasterDSP_create();

    fakehw_start_seq(f,"test_dac_gain");
    setup_masterdsp_and_fakehw(f, mdsp);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd100, 0x0200, 0x0032}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd200, 0x0200, 0x0019}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_run(f);

    fail_unless(MasterDSP_connect(mdsp, MASTER_DSP_IP));
    fail_unless(MasterDSP_initialize(mdsp));
    fail_unless(MasterDSP_set_dac_gain(mdsp, 0, 100));
    fail_unless(MasterDSP_set_dac_gain(mdsp, 1, 50));
    MasterDSP_destroy(mdsp);
    fail_unless(fakehw_check(f), f->failreason);
    fakehw_destroy(f);
}
END_TEST

START_TEST(test_dac_cutoff)
{
    fakehw *f = fakehw_create("10.1.2.1", NSPIKE_MESSAGE_PORT);
    MasterDSP *mdsp = MasterDSP_create();

    fakehw_start_seq(f,"test_dac_cutoff");
    setup_masterdsp_and_fakehw(f, mdsp);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd300, 0x0200, 0x6400}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd400, 0x0200, 0x9cff}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd500, 0x0200, 0x3200}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd600, 0x0200, 0xceff}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_run(f);

    fail_unless(MasterDSP_connect(mdsp, MASTER_DSP_IP));
    fail_unless(MasterDSP_initialize(mdsp));
    fail_unless(MasterDSP_set_dac_cutoff(mdsp, 0, 100));
    fail_unless(MasterDSP_set_dac_cutoff(mdsp, 1, 50));
    MasterDSP_destroy(mdsp);
    fail_unless(fakehw_check(f), f->failreason);
    fakehw_destroy(f);
}
END_TEST


START_TEST(test_dac_mute)
{
    fakehw *f = fakehw_create("10.1.2.1", NSPIKE_MESSAGE_PORT);
    MasterDSP *mdsp = MasterDSP_create();

    fakehw_start_seq(f,"test_dac_gain");
    setup_masterdsp_and_fakehw(f, mdsp);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd100, 0x0200, 0x0032}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd200, 0x0200, 0x0019}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd100, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd200, 0x0200, 0x0000}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd100, 0x0200, 0x0032}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xd200, 0x0200, 0x0019}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_run(f);

    fail_unless(MasterDSP_connect(mdsp, MASTER_DSP_IP));
    fail_unless(MasterDSP_initialize(mdsp));
    fail_unless(MasterDSP_set_dac_gain(mdsp, 0, 100));
    fail_unless(MasterDSP_set_dac_gain(mdsp, 1, 50));
    fail_unless(MasterDSP_set_dac_mute(mdsp, 0, true));
    fail_unless(MasterDSP_set_dac_mute(mdsp, 1, true));
    fail_unless(MasterDSP_set_dac_mute(mdsp, 0, false));
    fail_unless(MasterDSP_set_dac_mute(mdsp, 1, false));
    MasterDSP_destroy(mdsp);
    fail_unless(fakehw_check(f), f->failreason);
    fakehw_destroy(f);
}
END_TEST

START_TEST(test_dac_delay)
{
    fakehw *f = fakehw_create("10.1.2.1", NSPIKE_MESSAGE_PORT);
    MasterDSP *mdsp = MasterDSP_create();

    fakehw_start_seq(f,"test_dac_delay");
    setup_masterdsp_and_fakehw(f, mdsp);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xcc00, 0x0200, 0x3c00}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_recv(f,(unsigned short[]){0xfbff, 0x8200, 0x0040, 0xcd00, 0x0200, 0x5a00}, 12);
    fakehw_send(f,(unsigned short[]){0xfbff, 0x0100}, 4);
    fakehw_run(f);

    fail_unless(MasterDSP_connect(mdsp, MASTER_DSP_IP));
    fail_unless(MasterDSP_initialize(mdsp));
    fail_unless(MasterDSP_set_dac_delay(mdsp, 0, 2));
    fail_unless(MasterDSP_set_dac_delay(mdsp, 1, 3));
    MasterDSP_destroy(mdsp);
    fail_unless(fakehw_check(f), f->failreason);
    fakehw_destroy(f);
}
END_TEST

Suite *masterdsp_suite(void)
{
    Suite *s = suite_create("MasterDSP");

    TCase *tc_core = tcase_create("core tests");
    suite_add_tcase(s, tc_core);
    tcase_add_test(tc_core, test_create);
    tcase_add_test(tc_core, test_destroy);
    tcase_add_test(tc_core, test_connect);
    tcase_add_test(tc_core, test_initialize);
    tcase_add_test(tc_core, test_dac_gain);
    tcase_add_test(tc_core, test_dac_mute);
    tcase_add_test(tc_core, test_dac_cutoff);
    tcase_add_test(tc_core, test_dac_delay);

    return s;
}


