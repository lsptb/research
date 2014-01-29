
#include "libnspike.h"
#include <stdlib.h>
#include <check.h>
#include "suites.h"

Suite *libnspike_suite(void)
{
    Suite *s = suite_create("libnspike");

    TCase *tc_core = tcase_create("all tests");
    suite_add_tcase(s, tc_core);

    return s;
}

main()
{
    int number_failed;
    Suite *s = libnspike_suite();
    SRunner *sr = srunner_create(s);
    srunner_add_suite(sr, ringbuffer_suite());
    // srunner_add_suite(sr, channel_suite());
    // srunner_add_suite(sr, auxdsp_suite());
    srunner_add_suite(sr, nspike_suite());
    // srunner_add_suite(sr, masterdsp_suite());
    // srunner_add_suite(sr, daqengine_suite());
    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

