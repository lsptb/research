loader = test_loader();
runner = text_test_runner(1, 1);
suite = test_suite();
suite = set_name(suite, 'All tests');
suite = add_test(suite, load_tests_from_test_case(loader, 'test_nspike_to_ieeg'));
suite = add_test(suite, load_tests_from_test_case(loader, 'test_ieeg_to_epoch'));

run(runner, suite);